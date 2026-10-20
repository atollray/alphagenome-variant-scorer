"""Bounded parallel workers, failure draining and cross-run deduplication."""

from concurrent.futures import ThreadPoolExecutor
from contextlib import closing
import json
from pathlib import Path
import sqlite3
import tempfile
from threading import Barrier, Lock
import time
import unittest
from unittest.mock import Mock, patch

from filelock import FileLock

from alphagenome_variant_scorer.api import FetchError, ScoreClient
from alphagenome_variant_scorer.cache import SharedCache, identity
from alphagenome_variant_scorer.errors import UserError
from alphagenome_variant_scorer.prepare import prepare
from alphagenome_variant_scorer.reference import ensure_index
from alphagenome_variant_scorer.runner import run
from alphagenome_variant_scorer.scores import RESPONSE_FORMAT
from alphagenome_variant_scorer.store import run_status


class ParallelTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.fasta = self.root / "synthetic.fa"
        self.fasta.write_text(">chr1\n" + "ACGT" * 64 + "\n", encoding="ascii")
        self.index = self.root / "synthetic.fai"
        ensure_index(self.fasta, self.index, create=True)
        self.source = self.root / "input.tsv"
        self.source.write_text("#CHROM\tPOS\tREF\tALT\n" + "".join(
            f"chr1\t{p}\t{'ACGT'[(p-1)%4]}\t{'TGCA'[(p-1)%4]}\n"
            for p in range(1, 9)), encoding="ascii")
        self.cache = SharedCache(self.root / "cache")
        self.cache.initialize()

    def prepared(self, name):
        path = self.root / name
        prepare(self.source, path, self.fasta, self.index, input_format="literal-tsv")
        return path

    def execute(self, path, fetch, **kwargs):
        return run(path, self.cache.directory, fetch, min_interval=kwargs.pop("min_interval", 0),
                   retry_delay=0, **kwargs)

    def test_parallelism_bound_and_max_jobs(self):
        path = self.prepared("run")
        barrier, lock = Barrier(4), Lock()
        active = peak = calls = 0

        def fetch(request):
            nonlocal active, peak, calls
            with lock:
                active += 1
                calls += 1
                peak = max(peak, active)
            # Four calls must actually be in flight, not just queued serially.
            barrier.wait(timeout=5)
            with closing(sqlite3.connect(path / "run.private.sqlite")) as db:
                self.assertLessEqual(db.execute(
                    "SELECT count(*) FROM jobs WHERE state='running'").fetchone()[0], 4)
                self.assertGreaterEqual(db.execute(
                    "SELECT sum(attempts) FROM job_results").fetchone()[0], 4)
            with lock:
                active -= 1
            return b"synthetic", {}

        result = self.execute(path, fetch, workers=4, max_jobs=4)
        self.assertEqual((peak, calls), (4, 4))
        self.assertEqual(result["state"], "paused")
        self.assertEqual(result["scoring"]["fetched_jobs"], 4)
        resumed = self.execute(path, lambda _: (b"synthetic", {}), workers=2)
        self.assertEqual(resumed["state"], "completed")

    def test_rate_gate_includes_all_workers_and_retries(self):
        path = self.prepared("run")
        starts, counts, lock = [], {}, Lock()

        def fetch(request):
            pos = request["variant"]["pos"]
            with lock:
                starts.append(time.monotonic())
                counts[pos] = counts.get(pos, 0) + 1
                attempt = counts[pos]
            if attempt == 1:
                raise FetchError("unavailable")
            return b"synthetic", {}

        result = self.execute(path, fetch, workers=4, max_jobs=4, min_interval=0.025)
        self.assertEqual(result["scoring"]["fetch_attempts"], 8)
        self.assertEqual(len(starts), 8)
        self.assertTrue(all(b-a >= 0.020 for a, b in zip(starts, starts[1:])))

    def test_32_workers_can_finish_a_bounded_batch(self):
        self.source.write_text("#CHROM\tPOS\tREF\tALT\n" + "".join(
            f"chr1\t{p}\t{'ACGT'[(p-1)%4]}\t{'TGCA'[(p-1)%4]}\n"
            for p in range(1, 41)), encoding="ascii")
        path = self.prepared("run")
        barrier, lock = Barrier(32), Lock()
        calls = []

        def fetch(request):
            with lock:
                calls.append(1)
            barrier.wait(timeout=10)
            return b"synthetic", {}

        result = self.execute(path, fetch, workers=32, max_jobs=32)
        self.assertEqual(len(calls), 32)
        self.assertEqual(result["state"], "paused")
        self.assertEqual(result["scoring"]["fetched_jobs"], 32)
        self.assertEqual(result["runtime"]["workers"], 32)

    def failure_pair(self, interruption):
        path = self.prepared("run")
        barrier, lock = Barrier(2), Lock()
        calls = []

        def fetch(request):
            with lock:
                calls.append(1)
                number = len(calls)
            barrier.wait(timeout=5)
            if number == 1:
                raise interruption
            time.sleep(0.05)
            return b"saved-in-flight", {}

        return path, calls, fetch

    def test_fatal_error_stops_new_calls_and_keeps_inflight_success(self):
        path, calls, fetch = self.failure_pair(FetchError("authentication"))
        result = self.execute(path, fetch, workers=2)
        self.assertEqual(len(calls), 2)
        self.assertEqual(result["state"], "failed")
        self.assertEqual(result["scoring"]["fetched_jobs"], 1)
        self.assertEqual(self.cache.status()["responses"], 1)
        self.assertEqual(self.execute(path, lambda _: (b"resumed", {}), workers=3)["state"], "completed")

    def test_interrupt_drains_inflight_success_and_can_resume(self):
        path, calls, fetch = self.failure_pair(KeyboardInterrupt())
        with self.assertRaises(KeyboardInterrupt):
            self.execute(path, fetch, workers=2)
        result = run_status(path)
        self.assertEqual(result["state"], "paused")
        self.assertEqual(result["scoring"]["fetched_jobs"], 1)
        self.assertEqual(len(calls), 2)
        self.assertEqual(self.execute(path, lambda _: (b"resumed", {}), workers=2)["state"], "completed")

    def test_two_parallel_runs_fetch_each_request_once(self):
        a, b = self.prepared("a"), self.prepared("b")
        calls, lock = [], Lock()

        def fetch(request):
            with lock:
                calls.append(request["variant"]["pos"])
            time.sleep(0.01)
            return b"synthetic", {}

        with ThreadPoolExecutor(max_workers=2) as pool:
            results = list(pool.map(lambda path: self.execute(path, fetch, workers=3), [a, b]))
        self.assertEqual(len(calls), 8)
        self.assertEqual(len(set(calls)), 8)
        self.assertEqual(sum(r["scoring"]["cache_hits"] for r in results), 8)
        self.assertTrue(all(r["state"] == "completed" for r in results))

    def test_stop_cancels_wait_for_another_process_cache_lock(self):
        path = self.prepared("run")
        with closing(sqlite3.connect(path / "run.private.sqlite")) as db:
            request = json.loads(db.execute(
                "SELECT request_json FROM jobs ORDER BY id LIMIT 1").fetchone()[0])
        key, _ = identity(request, RESPONSE_FORMAT)
        def fail(request):
            raise FetchError("rate_limit")
        with FileLock(str(self.cache.directory / "locks" / (key + ".lock"))):
            with ThreadPoolExecutor(max_workers=1) as pool:
                result = pool.submit(self.execute, path, fail, workers=2).result(timeout=5)
        self.assertEqual(result["state"], "failed")
        self.assertEqual(result["scoring"]["fetch_attempts"], 1)

    def test_sdk_client_initializes_once_across_workers(self):
        provider = Mock(return_value="synthetic-credential")
        client = ScoreClient(key_provider=provider)
        with patch("grpc.secure_channel") as create, patch("grpc.channel_ready_future"), \
                patch("grpc.intercept_channel"):
            with ThreadPoolExecutor(max_workers=4) as pool:
                list(pool.map(lambda _: client._initialize(), range(4)))
            self.assertEqual(provider.call_count, 1)
            self.assertEqual(create.call_count, 1)
        client.close()

    def test_invalid_worker_counts_fail_before_fetch(self):
        path = self.prepared("run")
        for workers in (0, 129, -1, True):
            with self.assertRaises(UserError):
                self.execute(path, lambda _: self.fail("No fetch"), workers=workers)

    def test_program_limit_128_is_accepted(self):
        path = self.prepared("run")
        result = self.execute(path, lambda _: (b"synthetic", {}), workers=128, max_jobs=1)
        self.assertEqual(result["runtime"]["workers"], 128)
        self.assertEqual(result["scoring"]["fetched_jobs"], 1)
        self.assertEqual(result["state"], "paused")

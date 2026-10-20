"""Offline scoring integration: synthetic sequences, clients and matrices."""

from contextlib import closing
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest
from unittest.mock import patch

import anndata
import numpy as np
import pandas as pd

from alphagenome_variant_scorer.api import FetchError, ScoreClient
from alphagenome_variant_scorer.cache import SharedCache
from alphagenome_variant_scorer.errors import UserError
from alphagenome_variant_scorer.prepare import prepare
from alphagenome_variant_scorer.reference import ensure_index
from alphagenome_variant_scorer.runner import run
from alphagenome_variant_scorer.scores import encode_scores, decode_scores
from alphagenome_variant_scorer.store import run_status


def scores():
    return anndata.AnnData(
        X=np.array([[np.nan, 0.5], [np.inf, -np.inf]], dtype=np.float32),
        obs=pd.DataFrame({"gene_id": ["invented1", "invented2"]}, index=["row1", "row2"]),
        var=pd.DataFrame({"name": ["track A", "track B"]}, index=["t1", "t2"]),
        layers={"quantiles": np.array([[np.nan, 0.8], [0.0, 1.0]], dtype=np.float32)})


class RunnerTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.fasta = self.root / "synthetic.fa"
        self.fasta.write_text(">chr1\nACGTACGTACGTACGT\n", encoding="ascii")
        self.index = self.root / "synthetic.fai"
        ensure_index(self.fasta, self.index, create=True)
        self.input = self.root / "input.tsv"
        self.input.write_text("#CHROM\tPOS\tREF\tALT\nchr1\t2\tC\tT\nchr1\t3\tG\tA\n",
                              encoding="ascii")
        self.cache = SharedCache(self.root / "cache")
        self.cache.initialize()
        self.calls = 0

    def prepared(self, name):
        path = self.root / name
        prepare(self.input, path, self.fasta, self.index, input_format="literal-tsv")
        return path

    def fetch(self, request):
        self.calls += 1
        return encode_scores(scores()), {"source": "synthetic"}

    def execute(self, path, fetch=None, **kwargs):
        return run(path, self.cache.directory, fetch or self.fetch,
                   min_interval=0, retry_delay=0, **kwargs)

    def test_payload_roundtrip_preserves_values_labels_layers(self):
        meta, arrays = decode_scores(encode_scores(scores()))
        np.testing.assert_array_equal(arrays["raw"], scores().X)
        self.assertEqual(arrays["raw"].dtype, np.float32)
        np.testing.assert_array_equal(arrays["layer_0"], scores().layers["quantiles"])
        self.assertEqual(meta["obs"]["index"], ["row1", "row2"])
        self.assertEqual(meta["tracks"]["data"], [["track A"], ["track B"]])

    def test_two_runs_share_responses_without_key_or_reference(self):
        a, b = self.prepared("a"), self.prepared("b")
        self.assertEqual(self.execute(a)["state"], "completed")
        self.fasta.unlink()
        client = ScoreClient()
        self.assertEqual(self.execute(b, client)["scoring"]["cache_hits"], 2)
        self.assertEqual(self.calls, 2)
        self.assertEqual(self.execute(a, client)["state"], "completed")
        self.assertEqual(self.calls, 2)

    def test_bounded_run_and_resume(self):
        path = self.prepared("run")
        self.assertEqual(self.execute(path, max_jobs=1)["state"], "paused")
        self.assertEqual(self.execute(path)["state"], "completed")
        self.assertEqual(self.calls, 2)

    def test_absent_result_does_not_stop_other_jobs_or_retry_on_resume(self):
        for code in ("not_found", "not_returned"):
            with self.subTest(code=code):
                path = self.prepared(code)
                calls = []

                def fetch(request):
                    calls.append(1)
                    if len(calls) == 1:
                        raise FetchError(code)
                    return self.fetch(request)

                # Distinct cached successes from previous subtest are allowed;
                # force the first callback's request to a new cache namespace.
                cache = SharedCache(self.root / (code + "-cache"))
                cache.initialize()
                result = run(path, cache.directory, fetch, min_interval=0, retry_delay=0)
                self.assertEqual(result["state"], "completed")
                self.assertEqual(result["scoring"]["errors"], {code: 1})
                self.assertEqual(len(calls), 2)
                run(path, cache.directory, lambda _: self.fail("No retry on resume"))
                self.assertEqual(cache.status()["responses"], 1)

    def test_transient_retries_and_failure_sanitization(self):
        path = self.prepared("run")
        calls = []
        def retry(request):
            calls.append(1)
            if len(calls) < 3:
                raise FetchError("unavailable")
            return self.fetch(request)
        result = self.execute(path, retry, max_jobs=1)
        self.assertEqual(result["scoring"]["fetch_attempts"], 3)
        def fail(request):
            raise FetchError("authentication")
        result = self.execute(path, fail)
        self.assertEqual(result["state"], "failed")
        self.assertEqual(result["scoring"]["errors"], {"authentication": 1})
        self.assertEqual(result["scoring"]["fetch_attempts"], 4)
        self.assertEqual(self.execute(path)["state"], "completed")

    def test_post_cache_commit_interruption_resumes_without_refetch(self):
        path = self.prepared("run")
        original = SharedCache.get_or_fetch
        def interrupted(cache, *args, **kwargs):
            original(cache, *args, **kwargs)
            raise KeyboardInterrupt()
        with patch.object(SharedCache, "get_or_fetch", interrupted):
            with self.assertRaises(KeyboardInterrupt):
                self.execute(path)
        self.assertEqual(run_status(path)["state"], "paused")
        result = self.execute(path)
        self.assertEqual(result["state"], "completed")
        self.assertEqual(result["scoring"]["cache_hits"], 1)
        self.assertEqual(self.calls, 2)

    def test_modified_requests_and_switched_cache_rejected(self):
        path = self.prepared("run")
        self.execute(path, max_jobs=1)
        other = SharedCache(self.root / "other")
        other.initialize()
        with self.assertRaises(UserError):
            run(path, other.directory, self.fetch)
        with closing(sqlite3.connect(path / "run.private.sqlite")) as db, db:
            db.execute("UPDATE jobs SET request_json=replace(request_json,'\"C\"','\"A\"')")
        with self.assertRaises(UserError):
            self.execute(path)

    def test_run_lock_prevents_concurrent_runner(self):
        from filelock import FileLock, Timeout
        path = self.prepared("run")
        with FileLock(str(path / "run.lock")):
            with self.assertRaises(Timeout):
                self.execute(path)
        self.assertEqual(self.calls, 0)

    def test_interrupted_schema_upgrade_rolls_back_and_can_resume(self):
        path = self.prepared("run")
        connect = sqlite3.connect

        class InterruptedUpgrade(sqlite3.Connection):
            def execute(self, statement, *args, **kwargs):
                if statement == "PRAGMA user_version=2":
                    raise KeyboardInterrupt()
                return super().execute(statement, *args, **kwargs)

        def interrupted_connect(*args, **kwargs):
            kwargs["factory"] = InterruptedUpgrade
            return connect(*args, **kwargs)

        with patch("alphagenome_variant_scorer.runner.sqlite3.connect", interrupted_connect):
            with self.assertRaises(KeyboardInterrupt):
                self.execute(path)
        with closing(connect(path / "run.private.sqlite")) as db:
            self.assertEqual(db.execute("PRAGMA user_version").fetchone()[0], 1)
            self.assertIsNone(db.execute(
                "SELECT name FROM sqlite_master WHERE name='job_results'").fetchone())
        self.assertEqual(self.calls, 0)
        self.assertEqual(self.execute(path)["state"], "completed")

    def test_atlas_adapter_encodes_response_and_sanitizes_exception(self):
        from unittest.mock import Mock
        client = ScoreClient()
        client.channel = object()
        client.atlas = Mock()
        client.atlas.query_variant.return_value = {"RNA_SEQ": scores()}
        request = {"api": "Atlas", "assembly": "GRCh38", "organism": "HOMO_SAPIENS",
                   "requested_scorer": "RNA_SEQ",
                   "variant": {"chrom": "chr1", "pos": 2, "ref": "C", "alt": "T"}}
        payload, _ = client(request)
        self.assertEqual(decode_scores(payload)[1]["raw"].shape, (2, 2))
        client.atlas.query_variant.side_effect = RuntimeError("synthetic-secret-not-to-log")
        with self.assertRaises(FetchError) as raised:
            client(request)
        self.assertEqual(str(raised.exception), "api_error")

    def test_prediction_adapter_uses_saved_conditions_without_sdk_retry_wrapper(self):
        from alphagenome.models import dna_client
        from unittest.mock import Mock
        client = ScoreClient()
        client.channel = object()
        client.model = Mock()
        request = {
            "api": "Prediction", "assembly": "GRCh38", "organism": "HOMO_SAPIENS",
            "requested_scorer": "RNA_SEQ",
            "variant": {"chrom": "chr1", "pos": 200, "ref": "T", "alt": "A"},
            "interval": {"start0": 0, "end0": 2**20},
            "model_version": "ALL_FOLDS", "merge_stranded_gene_tracks": True,
            "scorer": {"type": "GeneMaskLFCScorer", "requested_output": "RNA_SEQ"},
        }
        with patch.object(dna_client.DnaClient.score_variant, "__wrapped__",
                          return_value=[scores()]) as rpc:
            client(request)
            self.assertEqual(rpc.call_count, 1)
            self.assertEqual(rpc.call_args.args[1].width, 2**20)
            self.assertEqual(rpc.call_args.args[2].alternate_bases, "A")
            self.assertIs(rpc.call_args.kwargs["merge_stranded_gene_tracks"], True)

    def test_deadline_applies_to_unary_and_streaming_rpc(self):
        from alphagenome_variant_scorer.api import Deadline, _Details
        from unittest.mock import Mock
        interceptor = Deadline(12)
        details = _Details("method", None, [], None, None, None)
        continuation = Mock()
        interceptor.intercept_unary_unary(continuation, details, "request")
        self.assertEqual(continuation.call_args.args[0].timeout, 12)
        details = details._replace(timeout=3)
        interceptor.intercept_stream_stream(continuation, details, iter([]))
        self.assertEqual(continuation.call_args.args[0].timeout, 3)

    def test_quota_failure_stops_without_retry_or_cache(self):
        path = self.prepared("run")
        def fail(request):
            raise FetchError("rate_limit")
        result = self.execute(path, fail)
        self.assertEqual(result["scoring"]["fetch_attempts"], 1)
        self.assertEqual(result["state"], "failed")
        self.assertEqual(self.cache.status()["responses"], 0)

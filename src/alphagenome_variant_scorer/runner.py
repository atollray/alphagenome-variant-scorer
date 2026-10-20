"""Bounded parallel, cache-first scoring; only the coordinator writes the run DB."""

from concurrent.futures import CancelledError, ThreadPoolExecutor
from contextlib import closing
import hashlib
import json
import math
import queue
from threading import Event, Lock
from pathlib import Path
import sqlite3
import time

from filelock import FileLock

from .api import FetchError, validate_request
from .cache import SharedCache
from .errors import UserError
from .prepare import canonical, save_manifest
from .scores import RESPONSE_FORMAT
from .store import run_status


def run(run_dir: Path, cache_dir: Path, fetch, *, max_jobs=None, max_attempts=3,
        min_interval=1.0, retry_delay=5.0, workers=1):
    if (type(workers) is not int or not 1 <= workers <= 128 or max_attempts < 1 or (max_jobs is not None and max_jobs < 1)
            or not math.isfinite(min_interval) or min_interval < 0
            or not math.isfinite(retry_delay) or retry_delay < 0):
        raise UserError("Invalid run limits.")
    path = run_dir / "run.private.sqlite"
    if not path.is_file():
        raise UserError("Prepare this run before scoring.")
    cache = SharedCache(cache_dir)
    cache.status()
    with FileLock(str(run_dir / "run.lock"), timeout=0), closing(
        sqlite3.connect(path.resolve().as_uri() + "?mode=rw", uri=True, timeout=30)
    ) as db:
        version = db.execute("PRAGMA user_version").fetchone()[0]
        if version not in {1, 2}:
            raise UserError("Unsupported run database schema.")
        manifest = json.loads(db.execute("SELECT value FROM meta WHERE key='manifest'").fetchone()[0])
        if manifest["state"] not in {"prepared", "running", "paused", "failed", "completed"}:
            raise UserError("Preparation did not complete successfully.")
        previous_cache = manifest.get("cache_directory")
        if previous_cache is not None and previous_cache != str(cache.directory):
            raise UserError("Resume requires the same shared cache directory as the original run.")
        digest = hashlib.sha256()
        for job_id, request_json in db.execute("SELECT id,request_json FROM jobs ORDER BY id"):
            validate_request(json.loads(request_json))
            digest.update((job_id + "\n" + request_json + "\n").encode())
        if digest.hexdigest() != manifest["prepared_requests_sha256"]:
            raise UserError("Prepared requests changed; refusing to run.")
        db.execute("PRAGMA journal_mode=WAL")
        db.execute("PRAGMA synchronous=FULL")
        db.execute("PRAGMA foreign_keys=ON")
        with db:
            if version == 1:
                # Python's sqlite3 does not begin a transaction for CREATE TABLE.
                # Keep the new table and version change atomic across interruptions.
                db.execute("BEGIN IMMEDIATE")
                db.execute("""CREATE TABLE job_results(
                    job_id TEXT PRIMARY KEY REFERENCES jobs(id),
                    response_id TEXT, cache_hit INTEGER, attempts INTEGER NOT NULL DEFAULT 0,
                    error_code TEXT)""")
                db.execute("PRAGMA user_version=2")
            manifest["schema_version"] = 2
            manifest["cache_directory"] = str(cache.directory)
            manifest["runtime"] = {
                "workers": workers, "min_interval": min_interval,
                "max_attempts": max_attempts,
            }
            # Preparation's zero is not a live RPC counter. Count adapter attempts separately.
            manifest.pop("api_calls", None)

        def save(state):
            manifest["state"] = state
            manifest["fetch_attempts"] = db.execute(
                "SELECT COALESCE(sum(attempts),0) FROM job_results").fetchone()[0]
            with db:
                db.execute("UPDATE meta SET value=? WHERE key='manifest'", (canonical(manifest),))
            save_manifest(run_dir / "manifest.private.json", manifest)

        # Previously completed jobs retain their exact response, even after a refresh.
        cache.verify_references(response_id for (response_id,) in db.execute(
            "SELECT response_id FROM job_results JOIN jobs ON jobs.id=job_results.job_id "
            "WHERE jobs.state='completed'"
        ))
        save("running")
        _execute_jobs(db, cache, fetch, save, max_jobs=max_jobs,
                      max_attempts=max_attempts, min_interval=min_interval,
                      retry_delay=retry_delay, workers=workers)
    return run_status(run_dir)


def _execute_jobs(db, cache, fetch, save, *, max_jobs, max_attempts,
                  min_interval, retry_delay, workers):
    stop = Event()
    wake = Event()
    attempts = queue.Queue()
    gate = Lock()
    last_call = None
    dispatched = 0
    last_job_id = ""
    exhausted = False
    failed = False
    interrupted = None
    futures = {}
    pool = ThreadPoolExecutor(max_workers=workers, thread_name_prefix="agvs-fetch")

    def begin_attempt(job_id):
        nonlocal last_call
        # One start-time gate for all workers, including retries and both APIs.
        with gate:
            if last_call is not None:
                # Python 3.12 on Windows uses a coarse clock for monotonic().
                # Use the high-resolution monotonic counter and recheck after
                # waiting so timer rounding cannot shorten the interval.
                deadline = last_call + min_interval
                while (delay := deadline - time.perf_counter()) > 0:
                    if stop.wait(delay):
                        raise CancelledError()
            if stop.is_set():
                raise CancelledError()
            acknowledged = Event()
            attempts.put((job_id, acknowledged))
            wake.set()
            while not acknowledged.wait(0.05):
                if stop.is_set():
                    raise CancelledError()
            if stop.is_set():
                raise CancelledError()
            last_call = time.perf_counter()

    def work(job_id, request):
        def fetch_with_retries():
            for attempt in range(max_attempts):
                begin_attempt(job_id)
                try:
                    return fetch(request)
                except FetchError as exc:
                    if (exc.code not in {"timeout", "unavailable"}
                            or attempt + 1 == max_attempts):
                        raise
                    if stop.wait(min(60, retry_delay * 2**attempt)):
                        raise CancelledError()
        try:
            return cache.get_or_fetch(request, RESPONSE_FORMAT, fetch_with_retries,
                                      stop_event=stop)
        except FetchError as exc:
            if exc.code not in {"not_found", "not_returned"}:
                stop.set()
            raise
        except BaseException:
            stop.set()
            raise

    def acknowledge_attempts():
        while True:
            try:
                job_id, acknowledged = attempts.get_nowait()
            except queue.Empty:
                return
            if not stop.is_set():
                with db:
                    db.execute("UPDATE job_results SET attempts=attempts+1 WHERE job_id=?",
                               (job_id,))
            acknowledged.set()

    def receive(future, job_id):
        nonlocal failed, interrupted
        try:
            response, hit = future.result()
        except CancelledError:
            with db:
                db.execute("UPDATE jobs SET state='pending' WHERE id=?", (job_id,))
        except FetchError as exc:
            terminal = exc.code in {"not_found", "not_returned"}
            failed = failed or not terminal
            with db:
                db.execute("UPDATE jobs SET state=? WHERE id=?",
                           (exc.code if terminal else "failed", job_id))
                db.execute("UPDATE job_results SET error_code=? WHERE job_id=?",
                           (exc.code, job_id))
        except BaseException as exc:
            stop.set()
            interrupted = interrupted or exc
        else:
            # In-flight successes are kept even after another worker fails.
            # Cache commit always precedes the atomic reference + completion.
            with db:
                db.execute("UPDATE job_results SET response_id=?,cache_hit=?,error_code=NULL "
                           "WHERE job_id=?", (response.response_id, int(hit), job_id))
                db.execute("UPDATE jobs SET state='completed' WHERE id=?", (job_id,))

    def drain():
        acknowledge_attempts()
        for future in list(futures):
            if future.done():
                job_id = futures.pop(future)
                receive(future, job_id)

    try:
        while True:
            wake.clear()
            drain()
            while (not stop.is_set() and not exhausted and len(futures) < workers
                   and (max_jobs is None or dispatched < max_jobs)):
                row = db.execute(
                    "SELECT id,request_json FROM jobs WHERE id>? "
                    "AND state NOT IN ('completed','not_found','not_returned') "
                    "ORDER BY id LIMIT 1", (last_job_id,)).fetchone()
                if row is None:
                    exhausted = True
                    break
                job_id, request_json = row
                last_job_id = job_id
                with db:
                    db.execute("INSERT OR IGNORE INTO job_results(job_id) VALUES(?)", (job_id,))
                    db.execute("UPDATE job_results SET error_code=NULL WHERE job_id=?", (job_id,))
                    db.execute("UPDATE jobs SET state='running' WHERE id=?", (job_id,))
                future = pool.submit(work, job_id, json.loads(request_json))
                futures[future] = job_id
                future.add_done_callback(lambda _: wake.set())
                dispatched += 1
            if not futures:
                break
            wake.wait(0.05)
    except BaseException as exc:
        interrupted = exc
        stop.set()
    finally:
        # RPCs already in flight finish under their deadlines; no new requests start.
        # Also collect results saved just before interruption for immediate progress.
        while futures:
            wake.clear()
            drain()
            if futures:
                wake.wait(0.05)
        pool.shutdown(wait=True)
    remaining = db.execute("SELECT count(*) FROM jobs WHERE state NOT IN "
                           "('completed','not_found','not_returned')").fetchone()[0]
    save("failed" if failed else "paused" if interrupted or remaining else "completed")
    if interrupted is not None:
        raise interrupted

"""Shared cache contracts, with synthetic responses and no API client."""

from concurrent.futures import ThreadPoolExecutor
from contextlib import closing, redirect_stdout
import io
from pathlib import Path
import sqlite3
import tempfile
import time
import unittest

from alphagenome_variant_scorer.cache import SharedCache
from alphagenome_variant_scorer.cli import main
from alphagenome_variant_scorer.config import load, save_setting
from alphagenome_variant_scorer.errors import UserError

REQUEST = {"api": "synthetic", "variant": {"position": 12, "ref": "A", "alt": "C"}}
FORMAT = "synthetic-json-v1"


class CacheTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.cache = SharedCache(self.root / "shared")
        self.cache.initialize()

    def test_cross_run_reuse_and_missing_scores(self):
        payload = b'{"scores":[null,1.5],"states":["P01","evaluated"]}'
        first, hit = self.cache.get_or_fetch(REQUEST, FORMAT, lambda: (payload, {"source": "test"}))
        self.assertFalse(hit)
        second, hit = SharedCache(self.cache.directory).get_or_fetch(
            dict(reversed(list(REQUEST.items()))), FORMAT,
            lambda: self.fail("A second run must reuse the saved response"))
        self.assertTrue(hit)
        self.assertEqual(first, second)
        self.assertEqual(second.payload, payload)

    def test_concurrent_misses_fetch_once(self):
        calls = []

        def fetch():
            calls.append(1)
            time.sleep(0.05)
            return b"success", {}

        def run(_):
            return SharedCache(self.cache.directory).get_or_fetch(REQUEST, FORMAT, fetch)

        with ThreadPoolExecutor(max_workers=4) as pool:
            results = list(pool.map(run, range(4)))
        self.assertEqual(len(calls), 1)
        self.assertEqual(sum(not hit for _, hit in results), 1)
        self.assertEqual(len({response.response_id for response, _ in results}), 1)

    def test_failure_is_not_cached_and_retry_works(self):
        def fail():
            raise RuntimeError("synthetic network failure")
        with self.assertRaises(RuntimeError):
            self.cache.get_or_fetch(REQUEST, FORMAT, fail)
        self.assertIsNone(self.cache.lookup(REQUEST, FORMAT))
        self.cache.get_or_fetch(REQUEST, FORMAT, lambda: (b"ok", {}))
        self.assertEqual(self.cache.status()["responses"], 1)

    def test_request_and_format_changes_separate_entries(self):
        for request, fmt in [(REQUEST, FORMAT), ({**REQUEST, "model": "different"}, FORMAT),
                             (REQUEST, "synthetic-json-v2")]:
            _, hit = self.cache.get_or_fetch(request, fmt, lambda: (b"ok", {}))
            self.assertFalse(hit)
        self.assertEqual(self.cache.status()["requests"], 3)

    def test_refresh_preserves_original_response(self):
        first, _ = self.cache.get_or_fetch(REQUEST, FORMAT, lambda: (b"old", {}))
        newer, hit = self.cache.get_or_fetch(REQUEST, FORMAT, lambda: (b"new", {}), refresh=True)
        self.assertFalse(hit)
        self.assertNotEqual(first.response_id, newer.response_id)
        self.assertEqual(self.cache.read_response(first.response_id), first)
        self.assertEqual(self.cache.lookup(REQUEST, FORMAT), newer)
        self.assertEqual(self.cache.status()["responses"], 2)

    def test_post_save_interruption_reuses_committed_response(self):
        self.cache.get_or_fetch(REQUEST, FORMAT, lambda: (b"saved", {}))
        # Simulate an interrupted caller that never recorded the returned ID.
        response, hit = SharedCache(self.cache.directory).get_or_fetch(
            REQUEST, FORMAT, lambda: self.fail("Committed data must survive caller restart"))
        self.assertTrue(hit)
        self.assertEqual(response.payload, b"saved")

    def test_foreign_database_and_newer_schema_are_not_modified(self):
        foreign = self.root / "foreign"
        foreign.mkdir()
        path = foreign / "scores.private.sqlite"
        with closing(sqlite3.connect(path)) as db:
            db.execute("PRAGMA user_version=99")
        before = path.read_bytes()
        with self.assertRaises(UserError):
            SharedCache(foreign).initialize()
        self.assertEqual(before, path.read_bytes())
        with closing(sqlite3.connect(self.cache.path)) as db:
            db.execute("PRAGMA user_version=99")
        with self.assertRaises(UserError):
            self.cache.lookup(REQUEST, FORMAT)

    def test_corrupt_payload_is_rejected(self):
        self.cache.get_or_fetch(REQUEST, FORMAT, lambda: (b"ok", {}))
        with closing(sqlite3.connect(self.cache.path)) as db, db:
            db.execute("UPDATE responses SET payload=?", (b"bad",))
        with self.assertRaises(UserError):
            self.cache.lookup(REQUEST, FORMAT)

    def test_config_and_cli_do_not_require_reference(self):
        config = self.root / "settings.toml"
        destination = self.root / "new-cache"
        save_setting(config, "cache.directory", str(destination))
        self.assertFalse(destination.exists())
        self.assertEqual(load(config)["cache"]["directory"], str(destination.resolve()))
        with redirect_stdout(io.StringIO()):
            self.assertEqual(main(["--config", str(config), "cache", "init"]), 0)
            self.assertEqual(main(["--config", str(config), "cache", "status"]), 0)
            self.assertEqual(main(["--config", str(config), "config", "show"]), 0)
        self.assertEqual(SharedCache(destination).status()["responses"], 0)

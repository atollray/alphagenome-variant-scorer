"""Synthetic credentials only; ensure precedence, redaction and lazy loading."""

from contextlib import redirect_stdout, redirect_stderr
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from alphagenome_variant_scorer.api import FetchError, ScoreClient
from alphagenome_variant_scorer.cli import main
from alphagenome_variant_scorer.config import load, save_setting
from alphagenome_variant_scorer.auth import resolve_key
from alphagenome_variant_scorer.errors import UserError


class CredentialTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.config = self.root / "config.toml"
        self.env = self.root / "credentials.env"

    def write_config(self, auth):
        self.config.write_text("schema_version=1\n[auth]\n" + auth, encoding="utf-8")

    def test_config_key_wins_over_stale_environment_and_show_redacts(self):
        self.write_config('api_key="synthetic-new-secret"\n')
        with patch.dict("os.environ", {"ALPHAGENOME_API_KEY": "synthetic-old-secret"}):
            self.assertEqual(resolve_key(config_file=self.config), "synthetic-new-secret")
            out, err = io.StringIO(), io.StringIO()
            with redirect_stdout(out), redirect_stderr(err):
                self.assertEqual(main(["--config", str(self.config), "config", "show"]), 0)
            self.assertNotIn("synthetic-new-secret", out.getvalue() + err.getvalue())
            self.assertEqual(json.loads(out.getvalue())["auth"]["api_key"], "[redacted]")

    def test_missing_config_never_uses_environment(self):
        with patch.dict("os.environ", {"ALPHAGENOME_API_KEY": "synthetic-old-secret"}):
            with self.assertRaises(UserError):
                resolve_key(config_file=self.config)

    def test_other_settings_preserve_key(self):
        self.write_config('api_key="synthetic-secret"\n')
        save_setting(self.config, "cache.directory", str(self.root / "cache"))
        self.assertEqual(load(self.config)["auth"]["api_key"], "synthetic-secret")

    def test_env_file_setting_rejected(self):
        self.write_config('env_file="credentials.env"\n')
        with self.assertRaises(UserError):
            load(self.config)

    def test_invalid_config_does_not_print_secret(self):
        self.write_config('api_key="synthetic-secret\n')
        out, err = io.StringIO(), io.StringIO()
        with redirect_stdout(out), redirect_stderr(err):
            self.assertEqual(main(["--config", str(self.config), "config", "show"]), 2)
        self.assertNotIn("synthetic-secret", out.getvalue() + err.getvalue())

    def test_client_resolves_key_lazily_and_reports_safe_error(self):
        calls = []
        def provider():
            calls.append(1)
            raise UserError("synthetic-secret")
        client = ScoreClient(key_provider=provider)
        self.assertEqual(calls, [])
        with self.assertRaises(FetchError) as raised:
            client._initialize()
        self.assertEqual(str(raised.exception), "credential_unavailable")
        self.assertEqual(calls, [1])

"""Regression checks for accidental publication of local artifacts."""

import importlib.util
from pathlib import Path
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("publication_check", ROOT / "tools" / "check_public_files.py")
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)


class PublicFilesTest(unittest.TestCase):
    def test_only_explicit_synthetic_data_are_allowed(self):
        for name in [
            "data/sample.tsv", "tests/fixtures/sample.fa", "config.toml",
            "docs/implementation-plan.md", "docs/publication.md",
            "results.sqlite", ".env", "src/sample.private.py", "tools/secret.json",
        ]:
            with self.subTest(name=name):
                self.assertFalse(checker.allowed(name))
        self.assertTrue(checker.allowed("tests/fixtures/synthetic.fa"))
        self.assertTrue(checker.allowed("src/alphagenome_variant_scorer/cli.py"))

    def test_secret_patterns_and_personal_paths(self):
        fake = ("ghp_" + "x" * 32).encode()
        self.assertIn("GitHub credential", checker.inspect_content("example", fake))
        fake_path = ("C:" + "/Users/" + "example/" + "private.fa").encode()
        self.assertIn("personal Windows path", checker.inspect_content("example", fake_path))

    @unittest.skipUnless((ROOT / ".git").exists(), "Git ignore checks require a checkout")
    def test_ignore_rules_cover_private_artifacts(self):
        paths = [
            ".env", "config.toml", "settings.local.toml", "docs/publication.md",
            "runs/example/batch.sqlite-wal", "cache/old.json",
            "input.private.json", "hg38.fa", "hg38.fa.fai", "variants.tsv",
            "scores.npz", "reference/chrSynthetic.gz",
        ]
        result = subprocess.run(
            ["git", "check-ignore", "-z", "--stdin"], cwd=ROOT,
            input=("\0".join(paths) + "\0").encode(), capture_output=True, check=False,
        )
        self.assertEqual(set(result.stdout.decode().split("\0")) - {""}, set(paths))
        public = subprocess.run(
            ["git", "check-ignore", "--no-index", "tests/fixtures/synthetic.fa",
             "examples/config.example.toml"], cwd=ROOT, capture_output=True,
        )
        self.assertEqual(public.returncode, 1)


if __name__ == "__main__":
    unittest.main()

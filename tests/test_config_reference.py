"""Offline behavior tests using only invented sequences and temporary paths."""

from contextlib import redirect_stdout, redirect_stderr
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from filelock import FileLock

from alphagenome_variant_scorer.cli import main
from alphagenome_variant_scorer.config import (
    default_config_path, load, resolved_reference, save_setting, select_reference,
)
from alphagenome_variant_scorer.errors import UserError
from alphagenome_variant_scorer.reference import Reference, check_reference, ensure_index

FIXTURE = Path(__file__).parent / "fixtures" / "synthetic.fa"


class WorkspaceTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.fasta = self.root / "reference space 日本語.fa"
        shutil.copyfile(FIXTURE, self.fasta)
        self.config = self.root / "settings" / "config.toml"
        self.index = self.root / "indices" / "synthetic.fai"

    def cli(self, *args):
        out, err = io.StringIO(), io.StringIO()
        with redirect_stdout(out), redirect_stderr(err):
            code = main(["--config", str(self.config), *args])
        return code, out.getvalue(), err.getvalue()

    def test_roundtrip_does_not_index_and_preserves_other_settings(self):
        save_setting(self.config, "reference.fasta", str(self.fasta))
        save_setting(self.config, "reference.index", str(self.index))
        save_setting(self.config, "reference.assembly", "GRCh38")
        data = load(self.config)
        self.assertEqual(data["reference"]["fasta"], str(self.fasta.resolve()))
        self.assertEqual(data["reference"]["index"], str(self.index.resolve()))
        self.assertFalse(self.index.exists())
        self.assertFalse(Path(str(self.fasta) + ".fai").exists())

    def test_changing_fasta_clears_old_index_without_deleting_it(self):
        save_setting(self.config, "reference.fasta", str(self.fasta))
        save_setting(self.config, "reference.index", str(self.index))
        ensure_index(self.fasta, self.index, create=True)
        old = self.index.read_bytes()
        other = self.root / "other.fa"
        shutil.copyfile(self.fasta, other)
        save_setting(self.config, "reference.fasta", str(other))
        self.assertNotIn("index", load(self.config)["reference"])
        self.assertEqual(self.index.read_bytes(), old)

    def test_relative_config_paths_use_config_directory(self):
        self.config.parent.mkdir()
        self.config.write_text(
            'schema_version=1\n[reference]\nfasta="../reference space 日本語.fa"\n',
            encoding="utf-8",
        )
        ref = resolved_reference(load(self.config), self.config)
        self.assertEqual(ref["fasta"], str(self.fasta.resolve()))

    def test_cli_override_does_not_inherit_index_or_change_config(self):
        save_setting(self.config, "reference.fasta", str(self.fasta))
        save_setting(self.config, "reference.index", str(self.index))
        original = self.config.read_bytes()
        other = self.root / "other.fa"
        shutil.copyfile(self.fasta, other)
        chosen, index, sources = select_reference(
            self.config, explicit_config=True, fasta=str(other), index=None,
        )
        self.assertEqual(chosen, other)
        self.assertEqual(index, Path(str(other) + ".fai"))
        self.assertEqual(sources, {"fasta": "argument", "index": "adjacent"})
        self.assertEqual(original, self.config.read_bytes())

    def test_missing_explicit_config_is_not_silently_ignored(self):
        with self.assertRaises(UserError):
            select_reference(self.config, explicit_config=True, fasta=str(self.fasta), index=None)

    def test_invalid_config_is_not_overwritten_or_echoed(self):
        self.config.parent.mkdir()
        for content in [
            'schema_version=2\n',
            'schema_version=true\n',
            'schema_version=1\n[reference]\nfasta=42\n',
            'schema_version=1\n[reference]\nassembly="GRCh37"\n',
            'schema_version=1\n[credentials]\npassword="do-not-echo-this"\n',
            'schema_version=1\nmalformed="do-not-echo-this\n',
        ]:
            with self.subTest(content=content):
                self.config.write_text(content, encoding="utf-8")
                original = self.config.read_bytes()
                code, out, err = self.cli("config", "set", "reference.fasta", str(self.fasta))
                self.assertEqual(code, 2)
                self.assertEqual(out, "")
                self.assertNotIn("do-not-echo-this", err)
                self.assertEqual(self.config.read_bytes(), original)

    def test_config_lock_blocks_a_second_writer(self):
        save_setting(self.config, "reference.fasta", str(self.fasta))
        original = self.config.read_bytes()
        with FileLock(str(self.config) + ".lock"):
            with self.assertRaises(UserError):
                save_setting(self.config, "reference.assembly", "GRCh38")
        self.assertEqual(self.config.read_bytes(), original)

    def test_failed_atomic_replace_preserves_config(self):
        save_setting(self.config, "reference.fasta", str(self.fasta))
        original = self.config.read_bytes()
        with patch("alphagenome_variant_scorer.config.os.replace", side_effect=OSError):
            with self.assertRaises(OSError):
                save_setting(self.config, "reference.index", str(self.index))
        self.assertEqual(self.config.read_bytes(), original)
        self.assertEqual(list(self.config.parent.glob("*.tmp")), [])

    def test_default_locations(self):
        expected = self.root / "alphagenome-variant-scorer" / "config.toml"
        self.assertEqual(default_config_path(
            platform="nt", environ={"LOCALAPPDATA": str(self.root)}), expected)
        self.assertEqual(default_config_path(
            platform="posix", environ={"XDG_CONFIG_HOME": str(self.root)}), expected)
        self.assertEqual(default_config_path(
            platform="posix", environ={}, home=self.root),
            self.root / ".config" / "alphagenome-variant-scorer" / "config.toml")
        with self.assertRaises(UserError):
            default_config_path(platform="posix", environ={"XDG_CONFIG_HOME": "relative"})
        with self.assertRaises(UserError):
            default_config_path(platform="nt", environ={})

    def test_reference_check_is_read_only_by_default(self):
        before = self.fasta.read_bytes()
        with self.assertRaises(UserError):
            ensure_index(self.fasta, self.index)
        self.assertFalse(self.index.exists())
        ensure_index(self.fasta, self.index, create=True)
        initial_index = self.index.read_bytes()
        initial_mtime = self.index.stat().st_mtime_ns
        report = check_reference(self.fasta, self.index, create=True, hash_file=True)
        self.assertEqual(report["sha256"], hashlib.sha256(before).hexdigest())
        self.assertEqual(report["contigs"], {"chrSynthetic1": 20, "chrSynthetic2": 8})
        self.assertFalse(report["assembly_verified"])
        self.assertFalse(report["api_reference_match_verified"])
        self.assertEqual(self.fasta.read_bytes(), before)
        self.assertEqual(self.index.read_bytes(), initial_index)
        self.assertEqual(self.index.stat().st_mtime_ns, initial_mtime)

    def test_extraction_coordinates_softmask_and_unknown_contigs(self):
        ensure_index(self.fasta, self.index, create=True)
        with Reference(self.fasta, self.index) as ref:
            self.assertEqual(ref.sequence("chrSynthetic1", 0, 5), "ACGTA")
            self.assertEqual(ref.sequence("chrSynthetic1", 8, 12), "NNAA")
            self.assertEqual(ref.sequence("chrSynthetic1", 20, 20), "")
            self.assertEqual(ref.extract(SimpleNamespace(
                chromosome="chrSynthetic2", start=1, end=4)), "TGC")
            for start, end in [(-1, 3), (0, 21), (5, 4), (True, 2)]:
                with self.assertRaises(UserError):
                    ref.sequence("chrSynthetic1", start, end)
            with self.assertRaises(UserError):
                ref.sequence("Synthetic1", 0, 1)

    def test_stale_index_is_not_replaced(self):
        ensure_index(self.fasta, self.index, create=True)
        before = self.index.read_bytes()
        timestamp = self.index.stat().st_mtime_ns + 2_000_000_000
        os.utime(self.fasta, ns=(timestamp, timestamp))
        with self.assertRaises(UserError):
            ensure_index(self.fasta, self.index, create=True)
        self.assertEqual(self.index.read_bytes(), before)

    def test_invalid_and_duplicate_fasta_do_not_publish_index(self):
        for payload in [b"", b"not fasta\n", b">x\nACGT\n>x\nTTTT\n"]:
            with self.subTest(payload=payload):
                self.fasta.write_bytes(payload)
                with self.assertRaises(UserError):
                    ensure_index(self.fasta, self.index, create=True)
                self.assertFalse(self.index.exists())

    def test_compressed_fasta_and_fasta_as_index_are_rejected(self):
        zipped = self.root / "reference.fa.gz"
        zipped.write_bytes(b"\x1f\x8b")
        for fasta, index in [(zipped, self.index), (self.fasta, self.fasta)]:
            with self.assertRaises(UserError):
                ensure_index(fasta, index, create=True)

    def test_cli_end_to_end_without_network(self):
        # Fails immediately if a future implementation accidentally adds a network call.
        with patch("socket.socket.connect", side_effect=AssertionError("network forbidden")):
            code, _, err = self.cli("config", "set", "reference.fasta", str(self.fasta))
            self.assertEqual((code, err), (0, ""))
            code, out, err = self.cli("config", "show")
            self.assertEqual((code, err), (0, ""))
            self.assertEqual(json.loads(out)["reference"]["fasta"], str(self.fasta))
            code, out, err = self.cli(
                "reference", "check", "--index", str(self.index), "--create-index", "--sha256")
            self.assertEqual((code, err), (0, ""))
            self.assertEqual(json.loads(out)["contig_count"], 2)


if __name__ == "__main__":
    unittest.main()

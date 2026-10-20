"""All positions and sequences in this module are synthetic."""

from contextlib import closing, redirect_stdout, redirect_stderr
import csv
import io
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest
from unittest.mock import patch

from alphagenome_variant_scorer.cli import main
from alphagenome_variant_scorer.errors import UserError
from alphagenome_variant_scorer.prepare import prepare, file_hash, DB_NAME, MANIFEST_NAME
from alphagenome_variant_scorer.reference import ensure_index
from alphagenome_variant_scorer.store import run_status


class PrepareTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.shared = tempfile.TemporaryDirectory()
        cls.fasta = Path(cls.shared.name) / "synthetic.fa"
        cls.index = Path(cls.shared.name) / "synthetic.fai"
        bases = list("ACGT" * (2**18 + 128))
        bases[999:1008] = list("CAAAAGTGC")
        cls.sequence = "".join(bases)
        cls.fasta.write_text(
            ">chr1 invented sequence for tests\n" + cls.sequence +
            "\n>chr2 short synthetic contig\nACGTACGT\n>chr3 synthetic homopolymer\nAAAAA\n",
            encoding="ascii")
        ensure_index(cls.fasta, cls.index, create=True)

    @classmethod
    def tearDownClass(cls):
        cls.shared.cleanup()

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.input = self.root / "input.tsv"
        self.run = self.root / "run"

    def literal(self, rows):
        self.write(["#CHROM", "POS", "REF", "ALT", "note"], rows)

    def write(self, header, rows):
        with self.input.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(header)
            writer.writerows(rows)

    def prepare(self, **kwargs):
        return prepare(self.input, self.run, self.fasta, self.index,
                       input_format=kwargs.pop("input_format", "literal-tsv"), **kwargs)

    def query(self, sql):
        with closing(sqlite3.connect((self.run / DB_NAME).as_uri() + "?mode=ro", uri=True)) as db:
            return db.execute(sql).fetchall()

    def test_literal_rows_preserve_duplicates_and_normalize_without_network(self):
        self.literal([
            ["chr1", "1000", "c", "g", "first"],
            ["chr1", "1000", "C", "G", "second"],
            ["chr1", "1002", "AA", "A", "repeat deletion one"],
            ["chr1", "1003", "AA", "A", "repeat deletion two"],
            ["chr1", "1000", "CA", "TG", "MNV"],
        ])
        original = self.input.read_bytes()
        reference_hash = file_hash(self.fasta)
        with patch("socket.socket.connect", side_effect=AssertionError("network forbidden")):
            result = self.prepare()
        self.assertEqual(result["counts"], {"rows": {"ready": 5}, "jobs": {"Atlas": 1, "Prediction": 2}})
        self.assertEqual((self.run / "input.private.tsv").read_bytes(), original)
        self.assertEqual(self.input.read_bytes(), original)
        self.assertEqual(file_hash(self.fasta), reference_hash)
        rows = self.query("SELECT row_number,original_json,job_id,variant_kind,details_json FROM rows ORDER BY row_number")
        self.assertEqual(rows[0][2], rows[1][2])
        self.assertEqual(rows[2][2], rows[3][2])
        self.assertEqual(json.loads(rows[0][1])[2], "c")
        self.assertEqual(rows[4][3], "MNV")
        normalized = json.loads(rows[2][4])["normalized_variant"]
        self.assertEqual(normalized, {"chrom": "chr1", "pos": 1000, "ref": "CA", "alt": "C"})
        for (serialized,) in self.query("SELECT request_json FROM jobs"):
            request = json.loads(serialized)
            self.assertNotIn("note", serialized)
            if request["api"] == "Prediction":
                self.assertEqual(request["interval"]["end0"] - request["interval"]["start0"], 2**20)
                self.assertEqual(request["model_version"], "ALL_FOLDS")
        manifest = json.loads((self.run / MANIFEST_NAME).read_text(encoding="utf-8"))
        self.assertEqual(manifest["input_sha256"], file_hash(self.input))
        self.assertEqual(manifest["reference"]["sha256"], reference_hash)
        self.assertEqual(manifest["api_calls"], 0)

    def test_invalid_unsupported_and_blank_rows_remain_in_database(self):
        self.literal([
            ["chr1", "1000", "A", "G", "bad ref"],
            ["chr1", "0", "A", "G", "bad position"],
            ["chr1", "1000", "C", "C", "no change"],
            ["chr1", "1000", "C", "A,G", "multi"],
            ["chr1", "1000", "C", "<DEL>", "symbolic"],
            ["chr1", "1000", "C", "N", "ambiguous allele"],
            ["1", "1000", "C", "G", "no inferred alias"],
            [],
        ])
        result = self.prepare()
        self.assertEqual(result["counts"]["rows"], {"invalid": 6, "unsupported": 2})
        self.assertEqual(result["counts"]["jobs"], {})
        self.assertEqual(len(self.query("SELECT * FROM rows")), 8)
        self.assertIn("ref_mismatch", result["reasons"])

    def test_snv_prediction_is_explicit_and_short_contigs_are_classified(self):
        self.literal([
            ["chr1", "1000", "C", "G", ""],
            ["chr2", "1", "A", "C", ""],
        ])
        result = self.prepare(snv_source="prediction")
        self.assertEqual(result["counts"]["jobs"], {"Prediction": 1})
        self.assertEqual(result["reasons"], {"contig_shorter_than_prediction_window": 1})

    def test_sv_structures_and_haplotype_equivalence(self):
        self.write(
            ["SV_chrom", "SV_start", "SV_end", "SV_type", "SV_length", "note"],
            [["chr1", "1005", "1007", "DEL", "-3", "deletion"],
             ["chr1", "1005", "1007", "DUP", "3", "duplication"]],
        )
        result = self.prepare(input_format="sv-tsv", dup_model="tandem-extra-copy")
        self.assertEqual(result["counts"]["rows"], {"ready": 2})
        for kind, serialized in self.query("SELECT variant_kind,details_json FROM rows"):
            detail = json.loads(serialized)
            v = detail["normalized_variant"]
            observed = self.sequence[:v["pos"]-1] + v["alt"] + self.sequence[v["pos"]-1+len(v["ref"]):]
            if kind == "SV_DEL":
                expected = self.sequence[:1004] + self.sequence[1007:]
            else:
                expected = self.sequence[:1007] + self.sequence[1004:1007] + self.sequence[1007:]
            self.assertEqual(observed, expected)
            self.assertTrue(detail["haplotype_equality_verified"])
            self.assertIsNotNone(detail["structure_assumption"])

    def test_sv_missing_structure_length_and_anchor(self):
        self.write(
            ["SV_chrom", "SV_start", "SV_end", "SV_type", "SV_length"],
            [["chr1", "1005", "1007", "DUP", "3"],
             ["chr1", "1005", "1007", "DEL", "2"],
             ["chr1", "1", "3", "DEL", "3"],
             ["chr1", "1005", "1007", "INV", "3"]],
        )
        result = self.prepare(input_format="sv-tsv")
        self.assertEqual(result["reasons"], {
            "dup_structure_unspecified": 1, "sv_length_mismatch": 1,
            "left_anchor_unavailable": 1, "sv_type_not_supported": 1,
        })

    def test_near_chromosome_end_and_normalization_context(self):
        size = len(self.sequence)
        self.literal([
            ["chr1", str(size), self.sequence[-1], "A", "terminal SNV"],
            ["chr1", str(size), self.sequence[-1], self.sequence[-1] + "A", "terminal INS"],
        ])
        result = self.prepare(snv_source="prediction")
        self.assertEqual(result["counts"]["rows"], {"ready": 1, "unsupported": 1})
        request = json.loads(self.query("SELECT request_json FROM jobs")[0][0])
        self.assertEqual(request["interval"]["end0"], size)
        self.assertEqual(result["reasons"], {"normalization_context_out_of_bounds": 1})

    def test_insertion_and_complex_replacement_keep_literal_sequence(self):
        self.literal([
            ["chr1", "1000", "C", "CAA", "insertion"],
            ["chr1", "1000", "CA", "GTG", "replacement"],
        ])
        result = self.prepare()
        self.assertEqual(result["counts"]["rows"], {"ready": 2})
        self.assertEqual({r[0] for r in self.query("SELECT variant_kind FROM rows")},
                         {"INS", "COMPLEX_REPLACEMENT"})
        for (details,) in self.query("SELECT details_json FROM rows"):
            self.assertTrue(json.loads(details)["haplotype_equality_verified"])

    def test_left_edge_normalization_and_oversized_integer_are_classified(self):
        self.literal([
            ["chr3", "2", "AA", "A", "left shift reaches chromosome start"],
            ["chr1", "1" * 100, "A", "C", "invalid coordinate"],
        ])
        result = self.prepare()
        self.assertEqual(result["reasons"], {
            "normalization_reference_bounds_or_bases": 1, "invalid_integer": 1,
        })

    def test_existing_run_is_never_overwritten(self):
        self.literal([["chr1", "1000", "C", "G", ""]])
        self.prepare()
        original = (self.run / DB_NAME).read_bytes()
        with self.assertRaises(UserError):
            self.prepare()
        self.assertEqual((self.run / DB_NAME).read_bytes(), original)

    def test_bad_header_stops_before_creating_output(self):
        self.write(["#CHROM", "POS", "REF", "REF"], [["chr1", "1000", "C", "G"]])
        with self.assertRaises(UserError):
            self.prepare()
        self.assertFalse(self.run.exists())

    def test_failed_preparation_is_not_presented_as_ready(self):
        self.literal([["chr1", "1000", "C", "G", ""]])
        with patch("alphagenome_variant_scorer.prepare.prepare_variant", side_effect=UserError("synthetic failure")):
            with self.assertRaises(UserError):
                self.prepare()
        self.assertEqual(run_status(self.run)["state"], "preparation_failed")
        self.assertTrue((self.run / "input.private.tsv").exists())

    def test_status_ignores_broken_config_and_needs_no_reference(self):
        self.literal([["chr1", "1000", "C", "G", ""]])
        self.prepare()
        config = self.root / "broken.toml"
        config.write_text("invalid", encoding="utf-8")
        before = file_hash(self.run / DB_NAME)
        out = io.StringIO()
        with redirect_stdout(out), patch(
            "alphagenome_variant_scorer.config.config_path", side_effect=AssertionError("config forbidden")
        ), patch("socket.socket.connect", side_effect=AssertionError("network forbidden")):
            code = main(["--config", str(config), "status", "--run-dir", str(self.run)])
        self.assertEqual(code, 0)
        self.assertEqual(json.loads(out.getvalue())["state"], "prepared")
        self.assertEqual(file_hash(self.run / DB_NAME), before)

    def test_cli_prepare_with_explicit_reference(self):
        self.literal([["chr1", "1000", "C", "G", ""]])
        out, err = io.StringIO(), io.StringIO()
        missing_default = self.root / "no-default.toml"
        with patch("alphagenome_variant_scorer.cli.config_path", return_value=missing_default):
            with redirect_stdout(out), redirect_stderr(err):
                code = main([
                    "prepare", "--input", str(self.input), "--input-format", "literal-tsv",
                    "--run-dir", str(self.run), "--reference-fasta", str(self.fasta), "--index", str(self.index),
                ])
        self.assertEqual((code, err.getvalue()), (0, ""))
        self.assertEqual(json.loads(out.getvalue())["counts"]["rows"], {"ready": 1})


if __name__ == "__main__":
    unittest.main()

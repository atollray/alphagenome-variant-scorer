"""Synthetic aggregation and read-only export contracts, without API access."""

from contextlib import redirect_stdout
import csv
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import anndata
import numpy as np
import pandas as pd

from alphagenome_variant_scorer.aggregate import summarize
from alphagenome_variant_scorer.api import FetchError
from alphagenome_variant_scorer.cache import SharedCache
from alphagenome_variant_scorer.cli import main
from alphagenome_variant_scorer.errors import UserError
from alphagenome_variant_scorer.export import export
from alphagenome_variant_scorer.prepare import prepare
from alphagenome_variant_scorer.reference import ensure_index
from alphagenome_variant_scorer.runner import run
from alphagenome_variant_scorer.scores import encode_scores


def payload(values, *, ids=None):
    array = np.asarray(values, dtype=np.float64)
    n, m = array.shape
    obs = pd.DataFrame(
        {"gene_id": ids or [f"invented{i}" for i in range(n)],
         "gene_name": [f"GENE{i}" for i in range(n)]},
        index=[str(i) for i in range(n)])
    tracks = pd.DataFrame({"name": [f"track{i}" for i in range(m)]},
                          index=[str(i) for i in range(m)])
    return encode_scores(anndata.AnnData(X=array, obs=obs, var=tracks))


class AggregateTest(unittest.TestCase):
    def test_complete_both_signs_and_unrounded_boundary(self):
        result = summarize(payload([[-0.7, 0.8], [0.49999, -0.49999], [0.5, 0]]), 0.5)
        self.assertEqual(result["evaluation"], "complete")
        self.assertEqual(result["hit_genes"], 2)
        self.assertTrue(result["genes"][0]["positive_hit"])
        self.assertTrue(result["genes"][0]["negative_hit"])
        self.assertFalse(result["genes"][1]["positive_hit"])
        self.assertEqual(result["genes"][1]["max_raw"], 0.49999)

    def test_partial_no_hits_and_partial_hits_keep_p01(self):
        for value, expected in [(0.2, "no_hits"), (0.6, "hits")]:
            result = summarize(payload([[value, 0], [np.nan, np.nan]]), 0.5)
            self.assertEqual(result["code"], "P01")
            self.assertEqual(result["evaluation"], "partial")
            self.assertEqual(result["threshold_result"], expected)
            self.assertEqual(result["genes"][1]["missing_reason"], "all_nan")
            self.assertIsNone(result["genes"][1]["max_raw"])

    def test_partial_tracks_are_excluded_even_with_large_finite_value(self):
        result = summarize(payload([[0.1, 0.2], [100, np.nan], [np.inf, -np.inf]]), 0.5)
        self.assertEqual(result["valid_genes"], 1)
        self.assertEqual(result["hit_genes"], 0)
        self.assertEqual(result["genes"][1]["missing_reason"], "partial_nonfinite")
        self.assertEqual(result["genes"][2]["missing_reason"], "all_nonfinite")
        self.assertEqual(result["genes"][2]["positive_infinity_tracks"], 1)

    def test_all_invalid_empty_genes_and_empty_tracks_are_distinct(self):
        for values, state in [
            ([[np.nan, np.nan]], "unavailable"),
            (np.empty((0, 2)), "no_genes"), (np.empty((2, 0)), "no_tracks"),
        ]:
            result = summarize(payload(values), 0.5)
            self.assertEqual(result["evaluation"], state)
            self.assertEqual(result["code"], "E02")
            self.assertIsNone(result["hit_genes"])
            self.assertEqual(result["threshold_result"], "not_evaluated")

    def test_bad_gene_metadata_is_not_p01(self):
        with self.assertRaises(UserError):
            summarize(payload([[1], [2]], ids=["duplicate", "duplicate"]), 0.5)
        for threshold in [-1, float("nan"), float("inf")]:
            with self.assertRaises(UserError):
                summarize(payload([[1]]), threshold)

    def test_corrupt_payload_is_a_validation_error(self):
        with self.assertRaises(UserError):
            summarize(b"not-an-npz-file", 0.5)


class ExportTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        fasta = self.root / "synthetic.fa"
        fasta.write_text(">chr1\nACGTACGTACGTACGT\n", encoding="ascii")
        index = self.root / "synthetic.fai"
        ensure_index(fasta, index, create=True)
        source = self.root / "input.tsv"
        source.write_text(
            "#CHROM\tPOS\tREF\tALT\nchr1\t2\tC\tT\nchr1\t2\tC\tT\n"
            "chr1\t3\tG\tA\nchr1\tbad\tA\tC\nchr1\t4\tT\tG\textra\n",
            encoding="ascii")
        self.run_dir = self.root / "run"
        prepare(source, self.run_dir, fasta, index, input_format="literal-tsv")
        self.cache = SharedCache(self.root / "cache")
        self.cache.initialize()
        self.response = payload([[0.2, -0.3], [np.nan, np.nan], [np.inf, -np.inf]])

    def fetch(self, request):
        return self.response, {"source": "synthetic"}

    def finish(self, **kwargs):
        return run(self.run_dir, self.cache.directory, self.fetch,
                   min_interval=0, retry_delay=0, **kwargs)

    def test_json_preserves_rows_duplicates_and_no_config_required(self):
        self.finish()
        output = self.root / "results.json"
        broken = self.root / "broken.toml"
        broken.write_text("invalid", encoding="ascii")
        with patch("socket.socket.connect", side_effect=AssertionError("No network")), redirect_stdout(io.StringIO()):
            self.assertEqual(main(["--config", str(broken), "export", "--run-dir",
                                   str(self.run_dir), "--output", str(output), "--format", "json"]), 0)
        result = json.loads(output.read_text(encoding="utf-8"))
        self.assertEqual(len(result["rows"]), 5)
        self.assertEqual(result["rows"][0]["response_id"], result["rows"][1]["response_id"])
        self.assertEqual(result["rows"][0]["summary"]["code"], "P01")
        self.assertEqual(result["rows"][4]["original"][-1], "extra")
        self.assertEqual(result["rows"][3]["input_state"], "invalid")
        sidecar = json.loads(Path(str(output) + ".manifest.json").read_text())
        self.assertEqual(sidecar["output_sha256"], hashlib.sha256(output.read_bytes()).hexdigest())
        self.assertNotIn("NaN", output.read_text())
        self.assertNotIn("Infinity", output.read_text())

    def test_tsv_preserves_source_cells_and_rethresholds_without_fetch(self):
        self.finish()
        for threshold in [0.5, 0.1]:
            output = self.root / f"result-{threshold}.tsv"
            with patch("socket.socket.connect", side_effect=AssertionError("No network")):
                export(self.run_dir, output, threshold=threshold)
            with output.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 5)
            self.assertEqual(rows[0]["POS"], "2")
            self.assertEqual(rows[0]["AGVS_code"], "P01")
            self.assertEqual(rows[0]["AGVS_hit_genes"], "0" if threshold == 0.5 else "1")
            self.assertEqual(json.loads(rows[4]["AGVS_original_values_json"])[-1], "extra")

    def test_long_output_keeps_nonfinite_and_excluded_values(self):
        self.finish()
        output = self.root / "scores.tsv"
        export(self.run_dir, output, format="scores-tsv")
        with output.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 18)  # 3 ready input rows, 3 genes, 2 tracks.
        self.assertEqual(rows[2]["raw"], "")
        self.assertEqual(rows[2]["value_state"], "nan")
        self.assertEqual(rows[4]["value_state"], "positive_infinity")
        self.assertEqual(rows[5]["value_state"], "negative_infinity")
        self.assertEqual(rows[4]["gene_valid"], "False")

    def test_partial_requires_opt_in_and_supports_prepared_schema(self):
        output = self.root / "partial.json"
        with self.assertRaises(UserError):
            export(self.run_dir, output, format="json")
        export(self.run_dir, output, format="json", allow_partial=True)
        result = json.loads(output.read_text())
        self.assertTrue(result["manifest"]["partial_snapshot"])
        self.assertEqual(result["rows"][0]["job_state"], "pending")
        self.assertIsNone(result["rows"][0]["summary"])

    def test_existing_output_is_never_replaced(self):
        self.finish()
        output = self.root / "existing.tsv"
        output.write_text("original", encoding="ascii")
        with self.assertRaises(UserError):
            export(self.run_dir, output)
        self.assertEqual(output.read_text(), "original")

    def test_absent_atlas_result_is_exported_as_a_final_retrieval_outcome(self):
        calls = []

        def fetch(request):
            calls.append(1)
            if len(calls) == 1:
                raise FetchError("not_found")
            return self.fetch(request)

        result = run(self.run_dir, self.cache.directory, fetch, min_interval=0)
        self.assertEqual(result["state"], "completed")
        output = self.root / "absent.json"
        export(self.run_dir, output, format="json")
        data = json.loads(output.read_text())
        self.assertFalse(data["manifest"]["partial_snapshot"])
        missing = [row for row in data["rows"] if row["job_state"] == "not_found"]
        self.assertTrue(missing)
        self.assertTrue(all(row["summary"] is None and row["reason"] == "not_found"
                            for row in missing))

    def test_export_keeps_snapshot_while_runner_completes(self):
        self.finish(max_jobs=1)
        original = SharedCache.read_response
        ran = []
        def read(cache, response_id):
            if not ran:
                ran.append(True)
                self.finish()
            return original(cache, response_id)
        output = self.root / "snapshot.json"
        with patch.object(SharedCache, "read_response", read):
            export(self.run_dir, output, format="json", allow_partial=True)
        result = json.loads(output.read_text())
        self.assertTrue(result["manifest"]["partial_snapshot"])
        self.assertIn("pending", [row["job_state"] for row in result["rows"]])

    def test_failure_before_publication_leaves_no_final_output(self):
        self.finish()
        output = self.root / "failure.tsv"
        with patch("alphagenome_variant_scorer.export.os.link", side_effect=OSError("synthetic")):
            with self.assertRaises(OSError):
                export(self.run_dir, output)
        self.assertFalse(output.exists())
        self.assertFalse(Path(str(output) + ".manifest.json").exists())

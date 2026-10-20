"""Streaming exports from a read-only run snapshot and immutable cache responses."""

from collections import OrderedDict
from contextlib import closing
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import sqlite3
import tempfile

from filelock import FileLock

from .aggregate import AGGREGATION_VERSION, long_scores, summarize
from .cache import SharedCache
from .errors import UserError

FIELDS = ["input_state", "job_state", "reason", "api", "variant_kind", "response_id",
          "evaluation", "code", "gene_count", "track_count", "valid_genes", "invalid_genes",
          "hit_genes", "threshold_result", "genes_json", "original_values_json"]


def dump(value):
    return json.dumps(value, ensure_ascii=False, allow_nan=False, separators=(",", ":"))


def export(run_dir: Path, output: Path, *, format="tsv", threshold=0.5, allow_partial=False):
    if format not in {"tsv", "json", "scores-tsv"}:
        raise UserError("Unknown export format.")
    if not math.isfinite(threshold) or threshold < 0:
        raise UserError("Threshold must be finite and non-negative.")
    database = run_dir / "run.private.sqlite"
    if not database.is_file():
        raise UserError("No prepared run database found.")
    sidecar = Path(str(output) + ".manifest.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    with FileLock(str(output) + ".lock", timeout=0), closing(
        sqlite3.connect(database.resolve().as_uri() + "?mode=ro", uri=True)
    ) as db:
        db.execute("PRAGMA query_only=ON")
        db.execute("BEGIN")  # One consistent WAL snapshot while a runner may continue.
        version = db.execute("PRAGMA user_version").fetchone()[0]
        if version not in {1, 2}:
            raise UserError("Unsupported run database schema.")
        manifest = json.loads(db.execute("SELECT value FROM meta WHERE key='manifest'").fetchone()[0])
        if manifest["state"] in {"preparing", "preparation_failed"}:
            raise UserError("Preparation must complete before export.")
        incomplete = db.execute("SELECT count(*) FROM jobs WHERE state NOT IN "
                                "('completed','not_found','not_returned')").fetchone()[0]
        partial = manifest["state"] != "completed" or bool(incomplete)
        if partial and not allow_partial:
            raise UserError("Run is incomplete; use --allow-partial to export a snapshot.")
        if output.exists() or sidecar.exists():
            raise UserError("Export output or manifest already exists; choose a new path.")
        columns = manifest["input_columns"]
        extras = ["AGVS_" + field for field in FIELDS]
        if format == "tsv" and set(columns) & set(extras):
            raise UserError("Input column collides with an export column; use JSON output.")
        metadata = {
            "export_schema": 1, "aggregation_version": AGGREGATION_VERSION,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "run_state": manifest["state"], "partial_snapshot": partial,
            "incomplete_jobs": incomplete, "format": format, "threshold": threshold,
            "input_columns": columns, "input_sha256": manifest["input_sha256"],
            "prepared_requests_sha256": manifest["prepared_requests_sha256"],
            "input_rows": db.execute("SELECT count(*) FROM rows").fetchone()[0],
            "row_states": dict(db.execute("SELECT state,count(*) FROM rows GROUP BY state")),
            "job_states": dict(db.execute("SELECT state,count(*) FROM jobs GROUP BY state")),
            "scope": "returned genes and all returned tracks; raw scores only",
            "scores_tsv_scope": "completed response cells; other input rows have no score cells",
        }
        cache = SharedCache(Path(manifest["cache_directory"])) if "cache_directory" in manifest else None
        saved = OrderedDict()

        def payload_and_summary(response_id):
            if response_id not in saved:
                if cache is None:
                    raise UserError("Completed job has no shared cache location.")
                payload = cache.read_response(response_id).payload
                saved[response_id] = (payload, summarize(payload, threshold))
                # Bound memory even when one response is unusually large.
                while len(saved) > 1 and (len(saved) > 16 or sum(len(v[0]) for v in saved.values()) > 16_000_000):
                    saved.popitem(last=False)
            saved.move_to_end(response_id)
            return saved[response_id]

        if version == 2:
            query = """SELECT r.row_number,r.original_json,r.state,r.reason,r.variant_kind,
                j.state,j.api,x.response_id,x.error_code
                FROM rows r LEFT JOIN jobs j ON r.job_id=j.id
                LEFT JOIN job_results x ON j.id=x.job_id ORDER BY r.row_number"""
        else:
            query = """SELECT r.row_number,r.original_json,r.state,r.reason,r.variant_kind,
                j.state,j.api,NULL,NULL FROM rows r LEFT JOIN jobs j ON r.job_id=j.id
                ORDER BY r.row_number"""
        temp = temp_meta = None
        published_meta = False
        published_output = False
        written = 0
        try:
            with tempfile.NamedTemporaryFile(
                mode="w", encoding="utf-8", newline="", dir=output.parent,
                prefix="." + output.name + ".", suffix=".partial", delete=False
            ) as handle:
                temp = Path(handle.name)
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                if format == "tsv":
                    writer.writerow(columns + extras)
                elif format == "scores-tsv":
                    writer.writerow(["input_row", "response_id", "gene_index", "track_index",
                                     "gene_id", "gene_name", "gene_valid", "gene_missing_reason",
                                     "raw", "value_state", "gene_metadata_json", "track_metadata_json"])
                else:
                    handle.write('{"manifest":' + dump(metadata) + ',"rows":[')
                first = True
                for row_number, original_json, input_state, reason, kind, job_state, api, response_id, error in db.execute(query):
                    original = json.loads(original_json)
                    record = dict(input_row=row_number, original=original, input_state=input_state,
                                  job_state=job_state, reason=reason or error, api=api,
                                  variant_kind=kind, response_id=response_id, summary=None)
                    payload = summary = None
                    if job_state == "completed":
                        if not response_id:
                            raise UserError("Completed job is missing its response reference.")
                        payload, summary = payload_and_summary(response_id)
                        record["summary"] = summary
                    if format == "json":
                        if not first:
                            handle.write(",")
                        handle.write(dump(record))
                        first = False
                        written += 1
                    elif format == "tsv":
                        values = {**record, **(summary or {})}
                        values["genes_json"] = dump(summary["genes"]) if summary else ""
                        values["original_values_json"] = original_json
                        # Original JSON retains extra cells in malformed source rows.
                        cells = original[:len(columns)] + [""] * max(0, len(columns) - len(original))
                        writer.writerow(cells + [values.get(name) for name in FIELDS])
                        written += 1
                    elif payload is not None:
                        for i, j, gene, track, value, state in long_scores(payload):
                            detail = summary["genes"][i]
                            writer.writerow([row_number, response_id, i, j, gene["gene_id"],
                                             gene.get("gene_name"), detail["valid"], detail["missing_reason"],
                                             value, state, dump(gene), dump(track)])
                            written += 1
                if format == "json":
                    handle.write("]}\n")
                handle.flush()
                os.fsync(handle.fileno())
            digest = hashlib.sha256()
            with temp.open("rb") as handle:
                for block in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(block)
            with tempfile.NamedTemporaryFile(
                mode="w", encoding="utf-8", dir=output.parent,
                prefix="." + output.name + ".", suffix=".partial", delete=False
            ) as handle:
                temp_meta = Path(handle.name)
                handle.write(dump({**metadata, "written_records": written, "output_sha256": digest.hexdigest()}) + "\n")
                handle.flush()
                os.fsync(handle.fileno())
            # Publish the manifest first; the data path is the final completion marker.
            # Hard links provide atomic no-clobber publication on the same filesystem.
            os.link(temp_meta, sidecar)
            published_meta = True
            os.link(temp, output)
            published_output = True
        except BaseException:
            if published_meta and not published_output:
                sidecar.unlink()
            raise
        finally:
            if temp is not None:
                temp.unlink(missing_ok=True)
            if temp_meta is not None:
                temp_meta.unlink(missing_ok=True)
        return {"output": str(output), "manifest": str(sidecar), "format": format,
                "written_records": written, "partial_snapshot": partial, "api_calls": 0}

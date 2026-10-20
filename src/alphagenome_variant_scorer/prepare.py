"""Streaming offline preparation into a new, private run directory."""

from collections import Counter
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import sqlite3

from . import __version__
from .errors import UserError
from .inputs import read_tsv
from .reference import Reference, check_reference
from .variants import PREPARATION_VERSION, RowIssue, prepare_variant

DB_NAME = "run.private.sqlite"
MANIFEST_NAME = "manifest.private.json"
SCHEMA_VERSION = 1


def canonical(data):
    return json.dumps(data, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False)


def digest(data):
    return hashlib.sha256(canonical(data).encode("utf-8")).hexdigest()


def file_hash(path):
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def signature(path):
    st = path.stat()
    return st.st_size, st.st_mtime_ns


def save_manifest(path, data):
    temp = path.with_suffix(path.suffix + ".partial")
    with temp.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(canonical(data) + "\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temp, path)


def counts(db):
    return {
        "rows": dict(db.execute("SELECT state,count(*) FROM rows GROUP BY state")),
        "jobs": dict(db.execute("SELECT api,count(*) FROM jobs GROUP BY api")),
    }


def prepare(
    input_path: Path, run_dir: Path, fasta: Path, index: Path, *,
    input_format: str, dup_model: str | None = None, snv_source: str = "atlas",
    setting_sources: dict | None = None,
):
    if input_format not in {"literal-tsv", "sv-tsv"} or snv_source not in {"atlas", "prediction"}:
        raise UserError("Invalid preparation options.")
    if dup_model not in {None, "tandem-extra-copy"}:
        raise UserError("Unknown duplication model.")
    if not input_path.is_file():
        raise UserError("TSV input is not a regular file.")
    sdk_version = importlib.metadata.version("alphagenome")
    if sdk_version != "0.9.0":
        raise UserError("Preparation requires AlphaGenome SDK 0.9.0.")
    if run_dir.exists():
        raise UserError("Run directory already exists; choose a new directory.")
    # Validate headers before creating output or hashing a large reference.
    with read_tsv(input_path, input_format):
        pass
    reference_report = check_reference(fasta, index, create=False, hash_file=True)
    reference_signature = signature(fasta)
    if reference_signature != (reference_report["size_bytes"], reference_report["mtime_ns"]):
        raise UserError("FASTA changed after reference validation.")
    index_signature = signature(index)
    index_sha = file_hash(index)
    try:
        run_dir.mkdir(parents=True, exist_ok=False)
    except FileExistsError as exc:
        raise UserError("Run directory already exists; choose a new directory.") from exc
    manifest = {
        "schema_version": SCHEMA_VERSION, "preparation_version": PREPARATION_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(), "state": "preparing",
        "program_version": __version__, "sdk_version": sdk_version,
        "input_format": input_format, "input_path": str(input_path),
        "reference": reference_report | {"index_sha256": index_sha},
        "setting_sources": setting_sources or {},
        "options": {"snv_source": snv_source, "dup_model": dup_model},
        "api_calls": 0,
    }
    db = None
    try:
        save_manifest(run_dir / MANIFEST_NAME, manifest)
        # Copy only this explicitly selected input, never any surrounding folder.
        source_signature = signature(input_path)
        input_digest = hashlib.sha256()
        snapshot = run_dir / "input.private.tsv"
        with input_path.open("rb") as src, snapshot.open("xb") as dst:
            for block in iter(lambda: src.read(1024 * 1024), b""):
                dst.write(block)
                input_digest.update(block)
        if source_signature != signature(input_path):
            raise UserError("TSV input changed during copying.")
        manifest["input_sha256"] = input_digest.hexdigest()
        db = sqlite3.connect(run_dir / DB_NAME)
        # Preparation has a single writer and needs no concurrent external readers.
        db.execute("PRAGMA foreign_keys=ON")
        db.execute("PRAGMA synchronous=FULL")
        db.executescript("""
            PRAGMA user_version=1;
            CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT NOT NULL);
            CREATE TABLE jobs(
                id TEXT PRIMARY KEY, api TEXT NOT NULL, request_json TEXT NOT NULL,
                state TEXT NOT NULL DEFAULT 'pending');
            CREATE TABLE rows(
                row_number INTEGER PRIMARY KEY, source_line INTEGER NOT NULL,
                original_json TEXT NOT NULL, state TEXT NOT NULL,
                reason TEXT, variant_kind TEXT, job_id TEXT REFERENCES jobs(id),
                details_json TEXT NOT NULL);
            CREATE INDEX row_state ON rows(state);
        """)
        db.execute("INSERT INTO meta VALUES('manifest',?)", (canonical(manifest),))
        db.commit()
        reasons = Counter()
        reference_key = {
            "assembly": "GRCh38", "fasta_sha256": reference_report["sha256"],
            "sdk_version": sdk_version, "preparation_version": PREPARATION_VERSION,
        }
        with Reference(fasta, index) as reference, read_tsv(snapshot, input_format) as (header, reader):
            manifest["input_columns"] = header
            for row_number, values in enumerate(reader, start=2):
                kind = None
                job_id = None
                details = {}
                reason = None
                state = "ready"
                try:
                    if len(values) != len(header):
                        raise RowIssue("invalid", "column_count_mismatch")
                    kind, request, details = prepare_variant(
                        dict(zip(header, values)), input_format, reference,
                        dup_model=dup_model, snv_source=snv_source,
                    )
                    job_id = digest({
                        "request": request, "reference": reference_key,
                        "structure_assumption": details["structure_assumption"],
                    })
                    db.execute("INSERT OR IGNORE INTO jobs(id,api,request_json) VALUES(?,?,?)",
                               (job_id, request["api"], canonical(request)))
                except RowIssue as issue:
                    state, reason = issue.state, issue.reason
                    reasons[reason] += 1
                db.execute(
                    "INSERT INTO rows VALUES(?,?,?,?,?,?,?,?)",
                    (row_number, reader.line_num, canonical(values), state, reason,
                     kind, job_id, canonical(details)),
                )
                if row_number % 1000 == 0:
                    db.commit()
        if reference_signature != signature(fasta) or index_signature != signature(index):
            raise UserError("Reference or index changed during preparation.")
        if file_hash(snapshot) != manifest["input_sha256"]:
            raise UserError("Input snapshot changed during preparation.")
        manifest.update(state="prepared", counts=counts(db), reasons=dict(reasons))
        request_digest = hashlib.sha256()
        for job_id, request_json in db.execute("SELECT id,request_json FROM jobs ORDER BY id"):
            request_digest.update((job_id + "\n" + request_json + "\n").encode("utf-8"))
        manifest["prepared_requests_sha256"] = request_digest.hexdigest()
        db.execute("UPDATE meta SET value=? WHERE key='manifest'", (canonical(manifest),))
        db.commit()
        save_manifest(run_dir / MANIFEST_NAME, manifest)
        return {
            "run_dir": str(run_dir), "state": "prepared",
            "counts": manifest["counts"], "reasons": manifest["reasons"], "api_calls": 0,
        }
    except BaseException:
        manifest["state"] = "preparation_failed"
        if db is not None:
            try:
                db.rollback()
                db.execute("UPDATE meta SET value=? WHERE key='manifest'", (canonical(manifest),))
                db.commit()
            except sqlite3.Error:
                pass  # Preserve the original failure if storage itself is unavailable.
        try:
            save_manifest(run_dir / MANIFEST_NAME, manifest)
        except OSError:
            pass
        raise
    finally:
        if db is not None:
            db.close()

"""Command line entry point; only run/resume can call the scoring APIs."""

from __future__ import annotations

import argparse
import json
import math
import sqlite3
import sys

from filelock import Timeout

from . import __version__
from .config import absolute_path, config_path, load, resolved_reference, save_setting, select_reference
from .errors import UserError
from .reference import check_reference


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(
        prog="alphagenome-variant-scorer",
        description="Reference checks, variant preparation and cached RNA-seq scoring.",
    )
    root.add_argument("--version", action="version", version=__version__)
    root.add_argument("--config", metavar="FILE", help="Use this config instead of the user config.")
    commands = root.add_subparsers(dest="command", required=True)
    config = commands.add_parser("config", help="Persist and inspect local settings.")
    settings = config.add_subparsers(dest="operation", required=True)
    change = settings.add_parser("set", help="Save a setting; does not create a FASTA index.")
    change.add_argument("key", choices=["reference.fasta", "reference.index", "reference.assembly",
                                        "cache.directory"])
    change.add_argument("value")
    settings.add_parser("show", help="Show supported settings and their source as JSON.")
    cache = commands.add_parser("cache", help="Initialize or inspect a local shared cache.")
    cache_ops = cache.add_subparsers(dest="operation", required=True)
    for operation in ("init", "status"):
        item = cache_ops.add_parser(operation)
        item.add_argument("--cache-dir", metavar="DIR")
    reference = commands.add_parser("reference", help="Inspect an external reference FASTA.")
    operations = reference.add_subparsers(dest="operation", required=True)
    check = operations.add_parser("check")
    check.add_argument("--reference-fasta", metavar="PATH")
    check.add_argument("--index", metavar="PATH", help="Override the index path for this check.")
    check.add_argument("--create-index", action="store_true", help="Create a missing index; never replace one.")
    check.add_argument("--sha256", action="store_true", help="Read the full FASTA to compute its SHA-256.")
    prepare = commands.add_parser("prepare", help="Validate TSV and save API requests without sending them.")
    prepare.add_argument("--input", required=True, metavar="TSV")
    prepare.add_argument("--input-format", required=True, choices=["literal-tsv", "sv-tsv"])
    prepare.add_argument("--run-dir", required=True, metavar="DIR")
    prepare.add_argument("--reference-fasta", metavar="PATH")
    prepare.add_argument("--index", metavar="PATH")
    prepare.add_argument("--snv-source", choices=["atlas", "prediction"], default="atlas")
    prepare.add_argument("--dup-model", choices=["tandem-extra-copy"])
    status = commands.add_parser("status", help="Read a run's status without config, FASTA or API keys.")
    status.add_argument("--run-dir", required=True, metavar="DIR")
    exporting = commands.add_parser("export", help="Export saved scores without config, FASTA or API keys.")
    exporting.add_argument("--run-dir", required=True, metavar="DIR")
    exporting.add_argument("--output", required=True, metavar="FILE")
    exporting.add_argument("--format", choices=["tsv", "json", "scores-tsv"], default="tsv")
    exporting.add_argument("--threshold", type=float, default=0.5)
    exporting.add_argument("--allow-partial", action="store_true")
    for name in ("run", "resume"):
        scoring = commands.add_parser(name, help="Score prepared requests, reusing shared responses.")
        scoring.add_argument("--run-dir", required=True, metavar="DIR")
        scoring.add_argument("--cache-dir", metavar="DIR")
        scoring.add_argument("--max-jobs", type=int)
        scoring.add_argument("--max-attempts", type=int, default=3)
        scoring.add_argument("--workers", type=int, default=1,
                             help="Concurrent cache/score workers (1-128; default 1). "
                                  "128 is a program limit, not a recommendation.")
        scoring.add_argument("--min-interval", type=float, default=1.0)
        scoring.add_argument("--timeout", type=float, default=120.0)
    return root


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        if args.command == "export":
            from .export import export

            result = export(
                absolute_path(args.run_dir), absolute_path(args.output), format=args.format,
                threshold=args.threshold, allow_partial=args.allow_partial)
        elif args.command == "status":
            from .store import run_status

            result = run_status(absolute_path(args.run_dir))
        elif args.command in {"cache", "run", "resume"}:
            from .cache import SharedCache

            if args.cache_dir:
                directory = absolute_path(args.cache_dir)
            else:
                path = config_path(args.config)
                data = load(path, missing_ok=args.config is None)
                value = data.get("cache", {}).get("directory")
                if value is None:
                    raise UserError("Set cache.directory or supply --cache-dir DIR.")
                directory = absolute_path(value, base=path.parent)
            cache = SharedCache(directory)
            if args.command == "cache":
                result = cache.initialize() if args.operation == "init" else cache.status()
            else:
                from .api import ScoreClient
                from .auth import resolve_key
                from .runner import run

                if not math.isfinite(args.timeout) or args.timeout <= 0:
                    raise UserError("timeout must be a positive finite number.")
                client = ScoreClient(
                    key_provider=lambda: resolve_key(
                        config_file=config_path(args.config),
                        explicit_config=args.config is not None),
                    timeout=args.timeout)
                try:
                    result = run(absolute_path(args.run_dir), directory, client,
                                 max_jobs=args.max_jobs, max_attempts=args.max_attempts,
                                 min_interval=args.min_interval, workers=args.workers)
                finally:
                    client.close()
        elif args.command == "config":
            path = config_path(args.config)
            if args.operation == "set":
                save_setting(path, args.key, args.value)
                result = {"saved": args.key, "config_file": str(path)}
            else:
                data = load(path, missing_ok=args.config is None)
                result = {
                    "config_file": str(path), "exists": path.exists(),
                    "schema_version": data["schema_version"],
                    "reference": resolved_reference(data, path),
                    "cache": {
                        "directory": str(absolute_path(data["cache"]["directory"], base=path.parent))
                    } if data.get("cache", {}).get("directory") else {},
                    "auth": {
                        "api_key": "[redacted]" if data.get("auth", {}).get("api_key", "").strip() else None,
                    },
                }
        else:
            path = config_path(args.config)
            fasta, index, sources = select_reference(
                path, explicit_config=args.config is not None,
                fasta=args.reference_fasta, index=args.index,
            )
            if args.command == "prepare":
                from .prepare import prepare

                result = prepare(
                    absolute_path(args.input), absolute_path(args.run_dir), fasta, index,
                    input_format=args.input_format, dup_model=args.dup_model,
                    snv_source=args.snv_source, setting_sources=sources,
                )
            else:
                result = check_reference(
                    fasta, index, create=args.create_index, hash_file=args.sha256,
                )
                result["setting_sources"] = sources
                result["config_file"] = str(path)
        print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
        return 1 if result.get("state") == "failed" else 0
    except UserError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except Timeout:
        print("error: another process holds the operation lock; retry later.", file=sys.stderr)
        return 2
    except OSError as exc:
        # Do not echo arbitrary file contents or OS exception payloads.
        print(f"error: filesystem operation failed ({type(exc).__name__}). Check paths and permissions.",
              file=sys.stderr)
        return 2
    except sqlite3.Error:
        print("error: database operation failed. Check database existence, disk space and permissions.",
              file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print("Interrupted.", file=sys.stderr)
        return 130

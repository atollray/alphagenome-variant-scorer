"""Read-only status. No SDK, reference, config, or credentials are needed."""

import json
from contextlib import closing
from pathlib import Path
import sqlite3

from .errors import UserError


def run_status(run_dir: Path) -> dict:
    path = run_dir / "run.private.sqlite"
    if not path.is_file():
        raise UserError("No prepared database found in this run directory.")
    try:
        with closing(sqlite3.connect(path.resolve().as_uri() + "?mode=ro", uri=True)) as db:
            db.execute("PRAGMA query_only=ON")
            version = db.execute("PRAGMA user_version").fetchone()[0]
            if version not in {1, 2}:
                raise UserError("Unsupported run database schema.")
            saved = db.execute("SELECT value FROM meta WHERE key='manifest'").fetchone()
            if saved is None:
                raise UserError("Run manifest is missing.")
            manifest = json.loads(saved[0])
            result = {
                "run_dir": str(run_dir), "state": manifest["state"],
                "rows": dict(db.execute("SELECT state,count(*) FROM rows GROUP BY state")),
                "jobs": [dict(api=api, state=state, count=n) for api, state, n in db.execute(
                    "SELECT api,state,count(*) FROM jobs GROUP BY api,state")],
                "reasons": dict(db.execute(
                    "SELECT reason,count(*) FROM rows WHERE reason IS NOT NULL GROUP BY reason")),
                "sdk_version": manifest["sdk_version"],
                "runtime": manifest.get("runtime"),
            }
            if version == 2:
                attempts, hits, fetched = db.execute(
                    "SELECT COALESCE(sum(attempts),0),"
                    "COALESCE(sum(cache_hit=1),0),COALESCE(sum(cache_hit=0),0) FROM job_results"
                ).fetchone()
                result["scoring"] = {
                    "fetch_attempts": attempts, "cache_hits": hits, "fetched_jobs": fetched,
                    "errors": dict(db.execute(
                        "SELECT error_code,count(*) FROM job_results "
                        "WHERE error_code IS NOT NULL GROUP BY error_code")),
                }
            return result
    except (sqlite3.Error, ValueError, KeyError) as exc:
        raise UserError("Invalid or incomplete run database.") from exc

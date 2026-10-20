"""Local shared response store. No legacy import, SDK client, or credentials."""

from contextlib import closing, contextmanager
from concurrent.futures import CancelledError
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sqlite3
import uuid

from filelock import FileLock, Timeout

from .errors import UserError

APPLICATION_ID = 0x41475653
SCHEMA = 1


@contextmanager
def request_lock(path, stop_event=None):
    lock = FileLock(str(path))
    while True:
        if stop_event is not None and stop_event.is_set():
            raise CancelledError()
        try:
            lock.acquire(timeout=0.1)
            break
        except Timeout:
            continue
    try:
        yield
    finally:
        lock.release()


def canonical(value: dict) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
                      allow_nan=False)


def identity(request: dict, response_format: str) -> tuple[str, str]:
    """Only pass effective API request fields, never run IDs or local provenance."""
    if not request or not isinstance(request, dict) or not response_format.strip():
        raise UserError("Cache identity requires a request and versioned response format.")
    encoded = canonical({"key_version": 1, "request": request,
                         "response_format": response_format})
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest(), encoded


@dataclass(frozen=True)
class Response:
    response_id: str
    payload: bytes
    provenance: dict
    created_at: str


class SharedCache:
    def __init__(self, directory: Path):
        self.directory = directory.resolve()
        self.path = self.directory / "scores.private.sqlite"

    def _connect(self, *, readonly=False):
        mode = "ro" if readonly else "rw"
        db = sqlite3.connect(self.path.as_uri() + f"?mode={mode}", uri=True, timeout=30)
        if (db.execute("PRAGMA application_id").fetchone()[0] != APPLICATION_ID
                or db.execute("PRAGMA user_version").fetchone()[0] != SCHEMA):
            db.close()
            raise UserError("Unknown shared cache database or unsupported schema.")
        if readonly:
            db.execute("PRAGMA query_only=ON")
        return db

    def initialize(self):
        """Create only a new dedicated directory, or validate an existing cache."""
        self.directory.parent.mkdir(parents=True, exist_ok=True)
        with FileLock(str(self.directory) + ".init.lock", timeout=30):
            if self.directory.exists():
                # Never initialize an arbitrary existing folder or database.
                return self.status()
            self.directory.mkdir()
            with closing(sqlite3.connect(self.path)) as db:
                db.execute("PRAGMA journal_mode=WAL")
                db.execute("PRAGMA synchronous=FULL")
                with db:
                    db.execute("""CREATE TABLE responses (
                        response_id TEXT PRIMARY KEY, cache_key TEXT NOT NULL,
                        identity TEXT NOT NULL, payload BLOB NOT NULL,
                        payload_sha256 TEXT NOT NULL, provenance TEXT NOT NULL,
                        created_at TEXT NOT NULL, sequence INTEGER NOT NULL UNIQUE)""")
                    db.execute("CREATE INDEX response_key ON responses(cache_key, sequence)")
                    db.execute(f"PRAGMA application_id={APPLICATION_ID}")
                    db.execute(f"PRAGMA user_version={SCHEMA}")
            (self.directory / "locks").mkdir()
        return self.status()

    def status(self):
        with closing(self._connect(readonly=True)) as db:
            count, keys = db.execute(
                "SELECT count(*), count(DISTINCT cache_key) FROM responses").fetchone()
        return {"directory": str(self.directory), "schema_version": SCHEMA,
                "responses": count, "requests": keys}

    def lookup(self, request: dict, response_format: str) -> Response | None:
        key, encoded = identity(request, response_format)
        with closing(self._connect(readonly=True)) as db:
            row = db.execute(
                """SELECT response_id,payload,provenance,created_at,payload_sha256,identity
                   FROM responses WHERE cache_key=? ORDER BY sequence DESC LIMIT 1""",
                (key,)).fetchone()
        if row is None:
            return None
        if row[5] != encoded or hashlib.sha256(row[1]).hexdigest() != row[4]:
            raise UserError("Shared cache integrity check failed; response was not reused.")
        return Response(row[0], row[1], json.loads(row[2]), row[3])

    def read_response(self, response_id: str) -> Response:
        """Read the exact immutable response referenced by a run."""
        with closing(self._connect(readonly=True)) as db:
            row = db.execute(
                """SELECT response_id,payload,provenance,created_at,payload_sha256
                   FROM responses WHERE response_id=?""", (response_id,)).fetchone()
        if row is None:
            raise UserError("Referenced cache response does not exist.")
        if hashlib.sha256(row[1]).hexdigest() != row[4]:
            raise UserError("Shared cache integrity check failed.")
        return Response(row[0], row[1], json.loads(row[2]), row[3])

    def verify_references(self, response_ids):
        """Check run references without loading every historical matrix on resume."""
        with closing(self._connect(readonly=True)) as db:
            for response_id in response_ids:
                if db.execute("SELECT 1 FROM responses WHERE response_id=?",
                              (response_id,)).fetchone() is None:
                    raise UserError("Referenced cache response does not exist.")

    def get_or_fetch(self, request: dict, response_format: str, fetch, *, refresh=False,
                     stop_event=None):
        """fetch returns (encoded success payload, provenance); errors must raise.

        Payload codec must preserve gene/track labels and missing-value states.
        A successful response with missing scores is still a success. A failed
        call is never persisted. Returns (response, cache_hit).
        """
        key, encoded = identity(request, response_format)
        # Validate before creating lock files; no DB transaction spans the RPC.
        with closing(self._connect(readonly=True)):
            pass
        with request_lock(self.directory / "locks" / (key + ".lock"), stop_event):
            if not refresh:
                found = self.lookup(request, response_format)
                if found is not None:
                    return found, True
            payload, provenance = fetch()
            if not isinstance(payload, bytes) or not isinstance(provenance, dict):
                raise UserError("Cache fetch must return bytes and provenance metadata.")
            metadata = canonical(provenance)
            result = Response(uuid.uuid4().hex, payload, json.loads(metadata),
                              datetime.now(timezone.utc).isoformat())
            with closing(self._connect()) as db:
                db.execute("PRAGMA synchronous=FULL")
                with db:
                    db.execute("""INSERT INTO responses
                        SELECT ?,?,?,?,?,?,?,COALESCE(MAX(sequence),0)+1 FROM responses""",
                        (result.response_id, key, encoded, payload,
                         hashlib.sha256(payload).hexdigest(), metadata, result.created_at))
            return result, False

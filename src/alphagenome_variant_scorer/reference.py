"""Read-only FASTA extraction and explicitly requested index creation."""

from __future__ import annotations

import hashlib
from functools import cached_property
import os
from pathlib import Path
import tempfile
import warnings

from filelock import FileLock, Timeout
from pyfaidx import Fasta, FastaIndexingError, FetchError

from .errors import UserError


def inspect_header(path: Path) -> None:
    """Cheap format check, without creating an index or reading the whole genome."""
    if not path.is_file():
        raise UserError("FASTA file does not exist or is not a regular file.")
    if path.suffix.lower() in {".gz", ".bgz", ".bgzf", ".zip", ".bz2"}:
        raise UserError("This version supports uncompressed FASTA only.")
    with path.open("rb") as handle:
        header = handle.readline(65537)
        if len(header) > 65536 or not header.startswith(b">") or not header[1:].strip():
            raise UserError("Expected an uncompressed FASTA with a named first record.")
        try:
            header.decode("ascii")
        except UnicodeError as exc:
            raise UserError("FASTA record names must be ASCII.") from exc


def _index_rows(index: Path, fasta_size: int) -> dict[str, int]:
    rows: dict[str, int] = {}
    try:
        with index.open(encoding="ascii") as handle:
            for line in handle:
                parts = line.rstrip("\r\n").split("\t")
                if len(parts) != 5 or not parts[0] or parts[0] in rows:
                    raise ValueError
                length, offset, bases, width = map(int, parts[1:])
                if length <= 0 or offset <= 0 or bases <= 0 or width < bases:
                    raise ValueError
                last_byte = offset + (length - 1) // bases * width + (length - 1) % bases
                if last_byte >= fasta_size:
                    raise ValueError
                rows[parts[0]] = length
    except (ValueError, UnicodeError) as exc:
        raise UserError("Invalid or incompatible FASTA index.") from exc
    if not rows:
        raise UserError("FASTA index has no sequences.")
    return rows


def ensure_index(fasta: Path, index: Path, *, create: bool = False) -> None:
    inspect_header(fasta)
    if fasta.resolve() == index.resolve():
        raise UserError("Index path must differ from the FASTA path.")
    if index.exists():
        if not index.is_file():
            raise UserError("Index path is not a regular file.")
        if fasta.stat().st_mtime_ns > index.stat().st_mtime_ns:
            raise UserError("Index is older than FASTA; create an index at a new path.")
        _index_rows(index, fasta.stat().st_size)
        return
    if not create:
        raise UserError("FASTA index is missing; use --create-index or --index PATH.")
    index.parent.mkdir(parents=True, exist_ok=True)
    try:
        with FileLock(str(index) + ".lock", timeout=0):
            if index.exists():
                ensure_index(fasta, index)
                return
            initial = fasta.stat()
            # Generate in the destination filesystem, never in the FASTA itself.
            with tempfile.TemporaryDirectory(prefix=".fai-", dir=index.parent) as temp:
                candidate = Path(temp) / "reference.fai"
                try:
                    with Fasta(
                        str(fasta), indexname=str(candidate), mutable=False,
                        build_index=True, rebuild=True, duplicate_action="stop",
                    ):
                        pass
                except (ValueError, RuntimeError, FastaIndexingError) as exc:
                    raise UserError("Unable to build a FASTA index; check the FASTA format.") from exc
                current = fasta.stat()
                if (initial.st_size, initial.st_mtime_ns) != (current.st_size, current.st_mtime_ns):
                    raise UserError("FASTA changed while indexing; retry with an unchanged reference.")
                _index_rows(candidate, current.st_size)
                # Atomic no-clobber publication. Existing indices are never replaced.
                try:
                    os.link(candidate, index)
                except FileExistsError as exc:
                    raise UserError("Index appeared during creation; inspect it and retry.") from exc
    except Timeout as exc:
        raise UserError("Another process is creating this index; retry later.") from exc


class Reference:
    """One handle per worker; coordinates are 0-based, half-open."""

    def __init__(self, fasta: Path, index: Path):
        ensure_index(fasta, index)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("error")
                self._fasta = Fasta(
                    str(fasta), indexname=str(index), as_raw=True,
                    strict_bounds=True, sequence_always_upper=True,
                    mutable=False, build_index=False, rebuild=False,
                    duplicate_action="stop",
                )
        except (ValueError, RuntimeError, Warning) as exc:
            raise UserError("Unable to open the FASTA with this index.") from exc

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def close(self) -> None:
        self._fasta.close()

    @cached_property
    def contigs(self) -> dict[str, int]:
        return {name: len(self._fasta[name]) for name in self._fasta.keys()}

    def sequence(self, chromosome: str, start: int, end: int) -> str:
        if type(start) is not int or type(end) is not int:
            raise UserError("Coordinates must be integers.")
        if chromosome not in self._fasta:
            raise UserError("Unknown contig; contig aliases are not inferred.")
        length = len(self._fasta[chromosome])
        if not 0 <= start <= end <= length:
            raise UserError("Interval is outside the contig (expected 0-based, half-open bounds).")
        if start == end:
            return ""
        try:
            value = str(self._fasta[chromosome][start:end])
        except (ValueError, FetchError) as exc:
            raise UserError("Unable to read this interval with the supplied index.") from exc
        if len(value) != end - start or set(value) - set("ACGTN"):
            raise UserError("Reference region has invalid bases or an incompatible index.")
        return value

    def extract(self, interval) -> str:
        """Adapter for an SDK-like interval without importing the SDK."""
        return self.sequence(interval.chromosome, interval.start, interval.end)


def check_reference(fasta: Path, index: Path, *, create: bool, hash_file: bool) -> dict:
    ensure_index(fasta, index, create=create)
    before = fasta.stat()
    with Reference(fasta, index) as ref:
        contigs = ref.contigs
        for name, length in contigs.items():
            ref.sequence(name, 0, min(32, length))
            ref.sequence(name, max(0, length - 32), length)
    result = {
        "fasta": str(fasta),
        "index": str(index),
        "size_bytes": before.st_size,
        "mtime_ns": before.st_mtime_ns,
        "contig_count": len(contigs),
        "contigs": contigs,
        "validation": "index structure and contig endpoint reads",
        "assembly_verified": False,
        "api_reference_match_verified": False,
    }
    if hash_file:
        digest = hashlib.sha256()
        with fasta.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        result["sha256"] = digest.hexdigest()
    after = fasta.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise UserError("FASTA changed during validation.")
    return result

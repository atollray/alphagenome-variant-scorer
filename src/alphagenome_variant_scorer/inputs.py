"""Streaming TSV input contracts. Original field values are never rewritten."""

import csv
from contextlib import contextmanager

from .errors import UserError

REQUIRED = {
    "literal-tsv": ("#CHROM", "POS", "REF", "ALT"),
    "sv-tsv": ("SV_chrom", "SV_start", "SV_end", "SV_type", "SV_length"),
}


@contextmanager
def read_tsv(path, input_format):
    if input_format not in REQUIRED:
        raise UserError("Unsupported input format.")
    previous_limit = csv.field_size_limit(2**21)
    try:
        with path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t", strict=True)
            header = next(reader, None)
            if not header or any(not name for name in header) or len(set(header)) != len(header):
                raise UserError("TSV must have a non-empty header with unique column names.")
            if not set(REQUIRED[input_format]).issubset(header):
                raise UserError("TSV is missing required columns for the selected input format.")
            yield header, reader
    except (UnicodeError, csv.Error) as exc:
        raise UserError("Invalid UTF-8 TSV or a field exceeds the supported size.") from exc
    finally:
        csv.field_size_limit(previous_limit)

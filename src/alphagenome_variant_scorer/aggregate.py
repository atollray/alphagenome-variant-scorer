"""RNA-seq summaries over returned genes/tracks; no pathogenicity interpretation."""

import math
import numpy as np

from .errors import UserError
from .scores import decode_scores

AGGREGATION_VERSION = "all-tracks-finite-v1"


def records(frame):
    columns, rows, index = frame["columns"], frame["data"], frame["index"]
    if len(columns) != len(set(columns)) or len(index) != len(rows):
        raise UserError("Invalid score metadata.")
    if any(len(row) != len(columns) for row in rows):
        raise UserError("Invalid score metadata row width.")
    return [dict(zip(columns, row)) for row in rows]


def summarize(payload, threshold):
    if not math.isfinite(threshold) or threshold < 0:
        raise UserError("Threshold must be finite and non-negative.")
    meta, arrays = decode_scores(payload)
    raw = arrays["raw"]
    obs, tracks = records(meta["obs"]), records(meta["tracks"])
    ids = [row.get("gene_id") for row in obs]
    if any(not isinstance(gene, str) or not gene for gene in ids) or len(set(ids)) != len(ids):
        raise UserError("Expected unique, non-empty gene_id metadata for one variant.")
    n_genes, n_tracks = raw.shape
    genes = []
    for i, row in enumerate(obs):
        values = raw[i]
        finite = np.isfinite(values)
        valid = bool(n_tracks and finite.all())
        positive = negative = False
        minimum = maximum = None
        min_track = max_track = None
        if valid:
            low, high = int(np.argmin(values)), int(np.argmax(values))
            minimum, maximum = float(values[low]), float(values[high])
            min_track, max_track = low, high
            positive, negative = maximum >= threshold, minimum <= -threshold
            missing = None
        elif not n_tracks:
            missing = "no_tracks"
        elif np.isnan(values).all():
            missing = "all_nan"
        elif finite.any():
            missing = "partial_nonfinite"
        else:
            missing = "all_nonfinite"
        genes.append({
            "gene_id": row["gene_id"], "gene_name": row.get("gene_name"),
            "valid": valid, "missing_reason": missing,
            "finite_tracks": int(finite.sum()), "nan_tracks": int(np.isnan(values).sum()),
            "positive_infinity_tracks": int(np.isposinf(values).sum()),
            "negative_infinity_tracks": int(np.isneginf(values).sum()),
            "min_raw": minimum, "max_raw": maximum,
            "min_track_index": min_track, "max_track_index": max_track,
            "positive_hit": positive, "negative_hit": negative,
        })
    valid_count = sum(gene["valid"] for gene in genes)
    hit_count = sum(gene["positive_hit"] or gene["negative_hit"] for gene in genes)
    if not n_genes:
        evaluation = "no_genes"
    elif not n_tracks:
        evaluation = "no_tracks"
    elif not valid_count:
        evaluation = "unavailable"
    elif valid_count < n_genes:
        evaluation = "partial"
    else:
        evaluation = "complete"
    return {
        "aggregation_version": AGGREGATION_VERSION, "threshold": threshold,
        "evaluation": evaluation,
        "code": "P01" if evaluation == "partial" else ("E02" if not valid_count else None),
        "gene_count": n_genes, "track_count": n_tracks,
        "valid_genes": valid_count, "invalid_genes": n_genes - valid_count,
        "hit_genes": hit_count if valid_count else None,
        "threshold_result": ("hits" if hit_count else "no_hits") if valid_count else "not_evaluated",
        "genes": genes,
    }


def long_scores(payload):
    """All raw cells, including excluded genes, with explicit nonfinite kinds."""
    meta, arrays = decode_scores(payload)
    obs, tracks = records(meta["obs"]), records(meta["tracks"])
    for i, row in enumerate(obs):
        for j, track in enumerate(tracks):
            value = float(arrays["raw"][i, j])
            state = ("finite" if math.isfinite(value) else
                     "nan" if math.isnan(value) else
                     "positive_infinity" if value > 0 else "negative_infinity")
            yield i, j, row, track, value if state == "finite" else None, state

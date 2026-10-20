"""Lossless numerical score payloads, with JSON labels and no pickle."""

import io
import json
import zipfile

import numpy as np

from .errors import UserError

RESPONSE_FORMAT = "rna-seq-npz-v1"


def encode_scores(adata) -> bytes:
    arrays = {"raw": np.asarray(adata.X)}
    names = {}
    for number, (name, values) in enumerate(adata.layers.items()):
        key = f"layer_{number}"
        names[key] = name
        arrays[key] = np.asarray(values)
    for value in arrays.values():
        if value.ndim != 2 or value.shape != adata.shape or value.dtype.kind not in "fiu":
            raise UserError("Unsupported score matrix shape or dtype.")
    obs = adata.obs.copy()
    if "variant" in obs:
        obs["variant"] = obs["variant"].map(str)
    metadata = {
        "format": RESPONSE_FORMAT, "shape": list(adata.shape), "layers": names,
        "obs": json.loads(obs.to_json(orient="split", date_format="iso")),
        "tracks": json.loads(adata.var.to_json(orient="split", date_format="iso")),
        # NaN and +/-Infinity are preserved in NPZ, never replaced with zero.
        "nonfinite": "preserved-in-numeric-arrays",
    }
    arrays["metadata"] = np.frombuffer(
        json.dumps(metadata, ensure_ascii=False, allow_nan=False).encode("utf-8"), dtype=np.uint8)
    output = io.BytesIO()
    np.savez_compressed(output, **arrays)
    return output.getvalue()


def decode_scores(payload: bytes) -> tuple[dict, dict]:
    try:
        return _decode_scores(payload)
    except (ValueError, KeyError, TypeError, IndexError, AttributeError,
            EOFError, UnicodeError, zipfile.BadZipFile):
        raise UserError("Invalid score payload or metadata.") from None


def _decode_scores(payload: bytes) -> tuple[dict, dict]:
    with np.load(io.BytesIO(payload), allow_pickle=False) as archive:
        metadata = json.loads(archive["metadata"].tobytes())
        if metadata["format"] != RESPONSE_FORMAT:
            raise UserError("Unsupported score payload format.")
        names = {"raw", *metadata["layers"]}
        arrays = {name: archive[name].copy() for name in names}
    shape = tuple(metadata["shape"])
    if len(shape) != 2 or any(value.shape != shape or value.dtype.kind not in "fiu" for value in arrays.values()):
        raise UserError("Invalid score payload shape or dtype.")
    if len(metadata["obs"]["data"]) != shape[0] or len(metadata["tracks"]["data"]) != shape[1]:
        raise UserError("Score labels do not match the matrix.")
    for name in ("obs", "tracks"):
        frame = metadata[name]
        columns, rows, index = frame["columns"], frame["data"], frame["index"]
        if (not all(isinstance(column, str) for column in columns)
                or len(columns) != len(set(columns)) or len(index) != len(rows)
                or any(len(row) != len(columns) for row in rows)):
            raise UserError("Invalid score metadata labels.")
    return metadata, arrays

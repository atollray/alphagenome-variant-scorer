"""Versioned local configuration; callers must never print unfiltered data."""

from __future__ import annotations

import os
from pathlib import Path
import tempfile
import tomllib

from filelock import FileLock, Timeout
import tomli_w

from .errors import UserError

SCHEMA_VERSION = 1
APP_NAME = "alphagenome-variant-scorer"
REFERENCE_KEYS = {"fasta", "index", "assembly"}


def absolute_path(value: str | Path, *, base: Path | None = None) -> Path:
    if not str(value).strip():
        raise UserError("Path must not be empty.")
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = (base if base is not None else Path.cwd()) / path
    return path.resolve()


def default_config_path(
    *, platform: str | None = None, environ: dict | None = None,
    home: Path | None = None,
) -> Path:
    platform = os.name if platform is None else platform
    env = os.environ if environ is None else environ
    if platform == "nt":
        root = env.get("LOCALAPPDATA")
        if not root or not Path(root).is_absolute():
            raise UserError("LOCALAPPDATA is unavailable; use --config FILE.")
        return Path(root) / APP_NAME / "config.toml"
    root = env.get("XDG_CONFIG_HOME")
    if root:
        if not Path(root).is_absolute():
            raise UserError("XDG_CONFIG_HOME must be absolute; use --config FILE.")
        folder = Path(root)
    else:
        folder = (home if home is not None else Path.home()) / ".config"
    return folder / APP_NAME / "config.toml"


def config_path(explicit: str | None) -> Path:
    # Packaged Windows applications can redirect LOCALAPPDATA. Expose the real
    # location so the same file can be opened outside the app or from WSL.
    return absolute_path(explicit) if explicit else default_config_path().resolve()


def validate(data: dict) -> dict:
    if set(data) - {"schema_version", "reference", "cache", "auth"}:
        raise UserError("Unknown config section. Supported: schema_version, reference, cache, auth.")
    if type(data.get("schema_version")) is not int or data["schema_version"] != SCHEMA_VERSION:
        raise UserError("Unsupported config schema_version; expected 1.")
    ref = data.get("reference", {})
    if not isinstance(ref, dict) or set(ref) - REFERENCE_KEYS:
        raise UserError("Unknown reference setting; supported keys: fasta, index, assembly.")
    for key, value in ref.items():
        if not isinstance(value, str) or not value.strip():
            raise UserError(f"reference.{key} must be a non-empty string.")
    if ref.get("assembly", "GRCh38") != "GRCh38":
        raise UserError("Only assembly GRCh38 is supported. No coordinate conversion is performed.")
    cache = data.get("cache", {})
    if not isinstance(cache, dict) or set(cache) - {"directory"}:
        raise UserError("Unknown cache setting; supported key: directory.")
    if "directory" in cache and (
        not isinstance(cache["directory"], str) or not cache["directory"].strip()
    ):
        raise UserError("cache.directory must be a non-empty string.")
    auth = data.get("auth", {})
    if not isinstance(auth, dict) or set(auth) - {"api_key"}:
        raise UserError("Unknown auth setting; supported key: api_key.")
    for value in auth.values():
        if not isinstance(value, str):
            raise UserError("auth.api_key must be a string.")
    return data


def load(path: Path, *, missing_ok: bool = False) -> dict:
    try:
        with path.open("rb") as handle:
            data = tomllib.load(handle)
    except FileNotFoundError as exc:
        if missing_ok:
            return {"schema_version": SCHEMA_VERSION, "reference": {"assembly": "GRCh38"}}
        raise UserError("Config file not found; use config set reference.fasta PATH.") from exc
    except (tomllib.TOMLDecodeError, UnicodeError) as exc:
        # Parser error strings can include a line containing a pasted secret.
        raise UserError("Invalid TOML config. Check its syntax and UTF-8 encoding.") from exc
    return validate(data)


def resolved_reference(data: dict, path: Path) -> dict:
    result = dict(data.get("reference", {}))
    result.setdefault("assembly", "GRCh38")
    for key in ("fasta", "index"):
        if key in result:
            result[key] = str(absolute_path(result[key], base=path.parent))
    return result


def save_setting(path: Path, key: str, value: str) -> None:
    if key not in {f"reference.{k}" for k in REFERENCE_KEYS} | {"cache.directory"}:
        raise UserError("Unsupported setting. Edit auth.api_key in the private config file.")
    section = key.split(".")[0]
    name = key.split(".")[1]
    if name == "fasta":
        from .reference import inspect_header

        target = absolute_path(value)
        inspect_header(target)
        value = str(target)
    elif name in {"index", "directory"}:
        value = str(absolute_path(value))
    elif value != "GRCh38":
        raise UserError("Only assembly GRCh38 is supported.")
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        # OS-backed locking releases on process exit; an empty lock file may remain.
        with FileLock(str(path) + ".lock", timeout=0):
            data = load(path, missing_ok=True)
            previous = resolved_reference(data, path)
            if name == "fasta" and previous.get("fasta") != value:
                # The previous index belongs to the previous reference.
                data.setdefault("reference", {}).pop("index", None)
            data.setdefault(section, {})[name] = value
            ref = resolved_reference(data, path)
            if ref.get("fasta") and ref.get("index") == ref["fasta"]:
                raise UserError("Index path must differ from the FASTA path.")
            validate(data)
            payload = tomli_w.dumps(data)
            validate(tomllib.loads(payload))
            temp: Path | None = None
            try:
                with tempfile.NamedTemporaryFile(
                    mode="w", encoding="utf-8", newline="\n", dir=path.parent,
                    prefix="." + path.name + ".", suffix=".tmp", delete=False,
                ) as handle:
                    temp = Path(handle.name)
                    handle.write(payload)
                    handle.flush()
                    os.fsync(handle.fileno())
                os.replace(temp, path)
            finally:
                if temp is not None:
                    temp.unlink(missing_ok=True)
    except Timeout as exc:
        raise UserError("Config is being changed by another process; retry later.") from exc


def select_reference(
    path: Path, *, explicit_config: bool, fasta: str | None, index: str | None,
) -> tuple[Path, Path, dict]:
    data = load(path, missing_ok=not explicit_config)
    ref = resolved_reference(data, path)
    chosen = absolute_path(fasta) if fasta is not None else (
        Path(ref["fasta"]) if "fasta" in ref else None
    )
    if chosen is None:
        raise UserError("Set reference.fasta or supply --reference-fasta PATH.")
    # A one-off reference must not inherit the index for a different reference.
    if index is not None:
        index_path = absolute_path(index)
        index_source = "argument"
    elif fasta is None and "index" in ref:
        index_path = Path(ref["index"])
        index_source = "config"
    else:
        index_path = Path(str(chosen) + ".fai")
        index_source = "adjacent"
    return chosen, index_path, {
        "fasta": "argument" if fasta is not None else "config",
        "index": index_source,
    }

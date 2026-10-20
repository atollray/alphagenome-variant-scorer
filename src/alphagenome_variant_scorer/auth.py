"""Credentials come only from the selected private TOML config."""

from pathlib import Path
from .config import load
from .errors import UserError


def resolve_key(*, config_file: Path, explicit_config=False):
    data = load(config_file, missing_ok=not explicit_config)
    value = data.get("auth", {}).get("api_key", "").strip()
    if not value or any(c.isspace() for c in value):
        raise UserError("Set auth.api_key in your private config file.")
    return value

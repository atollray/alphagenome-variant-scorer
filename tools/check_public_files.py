"""Audit Git-visible file paths and obvious secrets without printing their values."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import subprocess
import sys

ROOT_FILES = {".gitignore", ".gitattributes", "pyproject.toml", "README.md", "LICENSE", "MANIFEST.in"}
EXAMPLES = {"tests/fixtures/synthetic.fa", "examples/config.example.toml"}
PATTERNS = {
    "GitHub credential": re.compile(r"\b(?:gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,})"),
    "Google API key": re.compile(r"\bAIza[0-9A-Za-z_-]{30,}"),
    "private key": re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----"),
    "personal Windows path": re.compile(r"[A-Za-z]:[/\\]Users[/\\][^/\\\s]+"),
}


def allowed(name: str) -> bool:
    path = Path(name)
    if name in ROOT_FILES or name in EXAMPLES:
        return True
    if ".private" in name.lower():
        return False
    if path.suffix == ".py" and path.parts[0] in {"src", "tests", "tools"}:
        return True
    return False


def inspect_content(name: str, content: bytes) -> list[str]:
    issues = []
    if len(content) > 1_000_000:
        return ["file exceeds 1 MB review limit"]
    try:
        text = content.decode("utf-8")
    except UnicodeError:
        return ["unexpected binary file"]
    for label, pattern in PATTERNS.items():
        if pattern.search(text):
            issues.append(label)
    return issues


def audit(root: Path, *, staged: bool = False) -> tuple[list[str], list[str]]:
    command = ["git", "ls-files", "-z", "--cached"]
    if not staged:
        command += ["--others", "--exclude-standard"]
    names = sorted(set(subprocess.check_output(command, cwd=root).decode("utf-8").split("\0")) - {""})
    issues = []
    for name in names:
        if not allowed(name):
            issues.append(f"{name}: not on the public source allowlist")
        path = root / name
        if staged:
            entry = subprocess.check_output(["git", "ls-files", "-s", "--", name], cwd=root)
            if entry.startswith((b"120000", b"160000")):
                issues.append(f"{name}: symlink or submodule is not allowed")
                continue
            content = subprocess.check_output(["git", "show", ":" + name], cwd=root)
        else:
            if path.is_symlink() or not path.is_file():
                issues.append(f"{name}: not a regular file")
                continue
            content = path.read_bytes()
        issues.extend(f"{name}: {issue}" for issue in inspect_content(name, content))
    return names, issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--staged", action="store_true", help="Inspect the actual Git index contents.")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    names, issues = audit(root, staged=args.staged)
    for issue in issues:
        print(issue, file=sys.stderr)
    if not names:
        print("No files to audit.", file=sys.stderr)
        return 1
    print(f"Checked {len(names)} public candidate files; findings: {len(issues)}.")
    return 1 if issues else 0


if __name__ == "__main__":
    raise SystemExit(main())

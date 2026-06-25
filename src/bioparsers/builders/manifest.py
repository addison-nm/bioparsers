"""Reproducibility manifests for builder outputs.

A *build manifest* is a small JSON sidecar that records **how** a curated
dataset was produced: the bioparsers version and git state, the builder's
``name`` and ``description``, the environment, and optional run details
(output path, record count, a user-supplied description). It accompanies the
output file so a dataset can be traced back to the exact builder and code
that made it.

Modeled on the BioM3 ``build_manifest.json`` convention. Use it from a
recipe or your own script::

    n = write_jsonl(builder.build(records), "out.jsonl")
    write_manifest(builder, "out.jsonl.manifest.json",
                   description="SH3 domain subset", output="out.jsonl",
                   record_count=n)
"""

from __future__ import annotations

import json
import platform
import re
import subprocess
import sys
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version


def _bioparsers_version() -> str:
    try:
        return version("bioparsers")
    except PackageNotFoundError:
        return "unknown"


def _git(*args: str) -> str | None:
    """Run a git command, returning stripped stdout or None on any failure."""
    try:
        result = subprocess.run(
            ["git", *args], capture_output=True, text=True, timeout=5,
        )
    except Exception:
        return None
    return result.stdout.strip() if result.returncode == 0 else None


def _git_info() -> dict:
    """Best-effort git provenance; each value is None if git is unavailable."""
    remote = _git("remote", "get-url", "origin")
    if remote:
        # Strip embedded credentials: https://user:token@host/... -> https://host/...
        remote = re.sub(r"(https?://)[^@/]+@", r"\1", remote)
    porcelain = _git("status", "--porcelain")
    return {
        "hash": _git("rev-parse", "--short", "HEAD"),
        "dirty": (porcelain != "") if porcelain is not None else None,
        "branch": _git("rev-parse", "--abbrev-ref", "HEAD"),
        "remote": remote,
    }


def _build_string(ver: str, git: dict) -> str:
    """PEP 440-style local version, e.g. ``'0.0.1+gb310e6f.dirty'``."""
    out = ver
    if git.get("hash"):
        out += f"+g{git['hash']}"
        if git.get("dirty"):
            out += ".dirty"
    elif git.get("dirty"):
        out += "+dirty"
    return out


def generate_manifest(
    builder,
    *,
    description: str | None = None,
    output: str | None = None,
    record_count: int | None = None,
    extra: dict | None = None,
) -> dict:
    """Build a reproducibility manifest dict for a dataset made by *builder*.

    Captures auto-collected provenance — bioparsers version + git state,
    Python/platform, timestamp, and the invoking command — alongside the
    builder's ``name`` and ``description``. *description* is an optional
    free-text note from the caller; *output* and *record_count* describe the
    produced file; *extra* is merged in at the top level for any extra
    run-specific keys (e.g. the recipe's Pfam IDs).
    """
    ver = _bioparsers_version()
    git = _git_info()
    manifest = {
        "bioparsers_version": ver,
        "bioparsers_build": _build_string(ver, git),
        "git": git,
        "builder": {"name": builder.name, "description": builder.description},
        "description": description,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "command": " ".join(sys.argv),
        "environment": {"python": sys.version, "platform": platform.platform()},
        "output": output,
        "record_count": record_count,
    }
    if extra:
        manifest.update(extra)
    return manifest


def write_manifest(builder, path: str, *, indent: int = 2, **kwargs) -> str:
    """Write :func:`generate_manifest` for *builder* to *path* as JSON.

    Extra keyword args (``description``, ``output``, ``record_count``,
    ``extra``) are forwarded to :func:`generate_manifest`. Returns *path*.
    """
    manifest = generate_manifest(builder, **kwargs)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=indent, default=str)
    return path

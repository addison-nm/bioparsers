"""Shared orchestration for the `build_uniprot_by_pfam_*` recipes.

Filters a parsed UniProt JSONL stream to entries carrying one or more
specified Pfam domains, runs them through a builder, and writes JSONL —
either a single unioned file (``--join``) or one file per Pfam ID. Each
output file gets a ``<output>.manifest.json`` reproducibility sidecar.

The per-ID mode makes a single streaming pass over the (large) input,
routing each record to every matching ID's output file, so an entry that
carries several of the requested domains lands in each of their files.
"""

import os
import sys
from contextlib import ExitStack

from bioparsers.builders import jsonl_writer, load_jsonl, write_jsonl, write_manifest
from bioparsers.builders.uniprot import helpers


def with_pfam_suffix(path: str, pfam_id: str) -> str:
    """Insert *pfam_id* before the extension of *path*:
    ``data/sprot_flat.jsonl`` + ``PF00069`` -> ``data/sprot_flat.PF00069.jsonl``
    (``.jsonl.gz`` is preserved as a unit).
    """
    head, base = os.path.split(path)
    stem, dot, rest = base.partition(".")
    name = f"{stem}.{pfam_id}{dot}{rest}" if dot else f"{stem}.{pfam_id}"
    return os.path.join(head, name)


def _ensure_parent(path: str) -> None:
    parent = os.path.dirname(os.path.abspath(path))
    os.makedirs(parent, exist_ok=True)


def _write_sidecar(builder, out_path, *, count, description, extra):
    """Write a ``<out_path>.manifest.json`` reproducibility sidecar."""
    mpath = write_manifest(builder, out_path + ".manifest.json",
                           description=description, output=out_path,
                           record_count=count, extra=extra)
    print(f"{count} records -> {out_path}  (manifest: {mpath})", file=sys.stderr)


def run_by_pfam(builder, input_path, pfam_ids, output, *,
                join=False, gzip=False, description=None):
    """Build a Pfam-filtered dataset from *input_path* with *builder*.

    *join* True  -> one *output* file: the union of entries matching ANY ID.
    *join* False -> one file per ID, named via :func:`with_pfam_suffix`.

    Each output file gets a ``<output>.manifest.json`` sidecar recording the
    builder, environment, Pfam IDs, record count, and *description*.
    """
    _ensure_parent(output)

    if join:
        targets = set(pfam_ids)
        selected = (r for r in load_jsonl(input_path)
                    if not targets.isdisjoint(helpers.pfam_ids(r)))
        n = write_jsonl(builder.build(selected), output, gzip=gzip)
        _write_sidecar(builder, output, count=n, description=description,
                       extra={"pfam_ids": list(pfam_ids), "join": True})
        return

    paths = {pid: with_pfam_suffix(output, pid) for pid in pfam_ids}
    counts = {pid: 0 for pid in pfam_ids}
    with ExitStack() as stack:
        writers = {pid: stack.enter_context(jsonl_writer(paths[pid], gzip=gzip))
                   for pid in pfam_ids}
        for rec in load_jsonl(input_path):
            rec_pfams = set(helpers.pfam_ids(rec))
            matched = [pid for pid in pfam_ids if pid in rec_pfams]
            if not matched:
                continue
            # builder.build may drop the record (its own filters) -> 0 or 1 out.
            for out in builder.build([rec]):
                for pid in matched:
                    writers[pid](out)
                    counts[pid] += 1
    for pid in pfam_ids:
        _write_sidecar(builder, paths[pid], count=counts[pid],
                       description=description, extra={"pfam_ids": [pid], "join": False})

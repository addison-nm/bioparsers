"""Pfam-partitioned builder runner with a UniProt annotation join.

The Pfam counterpart to ``bioparsers.builders.uniprot.run_by_pfam``. Where the
UniProt runner makes a single streaming pass over one UniProt file, a Pfam
dataset is a **join of three sources**:

1. **members** — a parsed Pfam member FASTA JSONL (``bioparsers pfam-fasta``):
   one redundancy-reduced domain sequence per family member, carrying the
   member's UniProt ``accession``.
2. **families** — ``{pfam_accession: {name, description}}`` family metadata
   (loaded by the recipe from the small ``bioparsers pfam`` family JSONL).
3. **uniprot** — one or more parsed UniProt JSONL files (Swiss-Prot, TrEMBL),
   from which annotation (protein name, function, GO, lineage, ...) is pulled
   for each member by its accession.

The join strategy (chosen for the SH3 proof-of-concept, ~24k target
accessions) is an **in-memory accession index**: the runner collects the
member accessions for the requested families, then streams the UniProt
file(s) keeping only matching records in a ``{accession: record}`` dict
(stopping early once every target is found). Each member is then composed into
a single record — ``{...member, family, uniprot}`` — and handed to the
:class:`~bioparsers.builders.base.Builder`. A member whose accession is empty
or absent from UniProt gets ``uniprot=None`` (family-only annotation).

The UniProt scan is the dominant cost (a full TrEMBL pass is hundreds of GB),
so three optimizations apply:

- **member early-stop** — the member FASTA is grouped by family (contiguous
  blocks, in ``Pfam-A.full`` order), so reading stops once the requested
  families' blocks have passed rather than scanning the whole file.
- **accession prefilter** — each UniProt line is cheaply tested for a target
  accession (a substring/field check) *before* ``json.loads``, so the ~all
  non-matching lines skip JSON parsing. Decompression uses ``pigz`` when
  available (see :func:`bioparsers.builders.io.iter_text_lines`).
- **annotation cache** (opt-in, ``cache_path``) — the resolved
  ``{accession: record}`` subset is written to disk with the requested
  accession set + UniProt sources it was built from. A later run whose
  requested accessions are a subset (same sources) reuses it and skips the
  UniProt scan entirely.

Output mirrors ``run_by_pfam``: one file per Pfam ID (``with_pfam_suffix``) or
a single unioned file (``join=True``), each with a ``<output>.manifest.json``
sidecar.
"""

from __future__ import annotations

import json
import os
import sys

from bioparsers.builders import (
    iter_text_lines,
    load_jsonl,
    write_jsonl,
    write_manifest,
)
from bioparsers.builders.uniprot.runner import with_pfam_suffix, _ensure_parent


# ===========================================================================
# Member read (early-stop) & UniProt index (prefilter + cache)
# ===========================================================================

def _line_pfam(line: str):
    """Cheaply read the ``pfam_accession`` value from a member JSONL *line*
    without a full ``json.loads`` — the prefilter for the member scan."""
    i = line.find('"pfam_accession"')
    if i == -1:
        return None
    q1 = line.find('"', line.find(":", i) + 1)
    q2 = line.find('"', q1 + 1) if q1 != -1 else -1
    return line[q1 + 1:q2] if q1 != -1 and q2 != -1 else None


def _iter_target_members(members_path, targets):
    """Yield members of the requested *targets* families, stopping early.

    The member FASTA (and the JSONL parsed from it) is grouped by family in
    ``Pfam-A.full`` order — but that order is *not* accession-sorted, so a
    family can sit deep in the (multi-GB) file. Each line's ``pfam_accession``
    is read with a cheap substring check (decompressing via ``pigz`` when
    available) and only target-family lines are ``json.loads``-ed. Once every
    requested family's contiguous block has passed, the scan returns without
    reading the rest of the file. A requested family that never appears simply
    means the file is read to the end (nothing to stop on).
    """
    remaining = set(targets)
    last = None
    for line in iter_text_lines(members_path):
        if not line.strip():
            continue
        pfam = _line_pfam(line)
        if pfam != last and last in remaining:
            remaining.discard(last)
            if not remaining:
                return  # left the last requested family's block
        last = pfam
        if pfam in targets:
            yield json.loads(line)


def _candidate_accessions(line: str):
    """Cheaply pull the accession tokens out of a parsed-UniProt JSONL *line*
    without a full ``json.loads`` — the prefilter for the big scan.

    Reads the ``primary_accession`` value and the ``accessions`` array via
    ``str.find`` (locating each field's exact span, so a long accession list is
    never truncated). Returns the tokens; the caller tests them against the
    target set and only parses the line on a hit.
    """
    cands = []
    i = line.find('"primary_accession"')
    if i != -1:
        q1 = line.find('"', line.find(":", i) + 1)
        q2 = line.find('"', q1 + 1) if q1 != -1 else -1
        if q1 != -1 and q2 != -1:
            cands.append(line[q1 + 1:q2])
    i = line.find('"accessions"')
    if i != -1:
        lb = line.find("[", i)
        rb = line.find("]", lb) if lb != -1 else -1
        if lb != -1 and rb != -1:
            for tok in line[lb + 1:rb].split(","):
                tok = tok.strip().strip('"')
                if tok:
                    cands.append(tok)
    return cands


def _build_uniprot_index(uniprot_paths, accession_set):
    """Stream *uniprot_paths*, returning ``{accession: record}`` for every
    record matching any accession in *accession_set*.

    Each line is prefiltered on its accession tokens before ``json.loads`` so
    non-matching records (the overwhelming majority) skip JSON parsing. A
    matched UniProt entry is indexed under *all* of its accessions (primary +
    secondary), so a member referring to an old/secondary accession resolves.
    Scanning a file stops once every requested accession is found; later files
    are skipped entirely once nothing remains to look up.
    """
    index: dict = {}
    remaining = set(accession_set)
    if not remaining:
        return index
    for path in uniprot_paths:
        if not remaining:
            break
        print(f"[pfam-join] scanning {path} for {len(remaining)} accessions...",
              file=sys.stderr)
        for line in iter_text_lines(path):
            if not line.strip():
                continue
            hit = remaining.intersection(_candidate_accessions(line))
            if not hit:
                continue
            rec = json.loads(line)
            accs = set(rec.get("accessions") or [])
            primary = rec.get("primary_accession")
            if primary:
                accs.add(primary)
            for acc in accs:
                index[acc] = rec
            remaining -= hit
            if not remaining:
                break
    print(f"[pfam-join] resolved {len(accession_set) - len(remaining)}/"
          f"{len(accession_set)} accessions", file=sys.stderr)
    return index


def _load_uniprot_cache(cache_path, requested, uniprot_paths):
    """Return a ``{accession: record}`` index from a prior cache, or ``None``.

    Reuses the cache only when it was built from the *same* UniProt sources and
    its requested accession set is a **superset** of *requested* — so a miss in
    the cache is authoritative (that accession was already searched and not
    found) and no rescan is needed.
    """
    meta_path = cache_path + ".meta.json"
    if not (os.path.exists(cache_path) and os.path.exists(meta_path)):
        return None
    with open(meta_path) as f:
        meta = json.load(f)
    if meta.get("uniprot_sources") != list(uniprot_paths):
        return None
    if not requested <= set(meta.get("requested_accessions") or []):
        return None
    index: dict = {}
    for rec in load_jsonl(cache_path):
        accs = set(rec.get("accessions") or [])
        primary = rec.get("primary_accession")
        if primary:
            accs.add(primary)
        for acc in accs:
            index[acc] = rec
    print(f"[pfam-join] reusing UniProt cache {cache_path} "
          f"({len(index)} accessions)", file=sys.stderr)
    return index


def _save_uniprot_cache(cache_path, requested, uniprot_paths, index):
    """Write the resolved subset + a sidecar recording the requested accession
    set and UniProt sources, for reuse by :func:`_load_uniprot_cache`. The data
    file is gzip-compressed when *cache_path* ends with ``.gz``."""
    _ensure_parent(cache_path)
    seen = set()

    def _unique_records():
        for rec in index.values():
            key = rec.get("primary_accession") or id(rec)
            if key in seen:
                continue
            seen.add(key)
            yield rec

    n = write_jsonl(_unique_records(), cache_path, gzip=cache_path.endswith(".gz"))
    with open(cache_path + ".meta.json", "w") as f:
        json.dump({"uniprot_sources": list(uniprot_paths),
                   "requested_accessions": sorted(requested),
                   "resolved_records": n}, f)
    print(f"[pfam-join] wrote UniProt cache {cache_path} ({n} records)",
          file=sys.stderr)


def _is_dir_cache(cache_path: str) -> bool:
    """A *folder* cache (one file per Pfam ID) vs a single-file cache: true when
    the path ends with a path separator or names an existing directory."""
    return (cache_path.endswith(os.sep) or cache_path.endswith("/")
            or os.path.isdir(cache_path))


def _resolve_index(uniprot_paths, accs_by_pfam, cache_path):
    """Resolve ``{accession: record}`` for every requested family.

    No *cache_path*: one scan over the union of all families' accessions. A
    single-file *cache_path*: cache that union. A **directory** *cache_path*:
    keep one cache file per Pfam ID (``<dir>/<PF>.jsonl.gz``) — already-cached
    families load instantly and only the rest share a single UniProt scan, so
    each family's resolved subset is independently reusable across runs and
    recipes.
    """
    all_accs = set().union(*accs_by_pfam.values()) if accs_by_pfam else set()
    if not cache_path:
        return _build_uniprot_index(uniprot_paths, all_accs)
    if _is_dir_cache(cache_path):
        return _resolve_with_dir_cache(uniprot_paths, accs_by_pfam, cache_path)
    return _uniprot_index(uniprot_paths, all_accs, cache_path)


def _resolve_with_dir_cache(uniprot_paths, accs_by_pfam, cache_dir):
    """Per-family folder cache: load what is cached, scan the rest once, and
    write each newly-resolved family's subset to its own ``<PF>.jsonl.gz``."""
    os.makedirs(cache_dir, exist_ok=True)
    index: dict = {}
    to_scan, scan_accs = [], set()
    for pid, accs in accs_by_pfam.items():
        fam = _load_uniprot_cache(os.path.join(cache_dir, f"{pid}.jsonl.gz"),
                                  accs, uniprot_paths)
        if fam is not None:
            index.update(fam)
        else:
            to_scan.append(pid)
            scan_accs |= accs
    if scan_accs:
        scanned = _build_uniprot_index(uniprot_paths, scan_accs)
        index.update(scanned)
        for pid in to_scan:
            sub: dict = {}
            for acc in accs_by_pfam[pid]:
                rec = scanned.get(acc)
                if rec is None:
                    continue
                for key in (rec.get("accessions") or []):
                    sub[key] = rec
                if rec.get("primary_accession"):
                    sub[rec["primary_accession"]] = rec
            _save_uniprot_cache(os.path.join(cache_dir, f"{pid}.jsonl.gz"),
                                accs_by_pfam[pid], uniprot_paths, sub)
    return index


def _uniprot_index(uniprot_paths, accession_set, cache_path):
    """Resolve the accession index from a single-file *cache_path* (or build
    and write it). Used for the union scan and folder-mode's per-family files."""
    if cache_path:
        index = _load_uniprot_cache(cache_path, accession_set, uniprot_paths)
        if index is not None:
            return index
        index = _build_uniprot_index(uniprot_paths, accession_set)
        _save_uniprot_cache(cache_path, accession_set, uniprot_paths, index)
        return index
    return _build_uniprot_index(uniprot_paths, accession_set)


# ===========================================================================
# Runner
# ===========================================================================

def _composite(member: dict, families: dict, uniprot_index: dict) -> dict:
    """Join one member with its family metadata and UniProt record."""
    return {
        **member,
        "family": families.get(member.get("pfam_accession")) or {},
        "uniprot": uniprot_index.get(member.get("accession")),
    }


def _write_sidecar(builder, out_path, *, count, description, extra):
    mpath = write_manifest(builder, out_path + ".manifest.json",
                           description=description, output=out_path,
                           record_count=count, extra=extra)
    print(f"{count} records -> {out_path}  (manifest: {mpath})", file=sys.stderr)


def run_pfam_join(builder, members_path, pfam_filter, *, families,
                  uniprot_paths, output, join=False, gzip=False,
                  description=None, cache_path=None):
    """Build a Pfam dataset from *members_path*, joining family + UniProt data.

    *pfam_filter* is the iterable of Pfam accessions to select. *families* is
    the ``{pfam_accession: {name, description}}`` metadata map; *uniprot_paths*
    is the list of parsed UniProt JSONL files to draw annotation from.
    *cache_path*, when given, caches the resolved UniProt subset for reuse
    across runs: a single-file path caches the union of all requested families'
    accessions; a path that is a **directory** (or ends with a path separator)
    caches one file per Pfam ID, so each family is independently reusable (see
    :func:`_resolve_index`).

    *join* True  -> one *output* file: every selected member, in one file.
    *join* False -> one file per Pfam ID, named via :func:`with_pfam_suffix`.

    Each output gets a ``<output>.manifest.json`` sidecar recording the
    builder, environment, Pfam IDs, UniProt sources, record count, and
    *description*.
    """
    pfam_list = list(pfam_filter)
    targets = set(pfam_list)
    _ensure_parent(output)

    # Materialize the selected members (one short domain record each) so the
    # member stream can be walked twice: once to gather accessions, once to emit.
    members = list(_iter_target_members(members_path, targets))
    accs_by_pfam: dict = {pid: set() for pid in pfam_list}
    for m in members:
        acc = m.get("accession")
        if acc:
            accs_by_pfam.setdefault(m.get("pfam_accession"), set()).add(acc)
    uniprot_index = _resolve_index(uniprot_paths, accs_by_pfam, cache_path)

    extra_base = {"pfam_ids": pfam_list, "uniprot_sources": list(uniprot_paths)}

    if join:
        comps = (_composite(m, families, uniprot_index) for m in members)
        n = write_jsonl(builder.build(comps), output, gzip=gzip)
        _write_sidecar(builder, output, count=n, description=description,
                       extra={**extra_base, "join": True})
        return

    for pid in pfam_list:
        path = with_pfam_suffix(output, pid)
        comps = (_composite(m, families, uniprot_index)
                 for m in members if m.get("pfam_accession") == pid)
        n = write_jsonl(builder.build(comps), path, gzip=gzip)
        _write_sidecar(builder, path, count=n, description=description,
                       extra={**extra_base, "pfam_ids": [pid], "join": False})

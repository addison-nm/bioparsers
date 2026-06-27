#!/usr/bin/env python
"""Build Recipe: Concatenate labeled JSONL datasets into one source-tagged JSONL.

A generic utility recipe (no hardcoded sources): given any number of
``NAME=PATH`` arguments, concatenate their records in the given order, adding the
source ``NAME`` at the root of each record (key ``source`` by default). The
inputs are any JSONL[.gz] files — e.g. combine the three SH3 ``caption_fields``
datasets so a trainer gets one stream that knows where each entry came from:

    python recipes/concatenate_datasets.py \\
        supplemental=outputs/supplement_caption_fields_SH3.jsonl \\
        swissprot=outputs/swissprot_caption_fields_SH3.PF00018.jsonl \\
        pfam=outputs/pfam_caption_fields_SH3.PF00018.jsonl \\
        -o outputs/SH3_caption_fields_all.jsonl

Each output record is the input record with ``source`` added at the root, e.g.
``{source, accession, sequence, region?, pfam_ids, fields, caption_fields}``.
The core logic is :func:`bioparsers.builders.concatenate`.
"""

import argparse
import os
import sys

from bioparsers.builders import (
    ConcatenatedDataset,
    concatenate,
    jsonl_writer,
    write_manifest,
)


def _pair(arg: str):
    """Parse a ``NAME=PATH`` source argument into ``(name, path)``."""
    name, sep, path = arg.partition("=")
    if not sep or not name or not path:
        raise argparse.ArgumentTypeError(f"expected NAME=PATH, got {arg!r}")
    return (name, path)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=ConcatenatedDataset.description)
    p.add_argument("sources", nargs="+", type=_pair, metavar="NAME=PATH",
                   help="labeled JSONL[.gz] source(s) to concatenate, in order")
    p.add_argument("-o", "--output", required=True, help="output JSONL path")
    p.add_argument("--gzip", action="store_true", help="gzip the output")
    p.add_argument("--source-key", default="source",
                   help="root key holding the source name (default: source)")
    p.add_argument("--description", default=None,
                   help="free-text note recorded in the build manifest")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    counts: dict = {}
    with jsonl_writer(args.output, gzip=args.gzip) as write:
        for name, path in args.sources:
            n = 0
            for rec in concatenate([(name, path)], source_key=args.source_key):
                write(rec)
                n += 1
            counts[name] = counts.get(name, 0) + n

    total = sum(counts.values())
    mpath = write_manifest(
        ConcatenatedDataset(), args.output + ".manifest.json",
        description=args.description, output=args.output, record_count=total,
        extra={"source_key": args.source_key, "counts": counts,
               "sources": [{"name": n, "path": p} for n, p in args.sources]},
    )
    summary = ", ".join(f"{n}={counts[n]}" for n, _ in args.sources)
    print(f"{total} records -> {args.output}  ({summary})  (manifest: {mpath})",
          file=sys.stderr)


if __name__ == "__main__":
    main()

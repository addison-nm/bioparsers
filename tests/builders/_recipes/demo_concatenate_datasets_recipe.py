#!/usr/bin/env python
"""Test-fixture recipe — a self-contained copy of ``concatenate_datasets.py``."""

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
    name, sep, path = arg.partition("=")
    if not sep or not name or not path:
        raise argparse.ArgumentTypeError(f"expected NAME=PATH, got {arg!r}")
    return (name, path)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=ConcatenatedDataset.description)
    p.add_argument("sources", nargs="+", type=_pair, metavar="NAME=PATH")
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--gzip", action="store_true")
    p.add_argument("--source-key", default="source")
    p.add_argument("--description", default=None)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    counts = {}
    with jsonl_writer(args.output, gzip=args.gzip) as write:
        for name, path in args.sources:
            n = 0
            for rec in concatenate([(name, path)], source_key=args.source_key):
                write(rec)
                n += 1
            counts[name] = counts.get(name, 0) + n
    write_manifest(
        ConcatenatedDataset(), args.output + ".manifest.json",
        description=args.description, output=args.output,
        record_count=sum(counts.values()),
        extra={"source_key": args.source_key, "counts": counts,
               "sources": [{"name": n, "path": p} for n, p in args.sources]},
    )


if __name__ == "__main__":
    main()

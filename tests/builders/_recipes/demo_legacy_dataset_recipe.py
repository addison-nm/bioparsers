#!/usr/bin/env python
"""Test-fixture recipe — a self-contained copy of ``build_legacy_dataset.py``."""

import argparse
import csv
import gzip
import os
import sys
from typing import Iterable, Iterator

from bioparsers.builders import Builder, load_jsonl, write_manifest

COLUMNS = ["primary_Accession", "protein_sequence", "[final]text_caption", "pfam_label"]
SECTION_ORDER = ["supplement", "swissprot", "pfam"]


class LegacyDatasetBuilder(Builder):
    """Maps a section's legacy records to the legacy dataset CSV columns
    (test-fixture copy)."""

    name = "legacy_sh3_dataset"

    def build(self, records: Iterable[dict]) -> Iterator[dict]:
        for rec in records:
            yield {
                "primary_Accession": rec.get("accession") or "",
                "protein_sequence": rec.get("sequence") or "",
                "[final]text_caption": rec.get("caption") or "",
                "pfam_label": _pfam_label(rec.get("_section"), rec.get("pfam_ids") or []),
            }


def _pfam_label(section: str, pfam_ids: list) -> str:
    if section == "swissprot":
        return repr(list(pfam_ids))
    if section == "pfam":
        return pfam_ids[0] if pfam_ids else ""
    return ""


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=LegacyDatasetBuilder.description)
    p.add_argument("--supplement", default=None)
    p.add_argument("--swissprot", default=None)
    p.add_argument("--pfam", default=None)
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--gzip", action="store_true")
    p.add_argument("--description", default=None)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    sources = {"supplement": args.supplement, "swissprot": args.swissprot,
               "pfam": args.pfam}
    if not any(sources.values()):
        sys.exit("at least one of --supplement / --swissprot / --pfam is required")
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    builder = LegacyDatasetBuilder()
    opener = gzip.open if args.gzip else open
    counts = {}
    with opener(args.output, "wt", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS, lineterminator="\n")
        writer.writeheader()
        for section in SECTION_ORDER:
            path = sources[section]
            if not path:
                continue
            tagged = (dict(rec, _section=section) for rec in load_jsonl(path))
            n = 0
            for row in builder.build(tagged):
                writer.writerow(row)
                n += 1
            counts[section] = n
    write_manifest(builder, args.output + ".manifest.json",
                   description=args.description, output=args.output,
                   record_count=sum(counts.values()),
                   extra={"sections": counts, "section_order": SECTION_ORDER,
                          "inputs": {k: v for k, v in sources.items() if v}})


if __name__ == "__main__":
    main()

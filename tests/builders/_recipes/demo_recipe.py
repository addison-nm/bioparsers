#!/usr/bin/env python
"""Test-fixture recipe — a self-contained copy of a builder recipe.

Kept under ``tests/`` so the suite exercises a representative ``Builder`` and
the ``run_by_pfam`` CLI **without depending on the live scripts in
``recipes/``** (which are user-editable and get renamed/reorganized). Mirror
the recipe you actually want to cover here; this copy may drift from it.
"""

import argparse
from typing import Iterable, Iterator

from bioparsers.builders import Builder
from bioparsers.builders.uniprot import filters, helpers, run_by_pfam


class SwissProtDemoFields(Builder):
    """Curated records pairing a sequence with a nested block of text fields.

    Output record::

        {
          "accession": str,    # primary accession (always present)
          "sequence":  str,    # amino-acid sequence (always present)
          "fields": {          # always present; contains only populated keys
            "name":     str,   # RecName/SubName full name
            "function": str,   # FUNCTION comment, cleaned
            "domains":  str     # DOMAIN comment, cleaned
          }
        }

    Empty fields are omitted, so ``fields`` may be empty for a sparse entry.

    Options:
      - ``reviewed_only``: keep only Swiss-Prot (Reviewed) entries.
      - ``min_length``: drop entries shorter than this many residues.
    """

    name = "swissprot_demo_fields"

    def __init__(self, *, reviewed_only: bool = False, min_length: int = 0):
        self.reviewed_only = reviewed_only
        self.min_length = min_length

    def build(self, records: Iterable[dict]) -> Iterator[dict]:
        keep_len = filters.min_length(self.min_length)
        for rec in records:
            if self.reviewed_only and not filters.is_reviewed(rec):
                continue
            if not keep_len(rec):
                continue

            fields = {}
            name = helpers.full_name(rec)
            if name:
                fields["name"] = name
            function = helpers.joined_comment(rec, "FUNCTION")
            if function:
                fields["function"] = function
            domains = helpers.joined_comment(rec, "DOMAIN")
            if domains:
                fields["domains"] = domains

            yield {
                "accession": rec.get("primary_accession"),
                "sequence": rec.get("sequence"),
                "fields": fields,
            }


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=SwissProtDemoFields.description)
    p.add_argument("input", help="parsed Swiss-Prot JSONL (plain or .gz)")
    p.add_argument("--pfam-ids", nargs="+", required=True, metavar="PFAM_ID",
                   dest="pfam_ids", help="one or more Pfam accessions")
    p.add_argument("-o", "--output", required=True,
                   help="output JSONL path; per-ID mode inserts the Pfam ID")
    p.add_argument("--join", action="store_true",
                   help="write a single union file instead of one per Pfam ID")
    p.add_argument("--gzip", action="store_true", help="gzip the output")
    p.add_argument("--description", default=None,
                   help="free-text note recorded in each output's build manifest")
    p.add_argument("--reviewed-only", action="store_true")
    p.add_argument("--min-length", type=int, default=0)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    builder = SwissProtDemoFields(reviewed_only=args.reviewed_only,
                                  min_length=args.min_length)
    run_by_pfam(builder, args.input, args.pfam_ids, args.output,
                join=args.join, gzip=args.gzip, description=args.description)


if __name__ == "__main__":
    main()

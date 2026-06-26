#!/usr/bin/env python
"""Recipe: build ``uniprot_flat_demo`` records, filtered by Pfam domain.

A worked example of the `bioparsers.builders` framework: define a custom
``Builder`` (below), then use it with the framework's filters, helpers, and
streaming I/O to compose a dataset from parsed UniProt JSONL.

Keeps only entries carrying one or more specified Pfam domains, then runs
them through ``UniprotFlat``. By default writes one output file per Pfam
ID; with ``--join`` writes a single file containing the union of entries
matching any of the IDs (no duplication).

    # one file per Pfam ID -> sprot_flat.PF00069.jsonl, sprot_flat.PF00027.jsonl
    python recipes/build_uniprot_by_pfam_flat_demo.py outputs/uniprot_sprot.jsonl \\
        --pfam-ids PF00069 PF00027 -o outputs/sprot_flat.jsonl

    # union of both into one file
    python recipes/build_uniprot_by_pfam_flat_demo.py outputs/uniprot_sprot.jsonl \\
        --pfam-ids PF00069 PF00027 --join -o outputs/sprot_kinases.jsonl
"""

import argparse
from typing import Iterable, Iterator

from bioparsers.builders import Builder
from bioparsers.builders.uniprot import filters, helpers, run_by_pfam



class UniprotFlat(Builder):
    """Flat curated records, one per UniProt entry.

    Output record::

        {
          "accession":  str,   # primary accession (always present)
          "entry_name": str,   # ID-line mnemonic (always present)
          "length":     int,   # sequence length  (always present)
          "sequence":   str,   # amino-acid sequence (always present)
          "name":       str,   # RecName/SubName full name (omitted if absent)
          "function":   str,   # FUNCTION comment, cleaned (omitted if absent)
          "keywords":   [str]  # only when include_keywords=True and non-empty
        }

    ``name`` is the protein's full recommended name (falling back to the
    TrEMBL SubName); ``function`` is the entry's FUNCTION comment(s) with
    ``{ECO:...}`` evidence and ``(PubMed:...)`` citations removed. Optional
    text fields are omitted when the source has no value (variable schema),
    so every present key carries real data.

    Options:
      - ``reviewed_only``: keep only Swiss-Prot (Reviewed) entries.
      - ``min_length``: drop entries shorter than this many residues.
      - ``require_function``: skip entries with no FUNCTION comment.
      - ``include_keywords``: add an evidence-stripped ``keywords`` list.
    """

    name = "uniprot_flat_demo"

    def __init__(
        self,
        *,
        reviewed_only: bool = False,
        min_length: int = 0,
        require_function: bool = True,
        include_keywords: bool = False,
    ):
        self.reviewed_only = reviewed_only
        self.min_length = min_length
        self.require_function = require_function
        self.include_keywords = include_keywords

    def build(self, records: Iterable[dict]) -> Iterator[dict]:
        keep_len = filters.min_length(self.min_length)
        for rec in records:
            if self.reviewed_only and not filters.is_reviewed(rec):
                continue
            if not keep_len(rec):
                continue
            function = helpers.joined_comment(rec, "FUNCTION")
            if self.require_function and not function:
                continue

            out = {
                "accession": rec.get("primary_accession"),
                "entry_name": rec.get("entry_name"),
                "length": rec.get("sequence_length"),
                "sequence": rec.get("sequence"),
            }
            name = helpers.full_name(rec)
            if name:
                out["name"] = name
            if function:
                out["function"] = function
            if self.include_keywords:
                kws = helpers.keywords(rec)
                if kws:
                    out["keywords"] = kws
            yield out


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=UniprotFlat.description)
    p.add_argument("input", help="parsed UniProt JSONL (plain or .gz)")
    p.add_argument("--pfam-ids", nargs="+", required=True, metavar="PFAM_ID",
                   dest="pfam_ids", help="one or more Pfam accessions, e.g. PF00069")
    p.add_argument("-o", "--output", required=True,
                   help="output JSONL path; in per-ID mode the Pfam ID is "
                        "inserted before the extension")
    p.add_argument("--join", action="store_true",
                   help="write a single union file instead of one file per Pfam ID")
    p.add_argument("--gzip", action="store_true", help="gzip the output")
    p.add_argument("--description", default=None,
                   help="free-text note recorded in each output's build manifest")
    p.add_argument("--reviewed-only", action="store_true")
    p.add_argument("--min-length", type=int, default=0)
    p.add_argument("--include-keywords", action="store_true")
    p.add_argument("--no-require-function", action="store_true",
                   help="keep entries that have no FUNCTION comment")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    builder = UniprotFlat(
        reviewed_only=args.reviewed_only,
        min_length=args.min_length,
        require_function=not args.no_require_function,
        include_keywords=args.include_keywords,
    )
    run_by_pfam(builder, args.input, args.pfam_ids, args.output,
                join=args.join, gzip=args.gzip, description=args.description)


if __name__ == "__main__":
    main()

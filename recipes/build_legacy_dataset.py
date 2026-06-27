#!/usr/bin/env python
"""Build Recipe: Assemble the complete legacy BioM3 SH3 dataset CSV.

Concatenates the three section outputs — Supplemental, Swiss-Prot, Pfam, **in
that order** — into a single 4-column CSV matching
``FINAL_SH3_all_dataset_with_prompts.csv`` (columns: primary_Accession,
protein_sequence, [final]text_caption, pfam_label).

Each section is produced independently by its own recipe (so each keeps its own
inputs / Pfam filter / UniProt cache), and this recipe just stitches their
legacy JSONL outputs together and maps them to the legacy CSV columns:

  - ``primary_Accession`` <- ``accession`` (empty string when absent)
  - ``protein_sequence``  <- ``sequence``
  - ``[final]text_caption`` <- ``caption``
  - ``pfam_label`` is **section-specific** (a legacy quirk):
      * supplement -> ``""`` (no Pfam label)
      * swissprot  -> the Python list-repr of the entry's Pfam IDs,
                      e.g. ``"['PF07714', 'PF00017', 'PF00018']"``
      * pfam       -> the single family accession, e.g. ``"PF00018"``

Fidelity follows the sections: the Supplemental section reproduces exactly, while
the Swiss-Prot and Pfam sections are approximate (Pfam-release drift), so the
assembled file has the same form and ordering as the legacy CSV but not an
identical row count.

Usage (after building each section):
    python recipes/build_legacy_dataset.py \\
        --supplement outputs/supplement_legacy_SH3.jsonl \\
        --swissprot  outputs/swissprot_legacy_SH3.PF00018.jsonl \\
        --pfam       outputs/pfam_legacy_SH3.PF00018.jsonl \\
        -o outputs/FINAL_SH3_all_reproduced.csv
"""

import argparse
import csv
import gzip
import os
import sys
from typing import Iterable, Iterator

from bioparsers.builders import Builder, load_jsonl, write_manifest

#: The legacy dataset's columns, in order.
COLUMNS = ["primary_Accession", "protein_sequence", "[final]text_caption", "pfam_label"]

#: Section name -> the recipe's source artifact, in legacy concatenation order.
SECTION_ORDER = ["supplement", "swissprot", "pfam"]


class LegacyDatasetBuilder(Builder):
    """Maps a section's legacy records to the legacy dataset's CSV columns.

    Input records are the legacy JSONL emitted by the section recipes
    (``{accession, sequence, caption, pfam_ids, ...}``), each tagged with its
    ``_section`` (``supplement`` / ``swissprot`` / ``pfam``). Output record::

        {primary_Accession, protein_sequence, [final]text_caption, pfam_label}

    where ``pfam_label`` follows the section-specific legacy convention (empty /
    list-repr / single accession).
    """

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
    """The legacy ``pfam_label`` for *section*: empty (supplement), the
    list-repr of all IDs (swissprot), or the single family accession (pfam)."""
    if section == "swissprot":
        return repr(list(pfam_ids))
    if section == "pfam":
        return pfam_ids[0] if pfam_ids else ""
    return ""  # supplement (and any unlabeled section)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=LegacyDatasetBuilder.description)
    p.add_argument("--supplement", default=None, metavar="JSONL",
                   help="supplement_legacy output JSONL")
    p.add_argument("--swissprot", default=None, metavar="JSONL",
                   help="swissprot_legacy output JSONL")
    p.add_argument("--pfam", default=None, metavar="JSONL",
                   help="pfam_legacy output JSONL")
    p.add_argument("-o", "--output", required=True, help="output CSV path")
    p.add_argument("--gzip", action="store_true", help="gzip the output CSV")
    p.add_argument("--description", default=None,
                   help="free-text note recorded in the build manifest")
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

    total = sum(counts.values())
    mpath = write_manifest(builder, args.output + ".manifest.json",
                           description=args.description, output=args.output,
                           record_count=total,
                           extra={"sections": counts, "section_order": SECTION_ORDER,
                                  "inputs": {k: v for k, v in sources.items() if v}})
    summary = ", ".join(f"{s}={counts.get(s, 0)}" for s in SECTION_ORDER)
    print(f"{total} rows -> {args.output}  ({summary})  (manifest: {mpath})",
          file=sys.stderr)


if __name__ == "__main__":
    main()

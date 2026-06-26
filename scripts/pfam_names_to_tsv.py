#!/usr/bin/env python
"""Helper script. Write the Pfam accession -> family-name table as a two-column 
TSV (PF#####<TAB>name), sorted by accession.

Usage: pfam_names_to_tsv.py INPUT OUTPUT
  INPUT   Pfam-A.hmm.gz (or any Pfam Stockholm/HMM file)
  OUTPUT  destination .tsv path
"""

import sys

from bioparsers.parsers.pfam_stockholm import family_name_map


def main(argv):
    if len(argv) != 2:
        sys.exit("usage: pfam_names_to_tsv.py INPUT OUTPUT")
    src, out = argv
    names = family_name_map(src)
    with open(out, "w") as fh:
        for accession, family in names.items():
            fh.write(f"{accession}\t{family}\n")
    print(f"{len(names)} families -> {out}")


if __name__ == "__main__":
    main(sys.argv[1:])

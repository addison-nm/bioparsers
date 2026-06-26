#!/usr/bin/env bash
# Parse the TrEMBL .dat release to gzipped JSONL. Run from the repo root.
# NOTE: the TrEMBL input is ~160 GB gzipped; plain JSONL output would be on
# the order of a terabyte, so this defaults to gzip-compressed output.
set -euo pipefail

mkdir -p outputs

bioparsers uniprot databases/trembl/uniprot_trembl.dat.gz \
    -o /media/nm-data/data/trembl_json/uniprot_trembl.jsonl.gz --progress 1000000 --gzip

# bioparsers uniprot databases/trembl/uniprot_trembl.dat.gz \
#     -o outputs/uniprot_trembl.jsonl --progress

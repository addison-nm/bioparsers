#!/usr/bin/env bash
# Parse the Swiss-Prot .dat release to gzipped JSONL. Run from the repo root.
set -euo pipefail

mkdir -p outputs

bioparsers uniprot data/swissprot/uniprot_sprot.dat.gz \
    -o outputs/uniprot_sprot.jsonl --progress
    
# bioparsers uniprot data/swissprot/uniprot_sprot.dat.gz \
#     -o outputs/uniprot_sprot.jsonl.gz --progress --gzip


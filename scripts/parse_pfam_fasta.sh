#!/usr/bin/env bash
#
# Parse the Pfam-A member FASTA (Pfam-A.fasta.gz) to gzipped JSONL, one record 
# per member sequence (the redundancy-reduced set).
#
# Usage: scripts/parse_pfam_fasta.sh [OUTPUT] [--link]
#   With no OUTPUT, writes data/pfam_fasta.jsonl.gz. With an OUTPUT path,
#   writes there; add --link to also symlink it under data/.
#
# Examples:
#   scripts/parse_pfam_fasta.sh
#   scripts/parse_pfam_fasta.sh ${outdir}/pfam_fasta.jsonl.gz --link
#
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
name="pfam_fasta.jsonl.gz"
out=""; link=0
for a in "$@"; do
    if [[ "$a" == "--link" ]]; then link=1; else out="$a"; fi
done
[[ -n "$out" ]] || { out="data/$name"; link=0; }
mkdir -p data "$(dirname "$out")"
# ─────────────────────────────────────────────────────────────────────────────

bioparsers pfam-fasta databases/pfam/Pfam-A.fasta.gz -o "$out" --gzip --progress 1000000

if [[ "$link" -eq 1 ]]; then
    ln -sfn "$(readlink -f "$out")" "data/$name"
    echo "linked data/$name -> $out"
fi

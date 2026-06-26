#!/usr/bin/env bash
#
# Parse the full Pfam-A Stockholm release (Pfam-A.full.gz) to gzipped JSONL, one 
# record per family (family-level metadata only; pass --with-member-accessions 
# and/or --with-member-sequences to bioparsers to include members).
#
# Usage: scripts/parse_pfam_full.sh [OUTPUT] [--link]
#   With no OUTPUT, writes data/pfam.jsonl.gz. With an OUTPUT path, writes
#   there; add --link to also symlink it under data/.
#
# Examples:
#   scripts/parse_pfam_full.sh
#   scripts/parse_pfam_full.sh ${outdir}/pfam.jsonl.gz --link
#   scripts/parse_pfam_full.sh pfam_with_sequences.jsonl.gz --with-member-sequences
#
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
name="pfam.jsonl.gz"
out=""; link=0
for a in "$@"; do
    if [[ "$a" == "--link" ]]; then link=1; else out="$a"; fi
done
[[ -n "$out" ]] || { out="data/$name"; link=0; }
mkdir -p data "$(dirname "$out")"
# ─────────────────────────────────────────────────────────────────────────────

bioparsers pfam databases/pfam/Pfam-A.full.gz -o "$out" --gzip --progress 1000

if [[ "$link" -eq 1 ]]; then
    ln -sfn "$(readlink -f "$out")" "data/$name"
    echo "linked data/$name -> $out"
fi

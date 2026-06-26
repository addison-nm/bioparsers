#!/usr/bin/env bash
#
# Parse the TrEMBL release (uniprot_trembl.dat.gz) to gzipped JSONL.
#
# Usage: scripts/parse_trembl.sh [OUTPUT] [--link]
#   With no OUTPUT, writes data/uniprot_trembl.jsonl.gz. With an OUTPUT path,
#   writes there; add --link to also symlink it under data/.
#
# Examples:
#   scripts/parse_trembl.sh
#   scripts/parse_trembl.sh ${outdir}/trembl.jsonl.gz --link
#
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
name="uniprot_trembl.jsonl.gz"
out=""; link=0
for a in "$@"; do
    if [[ "$a" == "--link" ]]; then link=1; else out="$a"; fi
done
[[ -n "$out" ]] || { out="data/$name"; link=0; }
mkdir -p data "$(dirname "$out")"
# ─────────────────────────────────────────────────────────────────────────────

bioparsers uniprot databases/trembl/uniprot_trembl.dat.gz -o "$out" --gzip --progress 1000000

if [[ "$link" -eq 1 ]]; then
    ln -sfn "$(readlink -f "$out")" "data/$name"
    echo "linked data/$name -> $out"
fi

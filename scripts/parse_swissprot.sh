#!/usr/bin/env bash
#
# Parse the Swiss-Prot release (uniprot_sprot.dat.gz) to gzipped JSONL.
#
# Usage: scripts/parse_swissprot.sh [OUTPUT] [--link]
#   With no OUTPUT, writes data/uniprot_sprot.jsonl.gz. With an OUTPUT path,
#   writes there; add --link to also symlink it under data/.
#
# Examples:
#   scripts/parse_swissprot.sh
#   scripts/parse_swissprot.sh ${outdir}/sprot.jsonl.gz --link
#
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
name="uniprot_sprot.jsonl.gz"
out=""; link=0
for a in "$@"; do
    if [[ "$a" == "--link" ]]; then link=1; else out="$a"; fi
done
[[ -n "$out" ]] || { out="data/$name"; link=0; }
mkdir -p data "$(dirname "$out")"
# ─────────────────────────────────────────────────────────────────────────────

bioparsers uniprot databases/swissprot/uniprot_sprot.dat.gz -o "$out" --gzip --progress

if [[ "$link" -eq 1 ]]; then
    ln -sfn "$(readlink -f "$out")" "data/$name"
    echo "linked data/$name -> $out"
fi

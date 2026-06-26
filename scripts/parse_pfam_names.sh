#!/usr/bin/env bash
#
# Write the Pfam accession -> family-name table as a two-column TSV
# (PF#####<TAB>name), sorted by accession, read from Pfam-A.hmm.gz.
#
# Usage: scripts/parse_pfam_names.sh [OUTPUT] [--link]
#   With no OUTPUT, writes data/pfam_names.tsv. With an OUTPUT path, writes
#   there; add --link to also symlink it under data/.
#
# Examples:
#   scripts/parse_pfam_names.sh
#   scripts/parse_pfam_names.sh ${outdir}/pfam_names.tsv --link
#
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
name="pfam_names.tsv"
out=""; link=0
for a in "$@"; do
    if [[ "$a" == "--link" ]]; then link=1; else out="$a"; fi
done
[[ -n "$out" ]] || { out="data/$name"; link=0; }
mkdir -p data "$(dirname "$out")"
# ─────────────────────────────────────────────────────────────────────────────

python scripts/pfam_names_to_tsv.py databases/pfam/Pfam-A.hmm.gz "$out"

if [[ "$link" -eq 1 ]]; then
    ln -sfn "$(readlink -f "$out")" "data/$name"
    echo "linked data/$name -> $out"
fi

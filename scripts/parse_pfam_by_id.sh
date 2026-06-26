#!/usr/bin/env bash
#
# Parse selected Pfam families from the full Stockholm release (Pfam-A.full.gz)
# to JSONL, with each family's members and their sequences attached. Default:
# one file per Pfam ID. With --join: a single unioned file instead.
#
# Usage: scripts/parse_pfam_by_id.sh PFXXXXX [PFYYYYY ...] [OUTPUT] [--join] [--link]
#   Default (per ID): OUTPUT is a directory (default data/); writes
#     <OUTPUT>/pfam_<ID>.jsonl per family.
#   --join: OUTPUT is a file (default data/pfam_<ids>.jsonl); writes one union.
#   --link: also symlink the output(s) under data/.
#
# Examples:
#   scripts/parse_pfam_by_id.sh PF00018 PF07714
#   scripts/parse_pfam_by_id.sh PF00018 PF07714 --join ${outdir}/sh3_kinase.jsonl --link
#
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
out=""; link=0; join=0; ids=()
for a in "$@"; do
    case "$a" in
        --link)    link=1 ;;
        --join)    join=1 ;;
        PF[0-9]*)  ids+=("$a") ;;
        *)         out="$a" ;;
    esac
done
[[ ${#ids[@]} -gt 0 ]] || { echo "usage: $0 PFXXXXX [PFYYYYY ...] [OUTPUT] [--join] [--link]" >&2; exit 2; }
pfam_args=(); for id in "${ids[@]}"; do pfam_args+=(--pfam-id "$id"); done
db="databases/pfam/Pfam-A.full.gz"
# ─────────────────────────────────────────────────────────────────────────────

if [[ "$join" -eq 1 ]]; then
    name="pfam_$(IFS=_; echo "${ids[*]}").jsonl"
    [[ -n "$out" ]] || { out="data/$name"; link=0; }
    mkdir -p data "$(dirname "$out")"
    bioparsers pfam "$db" "${pfam_args[@]}" --join --with-member-sequences \
        -o "$out" --progress 1000
    if [[ "$link" -eq 1 ]]; then
        ln -sfn "$(readlink -f "$out")" "data/$name"
        echo "linked data/$name -> $out"
    fi
else
    dir="${out:-data}"
    [[ -n "$out" ]] || link=0
    mkdir -p data "$dir"
    bioparsers pfam "$db" "${pfam_args[@]}" --with-member-sequences \
        -o "$dir" --progress 1000
    if [[ "$link" -eq 1 ]]; then
        for id in "${ids[@]}"; do
            base="${id%%.*}"
            f="$dir/pfam_$base.jsonl"
            if [[ -e "$f" ]]; then
                ln -sfn "$(readlink -f "$f")" "data/pfam_$base.jsonl"
                echo "linked data/pfam_$base.jsonl -> $f"
            fi
        done
    fi
fi

#!/bin/bash
# Unsorted, untruncated query output with every metric column, for offline sort/filter selection.
# Held-out M-CSA queries (q755 minus q250). Usage: collect_heldout.sh <binary>
# Output: build/raw/<run>/<query>.tsv
set -uo pipefail
source "$(dirname "$0")/../../scripts/lib.sh"
D=$(cd "$(dirname "$0")/.." && pwd)
export XBIN=$1 FD_RAW=1 RAW=$D/build/raw_heldout IDX_MCSA IDX_HSAPIENS REPO
export PM=tid,node_count,idf,rmsd,e_value,tm_score,gdt_ts,gdt_ha,chamfer_distance,hausdorff_distance,drmsd,max_dist_deviation
export PS=tid,idf,total_match_count,node_count,edge_count,max_node_cov,min_rmsd,min_drmsd,nres,plddt
MUT=$D/build/heldout/mutant_spec.tsv
SPEC=$D/build/heldout/spec.tsv

mcsa_one() { # <spec line>   env: RUN KIND EXTRA
  IFS=$'\t' read -r m p r res _ <<< "$1"
  local out=$RAW/$RUN; [ -f "$out/$m.ok" ] && return 0
  [ "$KIND" = star ] && r=$(echo ",$r," | sed "s/,$res,/,$res:*,/; s/^,//; s/,\$//")
  /usr/bin/time -f "%e" -o "$out/$m.time" timeout 900 \
    $XBIN query -i "$IDX_MCSA" -p "$p" -q "$r" -t 4 --top 6000 $EXTRA > "$out/$m.tsv" 2>/dev/null && touch "$out/$m.ok"
}
export -f mcsa_one

mcsa_run() { # <run> <spec> <kind> <extra...>
  export RUN=$1 KIND=$3 EXTRA="${*:4}"; mkdir -p "$RAW/$RUN"
  local t0=$(date +%s)
  xargs -d'\n' -a "$2" -P 5 -I{} bash -c 'mcsa_one "$@"' _ {}
  echo "$RUN wall=$(( $(date +%s)-t0 ))s ok=$(ls $RAW/$RUN/*.ok 2>/dev/null | wc -l)"
}

motif_run() { # <run> <extra...>
  local run=$1; shift; mkdir -p "$RAW/$run"
  local t0=$(date +%s)
  while IFS=$'\t' read -r qid label motif mode pdb res threads flags; do
    [ -f "$RAW/$run/$qid.ok" ] && continue
    ( cd "$REPO" && /usr/bin/time -f "%e" -o "$RAW/$run/$qid.time" \
        $XBIN query -p "$pdb" -q "$res" -i "$IDX_HSAPIENS" -t 12 "$@" ) > "$RAW/$run/$qid.tsv" 2>/dev/null \
      && touch "$RAW/$run/$qid.ok"
  done < <(awk -F'\t' '!seen[$5 FS $6]++' "$SPEC_MOTIF")
  echo "$run wall=$(( $(date +%s)-t0 ))s"
}

mcsa_run mc_pm_def   $SPEC orig --format-output $PM
mcsa_run mc_ps_def   $SPEC orig --per-structure --format-output $PS
mcsa_run mc_pm_sens  $SPEC orig --format-output $PM --sensitive
mcsa_run mc_pm_star  $MUT  star --format-output $PM
echo done

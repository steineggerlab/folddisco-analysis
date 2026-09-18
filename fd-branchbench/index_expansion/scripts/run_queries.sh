#!/bin/bash
# Site queries against the motif indices. Usage: run_queries.sh [config ...]  (default: all)
# Output: build/out/<config>/<id>.tsv (tid, node_count, idf, rmsd; idf:desc,rmsd:asc) and .time
set -uo pipefail
D=$(cd "$(dirname "$0")/.." && pwd)
export BIN=${BIN:-$D/../build/bin/folddisco-2cd6e60}
export BUILD=$D/build TOP=${TOP:-500} NT=${NT:-2}
NP=${NP:-4}

# name<TAB>query set (orig|mut|star)<TAB>index<TAB>extra flags
configs() {
  for set in orig mut; do
    printf '%s\t%s\t%s\t%s\n' "${set}_I0"       $set I0 ""
    printf '%s\t%s\t%s\t%s\n' "${set}_I0_sens"  $set I0 "--sensitive"
    printf '%s\t%s\t%s\t%s\n' "${set}_I0_aa"    $set I0 "--aa-subst blosum62"
    printf '%s\t%s\t%s\t%s\n' "${set}_I1"       $set I1 ""
    printf '%s\t%s\t%s\t%s\n' "${set}_I2"       $set I2 ""
    printf '%s\t%s\t%s\t%s\n' "${set}_I3"       $set I3 ""
    printf '%s\t%s\t%s\t%s\n' "${set}_I4"       $set I4 ""
  done
  printf '%s\t%s\t%s\t%s\n' "mut_I0_star" star I0 ""
}

run_one() { # <site query line>
  IFS=$'\t' read -r id pdb mut_pdb site _ _ mres <<< "$1"
  local out=$BUILD/out/$CFG
  [ -f "$out/$id.ok" ] && return 0
  case $SET in
    orig) p=$pdb ;;
    mut)  p=$mut_pdb ;;
    star) p=$mut_pdb; site=$(echo ",$site," | sed "s/,$mres,/,$mres:*,/; s/^,//; s/,\$//") ;;
  esac
  ( cd "$BUILD" && /usr/bin/time -f "%e %M" -o "$out/$id.time" \
      $BIN query -i index/$IDX -p "$p" -q "$site" -t "$NT" --top "$TOP" \
      --format-output tid,node_count,idf,rmsd $EXTRA > "$out/$id.tsv" 2> "$out/$id.err" ) \
    && touch "$out/$id.ok"
}
export -f run_one

want=" $* "
while IFS=$'\t' read -r cfg set idx extra; do
  [ $# -gt 0 ] && [[ "$want" != *" $cfg "* ]] && continue
  export CFG=$cfg SET=$set IDX=$idx EXTRA=$extra
  mkdir -p "$BUILD/out/$cfg"
  t0=$(date +%s)
  xargs -d'\n' -a "${SPEC:-$BUILD/site_queries.tsv}" -P "$NP" -I{} bash -c 'run_one "$@"' _ {}
  echo "$cfg wall=$(( $(date +%s) - t0 ))s ok=$(ls "$BUILD/out/$cfg"/*.ok 2>/dev/null | wc -l)"
done < <(configs)

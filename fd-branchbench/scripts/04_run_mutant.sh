#!/bin/bash
# Substitution benchmark: M-CSA queries with one residue renamed to its closest BLOSUM62
# alternative, same answer sets. Exact search loses answers; substitution should win them back.
# Usage: 04_run_mutant.sh [60|250|755]   (default 250)
# Env: NP (5), NT (4), TOP (6000), TIMEOUT (600s per query).
set -uo pipefail
source "$(dirname "$0")/lib.sh"
SUBSET=${1:-250}; NP=${NP:-5}; NT=${NT:-4}; TOP=${TOP:-6000}; TIMEOUT=${TIMEOUT:-600}
SPEC=$MCSA_DIR/queryspec_${SUBSET}.tsv
MUT_DIR=$BENCH/build/mutant_$SUBSET
MUT_SPEC=$MUT_DIR/spec.tsv
OUT=$RESULT/mutant_$SUBSET; mkdir -p "$OUT"
F=$RESULT/mutant_fp1_${SUBSET}.tsv

[ -s "$MUT_SPEC" ] || "$VENV" "$BENCH/scripts/make_mutant_queries.py" "$SPEC" "$MUT_DIR/pdb" "$MUT_SPEC"
echo "mutant queries: $(wc -l < "$MUT_SPEC") (one renamed residue each)"
cut -f5 "$MUT_SPEC" | sort | uniq -c | sort -nr | head -5

run_one() { # <mutant-spec-line>
  IFS=$'\t' read -r m mut_pdb r res _ <<< "$1"
  local d=$OUT/$CFG pdb=$mut_pdb
  [ -f "$d/$m.ok" ] && return 0
  case $KIND in
    orig) pdb=$(awk -F'\t' -v m="$m" '$1==m{print $2}' "$SPEC") ;;
    star) r=$(echo ",$r," | sed "s/,$res,/,$res:*,/; s/^,//; s/,\$//") ;;
  esac
  /usr/bin/time -f "%e %M" -o "$d/$m.time" timeout "$TIMEOUT" \
    $BIN_BRANCH query -i "$IDX_MCSA" -p "$pdb" -q "$r" -t "$NT" --top "$TOP" $EXTRA > "$d/$m.tsv" 2>"$d/$m.err" \
    && touch "$d/$m.ok" || rm -f "$d/$m.tsv"
}
export -f run_one

score_one() { # <mutant-spec-line>
  IFS=$'\t' read -r m _ _ _ change <<< "$1"
  local d=$OUT/$CFG
  [ -f "$d/$m.ok" ] || return 0
  "$BIN_BRANCH" benchmark -i "$IDX_MCSA" -r "$d/$m.tsv" -a "$ANS_MCSA/${m}_answer.tsv" --fp 1 2>/dev/null \
    | awk -F'\t' -v c="$CFG" -v m="$m" -v ch="$change" -v rt="$(awk 'NR==1{print $1}' "$d/$m.time")" \
        'NF>=17{print c"\t"m"\t"ch"\t"$6"\t"$10"\t"$15"\t"rt; exit}'
}
export -f score_one

printf 'config\tquery\tchange\tanswer_len\ttp\tsens_at_1fp\truntime\n' > "$F"
while IFS=$'\t' read -r cfg kind extra; do
  export CFG=$cfg KIND=$kind EXTRA=$extra OUT SPEC IDX_MCSA NT TOP TIMEOUT ANS_MCSA BIN_BRANCH
  mkdir -p "$OUT/$cfg"
  t0=$(date +%s)
  xargs -d'\n' -a "$MUT_SPEC" -P "$NP" -I{} bash -c 'run_one "$@"' _ {}
  t1=$(date +%s)
  xargs -d'\n' -a "$MUT_SPEC" -P "$NP" -I{} bash -c 'score_one "$@"' _ {} | sort >> "$F"
  echo "$cfg: wall=$((t1-t0))s scored=$(grep -c "^$cfg	" "$F")"
done < <(mutant_configs)
echo "wrote $F"

#!/bin/bash
# M-CSA catalytic-site benchmark: Sens@1FP per query, all three configs.
# Usage: 02_run_mcsa.sh [60|250|755]   (default 250)
# Env: NP (parallel queries, 5), NT (threads each, 4), TOP (6000), TIMEOUT (600s per query).
# Scored with the BRANCH binary throughout so the metric code is constant.
set -uo pipefail
source "$(dirname "$0")/lib.sh"
SUBSET=${1:-250}; NP=${NP:-5}; NT=${NT:-4}; TOP=${TOP:-6000}; TIMEOUT=${TIMEOUT:-600}
SPEC=$MCSA_DIR/queryspec_${SUBSET}.tsv
[ -f "$SPEC" ] || { echo "no such subset spec: $SPEC"; exit 1; }
OUT=$RESULT/mcsa_$SUBSET; mkdir -p "$OUT"
F=$RESULT/mcsa_fp1_${SUBSET}.tsv
printf 'config\tquery\thits\tresult_len\tanswer_len\ttp\tfp\tfn\tprecision\tsens_at_1fp\truntime\tmax_rss_kb\n' > "$F"

run_one() { # <query-spec-line>
  IFS=$'\t' read -r m p r _ <<< "$1"
  local d=$OUT/$CFG
  [ -f "$d/$m.ok" ] && return 0
  /usr/bin/time -f "%e %M" -o "$d/$m.time" timeout "$TIMEOUT" \
    $BIN query -i "$IDX_MCSA" -p "$p" -q "$r" -t "$NT" --top "$TOP" $EXTRA > "$d/$m.tsv" 2>"$d/$m.err"
  rc=$?
  # A timeout (124) leaves a truncated result: discard it rather than score a partial list.
  if [ $rc -ne 0 ]; then
    rm -f "$d/$m.tsv"
    printf '%s\t%s\n' "$m" "$([ $rc -eq 124 ] && echo timeout || echo "rc=$rc")" >> "$OUT/${CFG}_failed.txt"
    return 0
  fi
  touch "$d/$m.ok"
}
export -f run_one

score_one() { # <query-spec-line>
  IFS=$'\t' read -r m p r _ <<< "$1"
  local d=$OUT/$CFG
  [ -f "$d/$m.ok" ] || return 0
  local hits rt rss
  hits=$(wc -l < "$d/$m.tsv")
  rt=$(awk 'NR==1{print $1}' "$d/$m.time" 2>/dev/null); rss=$(awk 'NR==1{print $2}' "$d/$m.time" 2>/dev/null)
  "$BIN_BRANCH" benchmark -i "$IDX_MCSA" -r "$d/$m.tsv" -a "$ANS_MCSA/${m}_answer.tsv" --fp 1 2>/dev/null \
    | awk -F'\t' -v c="$CFG" -v m="$m" -v h="$hits" -v rt="$rt" -v rss="$rss" \
        'NF>=17{print c"\t"m"\t"h"\t"$5"\t"$6"\t"$10"\t"$12"\t"$13"\t"$14"\t"$15"\t"rt"\t"rss; exit}'
}
export -f score_one

while IFS=$'\t' read -r cfg bin extra; do
  export CFG=$cfg BIN=$bin EXTRA=$extra OUT IDX_MCSA NT TOP TIMEOUT ANS_MCSA BIN_BRANCH
  mkdir -p "$OUT/$cfg"
  t0=$(date +%s)
  xargs -d'\n' -a "$SPEC" -P "$NP" -I{} bash -c 'run_one "$@"' _ {}
  t1=$(date +%s)
  xargs -d'\n' -a "$SPEC" -P "$NP" -I{} bash -c 'score_one "$@"' _ {} | sort >> "$F"
  nfail=0; [ -f "$OUT/${cfg}_failed.txt" ] && nfail=$(wc -l < "$OUT/${cfg}_failed.txt")
  echo "$cfg: wall=$((t1-t0))s scored=$(grep -c "^$cfg	" "$F") failed=$nfail"
done < <(configs)

# Two identity measures. Byte-identity also covers the matched-residue chain labels,
# which the branch's multichar-chainID work changes on multi-chain entries; scoring keys
# on field 1 only, so fields 1-4 (id, node_count, idf, rmsd) are what the metrics see.
echo "== default-vs-master identity ($SUBSET) =="
byte=0; rank=0; n=0
for m in $(cut -f1 "$SPEC"); do
  [ -f "$OUT/master/$m.tsv" ] && [ -f "$OUT/default/$m.tsv" ] || continue
  n=$((n+1))
  cmp -s "$OUT/master/$m.tsv" "$OUT/default/$m.tsv" && byte=$((byte+1))
  cmp -s <(cut -f1-4 "$OUT/master/$m.tsv") <(cut -f1-4 "$OUT/default/$m.tsv") && rank=$((rank+1))
done
printf 'byte-identical:       %s/%s\n' "$byte" "$n"
printf 'rank+score identical: %s/%s\n' "$rank" "$n"
printf 'subset\tn\tbyte_identical\trank_identical\n%s\t%s\t%s\t%s\n' "$SUBSET" "$n" "$byte" "$rank" > "$RESULT/mcsa_identity_${SUBSET}.tsv"
echo "wrote $F"

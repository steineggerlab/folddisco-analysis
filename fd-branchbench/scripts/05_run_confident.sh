#!/bin/bash
# What --confident returns, against no filter and against the protocol's per-query filters.
# Precision/recall/F1 of the returned set (not a ranked metric). Usage: 05_run_confident.sh [250]
set -uo pipefail
source "$(dirname "$0")/lib.sh"
SUBSET=${1:-250}; NP=${NP:-5}; NT=${NT:-4}; TOP=${TOP:-6000}
SPEC=$MCSA_DIR/queryspec_${SUBSET}.tsv
OUT=$RESULT/confident_$SUBSET; mkdir -p "$OUT"
M=$RESULT/confident_motif.tsv
F=$RESULT/confident_mcsa_${SUBSET}.tsv

# name<TAB>flags; "protocol" is filled per query from the motif spec
motif_configs() {
  printf '%s\t%s\n' unfiltered     ""
  printf '%s\t%s\n' protocol       "PROTOCOL"
  printf '%s\t%s\n' confident      "--confident"
  printf '%s\t%s\n' confident_sens "--confident --sensitive"
}
mcsa_configs() {
  printf '%s\t%s\n' unfiltered     ""
  printf '%s\t%s\n' confident      "--confident"
  printf '%s\t%s\n' confident_sens "--confident --sensitive"
}

echo "== motif"
printf 'config\tqid\tmotif\thits\tanswer_len\ttp\tprecision\trecall\tf1\truntime\n' > "$M"
while IFS=$'\t' read -r qid label motif mode pdb res threads matched_flags; do
  case $motif in zinc) ans=$ANS_ZINC;; serine) ans=$ANS_SERINE;; esac
  while IFS=$'\t' read -r cfg flags; do
    [ "$flags" = PROTOCOL ] && flags=$matched_flags
    t=$OUT/motif_${cfg}_${qid}.tsv
    s=$( { /usr/bin/time -f "%e" -o "$OUT/motif_${cfg}_${qid}.time" \
        bash -c "cd $REPO && $BIN_BRANCH query -p $pdb -q $res -i $IDX_HSAPIENS -t 12 $flags" > "$t" 2>/dev/null; } )
    read -r alen tp prec rec f1 <<< "$("$BIN_BRANCH" benchmark -i "$IDX_HSAPIENS" -r "$t" -a "$ans" --afdb-to-uniprot 2>/dev/null \
      | awk -F'\t' 'NF>=17{print $6,$10,$14,$15,$17; exit}')"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$cfg" "$qid" "$motif" "$(wc -l < "$t")" \
      "$alen" "$tp" "$prec" "$rec" "$f1" "$(cat "$OUT/motif_${cfg}_${qid}.time")" >> "$M"
  done < <(motif_configs)
done < <(awk -F'\t' '$4=="matched"' "$SPEC_MOTIF")
echo "wrote $M"

echo "== M-CSA q$SUBSET"
run_one() {
  IFS=$'\t' read -r m p r _ <<< "$1"
  local d=$OUT/$CFG
  [ -f "$d/$m.ok" ] && return 0
  /usr/bin/time -f "%e" -o "$d/$m.time" timeout 900 \
    $BIN_BRANCH query -i "$IDX_MCSA" -p "$p" -q "$r" -t "$NT" --top "$TOP" $EXTRA > "$d/$m.tsv" 2>/dev/null \
    && touch "$d/$m.ok"
}
score_one() {
  IFS=$'\t' read -r m p r _ <<< "$1"
  local d=$OUT/$CFG
  [ -f "$d/$m.ok" ] || return 0
  "$BIN_BRANCH" benchmark -i "$IDX_MCSA" -r "$d/$m.tsv" -a "$ANS_MCSA/${m}_answer.tsv" 2>/dev/null \
    | awk -F'\t' -v c="$CFG" -v m="$m" -v h="$(wc -l < "$d/$m.tsv")" -v rt="$(cat "$d/$m.time")" \
        'NF>=17{print c"\t"m"\t"h"\t"$6"\t"$10"\t"$14"\t"$15"\t"$17"\t"rt; exit}'
}
export -f run_one score_one
printf 'config\tquery\thits\tanswer_len\ttp\tprecision\trecall\tf1\truntime\n' > "$F"
while IFS=$'\t' read -r cfg flags; do
  export CFG=$cfg EXTRA=$flags OUT IDX_MCSA NT TOP ANS_MCSA BIN_BRANCH
  mkdir -p "$OUT/$cfg"
  t0=$(date +%s)
  xargs -d'\n' -a "$SPEC" -P "$NP" -I{} bash -c 'run_one "$@"' _ {}
  t1=$(date +%s)
  xargs -d'\n' -a "$SPEC" -P "$NP" -I{} bash -c 'score_one "$@"' _ {} | sort >> "$F"
  echo "$cfg: wall=$((t1-t0))s scored=$(grep -c "^$cfg	" "$F")"
done < <(mcsa_configs)
echo "wrote $F"

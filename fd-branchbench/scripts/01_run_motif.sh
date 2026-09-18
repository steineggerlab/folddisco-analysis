#!/bin/bash
# Serine-peptidase + zinc-finger motif search on the human index: P/R/F1, Sens@1FP, runtime.
# Every config is scored with the BRANCH binary so the metric code is constant across configs.
# Env: REPEATS (default 13), WARMUP (default 2).
set -uo pipefail
source "$(dirname "$0")/lib.sh"
REPEATS=${REPEATS:-13}; WARMUP=${WARMUP:-2}; RESUME=${RESUME:-0}
OUT=$RESULT/motif; mkdir -p "$OUT"
M=$RESULT/motif_metrics.tsv

# RESUME=1 keeps rows already in the metrics file and fills only what is missing,
# so a run killed part-way (e.g. by the OOM reaper) does not have to start over.
if [ "$RESUME" = 1 ] && [ -s "$M" ]; then
  echo "resuming: $(( $(wc -l < "$M") - 1 )) rows already present"
else
  printf 'config\tqid\tlabel\tmotif\tmode\thits\tanswer_len\ttp\tfp\tfn\tprecision\trecall\tf1\ttp_at_1fp\tsens_at_1fp\truntime_med\truntime_min\truntime_max\tn_repeat\n' > "$M"
fi

score() { # score <binary> <result> <answer> <fp-args...>
  "$1" benchmark -i "$IDX_HSAPIENS" -r "$2" -a "$3" --afdb-to-uniprot "${@:4}" 2>/dev/null \
    | awk -F'\t' 'NF>=17{print; exit}'
}

while IFS=$'\t' read -r qid label motif mode pdb res threads flags; do
  case $motif in zinc) ans=$ANS_ZINC;; serine) ans=$ANS_SERINE;; esac
  while IFS=$'\t' read -r cfg bin extra; do
    if [ "$RESUME" = 1 ] && grep -qP "^${cfg}\t${qid}\t" "$M" 2>/dev/null; then continue; fi
    res_tsv=$OUT/${cfg}_${qid}.tsv
    # correctness run
    ( cd "$REPO" && $bin query -p "$pdb" -q "$res" -i "$IDX_HSAPIENS" -t "$threads" $flags $extra ) > "$res_tsv" 2>"$OUT/${cfg}_${qid}.err"
    hits=$(wc -l < "$res_tsv")
    prf=$(score "$BIN_BRANCH" "$res_tsv" "$ans")
    fp1=$(score "$BIN_BRANCH" "$res_tsv" "$ans" --fp 1)
    read -r alen tp fp fn prec rec f1 <<< "$(echo "$prf" | awk -F'\t' '{print $6,$10,$12,$13,$14,$15,$17}')"
    read -r tp1 sens1      <<< "$(echo "$fp1" | awk -F'\t' '{print $10,$15}')"
    # timing
    for _ in $(seq "$WARMUP"); do ( cd "$REPO" && $bin query -p "$pdb" -q "$res" -i "$IDX_HSAPIENS" -t "$threads" $flags $extra ) >/dev/null 2>&1; done
    times=$(for _ in $(seq "$REPEATS"); do
      s=$(date +%s.%N)
      ( cd "$REPO" && $bin query -p "$pdb" -q "$res" -i "$IDX_HSAPIENS" -t "$threads" $flags $extra ) >/dev/null 2>&1
      e=$(date +%s.%N); awk -v s="$s" -v e="$e" 'BEGIN{printf "%.4f\n", e-s}'
    done)
    med=$(echo "$times" | median); mn=$(echo "$times" | sort -g | head -1); mx=$(echo "$times" | sort -g | tail -1)
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$cfg" "$qid" "$label" "$motif" "$mode" "$hits" "$alen" "$tp" "$fp" "$fn" "$prec" "$rec" "$f1" "$tp1" "$sens1" "$med" "$mn" "$mx" "$REPEATS" >> "$M"
    echo "  $cfg $qid hits=$hits f1=$f1 sens@1fp=$sens1 t=${med}s"
  done < <(configs)
done < "$SPEC_MOTIF"

# Identity control: branch default vs master. Byte-identity also covers matched-residue
# chain labels; fields 1-4 are what the metrics actually key on.
echo "== default-vs-master identity =="
byte=0; rank=0; n=0
for qid in $(cut -f1 "$SPEC_MOTIF"); do
  n=$((n+1))
  cmp -s "$OUT/master_${qid}.tsv" "$OUT/default_${qid}.tsv" && byte=$((byte+1)) || echo "  byte-differs: $qid"
  cmp -s <(cut -f1-4 "$OUT/master_${qid}.tsv") <(cut -f1-4 "$OUT/default_${qid}.tsv") && rank=$((rank+1))
done
printf 'byte-identical:       %s/%s\n' "$byte" "$n"
printf 'rank+score identical: %s/%s\n' "$rank" "$n"
printf 'n\tbyte_identical\trank_identical\n%s\t%s\t%s\n' "$n" "$byte" "$rank" > "$RESULT/motif_identity.tsv"
echo "wrote $M"

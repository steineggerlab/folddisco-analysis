#!/bin/bash
# End-to-end: setup, benchmarks, answer classes, plots. Usage: run_all.sh [60|250|755]  (default 250)
set -euo pipefail
D=$(dirname "$0"); SUBSET=${1:-250}
source "$D/lib.sh"
"$D/00_setup.sh"
"$D/01_run_motif.sh"
"$D/02_run_mcsa.sh" "$SUBSET"
"$D/04_run_mutant.sh" "$SUBSET"
"$D/05_run_confident.sh" "$SUBSET"
"$VENV" "$D/mcsa_answer_classes.py" "$MCSA_DIR/queryspec_$SUBSET.tsv" "$RESULT/mcsa_answer_classes_$SUBSET.tsv"
"$VENV" "$D/mcsa_class_metrics.py" "$RESULT/mcsa_answer_classes_$SUBSET.tsv" "$RESULT/mcsa_$SUBSET" "$RESULT/mcsa_classes_$SUBSET.tsv"
"$VENV" "$D/03_plot.py" --subset "$SUBSET"

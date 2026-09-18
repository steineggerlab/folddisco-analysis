#!/bin/bash
# Shared paths and config matrix for the master-vs-feature-integration comparison.
#
# Paths come from the environment so the suite runs off this machine. Override any of them
# in scripts/config.local.sh (untracked) or in the environment; defaults assume the layout
# that scripts/00_setup.sh checks and the README documents:
#   $FOLDDISCO_REPO          folddisco source checkout (cargo build --release)
#   $FOLDDISCO_BENCH_DATA    indices and answer sets (default: <bench>/data)
set -u

BENCH=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
[ -f "$BENCH/scripts/config.local.sh" ] && source "$BENCH/scripts/config.local.sh"

REPO=${FOLDDISCO_REPO:-$BENCH/../folddisco}                # feature-integration checkout
DATA=${FOLDDISCO_BENCH_DATA:-$BENCH/data}                  # benchmark inputs
WT_MASTER=${WT_MASTER:-$BENCH/build/wt-master}             # detached worktree at origin/master

BIN_MASTER=${BIN_MASTER:-$WT_MASTER/target/release/folddisco}
BIN_BRANCH=${BIN_BRANCH:-$REPO/target/release/folddisco}

# Indices
IDX_HSAPIENS=${IDX_HSAPIENS:-$DATA/index/h_sapiens_folddisco}   # 23,391 AFDB human structures
MCSA_DIR=${MCSA_DIR:-$DATA/mcsa_rebuild}
IDX_MCSA=${IDX_MCSA:-$MCSA_DIR/index/mcsa_subset_folddisco}     # 62,122 PDB entries

# Answer sets
ANS_ZINC=${ANS_ZINC:-$DATA/answers/zincfinger_answer.tsv}        # 761 accessions
ANS_SERINE=${ANS_SERINE:-$DATA/answers/serinepeptidase_answer.tsv} # 124 accessions
ANS_MCSA=${ANS_MCSA:-$DATA/mcsa/mcsa_answers}

SPEC_MOTIF=$BENCH/spec/motif_queries.tsv
RESULT=$BENCH/result
FIGURE=$BENCH/figure
VENV=$BENCH/.venv/bin/python

# Exported for the Python steps, which resolve the same inputs.
export FOLDDISCO_BENCH_DATA=$DATA MCSA_DIR IDX_MCSA ANS_MCSA ANS_ZINC ANS_SERINE

# Config matrix: name<TAB>binary<TAB>extra flags
# "default" is the branch with no new flags: the rank-identity control against master.
# --sensitive is geometric (--expand-radius 2); --aa-subst is the chemical axis.
# CONFIGS=core limits the matrix to the three master-vs-branch configs.
configs() {
  printf '%s\t%s\t%s\n' master       "$BIN_MASTER" ""
  printf '%s\t%s\t%s\n' default      "$BIN_BRANCH" ""
  printf '%s\t%s\t%s\n' sensitive    "$BIN_BRANCH" "--sensitive"
  [ "${CONFIGS:-all}" = core ] && return 0
  printf '%s\t%s\t%s\n' aa_blosum62  "$BIN_BRANCH" "--aa-subst blosum62"
  printf '%s\t%s\t%s\n' aa_group     "$BIN_BRANCH" "--aa-subst group"
  printf '%s\t%s\t%s\n' aa_size      "$BIN_BRANCH" "--aa-subst size"
  printf '%s\t%s\t%s\n' sensitive_aa "$BIN_BRANCH" "--sensitive --aa-subst blosum62"
}

# Mutant-query matrix (04_run_mutant.sh): name<TAB>query<TAB>extra flags.
# query is "orig" (the M-CSA motif) or "mut" (one residue renamed, see make_mutant_queries.py);
# "star" appends :* to the mutated residue.
mutant_configs() {
  printf '%s\t%s\t%s\n' orig           orig ""
  printf '%s\t%s\t%s\n' mut            mut  ""
  printf '%s\t%s\t%s\n' mut_star       star ""
  printf '%s\t%s\t%s\n' mut_blosum62   mut  "--aa-subst blosum62"
  printf '%s\t%s\t%s\n' mut_sens_star  star "--sensitive"
}

# median of stdin numbers
median() { sort -g | awk '{v[NR]=$1} END{if(NR==0){print "NA";exit} m=int(NR/2); print (NR%2)?v[m+1]:(v[m]+v[m+1])/2}'; }

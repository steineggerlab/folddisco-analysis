#!/bin/bash
# Build the binaries, create the plotting venv, verify inputs. Idempotent; safe to rerun.
#
#   ./scripts/00_setup.sh            build branch + master binaries, venv, check inputs
#   SKIP_MASTER=1 ./scripts/00_setup.sh   branch only (no origin/master comparison)
#
# Paths come from scripts/lib.sh; put host-specific ones in scripts/config.local.sh.
set -uo pipefail
source "$(dirname "$0")/lib.sh"

fail=0
note() { printf '%s\n' "$*" >&2; }

echo "== tools =="
for t in cargo python3 git; do
  if command -v "$t" >/dev/null; then echo "ok   $t  $("$t" --version 2>&1 | head -1)"
  else echo "MISS $t"; fail=1; fi
done
[ "$fail" = 0 ] || { note "Install the missing tools first (Rust toolchain: https://rustup.rs)."; exit 1; }

echo "== folddisco source =="
if [ -d "$REPO/.git" ]; then
  echo "ok   $REPO"
else
  echo "MISS $REPO"
  note "Set FOLDDISCO_REPO to a folddisco checkout, e.g."
  note "  git clone https://github.com/steineggerlab/folddisco \$HOME/folddisco"
  note "  FOLDDISCO_REPO=\$HOME/folddisco ./scripts/00_setup.sh"
  exit 1
fi

echo "== worktree =="
if [ "${SKIP_MASTER:-0}" = 1 ]; then
  echo "skip (SKIP_MASTER=1): master configs will be missing from the figures"
elif [ -d "$WT_MASTER" ]; then
  echo "exists: $WT_MASTER"
elif git -C "$REPO" rev-parse --verify -q origin/master >/dev/null; then
  git -C "$REPO" worktree add "$WT_MASTER" "$(git -C "$REPO" rev-parse origin/master)"
else
  echo "MISS origin/master in $REPO"
  note "Fetch it (git -C $REPO fetch origin master) or rerun with SKIP_MASTER=1."
  exit 1
fi

echo "== build =="
[ "${SKIP_MASTER:-0}" = 1 ] || ( cd "$WT_MASTER" && cargo build --release -q )
( cd "$REPO" && cargo build --release -q )
for b in "$BIN_BRANCH" $([ "${SKIP_MASTER:-0}" = 1 ] || echo "$BIN_MASTER"); do
  [ -x "$b" ] && echo "ok   $b" || { echo "MISS $b"; fail=1; }
done

echo "== venv =="
if [ ! -x "$VENV" ]; then
  python3 -m venv "$BENCH/.venv"
  "$BENCH/.venv/bin/pip" -q install seaborn matplotlib pandas numpy
fi
"$VENV" -c "import seaborn,matplotlib,pandas; print('plot deps ok')"
# Optional: report/index.html -> PDF (scripts/make_report_pdf.py)
"$VENV" -c "import weasyprint" 2>/dev/null && echo "ok   weasyprint (report PDF)" \
  || echo "note weasyprint not installed; make_report_pdf.py needs it"

echo "== inputs =="
for f in "$IDX_HSAPIENS.lookup" "$IDX_MCSA.lookup" "$ANS_ZINC" "$ANS_SERINE" "$SPEC_MOTIF"; do
  [ -e "$f" ] && echo "ok   $f" || { echo "MISS $f"; fail=1; }
done
[ -d "$ANS_MCSA" ] && echo "ok   $ANS_MCSA ($(ls "$ANS_MCSA" | wc -l) answers)" || { echo "MISS $ANS_MCSA"; fail=1; }
for s in 60 250 755; do
  f=$MCSA_DIR/queryspec_$s.tsv
  [ -e "$f" ] && echo "ok   $f ($(wc -l < "$f") queries)" || { echo "MISS $f"; fail=1; }
done
if [ "$fail" != 0 ]; then
  note ""
  note "Missing inputs. Fetch the dataset record (https://doi.org/10.5281/zenodo.16679607) and"
  note "point FOLDDISCO_BENCH_DATA at the unpacked copy (see README, Setup),"
  note "or override single paths (IDX_HSAPIENS, IDX_MCSA, MCSA_DIR, ANS_MCSA, ANS_ZINC,"
  note "ANS_SERINE) in scripts/config.local.sh. Current root: $DATA"
fi

echo "== provenance =="
printf 'master\t%s\n' "$(git -C "$REPO" rev-parse --short origin/master 2>/dev/null || echo NA)"
printf 'branch\t%s\t%s\n' "$(git -C "$REPO" rev-parse --abbrev-ref HEAD)" "$(git -C "$REPO" rev-parse --short HEAD)"
mkdir -p "$RESULT" "$FIGURE"
{ printf 'field\tvalue\n'
  printf 'master_commit\t%s\n' "$(git -C "$REPO" rev-parse origin/master 2>/dev/null || echo NA)"
  printf 'branch_name\t%s\n' "$(git -C "$REPO" rev-parse --abbrev-ref HEAD)"
  printf 'branch_commit\t%s\n' "$(git -C "$REPO" rev-parse HEAD)"
  printf 'host_cores\t%s\n' "$(nproc)"
  printf 'date\t%s\n' "$(date -Iseconds)"
} > "$RESULT/provenance.tsv"

exit $fail

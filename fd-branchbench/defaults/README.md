# Default sort order and `--confident`

Selects the shipped ranking keys and filter cutoffs offline: one raw, unsorted dump per query
is scored under every candidate ordering, so no candidate needs its own search.

```bash
scripts/collect.sh <folddisco binary>        # raw dumps for the 250 selection queries
scripts/collect_heldout.sh <folddisco binary> # the 492 held-out queries (q755 minus q250)
scripts/analyze.py                            # 281 lexicographic orders, filter grid
scripts/composite.py && scripts/composite2.py # 69 evidence x quality composites
scripts/plot.py                               # fig_sort, fig_confident
```

Paths come from `../scripts/lib.sh` (`FOLDDISCO_BENCH_DATA`, `MCSA_DIR`, `ANS_*`); run
`../scripts/00_setup.sh` first. `NP=<n>` sets the worker count.

Outputs in `result/`: `sort_per_query.tsv.gz`, `composite*_per_query.tsv.gz`,
`filter_per_query.tsv.gz` (per query and candidate), figures in `figure/`. The raw dumps under
`build/` are large and untracked.

Selected: `match_score` = idf x coverage² x TM-score per match, `structure_score` = matched² x
√idf / (1 + RMSD) per structure, `--confident` = coverage ≥ 0.8 with RMSD ≤ 1.0 Å. Every choice
was made on the 250 selection queries and confirmed on the 492 held-out ones; selections that
gained only on the selection set are recorded as such in `folddisco/docs/feature_evaluation.md` §15.

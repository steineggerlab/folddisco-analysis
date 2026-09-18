# Index-time expansion on an M-CSA motif-only database

Binary `../build/bin/folddisco-5ae7321` (feature-integration 5ae7321), 2026-09-17; 5 queries × 4
threads. The first run used 2cd6e60 (8 threads, shared host): Sens@1FP within ±0.01
(`result_2cd6e60/`).

## Method

- **Database**: `scripts/make_motif_db.py` writes one PDB per (M-CSA entry, answer PDB) for the
  250 entries of `queryspec_250.tsv`, holding only the annotated catalytic residues (all atoms,
  first model; first homologue row with no `-` and every residue present). 28,140 answers →
  24,762 motifs, 249 entries (median 28 motifs per entry, max 2,085). Files are named
  `M####__pdbXXXX.pdb`.
- **Indices** (`folddisco index -p build/motifs -i build/index/<I> -t 8`):

| index | flags | build (s) | postings (MB) | offsets (MB) |
| --- | --- | --- | --- | --- |
| I0 | — | 8.8 | 0.8 | 0.9 |
| I1 | `--expand-radius 1` | 9.3 | 3.0 | 2.7 |
| I2 | `--expand-radius 2` | 10.0 | 5.1 | 4.3 |
| I3 | `--aa-subst blosum62` | 10.6 | 7.8 | 7.1 |
| I4 | `--expand-radius 1 --aa-subst blosum62` | 11.7 | 27.0 | 18.9 |

  Peak RSS is 16.8 GB for every build (dense 2^30-entry hash table before pruning).
- **Queries**: whole query proteins took 3–34 s each against I0 (74k–293k hashes), so each query
  is the query-protein site instead: every residue with Cα ≤ 12 Å of a motif Cα
  (`scripts/make_site_queries.py`; 42–429 residues, median 100). *mut* uses the same site on the
  renamed-residue structure from `../build/mutant_250`; *star* adds `:*` to the renamed residue.
  `--top 500` (`scripts/run_queries.sh`; `BIN`, `NP`, `NT` override). On I1–I4 the query
  looks up exact hashes (the default on an expanded index).
- **Scoring** (`scripts/evaluate.py`): a hit is TP if its entry is the query's, neutral if it is
  the query's own PDB. Rankings: *idf* = output order (per-match IDF desc, RMSD asc); *coverage* =
  matched residues / motif size desc, then RMSD. 246 queries have ≥1 non-self motif.

## Results

Mean Sens@1FP (idf / coverage ranking), top-1 accuracy (idf) and seconds per query:

| config | published sites | mutant sites | top-1 pub / mut | s |
| --- | --- | --- | --- | --- |
| I0 | 0.471 / 0.301 | 0.339 / 0.076 | 0.842 / 0.626 | 0.32 |
| I0 `--sensitive` | **0.501** / 0.305 | 0.372 / 0.067 | 0.837 / 0.638 | 0.57 |
| I0 `--aa-subst blosum62` | 0.334 / 0.280 | 0.305 / 0.273 | 0.642 / 0.594 | 3.04 |
| I0, `:*` on renamed residue | — | **0.440 / 0.300** | — / 0.793 | 0.34 |
| I1 | 0.447 / 0.307 | 0.323 / 0.081 | 0.821 / 0.585 | 0.32 |
| I2 | 0.464 / 0.325 | 0.339 / 0.081 | 0.801 / 0.598 | 0.56 |
| I3 | 0.337 / 0.257 | 0.309 / 0.251 | 0.646 / 0.598 | 3.06 |
| I4 | 0.319 / 0.246 | 0.300 / 0.238 | 0.614 / 0.577 | 3.07 |

Paired Δ Sens@1FP, idf ranking [95% bootstrap CI] (win/loss); coverage ranking in brackets:

| comparison | Δ idf | Δ coverage |
| --- | --- | --- |
| I1 − I0 (index r1 vs query r1) | −0.025 [−0.040, −0.011] (23/57) | +0.005 |
| I2 − I0 `--sensitive` | −0.037 [−0.054, −0.023] (24/71) | +0.020 |
| I3 − I0, published | −0.134 [−0.161, −0.110] (17/142) | −0.044 |
| mut: I3 − I0 | −0.031 [−0.053, −0.009] (44/81) | +0.174 |
| mut: I0 `--aa-subst` − I0 | −0.035 [−0.057, −0.012] (40/83) | +0.197 |
| mut: I0 `:*` − I0 | **+0.101** [+0.079, +0.124] (117/10) | +0.224 |

Index-side vs query-side expansion (`result/equivalence.tsv`): 0/250 identical rankings; mean
Jaccard of the returned motifs 0.58 (r1), 0.43 (r2), 0.62 (blosum62). Candidate sets before
`--top` also differ (e.g. M0346: I0 2,639, I1 3,261, shared 2,500).

## Notes

- Query-side expansion perturbs the query feature, index-side the target feature; near a bin
  edge they reach different pairs, so the two are not interchangeable.
- On an expanded index the lookup is exact, so `count_query` scores index-generated
  substitutions as exact hits (weight 1). Per-match IDF uses the matching map and does weight them.
- `--aa-subst` here substitutes all ~100 site residues, not a 3–8 residue motif; the idf ranking
  loses on published sites while coverage ranking recovers mutant sites.
- Index size grows 3.5× (r1), 6× (r2), 9× (blosum62), 32× (r1 + blosum62) on motif-only data.

## Files

`result/`: `index_build.tsv`, `per_query.tsv`, `summary.tsv`, `paired.tsv`, `equivalence.tsv`
(`result_2cd6e60/`: first run);
`figure/index_expansion.{png,pdf}`; `build/` (motifs 128 MB, indices 88 MB, outputs 97 MB).

```bash
# folddisco = ../build/bin/folddisco-5ae7321
../.venv/bin/python scripts/make_motif_db.py $MCSA/queryspec_250.tsv build/motifs build/motif_manifest.tsv
for I in "I0|" "I1|--expand-radius 1" "I2|--expand-radius 2" "I3|--aa-subst blosum62" "I4|--expand-radius 1 --aa-subst blosum62"; do
  folddisco index -p build/motifs -i build/index/${I%%|*} -t 8 ${I#*|}; done
../.venv/bin/python scripts/make_site_queries.py $MCSA/queryspec_250.tsv ../build/mutant_250/spec.tsv build/site_queries.tsv 12
./scripts/run_queries.sh && ../.venv/bin/python scripts/evaluate.py && ../.venv/bin/python scripts/plot.py
```

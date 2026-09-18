# folddisco: `feature-integration` vs `origin/master`

Query runtime, motif precision/recall/F1, and M-CSA Sens@1FP for the branch against upstream
master, plus the amino-acid substitution axis. Plots follow ggpubr `theme_pubr` in seaborn
Run date and commits: `result/provenance.tsv` (branch b907df5, master 2a756d9).

Publication: [Kim H, Kim RS, Mirdita M, Yoon J, Steinegger M. Structural motif search across the
protein-universe with Folddisco. *Nature Biotechnology* (2026)](https://www.nature.com/articles/s41587-026-03162-9).

## Setup

`scripts/lib.sh` reads every path from the environment, so nothing is tied to one machine.
Set them once in `scripts/config.local.sh` (untracked) or export them:

| variable | default | what it points at |
| --- | --- | --- |
| `FOLDDISCO_REPO` | `../folddisco` | folddisco checkout to build and benchmark |
| `FOLDDISCO_BENCH_DATA` | `./data` | root of the inputs below ([Zenodo 10.5281/zenodo.16679607](https://doi.org/10.5281/zenodo.16679607)) |
| `IDX_HSAPIENS` | `$DATA/index/h_sapiens_folddisco` | human AFDB index (23,391) |
| `MCSA_DIR` | `$DATA/mcsa_rebuild` | M-CSA index, `store/`, `queryspec_{60,250,755}.tsv` |
| `ANS_MCSA` | `$DATA/mcsa/mcsa_answers` | M-CSA answer sets |
| `ANS_ZINC`, `ANS_SERINE` | `$DATA/answers/*.tsv` | zinc finger (761), MEROPS S01 (124) |
| `MCSA_HOMOLOGUES`, `MCSA_STORE` | under `$DATA` | inputs for `mcsa_answer_classes.py` |

```bash
./scripts/00_setup.sh          # tools, master worktree, cargo build, venv, input check
SKIP_MASTER=1 ./scripts/00_setup.sh   # branch only, no origin/master comparison
```

It is idempotent and names every missing input, so run it until it exits 0. Indices and answer
sets come from the folddisco dataset record on
[Zenodo](https://doi.org/10.5281/zenodo.16679607); point `FOLDDISCO_BENCH_DATA` at the
unpacked copy, or override the individual paths.

## Run

```bash
REPEATS=5 WARMUP=1 ./scripts/run_all.sh 250   # ~70 min on 20 cores
```

Steps: `00_setup.sh` → `01_run_motif.sh` → `02_run_mcsa.sh N` → `04_run_mutant.sh N` →
`05_run_confident.sh N` → answer classes → `03_plot.py`. `CONFIGS=core` limits to the first
three configs; `RESUME=1` resumes the motif runner; M-CSA runners resume via `.ok` files.
The ranking/filter selection (`defaults/`) and the motif-only index study (`index_expansion/`)
have their own scripts and READMEs.

## Configs

| name | flags |
| --- | --- |
| `master` | origin/master binary |
| `default` | branch, no flags |
| `sensitive` | `--sensitive` (`--expand-radius 2`) |
| `aa_blosum62` / `aa_group` / `aa_size` | `--aa-subst <scheme>` |
| `sensitive_aa` | `--sensitive --aa-subst blosum62` |

Mutant configs (`04_run_mutant.sh`): `orig` (published query), `mut` (one residue renamed),
`mut_star` (`:*` on it), `mut_blosum62`, `mut_sens_star` (`:*` + `--sensitive`).

## Data

- **Motif**: human AFDB index (23,391), 8 commands A1–D2 in `spec/motif_queries.tsv`; zinc
  finger (761) and MEROPS S01 (124) answer sets.
- **M-CSA**: 62,122-entry PDB index; nested subsets q60 ⊂ q250 ⊂ q755
  (`bench_data/mcsa_rebuild/queryspec_*.tsv`). q250 ≈ 18 min for the three core configs.
- **Mutant M-CSA**: `make_mutant_queries.py` renames the first query residue with a positive
  BLOSUM62 alternative to its best one (His→Tyr 60, Asp→Glu 39, Glu→Asp 30, Arg→Lys 26, …);
  coordinates and answer sets are unchanged.
- **Answer classes**: `mcsa_answer_classes.py` compares residue names at the aligned catalytic
  positions (M-CSA homologue table): q250 has 12,702 exact answers, 1,327 with one conservative
  (BLOSUM62 > 0) and 2,652 with one other substitution.

## Results (q250)

M-CSA, mean Sens@1FP and seconds per query (5 proc × 4 threads):

| config | Sens@1FP | Δ master (win/loss) | zero-TP | runtime |
| --- | --- | --- | --- | --- |
| master | 0.4475 | — | 10 | 5.35 |
| default | 0.4961 | +0.0486 (133/42) | 1 | 2.65 |
| sensitive | **0.5015** | +0.0540 (134/42) | 1 | 3.02 |
| aa_blosum62 | 0.4764 | +0.0289 (124/56) | 1 | 5.69 |
| aa_group | 0.4734 | +0.0259 (120/60) | 1 | 6.54 |
| aa_size | 0.4662 | +0.0188 (105/67) | 3 | 6.79 |
| sensitive_aa | 0.4807 | +0.0333 (120/61) | 2 | 5.55 |

Mutant M-CSA (renaming one residue costs 0.4961 → 0.3020):

| config | Sens@1FP | Δ vs `mut` [95% CI] (win/loss) | runtime |
| --- | --- | --- | --- |
| `mut_star` | 0.4794 | +0.177 [+0.148, +0.209] (172/22) | 3.37 |
| `mut_blosum62` | 0.4546 | +0.153 [+0.125, +0.183] (158/36) | 5.76 |
| `mut_sens_star` | **0.4830** | +0.181 [+0.150, +0.212] (176/23) | 3.83 |

Answers ranked before the first FP on exact M-CSA, by class (`result/mcsa_classes_250.tsv`):

| config | exact | conservative single | other single |
| --- | --- | --- | --- |
| master | 5,749 | 458 | 841 |
| default | 6,378 | 388 | 997 |
| sensitive | **6,413** | 391 | 887 |
| aa_blosum62 | 5,993 | **616** | 646 |
| sensitive_aa | 6,021 | 615 | 596 |

Motif (`result/motif_metrics.tsv`): default is byte-identical to master; `--sensitive` F1
+0.011 to +0.033 on A2/B1/B2/C1, −0.004 to −0.009 on A1/D1/D2. With substitution, Sens@1FP
equals default on the serine queries (C1/C2 0.2339; aa5542d had 0.0000) and D1/D2 F1 is
0.66–0.84 (aa5542d: 0.005–0.34). Set F1 on the prefilter commands still drops (C1 0.883 → 0.370
for blosum62): more structures pass `--covered-node`, and the answer sets hold exact residues only.

## Scoring changes

- **Substituted hashes** (78bb4a6): weight 0.75 per substituted residue, IDF capped at the query
  pair's, one substituted hit per query edge. At aa5542d they scored full IDF, so blosum62 fell
  to 0.3537 on q60 (default 0.4680) at 37 s per query. The 4096-hash cap never binds at default
  tolerances and was not the cause. Selection: `selection/scoring_sweep.txt`.
- **Neighbouring bins** (5ae7321): per-match IDF counts edges from expanded bins, like substituted
  ones, only between matched residues. Applying the per-edge cap to them too lowered zinc-finger
  Sens@50FP from 0.89 to 0.37, so candidate scoring still sums them. Selection:
  `selection/geometric_sweep.txt`.
- A code-based residue prefilter halves default M-CSA runtime (5.78 → 2.67 s), output unchanged.

Details: `folddisco/docs/feature_evaluation.md` §14. Index-time expansion: `index_expansion/`.

## Identity control

Motif commands: default vs master 8/8 byte-identical (they rank by candidate score). M-CSA
per-match output differs by design since 5ae7321 (59/250 rank-identical); at 78bb4a6 it was
249/250, the exception being M0093, where master returns nothing for the valid pair `B559,B866`.

## Outputs

`result/`: `motif_metrics.tsv`, `mcsa_fp1_250.tsv`, `mcsa_summary_250.tsv`,
`mutant_fp1_250.tsv`, `mutant_summary_250.tsv`, `mcsa_classes_250.tsv`, identity TSVs,
`provenance.tsv`. `figure/`: fig1 motif accuracy, fig2 runtime, fig3 M-CSA Sens@1FP, fig4 gain
vs master, fig5 substitution, fig6 mutant queries (PNG + PDF). Earlier results:
`../archive/fd-branchbench_aa5542d/`, `../archive/fd-branchbench_78bb4a6/`.

## Caveats

- Sens@1FP is a k=1 metric; single accessions move it. Read paired win/loss with the means.
- M-CSA timings are under 5-way parallelism; motif timings are serial medians of 5.
- Mutant queries rename a residue without changing its side chain; hashes use Cα/Cβ only.

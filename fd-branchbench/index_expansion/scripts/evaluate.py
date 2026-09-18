"""Score motif-index searches.

A hit (motif file M####__pdbXXXX.pdb) is TP if its entry equals the query's, neutral if it is
the query's own PDB, FP otherwise. Two rankings per query:
  idf       - output order (idf desc, rmsd asc), first row per motif
  coverage  - matched residues / motif size desc, then rmsd asc
Writes result/per_query.tsv, summary.tsv, paired.tsv, equivalence.tsv.
"""
import csv, os, re, sys
from collections import defaultdict
import numpy as np
import pandas as pd

D = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BUILD, RESULT = f"{D}/build", f"{D}/result"
NAME = re.compile(r"(M\d{4})__pdb(\w+)\.pdb$")

manifest = pd.read_csv(f"{BUILD}/motif_manifest.tsv", sep="\t", dtype=str)
size = {r.file: len(r.residues.split(",")) for r in manifest.itertuples()}
per_entry = manifest.groupby("mcsa").pdb.apply(set).to_dict()
queries = [l.rstrip("\n").split("\t") for l in open(f"{BUILD}/site_queries.tsv")]


def read_hits(path):
    best = {}
    order = []
    with open(path) as f:
        for row in csv.reader(f, delimiter="\t"):
            if len(row) < 4:
                continue
            name = os.path.basename(row[0])
            nodes, idf, rmsd = int(row[1]), float(row[2]), float(row[3])
            cov = min(nodes, size.get(name, nodes)) / size.get(name, max(nodes, 1))
            if name not in best:
                order.append(name)
                best[name] = (cov, rmsd)
            elif (cov, -rmsd) > (best[name][0], -best[name][1]):
                best[name] = (cov, rmsd)
    by_cov = sorted(best, key=lambda n: (-best[n][0], best[n][1]))
    return order, by_cov


def score(ranked, entry, self_pdb, n_true):
    labels = []
    for name in ranked:
        m = NAME.search(name)
        if not m or m.group(2) == self_pdb:
            continue
        labels.append(m.group(1) == entry)
    tp_before = 0
    for is_tp in labels:
        if not is_tp:
            break
        tp_before += 1
    return dict(top1=float(bool(labels) and labels[0]), sens1fp=tp_before / n_true if n_true else np.nan,
                any_tp=float(any(labels)), hits=len(labels))


rows = []
configs = sorted(os.listdir(f"{BUILD}/out"))
for cfg in configs:
    for q in queries:
        qid, pdb = q[0], q[1]
        path = f"{BUILD}/out/{cfg}/{qid}.tsv"
        if not os.path.exists(f"{BUILD}/out/{cfg}/{qid}.ok"):
            continue
        self_pdb = os.path.basename(pdb).removesuffix(".pdb").lower()
        n_true = len(per_entry.get(qid, set()) - {self_pdb})
        rt = float(open(f"{BUILD}/out/{cfg}/{qid}.time").read().split()[0])
        by_idf, by_cov = read_hits(path)
        for rank_name, ranked in (("idf", by_idf), ("coverage", by_cov)):
            rows.append(dict(config=cfg, query=qid, ranking=rank_name, n_true=n_true, runtime=rt,
                             **score(ranked, qid, self_pdb, n_true)))
df = pd.DataFrame(rows)
df = df[df.n_true > 0]
df.to_csv(f"{RESULT}/per_query.tsv", sep="\t", index=False, float_format="%.4f")

summary = (df.groupby(["config", "ranking"])
             .agg(n=("query", "nunique"), top1=("top1", "mean"), sens1fp=("sens1fp", "mean"),
                  any_tp=("any_tp", "mean"), hits=("hits", "mean"), runtime=("runtime", "mean"))
             .reset_index())
summary.to_csv(f"{RESULT}/summary.tsv", sep="\t", index=False, float_format="%.4f")

PAIRS = [("orig_I1", "orig_I0"), ("orig_I2", "orig_I0_sens"), ("orig_I0_sens", "orig_I0"),
         ("orig_I0_aa", "orig_I0"), ("orig_I3", "orig_I0"), ("orig_I4", "orig_I0"),
         ("mut_I0", "orig_I0"),
         ("mut_I0_star", "mut_I0"), ("mut_I0_aa", "mut_I0"), ("mut_I3", "mut_I0"), ("mut_I4", "mut_I0"),
         ("mut_I1", "mut_I0"), ("mut_I2", "mut_I0_sens")]
rng = np.random.default_rng(0)
paired = []
for a, b in PAIRS:
    for rank_name in ("idf", "coverage"):
        x = df[(df.config == a) & (df.ranking == rank_name)].set_index("query")
        y = df[(df.config == b) & (df.ranking == rank_name)].set_index("query")
        j = x.join(y, lsuffix="_a", rsuffix="_b", how="inner")
        if j.empty:
            continue
        d = (j.sens1fp_a - j.sens1fp_b).to_numpy()
        bs = [rng.choice(d, len(d)).mean() for _ in range(2000)]
        paired.append(dict(a=a, b=b, ranking=rank_name, n=len(j), d_sens1fp=d.mean(),
                           ci_lo=np.percentile(bs, 2.5), ci_hi=np.percentile(bs, 97.5),
                           win=int((d > 0).sum()), loss=int((d < 0).sum()),
                           d_top1=(j.top1_a - j.top1_b).mean(), runtime_ratio=j.runtime_a.mean() / j.runtime_b.mean()))
pd.DataFrame(paired).to_csv(f"{RESULT}/paired.tsv", sep="\t", index=False, float_format="%.4f")

# Index-side vs query-side expansion: same hits?
eq = []
for a, b in [("orig_I1", "orig_I0"), ("orig_I2", "orig_I0_sens"), ("mut_I1", "mut_I0"), ("mut_I2", "mut_I0_sens"),
             ("orig_I3", "orig_I0_aa"), ("mut_I3", "mut_I0_aa")]:
    same_order = jac = n = 0
    for q in queries:
        pa, pb = f"{BUILD}/out/{a}/{q[0]}.tsv", f"{BUILD}/out/{b}/{q[0]}.tsv"
        if not (os.path.exists(pa) and os.path.exists(pb)):
            continue
        ra, rb = read_hits(pa)[0], read_hits(pb)[0]
        n += 1
        same_order += ra == rb
        sa, sb = set(ra), set(rb)
        jac += len(sa & sb) / max(len(sa | sb), 1)
    if n:
        eq.append(dict(index_side=a, query_side=b, n=n, identical_ranking=same_order / n, mean_jaccard=jac / n))
pd.DataFrame(eq).to_csv(f"{RESULT}/equivalence.tsv", sep="\t", index=False, float_format="%.4f")
print(summary.to_string(index=False))

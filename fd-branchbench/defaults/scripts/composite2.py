"""Prior-selected composite rankings: evidence x quality, over every implemented metric.

Evidence is IDF scaled by the matched fraction of the query (coverage^p); quality is a
factor in (0, 1] built from one geometric metric. RMSD breaks ties. Structure-level
candidates use the same shape over max_node_cov / min_rmsd / min_drmsd.

Usage: composite2.py [raw_dir] [spec]      (spec: the M-CSA queryspec the raw dir was run on)
Writes result/composite2_<tag>.tsv.gz and prints the ranking.
"""
import os, sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze as A

RAW = sys.argv[1] if len(sys.argv) > 1 else f"{A.HERE}/build/raw"
if len(sys.argv) > 2:
    A.SPEC = sys.argv[2]
TAG = os.path.basename(RAW)
POWERS = (0.5, 1.0, 1.5, 2.0)


def quality(df):
    """name -> factor in (0, 1], higher is better."""
    g = {"none": None}
    for t in (0.5, 1.0, 2.0):
        g[f"rmsd/{t}"] = 1.0 / (1.0 + df.rmsd.to_numpy(float) / t)
    for t in (0.5, 1.0):
        g[f"drmsd/{t}"] = 1.0 / (1.0 + df.drmsd.to_numpy(float) / t)
    g["tm_score"] = df.tm_score.to_numpy(float)
    g["gdt_ha"] = df.gdt_ha.to_numpy(float)
    g["gdt_ts"] = df.gdt_ts.to_numpy(float)
    g["chamfer"] = 1.0 / (1.0 + df.chamfer_distance.to_numpy(float))
    g["hausdorff"] = 1.0 / (1.0 + df.hausdorff_distance.to_numpy(float))
    g["max_dist_dev/2"] = 1.0 / (1.0 + df.max_dist_deviation.to_numpy(float) / 2.0)
    return g


def match_orders(df, qlen):
    n, idf, rmsd = (df[c].to_numpy(float) for c in ("node_count", "idf", "rmsd"))
    cov = n / qlen
    out = {"idf,rmsd": np.lexsort((rmsd, -idf))}
    for p in POWERS:
        base = idf * cov ** p
        for name, q in quality(df).items():
            score = base if q is None else base * q
            out[f"idf*cov^{p}" + ("" if q is None else f"*{name}") + ",rmsd"] = np.lexsort((rmsd, -score))
    return out


def structure_orders(df, qlen):
    idf, mx, rmsd, drmsd = (df[c].to_numpy(float)
                            for c in ("idf", "max_node_count", "min_rmsd", "min_drmsd"))
    cov = mx / qlen
    out = {"idf,min_rmsd": np.lexsort((rmsd, -idf)),
           "max_node_count,min_rmsd": np.lexsort((rmsd, -mx))}
    for p in (1.0, 2.0):
        for q in (0.0, 0.5, 1.0):
            base = cov ** p * idf ** q
            for gname, g in (("none", None), ("rmsd/1", 1 / (1 + rmsd)), ("drmsd/1", 1 / (1 + drmsd))):
                score = base if g is None else base * g
                out[f"cov^{p}*idf^{q}" + ("" if g is None else f"*{gname}") + ",min_rmsd"] = \
                    np.lexsort((rmsd, -score))
    return out


def one(task):
    run, qid, meta = task
    answers, qlen, afdb, trunc = meta
    path = f"{RAW}/{run}/{qid}.tsv"
    if not os.path.exists(path.replace(".tsv", ".ok")):
        return []
    per_match = run.split("_")[1] == "pm"
    df = pd.read_csv(path, sep="\t", header=None, names=A.PM if per_match else A.PS)
    if df.empty:
        return []
    ids = np.array([A.target_id(t, afdb) for t in df.tid], dtype=object)
    orders = match_orders(df, qlen) if per_match else structure_orders(df, qlen)
    return [(run, qid, name) + A.rank_metrics(A.ranked_ids(ids, idx, trunc if per_match else None), answers)
            for name, idx in orders.items()]


if __name__ == "__main__":
    motif, mcsa = A.queries()
    runs = [d for d in sorted(os.listdir(RAW)) if not d.endswith("skip")]
    tasks = [(r, q, m) for r in runs for q, m in (motif if r.startswith("mo_") else mcsa).items()]
    rows = []
    with ProcessPoolExecutor(int(os.environ.get("NP", "16"))) as pool:
        for r in pool.map(one, tasks, chunksize=4):
            rows += r
    d = pd.DataFrame(rows, columns=["run", "query", "sort", "sens1", "sens5", "ap"])
    d.to_csv(f"{A.HERE}/result/composite2_{TAG}.tsv.gz", sep="\t", index=False)
    for kind in ("pm", "ps"):
        sub = d[d.run.str.contains(f"_{kind}_")]
        if sub.empty:
            continue
        g = sub.groupby(["run", "sort"])[["sens1", "ap"]].mean().reset_index()
        piv = g.pivot(index="sort", columns="run")
        mc = [c for c in piv["sens1"].columns if c.startswith("mc_")]
        piv[("sens1", "mcsa")] = piv["sens1"][mc].mean(axis=1)
        piv[("ap", "mcsa")] = piv["ap"][mc].mean(axis=1)
        print(f"\n=== {kind}: top 12 by M-CSA mean AP (of {piv.shape[0]} candidates)")
        print(piv.sort_values(("ap", "mcsa"), ascending=False).head(12).round(4).to_string())

"""Composite per-match orderings that plain --sort-by keys cannot express."""
import os, sys, numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze as A

RUNS = ["mc_pm_def", "mc_pm_sens", "mc_pm_star", "mc_pm_blosum", "mo_pm_def", "mo_pm_sens"]
RAW = sys.argv[1] if len(sys.argv) > 1 else f"{os.path.dirname(A.HERE)}/defaults/build/raw"


def orders(df, qlen):
    n, idf, rmsd = (df[c].to_numpy(dtype=float) for c in ("node_count", "idf", "rmsd"))
    cov = n / qlen
    out = {"idf,rmsd": np.lexsort((rmsd, -idf))}
    for p in (0.5, 1.0, 2.0, 3.0):
        out[f"idf*cov^{p},rmsd"] = np.lexsort((rmsd, -(idf * cov ** p)))
    out["full_first,idf,rmsd"] = np.lexsort((rmsd, -idf, -(cov >= 0.999).astype(float)))
    out["cov>=0.8_first,idf,rmsd"] = np.lexsort((rmsd, -idf, -(cov >= 0.8).astype(float)))
    out["idf/n,rmsd"] = np.lexsort((rmsd, -(idf / np.maximum(n, 1))))
    out["node_count,idf,rmsd"] = np.lexsort((rmsd, -idf, -n))
    return out


def one(task):
    run, qid, meta = task
    answers, qlen, afdb, trunc = meta
    p = f"{RAW}/{run}/{qid}.tsv"
    if not os.path.exists(p.replace(".tsv", ".ok")):
        return []
    df = pd.read_csv(p, sep="\t", header=None, names=A.PM)
    ids = np.array([A.target_id(t, afdb) for t in df.tid], dtype=object)
    rows = []
    for name, idx in orders(df, qlen).items():
        s1, s5, ap = A.rank_metrics(A.ranked_ids(ids, idx, trunc), answers)
        rows.append((run, qid, name, s1, s5, ap))
    return rows


if __name__ == "__main__":
    motif, mcsa = A.queries()
    tasks = [(r, q, m) for r in RUNS for q, m in (motif if r.startswith("mo_") else mcsa).items()]
    rows = []
    with ProcessPoolExecutor(int(os.environ.get("NP", "10"))) as pool:
        for r in pool.map(one, tasks, chunksize=4):
            rows += r
    d = pd.DataFrame(rows, columns=["run", "query", "sort", "sens1", "sens5", "ap"])
    d.to_csv(f"{A.HERE}/result/composite_per_query.tsv.gz", sep="\t", index=False)
    g = d.groupby(["run", "sort"])[["sens1", "sens5", "ap"]].mean().reset_index()
    for m in ("sens1", "ap"):
        piv = g.pivot(index="sort", columns="run", values=m)
        piv["mcsa"] = piv[[c for c in piv.columns if c.startswith("mc_")]].mean(axis=1)
        print(m); print(piv.sort_values("mcsa", ascending=False).round(4).to_string(), "\n")

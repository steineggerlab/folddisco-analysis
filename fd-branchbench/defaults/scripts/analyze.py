"""Offline selection of the default sort order and the --confident filter.

Reads build/raw/<run>/<query>.tsv (unsorted, untruncated; see collect.sh), re-sorts and
re-filters them the way the binary does, and writes result/sort_*.tsv and result/filter_*.tsv.
Usage: analyze.py [raw_dir] [result_dir] [queryspec]
"""
import itertools, os, sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = sys.argv[1] if len(sys.argv) > 1 else f"{HERE}/build/raw"
OUT = sys.argv[2] if len(sys.argv) > 2 else f"{HERE}/result"
# Inputs follow scripts/lib.sh: the environment wins, then the default layout under DATA.
DATA = os.environ.get("FOLDDISCO_BENCH_DATA", f"{os.path.dirname(HERE)}/data")
MCSA_DIR = os.environ.get("MCSA_DIR", f"{DATA}/mcsa_rebuild")
SPEC = sys.argv[3] if len(sys.argv) > 3 else f"{MCSA_DIR}/queryspec_250.tsv"
ANS_MCSA = os.environ.get("ANS_MCSA", f"{DATA}/mcsa/mcsa_answers")
ANS_ZINC = os.environ.get("ANS_ZINC", f"{DATA}/answers/zincfinger_answer.tsv")
ANS_SERINE = os.environ.get("ANS_SERINE", f"{DATA}/answers/serinepeptidase_answer.tsv")
TOP = 6000  # M-CSA runs use --top 6000; per-match output is cut after sorting

PM = ["tid", "node_count", "idf", "rmsd", "e_value", "tm_score", "gdt_ts", "gdt_ha",
      "chamfer_distance", "hausdorff_distance", "drmsd", "max_dist_deviation"]
PS = ["tid", "idf", "total_match_count", "node_count", "edge_count", "max_node_count",
      "min_rmsd", "min_drmsd", "nres", "plddt"]
DESC = {"node_count", "idf", "tm_score", "gdt_ts", "gdt_ha", "max_node_count",
        "total_match_count", "edge_count", "nres", "plddt"}
MATCH_KEYS = ["node_count", "idf", "rmsd", "drmsd", "tm_score", "gdt_ts", "gdt_ha",
              "chamfer_distance", "hausdorff_distance", "max_dist_deviation"]
STRUCT_KEYS = ["max_node_count", "node_count", "idf", "min_rmsd", "min_drmsd",
               "total_match_count", "edge_count"]
MOTIF = {"A1": ("zinc", 4), "B1": ("zinc", 3), "C1": ("serine", 3), "D1": ("zinc", 23)}


def target_id(path, afdb):
    name = os.path.basename(path)
    for ext in (".pdb.gz", ".cif.gz", ".fcz.gz", ".ent.gz", ".pdb", ".cif", ".fcz", ".ent"):
        if name.endswith(ext):
            name = name[:-len(ext)]
            break
    if afdb:
        parts = name.split("-")
        return parts[1] if len(parts) >= 2 else name
    return name


def first_column(path):
    return {l.split("\t")[0].strip() for l in open(path) if l.strip()}


def queries():
    """(run kind, query id) -> (answers, query length, afdb ids, truncate)"""
    zinc, serine = first_column(ANS_ZINC), first_column(ANS_SERINE)
    motif = {q: ({"zinc": zinc, "serine": serine}[a], n, True, None) for q, (a, n) in MOTIF.items()}
    mcsa = {}
    for line in open(SPEC):
        mid, _, res = line.rstrip("\n").split("\t")[:3]
        mcsa[mid] = (first_column(f"{ANS_MCSA}/{mid}_answer.tsv"), len(res.split(",")), False, TOP)
    return motif, mcsa


def sort_orders(keys, per_match):
    """Lexicographic key lists: 1-2 keys, and 3 keys when the first is a count."""
    counts = {"node_count", "max_node_count"}
    out = [[k] for k in keys] + [list(p) for p in itertools.permutations(keys, 2)]
    out += [list(p) for p in itertools.permutations(keys, 3) if p[0] in counts]
    return out


def order_by(df, keys):
    cols = []
    for k in reversed(keys):
        v = df[k].to_numpy(dtype=float)
        cols.append(-v if k in DESC else v)
    return np.lexsort(cols)  # stable, like rayon's par_sort_by


def ranked_ids(ids, idx, truncate):
    if truncate is not None:
        idx = idx[:truncate]
    seen, order = set(), []
    for t in ids[idx]:
        if t not in seen:
            seen.add(t)
            order.append(t)
    return order


def rank_metrics(order, answers):
    tp = fp = 0
    s1 = s5 = None
    ap = 0.0
    for t in order:
        if t in answers:
            tp += 1
            ap += tp / (tp + fp)
        else:
            fp += 1
            if fp == 1:
                s1 = tp
            if fp == 5:
                s5 = tp
                break
    n = len(answers)
    s1 = tp if s1 is None else s1
    s5 = tp if s5 is None else s5
    # AP over the whole list
    tp = fp = 0
    ap = 0.0
    for t in order:
        if t in answers:
            tp += 1
            ap += tp / (tp + fp)
        else:
            fp += 1
    return s1 / n, s5 / n, ap / n


FILTERS = [(cov, r, d) for cov in (1.0, 0.8, 0.0)
           for r in (0.5, 0.75, 1.0, 1.25, 1.5, 2.0, np.inf)
           for d in (0.25, 0.5, 0.75, 1.0, 1.5, np.inf)]


def evaluate(task):
    run, qid, meta = task
    answers, qlen, afdb, truncate = meta
    per_match = run.split("_")[1] == "pm"
    path = f"{RAW}/{run}/{qid}.tsv"
    if not os.path.exists(path.replace(".tsv", ".ok")):
        return [], []
    try:
        df = pd.read_csv(path, sep="\t", header=None, names=PM if per_match else PS)
    except pd.errors.EmptyDataError:
        df = pd.DataFrame(columns=PM if per_match else PS)
    ids = np.array([target_id(t, afdb) for t in df.tid], dtype=object)
    sorts = []
    skip = run.endswith("skip")
    if not skip:
        keys = MATCH_KEYS if per_match else STRUCT_KEYS
        for ks in sort_orders(keys, per_match):
            s1, s5, ap = rank_metrics(ranked_ids(ids, order_by(df, ks), truncate if per_match else None), answers)
            sorts.append((run, qid, ",".join(ks), s1, s5, ap))
    filters = []
    base = ["idf", "rmsd"] if per_match else ["idf", "min_rmsd"]
    if per_match:
        n, rmsd, drmsd = (df[c].to_numpy(dtype=float) for c in ("node_count", "rmsd", "drmsd"))
    elif skip:
        n, rmsd, drmsd = df.node_count.to_numpy(dtype=float), np.zeros(len(df)), np.zeros(len(df))
    else:
        n, rmsd, drmsd = (df[c].to_numpy(dtype=float) for c in ("max_node_count", "min_rmsd", "min_drmsd"))
    order = order_by(df, base)
    for cov, r, d in FILTERS:
        if skip and (r != np.inf or d != np.inf):
            continue
        mask = (n / qlen >= cov - 1e-9) & (rmsd <= r) & (drmsd <= d)
        idx = order[mask[order]]
        passed = ranked_ids(ids, idx, truncate if per_match else None)
        tp = sum(t in answers for t in passed)
        filters.append((run, qid, qlen, cov, r, d, len(passed), tp, len(answers)))
    return sorts, filters


def main():
    os.makedirs(OUT, exist_ok=True)
    motif, mcsa = queries()
    tasks = []
    for run in sorted(os.listdir(RAW)):
        qs = motif if run.startswith("mo_") else mcsa
        tasks += [(run, q, m) for q, m in qs.items()]
    sorts, filters = [], []
    with ProcessPoolExecutor(int(os.environ.get("NP", "16"))) as pool:
        for s, f in pool.map(evaluate, tasks, chunksize=4):
            sorts += s
            filters += f
    pd.DataFrame(sorts, columns=["run", "query", "sort", "sens1", "sens5", "ap"]).to_csv(
        f"{OUT}/sort_per_query.tsv.gz", sep="\t", index=False)
    pd.DataFrame(filters, columns=["run", "query", "qlen", "cov", "rmsd", "drmsd", "hits", "tp", "answers"]).to_csv(
        f"{OUT}/filter_per_query.tsv.gz", sep="\t", index=False)
    print("wrote", OUT, len(sorts), len(filters))


if __name__ == "__main__":
    main()

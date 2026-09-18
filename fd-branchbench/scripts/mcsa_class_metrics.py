"""Per-class M-CSA metrics: how many exact / single-substitution answers each config ranks
before its first false positive, and anywhere in the list.

Usage: mcsa_class_metrics.py <classes.tsv> <result_dir> <out.tsv> [config ...]
result_dir holds <config>/<query>.tsv.
"""
import sys, os
from collections import defaultdict

AA3 = "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split()
# BLOSUM62 positive off-diagonal pairs
CONSERVATIVE = {frozenset(p) for p in [
    ("ALA", "SER"), ("ARG", "GLN"), ("ARG", "LYS"), ("ASN", "ASP"), ("ASN", "HIS"), ("ASN", "SER"),
    ("ASP", "GLU"), ("GLN", "GLU"), ("GLN", "LYS"), ("GLU", "LYS"), ("HIS", "TYR"), ("ILE", "LEU"),
    ("ILE", "MET"), ("ILE", "VAL"), ("LEU", "MET"), ("LEU", "VAL"), ("MET", "VAL"), ("PHE", "TRP"),
    ("PHE", "TYR"), ("SER", "THR"), ("TRP", "TYR")]}


def answer_classes(path):
    """query -> answer -> class in {exact, sub1_cons, sub1_other, other}"""
    out = defaultdict(dict)
    with open(path) as f:
        next(f)
        for line in f:
            q, a, cls, n, detail = line.rstrip("\n").split("\t")
            if cls == "exact":
                c = "exact"
            elif cls == "substituted" and n == "1":
                src, dst = detail.split(">")
                if src in AA3 and dst in AA3:
                    c = "sub1_cons" if frozenset((src, dst)) in CONSERVATIVE else "sub1_other"
                else:
                    c = "other"
            else:
                c = "other"
            # an answer listed under several rows keeps its best class
            rank = ["exact", "sub1_cons", "sub1_other", "other"]
            if a not in out[q] or rank.index(c) < rank.index(out[q][a]):
                out[q][a] = c
    return out


def target_id(path):
    """Same id `folddisco benchmark` derives: basename without a structure extension."""
    name = os.path.basename(path)
    for ext in (".pdb.gz", ".cif.gz", ".fcz.gz", ".ent.gz", ".pdb", ".cif", ".fcz", ".ent"):
        if name.endswith(ext):
            return name[:-len(ext)]
    return name


def ranked_ids(path):
    seen, order = set(), []
    with open(path) as f:
        for line in f:
            tid = target_id(line.split("\t", 1)[0])
            if tid not in seen:
                seen.add(tid)
                order.append(tid)
    return order


CLASSES = ["exact", "sub1_cons", "sub1_other"]

if __name__ == "__main__":
    classes = answer_classes(sys.argv[1])
    result_dir, out_path = sys.argv[2], sys.argv[3]
    configs = sys.argv[4:] or sorted(os.listdir(result_dir))
    with open(out_path, "w") as out:
        cols = ["config", "query", "answers", "sens_at_1fp"]
        for c in CLASSES:
            cols += [f"n_{c}", f"fp1_{c}", f"found_{c}"]
        out.write("\t".join(cols) + "\n")
        for cfg in configs:
            d = os.path.join(result_dir, cfg)
            if not os.path.isdir(d):
                continue
            for q, answers in sorted(classes.items()):
                p = os.path.join(d, f"{q}.tsv")
                if not os.path.exists(p):
                    continue
                order = ranked_ids(p)
                before_fp = set()
                for tid in order:
                    if tid not in answers:
                        break
                    before_fp.add(tid)
                found = set(order) & set(answers)
                row = [cfg, q, len(answers), f"{len(before_fp) / len(answers):.4f}"]
                for c in CLASSES:
                    members = {a for a, k in answers.items() if k == c}
                    row += [len(members), len(members & before_fp), len(members & found)]
                out.write("\t".join(map(str, row)) + "\n")

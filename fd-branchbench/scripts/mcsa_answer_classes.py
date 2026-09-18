"""Split each M-CSA answer set by residue identity at the aligned catalytic positions.

class per (query, answer): exact (some chain row maps every query residue to the same
amino acid), substituted (fully mapped, >=1 differs), partial (no fully mapped row).
Usage: mcsa_answer_classes.py <queryspec.tsv> <out.tsv>
"""
import sys, os, csv
from multiprocessing import Pool
from collections import defaultdict

# Inputs follow scripts/lib.sh: the environment wins, then the default layout under DATA.
DATA = os.environ.get("FOLDDISCO_BENCH_DATA",
                      f"{os.path.dirname(os.path.dirname(os.path.abspath(__file__)))}/data")
MCSA_DIR = os.environ.get("MCSA_DIR", f"{DATA}/mcsa_rebuild")
HOMOLOGUES = os.environ.get("MCSA_HOMOLOGUES",
                            f"{DATA}/mcsa/metadata/catalytic_residues_homologues_parsed.tsv")
STORE = os.environ.get("MCSA_STORE", f"{MCSA_DIR}/store")
ANSWERS = os.environ.get("ANS_MCSA", f"{DATA}/mcsa/mcsa_answers")


def residue_names(path):
    names = {}
    try:
        with open(path) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "CA":
                    key = f"{line[21].strip()}{line[22:26].strip()}"
                    names.setdefault(key, line[17:20].strip())
                elif line.startswith("ENDMDL"):
                    break
    except FileNotFoundError:
        return None
    return names


def rows_by_entry():
    rows = defaultdict(lambda: defaultdict(list))
    with open(HOMOLOGUES) as f:
        next(f)
        for mcsa, pdb, residues, _ in csv.reader(f, delimiter="\t"):
            rows[int(mcsa)][pdb.lower()].append(residues.split(","))
    return rows


ROWS = rows_by_entry()


def classify(spec):
    mid, pdb_path, query = spec
    mcsa = int(mid[1:])
    qpdb = os.path.basename(pdb_path).removesuffix(".pdb").lower()
    qres = query.split(",")
    qnames = residue_names(pdb_path)
    # position of each query residue inside the query structure's own rows
    pos = []
    for r in qres:
        p = next((row.index(r) for row in ROWS[mcsa][qpdb] if r in row), None)
        pos.append(p)
    out = []
    with open(f"{ANSWERS}/{mid}_answer.tsv") as f:
        answers = [l.strip() for l in f if l.strip()]
    for ans in answers:
        pdb = ans.removeprefix("pdb")
        names = residue_names(f"{STORE}/pdb{pdb}.pdb")
        best, cls, detail = -1, "partial", ""
        for row in ROWS[mcsa].get(pdb, []):
            mapped = [row[p] if p is not None and p < len(row) else "-" for p in pos]
            if names is None or "-" in mapped or any(m not in names for m in mapped):
                continue
            same = sum(names[m] == qnames.get(q) for m, q in zip(mapped, qres))
            if same > best:
                best = same
                cls = "exact" if same == len(qres) else "substituted"
                detail = ",".join(f"{qnames.get(q)}>{names[m]}" for m, q in zip(mapped, qres)
                                  if names[m] != qnames.get(q))
        out.append((mid, ans, cls, len(qres) - best if best >= 0 else "NA", detail))
    return out


if __name__ == "__main__":
    specs = [l.rstrip("\n").split("\t")[:3] for l in open(sys.argv[1])]
    with Pool(16) as pool, open(sys.argv[2], "w") as out:
        out.write("query\tanswer\tclass\tn_substituted\tsubstitutions\n")
        for rows in pool.imap(classify, specs):
            for row in rows:
                out.write("\t".join(map(str, row)) + "\n")

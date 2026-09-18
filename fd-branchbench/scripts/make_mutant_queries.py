"""Mutate one residue per M-CSA query to its closest BLOSUM62 alternative.

The first query residue with a positive-scoring alternative is renamed (coordinates kept);
the answer set is unchanged, so exact search loses what substitution should recover.
Usage: make_mutant_queries.py <queryspec.tsv> <out_dir> <out_spec.tsv>
out_spec columns: id, mutant pdb, query, residue, from>to
"""
import sys, os

AA = "ARNDCQEGHILKMFPSTWYV"
THREE = "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split()
# BLOSUM62 positive off-diagonal scores
POS = {("A","S"):1,("R","Q"):1,("R","K"):2,("N","D"):1,("N","H"):1,("N","S"):1,("D","E"):2,
       ("Q","E"):2,("Q","K"):1,("E","K"):1,("H","Y"):2,("I","L"):2,("I","M"):1,("I","V"):3,
       ("L","M"):2,("L","V"):1,("M","V"):1,("F","W"):1,("F","Y"):3,("S","T"):1,("W","Y"):2}


def best_alternative(one):
    alts = [(s, b if a == one else a) for (a, b), s in POS.items() if one in (a, b)]
    return max(alts, key=lambda x: (x[0], -AA.index(x[1])))[1] if alts else None


def residue_names(lines):
    names = {}
    for l in lines:
        if l.startswith(("ATOM", "HETATM")):
            names.setdefault(f"{l[21].strip()}{l[22:26].strip()}", l[17:20])
    return names


if __name__ == "__main__":
    spec, out_dir, out_spec = sys.argv[1:4]
    os.makedirs(out_dir, exist_ok=True)
    with open(out_spec, "w") as out:
        for line in open(spec):
            mid, pdb, query = line.rstrip("\n").split("\t")[:3]
            if not pdb.endswith(".pdb"):
                continue
            lines = open(pdb).readlines()
            names = residue_names(lines)
            for res in query.split(","):
                name = names.get(res)
                if name not in THREE:
                    continue
                alt = best_alternative(AA[THREE.index(name)])
                if alt:
                    break
            else:
                continue
            new = THREE[AA.index(alt)]
            chain, num = res[0], res[1:]
            path = os.path.join(out_dir, f"{mid}_{res}{alt}.pdb")
            with open(path, "w") as f:
                for l in lines:
                    if l.startswith(("ATOM", "HETATM")) and l[21] == chain and l[22:26].strip() == num:
                        l = l[:17] + new + l[20:]
                    f.write(l)
            out.write(f"{mid}\t{path}\t{query}\t{res}\t{name}>{new}\n")

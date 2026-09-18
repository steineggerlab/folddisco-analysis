"""Site queries: every residue with a CA within RADIUS of any annotated motif CA (first model).

Usage: make_site_queries.py <queryspec.tsv> <mutant_spec.tsv> <out.tsv> [radius=12]
out columns: id, pdb, mutant_pdb, site residues (-q), motif residues, n_site, mutated residue
"""
import math, sys


def ca_coords(path):
    out = {}
    with open(path) as f:
        for line in f:
            if line.startswith("ENDMDL"):
                break
            if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "CA" and line[26] == " ":
                key = f"{line[21].strip()}{line[22:26].strip()}"
                out.setdefault(key, (float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return out


if __name__ == "__main__":
    spec, mutant_spec, out_path = sys.argv[1:4]
    radius = float(sys.argv[4]) if len(sys.argv) > 4 else 12.0
    mutants = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(mutant_spec)}
    with open(out_path, "w") as out:
        for line in open(spec):
            mid, pdb, motif = line.rstrip("\n").split("\t")[:3]
            if mid not in mutants or not pdb.endswith(".pdb"):
                continue
            ca = ca_coords(pdb)
            centers = [ca[r] for r in motif.split(",") if r in ca]
            site = [r for r, xyz in ca.items() if len(r) > 1 and r[0].isalpha()
                    and any(math.dist(xyz, c) <= radius for c in centers)]
            out.write(f"{mid}\t{pdb}\t{mutants[mid][1]}\t{','.join(site)}\t{motif}\t{len(site)}\t{mutants[mid][3]}\n")

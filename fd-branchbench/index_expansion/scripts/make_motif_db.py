"""Write one motif-only PDB per (M-CSA entry, answer PDB): ATOM/HETATM records of the
annotated catalytic residues only (first model).

Usage: make_motif_db.py <queryspec.tsv> <out_dir> <out_manifest.tsv> [max_per_entry]
Row choice: the first homologue row for (entry, pdb) with no '-' whose residues all have a CA.
"""
import csv, os, sys
from collections import defaultdict
from multiprocessing import Pool

ROOT = os.environ.get("BENCH_DATA_ROOT", os.environ.get(
    "FOLDDISCO_BENCH_DATA",
    f"{os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))}/data"))
HOMOLOGUES = f"{ROOT}/zenodo/folddisco_data/v2/mcsa/data/metadata/catalytic_residues_homologues_parsed.tsv"
ANSWERS = f"{ROOT}/zenodo/folddisco_data/v2/mcsa/data/mcsa_answers"
STORE = f"{ROOT}/mcsa_rebuild/store"


def load_rows():
    rows = defaultdict(list)
    with open(HOMOLOGUES) as f:
        next(f)
        for mcsa, pdb, residues, _ in csv.reader(f, delimiter="\t"):
            rows[(int(mcsa), pdb.lower())].append(residues.split(","))
    return rows


ROWS = load_rows()


def residue_atoms(path):
    """(chain+resseq) -> [lines], first model, blank insertion code, no water."""
    atoms, ca = defaultdict(list), set()
    with open(path) as f:
        for line in f:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith(("ATOM", "HETATM")) or line[26] != " " or line[17:20] == "HOH":
                continue
            key = f"{line[21].strip()}{line[22:26].strip()}"
            atoms[key].append(line)
            if line[12:16].strip() == "CA":
                ca.add(key)
    return atoms, ca


def write_motif(task):
    mid, pdb, out_dir = task
    path = f"{STORE}/pdb{pdb}.pdb"
    if not os.path.exists(path):
        return None
    atoms, ca = residue_atoms(path)
    for row in ROWS.get((int(mid[1:]), pdb), []):
        if "-" in row:
            continue
        residues = list(dict.fromkeys(row))
        if len(residues) < 3 or any(r not in ca for r in residues):
            continue
        name = f"{mid}__pdb{pdb}.pdb"
        with open(os.path.join(out_dir, name), "w") as out:
            for r in residues:
                out.writelines(atoms[r])
            out.write("END\n")
        return (mid, pdb, name, ",".join(residues))
    return None


if __name__ == "__main__":
    spec, out_dir, manifest = sys.argv[1:4]
    cap = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    os.makedirs(out_dir, exist_ok=True)
    tasks = []
    for line in open(spec):
        mid = line.split("\t")[0]
        answers = [l.strip().removeprefix("pdb") for l in open(f"{ANSWERS}/{mid}_answer.tsv") if l.strip()]
        answers = sorted(dict.fromkeys(answers))
        if cap:
            answers = answers[:cap]
        tasks += [(mid, pdb, out_dir) for pdb in answers]
    with Pool(8) as pool, open(manifest, "w") as out:
        out.write("mcsa\tpdb\tfile\tresidues\n")
        n = 0
        for row in pool.imap(write_motif, tasks, chunksize=64):
            if row:
                out.write("\t".join(row) + "\n")
                n += 1
    print(f"tasks={len(tasks)} motifs={n}")

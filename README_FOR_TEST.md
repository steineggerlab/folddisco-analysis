# Folddisco test directory

## Command
### Querying
```bash
# default
./bin/folddisco query -i <INDEX> -p <QUERY_PDB> -q <QUERY_RESIDUES> -r -t <THREADS>
# -r flag is for residue matching & rmsd calculation
# -v flag is for verbose output (measures step-by-step time)

# Using the whole structure as a query
./bin/folddisco query -i <INDEX> -p <QUERY_PDB> -r -t <THREADS>

# Using query.txt file
./bin/folddisco query -i <INDEX> -q <QUERY_FILE> -r -t <THREADS>

# Using distance and angle threshold
./bin/folddisco query -i <INDEX> -p <QUERY_PDB> -q <QUERY_RESIDUES> -d <DISTANCE_THRESHOLD> -a <ANGLE_THRESHOLD> -r -t <THREADS>
```

```bash
# Example
# Zinc finger query to human proteome
./bin/folddisco query -i index/h_sapiens_folddisco -p query/1G2F.pdb -q F207,F212,F225,F229 -r -d 0.5 -a 5 -t 12
./bin/folddisco query -i index/h_sapiens_folddisco -q query/zinc_finger.txt -r -d 0.5 -a 5 -t 12

# Serine protease query to AFDB representatives
./bin/folddisco query -i index/afdb_rep_v4_folddisco -p query/4CHA.pdb -q B57,B102,C195 -r -t 12
./bin/folddisco query -i index/afdb_rep_v4_folddisco -q query/serine_protease.txt -r -t 12
```


## Index list
- `index/`
  - `afdb_rep_v4_folddisco`: AFDB representatives, 2.3M structures
  - `swissprot_folddisco`: SwissProt, 500K structures
  - `h_sapiens_folddisco`: Human proteome, 23K structures
  - `e_coli_folddisco`: E. coli proteome, 4K structures
  - `s_cerevisiae_folddisco`: Yeast proteome, 6K structures

## Query list
- `query/`
  - `1G2F.pdb`: Zinc finger protein
  - `4CHA.pdb`: Serine protease
  - `1LAP.pdb`: Aminopeptidase
  - `zinc_finger.txt`: 1G2F.pdb F207,F212,F225,F229
  - `serine_protease.txt`: 4CHA.pdb B57,B102,C195
  - `aminopeptidase.txt`: 1LAP.pdb 250,255,273,332,334
  - `knottin.txt`: 2N6N.pdb 3,10,15,16,21,23,28,30
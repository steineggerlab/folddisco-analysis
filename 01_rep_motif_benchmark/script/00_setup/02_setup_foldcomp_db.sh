# File: setup_foldcomp_db.sh
# Project: script
# Created: 2024-09-30 15:17:13
# Author: Hyunbin Kim (khb7840@gmail.com)
# Description:
#     This code is written as part of project "script".
# ---
# Last Modified: 2024-10-02 14:38:21
# Modified By: Hyunbin Kim (khb7840@gmail.com)
# ---
# Copyright © 2024 Hyunbin Kim, All rights reserved

# 1. Download foldcomp
# To download foldcomp database, use "download_foldcomp.py"
# cd foldcomp
# python ../download_foldcomp.py
# cd ..

# 2. Filter cifs in foldcomp
# Get the pdb list
grep ".pdb" foldcomp/h_sapiens.lookup | awk -F'\t' '{print $2}' > foldcomp/h_sapiens_pdb.list
grep ".pdb" foldcomp/e_coli.lookup | awk -F'\t' '{print $2}' > foldcomp/e_coli_pdb.list
grep ".pdb" foldcomp/s_cerevisiae.lookup | awk -F'\t' '{print $2}' > foldcomp/s_cerevisiae_pdb.list

# Run mmseqs2 createsubdb
./bin/mmseqs createsubdb foldcomp/h_sapiens_pdb.list foldcomp/h_sapiens foldcomp/h_sapiens_filtered --subdb-mode 1 --id-mode 1
./bin/mmseqs createsubdb foldcomp/e_coli_pdb.list foldcomp/e_coli foldcomp/e_coli_filtered --subdb-mode 1 --id-mode 1
./bin/mmseqs createsubdb foldcomp/s_cerevisiae_pdb.list foldcomp/s_cerevisiae foldcomp/s_cerevisiae_filtered --subdb-mode 1 --id-mode 1

# DONE
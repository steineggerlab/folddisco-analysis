#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: 01_run_foldmason_per_family.py
# Project: run
# Created: 2024-11-01 12:37:10
# Author: Hyunbin Kim (khb7840@gmail.com)
# Copyright © 2024 Hyunbin Kim, All rights reserved

import sys
import os

FOLDMASON_BIN = "~/Tools/foldmason/build/src/foldmason"
# FOLDMASON_BIN = "foldmason"
## WARNING: hard-coded paths. no ending '/' for the directories
input_dir = sys.argv[1]
pdb_dir = sys.argv[2] # /fast/hyunbin/motif/scop_benchmark/pdb/scope40-2.08/pdbstyle-2.08
output_dir = sys.argv[3]
temp_dir = sys.argv[4]
foldmason_command_file = sys.argv[5]
new_answer_dir = sys.argv[6]

# Get the list of domain files
domain_files = os.listdir(input_dir)
domain_files = [f for f in domain_files if f.endswith(".txt")]
print(f"Found {len(domain_files)} domain files.")
print(f"Writing foldmason commands to {foldmason_command_file}")


# Iterate over domain files
for domain_file in domain_files:
    print(f"Processing {domain_file}")
    command = f"{FOLDMASON_BIN} easy-msa"
    domain_id = domain_file.replace(".txt", "")
    pdb_ids = []
    with open(f"{input_dir}/{domain_file}", 'r') as f:
        for line in f:
            pdb_id = line.strip().split()[0]
            pdb_file = f"{pdb_dir}/{pdb_id[2:4]}/{pdb_id}.ent"
            # If pdb file contains multiple models, don't use it
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith("MODEL"):
                        print(f"Skipping {pdb_id} because it contains multiple models.")
                        break
            pdb_ids.append(pdb_id)

    # If len(pdb_ids) < 2, skip
    if len(pdb_ids) < 2:
        print(f"Skipping {domain_id} because the number of PDBs is less than 2.")
        continue
    else:
        with open(f"{new_answer_dir}/{domain_file}", 'w') as f:
            for pdb_id in pdb_ids:
                f.write(f"{pdb_id}\t{domain_id}\n")
        # Append pdb files to the command
        for pdb_id in pdb_ids:
            pdb_file = f"{pdb_dir}/{pdb_id[2:4]}/{pdb_id}.ent"
            command += f" {pdb_file}"
        command += f" {output_dir}/{domain_id} {temp_dir}"
        # Write the command to the foldmason command file
        with open(foldmason_command_file, 'a') as f:
            f.write(f"{command}\n")

# Execute the foldmason commands
print("Executing foldmason commands.")
os.system(f"bash {foldmason_command_file}")
print("Foldmason finished.")

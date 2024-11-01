#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: make_folddisco_query_from_foldmason_msa.py
# Description
#   This script gets foldmason MSA file and samples conserved columns to make FoldDisco queries.
#   Requires Biopython.
# Created: 2024-10-30 17:21:22
# Author: Hyunbin Kim (khb7840@gmail.com)
# Copyright © 2024 Hyunbin Kim, All rights reserved

# Usage: python make_folddisco_query_from_foldmason_msa.py foldmason_fasta_dir pdb_dir output_dir

from Bio import AlignIO
import numpy as np
import sys
import os

# QUERY_PDB_DIRECTORY="/fast/hyunbin/motif/scop_benchmark/pdb/pdbstyle-2.08/"
# WARNING: hard-coded paths. no ending '/' for the directories
foldmason_fasta_dir = sys.argv[1]
pdb_dir = sys.argv[2]
output_dir = sys.argv[3]

column_threshold = 0.8 
residue_threshold = 0.8
replicates = 3
query_output_prefix = "/fast/hyunbin/motif/scop_benchmark/temp/"
percentage_list = [0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
# print(f"/fast/hyunbin/motif/scop_benchmark/pdb/pdbstyle-2.08/gk/d2gkma_.ent\t{random_indices}\t/fast/hyunbin/motif/scop_benchmark/temp/d2gkma_{percentage}_{i}.tsv")

def choose_random_indices(indices, ratio, as_string=True):
    num_samples = int(len(indices) * ratio)
    if num_samples < 2:
        return None
    random_indices = np.random.choice(indices, num_samples, replace=False)
    random_indices = np.sort(random_indices)
    if as_string:
        random_indices = ",".join([str(i) for i in random_indices])
    return random_indices

def main():
    # Get the list of MSA files
    fasta_files = os.listdir(foldmason_fasta_dir)
    fasta_files = [f for f in fasta_files if f.endswith("_aa.fa")]

    for fasta_file in fasta_files:
        # Load the MSA file in A3M format
        msa = AlignIO.read(fasta_file, "fasta")
        # Get the domain ID
        base = os.path.basename(fasta_file)
        # domain_id is before the extension
        domain_id = os.path.splitext(base)[0]
        # Make a directory for the domain ID
        if not os.path.exists(f"{output_dir}/{domain_id}"):
            os.makedirs(f"{output_dir}/{domain_id}")

        # Set conservation threshold (e.g., 90% identity)
        num_seqs = len(msa)
        # Convert MSA to a matrix of characters, ignoring insertions (lowercase)
        alignment_matrix = []
        for record in msa:
            alignment_matrix.append([res for res in record.seq if res.isupper() or res == '-'])
        alignment_matrix = np.array(alignment_matrix)
        
        # Identify conserved columns
        conserved_columns = []
        for col_idx in range(alignment_matrix.shape[1]):
            column = alignment_matrix[:, col_idx]
            unique, counts = np.unique(column, return_counts=True)
            max_freq = max(counts) / num_seqs
            residue_freq = max(counts) / sum(counts)
            # If the most frequent residue meets the threshold, it's a conserved column
            if max_freq >= column_threshold and '-' not in unique and residue_freq >= residue_threshold:
                conserved_columns.append(col_idx)

        # For each sequence, translate the column indices to the corresponding positions in the sequence. Print the positions not the residues.
        for record in msa:
            positions = []
            residues = []
            start = 0
            for idx, res in enumerate(record.seq):
                # If dash, add to the start index
                if res == '-':
                    continue
                else:
                    # If the column index is in the conserved columns, add the position
                    if idx in conserved_columns:
                        positions.append(start)
                        residues.append(res)
                    start += 1
            # Sample from the conserved columns
            for percentage in percentage_list:
                if percentage == 1.0:
                    replicates_to_run = 1
                else:
                    replicates_to_run = replicates
                for i in range(replicates_to_run):
                    random_indices = choose_random_indices(positions, percentage)
                    if random_indices is None:
                        continue
                    pdb_id = record.id
                    pdb_file = f"{pdb_dir}/{pdb_id[2:4]}/{pdb_id}.ent"
                    print(f"{pdb_file}\t{random_indices}\t{output_dir}/{domain_id}/{record.id}_{percentage}_{i+1}.tsv")

if __name__ == "__main__":
    main()
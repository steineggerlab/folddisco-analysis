#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: 00_parse_scope40_fasta.py
# Project: run
# Created: 2024-11-01 13:27:45
# Author: Hyunbin Kim (khb7840@gmail.com)
# Copyright © 2024 Hyunbin Kim, All rights reserved

import sys
import os

scope_fasta = sys.argv[1]
output_prefix = sys.argv[2]

print(f"Preparing foldmason input from {scope_fasta}")

# Create a directory if output_prefix does not exist
output_dir = os.path.dirname(output_prefix)
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Just get the header. >id domain_id --> Return as tab-separated of id and domain_id
# domain_id_dict: key - domain_id, value - ids (list)
domain_id_dict = {}
with open(scope_fasta, 'r') as f:
    for line in f:
        if line.startswith('>'):
            entries = line[1:].strip().split()
            id = entries[0]
            domain_id = entries[1]
            if domain_id not in domain_id_dict:
                domain_id_dict[domain_id] = []
            domain_id_dict[domain_id].append(id)

# Write the domain_id_dict to each domain_id file in output_prefix
for domain_id, ids in domain_id_dict.items():
    with open(f"{output_prefix}/{domain_id}.txt", 'w') as f:
        # TSV: id, domain_id
        for id in ids:
            f.write(f"{id}\t{domain_id}\n")

print(f"Done. Output files are in {output_prefix}")
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: 03_run_benchmark.py
# Project: run
# Created: 2024-05-28 16:50:18
# Author: Hyunbin Kim (khb7840@gmail.com)
# Copyright © 2024 Hyunbin Kim, All rights reserved


# # Run folddisco query 

# FOLDDISCO_BIN=""
# DEFAULT_THREADS=128
# query_file=$1
# index_path=$2
# answer_dir=$3
# result_file=$4

# # Run folddisco query
# # ./folddisco query -q $query_file -i $index_path -d 0.5 -a 5 -t $DEFAULT_THREADS > 

# # Get all tsv files in the directory
# for file in temp/0.75/parsed/*.tsv
# do
#     # ./folddisco benchmark -r $file -a temptem -i scope_all >> temp/temp_result5.txt
#     ./folddisco benchmark -r $file -a temptem -i scope_all >> temp/temp_result3.txt
# done

import os
import sys
import multiprocessing as mp

# CONSTANTS
FOLDDISCO_BIN = ""
QUERY_THREADS = 128
BENCHMARK_THREADS = 16

# Get arguments
query_file = sys.argv[1]
index_path = sys.argv[2]
answer_dir = sys.argv[3]
result_dir = sys.argv[4]
benchmark_result_file = sys.argv[5]

# Run folddisco query
print("Running FoldDisco query")
os.system(f"{FOLDDISCO_BIN} query -q {query_file} -i {index_path} -d 0.5 -a 5 -t {QUERY_THREADS}")

# Get the domain list. Domain list can be obtained by getting the text files in the answer_dir
domain_list = os.listdir(answer_dir)
domain_list = [d for d in domain_list if d.endswith(".txt")]
domain_list = [d.split(".")[0] for d in domain_list]

# Run folddisco benchmark
for domain in domain_list:
    domain_answer = os.path.join(answer_dir, f"{domain}.txt")
    domain_result_dir = os.path.join(result_dir, domain)
    # List all tsv files in the directory
    tsv_files = os.listdir(domain_result_dir)
    for tsv_file in tsv_files:
        os.system(f"{FOLDDISCO_BIN} benchmark -r {tsv_file} -a {domain_answer} -i {index_path} >> {benchmark_result_file}")
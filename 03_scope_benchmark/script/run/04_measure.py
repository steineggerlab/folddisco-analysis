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
FOLDDISCO_BIN = "~/Projects/06_Motifsearch/motifsearch/target/release/folddisco"
THREADS = 64
# BENCHMARK_THREADS = 16

# Get arguments
query_file = sys.argv[1]
index_path = sys.argv[2]
answer_dir = sys.argv[3]
result_dir = sys.argv[4]
superfamily_dir = sys.argv[5]
benchmark_input_file = sys.argv[6]
benchmark_result_file = sys.argv[7]


# Run folddisco query
# print("Running FoldDisco query")
# os.system(f"{FOLDDISCO_BIN} query -q {query_file} -i {index_path} -d 0.5 -a 5 --serial -t {QUERY_THREADS}")

# Get the domain list. Domain list can be obtained by getting the text files in the answer_dir
family_list = os.listdir(answer_dir)
family_list = [d for d in family_list if d.endswith(".txt")]
family_list = [os.path.splitext(d)[0] for d in family_list]

benchmark_input = open(benchmark_input_file, 'w')

# Run folddisco benchmark
for family in family_list:
    family_answer = os.path.join(answer_dir, f"{family}.txt")
    family_result_dir = os.path.join(result_dir, family)
    family_split = family.split('.')
    superfamily = '.'.join(family_split[:3])
    superfamily_answer = os.path.join(superfamily_dir, f"{superfamily}.txt")
    # List all tsv files in the directory
    tsv_files = os.listdir(family_result_dir)
    for tsv_file in tsv_files:
        tsv_file = os.path.join(family_result_dir, tsv_file)
        print(f"{tsv_file}\t{family_answer}\t{superfamily_answer}", file=benchmark_input)
benchmark_input.close()

# Execute the benchmark
os.system(f"{FOLDDISCO_BIN} benchmark --input {benchmark_input_file} -i {index_path} -t {THREADS}> {benchmark_result_file}.tsv")
os.system(f"{FOLDDISCO_BIN} benchmark --input {benchmark_input_file} -i {index_path} -t {THREADS} --fp 5 > {benchmark_result_file}_fp5.tsv")
os.system(f"{FOLDDISCO_BIN} benchmark --input {benchmark_input_file} -i {index_path} -t {THREADS} --fp 1 > {benchmark_result_file}_fp1.tsv")
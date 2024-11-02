import os
import re
import sys

# 
file_to_parse = sys.argv[1]

header = "db\tresult\tanswer\tdb_length\tresult_length\tanswer_length\thashtype\tdistbin\tanglebin\tTP\tTN\tFP\tFN\tprecision\trecall\taccuracy\tf1\tratio\trepeat"

# Read the file content
with open(file_to_parse, "r") as file:
    with open(file_to_parse + ".parsed.tsv", "w") as output:
        output.write(header + "\n")
        for line in file:
            # Get the second column. 
            # /mnt/scratch/hyunbin/scop_result/c.62.1.0/d3zqka__0.8_3.tsv
            # Get 0.8 and 3
            entries = line.strip().split("\t")
            result = entries[1].split("/")[-1][:-4]
            result_split = result.split("_")
            repeat = result_split[-1]
            ratio = result_split[-2]
            print(f"{line.strip()}\t{ratio}\t{repeat}", file=output)

print("Done")
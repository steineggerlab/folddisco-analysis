import os
# 
directory = "temp/0.75/temp"
out_directory = "temp/0.75/parsed"
# Remove starting "./pdbstyle-2.08/**/" and ".ent" in the first column
# Save the first column at out_directory
for filename in os.listdir(directory):
    if filename.endswith(".tsv"):
        # get basename
        outfile = os.path.join(out_directory, filename)
        with open(os.path.join(directory, filename), "r") as file:
            for line in file:
                line = line.strip()
                line = line.split("\t")
                pdb_id = line[0].split("/")[-1].replace(".ent", "")
                score = line[1]
                total_matches = line[2]
                node_coverage = line[3]
                with open(outfile, "a") as out_file:
                    out_file.write(f"{pdb_id}\t{score}\t{total_matches}\t{node_coverage}\n")

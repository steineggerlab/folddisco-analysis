

data = "scope_benchmark/data/a1.1.1.temp.txt"

#prefix = "/fast/hyunbin/motif/scop_benchmark/pdb/pdbstyle-2.08/"
#id: d2gkma_
#group: gk
data_file = open(data, "r")
for line in data_file:
    line = line.strip()
    group = line[2:4]
    print(f"{group}/{line}.ent", end=" ")
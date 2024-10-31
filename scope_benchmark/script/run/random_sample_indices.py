
import numpy as np

temp_indices =[14, 16, 17, 18, 20, 21, 22, 26, 29, 32, 33, 36, 39, 40, 42, 44, 49, 50, 51, 52, 58, 61, 62, 64, 65, 67, 68, 69, 72, 74, 75, 76, 77, 78, 80, 81, 84, 86, 87, 92, 93, 95, 97, 98, 101, 107, 110, 112, 114, 118, 127]

percentage_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

replicates = 10

for percentage in percentage_list:
    num_samples = int(len(temp_indices) * percentage)
    for i in range(replicates):
        random_indices = np.random.choice(temp_indices, num_samples, replace=False)
        random_indices = np.sort(random_indices)
        random_indices = ",".join([str(i) for i in random_indices])
        print(f"/fast/hyunbin/motif/scop_benchmark/pdb/pdbstyle-2.08/gk/d2gkma_.ent\t{random_indices}\t/fast/hyunbin/motif/scop_benchmark/temp/d2gkma_{percentage}_{i}.tsv")
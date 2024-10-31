from Bio import AlignIO
import numpy as np

# Load the MSA file in A3M format
msa = AlignIO.read("scope40_a.1.1.1_result_aa.fa", "fasta")

# Set conservation threshold (e.g., 90% identity)
threshold = 0.5
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
    
    # If the most frequent residue meets the threshold, it's a conserved column
    if max_freq >= threshold and '-' not in unique:
        conserved_columns.append(col_idx)

# Output the conserved column indices
print("Conserved column indices:", conserved_columns)

# For each sequence, translate the column indices to the corresponding positions in the sequence. Print the positions not the residues.
for record in msa:
    positions = []
    residues = []
    start = 1
    for idx, res in enumerate(record.seq):
        # If dash, add to the start index
        if res == '-':
            continue
        else:
            # If the column index is in the conserved columns, add the position
            if idx in conserved_columns:
                positions.append(start + 1)
                residues.append(res)
            start += 1
    print(record.id, ":", positions)
    print(record.id, ":", residues)

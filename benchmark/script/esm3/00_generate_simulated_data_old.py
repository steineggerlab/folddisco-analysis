
# Commented out IPython magic to ensure Python compatibility.
# %set_env TOKENIZERS_PARALLELISM=false
# !pip install esm
import numpy as np
import torch
# !pip install py3Dmol
# import py3Dmol

from esm.utils.structure.protein_chain import ProteinChain
from esm.models.esm3 import ESM3
from esm.sdk import client
from esm.sdk.api import (
    ESMProtein,
    GenerationConfig,
)

from esm.utils.misc import huggingfacehub_login

huggingfacehub_login()  # will prompt you to get an API key and accept the ESM3 license.
model = ESM3.from_pretrained("esm3_sm_open_v1", device=torch.device("cuda"))

"""Alternatively, you could use the Forge API running the model remotely, and use the local `client` to call the API just like you're used to with the model running locally on your GPU:
"""

# from getpass import getpass
# token = getpass("Token from Forge console: ")
# model = client(
#     model="esm3-large-2024-03",
#     url="https://forge.evolutionaryscale.ai",
#     token=token,
# )

"""# Let's construct a prompt for ESM3, focusing on the task of scaffolding a motif from a natural protein
Zinc finger dataset generation
"""

pdb_id = "1G2F"  # PDB ID corresponding 
chain_id = "F"  # Chain ID corresponding
zinc_finger_chain = ProteinChain.from_rcsb(pdb_id, chain_id)
# Alternatively, we could have used ProteinChain.from_pdb() to load a protein structure from a local PDB file

motif_inds = np.arange(5, 28)
# `ProteinChain` objects can be indexed like numpy arrays to extract the sequence and atomic coordinates of a subset of residues
motif_sequence = zinc_finger_chain[motif_inds].sequence
motif_sequence = "C____C____________H___H"
motif_atom37_positions = zinc_finger_chain[motif_inds].atom37_positions

# TODO: partial match?
# cc_inds = np.arange(5, 11)
# cc_sequence = zinc_finger_chain[cc_inds].sequence
# num_partial_targets_per_length = 2

print("Motif sequence: ", motif_sequence)
print("Motif atom37_positions shape: ", motif_atom37_positions.shape)

length_list = [50, 100, 200, 400]
temperature_list = [0.2, 0.5, 0.8, 1.0]

setting_list = [(50, 0.2, 1), (50, 0.5, 1), (50, 0.8, 1), (50, 1.0, 1),
 (100, 0.2, 1), (100, 0.5, 1), (100, 0.8, 1), (100, 1.0, 1),
 (100, 0.2, 2), (100, 0.5, 2), (100, 0.8, 2), (100, 1.0, 2),
 (200, 0.2, 1), (200, 0.5, 1), (200, 0.8, 1), (200, 1.0, 1), 
 (200, 0.2, 2), (200, 0.5, 2), (200, 0.8, 2), (200, 1.0, 2), 
 (400, 0.2, 1), (400, 0.5, 1), (400, 0.8, 1), (400, 1.0, 1),
 (400, 0.2, 2), (400, 0.5, 2), (400, 0.8, 2), (400, 1.0, 2)]

false_repeat = 5

log_file = open("log.tsv", "w")

def generate_motif_prompt(prompt_length, temperature, motif_count):
    motif_starting_indices = []
    # First, we can construct a sequence prompt of all masks
    sequence_prompt = ["_"] * prompt_length
    # Then, randomly insert the motif sequence into the prompt 
    # Consider the motif count. Non-overlapping.
    if motif_count == 1:
        motif_start = np.random.randint(0, prompt_length - len(motif_sequence))
        sequence_prompt[motif_start : motif_start + len(motif_sequence)] = list(motif_sequence)
        motif_starting_indices.append(motif_start)
    elif motif_count == 2:
        # Get one at the first half, one at the second half
        motif_start = np.random.randint(0, prompt_length // 2 - len(motif_sequence))
        sequence_prompt[motif_start : motif_start + len(motif_sequence)] = list(motif_sequence)
        motif_starting_indices.append(motif_start)
        motif_start = np.random.randint(prompt_length // 2, prompt_length - len(motif_sequence))
        sequence_prompt[motif_start : motif_start + len(motif_sequence)] = list(motif_sequence)
        motif_starting_indices.append(motif_start)
    sequence_prompt = "".join(sequence_prompt)
    print("Sequence prompt: ", sequence_prompt)
    # Next, we can construct a structure prompt of all nan coordinates
    structure_prompt = torch.full((prompt_length, 37, 3), np.nan)
    # Then, we can insert the motif atomic coordinates into the prompt
    for motif_start in motif_starting_indices:
        structure_prompt[motif_start : motif_start + len(motif_atom37_positions)] = torch.tensor(
            motif_atom37_positions
        )
    protein_prompt = ESMProtein(sequence=sequence_prompt, coordinates=structure_prompt)
    
    sequence_generation_config = GenerationConfig(
        track="sequence",  # We want ESM3 to generate tokens for the sequence track
        num_steps=sequence_prompt.count("_") // 2,  # We'll use num(mask tokens) // 2 steps to decode the sequence
        temperature=temperature,  # We'll use a temperature of 0.5 to control the randomness of the decoding process
    )
    sequence_generation = model.generate(protein_prompt, sequence_generation_config)
    print("Sequence Prompt:\n\t", protein_prompt.sequence)
    print("Generated sequence:\n\t", sequence_generation.sequence)
    structure_prediction_config = GenerationConfig(
        track="structure",  # We want ESM3 to generate tokens for the structure track
        num_steps=len(sequence_generation) // 8,
        temperature=0.7,
    )
    structure_prediction_prompt = ESMProtein(sequence=sequence_generation.sequence)
    structure_prediction = model.generate(
        structure_prediction_prompt, structure_prediction_config
    )
    structure_prediction_chain = structure_prediction.to_protein_chain()
    # Print log
    pdb_name = f"zinc_finger_{prompt_length}_{temperature}_{motif_count}_{",".join(map(str, motif_starting_indices))}.pdb"
    log_file.write(f"{prompt_length}\t{temperature}\t{motif_count}\t{",".join(map(str, motif_starting_indices))}\t{sequence_generation.sequence}\t{pdb_name}\n")
    # Save the generated structure
    structure_prediction_chain.to_pdb(pdb_name)

def generate_false(prompt_length, temperature, motif_count):
    sequence_prompt = ["_"] * prompt_length
    sequence_prompt = "".join(sequence_prompt)
    structure_prompt = torch.full((prompt_length, 37, 3), np.nan)
    protein_prompt = ESMProtein(sequence=sequence_prompt, coordinates=structure_prompt)
    
    sequence_generation_config = GenerationConfig(
        track="sequence",  # We want ESM3 to generate tokens for the sequence track
        num_steps=sequence_prompt.count("_") // 2,  # We'll use num(mask tokens) // 2 steps to decode the sequence
        temperature=temperature,  # We'll use a temperature of 0.5 to control the randomness of the decoding process
    )
    sequence_generation = model.generate(protein_prompt, sequence_generation_config)
    print("Sequence Prompt:\n\t", protein_prompt.sequence)
    print("Generated sequence:\n\t", sequence_generation.sequence)
    structure_prediction_config = GenerationConfig(
        track="structure",  # We want ESM3 to generate tokens for the structure track
        num_steps=len(sequence_generation) // 8,
        temperature=0.7,
    )
    structure_prediction_prompt = ESMProtein(sequence=sequence_generation.sequence)
    structure_prediction = model.generate(
        structure_prediction_prompt, structure_prediction_config
    )
    structure_prediction_chain = structure_prediction.to_protein_chain()
    pdb_name = f"false_{prompt_length}_{temperature}_{motif_count}.pdb"
    log_file.write(f"{prompt_length}\t{temperature}\t{motif_count}\tNA\t{sequence_generation.sequence}\t{pdb_name}\n")
    # Save the generated structure
    structure_prediction_chain.to_pdb(pdb_name)
    
for setting in setting_list:
    for _ in range(false_repeat):
        generate_false(*setting)
    generate_motif_prompt(*setting)

# Read and Parse SCOP Hierarchy
# Generate Random Queries from Aligned SEQRES

import os
import sys
import random
import numpy as np
import pandas as pd

# Read SCOP Hierarchy
# FA-DOMID FA-PDBID FA-PDBREG FA-UNIID FA-UNIREG SF-DOMID SF-PDBID SF-PDBREG SF-UNIID SF-UNIREG SCOPCLA
def read_scop_hierarchy(scop_file):
    scop_hierarchy = pd.read_csv(
        scop_file, sep=' ', header=None, 
        names=[
            'FA-DOMID', 'FA-PDBID', 'FA-PDBREG', 'FA-UNIID', 'FA-UNIREG', 
            'SF-DOMID', 'SF-PDBID', 'SF-PDBREG', 'SF-UNIID', 'SF-UNIREG', 'SCOPCLA'
        ],
        comment='#'
    )
    return scop_hierarchy

def read_scope_hierarchy(scope_file):
    scope_hierarchy = pd.read_csv(
        scope_file, sep='\t', header=None, 
        names=[
            'id', 'pdb_id', 'range_str', 'classification', 'numeric_id', 'hieararchy'
        ],
        comment='#'
    )
    scope_hierarchy["range"] = scope_hierarchy["range_str"].apply(parse_res_index_range)
    return scope_hierarchy

def read_scope_description(scope_file):
    scope_description = pd.read_csv(
        scope_file, sep='\t', header=None, 
        names=[
            'numeric_id', 'type', 'classification', 'id', 'pdb_id', 'range_str'
        ],
        comment='#'
    )
    return scope_description

def retrieve_same_family(scope_description, fam):
    # Return list of ids that belong to the same family
    same_family = scope_description.loc[scope_description["classification"] == fam]
    same_family = same_family[same_family["type"] == "px"]
    return same_family["id"].tolist()

def retrieve_same_superfamily(scope_description, fam):
    # Get the superfamily id from the id. Remove the last digit after splitting with '.'
    fam = fam.split('.')
    superfamily = '.'.join(fam[:-1])
    # Return list of ids that belong to the same superfamily
    same_superfamily = scope_description.loc[scope_description["classification"].str.startswith(superfamily)]
    same_superfamily = same_superfamily[same_superfamily["type"] == "px"]
    return same_superfamily["id"].tolist()

def retrieve_same_fold(scope_description, fam):
    # Get the superfamily id from the id. Remove the last digit after splitting with '.'
    fam = fam.split('.')
    fold = '.'.join(fam[:-2])
    # Return list of ids that belong to the same superfamily
    same_fold = scope_description.loc[scope_description["classification"].str.startswith(fold)]
    same_fold = same_fold[scope_description["type"] == "px"]
    return same_fold["id"].tolist()

def retrieve_not_same_fold(scope_description, fam):
    # Get the superfamily id from the id. Remove the last digit after splitting with '.'
    fam = fam.split('.')
    fold = '.'.join(fam[:-2])
    # Return list of ids that belong to the same superfamily
    not_same_fold = scope_description.loc[~scope_description["classification"].str.startswith(fold)]
    not_same_fold = not_same_fold[scope_description["type"] == "px"]
    return not_same_fold["id"].tolist()

def parse_res_index_range(res_range):
    # Split with ',' first and then split with '-'
    res_range = res_range.split(',')
    res_index = []
    try:
        for r in res_range:
            chain = r.split(':')[0]
            r = r.split(':')[1]
            if '-' in r:
                r = r.split('-')
                for i in range(int(r[0]), int(r[1])+1):
                    res_index.append(chain + str(i))
            else:
                res_index.append(chain + r)
    except:
        print("Error in parsing res_index_range")
    return res_index

def pick_random_with_percent(seqres, percent):
    # Pick random residues with given percent
    res_index = []
    num_res = len(seqres)
    num_pick = int(num_res * percent)
    res_index = random.sample(seqres, num_pick)
    return res_index

def pick_random_with_num(seqres, num_pick):
    # Pick random residues with given number
    res_index = random.sample(seqres, num_pick)
    return res_index

def test():
    scope_description_file = "SCOPe/dir.des.scope.2.08-stable.txt"
    scope_classification_file = "SCOPe/dir.cla.scope.2.08-stable.txt"
    scope_description = read_scope_description(scope_description_file)
    print(f"Total {len(scope_description[scope_description['type'] == 'px'])} pdbs")
    test_fam = "a.1.1.1"
    same_family = retrieve_same_family(scope_description, test_fam)
    print(f"Same Family: Total {len(same_family)} pdbs")
    same_fold = retrieve_same_fold(scope_description, test_fam)
    print(f"Same Fold: Total {len(same_fold)} pdbs")
    not_same_fold = retrieve_not_same_fold(scope_description, test_fam)
    print(f"Not Same Fold: Total {len(not_same_fold)} pdbs")
    

def main():
    output = "scope_homeolike_query.txt"
    prefix = "/fast/hyunbin/motif/scop_benchmark/pdb/pdbstyle-2.08/"
    result_prefix = "result/scope_homeo/homeo_result_"
    intermediate = ""
    suffix = ".ent"
    output_file = open(output, 'w')
    repeat = 1
    # scope_file = "./SCOPe/homeo.subset"
    scope_file = "./SCOPe/homeolike.superfamily.txt"
    scope_hierarchy = read_scope_hierarchy(scope_file)
    # ratio_list = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
    ratio_list = [0.05]
    # count_list = [3, 4, 5, 7, 10]
    # Iterate over scope_hierarchy and make random queries
    for i in range(len(scope_hierarchy)):
        for j in range(repeat):
            for r in ratio_list:
            # for c in count_list:
                row = scope_hierarchy.iloc[i]
                intermediate = row["id"][2:4] + "/"
                res_index = pick_random_with_percent(row["range"], r)
                # res_index = pick_random_with_num(row["range"], c)
                output_file.write(f"{prefix}{intermediate}{row['id']}{suffix}\t{','.join(res_index)}\t{result_prefix}{row['id']}_{r}_{j}.tsv\n")
    output_file.close()

if __name__ == '__main__':
    test()
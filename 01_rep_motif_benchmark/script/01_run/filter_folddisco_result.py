"""
Q2L4Q9	69.4088	6	3	6	4	0	553	86.4346	A77,A128,A224:0.1457;A341,A382,A478:0.0892	A77,A128,A224:0.1457;A341,A382,A478:0.0892	B57,B102,C195	query/4CHA.pdb	index/h_sapiens_variation/h_sapiens_trrosetta
Q9BYE2	68.9071	6	3	6	4	0	586	73.6421	A366,A414,A511:0.0484	A366,A414,A511:0.0484	B57,B102,C195	query/4CHA.pdb	index/h_sapiens_variation/h_sapiens_trrosetta
Q8IU80	66.0943	6	3	6	4	0	811	85.2112	A617,A668,A762:0.0752	A617,A668,A762:0.0752	B57,B102,C195	query/4CHA.pdb	index/h_sapiens_variation/h_sapiens_trrosetta
P00749	50.2382	4	3	4	2	0	431	81.9607	A224,A275,A376:0.0511	A224,A275,A376:0.0511	B57,B102,C195	query/4CHA.pdb	index/h_sapiens_variation/h_sapiens_trrosetta
P37058	22.2790	2	2	2	2	0	310	93.8135	_,A106,A152:0.3277	_,A106,A152:0.3277	B57,B102,C195	query/4CHA.pdb	index/h_sapiens_variation/h_sapiens_trrosetta
Q9UI38	22.2095	2	2	2	2	0	385	74.3194	A153,A206,_:0.0437	A153,A206,_:0.0437	B57,B102,C195	query/4CHA.pdb	index/h_sapiens_variation/h_sapiens_trrosetta
Q8NGN1	22.1605	2	2	2	2	0	323	74.8819	_,A89,A11:0.4671	_,A89,A11:0.4671	B57,B102,C195	query/4CHA.pdb	index/h_sapiens_variation/h_sapiens_trrosetta
P34932	22.1430	2	2	2	0	0	840	85.9245	_,A32,A12:0.2386	_,A32,A12:0.2386	B57,B102,C195	query/4CHA.pdb	index/h_sapiens_variation/h_sapiens_trrosetta
Q68DD2	22.1122	2	2	2	0	0	849	83.6178	_,A151,A90:0.6314	_,A151,A90:0.6314	B57,B102,C195	query/4CHA.pdb	index/h_sapiens_variation/h_sapiens_trrosetta
"""

# Filter 9th column
# Filter condition
# 1. split by ";" first and then split by ":" -> get all matches and rmsd
# 2. calculate "_" in matches and get the maximum match count of all matches
# 3. check if the maximum match count is greater or equal to match_cutoff
# 4. get the rmsd of the maximum match count -> if rmsd is less than given cutoff, keep the line

# Usage: python filter_folddisco_result.py folddisco_result_file match_cutoff rmsd_cutoff

import sys

def filter_folddisco_result(folddisco_result_file, match_cutoff, rmsd_cutoff):
    with open(folddisco_result_file, "r") as f:
        for line in f:
            line = line.strip()
            if line == "":
                continue
            line_split = line.split("\t")
            
            # 9th column: matches and rmsd
            matches = line_split[10].split(";")
            max_match_count = 0
            max_rmsd = float('inf')

            # Loop over each match set
            for match in matches:
                match_split = match.split(":")
                match_residues = match_split[0].split(",")
                match_count = len([res for res in match_residues if res != "_"])  # Count valid residues
                
                rmsd = float(match_split[1]) if len(match_split) > 1 else 0  # Ensure RMSD is handled properly
                
                if match_count > max_match_count:
                    max_match_count = match_count
                    max_rmsd = rmsd
                if match_count == max_match_count and rmsd < max_rmsd:
                    max_rmsd = rmsd

            # Apply filters: max match count and RMSD cutoffs
            if max_match_count >= int(match_cutoff) and max_rmsd <= float(rmsd_cutoff):
                print(line)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python filter_folddisco_result.py folddisco_result_file match_cutoff rmsd_cutoff")
        sys.exit(1)

    folddisco_result_file = sys.argv[1]
    match_cutoff = sys.argv[2]
    rmsd_cutoff = sys.argv[3]
    filter_folddisco_result(folddisco_result_file, match_cutoff, rmsd_cutoff)

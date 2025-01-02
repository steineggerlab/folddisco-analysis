#!/bin/bash

# Output file for appending time and conditions
time_tsv_file="result/folddisco_timings.tsv"

# Write the header of the TSV file
echo -e "query_type\tquery_pdb\tquery_residue\tindex\tdist_bin\tangle_bin\tdist_threshold\tangle_threshold\tquery_num_threads\tmatching_residue_filter\tretrieve\thash_type\ttime\tnum_results\tresult_path" > $time_tsv_file

# Function to generate the index path based on conditions
generate_index_path() {
    local dist_bin=$1
    local angle_bin=$2
    local hash_type=$3
    local index_prefix=$4

    # Handle special case for non-h_sapiens indices
    if [[ "$index_prefix" != "h_sapiens" ]]; then
        echo "index/${index_prefix}_folddisco"
    # Determine the index path based on dist_bin, angle_bin, and hash_type for h_sapiens
    elif [[ "$dist_bin" == "16" && "$angle_bin" == "4" && "$hash_type" == "folddisco:pdbtr" ]]; then
        echo "index/h_sapiens_folddisco"
    elif [[ "$hash_type" == "folddisco:pdb" ]]; then
        echo "index/h_sapiens_variation/h_sapiens_pdb"
    elif [[ "$hash_type" == "folddisco:tr" ]]; then
        echo "index/h_sapiens_variation/h_sapiens_trrosetta"
    else
        echo "index/h_sapiens_variation/h_sapiens_d${dist_bin}a${angle_bin}"
    fi
}

# Function to run the command and append timing to TSV file
run_folddisco() {
    local query_pdb=$1
    local query_residue=$2
    local index_prefix=$3
    local dist_bin=$4
    local angle_bin=$5
    local dist_threshold=$6
    local angle_threshold=$7
    local query_num_threads=$8
    local matching_residue_filter=$9
    local retrieve=${10}
    local hash_type=${11}
    local query_type=${12}

    # Generate the correct index path based on dist_bin, angle_bin, and hash_type
    local index_path=$(generate_index_path "$dist_bin" "$angle_bin" "$hash_type" "$index_prefix")

    # Set the retrieve flag based on "yes" or "no"
    local retrieve_flag=""
    if [[ "$retrieve" == "yes" ]]; then
        retrieve_flag="--retrieve"
    fi

    # If index_prefix == pdb, use idtype pdb, else use idtype uniprot
    local idtype="uniprot"
    if [[ "$index_prefix" == "pdb" ]]; then
        idtype="pdb"
    fi

    # Create filenames for the output
    output_file="result/${query_type}_${index_prefix}_${hash_type}_${dist_bin}_${angle_bin}_${dist_threshold}_${angle_threshold}_${matching_residue_filter}_${retrieve}.tsv"

    # Measure time and capture output
    start_time=$(date +%s%N) # Record start time in nanoseconds
    /usr/bin/time -f "%e" -o temp_time ./bin/folddisco query \
        -p "query/$query_pdb" \
        -q "$query_residue" \
        -i "$index_path" \
        -d "$dist_threshold" \
        -a "$angle_threshold" \
        -t "$query_num_threads" \
        --node "$matching_residue_filter" $retrieve_flag \
        --id "$idtype" \
        > "$output_file"
    end_time=$(date +%s%N)   # Record end time in nanoseconds

    # Calculate time taken in seconds
    duration=$(echo "($end_time - $start_time) / 1000000000" | bc -l)
    # Format the time to 3 decimal places
    duration=$(printf "%.3f" $duration)

    # Count the number of results by counting the number of lines in the output file
    num_results=$(wc -l < "$output_file")

    # Append the input parameters and time to the TSV file
    echo -e "$query_type\t$query_pdb\t$query_residue\t$index_path\t$dist_bin\t$angle_bin\t$dist_threshold\t$angle_threshold\t$query_num_threads\t$matching_residue_filter\t$retrieve\t$hash_type\t$duration\t$num_results\t$output_file" >> $time_tsv_file
}

# List of all queries to run
run_folddisco "4CHA.pdb" "B57,B102,C195" "h_sapiens" 16 4 0.5 5 12 0 "yes" "folddisco:pdbtr" "serine_peptidase"
run_folddisco "4CHA.pdb" "B57,B102,C195" "h_sapiens" 16 4 0.5 5 12 3 "yes" "folddisco:pdbtr" "serine_peptidase"
run_folddisco "4CHA.pdb" "B57,B102,C195" "h_sapiens" 16 3 0.5 5 12 0 "yes" "folddisco:pdbtr" "serine_peptidase"
run_folddisco "4CHA.pdb" "B57,B102,C195" "h_sapiens" 16 3 0.5 5 12 3 "yes" "folddisco:pdbtr" "serine_peptidase"
run_folddisco "4CHA.pdb" "B57,B102,C195" "h_sapiens" 16 2 0.5 5 12 0 "yes" "folddisco:pdbtr" "serine_peptidase"
run_folddisco "4CHA.pdb" "B57,B102,C195" "h_sapiens" 16 2 0.5 5 12 3 "yes" "folddisco:pdbtr" "serine_peptidase"
run_folddisco "4CHA.pdb" "B57,B102,C195" "h_sapiens" 8 4 0.5 5 12 0 "yes" "folddisco:pdbtr" "serine_peptidase"
run_folddisco "4CHA.pdb" "B57,B102,C195" "h_sapiens" 8 4 0.5 5 12 3 "yes" "folddisco:pdbtr" "serine_peptidase"
run_folddisco "4CHA.pdb" "B57,B102,C195" "h_sapiens" 12 4 0.5 5 12 0 "yes" "folddisco:pdbtr" "serine_peptidase"
run_folddisco "4CHA.pdb" "B57,B102,C195" "h_sapiens" 12 4 0.5 5 12 3 "yes" "folddisco:pdbtr" "serine_peptidase"
run_folddisco "4CHA.pdb" "B57,B102,C195" "h_sapiens" 16 4 0.5 5 12 0 "yes" "folddisco:pdb" "serine_peptidase"
run_folddisco "4CHA.pdb" "B57,B102,C195" "h_sapiens" 16 4 0.5 5 12 3 "yes" "folddisco:pdb" "serine_peptidase"
run_folddisco "4CHA.pdb" "B57,B102,C195" "h_sapiens" 16 4 0.5 5 12 0 "yes" "folddisco:tr" "serine_peptidase"
run_folddisco "4CHA.pdb" "B57,B102,C195" "h_sapiens" 16 4 0.5 5 12 3 "yes" "folddisco:tr" "serine_peptidase"
run_folddisco "4CHA.pdb" "B57,B102,C195" "h_sapiens" 16 4 0.5 5 12 0 "no" "folddisco:pdbtr" "serine_peptidase"
run_folddisco "4CHA.pdb" "B57,B102,C195" "h_sapiens" 16 4 0.5 5 12 3 "no" "folddisco:pdbtr" "serine_peptidase"
run_folddisco "4CHA.pdb" "B57,B102,C195" "h_sapiens" 16 4 0 0 12 0 "yes" "folddisco:pdbtr" "serine_peptidase"
run_folddisco "4CHA.pdb" "B57,B102,C195" "h_sapiens" 16 4 0 0 12 0 "no" "folddisco:pdbtr" "serine_peptidase"
run_folddisco "4CHA.pdb" "B57,B102,C195" "e_coli" 16 4 0.5 5 12 3 "yes" "folddisco:pdbtr" "serine_peptidase"
run_folddisco "4CHA.pdb" "B57,B102,C195" "s_cerevisiae" 16 4 0.5 5 12 3 "yes" "folddisco:pdbtr" "serine_peptidase"
run_folddisco "4CHA.pdb" "B57,B102,C195" "e_coli" 16 4 0.5 5 12 0 "yes" "folddisco:pdbtr" "serine_peptidase"
run_folddisco "4CHA.pdb" "B57,B102,C195" "s_cerevisiae" 16 4 0.5 5 12 0 "yes" "folddisco:pdbtr" "serine_peptidase"
run_folddisco "4CHA.pdb" "B57,B102,C195" "pdb" 16 4 0.5 5 12 0 "yes" "folddisco:pdbtr" "serine_peptidase"
run_folddisco "4CHA.pdb" "B57,B102,C195" "swissprot" 16 4 0.5 5 12 0 "yes" "folddisco:pdbtr" "serine_peptidase"
run_folddisco "4CHA.pdb" "B57,B102,C195" "afdb_rep_v4" 16 4 0.5 5 12 0 "yes" "folddisco:pdbtr" "serine_peptidase"
run_folddisco "4CHA.pdb" "B57,B102,C195" "pdb" 16 4 0.5 5 12 3 "yes" "folddisco:pdbtr" "serine_peptidase"
run_folddisco "4CHA.pdb" "B57,B102,C195" "swissprot" 16 4 0.5 5 12 3 "yes" "folddisco:pdbtr" "serine_peptidase"
run_folddisco "4CHA.pdb" "B57,B102,C195" "afdb_rep_v4" 16 4 0.5 5 12 3 "yes" "folddisco:pdbtr" "serine_peptidase"
run_folddisco "1G2F.pdb" "F207,F212,F225,F229" "h_sapiens" 16 4 0.5 5 12 0 "yes" "folddisco:pdbtr" "zinc_finger"
run_folddisco "1G2F.pdb" "F207,F212,F225,F229" "h_sapiens" 16 4 0.5 5 12 3 "yes" "folddisco:pdbtr" "zinc_finger"
run_folddisco "1G2F.pdb" "F207,F212,F225,F229" "h_sapiens" 16 3 0.5 5 12 0 "yes" "folddisco:pdbtr" "zinc_finger"
run_folddisco "1G2F.pdb" "F207,F212,F225,F229" "h_sapiens" 16 3 0.5 5 12 3 "yes" "folddisco:pdbtr" "zinc_finger"
run_folddisco "1G2F.pdb" "F207,F212,F225,F229" "h_sapiens" 16 2 0.5 5 12 0 "yes" "folddisco:pdbtr" "zinc_finger"
run_folddisco "1G2F.pdb" "F207,F212,F225,F229" "h_sapiens" 16 2 0.5 5 12 3 "yes" "folddisco:pdbtr" "zinc_finger"
run_folddisco "1G2F.pdb" "F207,F212,F225,F229" "h_sapiens" 8 4 0.5 5 12 0 "yes" "folddisco:pdbtr" "zinc_finger"
run_folddisco "1G2F.pdb" "F207,F212,F225,F229" "h_sapiens" 8 4 0.5 5 12 3 "yes" "folddisco:pdbtr" "zinc_finger"
run_folddisco "1G2F.pdb" "F207,F212,F225,F229" "h_sapiens" 12 4 0.5 5 12 0 "yes" "folddisco:pdbtr" "zinc_finger"
run_folddisco "1G2F.pdb" "F207,F212,F225,F229" "h_sapiens" 12 4 0.5 5 12 3 "yes" "folddisco:pdbtr" "zinc_finger"
run_folddisco "1G2F.pdb" "F207,F212,F225,F229" "h_sapiens" 16 4 0.5 5 12 0 "yes" "folddisco:pdb" "zinc_finger"
run_folddisco "1G2F.pdb" "F207,F212,F225,F229" "h_sapiens" 16 4 0.5 5 12 3 "yes" "folddisco:pdb" "zinc_finger"
run_folddisco "1G2F.pdb" "F207,F212,F225,F229" "h_sapiens" 16 4 0.5 5 12 0 "yes" "folddisco:tr" "zinc_finger"
run_folddisco "1G2F.pdb" "F207,F212,F225,F229" "h_sapiens" 16 4 0.5 5 12 3 "yes" "folddisco:tr" "zinc_finger"
run_folddisco "1G2F.pdb" "F207,F212,F225,F229" "h_sapiens" 16 4 0.5 5 12 0 "no" "folddisco:pdbtr" "zinc_finger"
run_folddisco "1G2F.pdb" "F207,F212,F225,F229" "h_sapiens" 16 4 0.5 5 12 3 "no" "folddisco:pdbtr" "zinc_finger"
run_folddisco "1G2F.pdb" "F207,F212,F225,F229" "h_sapiens" 16 4 0 0 12 0 "yes" "folddisco:pdbtr" "zinc_finger"
run_folddisco "1G2F.pdb" "F207,F212,F225,F229" "h_sapiens" 16 4 0 0 12 0 "no" "folddisco:pdbtr" "zinc_finger"
run_folddisco "1G2F.pdb" "F207,F212,F225,F229" "e_coli" 16 4 0.5 5 12 3 "yes" "folddisco:pdbtr" "zinc_finger"
run_folddisco "1G2F.pdb" "F207,F212,F225,F229" "s_cerevisiae" 16 4 0.5 5 12 3 "yes" "folddisco:pdbtr" "zinc_finger"
run_folddisco "1G2F.pdb" "F207,F212,F225,F229" "pdb" 16 4 0.5 5 12 3 "yes" "folddisco:pdbtr" "zinc_finger"
run_folddisco "1G2F.pdb" "F207,F212,F225,F229" "swissprot" 16 4 0.5 5 12 3 "yes" "folddisco:pdbtr" "zinc_finger"
run_folddisco "1G2F.pdb" "F207,F212,F225,F229" "afdb_rep_v4" 16 4 0.5 5 12 3 "yes" "folddisco:pdbtr" "zinc_finger"
run_folddisco "1G2F.pdb" "F207,F212,F225,F229" "e_coli" 16 4 0.5 5 12 0 "yes" "folddisco:pdbtr" "zinc_finger"
run_folddisco "1G2F.pdb" "F207,F212,F225,F229" "s_cerevisiae" 16 4 0.5 5 12 0 "yes" "folddisco:pdbtr" "zinc_finger"
run_folddisco "1G2F.pdb" "F207,F212,F225,F229" "pdb" 16 4 0.5 5 12 0 "yes" "folddisco:pdbtr" "zinc_finger"
run_folddisco "1G2F.pdb" "F207,F212,F225,F229" "swissprot" 16 4 0.5 5 12 0 "yes" "folddisco:pdbtr" "zinc_finger"
run_folddisco "1G2F.pdb" "F207,F212,F225,F229" "afdb_rep_v4" 16 4 0.5 5 12 0 "yes" "folddisco:pdbtr" "zinc_finger"
run_folddisco "1LAP.pdb" "250,255,273,332,334" "h_sapiens" 16 4 0.5 5 12 0 "yes" "folddisco:pdbtr" "aminopeptidase"
run_folddisco "1LAP.pdb" "250,255,273,332,334" "h_sapiens" 16 4 0.5 5 12 5 "yes" "folddisco:pdbtr" "aminopeptidase"
run_folddisco "1LAP.pdb" "250,255,273,332,334" "e_coli" 16 4 0.5 5 12 5 "yes" "folddisco:pdbtr" "aminopeptidase"
run_folddisco "1LAP.pdb" "250,255,273,332,334" "s_cerevisiae" 16 4 0.5 5 12 5 "yes" "folddisco:pdbtr" "aminopeptidase"
run_folddisco "1LAP.pdb" "250,255,273,332,334" "pdb" 16 4 0.5 5 12 5 "yes" "folddisco:pdbtr" "aminopeptidase"
run_folddisco "1LAP.pdb" "250,255,273,332,334" "swissprot" 16 4 0.5 5 12 5 "yes" "folddisco:pdbtr" "aminopeptidase"
run_folddisco "1LAP.pdb" "250,255,273,332,334" "afdb_rep_v4" 16 4 0.5 5 12 5 "yes" "folddisco:pdbtr" "aminopeptidase"
run_folddisco "2N6N.pdb" "3,10,15,16,21,23,28,30" "h_sapiens" 16 4 0.5 5 12 0 "yes" "folddisco:pdbtr" "knottin"
run_folddisco "2N6N.pdb" "3,10,15,16,21,23,28,30" "h_sapiens" 16 4 0.5 5 12 6 "yes" "folddisco:pdbtr" "knottin"
run_folddisco "2N6N.pdb" "3,10,15,16,21,23,28,30" "e_coli" 16 4 0.5 5 12 6 "yes" "folddisco:pdbtr" "knottin"
run_folddisco "2N6N.pdb" "3,10,15,16,21,23,28,30" "s_cerevisiae" 16 4 0.5 5 12 6 "yes" "folddisco:pdbtr" "knottin"
run_folddisco "2N6N.pdb" "3,10,15,16,21,23,28,30" "pdb" 16 4 0.5 5 12 6 "yes" "folddisco:pdbtr" "knottin"
run_folddisco "2N6N.pdb" "3,10,15,16,21,23,28,30" "swissprot" 16 4 0.5 5 12 6 "yes" "folddisco:pdbtr" "knottin"
run_folddisco "2N6N.pdb" "3,10,15,16,21,23,28,30" "afdb_rep_v4" 16 4 0.5 5 12 6 "yes" "folddisco:pdbtr" "knottin"


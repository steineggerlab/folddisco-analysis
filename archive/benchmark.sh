

FOLDDISCO_PATH=~/Projects/06_Motifsearch/motifsearch/target/release/motifsearch
RESULTS_DIR=result/scope_homeo
ANSWER=result/homeo_answer.tsv
INDEX=scop_benchmark/ver2_08/pdbtr/merged

# ~/Projects/06_Motifsearch/motifsearch/target/release/motifsearch benchmark -r result/scope_homeo/homeo_result_*.tsv -a result/homeo_answer.tsv -i scop_benchmark/ver2_08/pdbtr/merged  --format tsv > ~/homeo.tsv

# Iterate over all tsv files in the result directory
# for file in result/scope_homeo/homeo_result_*.tsv
# do
#     # Run the benchmark
#     $FOLDDISCO_PATH benchmark -r $file -a $ANSWER -i $INDEX --format tsv >> ~/
# done

# Use GNU parallel to run the benchmark in parallel
ls result/scope_homeo/homeo_result_*.tsv | parallel "$FOLDDISCO_PATH benchmark -r {} -a $ANSWER -i $INDEX --format tsv >> ~/homeo.tsv"
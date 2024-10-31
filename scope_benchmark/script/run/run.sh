#$ Temp

# ratio_list=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)
# rep_list=(0 1 2 3 4 5 6 7 8 9)

# for ratio in "${ratio_list[@]}"
# do
#     for rep in "${rep_list[@]}"
#     do
#         ./folddisco_fp1 benchmark -r temp/0.5/parsed/d2gkma_${ratio}_${rep}.tsv -a a.1.1.1.answer.txt -i scope_all >> temp/temp_result.txt
#     done
# done


# Get all tsv files in the directory
for file in temp/0.75/parsed/*.tsv
do
    # ./folddisco benchmark -r $file -a temptem -i scope_all >> temp/temp_result5.txt
    ./folddisco benchmark -r $file -a temptem -i scope_all >> temp/temp_result3.txt
done
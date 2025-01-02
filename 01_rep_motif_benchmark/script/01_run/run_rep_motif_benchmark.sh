#!
# Script to run benchmark 

#\time ./target/release/folddisco query -p analysis/query/betaturn.pdb -i analysis/h_sapiens/d16a4/index -t 6 -d 0.5 -a 5.0 -v --header
# NUM_THREADS=32
NUM_THREADS=8
DEFAULT_FLAGS="-v"
#FOLDDISCO_INDEX=$1
FOLDDISCO_INDEX=index/pdbtr_d16a4/index
# Run Folddisco
./folddisco query -q query/folddisco_query.tsv -i $FOLDDISCO_INDEX -t $NUM_THREADS $DEFAULT_FLAGS

# # Retrieve
\time ./folddisco query -p query/betaturn.pdb -r -i $FOLDDISCO_INDEX -t $NUM_THREADS $DEFAULT_FLAGS > result/folddisco/res_mapped/betaturn.tsv
\time ./folddisco query -p query/4CHA.pdb -q B57,B102,C195 -r -i $FOLDDISCO_INDEX -t $NUM_THREADS $DEFAULT_FLAGS --node 3 > result/folddisco/res_mapped/serine_protease.tsv
\time ./folddisco query -p query/1G2F.pdb -q F207,F212,F225,F229 -r -i $FOLDDISCO_INDEX -t $NUM_THREADS $DEFAULT_FLAGS --node 3 > result/folddisco/res_mapped/zinc_finger.tsv
\time ./folddisco query -p query/1LAP.pdb -q 250,255,273,332,334 -r -i $FOLDDISCO_INDEX -t $NUM_THREADS $DEFAULT_FLAGS --node 4 > result/folddisco/res_mapped/aminopeptidase.tsv
\time ./folddisco query -p query/2N6N.pdb -q 3,10,15,16,21,23,28,30 -r -i $FOLDDISCO_INDEX -t $NUM_THREADS $DEFAULT_FLAGS --node 6 > result/folddisco/res_mapped/knottin.tsv
\time ./folddisco query -p query/homeobox.pdb -r -i $FOLDDISCO_INDEX -t $NUM_THREADS $DEFAULT_FLAGS --node 15 > result/folddisco/res_mapped/homeobox.tsv

# Run Pyscomotif
# TODO:


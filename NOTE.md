# TODOs
- [ ] Define parameters to benchmark
- [ ] Reorganize code & files from archive

Figure 1. 
Comparison with other tools


Figure 2.
Comparison of features within folddisco


---

Index: analysis/h_sapiens/pdb_d16a4.lookup
Result: temp_zinc_pdb
Answer: data/zinc_answer.tsv
Total ids: 23391
Result length: 4983
Answer length: 1817
Hash type: PDBMotifSinCos
Number of distance bins: 16
Number of angle bins: 4
TP: 1412
TN: 23391
FP: 3571
FN: 405
Precision: 0.2834
Recall: 0.7771
Accuracy: 0.8618
F1 score: 0.4153

Index: analysis/h_sapiens/pdb_d16a4.lookup
Result: temp_zinc_pdb_node3
Answer: data/zinc_answer.tsv
Total ids: 23391
Result length: 1029
Answer length: 1817
Hash type: PDBMotifSinCos
Number of distance bins: 16
Number of angle bins: 4
TP: 764
TN: 23391
FP: 265
FN: 1053
Precision: 0.7425
Recall: 0.4205
Accuracy: 0.9483
F1 score: 0.5369

Index: index/pdbtr_d16a4/id.lookup
Result: result/folddisco/inverted_index/zinc_finger_id.tsv
Answer: result/zinc_answer.tsv
Total ids: 20504
Result length: 1047
Answer length: 1817
Hash type: PDBTrRosetta
Number of distance bins: 16
Number of angle bins: 4
TP: 822
TN: 18471
FP: 225
FN: 995
Precision: 0.7851
Recall: 0.4524
Accuracy: 0.9405
F1 score: 0.5740
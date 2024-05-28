# Use seaborn or matplotlib to plot the data

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Data looks like this
# index	query	fraction	rep	answer	total	result_len	answer_len	type	nbin_dist	nbin_angle	tp	tn	fp	fn	precision	recall	accuracy	f1_score
# scop_benchmark/ver2_08/pdbtr/merged	result/scope_homeo/homeo_result_d1ui5a1	0.05	0	result/homeo_answer	344850	0	269	PDBTrRosetta	16	4	0	344581	0	269	NaN	0	0.9992	NaN
# scop_benchmark/ver2_08/pdbtr/merged	result/scope_homeo/homeo_result_d1ui5a1	0.05	1	result/homeo_answer	344850	0	269	PDBTrRosetta	16	4	0	344581	0	269	NaN	0	0.9992	NaN
# Making boxplots for each fraction and showing the distribution of precision, recall, accuracy, f1_score
# data = pd.read_csv("homeo.9.tsv", sep='\t', header=0)
data = pd.read_csv("homeolike.fp5.tsv", sep='\t', header=0)
# Handle NaN --> change to 0
data = data.fillna(0)

# Use font "Inter"
plt.rcParams['font.family'] = 'Inter'
flierprops = dict(marker='o', markersize=1.5, linestyle='none')

# Make boxplots of precision, recall, accuracy, f1_score
sns.set(style="white")
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
sns.boxplot(y="precision", x="fraction", data=data, ax=axes[0], orient='v', hue='fraction', palette="vlag", width=0.6, dodge=False, flierprops=flierprops)
sns.boxplot(y="recall", x="fraction", data=data, ax=axes[1], orient='v', hue='fraction', palette="vlag", width=0.6, dodge=False, flierprops=flierprops)
sns.boxplot(y="f1_score", x="fraction", data=data, ax=axes[2], orient='v', hue='fraction', palette="vlag", width=0.6, dodge=False, flierprops=flierprops)
# Remove the legend
axes[0].get_legend().remove()
axes[1].get_legend().remove()
axes[2].get_legend().remove()
plt.tight_layout()
# plt.show()
# Save the plot
plt.savefig("homeolike.fp5.png", dpi=300)

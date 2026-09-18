"""Sens@1FP and runtime per config for the motif-index test (theme_pubr look)."""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(D, "..", "scripts"))
import importlib.util
spec = importlib.util.spec_from_file_location("p3", os.path.join(D, "..", "scripts", "03_plot.py"))
p3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(p3)

s = pd.read_csv(f"{D}/result/summary.tsv", sep="\t")
ORDER = [("I0", "plain index"), ("I0_sens", "plain, --sensitive"), ("I0_aa", "plain, --aa-subst"),
         ("I0_star", "plain, :* on mutated"), ("I1", "index r1"), ("I2", "index r2"),
         ("I3", "index blosum62"), ("I4", "index r1+blosum62")]
RANK = {"idf": ("#E94B8B", "rank by IDF (default)"), "coverage": ("#655AD0", "rank by coverage, RMSD")}

p3.pubr()
fig, axes = plt.subplots(1, 3, figsize=(15, 5.2), gridspec_kw={"width_ratios": [1, 1.15, 1]})
for ax, qset, title in [(axes[0], "orig", "Published query sites"), (axes[1], "mut", "Mutant query sites")]:
    keys = [(k, lab) for k, lab in ORDER if f"{qset}_{k}" in set(s.config)]
    for j, (rank, (color, _)) in enumerate(RANK.items()):
        v = [s[(s.config == f"{qset}_{k}") & (s.ranking == rank)].sens1fp.iloc[0] for k, _ in keys]
        x = np.arange(len(keys)) + (j - 0.5) * 0.4
        ax.bar(x, v, width=0.38, color=color, edgecolor="white", linewidth=0.8, zorder=3)
        for xi, vi in zip(x, v):
            ax.annotate(f"{vi:.3f}", (xi, vi), textcoords="offset points", xytext=(0, 2), ha="center",
                        va="bottom", fontsize=6.5, color=p3.MUTED, rotation=90)
    ax.set_xticks(range(len(keys))); ax.set_xticklabels([lab for _, lab in keys], rotation=35, ha="right", fontsize=8.5)
    ax.set_ylim(0, 0.62); ax.set_ylabel("Mean Sens@1FP"); ax.set_title(title)

ax = axes[2]
r = s[s.ranking == "idf"].set_index("config").runtime
keys = [(k, lab) for k, lab in ORDER if f"orig_{k}" in r.index]
v = [r[f"orig_{k}"] for k, _ in keys]
ax.bar(range(len(keys)), v, width=0.6, color="#4FBDF2", edgecolor="white", linewidth=0.8, zorder=3)
for i, vi in enumerate(v):
    ax.annotate(f"{vi:.1f}s", (i, vi), textcoords="offset points", xytext=(0, 2), ha="center", fontsize=7.5, color=p3.MUTED)
ax.set_xticks(range(len(keys))); ax.set_xticklabels([lab for _, lab in keys], rotation=35, ha="right", fontsize=8.5)
ax.set_ylabel("Mean wall time per query (s, 2 threads, 4 parallel)"); ax.set_title("Runtime, published sites")
ax.set_ylim(0, max(v) * 1.2)

fig.suptitle("M-CSA motif-only index (24,762 motifs): query-side vs index-side expansion (n=246)",
             fontsize=12.5, fontweight="bold", y=0.995)
fig.tight_layout(rect=[0, 0.06, 1, 0.95])
h = [plt.Rectangle((0, 0), 1, 1, fc=c, ec="none") for c, _ in RANK.values()]
fig.legend(h, [lab for _, lab in RANK.values()], loc="lower center", ncol=2, bbox_to_anchor=(0.36, 0.0), fontsize=9.5)
os.makedirs(f"{D}/figure", exist_ok=True)
for e in ("png", "pdf"):
    fig.savefig(f"{D}/figure/index_expansion.{e}", dpi=180)
print("ok")

#!/usr/bin/env python
"""Figures for the default sort order and the --confident filter (see analyze.py)."""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, f"{os.path.dirname(HERE)}/scripts")
from importlib import import_module
pubr = import_module("03_plot").pubr  # same ggpubr-style rcParams as the main figures

# Paper palette: gray for the 2.x baseline, folddisco pink for what 3.0 ships.
C = ["#8E8E93", "#655AD0", "#4FBDF2", "#B04BB6", "#E94B8B"]
INK, MUTED = "#0b0b0b", "#52514e"
RES = f"{HERE}/result"
FIG = f"{HERE}/figure"


def bars(ax, labels, values, colors, fmt="{:.4f}", ylabel="", title="", ymax=None):
    for i, (v, c) in enumerate(zip(values, colors)):
        ax.bar(i, v, width=0.7, color=c, edgecolor="white", linewidth=0.8, zorder=3)
        ax.annotate(fmt.format(v), (i, v), textcoords="offset points", xytext=(0, 2),
                    ha="center", va="bottom", fontsize=8, color=MUTED)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=8, rotation=30, ha="right")
    ax.set_ylabel(ylabel); ax.set_title(title, fontsize=10)
    ax.set_ylim(0, (ymax or max(values)) * 1.18)


def dots(ax, labels, values, colors, fmt="{:.4f}", ylabel="", title="", pad=0.012):
    """Candidate orderings differ by ~0.02-0.08, too little to read as bars from zero."""
    x = np.arange(len(labels))
    ax.plot(x, values, color="#c9c7c4", lw=1.4, zorder=2)
    for i, (v, c) in enumerate(zip(values, colors)):
        ax.plot(i, v, marker="o", ms=9 if c == C[-1] else 7.5, color=c, mec="white", mew=1.2,
                zorder=3)
        ax.annotate(fmt.format(v), (i, v), textcoords="offset points", xytext=(0, 9),
                    ha="center", fontsize=8, color=c, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8, rotation=30, ha="right")
    ax.set_ylabel(ylabel); ax.set_title(title, fontsize=10)
    ax.set_xlim(-0.4, len(labels) - 0.6)
    ax.set_ylim(min(values) - pad, max(values) + pad * 2.2)
    ax.grid(axis="y", color="#e6e4e1", lw=0.7, zorder=0)


def fig_sort():
    s = pd.read_csv(f"{RES}/sort_per_query.tsv.gz", sep="\t")
    comp = pd.read_csv(f"{RES}/composite_per_query.tsv.gz", sep="\t")
    d = pd.concat([s, comp[~comp["sort"].isin(s["sort"].unique())]])
    pm = d[d.run.str.startswith("mc_pm")].groupby("sort")[["sens1", "ap"]].mean()
    ps = d[d.run.str.startswith("mc_ps")].groupby("sort")[["sens1", "ap"]].mean()
    c2 = pd.read_csv(f"{RES}/composite2_raw.tsv.gz", sep="\t")
    d = pd.concat([d, c2[~c2["sort"].isin(d["sort"].unique())]])
    pm = d[d.run.str.startswith("mc_pm")].groupby("sort")[["sens1", "ap"]].mean()
    ps = d[d.run.str.startswith("mc_ps")].groupby("sort")[["sens1", "ap"]].mean()
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 5.0))
    order = ["idf,rmsd", "node_count,idf,rmsd", "idf*cov^1.0,rmsd", "idf*cov^2.0*rmsd/0.5,rmsd", "idf*cov^2.0*tm_score,rmsd"]
    lab = ["idf, rmsd\n(2.x default)", "node_count,\nidf, rmsd", "idf x cov, rmsd", "idf x cov²\nx rmsd term", "match_score\n(3.0)"]
    dots(axes[0], lab, [pm.loc[o, "sens1"] for o in order], C, ylabel="Mean Sens@1FP",
         title="Per match, M-CSA q250 (4 configs)")
    dots(axes[1], lab, [pm.loc[o, "ap"] for o in order], C, ylabel="Mean average precision",
         title="Per match, average precision")
    order2 = ["idf,min_rmsd", "max_node_count,min_rmsd", "cov^2.0*idf^0.0*rmsd/1,min_rmsd", "cov^2.0*idf^0.5*rmsd/1,min_rmsd"]
    lab2 = ["idf, min_rmsd\n(2.x default)", "max_node_count,\nmin_rmsd", "matched²\n/(1+rmsd)", "structure_score\n(3.0)"]
    dots(axes[2], lab2, [ps.loc[o, "sens1"] for o in order2], [C[0], C[1], C[2], C[4]],
         ylabel="Mean Sens@1FP", title="Per structure, M-CSA q250")
    fig.suptitle("Default ranking: 281 lexicographic orders and 69 evidence x quality scores",
                 fontsize=12.5, fontweight="bold", y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    save(fig, "fig_sort")


def fig_confident():
    f = pd.read_csv(f"{RES}/filter_per_query.tsv.gz", sep="\t")
    f["prec"] = np.where(f.hits > 0, f.tp / f.hits, np.nan)
    f["rec"] = f.tp / f.answers
    sel = [(0.0, np.inf, "no filter"), (0.8, np.inf, "cov ≥ 0.8"), (1.0, 1.0, "cov = 1.0\nRMSD ≤ 1.0"),
           (0.8, 0.75, "cov ≥ 0.8\nRMSD ≤ 0.75"), (0.8, 1.0, "cov ≥ 0.8\nRMSD ≤ 1.0 (3.0)")]
    mc = f[f.run == "mc_pm_def"]
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.6))
    rows = [mc[(mc["cov"] == c) & (mc.rmsd == r) & (mc.drmsd == np.inf)] for c, r, _ in sel]
    labels = [s[2] for s in sel]
    bars(axes[0], labels, [g.prec.mean() for g in rows], C, fmt="{:.3f}",
         ylabel="Mean precision", title="M-CSA q250: precision of the returned set")
    bars(axes[1], labels, [g.rec.mean() for g in rows], C, fmt="{:.3f}",
         ylabel="Mean recall", title="Recall")
    ax = axes[2]
    for i, (g, c) in enumerate(zip(rows, C)):
        v = g.hits.to_numpy(dtype=float) + 0.5
        ax.bar(i, np.median(v), width=0.7, color=c, edgecolor="white", linewidth=0.8, zorder=3)
        ax.annotate(f"{np.median(g.hits):.0f}", (i, np.median(v)), textcoords="offset points",
                    xytext=(0, 2), ha="center", va="bottom", fontsize=8, color=MUTED)
    ax.set_yscale("log"); ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=8, rotation=30, ha="right")
    ax.set_ylabel("Median hits per query (log)"); ax.set_title("Hit list size", fontsize=10)
    fig.suptitle("--confident: coverage and RMSD cutoffs on the returned hit list",
                 fontsize=12.5, fontweight="bold", y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    save(fig, "fig_confident")


def save(fig, stem):
    os.makedirs(FIG, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(f"{FIG}/{stem}.{ext}", dpi=200)
    plt.close(fig)
    print("wrote", f"{FIG}/{stem}.png")


if __name__ == "__main__":
    pubr()
    fig_sort()
    fig_confident()

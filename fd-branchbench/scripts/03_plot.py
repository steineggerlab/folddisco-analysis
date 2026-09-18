#!/usr/bin/env python
"""Plot master vs feature-integration: motif P/R/F1, runtime, M-CSA Sens@1FP, deltas."""
import argparse, pathlib
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns

# Paper palette (folddisco pink is the highlight, gray is the 2.x master). Colour follows
# the config entity in every figure, so gray is always master wherever it appears.
PINK, GRAY = "#E94B8B", "#8E8E93"
COLORS = {"master": GRAY, "default": PINK, "sensitive": "#655AD0",
          "aa_blosum62": "#4FBDF2", "aa_group": "#B04BB6", "aa_size": "#F7B4A7",
          "sensitive_aa": "#3B3490"}
NAMES  = {"master": "master", "default": "branch default", "sensitive": "branch --sensitive",
          "aa_blosum62": "--aa-subst blosum62", "aa_group": "--aa-subst group",
          "aa_size": "--aa-subst size", "sensitive_aa": "--sensitive + blosum62"}
# Mutant benchmark (fig6) is its own entity set; the exact-search baseline takes the gray.
MUT = ["orig", "mut", "mut_star", "mut_blosum62", "mut_sens_star"]
MUT_COLORS = dict(zip(MUT, ["#4FBDF2", GRAY, PINK, "#655AD0", "#B04BB6"]))
MUT_NAMES = {"orig": "original query", "mut": "mutant, exact", "mut_star": "mutant, :* on it",
             "mut_blosum62": "mutant, --aa-subst blosum62", "mut_sens_star": "mutant, :* + --sensitive"}
# branch default is byte/rank-identical to master by design: dashed so the hidden line reads.
DASH   = {k: "-" for k in COLORS}; DASH["default"] = (0, (4, 2))
WIDTH  = {k: 2.0 for k in COLORS}; WIDTH["master"] = 3.0; WIDTH["default"] = 1.8
# fig1-4 answer "branch vs master" and stay at three series; fig5 is the aa-subst axis.
ORDER  = ["master", "default", "sensitive"]
AA     = ["default", "sensitive", "aa_blosum62", "aa_group", "aa_size", "sensitive_aa"]
INK, MUTED = "#0b0b0b", "#52514e"
FAMILY = {"zinc4": "Zinc finger (4 res)", "zinc3": "Zinc finger (3 res)",
          "serine": "Serine peptidase", "zincseg": "Zinc finger (segments)"}


def pubr():
    """ggpubr theme_pubr(): no grid, black axis lines, outward ticks."""
    sns.set_theme(style="ticks")
    plt.rcParams.update({
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
        "axes.grid": False, "axes.edgecolor": "black", "axes.linewidth": 0.8,
        "axes.spines.top": False, "axes.spines.right": False,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "axes.titlesize": 10, "axes.titleweight": "bold", "axes.labelsize": 10,
        "legend.frameon": False, "pdf.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Inter", "Liberation Sans", "Arial", "DejaVu Sans"],
    })


def legend_on(fig, keys=ORDER, y=0.005, kind="bar"):
    if kind in ("line", "dots"):
        h = [Line2D([], [], color=COLORS[k], lw=WIDTH[k], ls=DASH[k], ms=5.5, mec="white",
                    mew=0.9, marker="o" if kind == "dots" else "") for k in keys]
    else:
        h = [plt.Rectangle((0, 0), 1, 1, fc=COLORS[k], ec="none") for k in keys]
    fig.legend(h, [NAMES[k] for k in keys], loc="lower center", ncol=len(keys),
               bbox_to_anchor=(0.5, y), fontsize=10, handlelength=1.6, handleheight=1.0)


def bars(ax, sub, xcats, xcol, ycol, fmt, keys=ORDER, base=0.0):
    """Grouped bars: white spacer between adjacent fills, direct value labels."""
    n, w = len(keys), 0.8 / len(keys)
    for i, cfg in enumerate(keys):
        d = sub[sub.config == cfg].set_index(xcol).reindex(xcats)
        x = np.arange(len(xcats)) + (i - (n - 1) / 2) * w
        v = d[ycol].to_numpy(dtype=float)
        ax.bar(x, v - base, width=w - 0.02, bottom=base, color=COLORS[cfg],
               edgecolor="white", linewidth=0.8, zorder=3)
        for xi, vi in zip(x, v):
            if np.isfinite(vi):
                up = vi >= base
                ax.annotate(fmt(vi), (xi, vi), textcoords="offset points",
                            xytext=(0, 2 if up else -2), ha="center",
                            va="bottom" if up else "top", fontsize=6.2, color=MUTED, rotation=90)
    ax.set_xticks(np.arange(len(xcats))); ax.set_xticklabels(xcats)


def dots(ax, sub, xcats, xcol, ycol, fmt, keys=ORDER, pad=0.02):
    """Connected dots on a zoomed scale: for panels where bars differ too little to read."""
    x, vals = np.arange(len(xcats)), []
    for i, cfg in enumerate(keys):
        v = sub[sub.config == cfg].set_index(xcol).reindex(xcats)[ycol].to_numpy(dtype=float)
        vals.append(v)
        ax.plot(x, v, color=COLORS[cfg], lw=WIDTH[cfg], ls=DASH[cfg], marker="o",
                ms=7 if cfg == "master" else 5.5, mec="white", mew=0.9,
                solid_capstyle="round", zorder=3 + i)
    for xi in range(len(xcats)):
        col = sorted(((vals[i][xi], keys[i]) for i in range(len(keys))
                      if np.isfinite(vals[i][xi])), reverse=True)
        shown, place = {}, [(0, 8, "bottom"), (0, -9, "top"), (26, 0, "center")]
        for v, k in col:
            r = round(v, 3)
            if r in shown:                      # identical series (master == default): label once
                shown[r].append(k); continue
            shown[r] = [k]
        for j, (r, ks) in enumerate(shown.items()):
            dx, dy, va = place[min(j, 2)]
            ax.annotate(fmt(r), (xi, r), textcoords="offset points", xytext=(dx, dy),
                        ha="center" if not dx else "left", va=va, fontsize=6.4,
                        color=MUTED if len(ks) > 1 else COLORS[ks[0]], fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(xcats); ax.set_xlim(-0.4, len(xcats) - 0.6)
    lo = min(np.nanmin(v) for v in vals); hi = max(np.nanmax(v) for v in vals)
    ax.set_ylim(lo - pad, hi + pad)


def ecdf_steps(v, reverse=False):
    v = np.sort(np.asarray(v, dtype=float))
    if reverse:
        return np.concatenate([[0], v]), np.concatenate([[1], 1 - np.arange(len(v)) / len(v)])
    return v, np.arange(1, len(v) + 1) / len(v)


def fig_accuracy(m, out):
    """Precision / recall / F1 per motif query, master vs branch. All values sit in 0.82-0.98,
    so the panels are connected dots on a zoomed scale instead of bars from zero."""
    mt = m.melt(id_vars=["config", "qid", "label", "motif", "mode"],
                value_vars=["precision", "recall", "f1"], var_name="metric", value_name="value")
    mt["family"] = mt.label.str.split("_").str[0]
    fams, modes = ["zinc4", "zinc3", "serine", "zincseg"], ["prefilter", "matched"]
    fig, axes = plt.subplots(2, 4, figsize=(13, 6.4), sharey=True)
    for r, mode in enumerate(modes):
        for c, fam in enumerate(fams):
            ax = axes[r, c]
            sub = mt[(mt.family == fam) & (mt["mode"] == mode)]
            dots(ax, sub, ["precision", "recall", "f1"], "metric", "value", lambda v: f"{v:.3f}")
            ax.set_xticklabels(["Precision", "Recall", "F1"], fontsize=9)
            ax.set_ylim(0.80, 1.005); ax.set_yticks([0.80, 0.85, 0.90, 0.95, 1.0])
            ax.grid(axis="y", color="#e6e4e1", lw=0.7, zorder=0)
            qid = sub.qid.iloc[0] if len(sub) else ""
            ax.set_title(f"{FAMILY[fam]}\n{mode} ({qid})", fontsize=9.5)
            if c == 0: ax.set_ylabel("Metric value (axis starts at 0.80)")
    fig.suptitle("Motif search accuracy on the human index (23,391 structures)",
                 fontsize=13, fontweight="bold", y=0.995)
    fig.tight_layout(rect=[0, 0.07, 1, 0.965]); legend_on(fig, kind="dots")
    save(fig, out, "fig1_motif_accuracy")


def fig_runtime(m, mc, out, subset):
    """Query runtime: motif queries (log bars) + M-CSA distribution (ECDF)."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.0), gridspec_kw={"width_ratios": [1.45, 1]})
    ax = axes[0]
    qids = sorted(m.qid.unique())
    bars(ax, m, qids, "qid", "runtime_med", lambda v: f"{v*1000:.0f}ms" if v < 1 else f"{v:.2f}s")
    ax.set_yscale("log"); ax.set_ylabel("Median wall time (s, log)"); ax.set_xlabel("Motif query")
    ax.set_title(f"Motif queries (median of {int(m.n_repeat.iloc[0])} repeats)")
    lbl = m.drop_duplicates("qid").set_index("qid").reindex(qids)
    ax.set_xticklabels([f"{q}\n{lbl.loc[q,'label'].replace('_',chr(10))}" for q in qids], fontsize=7.5)
    ax.set_ylim(top=ax.get_ylim()[1] * 3.5)

    ax = axes[1]
    for i, cfg in enumerate(ORDER):
        v = mc[mc.config == cfg].runtime.dropna()
        if not len(v): continue
        x, y = ecdf_steps(v)
        ax.step(x, y, where="post", color=COLORS[cfg], lw=WIDTH[cfg], ls=DASH[cfg],
                solid_capstyle="round", zorder=3 + i)
        ax.annotate(f"{NAMES[cfg]}: median {np.median(v):.2f}s", (0.03, 0.94 - 0.075 * i),
                    xycoords="axes fraction", fontsize=8.5, color=COLORS[cfg], fontweight="bold")
    ax.set_xscale("log"); ax.set_xlabel("Wall time per query (s, log)")
    ax.set_ylabel("Fraction of M-CSA queries"); ax.set_ylim(0, 1.02)
    ax.set_title(f"M-CSA queries (n={mc['query'].nunique()}, --top 6000, 5 proc x 4 threads)")
    fig.suptitle("Query runtime", fontsize=13, fontweight="bold", y=0.995)
    fig.tight_layout(rect=[0, 0.085, 1, 0.945]); legend_on(fig, kind="line")
    save(fig, out, f"fig2_runtime_{subset}")


def fig_mcsa(mc, out, subset):
    """Sens@1FP on M-CSA: reverse ECDF, summary, paired delta."""
    fig, axes = plt.subplots(1, 3, figsize=(13, 5.0), gridspec_kw={"width_ratios": [1.25, .8, 1.05]})
    ax = axes[0]
    for i, cfg in enumerate(ORDER):
        v = mc[mc.config == cfg].sens_at_1fp.dropna()
        if not len(v): continue
        x, y = ecdf_steps(v, reverse=True)
        ax.step(x, y, where="post", color=COLORS[cfg], lw=WIDTH[cfg], ls=DASH[cfg],
                solid_capstyle="round", zorder=3 + i)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
    ax.set_xlabel("Sensitivity at 1st FP"); ax.set_ylabel("Fraction of M-CSA queries")
    ax.set_title("Reverse ECDF (higher/right is better)")

    ax = axes[1]
    s = mc.groupby("config").sens_at_1fp.agg(["mean", "median"]).reindex(ORDER).reset_index()
    sm = s.melt(id_vars="config", var_name="stat", value_name="value")
    bars(ax, sm, ["mean", "median"], "stat", "value", lambda v: f"{v:.4f}")
    ax.set_xticklabels(["Mean", "Median"]); ax.set_ylabel("Sensitivity at 1st FP")
    ax.set_ylim(0, sm.value.max() * 1.35); ax.set_title("Summary")

    ax = axes[2]
    w = mc.pivot_table(index="query", columns="config", values="sens_at_1fp")
    if {"master", "sensitive"} <= set(w.columns):
        d = (w["sensitive"] - w["master"]).dropna()
        nz = d[d != 0]
        ax.hist(nz, bins=31, color=COLORS["sensitive"], edgecolor="white", linewidth=0.5, zorder=3)
        ax.set_ylim(top=ax.get_ylim()[1] * 1.32)
        ax.axvline(0, color=MUTED, lw=1, ls="--", zorder=4)
        ax.set_title(f"--sensitive − master: win {(d>0).sum()} / loss {(d<0).sum()} / tie {(d==0).sum()}")
        ax.annotate(f"mean Δ {d.mean():+.4f}   median Δ {d.median():+.4f}\n"
                    f"({len(nz)} queries changed; ties excluded from bars)",
                    (0.03, 0.88), xycoords="axes fraction", fontsize=8.5, color=INK, fontweight="bold")
    ax.set_xlabel("Δ Sensitivity at 1st FP"); ax.set_ylabel("Queries")
    fig.suptitle(f"M-CSA catalytic-site benchmark, Sens@1FP "
                 f"(n={mc['query'].nunique()}, 62,122-entry PDB index)",
                 fontsize=13, fontweight="bold", y=0.995)
    fig.tight_layout(rect=[0, 0.085, 1, 0.945]); legend_on(fig, kind="line")
    save(fig, out, f"fig3_mcsa_sens_at_1fp_{subset}")


def fig_delta(m, mc, out, subset):
    """Gain over master: per-metric deltas, runtime ratio, M-CSA delta."""
    qids = sorted(m.qid.unique())
    base = m[m.config == "master"].set_index("qid")
    d = m[m.config.isin([c for c in ORDER if c != "master"])].copy()
    for c in ["precision", "recall", "f1", "sens_at_1fp"]:
        d["d_" + c] = d[c].to_numpy() - base.loc[d.qid, c].to_numpy()
    d["ratio"] = d.runtime_med.to_numpy() / base.loc[d.qid, "runtime_med"].to_numpy()
    keys = ["default", "sensitive"]

    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.4))
    panels = [("d_precision", "Δ Precision"), ("d_recall", "Δ Recall"),
              ("d_f1", "Δ F1"), ("d_sens_at_1fp", "Δ Sens@1FP")]
    for i, (col, lab) in enumerate(panels):
        ax = axes[i // 3, i % 3]
        bars(ax, d, qids, "qid", col, lambda v: f"{v:+.3f}", keys=keys)
        ax.axhline(0, color=INK, lw=0.9, zorder=4)
        ax.set_ylabel(lab); ax.set_title(f"{lab} vs master")
        ax.set_xticklabels(qids, fontsize=8.5)
        lim = (np.nanmax(np.abs(d[col].to_numpy())) or 0.005) * 2.0
        ax.set_ylim(-lim, lim)

    ax = axes[1, 1]
    bars(ax, d, qids, "qid", "ratio", lambda v: f"{v:.2f}x", keys=keys, base=1.0)
    ax.axhline(1, color=INK, lw=0.9, zorder=4)
    ax.set_ylabel("Runtime / master"); ax.set_title("Runtime ratio vs master (1.0 = no change)")
    ax.set_xticklabels(qids, fontsize=8.5)
    lo, hi = d.ratio.min(), d.ratio.max()
    ax.set_ylim(min(0.5, lo - 0.15), hi + 0.35)

    ax = axes[1, 2]
    w = mc.pivot_table(index="query", columns="config", values="sens_at_1fp")
    rows = [(k, (w[k] - w["master"]).dropna()) for k in keys if k in w.columns]
    x = np.arange(len(rows))
    for i, (k, dd) in enumerate(rows):
        ci = 1.96 * dd.std(ddof=1) / np.sqrt(len(dd))
        ax.bar(i, dd.mean(), width=0.5, color=COLORS[k], edgecolor="white", linewidth=0.8, zorder=3)
        ax.errorbar(i, dd.mean(), yerr=ci, color=INK, capsize=4, lw=1, zorder=4)
        ax.annotate(f"{dd.mean():+.4f}\nwin {(dd>0).sum()} / loss {(dd<0).sum()}",
                    (i, max(dd.mean(), 0) + ci), textcoords="offset points", xytext=(0, 3),
                    ha="center", va="bottom", fontsize=8, color=MUTED)
    ax.axhline(0, color=INK, lw=0.9, zorder=4)
    ax.set_xticks(x); ax.set_xticklabels([NAMES[k] for k in keys], fontsize=9)
    ax.set_ylabel("Mean Δ Sens@1FP"); ax.set_title(f"M-CSA mean Δ (n={mc['query'].nunique()}, 95% CI)")
    ax.set_ylim(top=max(1e-3, max(dd.mean() + 1.96 * dd.std(ddof=1) / np.sqrt(len(dd))
                                  for _, dd in rows) * 1.6))

    fig.suptitle("Gain over origin/master (motif commands: branch default is byte-identical)",
                 fontsize=13, fontweight="bold", y=0.995)
    fig.tight_layout(rect=[0, 0.06, 1, 0.955]); legend_on(fig, keys=keys)
    save(fig, out, f"fig4_gain_vs_master_{subset}")


def fig_aa(m, out, subset, mc_aa=None):
    """The amino-acid substitution axis: --sensitive is geometric, --aa-subst is chemical."""
    keys = [k for k in AA if k in set(m.config)]
    mt = m.melt(id_vars=["config", "qid"], value_vars=["precision", "recall", "f1"],
                var_name="metric", value_name="value")
    fig, axes = plt.subplots(2, 3, figsize=(14.5, 8.0))

    for ax, (qid, title) in zip(axes[0], [("C2", "Serine peptidase, matched (C2)"),
                                          ("C1", "Serine peptidase, prefilter (C1)"),
                                          ("A2", "Zinc finger 4-res, matched (A2)")]):
        bars(ax, mt[mt.qid == qid], ["precision", "recall", "f1"], "metric", "value",
             lambda v: f"{v:.3f}", keys=keys)
        ax.set_xticklabels(["Precision", "Recall", "F1"], fontsize=9)
        ax.set_ylim(0, 1.22); ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_ylabel("Metric value"); ax.set_title(title, fontsize=9.5)

    ax = axes[1, 0]
    qids = sorted(m.qid.unique())
    bars(ax, m, qids, "qid", "hits", lambda v: f"{v:.0f}", keys=keys)
    ax.set_yscale("log"); ax.set_ylim(bottom=50, top=m[m.config.isin(keys)].hits.max() * 9)
    ax.set_ylabel("Result rows (log)"); ax.set_xlabel("Motif query")
    ax.set_title("Hits", fontsize=9.5)
    ax.set_xticklabels(qids, fontsize=8.5)

    ax = axes[1, 1]
    bars(ax, m, qids, "qid", "runtime_med", lambda v: f"{v*1000:.0f}ms" if v < 1 else f"{v:.1f}s",
         keys=keys)
    ax.set_yscale("log"); ax.set_ylabel("Median wall time (s, log)"); ax.set_xlabel("Motif query")
    ax.set_title("Runtime cost", fontsize=9.5)
    ax.set_xticklabels(qids, fontsize=8.5); ax.set_ylim(top=ax.get_ylim()[1] * 5)

    ax = axes[1, 2]
    if mc_aa is not None and len(mc_aa):
        k2 = [k for k in keys if k in set(mc_aa.config)]
        sm = (mc_aa.groupby("config").sens_at_1fp.mean().reindex(k2).reset_index()
              .assign(stat="mean").rename(columns={"sens_at_1fp": "value"}))
        bars(ax, sm, ["mean"], "stat", "value", lambda v: f"{v:.4f}", keys=k2)
        ax.set_xticklabels([f"M-CSA q{mc_aa.attrs.get('subset','?')}  (n={mc_aa['query'].nunique()})"])
        ax.set_ylabel("Mean Sens@1FP"); ax.set_ylim(0, sm.value.max() * 1.45)
        ax.set_title("M-CSA mean Sens@1FP", fontsize=9.5)
    else:
        ax.axis("off")
        ax.annotate("M-CSA aa run not present\n(run 02_run_mcsa.sh)", (0.5, 0.5),
                    xycoords="axes fraction", ha="center", va="center", fontsize=10, color=MUTED)

    fig.suptitle("Amino-acid substitution vs geometric expansion  "
                 "(substituted residues score below exact ones)",
                 fontsize=12.5, fontweight="bold", y=0.995)
    fig.tight_layout(rect=[0, 0.075, 1, 0.955])
    legend_on(fig, keys=keys, y=0.005)
    save(fig, out, f"fig5_aa_substitution_{subset}")


def fig_mutant(mu, out, res, subset):
    """One residue renamed per M-CSA query: does substitution win back the answers?"""
    keys = [k for k in MUT if k in set(mu.config)]
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.8), gridspec_kw={"width_ratios": [1, 1.1, 1]})
    s = mu.groupby("config").agg(mean=("sens_at_1fp", "mean"), runtime=("runtime", "mean"))
    for ax, col, lab, fmt in [(axes[0], "mean", "Mean Sens@1FP", lambda v: f"{v:.4f}"),
                              (axes[2], "runtime", "Mean wall time per query (s)", lambda v: f"{v:.1f}s")]:
        for i, k in enumerate(keys):
            v = s.loc[k, col]
            ax.bar(i, v, width=0.7, color=MUT_COLORS[k], edgecolor="white", linewidth=0.8, zorder=3)
            ax.annotate(fmt(v), (i, v), textcoords="offset points", xytext=(0, 2),
                        ha="center", va="bottom", fontsize=8, color=MUTED)
        ax.set_xticks(range(len(keys))); ax.set_xticklabels([])
        ax.set_ylabel(lab); ax.set_ylim(0, s.loc[keys, col].max() * 1.2)
    axes[0].set_title("Sensitivity at 1st FP")
    axes[2].set_title("Runtime (5 proc x 4 threads)")

    ax = axes[1]
    w = mu.pivot_table(index="query", columns="config", values="sens_at_1fp")
    rows = [k for k in keys if k != "mut" and k in w.columns]
    for i, k in enumerate(rows):
        dd = (w[k] - w["mut"]).dropna()
        ci = 1.96 * dd.std(ddof=1) / np.sqrt(len(dd))
        ax.bar(i, dd.mean(), width=0.6, color=MUT_COLORS[k], edgecolor="white", linewidth=0.8, zorder=3)
        ax.errorbar(i, dd.mean(), yerr=ci, color=INK, capsize=4, lw=1, zorder=4)
        ax.annotate(f"{dd.mean():+.4f}\nwin {(dd>0).sum()} / loss {(dd<0).sum()}",
                    (i, max(dd.mean(), 0) + ci), textcoords="offset points", xytext=(0, 3),
                    ha="center", va="bottom", fontsize=8, color=MUTED)
    ax.axhline(0, color=INK, lw=0.9, zorder=4)
    ax.set_xticks(range(len(rows))); ax.set_xticklabels([])
    ax.set_ylabel("Mean Δ Sens@1FP vs mutant, exact")
    ax.set_title(f"Paired gain over exact search (95% CI)")
    ax.set_ylim(top=ax.get_ylim()[1] * 1.5)
    fig.suptitle(f"Mutant M-CSA queries: one residue renamed to its closest BLOSUM62 alternative "
                 f"(n={mu['query'].nunique()}, same answer sets)", fontsize=12, fontweight="bold", y=0.995)
    fig.tight_layout(rect=[0, 0.08, 1, 0.95])
    h = [plt.Rectangle((0, 0), 1, 1, fc=MUT_COLORS[k], ec="none") for k in keys]
    fig.legend(h, [MUT_NAMES[k] for k in keys], loc="lower center", ncol=len(keys),
               bbox_to_anchor=(0.5, 0.005), fontsize=9.5, handlelength=1.6, handleheight=1.0)
    save(fig, out, f"fig6_mutant_{subset}")
    s.reindex(keys).to_csv(res / f"mutant_summary_{subset}.tsv", sep="\t",
                           float_format="%.4f")


def save(fig, out, stem):
    for e in ("png", "pdf"):
        fig.savefig(out / f"{stem}.{e}", dpi=200)
    plt.close(fig)


def main():
    root = pathlib.Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--result", default=str(root / "result"))
    ap.add_argument("--figure", default=str(root / "figure"))
    ap.add_argument("--subset", default="250")
    ap.add_argument("--aa-subset", dest="aa_subset", default=None,
                    help="M-CSA subset carrying the aa-subst configs (fig5) [--subset]")
    a = ap.parse_args()
    res, out = pathlib.Path(a.result), pathlib.Path(a.figure); out.mkdir(parents=True, exist_ok=True)
    pubr()
    m = pd.read_csv(res / "motif_metrics.tsv", sep="\t")
    mc = pd.read_csv(res / f"mcsa_fp1_{a.subset}.tsv", sep="\t")
    fig_accuracy(m, out); fig_runtime(m, mc, out, a.subset)
    fig_mcsa(mc, out, a.subset); fig_delta(m, mc, out, a.subset)
    mc_aa = None
    aa_subset = a.aa_subset or a.subset
    aa_path = res / f"mcsa_fp1_{aa_subset}.tsv"
    if aa_path.is_file():
        mc_aa = pd.read_csv(aa_path, sep="\t"); mc_aa.attrs["subset"] = aa_subset
    fig_aa(m, out, a.subset, mc_aa)
    mu_path = res / f"mutant_fp1_{a.subset}.tsv"
    if mu_path.is_file():
        fig_mutant(pd.read_csv(mu_path, sep="\t"), out, res, a.subset)

    # Table view (relief rule: aqua sits below 3:1 contrast on the light surface).
    (mc.groupby("config").agg(n=("query", "nunique"), mean_sens=("sens_at_1fp", "mean"),
                              median_sens=("sens_at_1fp", "median"),
                              zero_tp=("tp", lambda s: (s == 0).sum()),
                              mean_runtime=("runtime", "mean"), median_runtime=("runtime", "median"))
       .reindex(ORDER).to_csv(res / f"mcsa_summary_{a.subset}.tsv", sep="\t", float_format="%.4f"))
    print("figures ->", out)
    for f in sorted(out.glob("*.png")): print("  ", f.name)


if __name__ == "__main__":
    main()

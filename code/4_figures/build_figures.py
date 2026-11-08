# -*- coding: utf-8 -*-
"""Build Figures 1 and 2 for the IJE v2 package directly from the MR output files.

Usage: python code/4_figures/build_figures.py results <out dir>

Within-family BMI (Howe et al. 2022, IEU OpenGWAS ieu-b-4815) is distributed on
an unstandardised (kg/m^2) scale; its IVW log-OR and SE are multiplied by
BMI_SD so that it shares the per-1-SD axis with the population BMI estimates.
P values are scale-invariant and are taken from the files unchanged.
"""
import sys, os, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullLocator, FixedFormatter

BMI_SD = 4.8          # kg/m^2, UK Biobank (Vithayathil et al. 2021, Table 1)
Z = 1.959963984540054

SRC, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7,
    "xtick.labelsize": 6.5, "ytick.labelsize": 7, "legend.fontsize": 6.5,
    "pdf.fonttype": 42, "ps.fonttype": 42, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "lines.linewidth": 0.9,
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold", "axes.spines.top": False, "axes.spines.right": False,
})
MM = 1 / 25.4
BLACK, GREEN, PURPLE, ORANGE, BLUE = "#000000", "#009E73", "#7B4F9E", "#D55E00", "#0072B2"

OUTCOMES = ["AD", "CKD", "COLORECTAL", "HF", "PROSTATE", "STR"]
MORT_LAB = {"AD": "Alzheimer's disease", "CKD": "Chronic kidney disease", "COLORECTAL": "Colorectal cancer",
            "HF": "Heart failure", "PROSTATE": "Prostate cancer", "STR": "Stroke"}
CC_LAB = dict(MORT_LAB, STR="Ischaemic stroke")

DIRS = {"mort": "mortality", "cc_bmi": "case_control", "cc_wf": "case_control", "cc_inc": "case_control"}
EXPS = {"cc_bmi": "BMI", "cc_wf": "BMI_within_family", "cc_inc": "Household"}


def load():
    rows = []
    for k, d in DIRS.items():
        pattern = os.path.join(SRC, d, EXPS.get(k, "*"), "*.or.txt")
        for f in glob.glob(pattern):
            exp, out = os.path.basename(f)[:-7].split("--")
            t = pd.read_csv(f, sep="\t")
            t = t[t.method == "Inverse variance weighted"].copy()
            t["set"], t["exp"], t["out"] = k, exp, out
            rows.append(t)
    d = pd.concat(rows, ignore_index=True)
    wf = d.exp == "BMI_within_family"
    d.loc[wf, "b"] *= BMI_SD
    d.loc[wf, "se"] *= BMI_SD
    d["or"], d["lo"], d["hi"] = np.exp(d.b), np.exp(d.b - Z * d.se), np.exp(d.b + Z * d.se)
    return d


D = load()


def get(set_, exp):
    t = D[(D.set == set_) & (D.exp == exp)].set_index("out")
    missing = [o for o in OUTCOMES if o not in t.index]
    if missing:
        raise SystemExit(f"missing {set_}/{exp}: {missing}")
    return t.loc[OUTCOMES]


def fmt_p(p):
    """P in manuscript style: 2 significant figures; a×10^b below 0.001."""
    if p >= 0.001:
        s = f"{p:.2g}"
        if "e" in s:
            s = f"{p:.3f}"
        if len(s.replace("0.", "").lstrip("0")) == 1:   # 0.04 -> 0.040
            s += "0"
        return f"P = {s}"
    m, e = f"{p:.1e}".split("e")
    return rf"P = {m}$\times$10$^{{{int(e)}}}$"


def log_axis(ax, lo, hi, ticks):
    ax.set_xscale("log")
    ax.set_xlim(lo, hi)
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(FixedFormatter([f"{t:g}" for t in ticks]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.axvline(1, color="0.45", lw=0.6, ls=(0, (3, 2)), zorder=0)


title_pad = 4


def forest(ax, series, labels, xlim, ticks, xlabel, title, annot_x, ylabels=True):
    """series: list of (table, colour, marker, legend label, filled)."""
    n = len(OUTCOMES)
    k = len(series)
    offs = np.linspace(0.17, -0.17, k) if k > 1 else [0.0]
    for (t, col, mk, lab, filled), off in zip(series, offs):
        y = np.arange(n)[::-1] + off
        for yi, (_, r) in zip(y, t.iterrows()):
            if r.lo < xlim[0] or r.hi > xlim[1]:
                raise SystemExit(f"CI outside axis in '{title}': {r.out} {r.lo:.3f}-{r.hi:.3f}")
            ax.plot([r.lo, r.hi], [yi, yi], color=col, lw=0.9, solid_capstyle="butt")
            ax.plot(r["or"], yi, marker=mk, ms=3.6, color=col, mfc=col if filled else "white",
                    mew=0.8, ls="none", label=lab if yi == y[0] else None)
            ax.text(annot_x, yi, f"{int(r.nsnp)} SNPs   {fmt_p(r.pval)}", transform=ax.get_yaxis_transform(),
                    ha="left", va="center", fontsize=6.2, color=col,
                    fontweight="bold" if r.pval < 0.05 else "normal")
    ax.set_yticks(np.arange(n)[::-1])
    ax.set_yticklabels([labels[o] for o in OUTCOMES] if ylabels else [])
    ax.set_ylim(-0.6, n - 0.4)
    ax.tick_params(axis="y", length=0)
    log_axis(ax, *xlim, ticks)
    ax.set_xlabel(xlabel)
    ax.set_title(title, loc="left", pad=title_pad)


def panel_letter(fig, ax, letter, dx=-0.0):
    bb = ax.get_position()
    fig.text(0.01 + dx, bb.y1 + 0.012, letter, fontsize=9, fontweight="bold", va="bottom")


def save(fig, name):
    for ext, kw in (("pdf", {}), ("png", {"dpi": 600}), ("tiff", {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}})):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), **kw)
    print("wrote", name)


BMI_X, BMI_T = (0.45, 2.8), [0.5, 0.7, 1, 1.5, 2, 2.5]
INC_X, INC_T = (0.2, 3.5), [0.25, 0.5, 1, 2, 3]

# ------------------------------------------------------------------ Figure 1
fig = plt.figure(figsize=(180 * MM, 150 * MM))
axA1 = fig.add_axes([0.17, 0.60, 0.20, 0.34])
axA2 = fig.add_axes([0.62, 0.60, 0.20, 0.34])
axB = fig.add_axes([0.10, 0.15, 0.86, 0.32])

forest(axA1, [(get("mort", "BMI"), BLACK, "o", None, True)], MORT_LAB, BMI_X, BMI_T,
       "IVW OR per 1-SD higher BMI (log scale)", "BMI (population GWAS)", 1.03)
forest(axA2, [(get("mort", "Household"), PURPLE, "s", None, True)], MORT_LAB, INC_X, INC_T,
       "IVW OR per 1-SD higher household income (log scale)", "Household income", 1.03, ylabels=False)
panel_letter(fig, axA1, "a")

STRATA = [("BMI", "Full set"), ("BMI_q1_lowest25", "Q1\n(strongest)"), ("BMI_q2_25to50", "Q2"),
          ("BMI_q3_50to75", "Q3"), ("BMI_q4_highest25", "Q4\n(weakest)")]
xpos = np.array([0, 1.5, 2.5, 3.5, 4.5])
for out, col, mk, ls, dx in (("HF", ORANGE, "o", "-", -0.07), ("CKD", BLUE, "s", (0, (4, 2)), 0.07)):
    t = D[(D.set == "mort") & (D.out == out)].set_index("exp").loc[[s for s, _ in STRATA]]
    x = xpos + dx
    axB.vlines(x, t.lo, t.hi, color=col, lw=0.9)
    axB.plot(x[1:], t["or"].values[1:], color=col, ls=ls, lw=0.9)
    axB.plot(x, t["or"], marker=mk, ms=4, color=col, ls="none", label=MORT_LAB[out] + " mortality")
    for xi, (_, r) in zip(x, t.iterrows()):
        row = -0.20 if out == "HF" else -0.27
        axB.text(xi - dx, row, f"{int(r.nsnp)} SNPs; {fmt_p(r.pval)}", transform=axB.get_xaxis_transform(),
                 ha="center", va="center", fontsize=6, color=col,
                 fontweight="bold" if r.pval < 0.05 else "normal")
    axB.text(-0.075, -0.20 if out == "HF" else -0.27, "HF" if out == "HF" else "CKD", transform=axB.transAxes,
             ha="left", va="center", fontsize=6, color=col, fontweight="bold")
lo = D[(D.set == "mort") & D.out.isin(["HF", "CKD"]) & D.exp.isin([s for s, _ in STRATA])]
assert lo.lo.min() > 0.5 and lo.hi.max() < 2.5, (lo.lo.min(), lo.hi.max())
axB.set_yscale("log")
axB.set_ylim(0.5, 2.5)
yt = [0.5, 0.7, 1, 1.5, 2, 2.5]
axB.yaxis.set_major_locator(FixedLocator(yt))
axB.yaxis.set_major_formatter(FixedFormatter([f"{t:g}" for t in yt]))
axB.yaxis.set_minor_locator(NullLocator())
axB.axhline(1, color="0.45", lw=0.6, ls=(0, (3, 2)), zorder=0)
axB.axvline(0.75, color="0.6", lw=0.5, ls=":")
axB.set_xticks(xpos)
axB.set_xticklabels([l for _, l in STRATA])
axB.set_xlim(-0.5, 5.0)
axB.set_ylabel("IVW OR per 1-SD higher BMI (log scale)")
axB.set_xlabel("BMI instrument set (quartiles of exposure-association P value)", labelpad=30)
axB.legend(loc="upper left", frameon=False, handlelength=1.5)
panel_letter(fig, axB, "b")
save(fig, "Figure1")
plt.close(fig)

# ------------------------------------------------------------------ Figure 2
title_pad = 15
fig = plt.figure(figsize=(180 * MM, 190 * MM))
axes = [fig.add_axes([0.17, y, 0.45, 0.215]) for y in (0.735, 0.41, 0.07)]
pop, wf = ("Population BMI (Yengo et al.)", "Within-family BMI (Howe et al.)")
forest(axes[0], [(get("cc_bmi", "BMI"), BLACK, "o", pop, True), (get("cc_wf", "BMI_within_family"), GREEN, "D", wf, False)],
       CC_LAB, BMI_X, BMI_T, "IVW OR per 1-SD higher BMI (log scale)",
       "BMI instruments and case-control GWAS of disease", 1.03)
forest(axes[1], [(get("mort", "BMI"), BLACK, "o", pop, True), (get("mort", "BMI_within_family"), GREEN, "D", wf, False)],
       MORT_LAB, BMI_X, BMI_T, "IVW OR per 1-SD higher BMI (log scale)",
       "BMI instruments and disease-specific mortality GWAS", 1.03)
forest(axes[2], [(get("cc_inc", "Household"), PURPLE, "s", None, True)], CC_LAB, INC_X, INC_T,
       "IVW OR per 1-SD higher household income (log scale)",
       "Household-income instruments and case-control GWAS of disease", 1.03)
for ax in axes[:2]:
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2, frameon=False, handlelength=1.2,
              columnspacing=1.2, borderaxespad=0)
for ax, L in zip(axes, "abc"):
    panel_letter(fig, ax, L)
save(fig, "Figure2")
plt.close(fig)

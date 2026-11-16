r"""
67_fig78.py -- Fig. 7 (cross-platform) and Fig. 8 (section replication).

Same specification as 43_figures_v2.py: 180 mm double column, Arial 7-9 pt, no
grid, vector PDF with a frozen timestamp plus a 300 dpi PNG. Every number is
read from NUMBERS_LOCK_v5.json; the point clouds come from the npz files that
64_crossplatform_robust.py wrote, i.e. the same arrays that produced the
reported correlations.

用法: python 67_fig78.py --outdir figs_v3
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullLocator

MM = 1 / 25.4
W2, W1 = 180 * MM, 90 * MM
COL = {"colon": "#1B3A5C", "lung": "#4A6E8F", "brain": "#B5651D",
       "stereo": "#6A4C93"}
GREY, LGREY, MGREY = "#8A8A8A", "#D6D6D6", "#4D4D4D"
DS = [2, 4, 8, 16, 32, 64, 128, 256]

plt.rcParams.update({
    "font.family": "Arial", "font.size": 7.5, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "xtick.direction": "out", "ytick.direction": "out",
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    "legend.frameon": False, "legend.handlelength": 1.6,
    "legend.borderpad": 0.2, "legend.labelspacing": 0.25,
})

ROOT = Path(__file__).resolve().parents[1]
N = json.loads((ROOT / "02_OUT/NUMBERS_LOCK_v5.json").read_text(encoding="utf-8"))
X = N["cross"]
MU = N["multi"]
CP = ROOT / "02_OUT/crossplatform"


def style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)
    ax.tick_params(which="both", top=False, right=False)


def panel(ax, s, dx=-0.16, dy=1.03):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="bottom", ha="left")


def save(fig, od, name):
    fig.savefig(od / f"{name}.pdf", metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(od / f"{name}.png", dpi=300)
    plt.close(fig)
    print("saved", name)


def xt(ax, vals=DS):
    ax.set_xscale("log", base=2)
    ax.set_xticks(vals)
    ax.set_xticklabels([f"{v:g}" for v in vals])
    ax.xaxis.set_minor_locator(NullLocator())


# ============================================================ Fig 7
def fig7(od):
    fig = plt.figure(figsize=(W2, 62 * MM))
    gs = fig.add_gridspec(1, 3, wspace=0.42, left=0.07, right=0.99, bottom=0.17, top=0.90)

    # ---- (a) normalised pooled alpha, Stereo vs Visium brain, colon reference
    ax = fig.add_subplot(gs[0, 0])
    a_st = np.array(X["stereo"]["alpha"]) / X["stereo"]["alpha"][0]
    a_vi = np.array(X["visium"]["alpha"]) / X["visium"]["alpha"][0]
    a_co = np.array([N["datasets"]["colon"]["alpha"][str(d)] for d in DS], float)
    a_co = a_co / a_co[0]
    ax.plot(DS, a_co, "-", color=COL["colon"], lw=1.0, label="Colon, Visium HD")
    ax.plot(DS, a_vi, "-", color=COL["brain"], lw=1.0, label="Brain, Visium HD")
    ax.plot(DS, a_st, "--", color=COL["stereo"], lw=1.2,
            label="Brain, Stereo-seq")
    ref = np.array(DS, float)
    ax.plot(DS, (ref / 2.0) ** -2, ":", color=GREY, lw=0.8)
    ax.text(90, 0.02, r"$d^{-2}$", color=GREY, fontsize=7, ha="left", va="bottom")
    ax.set_xlabel("Bin side $d$ (\u00b5m)")
    ax.set_ylabel(r"$\alpha_{\rm pool}(d)\,/\,\alpha_{\rm pool}(2\,\mu{\rm m})$")
    ax.set_yscale("log")
    xt(ax)
    ax.set_ylim(4e-3, 2.0)
    style(ax)
    ax.legend(loc="lower left", fontsize=6.2, handlelength=1.8)
    panel(ax, "a")

    # ---- (b) cumulative share of the alpha contribution against gene rank
    ax = fig.add_subplot(gs[0, 1])
    top100 = {}
    for tag, col, lab, ls in (
        ("vbar_stereo_primary", COL["stereo"], "Brain, Stereo-seq", "--"),
        ("vbar_visium", COL["brain"], "Brain, Visium HD", "-"),
    ):
        z = np.load(CP / f"{tag}.npz", allow_pickle=True)
        v = np.sort(z["v_bar"])[::-1]
        v = v[v > 0]
        cum = np.cumsum(v) / v.sum()
        top100[tag] = float(cum[99])
        ax.plot(np.arange(1, len(cum) + 1), cum, ls, color=col, lw=1.0, label=lab)
    ax.set_xscale("log")
    ax.set_xlabel("Genes, ranked by contribution")
    ax.set_ylabel("Cumulative share of $\\alpha_{\\rm pool}$")
    ax.set_xlim(0.9, 2e4)
    ax.set_ylim(0, 1.02)
    ax.text(0.04, 0.96,
            "Top 100 genes\n"
            f"Stereo-seq {top100['vbar_stereo_primary']*100:.0f}%\n"
            f"Visium HD {top100['vbar_visium']*100:.0f}%",
            transform=ax.transAxes, fontsize=6.4, va="top", ha="left")
    style(ax)
    ax.legend(loc="lower right", fontsize=6.2)
    panel(ax, "b")

    # ---- (c) per-gene exponent, one platform against the other
    ax = fig.add_subplot(gs[0, 2])
    z = np.load(CP / "pergene_shared.npz", allow_pickle=True)
    ga, gb = z["gamma_stereo"], z["gamma_visium"]
    hb = ax.hexbin(ga, gb, gridsize=48, bins="log", cmap="Blues",
                   mincnt=1, linewidths=0, extent=(-0.5, 2.5, -0.5, 2.5))
    lim = (-0.5, 2.5)
    ax.plot(lim, lim, "-", color=LGREY, lw=0.8)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("Per-gene $\\gamma$, Stereo-seq")
    ax.set_ylabel("Per-gene $\\gamma$, Visium HD")
    style(ax)
    cb = fig.colorbar(hb, ax=ax, fraction=0.045, pad=0.02)
    cb.set_label("Genes per cell", fontsize=6.5)
    cb.ax.tick_params(labelsize=6, width=0.5, length=2)
    cb.outline.set_linewidth(0.4)
    ax.text(0.04, 0.96,
            f"$\\rho$ = {X['replication']['rho_raw']:.2f}\n"
            f"$\\rho\\,|\\,$abundance = {X['replication']['rho_partial_abundance']:.2f}\n"
            f"n = {X['replication']['n']:,}",
            transform=ax.transAxes, fontsize=6.5, va="top", ha="left")
    panel(ax, "c")

    save(fig, od, "Fig7")


# ============================================================ Fig 8
def fig8(od):
    fig = plt.figure(figsize=(W2, 58 * MM))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.12, 0.95], wspace=0.46,
                          left=0.085, right=0.99, bottom=0.19, top=0.90)

    order = [("Visium_HD_Human_Colon_Cancer_P1", "Colon P1", COL["colon"], "o"),
             ("colon", "Colon P0", COL["colon"], "s"),
             ("lung", "Lung", COL["lung"], "o"),
             ("brain", "Brain", COL["brain"], "o")]

    # ---- (a) forest plot of gamma_pool with the superblock bootstrap CI
    ax = fig.add_subplot(gs[0, 0])
    ypos = np.arange(len(order))[::-1]
    for y, (k, lab, col, mk) in zip(ypos, order):
        s = MU["sections"][k]
        ax.plot([s["boot_lo"], s["boot_hi"]], [y, y], "-", color=col, lw=1.2)
        ax.plot([s["gamma_pool"]], [y], mk, ms=3.6, color=col, mfc=col)
    ax.set_yticks(ypos)
    ax.set_yticklabels([o[1] for o in order], fontsize=6.8)
    ax.set_xlabel(r"Pooled exponent $\gamma_{\rm pool}$")
    ax.set_ylim(-0.7, len(order) - 0.3)
    ax.set_xlim(0.335, 0.625)
    style(ax)
    ax.tick_params(axis="y", length=0)
    # the same-tissue difference and the tissue gap, annotated where no bar sits
    c1 = MU["contrasts"]["same_tissue_colon"]
    c2 = MU["contrasts"]["colon_vs_brain"]
    ax.text(0.985, 0.03,
            f"P1 $-$ P0: $\\Delta$ = {abs(c1['delta']):.3f}, $z$ = {abs(c1['z']):.1f}\n"
            f"colon $-$ brain: $\\Delta$ = {c2['delta']:.3f}, $z$ = {c2['z']:.1f}\n"
            f"MDE (2 sections, 80%): {MU['mde']['value']:.3f}",
            transform=ax.transAxes, fontsize=6.4, va="bottom", ha="right")
    panel(ax, "a")

    # ---- (b) the eight-point curves behind each estimate
    ax = fig.add_subplot(gs[0, 1])
    pool = json.loads((ROOT / "02_OUT/multisection/pool.json").read_text(encoding="utf-8"))
    rows = {r["label"]: r for r in pool["rows"]}
    for k, lab, col, mk in order:
        a = np.array(rows[k]["alpha"], float)
        ls = "--" if k.endswith("P1") else "-"
        ax.plot(DS, a, ls, color=col, lw=1.0, marker=mk, ms=2.0, mfc="none",
                label=f"{lab} ({rows[k]['gamma_pool']:.3f})")
    ax.set_xlabel("Bin side $d$ (\u00b5m)")
    ax.set_ylabel(r"$\alpha_{\rm pool}(d)$")
    ax.set_yscale("log")
    xt(ax)
    style(ax)
    ax.legend(loc="upper right", fontsize=6.2, ncol=1, handlelength=1.8)
    panel(ax, "b")

    # ---- (c) median per-gene exponent by abundance decile, both colon sections
    ax = fig.add_subplot(gs[0, 2])
    pg = json.loads((ROOT / "02_OUT/multisection/pergene_shift.json").read_text(encoding="utf-8"))
    dec = pg["deciles"]
    x = [d["decile"] for d in dec]
    m0 = [d["med0"] for d in dec]
    m1 = [d["med1"] for d in dec]
    ax.plot(x, m0, "-s", color=COL["colon"], ms=2.4, lw=1.0, label="Colon P0")
    ax.plot(x, m1, "--o", color=COL["stereo"], ms=2.4, lw=1.0, mfc="none", label="Colon P1")
    # the two medians over all genes, which is what Table 2 would report
    xa = max(x) + 1.15
    ax.plot([xa, xa], [MU["pergene"]["median_p0"], MU["pergene"]["median_p1"]],
            "-", color=MGREY, lw=0.8)
    ax.plot([xa], [MU["pergene"]["median_p0"]], "s", color=COL["colon"], ms=3.0)
    ax.plot([xa], [MU["pergene"]["median_p1"]], "o", color=COL["stereo"], ms=3.0, mfc="none")
    ax.annotate("", xy=(xa, MU["pergene"]["median_p1"]), xytext=(xa, MU["pergene"]["median_p0"]),
                arrowprops=dict(arrowstyle="<->", lw=0.7, color=MGREY, shrinkA=3, shrinkB=3))
    ax.text(xa - 0.35, (MU["pergene"]["median_p0"] + MU["pergene"]["median_p1"]) / 2,
            f"$\\Delta$ = {MU['pergene']['shift_overall']:.3f}\nall genes",
            fontsize=6.0, ha="right", va="center")
    ax.set_xlabel("Abundance decile (1, lowest)")
    ax.set_ylabel("Median per-gene $\\gamma$")
    ticks = x[::2] + [xa]
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(v) for v in x[::2]] + ["all"], fontsize=6.6)
    ax.set_xlim(0.4, xa + 0.5)
    ax.set_ylim(0.17, 0.99)
    ax.xaxis.set_minor_locator(NullLocator())
    style(ax)
    ax.legend(loc="upper right", fontsize=6.4)
    within = max(abs(d["shift"]) for d in dec)
    ax.text(0.03, 0.05, f"within a decile: max $|\\Delta|$ = {within:.3f}",
            transform=ax.transAxes, fontsize=6.2, va="bottom", ha="left")
    panel(ax, "c")

    save(fig, od, "Fig8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="figs_v3")
    a = ap.parse_args()
    od = ROOT / "02_OUT" / a.outdir if not Path(a.outdir).is_absolute() else Path(a.outdir)
    od.mkdir(parents=True, exist_ok=True)
    fig7(od)
    fig8(od)
    for p in sorted(od.glob("Fig[78].*")):
        print(f"  {p.name}  {p.stat().st_size/1024:.0f} KB")


if __name__ == "__main__":
    main()

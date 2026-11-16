r"""
43_figures_v2.py — 六张主图 + 五张补充图，全部现场从 NUMBERS_LOCK_v4.json 取数。

规格沿用 v1：双栏 180 mm / 单栏 90 mm；Arial 7–9 pt；无网格、无 3D、无阴影、无图例边框；
矢量 PDF（冻结时间戳）+ 300 dpi PNG。
配色：#1B3A5C colon、#4A6E8F lung、#B5651D brain、#8A8A8A 辅助灰。

点云来自 41_gene_curves.py 的 npz（数据，不是数字），地图来自 40_layer_maps.py。

用法：python 43_figures_v2.py --outdir figs_v2
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.ticker import NullLocator

MM = 1 / 25.4
W2, W1 = 180 * MM, 90 * MM
COL = {"colon": "#1B3A5C", "lung": "#4A6E8F", "brain": "#B5651D"}
LAB = {"colon": "Colon", "lung": "Lung", "brain": "Brain"}
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
N = json.loads((ROOT / "02_OUT/NUMBERS_LOCK_v4.json").read_text(encoding="utf-8"))
D = N["datasets"]
GC = {t: np.load(ROOT / f"02_OUT/gene_curves/gene_curves_{t}.npz", allow_pickle=True)
      for t in ("colon", "lung", "brain")}


def style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)
    ax.tick_params(which="both", top=False, right=False)


def panel(ax, s, dx=-0.16, dy=1.03):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="bottom", ha="left")


# 含几百万散点的图用低一点的位图 dpi，矢量 PDF 仍然完整
DPI = {"FigS5": 200}


def save(fig, od, name):
    fig.savefig(od / f"{name}.pdf", metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(od / f"{name}.png", dpi=DPI.get(name, 300))
    plt.close(fig)
    print("saved", name)


def alpha_of(t, d):
    return D[t]["alpha"][str(d)]


def xt(ax, vals, labels=None):
    ax.set_xticks(vals)
    ax.set_xticklabels(labels or [f"{v:g}" for v in vals])
    ax.xaxis.set_minor_locator(NullLocator())


# ============================================================ Fig 1
def fig1(od):
    fig = plt.figure(figsize=(W2, 74 * MM))
    axs = fig.add_axes([0.075, 0.16, 0.245, 0.74])
    ax = fig.add_axes([0.42, 0.16, 0.545, 0.74])

    # (a) 精确重聚合示意
    lab_d = 32
    for i in range(lab_d):
        for j in range(lab_d):
            axs.add_patch(Rectangle((i, j), 1, 1, facecolor="#F2F5F8", edgecolor=LGREY, lw=0.15))
    for i in range(0, lab_d, 4):
        for j in range(0, lab_d, 4):
            fill = "#C9D6E2" if (i // 4 + j // 4) % 2 == 0 else "#DCE6EF"
            axs.add_patch(Rectangle((i, j), 4, 4, facecolor=fill, edgecolor="none", lw=0))
    for i in range(0, lab_d, 4):
        axs.plot([i, i], [0, lab_d], color="#8FA6BB", lw=0.35, zorder=3)
        axs.plot([0, lab_d], [i, i], color="#8FA6BB", lw=0.35, zorder=3)
    axs.add_patch(Rectangle((8, 8), 16, 16, facecolor="none", edgecolor=COL["colon"], lw=1.3, zorder=4))
    axs.annotate("", xy=(1, -2.2), xytext=(3, -2.2),
                 arrowprops=dict(arrowstyle="<->", color=MGREY, lw=0.5))
    axs.text(2, -5.2, "2 µm", ha="center", va="top", fontsize=6.5, color=MGREY)
    axs.annotate("", xy=(4, 34.8), xytext=(8, 34.8),
                 arrowprops=dict(arrowstyle="<->", color=MGREY, lw=0.5))
    axs.text(6, 35.6, "4 µm", ha="center", va="bottom", fontsize=6.5, color=MGREY)
    axs.annotate("", xy=(8, -2.2), xytext=(24, -2.2),
                 arrowprops=dict(arrowstyle="<->", color=COL["colon"], lw=0.6))
    axs.text(16, -5.2, "16 µm", ha="center", va="top", fontsize=6.5, color=COL["colon"])
    axs.set_xlim(-1, lab_d + 1)
    axs.set_ylim(-8, 38)
    axs.axis("off")
    axs.text(-1, 36.5, "a", fontsize=9, fontweight="bold", va="bottom")
    axs.text(16, 42.5, "exact integer re-aggregation\n8 scales, no interpolation",
             ha="center", va="bottom", fontsize=6.8, color=MGREY)

    # (b) 归一化 alpha(d)
    dd = np.logspace(np.log10(1.9), np.log10(270), 60)
    ax.plot(dd, (dd / 2.0) ** -2, color=GREY, lw=0.9, ls=(0, (4, 2)), zorder=1)
    ax.plot(dd, np.ones_like(dd), color=GREY, lw=0.8, ls=(0, (1, 2)), zorder=1)
    for t in ("colon", "lung", "brain"):
        a = np.array([alpha_of(t, d) for d in DS])
        ax.plot(DS, a / a[0], color=COL[t], lw=1.2, marker="o", ms=2.6,
                mfc="white", mew=0.8, zorder=3, label=LAB[t])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1.9, 270)
    ax.set_ylim(3.5e-5, 3.4)
    xt(ax, [2, 8, 32, 128])
    ax.set_yticks([1e-4, 1e-3, 1e-2, 1e-1, 1])
    ax.set_yticklabels(["$10^{-4}$", "$10^{-3}$", "$10^{-2}$", "$10^{-1}$", "1"])
    ax.yaxis.set_minor_locator(NullLocator())
    ax.set_xlabel("Bin size $d$ (µm)")
    ax.set_ylabel(r"$\alpha(d)\,/\,\alpha(2\,\mathrm{µm})$")
    style(ax)
    ax.legend(loc="lower left", fontsize=6.8, ncol=3, columnspacing=0.9,
              bbox_to_anchor=(0.0, 1.0))
    ax.annotate("counting statistics\n" + r"$\alpha \propto d^{-2}$",
                xy=(120, (120 / 2) ** -2), xytext=(26, 2.6e-3), fontsize=6.8,
                color=MGREY, va="top",
                arrowprops=dict(arrowstyle="-", color=GREY, lw=0.5, shrinkA=2, shrinkB=2))
    ax.text(60, 1.35, r"fully correlated  $\alpha \propto d^{0}$", fontsize=6.8, color=MGREY)
    panel(ax, "b")
    save(fig, od, "Fig1")


# ============================================================ Fig 2
def fig2(od):
    fig = plt.figure(figsize=(W2, 62 * MM))
    ax1 = fig.add_axes([0.065, 0.18, 0.27, 0.72])
    ax2 = fig.add_axes([0.42, 0.18, 0.235, 0.72])
    ax3 = fig.add_axes([0.735, 0.18, 0.245, 0.72])

    # (a) 局部斜率
    for t in ("colon", "lung", "brain"):
        ls = N["models"][t]["local_slope"]
        x = np.array(ls["seg_mid"], float)
        ax1.plot(x, ls["observed"], color=COL[t], lw=1.0, marker="o", ms=2.4, zorder=4)
        ax1.plot(x, ls["predicted_quad"], color=COL[t], lw=0.7, ls=(0, (3, 2)), alpha=0.9)
        ax1.plot(x, ls["predicted_xi"], color=COL[t], lw=0.7, ls=(0, (1, 1.6)), alpha=0.75)
    ax1.set_xscale("log")
    ax1.set_xlim(2.2, 260)
    ax1.set_ylim(-0.05, 1.75)
    xt(ax1, [2, 8, 32, 128])
    ax1.set_xlabel("Segment midpoint (µm)")
    ax1.set_ylabel(r"Local slope  $\gamma(d)$")
    ax1.axhline(2, color=GREY, lw=0.6, ls=":")
    ax1.text(2.6, 1.82, r"$\gamma=2$", fontsize=6.5, color=GREY)
    ax1.text(2.6, 1.60, "finite-$\\xi$ (dotted) · quadratic (dashed) · observed (solid)",
             fontsize=6.2, color=MGREY)
    style(ax1)
    panel(ax1, "a")

    # (b) dAIC
    mods = ["power", "quad", "exp", "gauss"]
    mlab = ["power law", "quadratic", "exponential $\\xi$", "Gaussian $\\xi$"]
    ypos = np.arange(len(mods))
    for k, t in enumerate(("colon", "lung", "brain")):
        v = [N["models"][t]["daic"][m] for m in mods]
        ax2.barh(ypos + (k - 1) * 0.26, v, height=0.24, color=COL[t], edgecolor="none")
    ax2.set_yticks(ypos)
    ax2.set_yticklabels(mlab)
    ax2.invert_yaxis()
    ax2.set_xlabel(r"$\Delta$AIC vs best model")
    ax2.set_xlim(0, 72)
    ax2.axvline(0, color=MGREY, lw=0.6)
    style(ax2)
    panel(ax2, "b", dx=-0.42)

    # (c) xi profile
    for t in ("colon", "lung", "brain"):
        pr = N["models"][t]["xi_profile"]
        g = np.array(pr["xi_grid"], float)
        r = np.array(pr["rss"], float)
        ax3.plot(g, r - pr["rss_min"], color=COL[t], lw=1.0)
        ax3.axvline(pr["xi_hat"], color=COL[t], lw=0.6, ls=":")
        ax3.plot([pr["ci_lo"], pr["ci_hi"]], [3.2e-3, 3.2e-3], color=COL[t], lw=2.2,
                 solid_capstyle="butt", clip_on=False)
    ax3.axhline(np.array(N["models"]["colon"]["xi_profile"]["threshold"])
                - N["models"]["colon"]["xi_profile"]["rss_min"],
                color=GREY, lw=0.6, ls="--")
    ax3.text(1.15, 2.4, "95% profile threshold", fontsize=6.3, color=GREY)
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax3.set_xlim(1, 700)
    ax3.set_ylim(2.0e-3, 16)
    ax3.set_yticks([1e-2, 1e-1, 1, 10])
    ax3.set_yticklabels(["$10^{-2}$", "$10^{-1}$", "1", "10"])
    ax3.yaxis.set_minor_locator(NullLocator())
    ax3.set_xlabel(r"Correlation length $\xi$ (µm)")
    ax3.set_ylabel(r"$\Delta$RSS")
    xt(ax3, [1, 10, 100])
    style(ax3)
    panel(ax3, "c")
    save(fig, od, "Fig2")


# ============================================================ Fig 3
def fig3(od):
    fig = plt.figure(figsize=(W2, 148 * MM))
    ax1 = fig.add_axes([0.075, 0.575, 0.40, 0.375])
    ax2 = fig.add_axes([0.565, 0.575, 0.375, 0.375])
    ax3 = fig.add_axes([0.075, 0.075, 0.375, 0.375])
    ax4 = fig.add_axes([0.565, 0.075, 0.36, 0.375])

    # (a) gamma_g vs mu32
    for t in ("colon", "lung", "brain"):
        z = GC[t]
        m = z["fitted"]
        mu = z["mu32"][m]
        g = z["gamma"][m]
        ax1.scatter(mu, g, s=0.7, color=COL[t], alpha=0.16, edgecolors="none", rasterized=True)
    for t in ("colon", "lung", "brain"):
        z = GC[t]
        m = z["fitted"]
        mu = z["mu32"][m]
        g = z["gamma"][m]
        edges = np.geomspace(0.25, 200, 26)
        idx = np.digitize(mu, edges)
        cx, cy = [], []
        for b in range(1, len(edges)):
            s = idx == b
            if s.sum() >= 30:
                cx.append(np.sqrt(edges[b - 1] * edges[b]))
                cy.append(np.median(g[s]))
        ax1.plot(cx, cy, color=COL[t], lw=1.5, zorder=5)
        ax1.plot([], [], color=COL[t], lw=1.5, label=LAB[t])
    ax1.set_xscale("log")
    ax1.set_xlim(0.25, 200)
    ax1.set_ylim(0, 1.5)
    ax1.set_xlabel(r"Mean counts per 32 µm bin  $\mu_{32}$")
    ax1.set_ylabel(r"Per-gene exponent  $\gamma_g$")
    ax1.legend(loc="upper right", fontsize=6.8)
    ax1.text(0.30, 1.46, "lines: binned median", fontsize=6.2, color=MGREY)
    style(ax1)
    panel(ax1, "a")

    # (b) 分布
    data = [GC[t]["gamma"][GC[t]["fitted"]] for t in ("colon", "lung", "brain")]
    parts = ax2.violinplot(data, positions=[0, 1, 2], widths=0.72, showextrema=False)
    for pc, t in zip(parts["bodies"], ("colon", "lung", "brain")):
        pc.set_facecolor(COL[t])
        pc.set_alpha(0.35)
        pc.set_edgecolor(COL[t])
        pc.set_linewidth(0.6)
    for i, t in enumerate(("colon", "lung", "brain")):
        g = GC[t]["gamma"][GC[t]["fitted"]]
        q25, med, q75 = np.percentile(g, [25, 50, 75])
        ax2.plot([i, i], [q25, q75], color=COL[t], lw=3.0, solid_capstyle="butt")
        ax2.plot(i, med, "o", color="white", mec=COL[t], ms=3.4, mew=0.9, zorder=5)
        pooled = N["gene_level"][t]["gamma_gene"]["pooled"]
        ax2.plot([i - 0.24, i + 0.24], [pooled, pooled], color=MGREY, lw=0.9, ls=(0, (3, 2)))
    ax2.set_xticks([0, 1, 2])
    ax2.set_xticklabels([LAB[t] for t in ("colon", "lung", "brain")])
    ax2.set_ylim(0, 1.6)
    ax2.set_ylabel(r"$\gamma_g$")
    ax2.text(1.62, 1.50, "dashed: pooled $\\gamma$", fontsize=6.3, color=MGREY)
    ax2.set_xlim(-0.6, 2.6)
    style(ax2)
    panel(ax2, "b", dx=-0.14)

    # (c) 两半可靠性：真实逐基因点
    rel = N["splithalf"]["reliability"]
    zz = np.load(ROOT / "02_OUT/splithalf_genes/splithalf_genes.npz", allow_pickle=True)
    both = zz["both"]
    ga = zz["gamma_a"][both].astype(float)
    gb = zz["gamma_b"][both].astype(float)
    ax3.scatter(ga, gb, s=0.7, color=COL["colon"], alpha=0.18, edgecolors="none",
                rasterized=True)
    ax3.plot([0.02, 2], [0.02, 2], color=GREY, lw=0.6, ls=":")
    ax3.set_xlim(0.05, 1.8)
    ax3.set_ylim(0.05, 1.8)
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax3.set_xlabel("$\\gamma_g$, first half")
    ax3.set_ylabel("$\\gamma_g$, second half")
    ax3.text(0.065, 1.35, f"r = {rel['pearson']:.2f}   n = {rel['n']:,}\nMAD = {rel['mad']:.3f}",
             fontsize=7, color=COL["colon"])
    ax3.text(0.55, 0.065, "dotted: identity", fontsize=6.3, color=GREY)
    style(ax3)
    panel(ax3, "c")

    # (d) 按数量级的中位 gamma + 合成对照
    dec_edges = np.array([0, 0.01, 0.1, 1, 10, 100])
    for t in ("colon", "lung", "brain"):
        by = N["jackknife"][t]["gamma_vs_mu"]["by_decade"]
        xs, ys = [], []
        for b in by:
            lo, hi = b["lo"], b["hi"]
            if b["med_gamma"] is None:
                continue
            xs.append(np.sqrt(max(lo, 1e-3) * min(hi, 100)))
            ys.append(b["med_gamma"])
        ax4.plot(xs, ys, color=COL[t], lw=1.2, marker="s", ms=2.6, mfc="white", mew=0.7)
    sy = N["synthetic"]["expression"]["by_layer"]
    sx = [np.sqrt(float(k.split("-")[0]) * float(k.split("-")[1])) for k in sy]
    vv = [sy[k]["med"] for k in sy]
    ax4.plot(sx, vv, color=MGREY, lw=1.0, ls=(0, (3, 2)), marker="^", ms=2.6, mfc="white", mew=0.7)
    ax4.axhline(N["synthetic"]["expression"]["beta_true"], color=MGREY, lw=0.6, ls=":")
    ax4.text(0.012, N["synthetic"]["expression"]["beta_true"] + 0.03,
             "synthetic truth 0.55", fontsize=6.3, color=MGREY)
    ax4.set_xscale("log")
    ax4.set_xlim(0.008, 200)
    ax4.set_ylim(0, 1.25)
    ax4.set_xlabel(r"$\mu_{32}$ (counts per bin)")
    ax4.set_ylabel(r"Median $\gamma$")
    ax4.text(0.012, 1.16, "real (solid) · synthetic (dashed)", fontsize=6.3, color=MGREY)
    style(ax4)
    panel(ax4, "d", dx=-0.14)
    save(fig, od, "Fig3")


# ============================================================ Fig 4
def fig4(od):
    fig = plt.figure(figsize=(W2, 60 * MM))
    ax1 = fig.add_axes([0.065, 0.19, 0.245, 0.71])
    ax2 = fig.add_axes([0.405, 0.19, 0.245, 0.71])
    ax3 = fig.add_axes([0.775, 0.19, 0.21, 0.71])

    # (a) 权重累积
    for t in ("colon", "lung", "brain"):
        z = GC[t]
        mu = z["mu32"][z["fitted"]]
        w = mu ** 2
        w = np.sort(w)[::-1]
        cw = np.cumsum(w) / w.sum()
        ax1.plot(np.arange(1, cw.size + 1), cw, color=COL[t], lw=1.2)
        n50 = int(np.searchsorted(cw, 0.5) + 1)
        ax1.plot(n50, 0.5, "o", color=COL[t], ms=3, mfc="white", mew=0.8)
    ax1.axhline(0.5, color=GREY, lw=0.6, ls=":")
    ax1.set_xscale("log")
    ax1.set_xlim(1, 2e4)
    ax1.set_ylim(0, 1.02)
    ax1.set_xlabel("Genes ranked by $\\mu_{32}^2$")
    ax1.set_ylabel(r"Cumulative share of $\Sigma\mu^2$")
    xt(ax1, [1, 10, 100, 1000, 10000])
    style(ax1)
    panel(ax1, "a")

    # (b) 累积 gamma
    for t in ("colon", "lung", "brain"):
        cu = N["jackknife"][t]["cumulative"]
        k = np.array(cu["k"], float)
        g = np.array(cu["gamma"], float)
        ax2.plot(k, g, color=COL[t], lw=1.2, marker="o", ms=2.0, mfc="white", mew=0.6)
        full = N["jackknife"][t]["gamma_full"]
        ax2.axhline(full, color=COL[t], lw=0.5, ls=(0, (4, 2)), alpha=0.8)
        rk = N["jackknife"][t]["random"]
        ks = [int(x) for x in rk if rk[x]["med"] == rk[x]["med"]]
        lo = [rk[str(x)]["lo"] for x in ks]
        hi = [rk[str(x)]["hi"] for x in ks]
        ax2.fill_between(ks, lo, hi, color=COL[t], alpha=0.10, lw=0)
    ax2.set_xscale("log")
    ax2.set_xlim(0.9, 2.2e4)
    ax2.set_ylim(0, 0.62)
    ax2.set_xlabel("Number of genes used (top-$k$)")
    ax2.set_ylabel(r"$\gamma$ recomputed")
    xt(ax2, [1, 10, 100, 1000, 10000])
    ax2.text(1.1, 0.585, "shaded: random-$k$ 95% range", fontsize=6.3, color=MGREY)
    style(ax2)
    panel(ax2, "b")

    # (c) 留一位移
    for k, t in enumerate(("colon", "lung", "brain")):
        jk = N["jackknife"][t]["loo"]
        d = abs(jk["max_abs_delta"])
        ax3.bar(k, d, width=0.55, color=COL[t], edgecolor="none")
        ax3.text(k, d + 0.003, f"{d:.3f}", ha="center", va="bottom", fontsize=7, color=COL[t])
        ax3.text(k, -0.008, jk["worst_gene"], ha="center", va="top", fontsize=6.5, color=COL[t])
        # 条内不再写字：小字号白字在 1% 高度里会糊成一团，改在顶部注记
    ax3.set_xticks([])
    ax3.set_ylim(0, 0.148)
    ax3.set_ylabel(r"Max $|\Delta\gamma|$ from dropping one gene")
    cnt = [N["jackknife"][t]["loo"]["n_delta_gt_0p05"] for t in ("colon", "lung", "brain")]
    ngn = [N["jackknife"][t]["loo"]["n"] for t in ("colon", "lung", "brain")]
    ax3.text(0.02, 0.141, "colon / lung / brain", fontsize=6.4, color=MGREY)
    ax3.text(0.02, 0.128,
             "genes with $|\\Delta\\gamma|>0.05$:  "
             + " \u00b7 ".join(f"{c:,}/{n:,}" for c, n in zip(cnt, ngn)),
             fontsize=5.4, color=MGREY)
    ax3.spines["bottom"].set_visible(True)
    style(ax3)
    panel(ax3, "c", dx=-0.22)
    save(fig, od, "Fig4")


# ============================================================ Fig 5
def fig5(od):
    fig = plt.figure(figsize=(W2, 118 * MM))
    ax1 = fig.add_axes([0.075, 0.60, 0.36, 0.34])
    ax2 = fig.add_axes([0.565, 0.60, 0.375, 0.34])
    ax3 = fig.add_axes([0.075, 0.09, 0.36, 0.34])
    ax4 = fig.add_axes([0.565, 0.09, 0.375, 0.34])

    bands = N["splithalf"]["bands"]
    order = ["q1", "q2", "q3", "q4", "halfA", "halfB"]
    vals = [bands[k]["gamma"] for k in order]
    ax1.bar(range(6), vals, width=0.62,
            color=[COL["colon"]] * 4 + ["#7F9CB5", "#7F9CB5"], edgecolor="none")
    full = N["jackknife"]["colon"]["gamma_full"]
    ax1.axhline(full, color=MGREY, lw=0.7, ls=(0, (4, 2)))
    ax1.text(5.55, full + 0.004, "whole section", fontsize=6.3, color=MGREY, ha="right")
    ax1.set_xticks(range(6))
    ax1.set_xticklabels(["q1", "q2", "q3", "q4", "A", "B"])
    ax1.set_ylim(0.42, 0.68)
    ax1.set_ylabel(r"$\gamma$ recomputed per region")
    ax1.set_xlabel("Column bands of one colon section")
    ax1.text(0.05, 0.655, f"range {N['splithalf']['quarter_range']:.3f}   "
                          f"SD {N['splithalf']['quarter_sd']:.3f}",
             fontsize=6.8, color=COL["colon"])
    style(ax1)
    panel(ax1, "a")

    # (b) 区域差 vs 组织差
    gaps = [N["derived3"]["tissue_gap_colon_lung"], N["derived3"]["tissue_gap_colon_brain"]]
    names = ["colon–lung", "colon–brain"]
    hw = N["derived3"]["halfwidth95_regional_2sec"]
    hw2 = 1.96 * np.sqrt(2) * N["datasets"]["colon"]["rms_gamma_cut"]
    ax2.barh([0, 1], gaps, height=0.42, color=[COL["lung"], COL["brain"]], edgecolor="none")
    ax2.errorbar(gaps, [0, 1], xerr=[[hw, hw], [hw, hw]], fmt="none", ecolor=MGREY,
                 elinewidth=0.8, capsize=2.5)
    for i, g in enumerate(gaps):
        ax2.text(g + 0.012, i + 0.02, f"{g:.3f}", va="center", fontsize=7,
                 color=[COL["lung"], COL["brain"]][i])
    ax2.axvline(N["splithalf"]["quarter_range"], color=GREY, lw=0.7, ls=(0, (4, 2)))
    ax2.text(N["splithalf"]["quarter_range"] - 0.004, 1.62,
             f"regional range {N['splithalf']['quarter_range']:.3f}",
             fontsize=6.3, color=GREY, ha="right")
    ax2.set_yticks([0, 1])
    ax2.set_yticklabels(names)
    ax2.set_xlim(0, 0.30)
    ax2.set_xlabel(r"Difference in $\gamma$   (bars: 95% half-width from regional SD)")
    style(ax2)
    panel(ax2, "b", dx=-0.20)

    # (c) 深度稀释
    runs = N["depth"]["runs"]
    ps = [runs[k]["p"] for k in ("full", "p0.46", "p0.25", "p0.1")]
    gs = [runs[k]["gamma"] for k in ("full", "p0.46", "p0.25", "p0.1")]
    ax3.plot([1, 0.46, 0.25, 0.10], gs, color=COL["colon"], lw=1.2, marker="o",
             ms=3.2, mfc="white", mew=0.8)
    ax3.set_xscale("log")
    ax3.set_xlim(0.085, 1.15)
    ax3.set_ylim(0.548, 0.562)
    ax3.invert_xaxis()
    xt(ax3, [1.0, 0.46, 0.25, 0.10], ["1.0", "0.46", "0.25", "0.10"])
    ax3.set_xlabel("Fraction of UMIs retained (binomial thinning)")
    ax3.set_ylabel(r"$\gamma$")
    ax3.text(1.05, 0.5595, f"max $|\\Delta\\gamma|$ = {N['derived3']['depth_delta_max']:.4f}",
             fontsize=6.8, color=COL["colon"])
    style(ax3)
    panel(ax3, "c")

    # (d) 块自助
    bt = (N.get("boot") or {}).get("by_tissue", {}).get("colon")
    if bt:
        bs = [int(k) for k in bt["results"]]
        sd = [bt["results"][str(b)]["sd"] for b in bs]
        ax4.plot(bs, sd, color=COL["colon"], lw=1.2, marker="o", ms=3.2, mfc="white", mew=0.8)
        for b, s_ in zip(bs, sd):
            ax4.text(b, s_ + 0.002, f"{s_:.3f}", ha="center", fontsize=6.2, color=COL["colon"])
    else:
        ax4.text(0.5, 0.5, "block bootstrap not run", transform=ax4.transAxes,
                 ha="center", va="center", fontsize=7, color=GREY)
    ax4.axhline(N["splithalf"]["quarter_sd"], color=GREY, lw=0.7, ls=(0, (4, 2)))
    ax4.text(520, N["splithalf"]["quarter_sd"] + 0.003,
             "SD across the four bands", fontsize=6.3, color=GREY)
    ax4.set_xscale("log")
    ax4.set_xlim(420, 2600)
    ax4.set_ylim(0, 0.085)
    xt(ax4, [512, 1024, 2048], ["512", "1024", "2048"])
    ax4.set_xlabel("Superblock side (µm)")
    ax4.set_ylabel(r"SD of $\gamma$ across 200 resamples")
    style(ax4)
    panel(ax4, "d", dx=-0.20)
    save(fig, od, "Fig5")


# ============================================================ Fig 6
def fig6(od):
    fig = plt.figure(figsize=(W2, 60 * MM))
    ax1 = fig.add_axes([0.065, 0.19, 0.245, 0.71])
    ax2 = fig.add_axes([0.405, 0.19, 0.245, 0.71])
    ax3 = fig.add_axes([0.775, 0.19, 0.21, 0.71])

    for t in ("colon", "lung", "brain"):
        obs = [D[t]["cv"][str(d)] for d in DS]
        poi = [D[t]["cv_poisson"][str(d)] for d in DS]
        ax1.plot(DS, obs, color=COL[t], lw=1.2, marker="o", ms=2.4, mfc="white", mew=0.7)
        ax1.plot(DS, poi, color=COL[t], lw=0.8, ls=(0, (2, 2)), alpha=0.85)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlim(1.8, 300)
    ax1.set_ylim(8e-4, 2.2)
    xt(ax1, [2, 8, 32, 128])
    ax1.set_xlabel("Bin size $d$ (µm)")
    ax1.set_ylabel("CV of counts per bin")
    ax1.text(2.2, 1.55, "solid: observed   dashed: Poisson", fontsize=6.3, color=MGREY)
    style(ax1)
    panel(ax1, "a")

    ns = np.array([1, 2, 3, 4, 5, 6], float)
    for t in ("colon", "lung", "brain"):
        rms = D[t]["rms_gamma_cut"]
        ax2.plot(ns, 2.802 * rms / np.sqrt(ns), color=COL[t], lw=1.2, marker="o",
                 ms=2.4, mfc="white", mew=0.7)
    ax2.axhline(N["derived3"]["tissue_gap_colon_lung"], color=COL["lung"], lw=0.6, ls=":")
    ax2.axhline(N["derived3"]["tissue_gap_colon_brain"], color=COL["brain"], lw=0.6, ls=":")
    ax2.text(6.25, N["derived3"]["tissue_gap_colon_brain"] * 1.18, "colon–brain gap",
             fontsize=6.2, color=COL["brain"], ha="right")
    ax2.text(6.25, N["derived3"]["tissue_gap_colon_lung"] * 1.22, "colon–lung gap",
             fontsize=6.2, color=COL["lung"], ha="right")
    ax2.set_yscale("log")
    ax2.set_xlim(0.8, 6.4)
    ax2.set_ylim(0.006, 0.5)
    ax2.set_yticks([0.01, 0.02, 0.05, 0.1, 0.2, 0.5])
    ax2.set_yticklabels(["0.01", "0.02", "0.05", "0.1", "0.2", "0.5"])
    ax2.yaxis.set_minor_locator(NullLocator())
    ax2.set_xlabel("Sections per group")
    ax2.set_ylabel(r"MDE of $\Delta\gamma$ (80% power)")
    style(ax2)
    panel(ax2, "b")

    items = [
        ("depth\nto 10%", N["derived3"]["depth_delta_max"], MGREY),
        ("drop top-1\ngene (colon)", N["jackknife"]["colon"]["loo"]["max_abs_delta"], COL["colon"]),
        ("regional band\n(colon)", N["splithalf"]["quarter_range"], "#7F9CB5"),
        ("fit window\n(lung)", N["derived3"]["lung_window_delta"], COL["lung"]),
        ("expression\nlayer", N["variance"]["pooled"]["range_layer"], "#B9A48C"),
    ]
    items.sort(key=lambda x: x[1])
    y = np.arange(len(items))
    ax3.barh(y, [v for _, v, _ in items], height=0.6, color=[c for _, _, c in items],
             edgecolor="none")
    ax3.set_yticks(y)
    ax3.set_yticklabels([n for n, _, _ in items], fontsize=6.2, linespacing=1.15)
    ax3.set_xscale("log")
    ax3.set_xlim(3e-4, 1.2)
    ax3.set_xlabel(r"$|\Delta\gamma|$ produced by one change")
    for i, (_, v, _) in enumerate(items):
        ax3.text(v * 1.15, i, f"{v:.3g}", va="center", fontsize=6.2, color=MGREY)
    style(ax3)
    panel(ax3, "c", dx=-0.46)
    save(fig, od, "Fig6")


# ============================================================ 补充图
def figS1(od):
    fig = plt.figure(figsize=(W1, 62 * MM))
    ax = fig.add_axes([0.17, 0.19, 0.78, 0.74])
    for t in ("colon", "lung", "brain"):
        w = N["legacy"]["window"][t]
        ax.plot(w["starts"], w["gamma_by_start"], color=COL[t], lw=1.1, marker="o",
                ms=2.6, mfc="white", mew=0.7, label=LAB[t])
    ax.set_xscale("log", base=2)
    ax.set_xticks([2, 4, 8, 16, 32, 64])
    ax.set_xticklabels(["2", "4", "8", "16", "32", "64"])
    ax.set_xlabel("First bin size in the fit (µm; last fixed at 256 µm)")
    ax.set_ylabel(r"$\gamma$")
    ax.legend(loc="lower right", fontsize=6.8)
    style(ax)
    save(fig, od, "FigS1")


def figS2(od):
    fig = plt.figure(figsize=(W1, 62 * MM))
    ax = fig.add_axes([0.17, 0.19, 0.78, 0.74])
    sy = N["synthetic"]["expression"]["by_layer"]
    ks = list(sy)
    x = np.arange(len(ks))
    med = [sy[k]["med"] for k in ks]
    lo = [sy[k]["q25"] for k in ks]
    hi = [sy[k]["q75"] for k in ks]
    ax.errorbar(x, med, yerr=[np.array(med) - lo, np.array(hi) - med], fmt="o",
                color=MGREY, ms=3, mfc="white", mew=0.8, elinewidth=0.8, capsize=2)
    ax.axhline(N["synthetic"]["expression"]["beta_true"], color=COL["colon"], lw=0.7, ls=":")
    ax.text(0.05, N["synthetic"]["expression"]["beta_true"] + 0.012,
            "input exponent 0.55", fontsize=6.6, color=COL["colon"])
    ax.set_xticks(x)
    ax.set_xticklabels([k.replace("-", "–") for k in ks], rotation=35, ha="right", fontsize=6.4)
    ax.set_xlabel("Expression layer (counts per 32 µm bin)")
    ax.set_ylabel(r"Recovered $\gamma_g$ (median, IQR)")
    ax.set_ylim(0.55, 0.72)
    style(ax)
    save(fig, od, "FigS2")


def figS3(od):
    fig = plt.figure(figsize=(W1, 62 * MM))
    ax = fig.add_axes([0.17, 0.19, 0.78, 0.74])
    for k, shade in zip(["q1", "q2", "q3", "q4"], ["#1B3A5C", "#43678A", "#7F9CB5", "#B9CBD9"]):
        a = N["splithalf"]["bands"][k]["alpha"]
        ax.plot(DS, a, color=shade, lw=1.0, marker="o", ms=2.2, mfc="white", mew=0.6)
    for k, ls in (("halfA", (0, (3, 2))), ("halfB", (0, (1, 1.6)))):
        a = N["splithalf"]["bands"][k]["alpha"]
        ax.plot(DS, a, color=GREY, lw=1.0, ls=ls)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1.8, 300)
    xt(ax, [2, 8, 32, 128])
    ax.set_xlabel("Bin size $d$ (µm)")
    ax.set_ylabel(r"Pooled $\alpha(d)$")
    ax.text(2.2, 60, "solid: four column bands\ndashed / dotted: two halves", fontsize=6.3,
            color=MGREY, va="top")
    style(ax)
    save(fig, od, "FigS3")


def figS4(od):
    fig = plt.figure(figsize=(W2, 58 * MM))
    axs = []
    for i, t in enumerate(("colon", "lung", "brain")):
        ax = fig.add_axes([0.065 + i * 0.325, 0.20, 0.25, 0.73])
        axs.append(ax)
        a = [alpha_of(t, d) for d in DS]
        al = [D[t]["alpha_local"][str(d)] for d in DS]
        ax.plot(DS, a, color=COL[t], lw=1.2, marker="o", ms=2.4, mfc="white", mew=0.7)
        ax.plot(DS, al, color=COL[t], lw=1.0, ls=(0, (4, 2)))
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(1.8, 300)
        xt(ax, [2, 8, 32, 128])
        ax.set_title(LAB[t], fontsize=7.5, color=COL[t], pad=2)
        if i == 0:
            ax.set_ylabel(r"Pooled $\alpha(d)$")
        ax.set_xlabel("Bin size $d$ (µm)")
        style(ax)
        ax.text(2.1, a[0] * 0.55,
                f"global {N['jackknife'][t]['gamma_full']:.2f}\nwithin-block "
                f"{abs(D[t]['gamma_local']):.2f}", fontsize=6.2, color=COL[t], va="top")
    save(fig, od, "FigS4")


def figS5(od):
    fig = plt.figure(figsize=(W2, 118 * MM))
    specs = [("colon", "val", "IGKC", 0), ("colon", "total", "all genes", 1),
             ("brain", "val", "Camk2n1", 2), ("brain", "total", "all genes", 3)]
    for tag, field, name, i in specs:
        ax = fig.add_axes([0.055 + (i % 2) * 0.52, 0.545 - (i // 2) * 0.45, 0.36, 0.36])
        z = np.load(ROOT / f"02_OUT/layer_maps/maps_{tag}.npz", allow_pickle=True)
        r, c, v = z["row"].astype(float), z["col"].astype(float), z[field].astype(float)
        vmax = np.percentile(v[v > 0], 99) if (v > 0).any() else 1
        sc = ax.scatter(c, -r, c=np.minimum(v, vmax), s=0.9, cmap="magma", vmin=0,
                        vmax=vmax, edgecolors="none", rasterized=True)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        cb = plt.colorbar(sc, ax=ax, fraction=0.040, pad=0.015)
        cb.ax.tick_params(labelsize=5.8, width=0.5, length=2)
        cb.outline.set_linewidth(0.4)
        cb.set_label("counts per 32 µm bin", fontsize=6)
        ax.set_title(f"{LAB[tag]} · {name}", fontsize=7, color=COL[tag], pad=2)
    save(fig, od, "FigS5")


def figS6(od):
    """(a) effective gene count under the two weight conventions.
       (b) stability of five candidate summary scalars across four bands."""
    dec = N["decomp"]["by_tissue"]
    rob = N["robust"]
    fig = plt.figure(figsize=(W2, 82 * MM))

    # ---- (a)
    ax = fig.add_axes([0.075, 0.20, 0.34, 0.70])
    tis = ["colon", "lung", "brain"]
    x = np.arange(len(tis))
    w = 0.34
    v_mu2 = [dec[t]["n_eff_mu2"] for t in tis]
    v_al = [dec[t]["n_eff_alpha"] for t in tis]
    ax.bar(x - w / 2, v_mu2, w, color=[COL[t] for t in tis], alpha=0.35,
           edgecolor=[COL[t] for t in tis], linewidth=0.6, label="by $\\mu^2$ share")
    ax.bar(x + w / 2, v_al, w, color=[COL[t] for t in tis],
           edgecolor=[COL[t] for t in tis], linewidth=0.6, label="by $\\alpha$ contribution")
    ax.set_yscale("log")
    ax.set_ylim(0.7, 200)
    ax.set_xticks(x)
    ax.set_xticklabels([LAB[t] for t in tis])
    ax.set_ylabel("effective number of genes")
    style(ax)
    for i, t in enumerate(tis):
        ax.text(i + w / 2, v_al[i] * 1.18, f"{v_al[i]:.1f}", ha="center", va="bottom", fontsize=6.5)
        ax.text(i - w / 2, v_mu2[i] * 1.18, f"{v_mu2[i]:.1f}", ha="center", va="bottom", fontsize=6.5)
        ax.text(i, 130, f"+{dec[t]['rel_diff'] * 100:.1f}%", ha="center", va="center", fontsize=6.2,
                color=MGREY)
    ax.text(0.02, 1.05, "reconstruction residual", transform=ax.transAxes, fontsize=5.6,
            color=MGREY, va="bottom", ha="left")
    ax.legend(loc="lower left", fontsize=6, handlelength=1.2, bbox_to_anchor=(-0.02, -0.30), ncol=2)
    panel(ax, "a")

    # ---- (b)
    ax = fig.add_axes([0.60, 0.20, 0.37, 0.70])
    keys = ["pooled", "median", "mean", "trimmed", "wmedian"]
    names = ["pooled $\\gamma$", "per-gene median", "per-gene mean", "trimmed mean",
             "$\\mu^2$-weighted median"]
    bands = [b for b in rob["bands"] if b["name"].startswith("q")]
    ypos = np.arange(len(keys))[::-1]
    for y, k in zip(ypos, keys):
        vals = np.array([b[k] for b in bands])
        ax.scatter(vals, np.full(vals.size, y), s=9, color=MGREY, zorder=3,
                   edgecolors="none", alpha=0.85)
        ax.plot([vals.mean(), vals.mean()], [y - 0.26, y + 0.26], color=COL["colon"], lw=1.1, zorder=4)
        sd = rob["summary"][k]["sd"]
        ax.text(0.775, y, f"{sd:.3f}", fontsize=6.2, va="center",
                ha="right", color=COL["colon"] if k == "mean" else MGREY)
    ax.set_yticks(ypos)
    ax.set_yticklabels(names, fontsize=6.6)
    ax.set_xlim(0.24, 0.80)
    ax.set_ylim(-0.6, len(keys) - 0.4)
    ax.set_xlabel("summary scalar within one band")
    ax.text(0.775, len(keys) - 0.35, "SD", fontsize=6.2, ha="right", va="bottom", color=MGREY)
    style(ax)
    panel(ax, "b", dx=-0.42)

    save(fig, od, "FigS6")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=str(ROOT / "07_SUBMISSION/_build_v2.0/figs"))
    a = ap.parse_args()
    od = Path(a.outdir)
    od.mkdir(parents=True, exist_ok=True)
    fig1(od)
    fig2(od)
    fig3(od)
    fig4(od)
    fig5(od)
    fig6(od)
    figS1(od)
    figS2(od)
    figS3(od)
    figS4(od)
    figS5(od)
    figS6(od)


if __name__ == "__main__":
    main()

r"""99_fig45_v4.py — v4 新增两图

Fig 4  两条独立路径给出同一个指数
  (a) 逐基因 gamma（聚合路径）vs beta（8um 变差函数路径），三套组织
  (b) 中位 rho(h) 的局部斜率 vs lag —— 结肠/肺平坦，脑单调上升（无单一 beta）
  (c) gamma 窗口 × beta lag 窗口的 6x6 相关系数矩阵，框出先验钉死的那格

Fig 5  移动转录本就会移动指数
  (a) 六个基因各自置换前后的 gamma（配对哑铃），95% 区间为置换重复
  (b) 剂量-反应：置换 top-N 后 pooled gamma

风格沿用 93_fig9_theory.py（Nature 双栏、空心点、深蓝/中蓝/橙、无上右脊、pdf.fonttype 42）。

用法
----
  python 99_fig45_v4.py --outdir 02_OUT/figs_v4
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[1]
C_COLON = "#2C5F7C"
C_LUNG = "#5B8FA8"
C_BRAIN = "#C1663D"
C_GREY = "#8C8C8C"
C_DARK = "#333333"
LAB = {"colon": "Colon", "lung": "Lung", "brain": "Brain"}
COL = {"colon": C_COLON, "lung": C_LUNG, "brain": C_BRAIN}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.linewidth": 0.7,
    "axes.labelsize": 7.5,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6.2,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.4,
    "ytick.major.size": 2.4,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "savefig.dpi": 400,
})

SETS = {
    "colon": ("02_OUT/gene_moments.npz", "02_OUT/variogram/beta_g_colon.npz"),
    "lung": ("D3_LungCancer/02_OUT/gene_moments_noMT.npz", "02_OUT/variogram/beta_g_lung.npz"),
    "brain": ("D2_MouseBrain/02_OUT/gene_moments_noMT.npz", "02_OUT/variogram/beta_g_brain.npz"),
}


def load_gb(tag):
    mo_p, bt_p = SETS[tag]
    mo = np.load(ROOT / mo_p, allow_pickle=True)
    ds = np.array([int(x) for x in mo["ds"]], float)
    genes_mo = np.array([str(g) for g in mo["genes"]])
    MU = np.column_stack([mo[f"mu_{int(d)}"] for d in ds])
    VA = np.column_stack([mo[f"var_{int(d)}"] for d in ds])
    with np.errstate(invalid="ignore", divide="ignore"):
        AL = (VA - MU) / np.where(MU > 0, MU ** 2, np.nan)
    bt = np.load(ROOT / bt_p, allow_pickle=True)
    bg = np.array([str(g) for g in bt["gene"]])
    beta = bt["beta"].astype(float)
    rho = bt["rho"].astype(float)
    lags = bt["lags"].astype(float)
    valid = bt["valid"].astype(bool)
    pos = {g: i for i, g in enumerate(genes_mo)}
    A = AL[np.array([pos[g] for g in bg])]
    ok = np.all(np.isfinite(A) & (A > 0), axis=1)
    x = np.log(ds); xm = x - x.mean()
    y = np.log(np.where(ok[:, None], A, 1.0))
    ym = y - y.mean(axis=1, keepdims=True)
    g = np.full(bg.size, np.nan)
    g[ok] = -np.sum(xm[None, :] * ym, axis=1)[ok] / np.sum(xm ** 2)
    m = valid & np.isfinite(g) & (beta > 0) & (g > 0)
    return g[m], beta[m], rho[valid], lags


def deming(x, y):
    sx2, sy2 = np.var(x, ddof=1), np.var(y, ddof=1)
    sxy = np.cov(x, y, ddof=1)[0, 1]
    d = sy2 - sx2
    b = (d + np.sqrt(d ** 2 + 4 * sxy ** 2)) / (2 * sxy)
    return float(b), float(np.mean(y) - b * np.mean(x))


def fig4(outdir):
    fig = plt.figure(figsize=(7.2, 2.35))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.0, 1.10], wspace=0.62,
                          left=0.072, right=0.985, top=0.90, bottom=0.17)

    # (a) 散点
    ax = fig.add_subplot(gs[0, 0])
    stats = {}
    for tag in ("colon", "lung", "brain"):
        g, b, _, _ = load_gb(tag)
        ax.scatter(b, g, s=1.6, c=COL[tag], alpha=0.30, linewidths=0, rasterized=True)
        bb, aa = deming(np.log(b), np.log(g))
        xs = np.array([b.min(), b.max()])
        ax.plot(xs, np.exp(aa) * xs ** bb, color=COL[tag], lw=0.9)
        stats[tag] = (np.corrcoef(g, b)[0, 1], len(g), np.median(g - b))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(0.15, 2.6); ax.set_ylim(0.16, 1.5)
    ax.plot([0.16, 1.45], [0.16, 1.45], ls=(0, (3, 2)), color=C_GREY, lw=0.7)
    tk = [0.2, 0.3, 0.5, 1.0, 2.0]
    ax.set_xticks(tk); ax.set_xticklabels([f"{v:g}" for v in tk])
    ax.set_yticks(tk); ax.set_yticklabels([f"{v:g}" for v in tk])
    ax.minorticks_off()
    ax.set_xlabel(r"covariogram exponent  $\beta$")
    ax.set_ylabel(r"aggregation exponent  $\gamma$")
    ax.text(0.045, 0.955, "a", transform=ax.transAxes, fontsize=9, va="top")
    txt = "\n".join(
        f"{LAB[t]}  r={stats[t][0]:.2f}  Δ={stats[t][2]:+.3f}" for t in ("colon", "lung", "brain"))
    ax.text(0.975, 0.045, txt, transform=ax.transAxes, ha="right", va="bottom", fontsize=5.7,
            linespacing=1.45, bbox=dict(fc="white", ec="none", alpha=0.82, pad=1.2))

    # (b) 局部斜率
    ax = fig.add_subplot(gs[0, 1])
    for tag in ("colon", "lung", "brain"):
        _, _, rho, lags = load_gb(tag)
        med = np.nanmedian(rho, axis=0)
        ls, hh = [], []
        for i in range(1, len(lags) - 1):
            if med[i] > 0 and med[i + 1] > 0:
                ls.append(-(np.log(med[i + 1]) - np.log(med[i])) / np.log(lags[i + 1] / lags[i]))
                hh.append(lags[i] * 8)
        ax.plot(hh, ls, "o-", color=COL[tag], ms=2.0, lw=0.9, mfc="white", mew=0.7,
                label=LAB[tag])
    ax.set_xscale("log", base=2)
    ax.set_xlabel("lag (µm)")
    ax.set_ylabel(r"$-\Delta\log\rho\,/\,\Delta\log h$")
    ax.set_ylim(0, 1.75)
    ax.set_xticks([8, 16, 32, 64, 128, 256])
    ax.set_xticklabels(["8", "16", "32", "64", "128", "256"])
    ax.minorticks_off()
    ax.legend(frameon=False, loc="upper left", handlelength=1.4)
    ax.text(0.045, 0.955, "b", transform=ax.transAxes, fontsize=9, va="top")

    # (c) 窗口矩阵
    ax = fig.add_subplot(gs[0, 2])
    W = json.loads((ROOT / "02_OUT/variogram/windowscan_colon.json").read_text(encoding="utf-8"))
    gw = sorted((tuple(w) for w in W["gamma_windows"]))   # 按 (下界, 上界) 升序，禁用文件里的乱序
    bw = sorted((tuple(w) for w in W["beta_windows"]))
    M = np.full((len(gw), len(bw)), np.nan)
    for r in W["grid"]:
        i = gw.index(tuple(r["gamma_window"]))
        j = bw.index(tuple(r["beta_window"]))
        M[i, j] = r["pearson"]
    im = ax.imshow(M, cmap="Blues", vmin=0, vmax=1, aspect="auto", origin="lower")
    for i in range(len(gw)):
        for j in range(len(bw)):
            v = M[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=5.0,
                    color="white" if v > 0.62 else C_DARK)
    ax.set_xticks(range(len(bw)))
    ax.set_xticklabels([f"{a}-{b}" for a, b in bw], rotation=45, ha="right", fontsize=5.4)
    ax.set_yticks(range(len(gw)))
    ax.set_yticklabels([f"{a}-{b}" for a, b in gw], fontsize=5.4)
    ax.set_xlabel(r"$\beta$ lag window (grid units)")
    ax.set_ylabel(r"$\gamma$ fit window (µm)")
    # 先验窗口 = 正文钉死的 γ 2–256 × β 2–8（网格单位），框必须落在这格
    ax.add_patch(Rectangle((bw.index((2, 8)) - 0.5, gw.index((2, 256)) - 0.5), 1, 1,
                           fill=False, ec="#B23A2E", lw=1.2))
    ax.text(-0.30, 1.055, "c", transform=ax.transAxes, fontsize=9, va="bottom", color=C_DARK)
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_linewidth(0.5)

    outdir.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(outdir / f"Fig4_gamma_beta.{ext}", bbox_inches="tight",
                    metadata={"CreationDate": None} if ext == "pdf" else None)
    plt.close(fig)
    print(f"写出 {outdir}/Fig4_gamma_beta.pdf|png")


def fig5(outdir):
    ex = json.loads((ROOT / "02_OUT/excise/excise_colon.json").read_text(encoding="utf-8"))
    mp = json.loads((ROOT / "02_OUT/excise/multigene_permute.json").read_text(encoding="utf-8"))

    fig = plt.figure(figsize=(7.2, 2.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.05], wspace=0.34,
                          left=0.115, right=0.975, top=0.90, bottom=0.17)

    # (a) 六基因哑铃
    ax = fig.add_subplot(gs[0, 0])
    genes = list(ex["genes"].keys())
    order = sorted(genes, key=lambda g: ex["genes"][g]["baseline"]["gamma"])
    for i, gname in enumerate(order):
        v = ex["genes"][gname]
        b = v["baseline"]["gamma"]
        p = v["permute"]["gamma_median"]
        lo, hi = v["permute"]["gamma_lo"], v["permute"]["gamma_hi"]
        ax.plot([b, p], [i, i], color=C_GREY, lw=0.7, zorder=1)
        ax.plot([lo, hi], [i, i], color=C_COLON, lw=2.6, alpha=0.30,
                solid_capstyle="butt", zorder=2)
        ax.plot(b, i, "o", ms=3.4, mfc="white", mec=C_DARK, mew=0.7, zorder=3)
        ax.plot(p, i, "o", ms=3.4, color=C_COLON, zorder=3)
        ax.text(1.72, i, gname, fontsize=6.0, va="center", ha="left")
    ax.set_yticks([])
    ax.set_xlim(0.2, 1.78)
    ax.set_ylim(-0.7, len(order) - 0.3)
    ax.set_xlabel(r"per-gene exponent  $\gamma_g$")
    ax.text(0.03, 0.965, "a", transform=ax.transAxes, fontsize=9, va="top")
    ax.plot([], [], "o", ms=3.4, mfc="white", mec=C_DARK, mew=0.7, label="observed")
    ax.plot([], [], "o", ms=3.4, color=C_COLON, label="positions permuted")
    ax.legend(frameon=False, loc="lower right", handletextpad=0.4, borderaxespad=0.2)

    # (b) 剂量-反应
    ax = fig.add_subplot(gs[0, 1])
    ns = [0, 10, 30, 100, 300, 1000]
    gam = [mp["curves"][str(n)]["gamma"] for n in ns]
    shr = [mp["curves"][str(n)]["mu2_share"] for n in ns]
    xx = np.array([max(n, 3) for n in ns], float)
    ax.plot(xx, gam, "o-", color=C_COLON, ms=3.4, lw=1.0, mfc="white", mew=0.8)
    for x, y, s in zip(xx, gam, shr):
        dy = 8 if s == 0 else -9
        va = "bottom" if s == 0 else "top"
        ax.annotate(f"{s:.0%}" if s > 0 else "0%", (x, y), textcoords="offset points",
                    xytext=(0, dy), ha="center", va=va, fontsize=5.4, color=C_GREY)
    ax.axhline(2.0, ls=(0, (3, 2)), color=C_GREY, lw=0.7)
    ax.axhline(gam[0], ls=(0, (1, 2)), color=C_GREY, lw=0.7)
    ax.text(3.2, 1.90, "uncorrelated limit  γ = 2", fontsize=5.8, color=C_GREY, va="top")
    ax.set_xscale("log")
    ax.set_xlim(2.4, 1600)
    ax.set_ylim(0.45, 2.12)
    ax.set_xlabel(r"genes permuted, top $N$ by $\mu^2$ weight")
    ax.set_ylabel(r"pooled exponent  $\gamma$")
    ax.text(0.03, 0.965, "b", transform=ax.transAxes, fontsize=9, va="top")
    ax.set_xticks([3, 10, 30, 100, 300, 1000])
    ax.set_xticklabels(["0", "10", "30", "100", "300", "1000"])

    outdir.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(outdir / f"Fig5_intervention.{ext}", bbox_inches="tight",
                    metadata={"CreationDate": None} if ext == "pdf" else None)
    plt.close(fig)
    print(f"写出 {outdir}/Fig5_intervention.pdf|png")



def fig3(outdir):
    """理论图：四个面板逐条对应正文的四条论断。"""
    th = json.loads((ROOT / "02_OUT/theory/theory_check.json").read_text(encoding="utf-8"))
    E = th["experiments"]
    fs = plt.rcParams["font.size"]
    fig = plt.figure(figsize=(7.2, 4.5))
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.32,
                          left=0.085, right=0.975, top=0.95, bottom=0.095)

    # (a) 恒等式：预测 vs 实测
    ax = fig.add_subplot(gs[0, 0])
    def collect(src):
        a, q = [], []
        for v in src.values():
            for row in v["rows"]:
                a.append(row["alpha"]); q.append(row["predicted"])
        return np.array(a), np.array(q)

    a1, p1 = collect(E["E1_identity"])
    a2, p2 = collect(E["E5_cv_sweep"])
    ax_ = np.concatenate([a1, a2]); pr_ = np.concatenate([p1, p2])
    ax.plot(p1, a1, "o", ms=1.8, color=C_COLON, alpha=0.55, mew=0, rasterized=True)
    ax.plot(p2, a2, "s", ms=1.8, color=C_BRAIN, alpha=0.55, mew=0, rasterized=True)
    lim = np.array([ax_.min() * 0.7, ax_.max() * 1.4])
    ax.plot(lim, lim, ls=(0, (3, 2)), color=C_GREY, lw=0.7)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lim); ax.set_ylim(lim)
    rel = np.abs(ax_ - pr_) / pr_
    ax.set_xlabel(r"predicted  $CV^2_0\,\bar{\rho}(A)$")
    ax.set_ylabel(r"measured  $\alpha(A)$")
    ax.text(0.045, 0.955, "a", transform=ax.transAxes, fontsize=9, va="top")
    ax.text(0.96, 0.06, f"max rel. err {rel.max():.1%}\nmedian {np.median(rel):.1%}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=5.9, linespacing=1.4)

    # (b) gamma vs beta_emp（幂律族）
    ax = fig.add_subplot(gs[0, 1])
    for tag, fam in (("E1", "beta_"), ("E2", "powerlaw_")):
        src = E["E1_identity"] if tag == "E1" else E["E2_families"]
        gs_, bs_ = [], []
        for k, v in src.items():
            if "powerlaw_" in k or tag == "E1":
                gs_.append(-(v["gamma_fit"] if "gamma_fit" in v else v["gamma_2_64"]))
                bs_.append(v["beta_emp"])
        mk = "o" if tag == "E1" else "s"
        ax.plot(bs_, gs_, mk, ms=3.2, mfc="white", mec=C_COLON, mew=0.8)
        r2 = np.corrcoef(bs_, gs_)[0, 1] ** 2
    allg = [-(v["gamma_fit"]) for v in E["E1_identity"].values()] + \
           [-(E["E2_families"]["powerlaw_0.55"]["gamma_2_64"]),
            -(E["E2_families"]["powerlaw_1.2"]["gamma_2_64"])]
    allb = [v["beta_emp"] for v in E["E1_identity"].values()] + \
           [E["E2_families"]["powerlaw_0.55"]["beta_emp"],
            E["E2_families"]["powerlaw_1.2"]["beta_emp"]]
    r2 = np.corrcoef(allb, allg)[0, 1] ** 2
    lim = np.array([0.25, 1.75])
    ax.plot(lim, lim, ls=(0, (3, 2)), color=C_GREY, lw=0.7)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel(r"measured covariogram exponent  $\beta_{emp}$")
    ax.set_ylabel(r"recovered  $\gamma$")
    ax.text(0.045, 0.955, "b", transform=ax.transAxes, fontsize=9, va="top")
    ax.text(0.96, 0.06, f"$r^2$ = {r2:.4f}", transform=ax.transAxes, ha="right",
            va="bottom", fontsize=6.4)

    # (c) 指数相关：拟合窗口一换，指数就变
    ax = fig.add_subplot(gs[1, 0])
    xis = ["exponential_2", "exponential_8", "exponential_32", "exponential_128"]
    lab = ["2", "8", "32", "128"]
    a64 = [-E["E2_families"][k]["gamma_2_64"] for k in xis]
    a16 = [-E["E2_families"][k]["gamma_16_64"] for k in xis]
    yy = np.arange(len(xis))
    for i in range(len(xis)):
        ax.annotate("", xy=(a16[i], yy[i]), xytext=(a64[i], yy[i]),
                    arrowprops=dict(arrowstyle="-|>", color=C_BRAIN, lw=0.8,
                                    shrinkA=0, shrinkB=0))
    ax.plot(a64, yy, "o", ms=3.4, mfc="white", mec=C_BRAIN, mew=0.8)
    ax.plot(a16, yy, "o", ms=3.4, color=C_BRAIN)
    ax.axvline(2.0, ls=(0, (3, 2)), color=C_GREY, lw=0.7)
    ax.text(1.96, len(xis) - 0.35, "uncorrelated\n$\gamma=2$", fontsize=5.8,
            color=C_GREY, ha="right", va="top", linespacing=1.3)
    ax.plot([-E["E2_families"]["white_0"]["gamma_2_64"]], [len(xis)],
            "D", ms=3.4, color=C_DARK)
    ax.set_yticks(list(yy) + [len(xis)])
    ax.set_yticklabels([f"$\\xi$ = {l}" for l in lab] + ["white"])
    ax.set_ylim(-0.6, len(xis) + 0.6)
    ax.set_xlim(0, 2.15)
    ax.set_xlabel(r"recovered  $\gamma$")
    ax.text(0.045, 0.955, "c", transform=ax.transAxes, fontsize=9, va="top")
    ax.text(0.03, 0.10, "○ fit 2–64 µm    ● fit 16–64 µm", transform=ax.transAxes,
            fontsize=5.8, color=C_BRAIN)

    # (d) 深度不变
    ax = fig.add_subplot(gs[1, 1])
    rows = E["E3_depth"]["rows"]
    mc = np.array([rows[k]["mean_count_at_2um"] for k in rows])
    gg = np.array([-rows[k]["gamma"] for k in rows])
    ax.plot(mc, gg, "o", ms=3.2, color=C_COLON, mfc="white", mew=0.8)
    ax.set_xscale("log")
    ax.set_ylim(gg.min() - 0.02, gg.max() + 0.02)
    ax.set_xlabel("mean counts per 2 µm cell")
    ax.set_ylabel(r"recovered  $\gamma$")
    ax.text(0.045, 0.955, "d", transform=ax.transAxes, fontsize=9, va="top")
    ax.text(0.04, 0.90, f"range {gg.max() - gg.min():.4f}\nover {mc.max() / mc.min():.0f}× depth",
            transform=ax.transAxes, ha="left", va="top", fontsize=5.9, linespacing=1.4)

    outdir.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(outdir / f"Fig3_theory.{ext}", bbox_inches="tight",
                    metadata={"CreationDate": None} if ext == "pdf" else None)
    plt.close(fig)
    print(f"写出 {outdir}/Fig3_theory.pdf|png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="02_OUT/figs_v4")
    a = ap.parse_args()
    od = ROOT / a.outdir
    fig3(od)
    fig4(od)
    fig5(od)


if __name__ == "__main__":
    main()

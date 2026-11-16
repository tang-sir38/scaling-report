r"""93_fig9_theory.py — T2c：理论图（新 Fig 2，投稿版升为 Fig 2）

三个面板讲一件事：alpha 曲线不是经验拟合，是相关函数。
  (a) 三种相关结构下 α(d) 的形状（合成场，真值已知）
      —— 白噪声给出精确的 d^-2；指数给出"先平台后 d^-2"；幂律给出平行直线
  (b) 归一化 α 曲线 = bin 内平均相关函数 rho_bar(d)
      —— 真实三套数据 + 合成幂律场叠加，说明"形状"承载的信息
  (c) beta_g（空间协方差独立测得）vs gamma_g（聚合路径）
      —— 两条独立路径的一致性，即 P2

风格沿用现稿 Fig1：Nature 双栏、空心圆点、深蓝/中蓝/橙、灰色虚线参考、无网格。

用法
----
  python 93_fig9_theory.py --theory 02_OUT/theory/theory_check.json \
      --betag 02_OUT/variogram/beta_g_colon.npz --outdir 02_OUT/figs_v4
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

C_COLON = "#2C5F7C"
C_LUNG = "#5B8FA8"
C_BRAIN = "#C1663D"
C_GREY = "#8C8C8C"
C_DARK = "#333333"

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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--theory", default="02_OUT/theory/theory_check.json")
    ap.add_argument("--betag", default="02_OUT/variogram/beta_g_colon.npz")
    ap.add_argument("--outdir", default="02_OUT/figs_v4")
    ap.add_argument("--ds", default="2,4,8,16,32,64,128,256")
    a = ap.parse_args()

    od = Path(a.outdir); od.mkdir(parents=True, exist_ok=True)
    th = json.load(open(a.theory, encoding="utf-8"))
    ds_real = [int(x) for x in a.ds.split(",")]

    fig = plt.figure(figsize=(7.2, 2.3))
    gs = fig.add_gridspec(1, 3, wspace=0.44, left=0.075, right=0.985,
                          bottom=0.205, top=0.955)

    # ---------------- (a) 三种相关结构的 alpha 形状 ----------------
    ax = fig.add_subplot(gs[0, 0])
    e2 = th["experiments"]["E2_families"]
    style = {
        "white_0": (C_GREY, "--", "white noise"),
        "exponential_8": (C_BRAIN, "-", "exponential $\\xi$=8"),
        "exponential_32": ("#D9A05B", "-", "exponential $\\xi$=32"),
        "powerlaw_0.55": (C_COLON, "-", "power law $\\beta$=0.55"),
    }
    for key, (col, ls, lab) in style.items():
        if key not in e2:
            continue
        ks = np.array(e2[key]["ks"], float)
        al = np.array(e2[key]["alphas"], float)
        m = np.isfinite(al) & (al > 0)
        if key == "white_0":
            m &= ks <= 32          # d^-2 掉得太快，截到 32 免得压扁其他曲线
        ref = float(np.interp(1.0, ks, al))
        ax.plot(ks[m] * 2, al[m] / ref, ls, color=col, lw=0.9,
                label=lab, zorder=2)
    xx = np.array([2, 64])
    ax.plot(xx, (xx / 2) ** -2, ":", color=C_DARK, lw=0.65, zorder=1)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Bin size $d$ (µm)")
    ax.set_ylabel("$\\alpha(d)\\,/\\,\\alpha(2\\,\\mathrm{µm})$")
    ax.set_xlim(1.8, 150); ax.set_ylim(8e-4, 1.8)
    ax.set_xticks([2, 8, 32, 128]); ax.set_xticklabels(["2", "8", "32", "128"])
    ax.legend(frameon=False, loc="lower left", handlelength=1.5,
              borderaxespad=0.15, labelspacing=0.22, fontsize=5.6)
    ax.annotate("$d^{-2}$", xy=(22, (22 / 2) ** -2), xytext=(34, 3.2e-3),
            fontsize=6, color=C_DARK, ha="left", va="center")
    ax.text(-0.22, 1.02, "a", transform=ax.transAxes, fontsize=9, fontweight="bold")

    # ---------------- (b) 真实数据归一化曲线 = rho_bar ----------------
    ax = fig.add_subplot(gs[0, 1])
    lock = json.load(open("02_OUT/NUMBERS_LOCK_v6.json", encoding="utf-8"))
    for tag, col, lab in [("colon", C_COLON, "colon"),
                          ("lung", C_LUNG, "lung"),
                          ("brain", C_BRAIN, "brain")]:
        al = lock["datasets"][tag]["alpha"]
        xs = np.array([d for d in ds_real if str(d) in al], float)
        ys = np.array([al[str(d)] for d in ds_real if str(d) in al], float)
        if not ys.size:
            continue
        ax.plot(xs, ys / ys[0], "o-", color=col, ms=2.5, lw=0.9,
                markerfacecolor="white", markeredgewidth=0.75, label=lab, zorder=3)
    e4 = th["experiments"]["E4_collapse"]
    for key, col, ls in [("beta_0.3", "#B8B8B8", (0, (4, 2))),
                         ("beta_0.55", "#787878", (0, (2, 1.6)))]:
        if key not in e4:
            continue
        ks = np.array(e4[key]["ks"], float) * 2
        yv = np.array(e4[key]["alpha_norm"], float)
        ax.plot(ks, yv, linestyle=ls, color=col, lw=0.75, zorder=1,
                label="synth " + key.split("_")[1])
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Bin size $d$ (µm)")
    ax.set_ylabel("$\\alpha(d)\\,/\\,\\alpha(2\\,\\mathrm{µm})$")
    ax.set_xlim(1.8, 320); ax.set_ylim(0.05, 1.35)
    ax.set_xticks([2, 8, 32, 128]); ax.set_xticklabels(["2", "8", "32", "128"])
    ax.legend(frameon=False, loc="upper right", handlelength=1.7,
              borderaxespad=0.15, labelspacing=0.2, fontsize=5.2, ncol=1)
    ax.text(-0.22, 1.02, "b", transform=ax.transAxes, fontsize=9, fontweight="bold")

    # ---------------- (c) beta_g vs gamma_g ----------------
    ax = fig.add_subplot(gs[0, 2])
    z = np.load(a.betag, allow_pickle=True)
    beta = z["beta"]; gamma = z["gamma"]
    r2 = z["beta_r2"]; ok = z["valid"] & np.isfinite(beta) & np.isfinite(gamma) & (r2 > 0.5)
    ax.scatter(gamma[ok], beta[ok], s=2.0, c=C_COLON, alpha=0.32,
               linewidths=0, rasterized=True)
    lim = [0.05, 1.0]
    ax.plot(lim, lim, "-", color=C_GREY, lw=0.8, zorder=1)
    if ok.sum() > 10:
        p = np.polyfit(gamma[ok], beta[ok], 1)
        xx = np.array(lim)
        ax.plot(xx, np.polyval(p, xx), "--", color=C_BRAIN, lw=0.9, zorder=2)
        r = float(np.corrcoef(gamma[ok], beta[ok])[0, 1])
        ax.text(0.05, 0.95, f"$r$ = {r:.2f}\n$n$ = {int(ok.sum()):,}",
                transform=ax.transAxes, fontsize=6, va="top")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("$\\gamma_g$  (aggregation path)")
    ax.set_ylabel("$\\beta_g$  (spatial covariance)")
    ax.text(-0.22, 1.02, "c", transform=ax.transAxes, fontsize=9, fontweight="bold")
    ax.annotate("1:1", xy=(0.72, 0.72), xytext=(0.86, 0.60), fontsize=5.8,
            color=C_GREY, ha="left", va="center",
            arrowprops=dict(arrowstyle="-", color=C_GREY, lw=0.5,
                            shrinkA=0, shrinkB=1))

    for ext in ("pdf", "png"):
        fig.savefig(od / f"Fig9_theory.{ext}", bbox_inches="tight",
                    metadata={"CreationDate": None} if ext == "pdf" else None)
    print(f"写出 {od}/Fig9_theory.pdf|png")


if __name__ == "__main__":
    main()

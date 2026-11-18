r"""123_fig1_v551.py — M1 的 Fig 1 重画（v5.5.1）。

两处改动，只改这两处：

1. panel b：三条曲线 → 七条曲线，并在 256 µm 端点标出每片的下降倍数。
   目的：标题的 "2- to 15-fold" 第一次在图上落地（v5.5 里它只在 Fig. S3 有曲线、
   任何一张图都没有数字）。

2. panel a：几何自洽。v5.5 的图画的是 32×32 个 1 µm 小格，4×4 小格合成一个
   4 µm 块（8×8 块），方框 16×16 小格；但三个标尺分别跨越 2 格（标 2 µm）、
   4 格（标 4 µm）、16 格（标 16 µm）——同一格同时是 0.5、1、1 µm，读不出来。
   新版：8×8 格，每格 2 µm（原生 bin），整块 16 µm，方框就是整块本身；
   只留两个标尺（2 µm 单格、16 µm 整块），256 µm 的换算写成文字。

数据全部来自 02_OUT/extra_tissues_20261116.json 的 alpha_by_d，现场算 fold，
不手抄任何常数。样式沿用 43_figures_v2.py（Arial 7.5 pt、180 mm 双栏、
PDF 冻结时间戳、PNG 300 dpi、同一配色）。

用法：python 123_fig1_v551.py --outdir <04_图件 或 _build 目录>
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.ticker import NullLocator

MM = 1 / 25.4
W2 = 180 * MM

REPRO = Path(__file__).resolve().parents[1]
SUMMARY = REPRO / "02_OUT" / "extra_tissues_20261116.json"
DS = [2, 4, 8, 16, 32, 64, 128, 256]

STYLE = {
    "colon": ("#1B3A5C", "Colon"),
    "lung": ("#4A6E8F", "Lung"),
    "brain": ("#B5651D", "Brain"),
    "mouse_small_intestine": ("#6E8B3D", "Small intestine"),
    "mouse_kidney": ("#8C4A7E", "Kidney"),
    "human_pancreas": ("#A8452F", "Pancreas"),
    "mouse_embryo": ("#3F7F6F", "Embryo"),
}
# 图例顺序按 fold 从大到小，读者从上往下就是标题区间的两端
ORDER = ["mouse_small_intestine", "colon", "lung", "brain",
         "mouse_kidney", "mouse_embryo", "human_pancreas"]

GREY, LGREY, MGREY = "#8A8A8A", "#D6D6D6", "#4D4D4D"

plt.rcParams.update({
    "font.family": "Arial", "font.size": 7.5, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "xtick.direction": "out", "ytick.direction": "out",
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    "legend.frameon": False, "legend.handlelength": 1.6,
    "legend.borderpad": 0.2, "legend.labelspacing": 0.22,
})


def style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)
    ax.tick_params(which="both", top=False, right=False)


def panel_a(axs):
    """8×8 原生格 = 16 µm；一块 = 8×8 格；只两个标尺。"""
    n = 8                      # 每边 8 个原生 bin
    for i in range(n):
        for j in range(n):
            axs.add_patch(Rectangle((i, j), 1, 1, facecolor="#EDF2F7",
                                    edgecolor=LGREY, lw=0.5))
    # 一个原生 bin 加深，作为 2 µm 标尺的载体
    axs.add_patch(Rectangle((0, 0), 1, 1, facecolor="#1B3A5C", edgecolor="none"))
    # 整块轮廓：这就是 16 µm 的方框，与网格严丝合缝，不再有内嵌的第二个框
    axs.add_patch(Rectangle((0, 0), n, n, facecolor="none",
                            edgecolor="#1B3A5C", lw=1.3, clip_on=False))

    # 标尺 1：单格 = 2 µm。一格只有 3.9 mm 宽，画成带端帽的短杠而不是箭头
    for x0 in (0.0, 1.0):
        axs.plot([x0, x0], [-0.70, -0.42], color=MGREY, lw=0.5, clip_on=False)
    axs.plot([0.0, 1.0], [-0.56, -0.56], color=MGREY, lw=0.5, clip_on=False)
    axs.text(0.5, -0.78, "2 µm", ha="center", va="top", fontsize=6.2, color=MGREY)

    # 标尺 2：整块 = 16 µm
    for x0 in (0.0, float(n)):
        axs.plot([x0, x0], [-2.30, -2.02], color="#1B3A5C", lw=0.6, clip_on=False)
    axs.plot([0.0, float(n)], [-2.16, -2.16], color="#1B3A5C", lw=0.6, clip_on=False)
    axs.text(n / 2, -2.38, "16 µm  (8 native bins)", ha="center", va="top",
             fontsize=6.2, color="#1B3A5C")

    axs.set_xlim(0, n)
    axs.set_ylim(0, n)
    axs.set_aspect("equal")
    axs.axis("off")
    axs.text(0.02, 1.02, "a", transform=axs.transAxes, fontsize=9, fontweight="bold",
             va="bottom", ha="left")
    axs.text(0.5, -0.44, "16,384 native bins\nmake one 256 µm bin",
             transform=axs.transAxes, ha="center", va="top",
             fontsize=6.0, color=MGREY, linespacing=1.3)



def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--summary", default=str(SUMMARY))
    a = ap.parse_args()

    summ = json.loads(Path(a.summary).read_text(encoding="utf-8"))
    sec = {s["key"]: s for s in summ["sections"]}
    C = {k: np.array([float(sec[k]["alpha_by_d"][str(d)]) for d in DS]) for k in ORDER}
    fold = {k: float(C[k][0] / C[k][-1]) for k in ORDER}
    for k in ORDER:                      # 与锁里的 drop_fold 对表，防止两处口径漂移
        assert abs(fold[k] - float(sec[k]["drop_fold"])) < 1e-9, (k, fold[k])

    fig = plt.figure(figsize=(W2, 74 * MM))
    axs = fig.add_axes([0.075, 0.20, 0.175, 0.60])
    ax = fig.add_axes([0.375, 0.155, 0.50, 0.745])

    panel_a(axs)

    # (b) 七条归一化曲线
    dd = np.logspace(np.log10(1.9), np.log10(270), 60)
    ax.plot(dd, (dd / 2.0) ** -2, color=GREY, lw=0.9, ls=(0, (4, 2)), zorder=1)
    ax.plot(dd, np.ones_like(dd), color=GREY, lw=0.8, ls=(0, (1, 2)), zorder=1)
    for k in ORDER:
        col, lab = STYLE[k]
        ax.plot(DS, C[k] / C[k][0], color=col, lw=1.0, marker="o", ms=2.3,
                mfc="white", mew=0.7, zorder=3, label=lab)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1.9, 420)
    ax.set_ylim(3.5e-5, 3.4)
    ax.set_xticks([2, 8, 32, 128, 256])
    ax.set_xticklabels(["2", "8", "32", "128", "256"])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_yticks([1e-4, 1e-3, 1e-2, 1e-1, 1])
    ax.set_yticklabels(["$10^{-4}$", "$10^{-3}$", "$10^{-2}$", "$10^{-1}$", "1"])
    ax.yaxis.set_minor_locator(NullLocator())
    ax.set_xlabel("Bin size $d$ (µm)")
    ax.set_ylabel(r"$\alpha(d)\,/\,\alpha(2\,\mathrm{µm})$")
    style(ax)
    ax.legend(loc="lower left", fontsize=6.2, ncol=4, columnspacing=0.8,
              handletextpad=0.4, bbox_to_anchor=(0.0, 1.0))
    ax.annotate("counting statistics\n" + r"$\alpha \propto d^{-2}$",
                xy=(120, (120 / 2) ** -2), xytext=(26, 2.6e-3), fontsize=6.6,
                color=MGREY, va="top",
                arrowprops=dict(arrowstyle="-", color=GREY, lw=0.5, shrinkA=2, shrinkB=2))
    ax.text(60, 1.35, r"fully correlated  $\alpha \propto d^{0}$", fontsize=6.6, color=MGREY)

    # 端点标 fold：先按显示坐标算实际 y，再逐个下推，避免三个陡曲线互相压字
    fig.canvas.draw()
    inv = ax.transData.inverted()
    x_lab = 330.0
    y_end = {k: C[k][-1] / C[k][0] for k in ORDER}
    order_desc = sorted(ORDER, key=lambda k: -y_end[k])
    px_per_pt = fig.dpi / 72.0
    min_px = 9.2 * px_per_pt                       # 一行 6.4 pt 文字的最小间距
    placed: list[float] = []
    for k in order_desc:
        y_disp = ax.transData.transform((x_lab, y_end[k]))[1]
        if placed and placed[-1] - y_disp < min_px:
            y_disp = placed[-1] - min_px
        placed.append(y_disp)
        y_txt = inv.transform((x_lab, y_disp))[1]
        col = STYLE[k][0]
        ax.plot([256, x_lab - 8], [y_end[k], y_txt], color=col, lw=0.4,
                solid_capstyle="butt", zorder=2, clip_on=False, alpha=0.75)
        ax.text(x_lab, y_txt, f"{fold[k]:.1f}×", color=col, fontsize=6.4,
                ha="left", va="center", clip_on=False)

    panel_letter = ax.text(-0.155, 1.03, "b", transform=ax.transAxes, fontsize=9,
                           fontweight="bold", va="bottom", ha="left")
    _ = panel_letter

    od = Path(a.outdir)
    od.mkdir(parents=True, exist_ok=True)
    fig.savefig(od / "Fig1.pdf", metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(od / "Fig1.png", dpi=300)
    plt.close(fig)

    rec = dict(
        figure="Fig1", version="v5.5.1",
        source="02_OUT/extra_tissues_20261116.json :: sections[].alpha_by_d",
        panel_a=dict(grid="8x8 native bins", cell_um=2, block_um=16,
                     box="outline of the whole 8x8 field, no inner box",
                     scale_bars=["2 µm (one cell)", "16 µm (eight cells)"],
                     bin256="128 x 128 native bins = 16,384"),
        panel_b=dict(curves=[STYLE[k][1] for k in ORDER],
                     fold={STYLE[k][1]: round(fold[k], 3) for k in ORDER},
                     labels_at="right end of each curve, 256 µm"),
        files={f.name: dict(bytes=f.stat().st_size,
                            sha256=hashlib.sha256(f.read_bytes()).hexdigest())
               for f in (od / "Fig1.pdf", od / "Fig1.png")},
    )
    audit = od.parent / "audit" if od.name != "audit" else od
    audit.mkdir(parents=True, exist_ok=True)
    (audit / "fig1_v551.json").write_text(
        json.dumps(rec, ensure_ascii=False, indent=1), encoding="utf-8")

    print(f"written {od / 'Fig1.pdf'}")
    for k in ORDER:
        print(f"  {STYLE[k][1]:16s} fold={fold[k]:6.2f}  endpoint={y_end[k]:.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

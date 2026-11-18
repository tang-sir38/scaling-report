r"""124_fig_v57.py — Fig 1 与 Fig S3 的 v5.7 版：七片 → 八片。

改动
----
1. 新增 Visium HD human colon cancer P5（第三片结肠癌，fold 19.86，全组最大）。
   它与 colon 主片（P2）同组织、同平台、同流程，是"同组织重复能不能复现"的
   第二个独立对照。图上用**同一颜色 + 实心/空心圆**表示同组织，避免读者以为
   它们是两个组织。
2. Fig 1b 与 Fig S3 用同一份数据、同一顺序、同一配色，逐点对齐；两图曲线数
   同为 8，图注里的 "seven" 一律改 "eight"。
3. 所有 fold / span 现场从 alpha 数组算，不手抄；并与两份源文件的 fold 对表。

数据
----
  七片：02_OUT/extra_tissues_20261116.json :: sections[].alpha_by_d
  P5 ：E:/A_spatial_20261116/OUT/human_colon_cancer_p5/pooled_alpha_noMT.json :: per_d[*].alpha_pool

用法：python 124_fig_v57.py --outdir <04_图件 或 _build 目录>
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
P5JSON = Path(r"E:\A_spatial_20261116\OUT\human_colon_cancer_p5\pooled_alpha_noMT.json")
DS = [2, 4, 8, 16, 32, 64, 128, 256]

# (颜色, 图例, 空心/实心)
STYLE = {
    "colon":      ("#1B3A5C", "Colon (P2)", False),
    "colon_p5":   ("#1B3A5C", "Colon (P5)", True),
    "lung":       ("#4A6E8F", "Lung", False),
    "brain":      ("#B5651D", "Brain", False),
    "mouse_small_intestine": ("#6E8B3D", "Small intestine", False),
    "mouse_kidney": ("#8C4A7E", "Kidney", False),
    "human_pancreas": ("#A8452F", "Pancreas", False),
    "mouse_embryo": ("#3F7F6F", "Embryo", False),
}
# Fig 1b 按 fold 从大到小；Fig S3 按 2 µm 绝对值从大到小
ORDER_B = ["colon_p5", "mouse_small_intestine", "colon", "lung", "brain",
           "mouse_kidney", "mouse_embryo", "human_pancreas"]
ORDER_S3 = ["colon_p5", "colon", "lung", "mouse_small_intestine", "mouse_kidney",
            "mouse_embryo", "brain", "human_pancreas"]

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


def load_curves():
    summ = json.loads(SUMMARY.read_text(encoding="utf-8"))
    sec = {s["key"]: s for s in summ["sections"]}
    C, fold_declared = {}, {}
    for k in STYLE:
        if k == "colon_p5":
            continue
        C[k] = np.array([float(sec[k]["alpha_by_d"][str(d)]) for d in DS])
        fold_declared[k] = float(sec[k]["drop_fold"])
    p5 = json.loads(P5JSON.read_text(encoding="utf-8"))
    C["colon_p5"] = np.array([float(p5["per_d"][str(d)]["alpha_pool"]) for d in DS])
    fold_declared["colon_p5"] = C["colon_p5"][0] / C["colon_p5"][-1]
    fold = {k: float(C[k][0] / C[k][-1]) for k in C}
    for k in C:                          # 与源文件对表，防两处口径漂移
        assert abs(fold[k] - fold_declared[k]) < 1e-6, (k, fold[k], fold_declared[k])
    return C, fold


def style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)
    ax.tick_params(which="both", top=False, right=False)


def panel_a(axs):
    """8×8 原生格 = 16 µm；一块 = 8×8 格；只两个标尺。"""
    n = 8
    for i in range(n):
        for j in range(n):
            axs.add_patch(Rectangle((i, j), 1, 1, facecolor="#EDF2F7",
                                    edgecolor=LGREY, lw=0.5))
    axs.add_patch(Rectangle((0, 0), 1, 1, facecolor="#1B3A5C", edgecolor="none"))
    axs.add_patch(Rectangle((0, 0), n, n, facecolor="none",
                            edgecolor="#1B3A5C", lw=1.3, clip_on=False))

    for x0 in (0.0, 1.0):
        axs.plot([x0, x0], [-0.70, -0.42], color=MGREY, lw=0.5, clip_on=False)
    axs.plot([0.0, 1.0], [-0.56, -0.56], color=MGREY, lw=0.5, clip_on=False)
    axs.text(0.5, -0.78, "2 µm", ha="center", va="top", fontsize=6.2, color=MGREY)

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


def draw_fig1(od: Path, C, fold):
    fig = plt.figure(figsize=(W2, 74 * MM))
    axs = fig.add_axes([0.075, 0.20, 0.175, 0.60])
    ax = fig.add_axes([0.375, 0.155, 0.50, 0.745])
    panel_a(axs)

    dd = np.logspace(np.log10(1.9), np.log10(270), 60)
    ax.plot(dd, (dd / 2.0) ** -2, color=GREY, lw=0.9, ls=(0, (4, 2)), zorder=1)
    ax.plot(dd, np.ones_like(dd), color=GREY, lw=0.8, ls=(0, (1, 2)), zorder=1)
    for k in ORDER_B:
        col, lab, solid = STYLE[k]
        ax.plot(DS, C[k] / C[k][0], color=col, lw=1.0, marker="o", ms=2.3,
                mfc=(col if solid else "white"), mew=0.7, zorder=3, label=lab)

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

    fig.canvas.draw()
    inv = ax.transData.inverted()
    x_lab = 330.0
    y_end = {k: C[k][-1] / C[k][0] for k in ORDER_B}
    order_desc = sorted(ORDER_B, key=lambda k: -y_end[k])
    px_per_pt = fig.dpi / 72.0
    min_px = 9.2 * px_per_pt
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

    ax.text(-0.155, 1.03, "b", transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="bottom", ha="left")

    fig.savefig(od / "Fig1.pdf", metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(od / "Fig1.png", dpi=300)
    plt.close(fig)
    return y_end


def draw_figs3(od: Path, C):
    fig, ax = plt.subplots(figsize=(W2, 62 * MM))
    for k in ORDER_S3:
        col, lab, solid = STYLE[k]
        ax.plot(DS, C[k], color=col, lw=1.0, marker="o", ms=2.4,
                mfc=(col if solid else "white"), mew=0.7, label=lab, zorder=3)

    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xticks(DS)
    ax.set_xticklabels([str(d) for d in DS])
    ax.minorticks_off()
    ax.set_yticks([0.5, 1, 2, 5, 10, 20, 50])
    ax.set_yticklabels(["0.5", "1", "2", "5", "10", "20", "50"])
    ax.yaxis.set_minor_locator(NullLocator())
    ax.set_xlabel("Bin size (µm)")
    ax.set_ylabel("overdispersion α(d)")
    ax.set_xlim(1.7, 300)
    ax.set_ylim(0.42, 55)

    lo2, hi2 = min(C[k][0] for k in ORDER_S3), max(C[k][0] for k in ORDER_S3)
    lo256, hi256 = min(C[k][-1] for k in ORDER_S3), max(C[k][-1] for k in ORDER_S3)
    ax.annotate("", xy=(2.0, lo2), xytext=(2.0, hi2),
                arrowprops=dict(arrowstyle="<->", color=GREY, lw=0.5, shrinkA=0, shrinkB=0))
    ax.text(2.02, hi2 * 1.18, f"{hi2 / lo2:.0f}-fold", fontsize=6.4,
            color=GREY, va="bottom", ha="left")
    ax.annotate("", xy=(256, lo256), xytext=(256, hi256),
                arrowprops=dict(arrowstyle="<->", color=GREY, lw=0.5, shrinkA=0, shrinkB=0))
    ax.text(217, np.sqrt(lo256 * hi256), f"{hi256 / lo256:.1f}-fold", fontsize=6.4,
            color=GREY, va="center", ha="right")

    style(ax)
    ax.legend(ncol=4, loc="upper right", fontsize=6.4, columnspacing=0.9,
              handletextpad=0.4)

    fig.savefig(od / "FigS3.pdf", metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(od / "FigS3.png", dpi=300)
    plt.close(fig)

    from PIL import Image
    w, h = Image.open(od / "FigS3.png").size
    return dict(lo2=lo2, hi2=hi2, ratio2=hi2 / lo2,
                lo256=lo256, hi256=hi256, ratio256=hi256 / lo256, px=[w, h])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    od = Path(a.outdir)
    od.mkdir(parents=True, exist_ok=True)

    C, fold = load_curves()
    y_end = draw_fig1(od, C, fold)
    s3 = draw_figs3(od, C)

    rec = dict(figure="Fig1+FigS3", version="v5.7",
               n_sections=len(C), sections=list(C),
               source=dict(seven=str(SUMMARY), p5=str(P5JSON)),
               fold={k: round(fold[k], 3) for k in C},
               endpoint_norm={k: round(y_end[k], 5) for k in C},
               fig_s3=dict(span_2um=[s3["lo2"], s3["hi2"]], ratio_2um=round(s3["ratio2"], 2),
                           span_256um=[s3["lo256"], s3["hi256"]],
                           ratio_256um=round(s3["ratio256"], 2), png_px=s3["px"]),
               files={f.name: dict(bytes=f.stat().st_size,
                                   sha256=hashlib.sha256(f.read_bytes()).hexdigest())
                      for f in (od / "Fig1.pdf", od / "Fig1.png",
                                od / "FigS3.pdf", od / "FigS3.png")})
    audit = od.parent / "audit" if od.name != "audit" else od
    audit.mkdir(parents=True, exist_ok=True)
    (audit / "fig_v57.json").write_text(
        json.dumps(rec, ensure_ascii=False, indent=1), encoding="utf-8")

    print(f"written {od}/Fig1.pdf, FigS3.pdf   ({len(C)} sections)")
    for k in ORDER_B:
        print(f"  {STYLE[k][1]:16s} fold={fold[k]:6.2f}  norm_end={y_end[k]:.4f}")
    print(f"  FigS3  2µm span {s3['lo2']:.3f}-{s3['hi2']:.3f} = {s3['ratio2']:.1f}x"
          f"   |  256µm span {s3['lo256']:.3f}-{s3['hi256']:.3f} = {s3['ratio256']:.2f}x")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

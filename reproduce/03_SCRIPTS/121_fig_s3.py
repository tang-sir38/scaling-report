#!/usr/bin/env python3
r"""
121_fig_s3.py — Figure S3：七个组织的 α(d) 曲线，画在与 Fig. 1b 相同的坐标上。

为什么需要这张图
----------------
正文的标题现在写的是 "5- to 15-fold"，这个区间来自七个组织；而 Fig. 1b 只画了
三片（colon / lung / brain，那是逐基因与 covariogram 分析的三片）。读者看到
"seven sections" 却没有一张画全的图，只能去 Table S1 里自己拼。这张图补上载体：
同样的归一化（α(d)/α(2 µm)）、同样的 log–log 轴、同样的两条参考线，
所以它能与 Fig. 1b 直接对照。

样式沿用 43_figures_v2.py：Arial 7.5 pt、180 mm 双栏、无网格、PDF 冻结时间戳、
PNG 300 dpi；前三片的颜色不变（#1B3A5C / #4A6E8F / #B5651D），新增四片各给一色。

用法
----
  python 121_fig_s3.py --outdir <包的 04_图件>
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

MM = 1 / 25.4
W2 = 180 * MM

REPRO = Path(__file__).resolve().parents[1]
LOCK = REPRO / "02_OUT" / "NUMBERS_LOCK_v7.json"
SUMMARY = REPRO / "02_OUT" / "extra_tissues_20261116.json"
WS = Path(r"E:\A_spatial_20261116")

DS = [2, 4, 8, 16, 32, 64, 128, 256]

# 颜色与标签：前三片沿用主图配色，新增四片按色相拉开、印刷可辨
STYLE = {
    "colon": ("#1B3A5C", "Colon"),
    "lung": ("#4A6E8F", "Lung"),
    "brain": ("#B5651D", "Brain"),
    "mouse_small_intestine": ("#6E8B3D", "Small intestine"),
    "mouse_kidney": ("#8C4A7E", "Kidney"),
    "human_pancreas": ("#A8452F", "Pancreas"),
    "mouse_embryo": ("#3F7F6F", "Embryo"),
}
GREY = "#8A8A8A"

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


def curves(path: Path) -> dict[str, list[float]]:
    """整条 α(d) 曲线只从汇总文件取，不回头读原始 OUT 目录。"""
    summary = json.loads(path.read_text(encoding="utf-8"))
    return {s["key"]: [float(s["alpha_by_d"][str(d)]) for d in DS]
            for s in summary["sections"] if s["key"] != "colon_p1"}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--summary", default=str(SUMMARY))
    a = ap.parse_args()

    C = curves(Path(a.summary))
    order = ["brain", "mouse_small_intestine", "lung", "colon",
             "mouse_kidney", "human_pancreas", "mouse_embryo"]
    assert set(order) == set(C), f"缺曲线: {set(order) ^ set(C)}"

    fig, ax = plt.subplots(figsize=(W2, 58 * MM))
    dd = np.array(DS, dtype=float)
    ax.plot(dd, (dd / 2.0) ** -2, color=GREY, lw=0.9, ls=(0, (4, 2)), zorder=1)
    ax.plot(dd, np.ones_like(dd), color=GREY, lw=0.8, ls=(0, (1, 2)), zorder=1)

    for k in order:
        col, lab = STYLE[k]
        v = np.array(C[k], dtype=float)
        ax.plot(DS, v / v[0], color=col, lw=1.1, marker="o", ms=2.4,
                mec="white", mew=0.3, label=lab, zorder=3)

    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xticks(DS)
    ax.set_xticklabels([str(d) for d in DS])
    ax.minorticks_off()
    ax.set_xlabel("Bin size (\u00b5m)")
    ax.set_ylabel("\u03b1(d) / \u03b1(2 \u00b5m)")
    ax.set_xlim(1.7, 300)
    ax.set_ylim(0.008, 1.6)

    # 参考线的标签要跟着线走，否则斜率和屏幕角度对不上；角度用显示坐标算
    fig.canvas.draw()
    p1 = ax.transData.transform((12.0, (12.0 / 2.0) ** -2))
    p2 = ax.transData.transform((48.0, (48.0 / 2.0) ** -2))
    ang = np.degrees(np.arctan2(p2[1] - p1[1], p2[0] - p1[0]))
    ax.text(9.5, (9.5 / 2.0) ** -2 * 1.35, r"$\alpha \propto d^{-2}$", fontsize=6.8,
            color=GREY, ha="left", va="bottom", rotation=ang,
            rotation_mode="anchor")
    ax.text(150, 1.06, r"$\alpha \propto d^{0}$", fontsize=6.8, color=GREY,
            ha="right", va="bottom")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(which="both", top=False, right=False)
    ax.legend(ncol=2, loc="lower left", fontsize=6.6, columnspacing=1.0)

    od = Path(a.outdir)
    od.mkdir(parents=True, exist_ok=True)
    fig.savefig(od / "FigS3.pdf", metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(od / "FigS3.png", dpi=300)

    import hashlib
    rec = dict(
        figure="FigS3",
        source="02_OUT/extra_tissues_20261116.json :: sections[].alpha_by_d",
        normalisation="alpha(d) / alpha(2um), log-log, same axes and guides as Fig. 1b",
        sections={k: dict(alpha_2um=C[k][0], alpha_256um=C[k][-1],
                          fall=C[k][0] / C[k][-1]) for k in order},
        files={f.name: dict(bytes=f.stat().st_size,
                            sha256=hashlib.sha256(f.read_bytes()).hexdigest())
               for f in (od / "FigS3.pdf", od / "FigS3.png")},
    )
    audit = od.parent / "audit"
    audit.mkdir(parents=True, exist_ok=True)
    (audit / "fig_s3.json").write_text(
        json.dumps(rec, ensure_ascii=False, indent=1), encoding="utf-8")

    print(f"written {od / 'FigS3.pdf'}")
    print(f"written {od / 'FigS3.png'}")
    for k in order:
        v = np.array(C[k], dtype=float)
        print(f"  {STYLE[k][1]:16s} alpha(2um)={v[0]:7.3f}  alpha(256)={v[-1]:6.3f}  "
              f"fall={v[0] / v[-1]:6.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

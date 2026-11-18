r"""126_fig_s3_v551.py — Figure S3 改版（v5.5.1）。

为什么要改
----------
v5.5 的 Fig S3 画的是七片 α(d)/α(2 µm)，用途是"Fig 1b 只画三片，七片在这里"。
v5.5.1 的 Fig 1b 改成七条曲线并标了倍数，S3 与它逐点重合，成了重复图。

改法：S3 改画**未归一化**的 α(d)。这不是换包装，是换问题：
  - Fig 1b 归一化后只能看斜率；
  - S3 看绝对水平——七片在 2 µm 处相差 21 倍，在 256 µm 处只差 3.8 倍，
    因为陡的几条落下来正好落在平的那几条附近。
归一化图看不出这一点，绝对图一眼看出来。这也正面支持正文那句
"level is not comparable between sections, the fall is"。

数据同 Fig 1b：02_OUT/extra_tissues_20261116.json :: alpha_by_d。
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
ORDER = ["colon", "lung", "brain", "mouse_small_intestine",
         "mouse_kidney", "human_pancreas", "mouse_embryo"]
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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--summary", default=str(SUMMARY))
    a = ap.parse_args()

    summ = json.loads(Path(a.summary).read_text(encoding="utf-8"))
    sec = {s["key"]: s for s in summ["sections"]}
    C = {k: np.array([float(sec[k]["alpha_by_d"][str(d)]) for d in DS]) for k in ORDER}

    fig, ax = plt.subplots(figsize=(W2, 62 * MM))
    for k in ORDER:
        col, lab = STYLE[k]
        ax.plot(DS, C[k], color=col, lw=1.0, marker="o", ms=2.4, mfc="white",
                mew=0.7, label=lab, zorder=3)

    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xticks(DS)
    ax.set_xticklabels([str(d) for d in DS])
    ax.minorticks_off()
    ax.set_yticks([0.5, 1, 2, 5, 10, 20])
    ax.set_yticklabels(["0.5", "1", "2", "5", "10", "20"])
    ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlabel("Bin size (µm)")
    ax.set_ylabel("overdispersion α(d)")
    ax.set_xlim(1.7, 300)
    ax.set_ylim(0.42, 42)

    # 两端各画一道括注：跨度 21 倍 → 3.8 倍，这是本图唯一要读者带走的数
    lo2, hi2 = min(C[k][0] for k in ORDER), max(C[k][0] for k in ORDER)
    lo256, hi256 = min(C[k][-1] for k in ORDER), max(C[k][-1] for k in ORDER)
    ax.annotate("", xy=(2.0, lo2), xytext=(2.0, hi2),
                arrowprops=dict(arrowstyle="<->", color=GREY, lw=0.5, shrinkA=0, shrinkB=0))
    ax.text(2.02, hi2 * 1.18, f"{hi2 / lo2:.0f}-fold", fontsize=6.4,
            color=GREY, va="bottom", ha="left")
    ax.annotate("", xy=(256, lo256), xytext=(256, hi256),
                arrowprops=dict(arrowstyle="<->", color=GREY, lw=0.5, shrinkA=0, shrinkB=0))
    ax.text(217, np.sqrt(lo256 * hi256), f"{hi256 / lo256:.1f}-fold", fontsize=6.4,
            color=GREY, va="center", ha="right")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(which="both", top=False, right=False)
    ax.legend(ncol=4, loc="upper right", fontsize=6.4, columnspacing=0.9,
              handletextpad=0.4)

    od = Path(a.outdir)
    od.mkdir(parents=True, exist_ok=True)
    fig.savefig(od / "FigS3.pdf", metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(od / "FigS3.png", dpi=300)
    plt.close(fig)

    from PIL import Image
    w, h = Image.open(od / "FigS3.png").size
    rec = dict(
        figure="FigS3", version="v5.5.1",
        source="02_OUT/extra_tissues_20261116.json :: sections[].alpha_by_d",
        scale="absolute alpha(d), log-log; v5.5 plotted alpha(d)/alpha(2um)",
        span_2um=dict(lo=lo2, hi=hi2, ratio=hi2 / lo2),
        span_256um=dict(lo=lo256, hi=hi256, ratio=hi256 / lo256),
        png_px=[w, h],
        files={f.name: dict(bytes=f.stat().st_size,
                            sha256=hashlib.sha256(f.read_bytes()).hexdigest())
               for f in (od / "FigS3.pdf", od / "FigS3.png")},
    )
    audit = od.parent / "audit" if od.name != "audit" else od
    audit.mkdir(parents=True, exist_ok=True)
    (audit / "fig_s3_v551.json").write_text(
        json.dumps(rec, ensure_ascii=False, indent=1), encoding="utf-8")

    print(f"written {od / 'FigS3.pdf'}  png {w}x{h}")
    print(f"  2 µm span   {lo2:.3f}-{hi2:.3f}  = {hi2 / lo2:.1f}x")
    print(f"  256 µm span {lo256:.3f}-{hi256:.3f}  = {hi256 / lo256:.2f}x")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

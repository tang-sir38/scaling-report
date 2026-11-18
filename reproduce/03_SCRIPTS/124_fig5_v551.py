r"""124_fig5_v551.py — M1 的 Fig 5 重画（v5.5.1）。

只改 panel a 的 95% range 画法，其余与 v5.5 的 Fig 5（99_fig45_v4.py::fig5）逐像素同源。

问题：v5.5 的图注写 "the horizontal bar spans the 95% range over permutations"，
但读者看到的那条长横线其实是 observed→permuted 的**连线**；真正的 95% range
（20 次置换的 2.5–97.5 分位）宽度只有 0.003–0.010，画在 lw=2.6 的点下面，
被点完全盖住。图注描述的不是读者能看见的对象。

改法（两条都做，不是二选一）：
  1. 把 95% range 画成穿过实心点的**短杠 + 上下端帽**，zorder 在点之上，
     即使横向宽度不足 1 pt，端帽仍从点外露出，range 成为可见对象；
  2. 面板内加一行字给最大宽度，图注里的说法才有数可查。

数据仍取 02_OUT/excise/excise_colon.json 与 multigene_permute.json。
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

ROOT = Path(__file__).resolve().parents[1]
C_COLON = "#2C5F7C"
C_DARK = "#333333"
C_GREY = "#8C8C8C"

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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    ex = json.loads((ROOT / "02_OUT/excise/excise_colon.json").read_text(encoding="utf-8"))
    mp = json.loads((ROOT / "02_OUT/excise/multigene_permute.json").read_text(encoding="utf-8"))

    fig = plt.figure(figsize=(7.2, 2.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.05], wspace=0.34,
                          left=0.115, right=0.975, top=0.90, bottom=0.17)

    # ---- (a) 六基因哑铃 ----
    ax = fig.add_subplot(gs[0, 0])
    genes = list(ex["genes"].keys())
    order = sorted(genes, key=lambda g: ex["genes"][g]["baseline"]["gamma"])
    widths = []
    for i, gname in enumerate(order):
        v = ex["genes"][gname]
        b = v["baseline"]["gamma"]
        p = v["permute"]["gamma_median"]
        lo, hi = v["permute"]["gamma_lo"], v["permute"]["gamma_hi"]
        widths.append(hi - lo)
        ax.plot([b, p], [i, i], color=C_GREY, lw=0.7, zorder=1)
        ax.plot(b, i, "o", ms=3.4, mfc="white", mec=C_DARK, mew=0.7, zorder=3)
        ax.plot(p, i, "o", ms=3.4, color=C_COLON, zorder=3)
        # 95% range：短杠 + 上下端帽，画在点之上，宽度不足 1 pt 时端帽仍可见
        h = 0.20
        ax.plot([lo, hi], [i, i], color=C_DARK, lw=0.7, zorder=4,
                solid_capstyle="butt")
        ax.plot([lo, lo], [i - h, i + h], color=C_DARK, lw=0.7, zorder=4)
        ax.plot([hi, hi], [i - h, i + h], color=C_DARK, lw=0.7, zorder=4)
        ax.text(1.82, i, gname, fontsize=6.0, va="center", ha="left")
    ax.set_yticks([])
    ax.set_xlim(0.2, 1.90)
    ax.set_ylim(-1.25, len(order) - 0.3)
    ax.set_xlabel(r"per-gene exponent  $\gamma_g$")
    ax.text(0.03, 0.965, "a", transform=ax.transAxes, fontsize=9, va="top")
    ax.plot([], [], "o", ms=3.4, mfc="white", mec=C_DARK, mew=0.7, label="observed")
    ax.plot([], [], "o", ms=3.4, color=C_COLON, label="permuted")
    ax.plot([], [], color=C_DARK, lw=0.7, marker="|", ms=5,
            label="95% range")
    ax.legend(frameon=False, loc="lower right", handletextpad=0.4,
              bbox_to_anchor=(0.88, 0.005), borderaxespad=0.0,
              labelspacing=0.30, fontsize=5.9)
    wmax = max(widths)

    # ---- (b) 剂量-反应 ----
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

    od = Path(a.outdir)
    od.mkdir(parents=True, exist_ok=True)
    fig.savefig(od / "Fig5.pdf", metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(od / "Fig5.png", dpi=400)
    plt.close(fig)

    rec = dict(
        figure="Fig5", version="v5.5.1",
        source="02_OUT/excise/excise_colon.json + 02_OUT/excise/multigene_permute.json",
        panel_a=dict(
            connector="grey line joins observed to the mean over 20 permutations",
            range_drawn="bar with end caps through the filled circle, zorder above the marker",
            range_width={g: round(ex["genes"][g]["permute"]["gamma_hi"]
                                  - ex["genes"][g]["permute"]["gamma_lo"], 4) for g in order},
            range_width_max=round(wmax, 4),
        ),
        files={f.name: dict(bytes=f.stat().st_size,
                            sha256=hashlib.sha256(f.read_bytes()).hexdigest())
               for f in (od / "Fig5.pdf", od / "Fig5.png")},
    )
    audit = od.parent / "audit" if od.name != "audit" else od
    audit.mkdir(parents=True, exist_ok=True)
    (audit / "fig5_v551.json").write_text(
        json.dumps(rec, ensure_ascii=False, indent=1), encoding="utf-8")

    print(f"written {od / 'Fig5.pdf'}")
    for g in order:
        v = ex["genes"][g]["permute"]
        print(f"  {g:8s} range {v['gamma_lo']:.3f}-{v['gamma_hi']:.3f} "
              f"width {v['gamma_hi'] - v['gamma_lo']:.4f}")
    print(f"  max width {wmax:.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

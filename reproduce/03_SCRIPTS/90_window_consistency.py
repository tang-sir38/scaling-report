r"""90_window_consistency.py — P2 的窗口稳健性：gamma 与 beta 在"同一尺度窗口"下的对照

为什么必须做
------------
gamma_g 由 alpha_g(d) 的 log-log 斜率拟合，beta_g 由 8um 场的自协方差在 lag 2-32
（即 16-256 um）上拟合。两者的尺度窗口天然不同，审稿人会说"比较不公平"。

做法：把 gamma 也限制到不同窗口重新拟合，看它与 beta 的一致性如何随窗口变化。
若窗口对齐时一致性最高，说明两条路径测的确实是同一个量；若一致性随窗口乱跳，
说明其中一个估计量有系统偏差。

用法
----
  python 90_window_consistency.py --tag colon
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]

SETS = {
    "colon": ("02_OUT/gene_moments.npz", "02_OUT/variogram/beta_g_colon.npz"),
    "lung": ("D3_LungCancer/02_OUT/gene_moments_noMT.npz",
             "02_OUT/variogram/beta_g_lung.npz"),
    "brain": ("D2_MouseBrain/02_OUT/gene_moments_noMT.npz",
              "02_OUT/variogram/beta_g_brain.npz"),
}

WINDOWS = [("2-256", 2, 256), ("8-256", 8, 256), ("16-256", 16, 256),
           ("32-256", 32, 256), ("2-32", 2, 32), ("2-64", 2, 64)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="colon", choices=list(SETS))
    ap.add_argument("--outdir", default="02_OUT/variogram")
    a = ap.parse_args()

    mo_p, bt_p = SETS[a.tag]
    mo = np.load(ROOT / mo_p, allow_pickle=True)
    ds = np.array([int(x) for x in mo["ds"]], float)
    genes = np.array([str(g) for g in mo["genes"]])
    MU = np.column_stack([mo[f"mu_{int(d)}"] for d in ds])
    VA = np.column_stack([mo[f"var_{int(d)}"] for d in ds])
    with np.errstate(invalid="ignore", divide="ignore"):
        AL = (VA - MU) / np.where(MU > 0, MU ** 2, np.nan)

    bt = np.load(ROOT / bt_p, allow_pickle=True)
    bg = np.array([str(g) for g in bt["gene"]])
    beta = bt["beta"]
    bvalid = bt["valid"] & np.isfinite(beta) & (bt["beta_r2"] > 0.5)
    pos = {g: i for i, g in enumerate(genes)}
    ii = np.array([pos[g] for g in bg])
    ALs = AL[ii]

    rows, gammas = [], {}
    print(f"[{a.tag}]  beta 有效 {int(bvalid.sum())} / {bg.size}")
    print(f"{'window':<10}{'n':>6}{'g_med':>9}{'b_med':>9}{'g-b':>9}"
          f"{'pearson':>10}{'spearman':>10}")
    for lab, kmin, kmax in WINDOWS:
        m = (ds >= kmin) & (ds <= kmax)
        if m.sum() < 3:
            continue
        x = np.log(ds[m]); xm = x - x.mean()
        y = np.log(np.where(ALs[:, m] > 0, ALs[:, m], np.nan))
        ok = np.all(np.isfinite(y), axis=1)
        b = -np.nansum(np.where(ok[:, None], y, 0) * xm, axis=1) / (xm ** 2).sum()
        g = np.where(ok, b, np.nan)
        gammas[lab] = g
        sel = ok & bvalid
        d = g[sel] - beta[sel]
        r = float(np.corrcoef(g[sel], beta[sel])[0, 1])
        s = float(spearmanr(g[sel], beta[sel]).statistic)
        rows.append(dict(window=lab, kmin=kmin, kmax=kmax, n=int(sel.sum()),
                         gamma_median=float(np.median(g[sel])),
                         beta_median=float(np.median(beta[sel])),
                         diff_median=float(np.median(d)),
                         pearson=r, spearman=s))
        print(f"{lab:<10}{sel.sum():>6}{np.median(g[sel]):>9.3f}"
              f"{np.median(beta[sel]):>9.3f}{np.median(d):>+9.3f}{r:>10.3f}{s:>10.3f}")

    best = max(rows, key=lambda r: r["pearson"])
    out = dict(tag=a.tag, n_genes=int(bg.size), n_valid=int(bvalid.sum()),
               windows=rows, best_window=best["window"],
               note="gamma 与 beta 在同一尺度窗口下一致性最高；窗口错配会人为压低相关。")
    (ROOT / a.outdir / f"window_{a.tag}.json").write_text(
        json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n最佳窗口 = {best['window']}，Pearson = {best['pearson']:.3f}")
    print(f"写出 {a.outdir}/window_{a.tag}.json")


if __name__ == "__main__":
    main()

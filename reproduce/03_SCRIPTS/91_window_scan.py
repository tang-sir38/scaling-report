r"""91_window_scan.py — 标度窗口扫描：gamma 与 beta 在哪个尺度区间内相等

背景
----
90_window_consistency.py 发现"窗口对齐后 gamma 与 beta 相关升到 0.9"。
但 87 的 beta 是在 8um 场的 lag 2-24 上拟合的，而 alpha 的 gamma 用的是 2-256um。
若两套数据的有效标度区不同，笼统比较会得出"某组织不一致"的错误结论。

本脚本对**两侧窗口同时扫描**：
  gamma 侧：alpha(d) 的拟合区间
  beta  侧：rho(h) 的拟合 lag 区间
输出一致性矩阵，并给出"该组织的有效标度区间"。

判据（先验，不看结果调参）
--------------------------
  有效标度区 := 使 |median(gamma) - median(beta)| 最小、且 Pearson 最高的窗口。
  若不同组织的有效区不同，说明幂律区间的宽度本身是组织属性——这是结论，不是问题。

物理背景
--------
README 已记录：总量场变差函数在 ~150 um 饱和，那是细胞密度尺度。
若单基因协方差在大 lag 处也受密度尺度污染，beta 会随 lag 增大而升高。
本扫描正是要检测这一点。

用法
----
  python 91_window_scan.py --tag brain
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

GAMMA_WINDOWS = [(2, 256), (8, 256), (16, 256), (2, 32), (2, 64), (2, 128)]
BETA_WINDOWS = [(2, 8), (2, 12), (2, 16), (2, 24), (2, 32), (3, 16)]


def fit_beta(rho, lags, kmin, kmax):
    m = np.isfinite(rho) & (rho > 0) & (lags >= kmin) & (lags <= kmax)
    if m.sum() < 3:
        return np.nan, np.nan
    x, y = np.log(lags[m]), np.log(rho[m])
    A = np.vstack([x, np.ones_like(x)]).T
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ c
    ss = ((y - pred) ** 2).sum()
    tot = ((y - y.mean()) ** 2).sum()
    return float(-c[0]), float(1 - ss / tot if tot > 0 else np.nan)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="brain", choices=list(SETS))
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
    rho = bt["rho"]; lags = bt["lags"].astype(float)
    gam0 = bt["gamma"]; valid = bt["valid"]
    pos = {g: i for i, g in enumerate(genes)}
    ii = np.array([pos[g] for g in bg])
    ALs = AL[ii]

    # ---- gamma 侧
    gam = {}
    for kmin, kmax in GAMMA_WINDOWS:
        m = (ds >= kmin) & (ds <= kmax)
        x = np.log(ds[m]); xm = x - x.mean()
        y = np.log(np.where(ALs[:, m] > 0, ALs[:, m], np.nan))
        ok = np.all(np.isfinite(y), axis=1)
        b = -np.nansum(np.where(ok[:, None], y, 0) * xm, axis=1) / (xm ** 2).sum()
        gam[(kmin, kmax)] = np.where(ok, b, np.nan)

    # ---- beta 侧
    bet = {}
    for kmin, kmax in BETA_WINDOWS:
        out = np.full(bg.size, np.nan)
        for i in range(bg.size):
            out[i], _ = fit_beta(rho[i], lags, kmin, kmax)
        bet[(kmin, kmax)] = out
    # 参考：87 主口径
    bet[(2, 24)] = bt["beta"]

    rows = []
    print(f"[{a.tag}]  n={bg.size}  一致性扫描（行=gamma 窗口，列=beta lag 窗口）")
    hdr = "  ".join(f"{a2}-{b2:<4}" for _, a2, b2 in
                    [("", w[0], w[1]) for w in GAMMA_WINDOWS])
    print(f"{'gamma\\beta':<12}" + "".join(f"{f'{w[0]}-{w[1]}':>12}" for w in BETA_WINDOWS))
    for gw in GAMMA_WINDOWS:
        line = f"{f'{gw[0]}-{gw[1]}':<12}"
        for bw in BETA_WINDOWS:
            m = valid & np.isfinite(gam[gw]) & np.isfinite(bet[bw])
            if m.sum() < 10:
                line += f"{'--':>12}"
                continue
            r = float(np.corrcoef(gam[gw][m], bet[bw][m])[0, 1])
            line += f"{r:>12.3f}"
            rows.append(dict(gamma_window=list(gw), beta_window=list(bw),
                             n=int(m.sum()),
                             gamma_median=float(np.median(gam[gw][m])),
                             beta_median=float(np.median(bet[bw][m])),
                             diff=float(np.median(gam[gw][m] - bet[bw][m])),
                             pearson=r,
                             spearman=float(spearmanr(gam[gw][m], bet[bw][m]).statistic)))
        print(line)

    # 判定：|diff| 最小 + Pearson 最高
    best_r = max(rows, key=lambda r: r["pearson"])
    best_d = min(rows, key=lambda r: abs(r["diff"]))
    print(f"\n最高相关: gamma {best_r['gamma_window']} vs beta {best_r['beta_window']}  "
          f"r={best_r['pearson']:.3f}  diff={best_r['diff']:+.3f}")
    print(f"最小偏差: gamma {best_d['gamma_window']} vs beta {best_d['beta_window']}  "
          f"r={best_d['pearson']:.3f}  diff={best_d['diff']:+.3f}")

    # beta 随 lag 上界的漂移（检验密度尺度污染）
    drift = []
    for bw in sorted([w for w in BETA_WINDOWS], key=lambda w: w[1]):
        m = valid & np.isfinite(bet[bw])
        if m.sum() > 10:
            drift.append(dict(window=list(bw), beta_median=float(np.median(bet[bw][m]))))
    print("\nbeta 随 lag 上界的漂移（单调上升 = 大 lag 被密度尺度污染）：")
    for d in drift:
        print(f"  lag {d['window'][0]}-{d['window'][1]:<3} beta={d['beta_median']:.3f}")

    out = dict(tag=a.tag, n=int(bg.size),
               gamma_windows=[list(w) for w in GAMMA_WINDOWS],
               beta_windows=[list(w) for w in BETA_WINDOWS],
               grid=rows, best_by_pearson=best_r, best_by_absdiff=best_d,
               beta_drift_with_lag=drift)
    (ROOT / a.outdir / f"windowscan_{a.tag}.json").write_text(
        json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n写出 {a.outdir}/windowscan_{a.tag}.json")


if __name__ == "__main__":
    main()

r"""94_window_reliability.py — 用两半可靠性选 gamma 的拟合窗口（与 beta 完全无关）

判据
----
候选窗口 W。对每个 W：
  两半各自算 gamma_g(W)，比较两半之间的 Pearson r = 可靠性。
可靠性只由估计噪声决定，与空间协方差结构 beta 无任何关系，
因此"取可靠性达标的窗口"不构成挑结果。

规则（先验，写死）
------------------
  取**最宽**（下界最小）的窗口，使三套组织的可靠性都 >= 0.70。
  宽窗口优先是因为它保留了最大的动态范围；
  可靠性是唯一门槛，不涉及 gamma 的数值、也不涉及 beta。

附：Spearman-Brown 校正 r_full = 2r/(1+r)，把"半样本可靠性"换成"全样本可靠性"。

用法
----
  python 94_window_reliability.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[1]
DS = [2, 4, 8, 16, 32, 64, 128, 256]
TAGS = ["colon", "lung", "brain"]
WINDOWS = [(2, 256), (2, 128), (2, 64), (8, 256), (8, 128), (16, 256), (16, 128), (32, 256)]
THRESH = 0.70


def fit_gamma_win(al, ds_win):
    j = np.array([DS.index(int(d)) for d in ds_win])
    sub = al[:, j]
    ok = np.all(np.isfinite(sub) & (sub > 0), axis=1)
    x = np.log(np.array(ds_win, float))
    xm = x - x.mean()
    y = np.log(np.where(ok[:, None], sub, 1.0))
    ym = y - y.mean(axis=1, keepdims=True)
    g = -np.sum(xm[None, :] * ym, axis=1) / np.sum(xm ** 2)
    out = np.full(al.shape[0], np.nan)
    out[ok] = g[ok]
    return out, ok


def main():
    allrows = {}
    for tag in TAGS:
        p = ROOT / f"02_OUT/variogram/splithalf_alpha_{tag}.npz"
        if not p.exists():
            print(f"缺 {p}，跳过 {tag}")
            continue
        z = np.load(p, allow_pickle=True)
        A, B = z["al_a"].astype(np.float64), z["al_b"].astype(np.float64)
        MA, MB = z["mu_a"].astype(np.float64), z["mu_b"].astype(np.float64)
        rows = []
        print("=" * 88)
        print(f"[{tag}]  两半可靠性 vs 拟合窗口（可靠性只看噪声，与 beta 无关）")
        print(f"{'window':>9} {'npts':>5} {'n_both':>7} {'r_half':>7} {'r_SB':>7} "
              f"{'med_g_a':>8} {'med_g_b':>8} {'MAD':>7} {'med_mu@lo':>10}")
        for w in WINDOWS:
            ga, oka = fit_gamma_win(A, w)
            gb, okb = fit_gamma_win(B, w)
            both = oka & okb
            if both.sum() < 100:
                continue
            r = float(pearsonr(ga[both], gb[both]).statistic)
            rho = float(spearmanr(ga[both], gb[both]).statistic)
            rsb = 2 * r / (1 + r) if r > -1 else np.nan
            lo_j = DS.index(int(w[0]))
            mum = float(np.median(0.5 * (MA[both, lo_j] + MB[both, lo_j])))
            rows.append(dict(window=list(w), n_points=len(w), n_both=int(both.sum()),
                             r_half=r, r_spearman=rho, r_spearman_brown=float(rsb),
                             gamma_median_a=float(np.median(ga[both])),
                             gamma_median_b=float(np.median(gb[both])),
                             mad=float(np.median(np.abs(ga[both] - gb[both]))),
                             median_mu_at_low_d=mum))
            print(f"{f'{w[0]}-{w[1]}':>9} {len(w):>5} {int(both.sum()):>7} {r:>7.3f} "
                  f"{rsb:>7.3f} {np.median(ga[both]):>8.3f} {np.median(gb[both]):>8.3f} "
                  f"{rows[-1]['mad']:>7.3f} {mum:>10.5f}")
        allrows[tag] = rows

    # ---- 规则：最宽且三套可靠性都 >= 阈值的窗口
    print("\n" + "=" * 88)
    print(f"判据：取最宽（下界最小、跨度最大）且三套组织 r_half 均 >= {THRESH} 的窗口")
    verdict = None
    for w in WINDOWS:
        vals = []
        for tag in TAGS:
            hit = [r for r in allrows.get(tag, []) if tuple(r["window"]) == w]
            if not hit:
                vals = None
                break
            vals.append(hit[0]["r_half"])
        if vals is None:
            continue
        okall = all(v >= THRESH for v in vals)
        span = np.log2(w[1] / w[0])
        print(f"  {str(w):>10}  跨度 {span:.0f} oct  r = "
              + ", ".join(f"{v:.3f}" for v in vals) + ("   <= 通过" if okall else ""))
        if okall and verdict is None:
            verdict = dict(window=list(w), span_octaves=float(span),
                           r_half=dict(zip(TAGS, [float(v) for v in vals])))

    print(f"\n结论：gamma 的可靠性合格窗口 = {verdict['window'] if verdict else '无'}")
    out = dict(threshold=THRESH, rule="widest window with r_half>=threshold in all three tissues",
               windows=[list(w) for w in WINDOWS], by_tissue=allrows, verdict=verdict)
    (ROOT / "02_OUT/variogram/window_reliability.json").write_text(
        json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print("写出 02_OUT/variogram/window_reliability.json")


if __name__ == "__main__":
    main()

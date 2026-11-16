r"""96_gamma_from_rhobar.py — 由 rho(h) 算块平均 rhō 的衰减率，与实际 gamma 比

为什么不是"gamma vs beta"
--------------------------
gamma 是 **alpha(d) 的衰减率**，而 alpha(d) = CV^2_0 * rhobar(d)，
所以 gamma = -dlog rhobar / dlog d。rhobar 是块内平均相关函数。
beta 是 **rho(h) 的局部衰减率**，即 -dlog rho / dlog h。

两者是不同的泛函：
  rhobar(k) = sum_{a,b<k} m(a,b) rho(sqrt(a^2+b^2)) / k^4
块内配对数近似正比于 k^2（每个间距都有一份），所以 rhobar 把 rho 在各个间距上
**摊平**。因此 rho 变陡时 rhobar 只轻微变陡。
只有当 rho 本身是纯幂律时，gamma 才等于 beta。

实测正是这样（91 + 本文档的局部斜率）：
  colon / lung：rho(h) 的局部斜率 8-256um 全程 0.20-0.50，无趋势 → 纯幂律 → gamma = beta
  brain：       局部斜率 0.40(8-16um) -> 0.65(72-88um) -> 0.9(128um) -> 1.5(232-248um) → 不是幂律
但 brain 的 alpha 曲线却是三套里最干净的幂律（pooled 残差 ±0.03）。
这不是矛盾，正是"摊平"的必然结果。

所以本脚本检验的是**更普适的那条**：
    gamma_pred = -slope[ log rhobar(k) ~ log d ]，k = d/8
对三套数据都应成立，包括脑。

近似
----
rho(h) 只在 8um 网格的整数 lag 上有值，且每个 lag 值本身已被 8um 块平滑过一次。
用它做块平均 = 二次平滑。对**斜率**的影响远小于对**绝对值**的影响
（二次平滑近似把有效尺度放大一个常数倍，只移动截距不改变斜率），
这也是本脚本只比斜率、不比绝对高度的原因。

用法
----
  python 96_gamma_from_rhobar.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
SETS = {
    "colon": ("02_OUT/gene_moments.npz", "02_OUT/variogram/beta_g_colon.npz"),
    "lung": ("D3_LungCancer/02_OUT/gene_moments_noMT.npz", "02_OUT/variogram/beta_g_lung.npz"),
    "brain": ("D2_MouseBrain/02_OUT/gene_moments_noMT.npz", "02_OUT/variogram/beta_g_brain.npz"),
}
# k=2,4,8,16 对应 d=16,32,64,128，块角最大 lag=sqrt(15^2+15^2)=21.2 < 32，安全
KS = [2, 4, 8, 16]
DS_FIT = [16, 32, 64, 128]


def rhobar(rho_row, k, max_lag):
    tot = 0.0
    cnt = 0
    for a in range(k):
        for b in range(k):
            m = 1 if (a == 0 and b == 0) else (2 if (a == 0 or b == 0) else 4)
            h = int(round(np.hypot(a, b)))
            if h > max_lag or not np.isfinite(rho_row[h]):
                return np.nan, cnt
            tot += m * rho_row[h]
            cnt += m
    return tot / cnt, cnt


def main():
    out = {}
    print("=" * 84)
    print(f"{'tag':>6} {'n':>6} {'g_obs':>7} {'g_pred':>7} {'diff':>7} "
          f"{'r_log':>7} {'rho_s':>7} {'slope':>7}")
    for tag, (mo_p, bt_p) in SETS.items():
        mo = np.load(ROOT / mo_p, allow_pickle=True)
        ds = np.array([int(x) for x in mo["ds"]], float)
        genes_mo = np.array([str(g) for g in mo["genes"]])
        MU = np.column_stack([mo[f"mu_{int(d)}"] for d in ds])
        VA = np.column_stack([mo[f"var_{int(d)}"] for d in ds])
        with np.errstate(invalid="ignore", divide="ignore"):
            AL = (VA - MU) / np.where(MU > 0, MU ** 2, np.nan)

        bt = np.load(ROOT / bt_p, allow_pickle=True)
        bg = np.array([str(g) for g in bt["gene"]])
        rho_all = bt["rho"].astype(np.float64)
        valid = bt["valid"].astype(bool)
        max_lag = int(bt["lags"].max())
        pos = {g: i for i, g in enumerate(genes_mo)}
        ii = np.array([pos[g] for g in bg])
        A = AL[ii]

        # 逐基因 rhobar
        RB = np.full((bg.size, len(KS)), np.nan)
        for i in range(bg.size):
            if not valid[i]:
                continue
            for j, k in enumerate(KS):
                RB[i, j], _ = rhobar(rho_all[i], k, max_lag)

        # 逐基因 gamma_pred（在 d=16..128 上拟合 rhobar 的衰减）
        okr = np.all(np.isfinite(RB) & (RB > 0), axis=1)
        x = np.log(np.array(DS_FIT, float))
        xm = x - x.mean()
        y = np.log(np.where(okr[:, None], RB, 1.0))
        ym = y - y.mean(axis=1, keepdims=True)
        g_pred = np.full(bg.size, np.nan)
        g_pred[okr] = -np.sum(xm[None, :] * ym, axis=1)[okr] / np.sum(xm ** 2)

        # 逐基因 gamma_obs——必须用同一区间 16-128，否则不可比
        jd = [int(np.where(ds == d)[0][0]) for d in DS_FIT]
        sub = A[:, jd]
        oko = np.all(np.isfinite(sub) & (sub > 0), axis=1)
        yo = np.log(np.where(oko[:, None], sub, 1.0))
        ymo = yo - yo.mean(axis=1, keepdims=True)
        g_obs = np.full(bg.size, np.nan)
        g_obs[oko] = -np.sum(xm[None, :] * ymo, axis=1)[oko] / np.sum(xm ** 2)

        # 同时给全窗 gamma（现稿口径 2-256）作参照
        okA = np.all(np.isfinite(A) & (A > 0), axis=1)
        xa = np.log(ds); xma = xa - xa.mean()
        ya = np.log(np.where(okA[:, None], A, 1.0))
        yma = ya - ya.mean(axis=1, keepdims=True)
        g_full = np.full(bg.size, np.nan)
        g_full[okA] = -np.sum(xma[None, :] * yma, axis=1)[okA] / np.sum(xma ** 2)

        m = okr & oko & valid
        go, gp = g_obs[m], g_pred[m]
        r_log = float(np.corrcoef(np.log(np.clip(go, 1e-6, None)),
                                  np.log(np.clip(gp, 1e-6, None)))[0, 1])
        rs = float(spearmanr(go, gp).statistic)
        slope = float(np.polyfit(gp, go, 1)[0])
        print(f"{tag:>6} {int(m.sum()):>6} {np.median(go):>7.3f} {np.median(gp):>7.3f} "
              f"{np.median(go - gp):>+7.3f} {r_log:>7.3f} {rs:>7.3f} {slope:>7.3f}")
        out[tag] = dict(n=int(m.sum()),
                        gamma_obs_median=float(np.median(go)),
                        gamma_pred_median=float(np.median(gp)),
                        diff_median=float(np.median(go - gp)),
                        diff_mean=float(np.mean(go - gp)),
                        diff_sd=float(np.std(go - gp)),
                        r_log=float(r_log), r_spearman=float(rs),
                        deming_like_slope_obs_on_pred=float(slope),
                        gamma_full_median=float(np.nanmedian(g_full[m])),
                        fit_range_um=DS_FIT, ks=KS,
                        rhobar_median={f"k{k}": float(np.nanmedian(RB[m][:, j]))
                                       for j, k in enumerate(KS)})

    # 合并三套
    allg = [(t, out[t]) for t in out]
    print("\n合并（三套加权前先看各自是否同号）:")
    print("  " + "  ".join(f"{t}: diff={v['diff_median']:+.3f} r={v['r_log']:.2f}"
                           for t, v in allg))
    (ROOT / "02_OUT/variogram/gamma_from_rhobar.json").write_text(
        json.dumps(dict(by_tissue=out,
                        note=("gamma_pred 由 rho(h) 的块平均 rhobar 在 d=16-128 上拟合；"
                              "gamma_obs 用同一区间以保可比；rho 的 8um 平滑使绝对值有偏，"
                              "故只比斜率与中位差")),
                   indent=2, ensure_ascii=False), encoding="utf-8")
    print("\n写出 02_OUT/variogram/gamma_from_rhobar.json")


if __name__ == "__main__":
    main()

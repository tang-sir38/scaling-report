r"""95_curve_shape.py — 用独立协方差预测整条 alpha 曲线（不只是斜率）

理论
----
对条件独立（Poisson）采样，
    Var(计数) - E(计数) = Var(Lambda_A)          （Poisson 项被精确约掉）
    Var(Lambda_A)       = |A|^2 sigma^2 rhobar(A) （rhobar = 块内平均相关函数）
    alpha(A) = (Var - E)/E^2 = (sigma^2/lambda^2) rhobar(A) = CV^2_0 * rhobar(A)
CV^2_0 是**点尺度**方差，与 A 无关。因此

    alpha(d) / alpha(d_ref) = rhobar(d) / rhobar(d_ref)

这是一条**整条曲线**的预测，不是一个斜率。斜率相等（gamma = beta）只是它的推论。

数据从哪来（两条独立路径）
--------------------------
  左式：聚合路径。gene_moments 的 (mu_d, var_d)。
  右式：空间路径。87 在 8um 网格上测的 rho_g(h)，再按 k x k 块做加权平均。
两条路径用的估计量完全不同，只有"同一个场"这一点是共同的。

量纲约定
--------
  网格间距 8um。d = 8k um，k = 1,2,4,8,16。
  rhobar(k) 由 rho(h) 各向同性近似算：
      rhobar(k) = sum_{a,b=0..k-1} m(a,b) rho(round(sqrt(a^2+b^2))) / k^4
      m = 1 (a=b=0), 2 (恰一个为 0), 4 (其余)
  k=1 时 rhobar=rho(0)=1，即把 8um 块当作点。
  这忽略了 8um 块内部的相关（网格分辨不出），
  表现为**一个与 d 无关的乘性常数**。因此：
      - 零参数版本：直接比 alpha(d)/alpha(8) 与 rhobar(k)
      - 单位置标定版本：在 d=16 处定常数，再预测 d=32/64/128
  只有后者能区分"常数偏移"与"形状错误"。

输出
----
02_OUT/variogram/curveshape_<tag>.json

用法
----
  python 95_curve_shape.py --tag colon
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

SETS = {
    "colon": ("02_OUT/gene_moments.npz", "02_OUT/variogram/beta_g_colon.npz"),
    "lung": ("D3_LungCancer/02_OUT/gene_moments_noMT.npz", "02_OUT/variogram/beta_g_lung.npz"),
    "brain": ("D2_MouseBrain/02_OUT/gene_moments_noMT.npz", "02_OUT/variogram/beta_g_brain.npz"),
}
KS = [1, 2, 4, 8, 16]


def rhobar_from_rho(rho_row, ks, max_lag):
    """rho_row: (max_lag+1,) 已归一化 rho(0)=1。返回 {k: rhobar(k)}。"""
    out = {}
    for k in ks:
        tot = 0.0
        cnt = 0
        for a0 in range(k):
            for b0 in range(k):
                m = 1 if (a0 == 0 and b0 == 0) else (2 if (a0 == 0 or b0 == 0) else 4)
                h = int(round(np.hypot(a0, b0)))
                if h > max_lag:
                    return out, False
                v = rho_row[h]
                if not np.isfinite(v):
                    return out, False
                tot += m * v
                cnt += m
        out[k] = tot / cnt if cnt else np.nan
    return out, True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="colon", choices=list(SETS))
    args = ap.parse_args()

    mo_p, bt_p = SETS[args.tag]
    mo = np.load(ROOT / mo_p, allow_pickle=True)
    ds = np.array([int(x) for x in mo["ds"]], float)
    genes_mo = np.array([str(g) for g in mo["genes"]])
    MU = np.column_stack([mo[f"mu_{int(d)}"] for d in ds])
    VA = np.column_stack([mo[f"var_{int(d)}"] for d in ds])
    with np.errstate(invalid="ignore", divide="ignore"):
        AL = (VA - MU) / np.where(MU > 0, MU ** 2, np.nan)

    bt = np.load(ROOT / bt_p, allow_pickle=True)
    bg = np.array([str(g) for g in bt["gene"]])
    rho = bt["rho"].astype(np.float64)
    lags = bt["lags"].astype(int)
    valid = bt["valid"].astype(bool)
    max_lag = int(lags.max())
    pos = {g: i for i, g in enumerate(genes_mo)}
    ii = np.array([pos[g] for g in bg])

    # 逐基因：alpha 在 d 网格上，rho 在 8um 网格上
    j8 = int(np.where(ds == 8)[0][0])
    jd = {8 * k: int(np.where(ds == 8 * k)[0][0]) for k in KS if 8 * k in ds}

    RB = np.full((bg.size, len(KS)), np.nan)
    okrow = np.zeros(bg.size, bool)
    for i in range(bg.size):
        rb, ok = rhobar_from_rho(rho[i], KS, max_lag)
        if ok and len(rb) == len(KS):
            RB[i] = [rb[k] for k in KS]
            okrow[i] = True

    A = AL[ii]
    okA = np.all(np.isfinite(A) & (A > 0), axis=1)
    use = valid & okrow & okA & (A[:, j8] > 0)
    print(f"[{args.tag}] 可用基因 {int(use.sum())} / {bg.size}  (rho 有效 {int(valid.sum())})")

    # ---------- 逐基因零参数对照 ----------
    dev = {}
    print(f"\n  逐基因  log[alpha(d)/alpha(8)]  -  log[rhobar(d/8)]")
    print(f"{'d':>6} {'k':>3} {'med_dev':>9} {'med|dev|':>9} {'q25':>8} {'q75':>8} {'r':>7}")
    for k in KS:
        d = 8 * k
        if d not in jd:
            continue
        obs = np.log(A[use, jd[d]] / A[use, j8])
        pred = np.log(RB[use, KS.index(k)])
        dv = obs - pred
        has_var = bool(np.std(pred) > 1e-12)
        dev[d] = dict(k=k, median=float(np.median(dv)),
                      median_abs=float(np.median(np.abs(dv))),
                      q25=float(np.percentile(dv, 25)), q75=float(np.percentile(dv, 75)),
                      pearson=float(np.corrcoef(obs, pred)[0, 1]) if has_var else None,
                      slope=float(np.polyfit(pred, obs, 1)[0]) if has_var else None,
                      rms=float(np.sqrt(np.mean(dv ** 2))))
        print(f"{d:>6} {k:>3} {dev[d]['median']:>+9.3f} {dev[d]['median_abs']:>9.3f} "
              f"{dev[d]['q25']:>+8.3f} {dev[d]['q75']:>+8.3f} "
              f"{(dev[d]['pearson'] if dev[d]['pearson'] is not None else float('nan')):>7.3f}")

    # ---------- 单位置标定后预测其余尺度 ----------
    print(f"\n  在 d=16 定标定常数，预测 d=32/64/128（每点都是样本外）")
    k16 = 2
    d16 = 16
    if d16 in jd:
        # 常数：c_g = log[alpha(16)/alpha(8)] - log rhobar(2)
        c = np.log(A[use, jd[d16]] / A[use, j8]) - np.log(RB[use, KS.index(k16)])
        cal = {}
        for k in KS:
            d = 8 * k
            if d not in jd or d <= d16:
                continue
            pred = np.log(RB[use, KS.index(k)]) + c
            obs = np.log(A[use, jd[d]] / A[use, j8])
            dv = obs - pred
            cal[d] = dict(k=k, median=float(np.median(dv)),
                          median_abs=float(np.median(np.abs(dv))),
                          rms=float(np.sqrt(np.mean(dv ** 2))),
                          pearson=float(np.corrcoef(obs, pred)[0, 1]),
                          slope=float(np.polyfit(pred, obs, 1)[0]))
            print(f"    d={d:>4}  med dev {cal[d]['median']:+.3f}  "
                  f"med|dev| {cal[d]['median_abs']:.3f}  r={cal[d]['pearson']:.3f}  "
                  f"slope={cal[d]['slope']:.3f}")

    # ---------- pooled（mu^2 加权，权重随 d 变，需逐 d 重算）----------
    def pooled(vals, cols):
        """vals: (G, len(cols)); 用各 d 自己的 mu^2 权重。"""
        out = []
        for j, d in zip(cols, ds[cols]):
            w = MU[ii][use, j] ** 2
            out.append(float(np.sum(w * vals[:, j]) / np.sum(w)))
        return np.array(out)

    cols = [jd[8 * k] for k in KS if 8 * k in jd]
    kk = [k for k in KS if 8 * k in jd]
    a_pool = pooled(A, cols)
    pred_pool = np.array([np.sum(MU[ii][use, j] ** 2 * RB[use, KS.index(k)])
                          for j, k in zip(cols, kk)]) / np.array(
        [np.sum(MU[ii][use, j] ** 2) for j in cols])
    print(f"\n  pooled: d, alpha_pool, 归一化观测, 归一化预测")
    obs_n = a_pool / a_pool[0]
    prd_n = pred_pool / pred_pool[0]
    pool_tab = []
    for k, d, o, p in zip(kk, [8 * k for k in kk], obs_n, prd_n):
        pool_tab.append(dict(d=int(d), alpha_pool=float(a_pool[kk.index(k)]),
                             obs_norm=float(o), pred_norm=float(p),
                             ratio=float(o / p)))
        print(f"    d={d:>4}  alpha={a_pool[kk.index(k)]:>10.4f}  obs={o:.4f}  "
              f"pred={p:.4f}  比值={o/p:.3f}")

    out = dict(tag=args.tag, n_genes=int(use.sum()), ks=kk,
               per_gene_zero_param=dev, per_gene_calibrated_at_16=cal,
               pooled=pool_tab,
               note=("rhobar 由 rho(h) 的 8um 环平均按各向同性近似做块平均；"
                     "8um 块内部相关不可分辨，故零参数版本含一个与 d 无关的常数偏移，"
                     "该偏移由 d=16 定标后消除"))
    (ROOT / "02_OUT/variogram" / f"curveshape_{args.tag}.json").write_text(
        json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n写出 02_OUT/variogram/curveshape_{args.tag}.json")


if __name__ == "__main__":
    main()

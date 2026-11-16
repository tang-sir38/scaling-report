r"""88_region_excise.py — T5：空间操控实验（因果证据）

目的
----
现稿已有的 "drop-IGKC" 只是把一个基因删掉（γ 从 0.555 降到 0.474），这是相关性的
削弱，不是因果。本脚本做的是**保持该基因总计数与强度边缘分布不变、只改其空间排布**，
因此任何 α/γ 的变化都只能归因于空间结构，不能归因于"基因身份"或"表达量"。

三种操作（都在 2 um 原生网格上做，再按同一套聚合器重算）
--------------------------------------------------------
  baseline   原样
  permute    把强度值在域内格子上随机置换 —— 边缘分布完全不变，空间相关被打散。
             理论预测：ρ 变成白噪声，rhobar(A) ∝ 1/A，故 α ∝ d^-2，γ -> 2。
  excise     把 IGKC 计数最高的 q 分位格（浆细胞聚集区）置零
  isolate    只保留这些高值格，其余置零

输出
----
02_OUT/excise/excise.json   每个操作、每个基因的 α(d) 与 γ
02_OUT/excise/pooled.json   IGKC 被改造后，按 μ² 权重重算的 pooled γ

为什么用权重重算 pooled 而不是重跑全矩阵
----------------------------------------
  α_pool(d) = Σ_g μ_g(d)^2 α_g(d) / Σ_g μ_g(d)^2
只改一个基因时，其余基因的 (μ, α) 不受影响（它们的场没动），所以
  α_pool_new(d) = [S_other(d) + μ_IGKC(d)^2 α_IGKC_new(d)] / [W_other + μ_IGKC(d)^2]
其中 S_other、W_other 从 gene_moments.npz 直接算。这在数学上与原式等价。

用法
----
  python 88_region_excise.py --tag colon --genes IGKC,IGHG1,JCHAIN,MZB1,EPCAM,COL1A1 \
      --nperm 20
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

ROOT = Path(__file__).resolve().parents[1]
KEY_STRIDE = 4096
DS = [2, 4, 8, 16, 32, 64, 128, 256]


def load_gene_fields(npz_path, gene_names):
    """返回 {gene: 2D float32 场}, 掩码 M, 网格信息。"""
    z = np.load(npz_path, allow_pickle=True)
    genes = np.array([str(g) for g in z["genes"]])
    shape = tuple(int(v) for v in z["shape"])
    coords = z["coords"].astype(np.int64)
    r0, c0 = int(coords[:, 0].min()), int(coords[:, 1].min())
    H = int(coords[:, 0].max()) - r0 + 1
    W = int(coords[:, 1].max()) - c0 + 1

    M = np.zeros((H, W), dtype=bool)
    keys = z["dom_any"].astype(np.int64) if "dom_any" in z.files else None
    if keys is not None:
        M[keys // KEY_STRIDE - r0, keys % KEY_STRIDE - c0] = True
    else:
        M[coords[:, 0] - r0, coords[:, 1] - c0] = True

    pos_r = coords[:, 0] - r0
    pos_c = coords[:, 1] - c0

    X = csr_matrix((z["data"], z["indices"], z["indptr"]), shape=shape)
    out = {}
    for g in gene_names:
        hit = np.where(genes == g)[0]
        if hit.size == 0:
            print(f"  [warn] {g} 不在矩阵中，跳过")
            continue
        row = X.getrow(int(hit[0]))
        F = np.zeros((H, W), dtype=np.float32)
        F[pos_r[row.indices], pos_c[row.indices]] = row.data
        out[g] = F
    return out, M, (r0, c0, H, W)


def alpha_curve_from_field(F, M, ds):
    """把 2um 场按 k x k 聚合（k = d/2），返回 {d: dict(mu,var,alpha)}。"""
    out = {}
    for d in ds:
        k = max(int(round(d / 2)), 1)
        H, W = F.shape
        nb_h, nb_w = H // k, W // k
        if nb_h < 4 or nb_w < 4:
            continue
        B = F[: nb_h * k, : nb_w * k].reshape(nb_h, k, nb_w, k).sum(axis=(1, 3))
        Mb = M[: nb_h * k, : nb_w * k].reshape(nb_h, k, nb_w, k).any(axis=(1, 3))
        v = B[Mb].astype(np.float64)
        if v.size < 16:
            continue
        mu = v.mean()
        var = v.var(ddof=1)
        out[d] = dict(mu=float(mu), var=float(var),
                      alpha=float((var - mu) / mu ** 2) if mu > 0 else np.nan,
                      n=int(v.size))
    return out


def fit_gamma(curve):
    """对 alpha(d) 做 log-log 拟合，返回 (gamma, se, r2)；只用 alpha>0。"""
    ds = np.array(sorted(curve))
    al = np.array([curve[d]["alpha"] for d in ds])
    m = np.isfinite(al) & (al > 0)
    if m.sum() < 3:
        return np.nan, np.nan, np.nan
    x, y = np.log(ds[m]), np.log(al[m])
    A = np.vstack([x, np.ones_like(x)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ coef
    r2 = 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    se = np.sqrt(((y - pred) ** 2).sum() / max(len(x) - 2, 1) /
                 ((x - x.mean()) ** 2).sum())
    return float(-coef[0]), float(se), float(r2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="colon")
    ap.add_argument("--bins", default="01_WORK/bins_002um.npz")
    ap.add_argument("--moments", default="02_OUT/gene_moments.npz")
    ap.add_argument("--genes", default="IGKC,IGHG1,JCHAIN,MZB1,EPCAM,COL1A1")
    ap.add_argument("--q", type=float, default=0.95, help="高值区阈值分位")
    ap.add_argument("--nperm", type=int, default=20)
    ap.add_argument("--outdir", default="02_OUT/excise")
    ap.add_argument("--seed", type=int, default=20261115)
    a = ap.parse_args()

    od = ROOT / a.outdir
    od.mkdir(parents=True, exist_ok=True)
    gene_list = [g.strip() for g in a.genes.split(",") if g.strip()]

    print(f"[{a.tag}] 读 2um 矩阵，取基因 {gene_list}")
    fields, M, geom = load_gene_fields(ROOT / a.bins, gene_list)
    H, W = geom[2], geom[3]
    print(f"  网格 {H}x{W}，域内 {int(M.sum()):,} 格")

    rng = np.random.default_rng(a.seed)
    res = {"tag": a.tag, "geom": dict(H=H, W=W, domain=int(M.sum())),
           "q": a.q, "nperm": a.nperm, "genes": {}}

    for g, F in fields.items():
        print(f"\n=== {g} ===")
        tot = float(F.sum())
        occ = float((F[M] > 0).mean())
        base = alpha_curve_from_field(F, M, DS)
        g_base, se_b, r2_b = fit_gamma(base)
        print(f"  baseline  total={tot:.3e}  occupancy={occ:.3f}  "
              f"γ={g_base:.3f}±{se_b:.3f} (r2={r2_b:.4f})  α(2um)={base.get(2,{}).get('alpha',np.nan):.4f}")

        # --- permute：边缘分布不变，只打乱空间位置 ---
        idx = np.flatnonzero(M)
        vals = F.ravel()[idx]
        g_perm, a2_perm = [], []
        for it in range(a.nperm):
            p = rng.permutation(vals.size)
            Fp = np.zeros_like(F)
            Fp.ravel()[idx] = vals[p]
            c = alpha_curve_from_field(Fp, M, DS)
            gg, _, r2p = fit_gamma(c)
            if np.isfinite(gg):
                g_perm.append(gg)
                a2_perm.append(c.get(2, {}).get("alpha", np.nan))
        g_perm = np.array(g_perm)
        print(f"  permute   γ={np.nanmedian(g_perm):.3f} "
              f"[{np.nanpercentile(g_perm,2.5):.3f}, {np.nanpercentile(g_perm,97.5):.3f}] "
              f"(n={g_perm.size})  α(2um)={np.nanmedian(a2_perm):.4f}")

        # --- excise / isolate ---
        vals_in = F[M]
        thr = np.quantile(vals_in[vals_in > 0], a.q) if (vals_in > 0).any() else 0
        hot = M & (F >= thr) & (F > 0)
        print(f"  高值区阈值={thr:.0f}  格数={int(hot.sum()):,} "
              f"({hot.sum()/M.sum():.2%})  占总计数 {F[hot].sum()/max(tot,1):.1%}")

        Fe = F.copy(); Fe[hot] = 0
        Fi = np.zeros_like(F); Fi[hot] = F[hot]
        ce = alpha_curve_from_field(Fe, M, DS)
        ci = alpha_curve_from_field(Fi, M, DS)
        g_ex, se_e, r2_e = fit_gamma(ce)
        g_is, se_i, r2_i = fit_gamma(ci)
        print(f"  excise    γ={g_ex:.3f} (r2={r2_e:.3f})  α(2um)={ce.get(2,{}).get('alpha',np.nan):.4f}")
        print(f"  isolate   γ={g_is:.3f} (r2={r2_i:.3f})  α(2um)={ci.get(2,{}).get('alpha',np.nan):.4f}")

        res["genes"][g] = dict(
            total=tot, occupancy=occ, hot_threshold=float(thr),
            hot_cells=int(hot.sum()), hot_share_of_counts=float(F[hot].sum() / max(tot, 1)),
            baseline=dict(gamma=g_base, se=se_b, r2=r2_b, curve=base),
            permute=dict(gamma_median=float(np.nanmedian(g_perm)),
                         gamma_lo=float(np.nanpercentile(g_perm, 2.5)),
                         gamma_hi=float(np.nanpercentile(g_perm, 97.5)),
                         alpha2_median=float(np.nanmedian(a2_perm)),
                         gammas=g_perm.tolist()),
            excise=dict(gamma=g_ex, se=se_e, r2=r2_e, curve=ce),
            isolate=dict(gamma=g_is, se=se_i, r2=r2_i, curve=ci),
        )

    # ---------- pooled 重算（只对第一个基因，默认 IGKC）----------
    z = np.load(ROOT / a.moments, allow_pickle=True)
    ds_all = [int(x) for x in z["ds"]]
    genes_mo = np.array([str(g) for g in z["genes"]])
    target = gene_list[0]
    ti = np.where(genes_mo == target)[0]
    pooled = {}
    if ti.size and target in res["genes"]:
        ti = int(ti[0])
        MU = np.column_stack([z[f"mu_{d}"] for d in ds_all])
        VA = np.column_stack([z[f"var_{d}"] for d in ds_all])
        with np.errstate(invalid="ignore", divide="ignore"):
            AL = (VA - MU) / np.where(MU > 0, MU ** 2, np.nan)
        W2 = MU ** 2
        S_all = np.nansum(W2 * AL, axis=0)
        W_all = np.nansum(W2, axis=0)
        S_other = S_all - W2[ti] * AL[ti]
        W_other = W_all - W2[ti]

        def pooled_gamma(alpha_new):
            num = S_other + W2[ti] * np.array(alpha_new)
            with np.errstate(invalid="ignore", divide="ignore"):
                ap = num / W_other
            ok = np.isfinite(ap) & (ap > 0)
            if ok.sum() < 3:
                return np.nan, np.nan, ap
            x, y = np.log(np.array(ds_all)[ok]), np.log(ap[ok])
            A = np.vstack([x, np.ones_like(x)]).T
            coef, *_ = np.linalg.lstsq(A, y, rcond=None)
            r2 = 1 - ((y - A @ coef) ** 2).sum() / ((y - y.mean()) ** 2).sum()
            return float(-coef[0]), float(r2), ap

        def curve_to_alpha(c):
            return [c[d]["alpha"] if d in c else np.nan for d in ds_all]

        g0, r20, ap0 = pooled_gamma([AL[ti][j] for j in range(len(ds_all))])
        cases = {}
        for op in ("baseline", "excise", "isolate"):
            c = res["genes"][target][op]["curve"]
            gg, rr, ap = pooled_gamma(curve_to_alpha(c))
            cases[op] = dict(gamma=gg, r2=rr, alpha_pool=ap.tolist())
        # permute 用中位 α(d)（需要逐次留存曲线，这里重跑少量）
        gp, rp, app = [], [], []
        F = fields[target]
        idx = np.flatnonzero(M); vals = F.ravel()[idx]
        for it in range(min(a.nperm, 8)):
            p = rng.permutation(vals.size)
            Fp = np.zeros_like(F); Fp.ravel()[idx] = vals[p]
            c = alpha_curve_from_field(Fp, M, DS)
            gg, rr, ap = pooled_gamma(curve_to_alpha(c))
            gp.append(gg); rp.append(rr); app.append(ap)
        cases["permute"] = dict(gamma=float(np.nanmedian(gp)), r2=float(np.nanmedian(rp)),
                                alpha_pool=np.nanmedian(np.array(app), axis=0).tolist())
        pooled = dict(target=target, ds=ds_all,
                      gamma_baseline=g0, r2_baseline=r20,
                      cases=cases,
                      delta_permute=cases["permute"]["gamma"] - g0,
                      delta_excise=cases["excise"]["gamma"] - g0,
                      delta_isolate=cases["isolate"]["gamma"] - g0)
        print("\n" + "=" * 72)
        print(f"pooled γ（只改造 {target}）")
        for k, v in cases.items():
            print(f"  {k:<10} γ={v['gamma']:.4f} (r2={v['r2']:.4f})  Δ vs baseline = "
                  f"{v['gamma']-g0:+.4f}")
        print("=" * 72)

    res["pooled"] = pooled
    (od / f"excise_{a.tag}.json").write_text(
        json.dumps(res, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n写出 {od}/excise_{a.tag}.json")


if __name__ == "__main__":
    main()

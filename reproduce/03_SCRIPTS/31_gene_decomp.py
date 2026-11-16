r"""
31_gene_decomp.py — 逐基因分解：pooled 指数是不是少数基因 + 混合的产物？

三个问题
--------
Q1 构成：pooled alpha 的 mu^2 权重有多少落在 top-100 / top-500 基因上？有效基因数 1/sum(w^2)？
Q2 异质性：逐基因的 gamma_g 分布（中位数、四分位）；未加权中位 gamma 与 pooled gamma 差多少？
Q3 单基因形状：每个基因单独拟合 "纯幂律" vs "有限相关长度 H(d/xi)"，谁的 AIC 更低？
   若多数基因偏好有限 xi，而 pooled 偏好幂律，则 pooled 的幂律是**多尺度混合的假象**，
   这正好解释 "gamma 随表达量变化"（不同 xi_g 的混合权重随表达层变化）。

方法
----
逐基因 alpha_g(d) = (var_g(d) - mu_g(d)) / mu_g(d)^2。
筛选：8 个尺度上 mu_g>0 且 alpha_g>0，且 mu_g(256um) >= 0.5（保证不是纯噪声）。
xi 拟合用 profile（固定 xi 网格，A 解析最优），避免 19k 次非线性优化。
"""
import argparse, json
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline

ROOT = Path(__file__).resolve().parents[1]

def build_H(n=200, wmax=400.0, nw=180):
    x, wq = np.polynomial.legendre.leggauss(n)
    s, t = x[:, None], x[None, :]
    wt = wq[:, None] * wq[None, :]
    taper = (1 - np.abs(s)) * (1 - np.abs(t))
    r = np.sqrt(s ** 2 + t ** 2)
    wg = np.concatenate([[0.0], np.geomspace(1e-3, wmax, nw)])
    H = np.array([float(np.sum(wt * taper * np.exp(-w * r))) for w in wg])
    H /= H[0]
    return CubicSpline(np.log(wg[1:]), np.log(H[1:]), extrapolate=True)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=str(ROOT / "02_OUT" / "gene_decomp"))
    ap.add_argument("--min-mu256", type=float, default=0.5)
    ap.add_argument("--nrand", type=int, default=200)
    ap.add_argument("--seed", type=int, default=11)
    a = ap.parse_args()
    od = Path(a.outdir); od.mkdir(parents=True, exist_ok=True)

    sets = {
        "colon": str(ROOT / "02_OUT" / "gene_moments.npz"),
        "lung": str(ROOT / "D3_LungCancer" / "02_OUT" / "gene_moments_noMT.npz"),
        "brain": str(ROOT / "D2_MouseBrain" / "02_OUT" / "gene_moments_noMT.npz"),
    }
    spl = build_H()
    xg = np.geomspace(0.5, 1e6, 120)          # xi 候选网格（um）
    lx = np.log(xg)
    rng = np.random.default_rng(a.seed)
    out = {}

    for tag, p in sets.items():
        if not Path(p).exists():
            print("缺", p); continue
        z = np.load(p, allow_pickle=True)
        ds = np.array([int(x) for x in z["ds"]], float)
        genes = np.array([str(g) for g in z["genes"]])
        MU = np.column_stack([z[f"mu_{int(d)}"] for d in ds])       # (G, 8)
        VA = np.column_stack([z[f"var_{int(d)}"] for d in ds])
        AL = (VA - MU) / np.where(MU > 0, MU ** 2, np.nan)

        # ---------- Q1 构成 ----------
        mu32 = MU[:, int(np.where(ds == 32)[0][0])]
        w = mu32 ** 2
        w = w / w.sum()
        order = np.argsort(w)[::-1]
        top100 = float(w[order[:100]].sum())
        top500 = float(w[order[:500]].sum())
        neff = float(1.0 / np.sum(w ** 2))
        # 累计到 50% / 90% 需要多少基因
        cum = np.cumsum(w[order])
        n50 = int(np.searchsorted(cum, 0.5) + 1)
        n90 = int(np.searchsorted(cum, 0.9) + 1)

        # ---------- Q3 单基因拟合 ----------
        ok = np.all(np.isfinite(AL) & (AL > 0), axis=1) & (MU[:, -1] >= a.min_mu256)
        idx = np.where(ok)[0]
        n_ok = idx.size
        gamma_g = np.full(genes.size, np.nan)
        xi_g = np.full(genes.size, np.nan)
        daic = np.full(genes.size, np.nan)
        pref_xi = np.zeros(genes.size, bool)

        if n_ok:
            y = np.log(AL[idx])                       # (n,8)
            xd = np.log(ds)
            # 纯幂律：解析最小二乘
            xm = xd - xd.mean()
            b = (y * xm).sum(1) / (xm ** 2).sum()
            a0 = y.mean(1) - b * xd.mean()
            rss_p = ((y - (a0[:, None] + b[:, None] * xd[None, :])) ** 2).sum(1)
            gamma_g[idx] = -b
            # 有限 xi：固定 xi 网格，解析 A
            lh = spl(np.log(ds[None, :] / xg[:, None]))          # (nxi,8)
            d2 = y[:, None, :] - lh[None, :, :]                   # (n_gene,nxi,8)
            rss_x = ((d2 - d2.mean(2, keepdims=True)) ** 2).sum(2)
            jbest = np.argmin(rss_x, axis=1)
            rss_xi = rss_x[np.arange(idx.size), jbest]
            xi_g[idx] = xg[jbest]
            n = ds.size
            aic_p = n * np.log(rss_p / n) + 2 * 2
            aic_x = n * np.log(np.maximum(rss_xi, 1e-12) / n) + 2 * 2
            daic[idx] = aic_p - aic_x        # >0 表示 xi 模型更好
            pref_xi[idx] = daic[idx] > 2

        # ---------- Q2 gamma 分布 ----------
        g_ok = gamma_g[idx]
        med = float(np.median(g_ok)) if n_ok else None
        q25, q75 = (float(np.percentile(g_ok, 25)), float(np.percentile(g_ok, 75))) if n_ok else (None, None)
        # pooled gamma（与主分析同口径）
        num = np.nansum(VA - MU, axis=0)
        den = np.nansum(MU ** 2, axis=0)
        al_pool = num / den
        pooled_gamma = float(-np.polyfit(np.log(ds), np.log(al_pool), 1)[0])

        # 随机 100 基因子集的 pooled gamma 分布
        rand_g = []
        if n_ok > 200:
            for _ in range(a.nrand):
                sel = rng.choice(idx, size=100, replace=False)
                num_r = np.nansum(VA[sel] - MU[sel], axis=0)
                den_r = np.nansum(MU[sel] ** 2, axis=0)
                if np.all(den_r > 0) and np.all(num_r > 0):
                    rand_g.append(float(-np.polyfit(np.log(ds), np.log(num_r / den_r), 1)[0]))

        # xi 与表达量的关系
        corr = None
        if n_ok > 50:
            from scipy.stats import spearmanr
            corr = float(spearmanr(np.log10(MU[idx, -1] + 1e-12), np.log10(xi_g[idx]))[0])

        rec = dict(
            n_genes=int(genes.size), n_fit=int(n_ok),
            frac_fit=float(n_ok / genes.size),
            composition=dict(top100_share=top100, top500_share=top500,
                             n_eff=neff, n_for_50pct=n50, n_for_90pct=n90),
            gamma_gene=dict(median=med, q25=q25, q75=q75,
                            iqr=float(q75 - q25) if q75 is not None else None,
                            pooled=pooled_gamma,
                            ratio_median_to_pooled=float(med / pooled_gamma) if med else None,
                            p05=float(np.percentile(g_ok, 5)) if n_ok else None,
                            p95=float(np.percentile(g_ok, 95)) if n_ok else None,
                            gmax=float(np.max(g_ok)) if n_ok else None,
                            frac_gt_1=float(np.mean(g_ok > 1.0)) if n_ok else None,
                            frac_gt_1p5=float(np.mean(g_ok > 1.5)) if n_ok else None),
            random100=dict(n=int(len(rand_g)),
                           med=float(np.median(rand_g)) if rand_g else None,
                           lo=float(np.percentile(rand_g, 2.5)) if rand_g else None,
                           hi=float(np.percentile(rand_g, 97.5)) if rand_g else None,
                           sd=float(np.std(rand_g)) if rand_g else None),
            xi_fit=dict(n_pref_xi=int(pref_xi.sum()),
                        frac_pref_xi=float(pref_xi.sum() / n_ok) if n_ok else None,
                        xi_median=float(np.median(xi_g[idx])) if n_ok else None,
                        xi_q25=float(np.percentile(xi_g[idx], 25)) if n_ok else None,
                        xi_q75=float(np.percentile(xi_g[idx], 75)) if n_ok else None,
                        xi_min=float(np.min(xi_g[idx])) if n_ok else None,
                        xi_max=float(np.max(xi_g[idx])) if n_ok else None,
                        frac_at_bound=float(np.mean((xi_g[idx] <= xg[1]) | (xi_g[idx] >= xg[-2]))) if n_ok else None,
                        spearman_xi_vs_mu=corr),
        )
        out[tag] = rec
        print(f"\n== {tag} ==")
        print(f"  可拟合基因 {n_ok}/{genes.size} ({100*n_ok/genes.size:.1f}%), 门槛 mu256>={a.min_mu256}")
        print(f"  top100 占 mu^2 权重 {100*top100:.1f}% ; top500 {100*top500:.1f}% ; "
              f"有效基因数 {neff:.1f} ; 50% 需 {n50} 个基因, 90% 需 {n90} 个")
        print(f"  逐基因 gamma 中位 {med:.3f}  IQR [{q25:.3f},{q75:.3f}]  vs pooled {pooled_gamma:.3f}"
              if med else "  gamma 不可估")
        gg = rec["gamma_gene"]
        print(f"    p05 {gg['p05']:.3f}  p95 {gg['p95']:.3f}  max {gg['gmax']:.3f}  "
              f">1.0 占 {100*gg['frac_gt_1']:.1f}%  >1.5 占 {100*gg['frac_gt_1p5']:.2f}%")
        if rand_g:
            print(f"  随机100基因 pooled gamma: {np.median(rand_g):.3f} "
                  f"[{np.percentile(rand_g,2.5):.3f},{np.percentile(rand_g,97.5):.3f}] sd={np.std(rand_g):.3f}")
        xf = rec["xi_fit"]
        print(f"  单基因更偏好有限 xi 的比例 {100*xf['frac_pref_xi']:.1f}%; "
              f"xi 中位 {xf['xi_median']:.1f} um, IQR [{xf['xi_q25']:.1f},{xf['xi_q75']:.1f}]")
        print(f"  xi 落在网格边界的比例 {100*xf['frac_at_bound']:.1f}% ; xi vs mu Spearman {xf['spearman_xi_vs_mu']}")

    json.dump(out, open(od / "gene_decomp.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(f"\n写出 {od/'gene_decomp.json'}")

if __name__ == "__main__":
    main()

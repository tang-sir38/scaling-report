r"""
36_jackknife.py — pooled 指数有多稳？留一基因刀切 + 基因累积曲线 + 分窗一致性。

问题
----
pooled alpha(d) = SUM_g (var_g - mu_g) / SUM_g mu_g^2 是逐基因 alpha 的 mu^2 加权平均。
若权重集中在极少数基因，pooled gamma 就只是那几个基因的 gamma。

三项检查
--------
J1 留一基因（LOGO）：每次剔除 1 个基因，重算 pooled alpha(8 档) 与 gamma。
   报 gamma 的刀切分布、最大位移、以及贡献最大的基因。
J2 累积曲线：按 mu^2 从大到小逐个加入基因，看 gamma 如何被"搭建"出来；
   另做随机子集（k = 1..all，每个 k 抽 nrand 次），给随机基因集的分布。
J3 分窗一致性：每个基因在细窗（2-32 um）与粗窗（32-256 um）各自的 gamma，
   若两者高度相关，说明逐基因 gamma 是基因属性而非拟合噪声。

输出 <outdir>/jackknife.json
"""
import argparse, json
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]

SETS = {
    "colon": str(ROOT / "02_OUT" / "gene_moments.npz"),
    "lung": str(ROOT / "D3_LungCancer" / "02_OUT" / "gene_moments_noMT.npz"),
    "brain": str(ROOT / "D2_MouseBrain" / "02_OUT" / "gene_moments_noMT.npz"),
}


def fit_gamma(ds, alpha):
    """log alpha = c - gamma log d 的最小二乘解。alpha 可为 (G,8) 或 (8,)。"""
    x = np.log(ds)
    y = np.log(alpha)
    xm = x - x.mean(axis=-1, keepdims=True) if y.ndim > 1 else x - x.mean()
    if y.ndim == 1:
        return -float(np.sum(xm * (y - y.mean())) / np.sum(xm ** 2))
    ym = y - y.mean(axis=-1, keepdims=True)
    return -np.sum(xm * ym, axis=-1) / np.sum(xm ** 2, axis=-1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=str(ROOT / "02_OUT" / "jackknife"))
    ap.add_argument("--min-mu256", type=float, default=0.5)
    ap.add_argument("--nrand", type=int, default=200)
    ap.add_argument("--seed", type=int, default=23)
    a = ap.parse_args()
    od = Path(a.outdir); od.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    out = {}

    for tag, p in SETS.items():
        if not Path(p).exists():
            print("缺", p); continue
        z = np.load(p, allow_pickle=True)
        ds = np.array([int(x) for x in z["ds"]], float)
        genes = np.array([str(g) for g in z["genes"]])
        MU = np.column_stack([z[f"mu_{int(d)}"] for d in ds])        # (G,8)
        VA = np.column_stack([z[f"var_{int(d)}"] for d in ds])
        EX = VA - MU                                                  # 过散量（可加）
        MU2 = MU ** 2

        S1 = EX.sum(0)                       # (8,)
        S2 = MU2.sum(0)
        g_full = fit_gamma(ds, S1 / S2)

        # ---------- J1 留一基因 ----------
        # 分母逐基因可加：alpha_{-g} = (S1 - ex_g) / (S2 - mu2_g)
        keep = S2[None, :] - MU2 > 0
        num = S1[None, :] - EX
        den = S2[None, :] - MU2
        al = np.where(keep, num / np.where(keep, den, 1.0), np.nan)
        okrow = np.all(np.isfinite(al) & (al > 0), axis=1)
        g_loo = np.full(genes.size, np.nan)
        g_loo[okrow] = fit_gamma(ds, al[okrow])
        d_loo = g_loo - g_full
        j = int(np.nanargmax(np.abs(d_loo)))
        # 权重最大的基因（32 um 口径）
        i32 = int(np.where(ds == 32)[0][0])
        w32 = MU2[:, i32] / MU2[:, i32].sum()
        order = np.argsort(w32)[::-1]
        top1 = int(order[0])

        # ---------- J2a 累积（按 mu^2 降序）----------
        ks = [1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, genes.size]
        ks = [k for k in ks if k <= genes.size]
        cum_gamma = []
        for k in ks:
            sel = order[:k]
            s1, s2 = EX[sel].sum(0), MU2[sel].sum(0)
            cum_gamma.append(float(fit_gamma(ds, s1 / s2)) if np.all(s2 > 0) else None)

        # ---------- J2b 随机子集 ----------
        rand = {}
        for k in (1, 2, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, genes.size):
            if k > genes.size:
                continue
            gs = []
            for _ in range(a.nrand if k < genes.size else 1):
                sel = rng.choice(genes.size, size=k, replace=False)
                s1, s2 = EX[sel].sum(0), MU2[sel].sum(0)
                if np.all(s2 > 0):
                    gs.append(float(fit_gamma(ds, s1 / s2)))
            if gs:
                rand[str(k)] = dict(n=len(gs), med=float(np.median(gs)),
                                    lo=float(np.percentile(gs, 5)),
                                    hi=float(np.percentile(gs, 95)))

        # ---------- J3 分窗一致性 ----------
        AL = EX / np.where(MU > 0, MU2, np.nan)
        okg = np.all(np.isfinite(AL) & (AL > 0), axis=1) & (MU[:, -1] >= a.min_mu256)
        gf = gc = None
        rho = None
        if okg.sum() > 50:
            fine = ds <= 32
            coarse = ds >= 32
            gf = fit_gamma(ds[fine], AL[okg][:, fine])
            gc = fit_gamma(ds[coarse], AL[okg][:, coarse])
            rho = float(spearmanr(gf, gc).statistic)

        # ---------- J4 逐基因 gamma 与表达量的关系 ----------
        gam = np.full(genes.size, np.nan)
        gam[okg] = fit_gamma(ds, AL[okg])
        rho_mu = float(spearmanr(gam[okg], MU[okg, i32]).statistic)
        # 按 mu32 数量级分箱：权重占比 + 中位 gamma
        wdec = []
        for lo, hi in ((0, .01), (.01, .1), (.1, 1), (1, 10), (10, 100), (100, 1e12)):
            m = (MU[:, i32] >= lo) & (MU[:, i32] < hi)
            if m.sum() == 0:
                continue
            wdec.append(dict(lo=lo, hi=hi, n=int(m.sum()),
                             w=float(w32[m].sum()),
                             med_gamma=float(np.nanmedian(gam[m])) if np.isfinite(gam[m]).any() else None))

        out[tag] = dict(
            gamma_full=float(g_full),
            n_genes=int(genes.size),
            top1_gene=str(genes[top1]), top1_weight=float(w32[top1]),
            loo=dict(n=int(okrow.sum()),
                     med=float(np.nanmedian(g_loo)),
                     p05=float(np.nanpercentile(g_loo, 5)),
                     p95=float(np.nanpercentile(g_loo, 95)),
                     min=float(np.nanmin(g_loo)), max=float(np.nanmax(g_loo)),
                     max_abs_delta=float(np.nanmax(np.abs(d_loo))),
                     worst_gene=str(genes[j]), worst_delta=float(d_loo[j]),
                     n_delta_gt_0p05=int(np.nansum(np.abs(d_loo) > 0.05)),
                     n_delta_gt_0p02=int(np.nansum(np.abs(d_loo) > 0.02))),
            cumulative=dict(k=ks, gamma=cum_gamma),
            random=rand,
            split=dict(n=int(okg.sum()), spearman=rho,
                       gamma_fine_med=float(np.median(gf)) if gf is not None else None,
                       gamma_coarse_med=float(np.median(gc)) if gc is not None else None,
                       mad=float(np.median(np.abs(gf - gc))) if gf is not None else None),
            gamma_vs_mu=dict(spearman=rho_mu, by_decade=wdec),
            source="36_jackknife.py")

        r = out[tag]
        print(f"\n== {tag} ==  gamma_full={g_full:.4f}")
        print(f"  spearman(gamma,mu32)={rho_mu:+.3f}")
        print("  按表达量数量级的权重与中位 gamma:")
        for x in wdec:
            print(f"    [{x['lo']:>6g},{x['hi']:<8g}) n={x['n']:5d} w={x['w']*100:6.2f}%"
                  f" med_gamma={x['med_gamma']}")
        print(f"  LOGO: med={r['loo']['med']:.4f} range=[{r['loo']['min']:.3f},{r['loo']['max']:.3f}]"
              f"  max|d|={r['loo']['max_abs_delta']:.4f} by {r['loo']['worst_gene']}"
              f"  n(|d|>0.05)={r['loo']['n_delta_gt_0p05']}")
        print(f"  top1 by mu2: {r['top1_gene']} {r['top1_weight']*100:.1f}%")
        print(f"  累积 gamma: " + " ".join(f"{k}:{v:.3f}" for k, v in zip(ks, cum_gamma) if v))
        print(f"  随机 {list(rand)[:5]} ... 全基因 {rand.get(str(genes.size),{})}")
        print(f"  分窗 spearman={rho:.3f} med_fine={r['split']['gamma_fine_med']:.3f}"
              f" med_coarse={r['split']['gamma_coarse_med']:.3f} MAD={r['split']['mad']:.3f}")

    json.dump(out, open(od / "jackknife.json", "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print(f"\n写出 {od/'jackknife.json'}")


if __name__ == "__main__":
    main()

r"""
07_pooled_alpha.py — 池化过散度：不依赖基因筛选的标度律估计。

为什么需要它
------------
逐基因的 alpha 要求"该基因在每个分辨率下都有足够计数"。Visium HD 的单基因计数
在 8-16 um 档只有 0.01-0.1 UMI/bin，能过阈值的基因只剩几十个，且都是极端高表达
基因（线粒体、核糖体），斜率由这个小而偏的基因子集决定，不可靠，也无法跨数据集比较。

池化估计
--------
矩估计：  E[ Var_g ] = mu_g + alpha_g * mu_g^2
对基因求和，alpha 用 mu^2 加权：

    alpha_pool(d) = SUM_g ( Var_g(d) - mu_g(d) ) / SUM_g mu_g(d)^2

分子是"总过散量"，分母是"mu^2 总量"。这是 alpha 的 mu^2 加权平均，
无需任何基因筛选，所有基因都参与，且对每个基因的估计噪声无偏
（Var_g 用 n-1 分母时 E[Var_g]=mu_g+alpha_g mu_g^2 严格成立）。

物理解释
--------
alpha ~ 1 / N_eff，N_eff 是一个 bin 内的"有效独立抽样单位数"。
  alpha ∝ d^-2  <=>  N_eff ∝ 面积（细胞各自独立）
  alpha ∝ d^0   <=>  N_eff 不随面积增长（表达在空间上完全相关）

输出
----
  <outdir>/pooled_alpha.json
  <outdir>/pooled_alpha_bybin.json   按表达量分层的版本（检查是否依赖表达水平）
"""
import argparse
import json
from pathlib import Path

import numpy as np


def load_bins(p):
    z = np.load(p, allow_pickle=False)
    shape = tuple(int(x) for x in z["shape"])
    from scipy.sparse import csr_matrix
    X = csr_matrix((z["data"], z["indices"], z["indptr"]), shape=shape)
    return X, z["coords"], z["genes"], z


def per_gene_moments(X, n_cells):
    """返回每基因的 s1、s2、mu、var。n_cells 是口径内的格子总数（含零格）。"""
    s1 = np.asarray(X.sum(axis=1)).ravel().astype(np.float64)
    s2 = np.asarray(X.multiply(X).sum(axis=1)).ravel().astype(np.float64)
    mu = s1 / n_cells
    var = (s2 - s1 ** 2 / n_cells) / (n_cells - 1)
    return s1, s2, mu, var


def pooled_alpha(mu, var, mask=None):
    """alpha_pool = SUM(var - mu) / SUM(mu^2)。"""
    if mask is not None:
        mu, var = mu[mask], var[mask]
    num = float(np.sum(var - mu))
    den = float(np.sum(mu ** 2))
    if den <= 0:
        return np.nan, num, den
    return num / den, num, den


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--ds", nargs="+", type=int, default=[2, 4, 8, 16, 32, 64])
    ap.add_argument("--domain", default="any", choices=["any", "full"])
    ap.add_argument("--nbins", type=int, default=6, help="按表达量分层数")
    ap.add_argument("--exclude-mt", action="store_true",
                    help="剔除线粒体基因（人 MT-，鼠 mt-）；它们是技术噪声的常见来源")
    a = ap.parse_args()

    wd, od = Path(a.workdir), Path(a.outdir)
    od.mkdir(parents=True, exist_ok=True)

    rows = {}
    MOM = {}        # d -> (mu, var)，缓存备用
    mu_ref = None   # 用最大分辨率的 mu 做分层依据
    keep = None     # 非线粒体基因掩码（--exclude-mt 时启用）
    for d in a.ds:
        p = wd / f"bins_{d:03d}um.npz"
        if not p.exists():
            print(f"  缺 {p.name}，跳过")
            continue
        X, coords, genes, z = load_bins(p)
        key = f"n_domain_{a.domain}"
        n_cells = int(z[key]) if key in z.files else int(X.shape[1])
        s1, s2, mu, var = per_gene_moments(X, n_cells)

        if a.exclude_mt:
            if keep is None:
                gs = [g.decode() if isinstance(g, bytes) else str(g) for g in genes]
                keep = ~np.array([s.upper().startswith("MT-") for s in gs])
                print(f"  剔除线粒体基因 {int((~keep).sum())} 个，保留 {int(keep.sum())} 个",
                      flush=True)
            mu, var, s1 = mu[keep], var[keep], s1[keep]

        MOM[d] = (mu, var)

        ap_val, num, den = pooled_alpha(mu, var)
        rows[d] = dict(n_cells=n_cells, n_obs=int(X.shape[1]),
                       alpha_pool=ap_val, num=num, den=den,
                       total_umi=float(s1.sum()),
                       mu_mean=float(mu.mean()), mu_median=float(np.median(mu)))
        print(f"  d={d:3d}um  bins={n_cells:>9,}  alpha_pool={ap_val:.5f}  "
              f"总UMI={s1.sum():.3e}", flush=True)
        if d == max(a.ds):
            mu_ref = mu.copy()
        del X, s1, s2

    ds = sorted(rows)
    if len(ds) < 3:
        print("分辨率不足 3 档，无法拟合")
        return

    # ---- 主拟合：log alpha_pool ~ log d ----
    y = np.array([np.log(rows[d]["alpha_pool"]) for d in ds])
    x = np.log(np.array(ds, float))
    slope = float(np.polyfit(x, y, 1)[0])
    # 逐段比值
    seg = []
    for i in range(len(ds) - 1):
        r = rows[ds[i + 1]]["alpha_pool"] / rows[ds[i]]["alpha_pool"]
        seg.append(dict(frm=ds[i], to=ds[i + 1], ratio=float(r),
                        slope=float(np.log(r) / np.log(ds[i + 1] / ds[i]))))
        print(f"  {ds[i]:>3}->{ds[i+1]:>3}um  比值 {r:.3f}  斜率 {seg[-1]['slope']:+.3f}")

    print(f"\n  == 池化斜率（全部 {len(ds)} 档 d={ds[0]}-{ds[-1]}）= {slope:+.3f} ==")

    # ---- 只用 d>=8 的子区间（与逐基因口径可比）----
    sub = [d for d in ds if d >= 8]
    slope_sub = None
    if len(sub) >= 3:
        ys = np.array([np.log(rows[d]["alpha_pool"]) for d in sub])
        xs = np.log(np.array(sub, float))
        slope_sub = float(np.polyfit(xs, ys, 1)[0])
        print(f"  == 池化斜率（d=8-{sub[-1]}，与逐基因口径同区间）= {slope_sub:+.3f} ==")

    # ---- 按表达量分层，检查斜率是否依赖表达水平 ----
    bybin = {}
    if mu_ref is not None:
        edges = np.quantile(mu_ref[mu_ref > 0], np.linspace(0, 1, a.nbins + 1))
        edges = np.unique(edges)
        print(f"\n  按 mu({max(a.ds)}um) 分 {len(edges)-1} 层：")
        for k in range(len(edges) - 1):
            lo, hi = edges[k], edges[k + 1]
            mask = (mu_ref >= lo) & (mu_ref <= hi if k == len(edges) - 2 else mu_ref < hi)
            if mask.sum() < 20:
                continue
            vals = []
            for d in ds:
                mu, var = MOM[d]
                av, _, _ = pooled_alpha(mu, var, mask)
                vals.append(av)
            lab = f"{lo:.3g}-{hi:.3g}"
            if np.all(np.isfinite(vals)) and min(vals) > 0:
                sl = float(np.polyfit(np.log(np.array(ds, float)), np.log(vals), 1)[0])
                bybin[lab] = dict(n_genes=int(mask.sum()), alpha=vals, slope=sl)
                print(f"    mu∈[{lab:>14s}]  {mask.sum():>6d} 基因  alpha(d)="
                      f"{['%.4g' % v for v in vals]}  斜率 {sl:+.3f}")
            else:
                bybin[lab] = dict(n_genes=int(mask.sum()), alpha=vals, slope=None)
                print(f"    mu∈[{lab:>14s}]  {mask.sum():>6d} 基因  alpha(d)="
                      f"{['%.4g' % v for v in vals]}  斜率 不可估")

    out = dict(domain=a.domain, ds=ds, per_d=rows, slope_all=slope,
               slope_d8plus=slope_sub, segments=seg, by_expression=bybin,
               exclude_mt=bool(a.exclude_mt))
    suffix = "_noMT" if a.exclude_mt else ""
    fname = f"pooled_alpha{suffix}.json"
    json.dump(out, open(od / fname, "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print(f"\n已写出 {od/fname}")


if __name__ == "__main__":
    main()

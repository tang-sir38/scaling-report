r"""
34_layer_variance.py — 分层重建 + 方差分解：gamma 由表达量决定，还是由组织决定？

背景
----
逐基因拟合显示 gamma 从低表达端的 ~0.9 降到高表达端的 ~0.3，三套数据一致。
但"分层后的 pooled gamma"仍受层内少数基因的 mu^2 加权影响，所以同时报告
逐基因 gamma 的中位数（不受加权影响），两者互为验证。

分层
----
用 32 um 下每 bin 平均计数（绝对表达量），三套数据共用同一组对数边界，
从 1e-2 到 1e2，每 0.5 个数量级一层。每层报告基因数与两种 gamma。

方差分解
--------
对 (5 层 x 3 组织) 的 gamma 做 two-way 分解（无重复，交互项并入残差），
报告表达量、组织、残差各占的平方和比例。

输出 <outdir>/layer_variance.json
"""
import argparse, json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

SETS = {
    "colon": str(ROOT / "02_OUT" / "gene_moments.npz"),
    "lung": str(ROOT / "D3_LungCancer" / "02_OUT" / "gene_moments_noMT.npz"),
    "brain": str(ROOT / "D2_MouseBrain" / "02_OUT" / "gene_moments_noMT.npz"),
}

def load(p):
    z = np.load(p, allow_pickle=True)
    ds = np.array([int(x) for x in z["ds"]], float)
    g = np.array([str(x) for x in z["genes"]])
    MU = np.column_stack([z[f"mu_{int(d)}"] for d in ds])
    VA = np.column_stack([z[f"var_{int(d)}"] for d in ds])
    return ds, g, MU, VA

def pooled_layer(ds, MU, VA, mask):
    num = np.sum(VA[mask] - MU[mask], axis=0)
    den = np.sum(MU[mask] ** 2, axis=0)
    if np.any(den <= 0):
        return None, None
    al = num / den
    if np.any(al <= 0) or not np.all(np.isfinite(al)):
        return al, None
    return al, float(-np.polyfit(np.log(ds), np.log(al), 1)[0])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=str(ROOT / "02_OUT" / "layer_variance"))
    ap.add_argument("--lo", type=float, default=-2.0)
    ap.add_argument("--hi", type=float, default=2.0)
    ap.add_argument("--step", type=float, default=0.5)
    ap.add_argument("--ref-d", type=int, default=32)
    a = ap.parse_args()
    od = Path(a.outdir); od.mkdir(parents=True, exist_ok=True)

    edges = 10 ** np.arange(a.lo, a.hi + 1e-9, a.step)
    out = dict(edges=edges.tolist(), ref_d=a.ref_d, by_tissue={})

    for tag, p in SETS.items():
        ds, g, MU, VA = load(p)
        j = int(np.where(ds == a.ref_d)[0][0])
        AL = (VA - MU) / np.where(MU > 0, MU ** 2, np.nan)
        ok_g = np.all(np.isfinite(AL) & (AL > 0), axis=1) & (MU[:, -1] >= 0.5)
        xd = np.log(ds); xm = xd - xd.mean()
        gg = np.full(g.size, np.nan)
        gg[ok_g] = -((np.log(AL[ok_g]) * xm).sum(1) / (xm ** 2).sum())

        m32 = MU[:, j]
        tot2 = float(np.sum(m32 ** 2))
        rows = []
        for k in range(len(edges) - 1):
            lo, hi = edges[k], edges[k + 1]
            m = (m32 >= lo) & (m32 < hi)
            al, gam_pool = pooled_layer(ds, MU, VA, m)
            sel = m & ok_g
            rows.append(dict(lo=float(lo), hi=float(hi), n=int(m.sum()), n_gamma=int(sel.sum()),
                             gamma_pooled=gam_pool,
                             gamma_median=float(np.median(gg[sel])) if sel.sum() >= 20 else None,
                             share=float(np.sum(m32[m] ** 2) / tot2),
                             alpha=[float(x) for x in al] if al is not None and np.all(np.isfinite(al)) else None))
        # 最高一层：把 >10 counts/bin 的所有基因并成一行，表中据此补齐 100% 的权重
        mtop = m32 >= float(edges[-3])
        sel_t = mtop & ok_g
        al_t, gp_t = pooled_layer(ds, MU, VA, mtop)
        rows.append(dict(lo=float(edges[-3]), hi=None, n=int(mtop.sum()), n_gamma=int(sel_t.sum()),
                         gamma_pooled=gp_t,
                         gamma_median=float(np.median(gg[sel_t])) if sel_t.sum() >= 20 else None,
                         share=float(np.sum(m32[mtop] ** 2) / tot2),
                         alpha=None, merged=True))
        # 表外的低端：<0.01 counts/bin，只为说明剩余权重可忽略
        mbelow = m32 < float(edges[0])
        out.setdefault("below", {})[tag] = dict(n=int(mbelow.sum()),
                                               share=float(np.sum(m32[mbelow] ** 2) / tot2))
        out["by_tissue"][tag] = rows
        print(f"\n== {tag} ==  (mu_{a.ref_d} 分层, 每层 {a.step} 个数量级)")
        print(f"  {'层':>18s} {'基因':>6s} {'可拟合':>6s} {'gamma_pooled':>13s} {'gamma_median':>13s}")
        for r in rows:
            gp = f"{r['gamma_pooled']:.3f}" if r["gamma_pooled"] is not None else "  n.a."
            gm = f"{r['gamma_median']:.3f}" if r["gamma_median"] is not None else "  n.a."
            hi_s = f"{r['hi']:8.3g}" if r["hi"] is not None else "     inf"
            print(f"  {r['lo']:8.3g}-{hi_s} {r['n']:6d} {r['n_gamma']:6d} {gp:>13s} {gm:>13s} "
                  f"{100 * r['share']:6.1f}%")

    # ---- 方差分解（只用各组织都有基因的层，gamma_pooled）----
    labs = list(SETS)
    common = []
    for k in range(len(edges) - 1):
        vals = [out["by_tissue"][t][k]["gamma_pooled"] for t in labs]
        if all(v is not None for v in vals) and all(out["by_tissue"][t][k]["n"] >= 20 for t in labs):
            common.append((k, vals))
    if len(common) >= 3:
        Y = np.array([v for _, v in common])          # (n_layer, n_tissue)
        nl, nt = Y.shape
        gm = Y.mean()
        lay = Y.mean(axis=1) - gm
        tis = Y.mean(axis=0) - gm
        ss_lay = nt * float(np.sum(lay ** 2))
        ss_tis = nl * float(np.sum(tis ** 2))
        ss_tot = float(np.sum((Y - gm) ** 2))
        ss_res = ss_tot - ss_lay - ss_tis
        out["variance"] = dict(
            n_layers=int(nl), n_tissues=int(nt),
            ss_total=ss_tot, ss_layer=ss_lay, ss_tissue=ss_tis, ss_resid=ss_res,
            frac_layer=ss_lay / ss_tot, frac_tissue=ss_tis / ss_tot, frac_resid=ss_res / ss_tot,
            range_layer=float(lay.max() - lay.min()), range_tissue=float(tis.max() - tis.min()),
            layer_labels=[f"{edges[k]:.3g}-{edges[k+1]:.3g}" for k, _ in common],
            values=Y.tolist(), tissues=labs)
        v = out["variance"]
        print(f"\n== 方差分解（{nl} 层 x {nt} 组织，gamma_pooled）==")
        print(f"  表达量(层)  {100*v['frac_layer']:.1f}%   组织 {100*v['frac_tissue']:.1f}%   "
              f"残差(层x组织) {100*v['frac_resid']:.1f}%")
        print(f"  层间极差 {v['range_layer']:.3f}   组织间极差 {v['range_tissue']:.3f}   "
              f"比值 {v['range_layer']/v['range_tissue']:.1f}x")

        # 逐基因中位数的同一分解（独立于加权）
        Y2 = []
        for k, _ in common:
            row = [out["by_tissue"][t][k]["gamma_median"] for t in labs]
            Y2.append(row)
        Y2 = np.array(Y2)
        if np.all(np.isfinite(Y2)):
            gm2 = Y2.mean(); l2 = Y2.mean(1) - gm2; t2 = Y2.mean(0) - gm2
            s2 = np.sum((Y2 - gm2) ** 2)
            out["variance_median"] = dict(
                frac_layer=float(nt * np.sum(l2 ** 2) / s2),
                frac_tissue=float(nl * np.sum(t2 ** 2) / s2),
                range_layer=float(l2.max() - l2.min()),
                range_tissue=float(t2.max() - t2.min()))
            print(f"  [逐基因中位口径] 表达量 {100*out['variance_median']['frac_layer']:.1f}%  "
                  f"组织 {100*out['variance_median']['frac_tissue']:.1f}%  "
                  f"极差比 {out['variance_median']['range_layer']/out['variance_median']['range_tissue']:.1f}x")

    json.dump(out, open(od / "layer_variance.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(f"\n写出 {od/'layer_variance.json'}")

if __name__ == "__main__":
    main()

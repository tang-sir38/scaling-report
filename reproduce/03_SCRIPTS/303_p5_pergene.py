# -*- coding: utf-8 -*-
"""303_p5_pergene.py — 三片结肠（P1/P2/P5）逐基因 γ 的同口径对比。

口径照抄 41_gene_curves.py：8 档 α 全有限且为正、μ(256 µm) ≥ 0.5、γ = −slope(log α ~ log d)。
自检：P1 应复现正文的 0.618、P2 应复现 0.428；对不上就说明口径漂了，直接报错退出。
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

REPRO = Path(__file__).resolve().parents[1]
MIN_MU256 = 0.5

SETS = {
    "P2": REPRO / "02_OUT" / "gene_moments.npz",
    "P1": REPRO / "MS_Visium_HD_Human_Colon_Cancer_P1" / "02_OUT" / "gene_moments_noMT.npz",
    "P5": Path(r"E:\A_spatial_20261116\OUT\human_colon_cancer_p5\gene_moments__noMT.npz"),
}
EXPECT = {"P1": 0.618, "P2": 0.428}          # 正文 [46] 段已发表的两个中位数


def gamma_of(p: Path):
    z = np.load(p, allow_pickle=True)
    ds = np.array([int(x) for x in z["ds"]], float)
    genes = np.array([str(g) for g in z["genes"]])
    MU = np.column_stack([z[f"mu_{int(d)}"] for d in ds])
    VA = np.column_stack([z[f"var_{int(d)}"] for d in ds])
    AL = (VA - MU) / np.where(MU > 0, MU ** 2, np.nan)

    ok = np.all(np.isfinite(AL) & (AL > 0), axis=1) & (MU[:, -1] >= MIN_MU256)
    idx = np.where(ok)[0]
    xd = np.log(ds)
    xm = xd - xd.mean()
    b = (np.log(AL[idx]) * xm).sum(1) / (xm ** 2).sum()
    gamma = np.full(genes.size, np.nan)
    gamma[idx] = -b
    return dict(genes=genes, gamma=gamma, idx=idx, n_genes=int(genes.size),
                ds=ds, mu256=MU[:, -1])


def main() -> int:
    out = {}
    G = {}
    for tag, p in SETS.items():
        if not p.exists():
            raise SystemExit(f"缺 {tag}: {p}")
        r = gamma_of(p)
        g = r["gamma"][r["idx"]]
        G[tag] = r
        out[tag] = dict(n_genes=r["n_genes"], n_fit=int(r["idx"].size),
                        median=float(np.median(g)),
                        q25=float(np.percentile(g, 25)),
                        q75=float(np.percentile(g, 75)),
                        mean=float(g.mean()))
        print(f"{tag}: {r['n_genes']:,} genes -> {r['idx'].size:,} fit  "
              f"median γ = {np.median(g):.4f}  "
              f"IQR [{np.percentile(g,25):.3f}, {np.percentile(g,75):.3f}]")

    # 自检
    for tag, exp in EXPECT.items():
        got = out[tag]["median"]
        if abs(got - exp) > 0.002:
            raise SystemExit(f"口径漂移：{tag} per-gene median {got:.4f} != 正文 {exp}")

    # 共同基因集上的顺序一致性（P1 vs P2 vs P5 两两）
    common = set(G["P1"]["genes"][G["P1"]["idx"]])
    for t in ("P2", "P5"):
        common &= set(G[t]["genes"][G[t]["idx"]])
    common = sorted(common)
    print(f"\n三片共同可拟合基因 {len(common):,}")
    pair = {}
    from scipy.stats import spearmanr, pearsonr
    for a, b in (("P1", "P2"), ("P1", "P5"), ("P2", "P5")):
        ia = {g: i for i, g in enumerate(G[a]["genes"])}
        ib = {g: i for i, g in enumerate(G[b]["genes"])}
        va = np.array([G[a]["gamma"][ia[g]] for g in common])
        vb = np.array([G[b]["gamma"][ib[g]] for g in common])
        ok = np.isfinite(va) & np.isfinite(vb)
        rho = float(spearmanr(va[ok], vb[ok]).statistic)
        r = float(pearsonr(va[ok], vb[ok]).statistic)
        med_a, med_b = float(np.median(va[ok])), float(np.median(vb[ok]))
        pair[f"{a}_{b}"] = dict(n=int(ok.sum()), spearman=rho, pearson=r,
                                median_a=med_a, median_b=med_b,
                                gap=med_b - med_a)
        print(f"  {a} vs {b}: n={ok.sum():,}  ρ={rho:.4f}  r={r:.4f}  "
              f"median {med_a:.4f} vs {med_b:.4f}  Δ={med_b-med_a:+.4f}")

    res = dict(per_section=out, common_genes=len(common), pairs=pair,
               rule="41_gene_curves.py 同口径：8 档 α 全正且有限、mu256>=0.5",
               source={k: str(v) for k, v in SETS.items()})
    dst = Path(r"E:\A_spatial_20261116\OUT\human_colon_cancer_p5\pergene_three.json")
    dst.write_text(json.dumps(res, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"\n写入 {dst}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

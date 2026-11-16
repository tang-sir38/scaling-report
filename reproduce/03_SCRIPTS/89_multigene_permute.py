r"""89_multigene_permute.py — T5b：剂量-反应（置换多少基因，pooled γ 升多少）

问题
----
T5 只置换了 IGKC 一个基因，pooled γ 从 0.560 升到 0.873。这是"下界"。
如果所有基因都没有空间结构，pooled γ 应当回到计数统计极限 2（解析必然：
独立格下 alpha_g(d) = CV_g^2 / (d/2)^2，加权和仍正比于 d^-2）。

本脚本按"权重从大到小"逐步置换基因，测 pooled γ 的上升曲线。
若单调上升并趋近 2，"γ 由空间结构设定"就从单点证据变成剂量-反应关系。

设计
----
  取 top-N（按 2um 下的 mu^2 权重排序），N = 0, 10, 30, 100, 300, 1000
  对每个 N：把这 N 个基因在 2um 域内**各自独立随机置换**（每个基因的边缘分布不变）
            重算它们的 alpha_g(d)，其余基因沿用原值，按 mu^2 权重重组 pooled alpha
  输出 pooled gamma vs N

为什么其余基因可以沿用原值
--------------------------
  pooled alpha(d) = sum_g mu_g(d)^2 alpha_g(d) / sum_g mu_g(d)^2
  未置换的基因场没动，其 (mu, alpha) 完全不变，直接取现成值即可。

用法
----
  python 89_multigene_permute.py --ns 0,10,30,100,300,1000 --outdir 02_OUT/excise
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

ROOT = Path(__file__).resolve().parents[1]
KEY_STRIDE = 4096


def load_target_rows(npz_path, want_names):
    """只取需要的基因行，返回 {gene: (pos_r, pos_c, vals)}, 掩码 M, 网格。"""
    z = np.load(npz_path, allow_pickle=True)
    genes = np.array([str(g) for g in z["genes"]])
    shape = tuple(int(v) for v in z["shape"])
    coords = z["coords"].astype(np.int64)

    keys = z["dom_any"].astype(np.int64) if "dom_any" in z.files else None
    if keys is not None:
        kr, kc = keys // KEY_STRIDE, keys % KEY_STRIDE
        r0 = int(min(coords[:, 0].min(), kr.min())); c0 = int(min(coords[:, 1].min(), kc.min()))
        r1 = int(max(coords[:, 0].max(), kr.max())); c1 = int(max(coords[:, 1].max(), kc.max()))
    else:
        r0, c0 = int(coords[:, 0].min()), int(coords[:, 1].min())
        r1, c1 = int(coords[:, 0].max()), int(coords[:, 1].max())
    H, W = r1 - r0 + 1, c1 - c0 + 1

    M = np.zeros((H, W), dtype=bool)
    if keys is not None:
        M[kr - r0, kc - c0] = True
    else:
        M[coords[:, 0] - r0, coords[:, 1] - c0] = True
    dom_idx = np.flatnonzero(M)

    pos_r = coords[:, 0] - r0
    pos_c = coords[:, 1] - c0
    X = csr_matrix((z["data"], z["indices"], z["indptr"]), shape=shape)

    out = {}
    for g in want_names:
        hit = np.where(genes == g)[0]
        if hit.size == 0:
            continue
        row = X.getrow(int(hit[0]))
        out[g] = (pos_r[row.indices], pos_c[row.indices], row.data.astype(np.float64))
    return out, M, dom_idx, (H, W), genes


def alpha_after_permute(H, W, M, dom_idx, rr, cc, vals, ds, rng):
    """在 2um 域内随机置换该基因的取值，返回 {d: alpha}。"""
    F = np.zeros(H * W, dtype=np.float64)
    F[dom_idx] = 0.0
    lin = rr * W + cc
    # 把非零值放到域内随机位置
    pos = rng.choice(dom_idx.size, size=lin.size, replace=False)
    F[dom_idx[pos]] = vals
    F = F.reshape(H, W)
    out = {}
    for d in ds:
        k = max(int(round(d / 2)), 1)
        nb_h, nb_w = H // k, W // k
        if nb_h < 4 or nb_w < 4:
            continue
        B = F[: nb_h * k, : nb_w * k].reshape(nb_h, k, nb_w, k).sum(axis=(1, 3))
        Mb = M[: nb_h * k, : nb_w * k].reshape(nb_h, k, nb_w, k).any(axis=(1, 3))
        v = B[Mb]
        if v.size < 16:
            continue
        mu = v.mean()
        if mu <= 0:
            continue
        out[d] = float((v.var(ddof=1) - mu) / mu ** 2)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bins", default="01_WORK/bins_002um.npz")
    ap.add_argument("--moments", default="02_OUT/gene_moments.npz")
    ap.add_argument("--ns", default="0,10,30,100,300,1000")
    ap.add_argument("--outdir", default="02_OUT/excise")
    ap.add_argument("--seed", type=int, default=20261115)
    a = ap.parse_args()

    od = ROOT / a.outdir; od.mkdir(parents=True, exist_ok=True)
    ns = [int(x) for x in a.ns.split(",")]

    z = np.load(ROOT / a.moments, allow_pickle=True)
    ds = [int(x) for x in z["ds"]]
    genes_mo = np.array([str(g) for g in z["genes"]])
    MU = np.column_stack([z[f"mu_{d}"] for d in ds])
    VA = np.column_stack([z[f"var_{d}"] for d in ds])
    with np.errstate(invalid="ignore", divide="ignore"):
        AL = (VA - MU) / np.where(MU > 0, MU ** 2, np.nan)

    # 按 2um 的 mu^2 排序（决定谁主导 pooled）
    w2 = MU[:, 0] ** 2
    order = np.argsort(-w2)
    w2_sorted = w2[order]
    cum = np.cumsum(w2_sorted) / w2_sorted.sum()

    nmax = max(ns)
    top = order[:nmax]
    target_names = [genes_mo[i] for i in top]
    print(f"选 top-{nmax} 基因，其 2um mu^2 累积占比 {cum[nmax-1]:.1%}")

    rows, M, dom_idx, (H, W), all_genes = load_target_rows(ROOT / a.bins, target_names)
    print(f"  网格 {H}x{W}，域内 {dom_idx.size:,} 格，取到 {len(rows)} 个基因")

    rng = np.random.default_rng(a.seed)
    # 预先算好每个基因置换后的 alpha 曲线（每个基因只置换一次，N 之间复用）
    perm_alpha = {}
    for i, g in enumerate(target_names):
        if g not in rows:
            continue
        rr, cc, vals = rows[g]
        perm_alpha[g] = alpha_after_permute(H, W, M, dom_idx, rr, cc, vals, ds, rng)
        if (i + 1) % 100 == 0:
            print(f"    已置换 {i+1}/{len(target_names)}", flush=True)

    res = {"ds": ds, "cum_share": {}, "curves": {}}
    W2_all = np.nansum(MU ** 2, axis=0)
    S_all = np.nansum(MU ** 2 * AL, axis=0)

    for N in ns:
        idx = top[:N]
        names = [genes_mo[i] for i in idx]
        S = S_all.copy()
        hit = 0
        for i, g in zip(idx, names):
            if g not in perm_alpha:
                continue
            hit += 1
            for j, d in enumerate(ds):
                if d in perm_alpha[g]:
                    S[j] += MU[i, j] ** 2 * (perm_alpha[g][d] - AL[i, j])
        with np.errstate(invalid="ignore", divide="ignore"):
            ap = S / W2_all
        ok = np.isfinite(ap) & (ap > 0)
        if ok.sum() >= 3:
            x, y = np.log(np.array(ds)[ok]), np.log(ap[ok])
            A = np.vstack([x, np.ones_like(x)]).T
            coef, *_ = np.linalg.lstsq(A, y, rcond=None)
            r2 = 1 - ((y - A @ coef) ** 2).sum() / ((y - y.mean()) ** 2).sum()
            gam = float(-coef[0])
        else:
            gam, r2 = np.nan, np.nan
        share = float(cum[N - 1]) if N > 0 else 0.0
        res["curves"][str(N)] = dict(n_permuted=int(hit), mu2_share=share,
                                    gamma=gam, r2=float(r2),
                                    alpha_pool=ap.tolist())
        print(f"  N={N:<5} 置换 {hit:<5} 个  μ²累积 {share:6.1%}  "
              f"pooled γ={gam:.4f} (r2={r2:.3f})")

    (od / "multigene_permute.json").write_text(
        json.dumps(res, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n写出 {od}/multigene_permute.json")


if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
r"""
50_decomp_exact.py -- is the pooled exponent a weighted mean of the per-gene ones?

Two questions, both answered arithmetically rather than by argument.

(1) Identity.  At any fixed scale d, alpha_pool(d) = sum_g w_g(d) alpha_g(d) with
    w_g = mu_g^2 / sum mu^2.  This is exact by construction; we verify it against
    the direct estimator (sum(var-mu)/sum(mu^2)) to machine precision.

(2) Weighted mean.  gamma_pool is a slope of a sum, and the weights do not sit
    still as d grows, so gamma_pool is only approximately sum_g v_g gamma_g.
    We quantify the residual, and split the distance between the pooled value and
    the per-gene median into a composition part (the mu^2 weights) and a drift
    part (the weights moving with scale).

    v_g(d) = w_g(d) alpha_g(d) / alpha_pool(d)  is gene g's share of the pooled
    numerator at scale d; v_g = mean_d v_g(d) is its share of the reported number.
    n_eff = 1 / sum_g v_g^2 is the effective number of genes behind gamma_pool.
    The paper reports 1/sum w^2 instead, which is the same functional form applied
    to the wrong weights: w is the weight of the *variance*, v is the weight of the
    *exponent*.

Output: 02_OUT/decomp_exact/decomp_exact.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "02_OUT" / "decomp_exact"
DS = np.array([2, 4, 8, 16, 32, 64, 128, 256], float)
LGD = np.log(DS)
FILES = {
    "colon": ROOT / "02_OUT" / "gene_moments.npz",
    "lung": ROOT / "D3_LungCancer" / "02_OUT" / "gene_moments_noMT.npz",
    "brain": ROOT / "D2_MouseBrain" / "02_OUT" / "gene_moments_noMT.npz",
}


def load(path):
    z = np.load(path, allow_pickle=True)
    mu = np.column_stack([z[f"mu_{int(d)}"] for d in DS]).astype(np.float64)
    var = np.column_stack([z[f"var_{int(d)}"] for d in DS]).astype(np.float64)
    genes = np.array([str(g) for g in z["genes"]])
    return mu, var, genes


def slope(y):
    return -float(np.polyfit(LGD, np.log(y), 1)[0])


def one(tis, path):
    mu, var, genes = load(path)
    keep = mu.sum(axis=1) > 0                      # drop all-zero genes
    mu, var, genes = mu[keep], var[keep], genes[keep]
    G = mu.shape[0]

    with np.errstate(invalid="ignore", divide="ignore"):
        alpha_g = (var - mu) / mu ** 2
    mu2 = mu ** 2
    S2 = mu2.sum(axis=0)
    w = mu2 / S2                                   # mu^2 share, sums to 1 at every d

    # (1) identity check against the direct estimator
    alpha_pool = (w * alpha_g).sum(axis=0)
    direct = (var - mu).sum(axis=0) / S2
    identity_err = float(np.max(np.abs(alpha_pool - direct) / np.abs(direct)))

    gamma_pool = slope(alpha_pool)

    # gene filter as in the manuscript
    ok = np.all(np.isfinite(alpha_g) & (alpha_g > 0), axis=1) & (mu[:, -1] >= 0.5)
    gamma_g = np.array([slope(alpha_g[i]) if ok[i] else np.nan for i in range(G)])

    # (2) contribution shares and the weighted mean
    with np.errstate(invalid="ignore", divide="ignore"):
        v = w * alpha_g / alpha_pool
    v_bar = np.where(ok, v.mean(axis=1), 0.0)
    v_bar = v_bar / v_bar.sum()
    gamma_weighted = float(np.nansum(v_bar * gamma_g))
    rel_diff = float((gamma_weighted - gamma_pool) / gamma_pool)

    # effective number of genes, two conventions
    n_eff_mu2 = 1.0 / float((w ** 2).sum(axis=0).mean())
    n_eff_alpha = 1.0 / float((v_bar ** 2).sum())

    # heaviest contributors
    order = np.argsort(-v_bar)
    top1, top3 = int(order[0]), order[:3]
    top1_mu2_share = float(w[top1].mean())

    # drift: gamma_pool = sum(w0 gamma) - drift, evaluated at the finest scale
    w0 = w[:, 0]
    gamma_mu2 = float(np.nansum(np.where(ok, w0, 0.0) / np.where(ok, w0, 0.0).sum() * gamma_g))
    med = float(np.nanmedian(gamma_g[ok]))
    drift = gamma_weighted - gamma_pool
    gap = gamma_mu2 - med
    drift_share = float(drift / gap) if gap != 0 else float("nan")

    return dict(
        n_genes=int(G), n_fit=int(ok.sum()),
        gamma_pool=gamma_pool,
        gamma_weighted=gamma_weighted,
        rel_diff=rel_diff,
        gamma_mu2_weighted=gamma_mu2,
        gamma_median=med,
        identity_err=identity_err,
        n_eff_mu2=n_eff_mu2,
        n_eff_alpha=n_eff_alpha,
        n_eff_ratio=float(n_eff_mu2 / n_eff_alpha),
        top1_gene=str(genes[top1]),
        top1_alpha_share=float(v_bar[top1]),
        top1_mu2_share=top1_mu2_share,
        top3_alpha_share=float(v_bar[top3].sum()),
        top3_genes=[str(genes[i]) for i in top3],
        v_sum_fitted=float(v_bar.sum()),
        drift=drift,
        drift_share=drift_share,
        v_by_scale_top1=[float(v[top1, j]) for j in range(len(DS))],
    )


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    res = {}
    for tis in ("colon", "lung", "brain"):
        res[tis] = one(tis, FILES[tis])
        r = res[tis]
        print(f"\n===== {tis} =====")
        print(f"  identity max rel err      {r['identity_err']:.2e}")
        print(f"  gamma_pool                {r['gamma_pool']:.4f}")
        print(f"  sum(v*gamma)              {r['gamma_weighted']:.4f}"
              f"   ({r['rel_diff']*100:+.2f}%)")
        print(f"  mu2-weighted mean         {r['gamma_mu2_weighted']:.4f}"
              f"   median {r['gamma_median']:.4f}")
        print(f"  n_eff  mu2 {r['n_eff_mu2']:.1f}   alpha {r['n_eff_alpha']:.1f}"
              f"   ratio {r['n_eff_ratio']:.2f}")
        print(f"  top1 {r['top1_gene']}: alpha share {r['top1_alpha_share']*100:.1f}%"
              f" vs mu2 share {r['top1_mu2_share']*100:.1f}%")
        print(f"  top3 alpha share          {r['top3_alpha_share']*100:.1f}%")
        print(f"  drift {r['drift']:+.4f}  drift_share {r['drift_share']*100:.0f}%")

    p = OUT / "decomp_exact.json"
    p.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print("\nwritten", p)


if __name__ == "__main__":
    main()

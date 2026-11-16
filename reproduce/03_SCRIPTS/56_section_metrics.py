# -*- coding: utf-8 -*-
r"""
56_section_metrics.py -- one row per section, for the cross-section test.

For a dataset's per-gene per-scale moments this returns the quantities the
manuscript compares between sections:

  gamma_pool        slope of the pooled overdispersion over 2-256 um
  gamma_median      median of the per-gene exponents
  n_eff_mu2         1 / sum w^2 with w = mu_32^2 / sum mu_32^2
  n_eff_alpha       1 / sum v^2 with v the share of the pooled numerator
  top1_gene         the transcript carrying the largest share of sum mu^2
  top1_weight       that share
  top1_alpha_share  the same transcript's share of the exponent
  n_genes, n_fit

Usage:
  python 56_section_metrics.py --moments <npz> --label P1 --out <json>
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

DS = np.array([2, 4, 8, 16, 32, 64, 128, 256], float)
LGD = np.log(DS)


def slope(y):
    return -float(np.polyfit(LGD, np.log(y), 1)[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--moments", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    z = np.load(a.moments, allow_pickle=True)
    mu = np.column_stack([z[f"mu_{int(d)}"] for d in DS]).astype(np.float64)
    var = np.column_stack([z[f"var_{int(d)}"] for d in DS]).astype(np.float64)
    genes = np.array([str(g) for g in z["genes"]])
    keep = mu.sum(axis=1) > 0
    mu, var, genes = mu[keep], var[keep], genes[keep]
    G = mu.shape[0]

    with np.errstate(invalid="ignore", divide="ignore"):
        alpha_g = (var - mu) / mu ** 2
    mu2 = mu ** 2
    w = mu2 / mu2.sum(axis=0)
    alpha_pool = (w * alpha_g).sum(axis=0)
    gamma_pool = slope(alpha_pool)

    ok = np.all(np.isfinite(alpha_g) & (alpha_g > 0), axis=1) & (mu[:, -1] >= 0.5)
    gamma_g = np.array([slope(alpha_g[i]) if ok[i] else np.nan for i in range(G)])

    with np.errstate(invalid="ignore", divide="ignore"):
        v = w * alpha_g / alpha_pool
    v_bar = np.where(ok, v.mean(axis=1), 0.0)
    v_bar = v_bar / v_bar.sum()

    top1 = int(np.argmax(w[:, 4]))
    row = dict(
        label=a.label,
        n_genes=int(G),
        n_fit=int(ok.sum()),
        gamma_pool=gamma_pool,
        gamma_median=float(np.nanmedian(gamma_g[ok])),
        gamma_q25=float(np.nanpercentile(gamma_g[ok], 25)),
        gamma_q75=float(np.nanpercentile(gamma_g[ok], 75)),
        r2=float(1 - np.sum((np.log(alpha_pool) - np.polyval(
            np.polyfit(LGD, np.log(alpha_pool), 1), LGD)) ** 2)
                 / np.sum((np.log(alpha_pool) - np.log(alpha_pool).mean()) ** 2)),
        n_eff_mu2=1.0 / float((w ** 2).sum(axis=0).mean()),
        n_eff_alpha=1.0 / float((v_bar ** 2).sum()),
        top1_gene=str(genes[top1]),
        top1_weight=float(w[top1, 4]),
        top1_alpha_share=float(v_bar[top1]),
        alpha=[float(x) for x in alpha_pool],
        umi_total=float(mu.sum()),
    )
    p = Path(a.out)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(row, indent=2), encoding="utf-8")
    print(f"{a.label}: gamma_pool={gamma_pool:.4f}  median={row['gamma_median']:.4f}  "
          f"n_eff {row['n_eff_mu2']:.1f}/{row['n_eff_alpha']:.1f}  "
          f"top1={row['top1_gene']} {row['top1_weight']*100:.1f}% -> "
          f"{row['top1_alpha_share']*100:.1f}%")
    print("written", p)


if __name__ == "__main__":
    main()

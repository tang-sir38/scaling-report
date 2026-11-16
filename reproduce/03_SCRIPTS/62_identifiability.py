# -*- coding: utf-8 -*-
r"""
62_identifiability.py -- what the pooled exponent does not identify.

gamma_pool is a v-weighted mean of the per-gene exponents, with
v_g = w_g alpha_g / alpha_pool and v_bar its mean over scales (the same
definitions as 50_decomp_exact.py, which this script mirrors exactly).

Question: holding the reported gamma_pool fixed, how far can the median of the
per-gene exponents move? If the answer is "almost anywhere", the pooled value
does not identify the typical transcript and cannot rank tissues by their
transcript-level scaling.

Construction (existence proof, no optimisation): let gene 1 be the heaviest
contributor with share v1. Put gamma_1 = a and every other gene at one constant c:

    gamma_pool = v1 a + (1 - v1) c   =>   c = (gamma_pool - v1 a) / (1 - v1)

The remaining genes are the overwhelming majority, so median(gamma) = c, and
d(median)/da = -v1/(1-v1): the more concentrated the weights, the less the
pooled value constrains the typical gene.

Usage:
  python 62_identifiability.py [--out 02_OUT/identifiability/identifiability.json]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DS = np.array([2, 4, 8, 16, 32, 64, 128, 256], float)
LGD = np.log(DS)
FILES = {
    "colon": ROOT / "02_OUT" / "gene_moments.npz",
    "lung": ROOT / "D3_LungCancer" / "02_OUT" / "gene_moments_noMT.npz",
    "brain": ROOT / "D2_MouseBrain" / "02_OUT" / "gene_moments_noMT.npz",
}
_REF = ROOT / "02_OUT" / "decomp_exact" / "decomp_exact.json"
REF = json.loads(_REF.read_text(encoding="utf-8")) if _REF.exists() else {}


def slope(y: np.ndarray) -> float:
    return -float(np.polyfit(LGD, np.log(y), 1)[0])


def per_gene(path: Path):
    z = np.load(path, allow_pickle=True)
    mu = np.column_stack([z[f"mu_{int(d)}"] for d in DS]).astype(np.float64)
    var = np.column_stack([z[f"var_{int(d)}"] for d in DS]).astype(np.float64)
    genes = np.array([str(g) for g in z["genes"]])
    keep = mu.sum(axis=1) > 0
    mu, var, genes = mu[keep], var[keep], genes[keep]

    with np.errstate(invalid="ignore", divide="ignore"):
        alpha_g = (var - mu) / mu ** 2
    mu2 = mu ** 2
    S2 = mu2.sum(axis=0)
    w = mu2 / S2
    alpha_pool = (w * alpha_g).sum(axis=0)
    gamma_pool = slope(alpha_pool)

    ok = np.all(np.isfinite(alpha_g) & (alpha_g > 0), axis=1) & (mu[:, -1] >= 0.5)
    gamma_g = np.array([slope(alpha_g[i]) if ok[i] else np.nan for i in range(len(genes))])

    with np.errstate(invalid="ignore", divide="ignore"):
        v = w * alpha_g / alpha_pool
    v_bar = np.where(ok, v.mean(axis=1), 0.0)
    v_bar = v_bar / v_bar.sum()
    return dict(genes=genes, ok=ok, gamma=gamma_g, v=v_bar,
                gamma_pool=float(gamma_pool),
                gamma_weighted=float(np.nansum(v_bar * gamma_g)))


def analyse(tissue: str) -> dict:
    d = per_gene(FILES[tissue])
    ok, gamma, v = d["ok"], d["gamma"], d["v"]
    gpool = float(REF.get(tissue, {}).get("gamma_pool", d["gamma_pool"]))
    order = np.argsort(-np.where(ok, v, 0.0))
    top = int(order[0])
    v1 = float(v[top])
    med_obs = float(np.nanmedian(gamma[ok]))

    lo_c, hi_c = 0.0, 1.0
    a_lo = max((gpool - (1 - v1) * hi_c) / v1, 0.0)
    a_hi = min((gpool - (1 - v1) * lo_c) / v1, 1.0)
    med_hi = (gpool - v1 * a_lo) / (1 - v1)
    med_lo = (gpool - v1 * a_hi) / (1 - v1)

    return dict(
        tissue=tissue, scales_um=[int(x) for x in DS], n_genes=int(ok.sum()),
        gamma_pool=gpool, gamma_weighted_sum=d["gamma_weighted"],
        rel_diff=float((d["gamma_weighted"] - gpool) / gpool),
        gamma_median_observed=med_obs,
        top1_gene=str(d["genes"][top]), top1_alpha_share=v1,
        top1_gamma=float(gamma[top]),
        median_at_fixed_pool_min=float(med_lo),
        median_at_fixed_pool_max=float(med_hi),
        median_span=float(med_hi - med_lo),
        slope_dmedian_danchor=float(-v1 / (1 - v1)),
        anchor_range=[float(a_lo), float(a_hi)],
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="02_OUT/identifiability/identifiability.json")
    a = ap.parse_args()

    out = {}
    for t in FILES:
        r = out[t] = analyse(t)
        print(f"{t:6s} n={r['n_genes']:>6,}  gamma_pool={r['gamma_pool']:.4f}  "
              f"sum(v*gamma)={r['gamma_weighted_sum']:.4f} ({r['rel_diff']*100:+.1f}%)  "
              f"top1={r['top1_gene']} {r['top1_alpha_share']*100:.1f}%", flush=True)
        print(f"       observed median {r['gamma_median_observed']:.3f} | same pooled value "
              f"admits median {r['median_at_fixed_pool_min']:.3f}"
              f"..{r['median_at_fixed_pool_max']:.3f} (span {r['median_span']:.3f})", flush=True)
        print(f"       d(median)/d(anchor) = {r['slope_dmedian_danchor']:+.2f}", flush=True)

    p = Path(a.out)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"wrote {p}")


if __name__ == "__main__":
    main()

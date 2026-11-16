r"""Build a moments table with a known answer, to check the estimator.

    python examples/make_synthetic.py --out examples/synthetic_moments.npz

Construction. Mean counts grow with the area of the bin, mu_g(d) = m_g (d/2)^2,
and the per-gene overdispersion follows a power law by construction,

    alpha_g(d) = a_g (d/2)^(-gamma_g),   var_g = mu_g + alpha_g mu_g^2.

So gamma_g is known exactly for every gene and the estimator must return it to
machine precision. The pooled exponent is *not* the plain mean of the gene
exponents: the weights move with scale, which is the whole point.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

DS = np.array([2, 4, 8, 16, 32, 64, 128, 256], float)


def build(n_genes: int = 400, seed: int = 20261115):
    rng = np.random.default_rng(seed)
    # log-normal abundances: a few transcripts dominate, as in tissue
    mu256 = np.exp(rng.normal(np.log(30.0), 1.6, n_genes))
    mu256[:3] = [40000.0, 12000.0, 3000.0]
    gamma = rng.uniform(0.2, 0.9, n_genes).round(4)
    a256 = rng.uniform(0.05, 0.6, n_genes)

    mu = np.empty((n_genes, DS.size))
    var = np.empty_like(mu)
    for k, d in enumerate(DS):
        scale = (d / 256.0) ** 2          # counts per bin scale with area
        mu[:, k] = mu256 * scale
        alpha = a256 * (d / 256.0) ** (-gamma)
        var[:, k] = mu[:, k] + alpha * mu[:, k] ** 2
    genes = np.array([f"G{i:04d}" for i in range(n_genes)])
    return DS, genes, mu, var, {"gamma": gamma, "mu256": mu256}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path("examples/synthetic_moments.npz"))
    ap.add_argument("--n-genes", type=int, default=400)
    a = ap.parse_args()
    ds, genes, mu, var, truth = build(a.n_genes)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(a.out, ds=ds, genes=genes,
                        **{f"mu_{int(d)}": mu[:, k] for k, d in enumerate(ds)},
                        **{f"var_{int(d)}": var[:, k] for k, d in enumerate(ds)})
    np.savez_compressed(a.out.with_name("synthetic_truth.npz"),
                        genes=genes, gamma=truth["gamma"])
    print(f"wrote {a.out} ({a.out.stat().st_size/1e3:.0f} kB) and synthetic_truth.npz")


if __name__ == "__main__":
    main()

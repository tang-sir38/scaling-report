"""135_calibration.py -- quasi-likelihood calibration for binned counts.

Three layers, in increasing order of assumption:

  A. analytic   : given alpha, the Wald test's null distribution is N(0, phi)
                  exactly, so FPR and power are closed-form.  This layer carries
                  the claim; it is a theorem, not a simulation.
  B. simulated  : draw NB2 counts, estimate phi two ways (single-scale Pearson
                  vs eight-scale moments), and measure how much the estimation
                  error costs.  This layer is the second-order correction.
  C. clustered  : add a block-level log-normal frailty, and separately assign
                  group labels by whole 512 um block.  This is the pseudoreplication
                  scenario, and the one where a dispersion correction should fail.

Outputs out/calibration.json and the numbers used by the manuscript.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

DS = np.array([2, 4, 8, 16, 32, 64, 128, 256])
TOP_N = 2000
SEED = 20261125
OUT = Path(__file__).resolve().parent.parent / "out"
DATA = Path(__file__).resolve().parent.parent / "data"
OUT.mkdir(parents=True, exist_ok=True)

Z = 1.959963984540054


# ---------------------------------------------------------------- normal cdf
def phi_cdf(x):
    from math import erf, sqrt
    x = np.asarray(x, dtype=float)
    v = np.vectorize(lambda t: 0.5 * (1.0 + erf(t / sqrt(2.0))))
    return v(x)


def two_sided_p(z):
    return 2.0 * (1.0 - phi_cdf(np.abs(z)))


def fpr_from_phi(phi, alpha=0.05):
    """Two-sided FPR of a test whose statistic has variance phi under the null."""
    return float(2.0 * (1.0 - phi_cdf(Z / np.sqrt(phi))))


def phi_for_fpr(target, alpha=0.05):
    """Largest variance ratio that keeps the two-sided FPR at `target`."""
    from math import erf, sqrt

    def cdf(t):
        return 0.5 * (1.0 + erf(t / sqrt(2.0)))

    def ppf(p, lo=-10.0, hi=10.0, it=200):
        for _ in range(it):
            mid = 0.5 * (lo + hi)
            if cdf(mid) < p:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    zc = ppf(1.0 - alpha / 2.0)
    zt = ppf(1.0 - target / 2.0)
    return float((zc / zt) ** 2)


# ---------------------------------------------------------------- data
def load(path):
    d = np.load(path)
    mu = np.stack([d[f"mu_{x}"] for x in DS], axis=1).astype(float)
    var = np.stack([d[f"var_{x}"] for x in DS], axis=1).astype(float)
    n = d["n_cells"].astype(float)
    genes = d["genes"]
    return mu, var, n, genes


def alpha_from_moments(mu, var):
    with np.errstate(invalid="ignore", divide="ignore"):
        a = (var - mu) / mu ** 2
    a[~np.isfinite(a)] = np.nan
    return a


def top_genes(mu, j, n=TOP_N):
    col = mu[:, j]
    order = np.argsort(-col, kind="stable")
    order = order[col[order] > 0]
    return order[:n]


def per_gene_alpha_8scale(alpha, mu, j):
    """alpha_g(d) for every gene, read straight from the moments (no fit).

    The manuscript's gamma is a fit to log alpha vs log d; for the calibration we
    only need alpha at the scale being tested, which is measured, not fitted.
    """
    return alpha[:, j]


# ---------------------------------------------------------------- layer A
def analytic_table(mu, alpha, n_cells, j, n1=None, n2=None, deltas=(0.25, 0.5, 1.0)):
    """Closed-form FPR and power for Poisson vs quasi-Poisson at one scale.

    Null:  the statistic z_poisson = d / se_poisson has variance phi under H0,
           so |z| > 1.96 at rate 2[1-Phi(1.96/sqrt(phi))].
    Alt :  E[z] = delta / se_poisson and Var[z] = phi, so the same shift divided
           by sqrt(phi) applies to the QL statistic.
    """
    g = top_genes(mu, j)
    m = mu[g, j]
    a = alpha[g, j]
    ok = np.isfinite(a) & (a > 0) & (m > 0)
    m, a = m[ok], a[ok]

    N = n_cells[j]
    if n1 is None:
        n1 = max(30.0, 0.3 * N)
    if n2 is None:
        n2 = n1
    phi = 1.0 + a * m
    fpr_g = np.array([fpr_from_phi(float(x)) for x in phi])

    # Poisson SE of log(mu1/mu2) under the null, per gene
    se0 = np.sqrt(1.0 / (n1 * m) + 1.0 / (n2 * m))
    row = {
        "bin_um": float(DS[j]),
        "n_bins_total": int(N),
        "n1": float(n1),
        "n2": float(n2),
        "n_genes": int(m.size),
        "phi_median": float(np.median(phi)),
        "phi_mean": float(np.mean(phi)),
        "phi_q25": float(np.percentile(phi, 25)),
        "phi_q75": float(np.percentile(phi, 75)),
        "phi_q95": float(np.percentile(phi, 95)),
        # two readings of the same curve: the median gene and the average gene
        "fpr_poisson": fpr_from_phi(float(np.median(phi))),
        "fpr_poisson_mean": float(np.mean(fpr_g)),
        "fpr_poisson_q75": fpr_from_phi(float(np.percentile(phi, 75))),
        # QL with a known phi is calibrated by construction: 5%
        "fpr_ql": 0.05,
        "power": {},
    }
    order = np.argsort(m)
    thirds = np.array_split(order, 3)
    for d in deltas:
        z0 = d / se0
        p_pois = 1.0 - phi_cdf(Z - z0) + phi_cdf(-Z - z0)
        p_ql = 1.0 - phi_cdf(Z - z0 / np.sqrt(phi)) + phi_cdf(-Z - z0 / np.sqrt(phi))
        by_ab = {}
        for name, ix in zip(("low", "mid", "high"), thirds):
            by_ab[name] = {
                "poisson": float(np.mean(p_pois[ix])),
                "ql": float(np.mean(p_ql[ix])),
            }
        row["power"][f"{d}"] = {
            "poisson_median": float(np.median(p_pois)),
            "ql_median": float(np.median(p_ql)),
            "poisson_mean": float(np.mean(p_pois)),
            "ql_mean": float(np.mean(p_ql)),
            "loss_pp_median": float(100 * (np.median(p_pois) - np.median(p_ql))),
            "by_abundance": by_ab,
        }
    # equivalents: a corrected test needs phi times the bins for the same power
    row["n_equivalent_fraction"] = float(np.median(1.0 / phi))
    row["n_equivalent_fraction_q25"] = float(np.percentile(1.0 / phi, 25))
    return row


# ---------------------------------------------------------------- layer B
def nb_draw(rng, mu, alpha, size):
    """NB2 with mean mu and Var = mu + alpha*mu^2, shape (size, len(mu))."""
    alpha = np.asarray(alpha, dtype=float)
    mu = np.asarray(mu, dtype=float)
    # numpy's negative_binomial(n, p) has mean n(1-p)/p
    # choose n = 1/alpha, p = 1/(1+alpha*mu)
    npar = np.where(alpha > 0, 1.0 / np.maximum(alpha, 1e-12), 1e12)
    p = np.where(alpha > 0, 1.0 / (1.0 + alpha * mu), 1.0 - 1e-12)
    return rng.negative_binomial(npar, p, size=(size, mu.size)).astype(float)


def simulate_scale(rng, mu, alpha, j, n_genes=300, n_bin=500, reps=200,
                   delta=0.0, block_sd=0.0, n_blocks=None, min_count=20.0):
    """Empirical FPR / power with phi estimated two ways."""
    g = top_genes(mu, j)[:n_genes]
    m = mu[g, j]
    a = alpha[g, j]
    ok = np.isfinite(a) & (a > 0) & (m > 0) & (m * n_bin >= min_count)
    m, a = m[ok], a[ok]
    G = m.size
    if G < 10:
        return None

    n1 = n_bin
    n2 = n_bin
    mu1 = m
    mu2 = m * np.exp(delta)

    # optional block frailty: bins share a multiplier
    if block_sd > 0:
        nb = n_blocks or max(2, n_bin // 64)
        lab1 = rng.integers(0, nb, n1)
        lab2 = rng.integers(0, nb, n2)
        f1 = np.exp(rng.normal(0, block_sd, nb) - 0.5 * block_sd ** 2)
        f2 = np.exp(rng.normal(0, block_sd, nb) - 0.5 * block_sd ** 2)
        m1_eff = mu1[None, :] * f1[lab1][:, None]
        m2_eff = mu2[None, :] * f2[lab2][:, None]
    else:
        m1_eff = np.broadcast_to(mu1, (n1, G))
        m2_eff = np.broadcast_to(mu2, (n2, G))

    rej_p = np.zeros(G)
    rej_q = np.zeros(G)
    rej_qs = np.zeros(G)
    zs = np.empty((reps, G))
    for it in range(reps):
        y1 = nb_draw(rng, m1_eff.ravel(), np.broadcast_to(a, (n1, G)).ravel(), 1).reshape(n1, G)
        y2 = nb_draw(rng, m2_eff.ravel(), np.broadcast_to(a, (n2, G)).ravel(), 1).reshape(n2, G)
        # clip: a frailty can push a bin's mean below zero only through numerical noise
        m1h = np.maximum(y1.mean(0), 1e-9)
        m2h = np.maximum(y2.mean(0), 1e-9)
        se = np.sqrt(1.0 / (n1 * m1h) + 1.0 / (n2 * m2h))
        z = np.log(m2h / m1h) / se
        zs[it] = z

        # phi, way 1: eight-scale moments (what the paper has)
        phi_mom = 1.0 + a * np.sqrt(m1h * m2h)
        # phi, way 2: single-scale Pearson, pooled over both arms
        v1 = y1.var(0, ddof=1)
        v2 = y2.var(0, ddof=1)
        phi_pear = (v1 / m1h + v2 / m2h) / 2.0
        phi_pear = np.clip(phi_pear, 1.0, None)

        # phi, way 3: Pearson shrunk to the median (simple empirical Bayes)
        med = np.median(phi_pear)
        phi_shr = 0.5 * phi_pear + 0.5 * med

        rej_p += np.abs(z) > Z
        rej_q += np.abs(z / np.sqrt(phi_mom)) > Z
        rej_qs += np.abs(z / np.sqrt(phi_shr)) > Z

    # stratified refusal rate, low / mid / high abundance
    order = np.argsort(m)
    thirds = np.array_split(order, 3)
    strat = {name: float((rej_p[ix] / reps).mean())
             for name, ix in zip(("low", "mid", "high"), thirds)}
    phi_true = 1.0 + a * m
    fpr_analytic = np.array([fpr_from_phi(float(x)) for x in phi_true])

    return {
        "bin_um": float(DS[j]),
        "n_bin": int(n_bin),
        "n_genes": int(G),
        "reps": int(reps),
        "delta": float(delta),
        "block_sd": float(block_sd),
        "fpr_poisson": float((rej_p / reps).mean()),
        "fpr_ql_moments": float((rej_q / reps).mean()),
        "fpr_ql_shrunk": float((rej_qs / reps).mean()),
        "phi_median_true": float(np.median(1.0 + a * m)),
        # the same genes, closed form: agreement here means the simulation is clean
        "fpr_analytic_same_genes": float(np.mean(fpr_analytic)),
        # diagnostics: under the null the variance of the Poisson z should equal phi
        "z_var_empirical": float(np.mean(zs.var(0))),
        "z_var_over_phi": float(np.mean(zs.var(0) / (1.0 + a * m))),
        "fpr_poisson_by_abundance": strat,
        "fpr_poisson_by_abundance_ql": {name: float((rej_q[ix] / reps).mean())
                                        for name, ix in zip(("low", "mid", "high"), thirds)},
    }


# ---------------------------------------------------------------- layer C
def block_assignment(rng, mu, alpha, j, n_bin=500, reps=100, delta=0.0,
                     block_sd=0.5, block_size=16):
    """Groups assigned by whole block: the effective n is the number of blocks.

    This is the scenario a dispersion correction cannot fix, because the
    inflation is in the design, not in the count distribution.
    """
    g = top_genes(mu, j)[:300]
    m = mu[g, j]
    a = alpha[g, j]
    ok = np.isfinite(a) & (a > 0) & (m > 0)
    m, a = m[ok], a[ok]
    G = m.size
    if G == 0:
        return None
    nb = max(4, (2 * n_bin) // block_size)
    lab = np.repeat(np.arange(nb), block_size)
    if lab.size < 2 * n_bin:                      # pad the tail into the last block
        lab = np.concatenate([lab, np.full(2 * n_bin - lab.size, nb - 1)])
    lab = lab[: 2 * n_bin]
    # whole blocks are assigned to one arm or the other
    grp = np.random.default_rng(rng.integers(1 << 30)).integers(0, 2, nb)[lab]
    n1, n2 = int((grp == 0).sum()), int((grp == 1).sum())

    rej_p = np.zeros(G)
    rej_q = np.zeros(G)
    rej_b = np.zeros(G)
    blk_grp = np.array([grp[lab == b][0] for b in range(nb)])
    n1b, n2b = int((blk_grp == 0).sum()), int((blk_grp == 1).sum())
    tcrit = 1.96 if min(n1b, n2b) > 30 else 2.05   # t(60) ~ 2.00; conservative at 50
    for _ in range(reps):
        f = np.exp(rng.normal(0, block_sd, nb) - 0.5 * block_sd ** 2)
        mean_bin = m[None, :] * f[lab][:, None]
        mean_bin[grp == 1] *= np.exp(delta)
        y = nb_draw(rng, mean_bin.ravel(), np.broadcast_to(a, mean_bin.shape).ravel(), 1).reshape(mean_bin.shape)
        m1h = np.maximum(y[grp == 0].mean(0), 1e-9)
        m2h = np.maximum(y[grp == 1].mean(0), 1e-9)
        se = np.sqrt(1.0 / (n1 * m1h) + 1.0 / (n2 * m2h))
        z = np.log(m2h / m1h) / se
        phi = 1.0 + a * np.sqrt(m1h * m2h)
        rej_p += np.abs(z) > Z
        rej_q += np.abs(z / np.sqrt(phi)) > Z

        # block as the unit: sum within block, then Welch on log block means
        tot = np.zeros((nb, G))
        np.add.at(tot, lab, y)
        cnt = np.bincount(lab, minlength=nb)[:, None]
        bmean = tot / np.maximum(cnt, 1)
        b1 = bmean[blk_grp == 0]
        b2 = bmean[blk_grp == 1]
        mu1b = np.maximum(b1.mean(0), 1e-9)
        mu2b = np.maximum(b2.mean(0), 1e-9)
        s1 = b1.std(0, ddof=1)
        s2 = b2.std(0, ddof=1)
        se_b = np.sqrt((s1 ** 2) / n1b / mu1b ** 2 + (s2 ** 2) / n2b / mu2b ** 2)
        se_b = np.maximum(se_b, 1e-12)
        zb = np.log(mu2b / mu1b) / se_b
        rej_b += np.abs(zb) > tcrit

    return {
        "bin_um": float(DS[j]),
        "n_blocks": int(nb),
        "n1": n1, "n2": n2,
        "n1_blocks": n1b, "n2_blocks": n2b,
        "block_sd": float(block_sd),
        "reps": int(reps),
        "delta": float(delta),
        "fpr_poisson": float((rej_p / reps).mean()),
        "fpr_ql_moments": float((rej_q / reps).mean()),
        "fpr_blocklevel": float((rej_b / reps).mean()),
    }


# ---------------------------------------------------------------- driver
def run(prefix, path, quick=False):
    mu, var, n, genes = load(path)
    alpha = alpha_from_moments(mu, var)
    out = {"section": prefix, "ds": DS.tolist(), "n_cells": n.tolist(),
           "analytic": [], "simulated": [], "blocked": []}

    for j in range(len(DS)):
        out["analytic"].append(analytic_table(mu, alpha, n, j))

    reps = 60 if quick else 200
    for j in [2, 3, 4, 5, 6]:
        rng = np.random.default_rng(SEED + 100 * j)
        nsim = int(min(n[j], 600))
        rows = []
        for delta in (0.0, 1.0):
            r = simulate_scale(rng, mu, alpha, j, n_bin=nsim, reps=reps, delta=delta)
            if r:
                rows.append(r)
        out["simulated"] += rows

    for j in [3, 4, 5]:
        rng = np.random.default_rng(SEED + 900 + j)
        for sd in (0.0, 0.3, 0.6):
            r = simulate_scale(rng, mu, alpha, j, n_bin=400, reps=60,
                               delta=0.0, block_sd=sd, n_blocks=25)
            if r:
                out["simulated"].append({**r, "kind": "frailty"})
        r = block_assignment(rng, mu, alpha, j, n_bin=400, reps=60, block_sd=0.5)
        if r:
            out["blocked"].append(r)

    return out


def main():
    quick = "--quick" in sys.argv
    res = {}
    for pre, fn in (("P1", "colon_P1_moments.npz"), ("P2", "colon_P2_moments.npz")):
        res[pre] = run(pre, DATA / fn, quick=quick)
        print(f"[{pre}] done", flush=True)

    (OUT / "calibration.json").write_text(
        json.dumps(res, indent=2, ensure_ascii=False), encoding="utf-8")

    print("\n=== analytic, phi_median and implied FPR ===")
    for pre in res:
        for row in res[pre]["analytic"]:
            print(f"{pre} {row['bin_um']:>4.0f}um  phi_med={row['phi_median']:8.3f} "
                  f"FPR_med={row['fpr_poisson']*100:5.1f}%  "
                  f"FPR_mean={row['fpr_poisson_mean']*100:5.1f}%  n={row['n_bins_total']:,}")

    print("\n=== power, median gene, delta=0.25 (1.28x) ===")
    for pre in res:
        for row in res[pre]["analytic"]:
            p = row["power"]["0.25"]
            print(f"{pre} {row['bin_um']:>4.0f}um  P={p['poisson_median']:.3f}->{p['ql_median']:.3f} "
                  f"loss={p['loss_pp_median']:+5.1f}pp  "
                  f"n_eff={row['n_equivalent_fraction']*100:4.1f}%")

    print("\n=== simulated, null (delta=0) ===")
    for pre in res:
        for r in res[pre]["simulated"]:
            if r["delta"] == 0 and r.get("kind") != "frailty":
                print(f"{pre} {r['bin_um']:>4.0f}um  P={r['fpr_poisson']*100:5.1f}%  "
                      f"QLmom={r['fpr_ql_moments']*100:5.1f}%  QLshr={r['fpr_ql_shrunk']*100:5.1f}%")

    print("\n=== block assignment (pseudoreplication) ===")
    for pre in res:
        for r in res[pre]["blocked"]:
            print(f"{pre} {r['bin_um']:>4.0f}um  P={r['fpr_poisson']*100:5.1f}%  "
                  f"QL={r['fpr_ql_moments']*100:5.1f}%  "
                  f"block-level={r['fpr_blocklevel']*100:5.1f}%  "
                  f"blocks={r['n1_blocks']}+{r['n2_blocks']}")
    print("\nwrote", OUT / "calibration.json")


if __name__ == "__main__":
    main()

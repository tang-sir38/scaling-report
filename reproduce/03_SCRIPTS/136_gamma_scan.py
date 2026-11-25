"""136_gamma_scan.py -- separate the two drivers of the cost, v2.

Why v2: v1 predicted phi(256) by extrapolating a per-gene power law fitted over
2-256 um, and got 9.8 against a measured 29.5.  That is the same failure the
manuscript reports in S2 Text (predicted exponent steeper than observed by a
factor of 2.4-8.2), and it means a closed-form d* must not be offered as a
prescription.  v2 therefore reports

  1. the measured phi at every scale (lookup, no extrapolation),
  2. a gamma scan in which only gamma moves and mu(256) is the measured value,
  3. where the power law itself stops being accurate (predicted vs measured at
     every scale), which is the operational form of "price the test from gamma",
  4. the depth scan that explains why P1 loses more power than P2 despite a
     smaller phi, anchored on the measured phi(256).
"""
from __future__ import annotations

import json
from math import erf, sqrt
from pathlib import Path

import numpy as np

DS = np.array([2, 4, 8, 16, 32, 64, 128, 256])
TOP_N = 2000
OUT = Path(__file__).resolve().parent.parent / "out"
DATA = Path(__file__).resolve().parent.parent / "data"
Z = 1.959963984540054


def cdf(x):
    return 0.5 * (1.0 + erf(x / sqrt(2.0)))


def cdfv(x):
    """Vectorised standard normal cdf (math.erf is scalar only)."""
    x = np.asarray(x, dtype=float)
    return 0.5 * (1.0 + np.vectorize(erf)(x / sqrt(2.0)))


def ppf(p, lo=-10.0, hi=10.0):
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if cdf(mid) < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def fpr(phi):
    return 2.0 * (1.0 - cdf(Z / sqrt(phi)))


def load(path):
    d = np.load(path)
    mu = np.stack([d[f"mu_{x}"] for x in DS], 1).astype(float)
    var = np.stack([d[f"var_{x}"] for x in DS], 1).astype(float)
    return mu, var, d["n_cells"].astype(float)


def median_gamma(mu, var):
    """Per-gene slope of log alpha on log d, median over genes fitted at all scales."""
    with np.errstate(invalid="ignore", divide="ignore"):
        a = (var - mu) / mu ** 2
    ok = np.all(np.isfinite(a) & (a > 0), axis=1)
    lg = np.log(a[ok])
    x = np.log(DS.astype(float))
    x = x - x.mean()
    g = (lg * x).sum(1) / (x ** 2).sum()
    return float(np.median(-g)), int(ok.sum())


def measured(mu, var, j, n_top=TOP_N):
    """Median mu and median phi over the n_top most abundant genes at scale j."""
    m = mu[:, j]
    with np.errstate(invalid="ignore", divide="ignore"):
        a = (var[:, j] - m) / m ** 2
    order = np.argsort(-m)
    order = order[m[order] > 0][:n_top]
    keep = order[np.isfinite(a[order]) & (a[order] > 0)]
    mm, aa = m[keep], a[keep]
    phi = 1.0 + aa * mm
    return float(np.median(mm)), float(np.median(phi))


def main():
    res = {}
    tab = {}
    for pre, fn in (("P1", "colon_P1_moments.npz"), ("P2", "colon_P2_moments.npz")):
        mu, var, n = load(DATA / fn)
        g_med, n_fit = median_gamma(mu, var)
        mu_by, phi_by = [], []
        for j in range(len(DS)):
            m, p = measured(mu, var, j)
            mu_by.append(m)
            phi_by.append(p)
        base = phi_by[0] - 1.0
        mu256 = mu_by[7]
        nn = max(30.0, 0.3 * n[7])
        se = sqrt(2.0 / (nn * mu256))          # log-ratio SE at 256 um, median gene

        res[pre] = {
            "gamma_median_per_gene": g_med,
            "n_genes_fitted": n_fit,
            "n_cells": n.tolist(),
            "mu_median_by_scale": mu_by,
            "phi_median_by_scale": phi_by,
            "phi_2_minus_1": base,
            "mu_median_256": mu256,
            "gamma_scan": [],
            "depth_scan": [],
            "effective_n_by_scale": [],
        }

        # how far does the power law carry?  predict phi(256) from phi(2) and the
        # median gamma, then compare with the measured value
        pred256 = 1.0 + base * (256.0 / 2.0) ** (2.0 - g_med)
        res[pre]["phi_256_predicted_from_gamma"] = float(pred256)
        res[pre]["phi_256_measured"] = phi_by[7]
        res[pre]["extrapolation_ratio"] = float(pred256 / phi_by[7])

        for gam in (0.30, 0.364, 0.45, 0.50, 0.554, 0.60, 0.70):
            phi256 = 1.0 + base * (256.0 / 2.0) ** (2.0 - gam)
            z0 = 0.25 / se
            p_p = 1 - cdf(Z - z0) + cdf(-Z - z0)
            p_q = 1 - cdf(Z - z0 / sqrt(phi256)) + cdf(-Z - z0 / sqrt(phi256))
            res[pre]["gamma_scan"].append({
                "gamma": gam,
                "phi_256_if_power_law_held": float(phi256),
                "power_poisson": float(p_p),
                "power_ql": float(p_q),
                "loss_pp": float(100 * (p_p - p_q)),
                "n_effective_fraction": float(1.0 / phi256),
            })

        # Depth scan.  Anchored on the MEASURED phi(256), not on the power-law
        # prediction: with the prediction as anchor the mult=1 row returned a
        # 2.4 pp loss where the observed curve shows 37.5 pp, i.e. the table
        # contradicted the figure printed one page later.
        #
        # Evaluated gene by gene and then summarised by the median, which is the
        # same estimator 135_calibration.py uses; at mult=1 this reproduces its
        # loss_pp_median to the last digit.  Scaling depth multiplies mu and
        # leaves alpha alone (alpha is a property of the field and the bin), so
        # phi - 1 scales with the depth multiplier.
        m7 = mu[:, 7]
        with np.errstate(invalid="ignore", divide="ignore"):
            a7 = (var[:, 7] - m7) / m7 ** 2
        order7 = np.argsort(-m7, kind="stable")
        order7 = order7[m7[order7] > 0][:TOP_N]
        keep7 = order7[np.isfinite(a7[order7]) & (a7[order7] > 0)]
        mm7, aa7 = m7[keep7], a7[keep7]
        n1 = max(30.0, 0.3 * n[7])
        for mult in (0.5, 1.0, 2.0, 4.0):
            m_eff = mm7 * mult
            phi_eff = 1.0 + aa7 * m_eff
            se0 = np.sqrt(2.0 / (n1 * m_eff))
            z0 = 0.25 / se0
            p_p = 1 - cdfv(Z - z0) + cdfv(-Z - z0)
            p_q = 1 - cdfv(Z - z0 / np.sqrt(phi_eff)) + cdfv(-Z - z0 / np.sqrt(phi_eff))
            res[pre]["depth_scan"].append({
                "depth_multiplier": mult,
                "phi_256_median": float(np.median(phi_eff)),
                "power_poisson_median": float(np.median(p_p)),
                "power_ql_median": float(np.median(p_q)),
                "loss_pp": float(100 * (np.median(p_p) - np.median(p_q))),
            })
        res[pre]["depth_scan_anchor"] = "measured phi(256), per-gene medians"

        # Where does the power law stop being usable?  Predicted against measured
        # phi at every scale.  This is the operational version of "you can price
        # the test from gamma": it holds for the first two doublings and has to be
        # abandoned after that.
        fid = []
        for j in range(len(DS)):
            pred = 1.0 + base * (DS[j] / 2.0) ** (2.0 - g_med)
            meas = phi_by[j]
            fid.append({
                "bin_um": float(DS[j]),
                "phi_predicted": float(pred),
                "phi_measured": float(meas),
                "ratio_pred_over_meas": float(pred / meas),
                "rel_error": float((pred - meas) / meas),
            })
        res[pre]["power_law_fidelity"] = fid

        # The prospective version of the same question.  A user who has measured
        # down to 16 um does not have phi at 256 um; they have to extrapolate.
        # Fit gamma on 2..U only and price the coarser scales from it.
        pros = []
        for U in (8.0, 16.0, 32.0):
            sel = DS <= U
            with np.errstate(invalid="ignore", divide="ignore"):
                a_all = (var - mu) / mu ** 2
            okU = np.all(np.isfinite(a_all[:, sel]) & (a_all[:, sel] > 0), axis=1)
            lg = np.log(a_all[okU][:, sel])
            xU = np.log(DS[sel].astype(float))
            xU = xU - xU.mean()
            gU = float(np.median(-(lg * xU).sum(1) / (xU ** 2).sum()))
            for j in range(len(DS)):
                if DS[j] <= U:
                    continue
                pred = 1.0 + base * (DS[j] / 2.0) ** (2.0 - gU)
                pros.append({
                    "fit_window_um": [2.0, float(U)],
                    "gamma_fitted": gU,
                    "target_um": float(DS[j]),
                    "phi_predicted": float(pred),
                    "phi_measured": float(phi_by[j]),
                    "factor_pred_over_meas": float(pred / phi_by[j]),
                })
        res[pre]["prospective_extrapolation"] = pros

        # The effective number of bins at every scale.  Power itself is not
        # computed here: 135_calibration.py owns that number, and two estimators
        # for the same quantity is how a paper ends up quoting two different
        # costs for one scale.
        for j in range(len(DS)):
            phi = phi_by[j]
            res[pre]["effective_n_by_scale"].append({
                "bin_um": float(DS[j]),
                "phi": float(phi),
                "n_bins": float(n[j]),
                "n_effective_fraction": float(1.0 / phi),
                "n_effective": float(n[j] / phi),
            })

        tab[pre] = list(zip(DS.tolist(), mu_by, phi_by))

    (OUT / "gamma_scan.json").write_text(json.dumps(res, indent=2), encoding="utf-8")

    for pre in res:
        r = res[pre]
        print(f"\n== {pre} ==")
        print(f"  median per-gene gamma = {r['gamma_median_per_gene']:.3f} "
              f"(n={r['n_genes_fitted']:,})")
        print(f"  measured phi by scale: " +
              " ".join(f"{d:g}:{p:.2f}" for d, p in zip(DS, r["phi_median_by_scale"])))
        print(f"  power law from gamma predicts phi(256) = "
              f"{r['phi_256_predicted_from_gamma']:.1f} against measured "
              f"{r['phi_256_measured']:.1f}  (ratio {r['extrapolation_ratio']:.2f})")
        print("  gamma scan at 256 um (mu held at its measured value):")
        print("    gamma  phi*   P(QL)  loss    n_eff")
        for s in r["gamma_scan"]:
            print(f"    {s['gamma']:.3f} {s['phi_256_if_power_law_held']:7.1f} "
                  f"{s['power_ql']:.3f} {s['loss_pp']:+5.1f} "
                  f"{s['n_effective_fraction']*100:5.1f}%")
        print("  effective number of bins with the measured phi:")
        print("    bin    phi     n_bin   n_eff   n_eff%")
        for s in r["effective_n_by_scale"]:
            print(f"    {s['bin_um']:4.0f} {s['phi']:7.2f} {s['n_bins']:9.0f} "
                  f"{s['n_effective']:7.1f} {s['n_effective_fraction']*100:6.1f}%")
        print("  depth scan (anchored on measured phi(256)):")
        for s in r["depth_scan"]:
            print(f"    x{s['depth_multiplier']:.1f} phi={s['phi_256_median']:7.2f} "
                  f"P(QL)={s['power_ql_median']:.3f} loss={s['loss_pp']:+5.1f}")
        print("  prospective: fit gamma on 2..U, then extrapolate")
        for s in r["prospective_extrapolation"]:
            print(f"    window 2-{s['fit_window_um'][1]:3.0f} (gamma={s['gamma_fitted']:.3f}) "
                  f"-> {s['target_um']:3.0f} um  factor={s['factor_pred_over_meas']:5.2f}")
        print("  where the power law stops working:")
        for s in r["power_law_fidelity"]:
            print(f"    {s['bin_um']:4.0f} pred={s['phi_predicted']:7.2f} "
                  f"meas={s['phi_measured']:7.2f} err={100*s['rel_error']:+6.1f}%")
    print("\nwrote", OUT / "gamma_scan.json")


if __name__ == "__main__":
    main()

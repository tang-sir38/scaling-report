"""Turn the exponent into a decision.

A reader who measures gamma = 0.55 has learned something true and, on its own,
useless.  The question behind the measurement is always one of:

    May I use this bin size with a Poisson test?
    If not, what does the mistake cost me?
    What is the largest bin I can use and still trust a p-value?

This module answers those three, from the same numbers the report already
prints.  Nothing here is a new estimate; it is arithmetic on alpha_pool(d) and
the per-gene moments, and every quantity it returns is stated with the scale it
belongs to.

The model
---------
A negative-binomial count with mean mu and overdispersion alpha has variance
mu + alpha*mu^2, so its variance-to-mean ratio is

    phi = 1 + alpha * mu

A test that assumes Poisson divides by sqrt(mu) where it should divide by
sqrt(phi*mu), so its z statistic is inflated by sqrt(phi).  A two-sided test at
level a then rejects at

    FPR(phi) = 2 * (1 - Phi(z_crit / sqrt(phi))),    z_crit = Phi^-1(1 - a/2)

which is 5% at phi = 1 and rises steeply: phi = 30 gives 72%, the value the
paper reports at 256 um.
"""
from __future__ import annotations

import math

import numpy as np

from .core import Fit, Moments

# --------------------------------------------------------------------------
# normal distribution, standard library only
# --------------------------------------------------------------------------


def _phi_cdf(z: float) -> float:
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def _phi_ppf(p: float, lo: float = -10.0, hi: float = 10.0,
             iters: int = 200) -> float:
    """Inverse normal CDF by bisection; p in (0, 1)."""
    if not 0.0 < p < 1.0:
        raise ValueError("p must be in (0, 1)")
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if _phi_cdf(mid) < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def fpr_under_poisson(phi: float, alpha: float = 0.05) -> float:
    """Two-sided false-positive rate of a Poisson test when the true ratio is phi."""
    if phi <= 0:
        raise ValueError("phi must be positive")
    z_crit = _phi_ppf(1.0 - alpha / 2.0)
    return float(2.0 * (1.0 - _phi_cdf(z_crit / math.sqrt(phi))))


def phi_for_fpr(fpr: float, alpha: float = 0.05) -> float:
    """The largest variance ratio that still keeps the false-positive rate at fpr."""
    if not 0.0 < fpr < 1.0:
        raise ValueError("fpr must be in (0, 1)")
    z_crit = _phi_ppf(1.0 - alpha / 2.0)
    z = _phi_ppf(1.0 - fpr / 2.0)
    if z <= 0:
        return float("inf")
    return float((z_crit / z) ** 2)


# --------------------------------------------------------------------------
# per-gene dispersion at each scale
# --------------------------------------------------------------------------


def per_gene_phi(m: Moments, f: Fit) -> np.ndarray:
    """phi_g(d) = 1 + alpha_g(d) * mu_g(d), shape (G, n_scales).

    Genes with no counts at a scale are left at 1 (Poisson), which is what a
    test on that gene would effectively see.
    """
    with np.errstate(invalid="ignore", divide="ignore"):
        alpha = (m.var - m.mu) / m.mu ** 2
        phi = 1.0 + alpha * m.mu
    phi = np.where(np.isfinite(phi), phi, 1.0)
    return np.maximum(phi, 1.0)


# The paper reports phi as the median over the 2,000 most abundant genes at each
# scale (S2 Fig legend), not over every fitted gene. The two differ by a factor
# of about five at the coarsest bin, because the fitted set is dominated by
# low-count genes whose Poisson term still hides the overdispersion. The default
# is the paper's definition so that a report from this tool is comparable with
# the published numbers; the fitted-gene median is carried alongside it.
PAPER_TOP_N = 2000


def scale_table(m: Moments, f: Fit, *, fpr_target: float = 0.10,
                top_n: int = PAPER_TOP_N) -> list[dict]:
    """One row per bin size: what a Poisson test would actually do there.

    ``phi_median`` is the median over the ``top_n`` most abundant genes, which is
    the quantity the paper reports. ``phi_median_all_fitted`` is the median over
    every gene that entered the per-gene fit, kept because it is the more
    conservative of the two and a reader may want it.
    """
    phi = per_gene_phi(m, f)
    ok = f.ok if f.ok.any() else np.ones(m.n_genes, bool)
    rows = []
    for j, d in enumerate(m.ds):
        col = phi[:, j]
        expressed = m.mu[:, j] > 0
        if not expressed.any():
            continue
        # the paper's set: the most abundant genes at this scale
        order = np.argsort(-m.mu[:, j], kind="stable")
        order = order[m.mu[order, j] > 0][:top_n]
        top = col[order]
        top = top[np.isfinite(top) & (top > 0)]
        if top.size == 0:
            continue
        fitted = col[ok]
        fitted = fitted[np.isfinite(fitted) & (fitted > 0)]
        med = float(np.median(top))
        rows.append(dict(
            bin_um=float(d),
            phi_median=med,
            phi_q25=float(np.percentile(top, 25)),
            phi_q75=float(np.percentile(top, 75)),
            phi_median_all_fitted=float(np.median(fitted)) if fitted.size else None,
            fpr_median=fpr_under_poisson(med),
            fpr_q25=fpr_under_poisson(float(np.percentile(top, 25))),
            fpr_q75=fpr_under_poisson(float(np.percentile(top, 75))),
            alpha_pool=float(f.alpha_pool[j]),
            n_genes_in_top=int(top.size),
            sample_multiplier=med,
            acceptable=bool(fpr_under_poisson(med) <= fpr_target),
        ))
    return rows


def max_usable_bin(table: list[dict], *, fpr_target: float = 0.10) -> float | None:
    """Largest bin size at which the median gene still passes the target."""
    ok = [r["bin_um"] for r in table if r["acceptable"]]
    return max(ok) if ok else None


# --------------------------------------------------------------------------
# the advice, in sentences a reader can paste into a methods section
# --------------------------------------------------------------------------


def verdict(m: Moments, f: Fit, *, fpr_target: float = 0.10,
            top_n: int = PAPER_TOP_N) -> dict:
    table = scale_table(m, f, fpr_target=fpr_target, top_n=top_n)
    if not table:
        return dict(available=False, reason="no usable scale")

    finest, coarsest = table[0], table[-1]
    usable = max_usable_bin(table, fpr_target=fpr_target)

    lines: list[str] = []
    lines.append(
        f"On this section a Poisson test at {coarsest['bin_um']:.0f} um would "
        f"reject {coarsest['fpr_median']*100:.0f}% of true nulls (nominal 5%), "
        f"against {finest['fpr_median']*100:.0f}% at {finest['bin_um']:.0f} um; "
        f"this is the median over the {coarsest['n_genes_in_top']:,} most "
        f"abundant genes.")

    if usable is None:
        lines.append(
            f"No bin size on this ladder keeps the median gene under "
            f"{fpr_target*100:.0f}% false positives. Use a negative-binomial or "
            f"quasi-likelihood test at every scale, or report the inflation "
            f"alongside the p-value.")
    elif usable >= coarsest["bin_um"]:
        lines.append(
            f"Every bin size on this ladder stays under {fpr_target*100:.0f}% "
            f"false positives for the median gene, so a Poisson test is defensible "
            f"throughout this section.")
    else:
        lines.append(
            f"To keep the median gene under {fpr_target*100:.0f}% false positives, "
            f"do not bin coarser than {usable:.0f} um; beyond that, use a "
            f"negative-binomial or quasi-likelihood test.")

    # the concentration warning belongs in the decision, not only in the report
    try:
        top = float(f.top1_share)
        name = str(m.genes[f.top1])
    except Exception:
        top, name = float("nan"), "?"
    if np.isfinite(top) and top > 0.5:
        lines.append(
            f"The exponent is set by {name} ({top*100:.0f}% of it). Quote that share "
            f"next to the exponent, or the number cannot be compared with another "
            f"tissue.")

    if coarsest["fpr_median"] > 0.2:
        lines.append(
            f"At {coarsest['bin_um']:.0f} um a fixed effect needs "
            f"{coarsest['sample_multiplier']:.0f}x the sample size to reach the "
            f"same power it would have had under Poisson.")

    checklist = [
        "bin ladder: " + ", ".join(f"{r['bin_um']:.0f}" for r in table) + " um",
        f"gene filter: mean count >= {f.mu_min:g} at the largest bin "
        f"({f.n_fitted:,} of {m.n_genes:,} genes fitted)",
        f"platform and binning: state the instrument and the bin edge length",
        f"exponent: {f.gamma_pool:.3f} with its scale range",
        f"variance inflation: phi = 1 + mu*alpha, median over the {top_n} most "
        f"abundant genes, with the bin size it belongs to",
    ]
    if np.isfinite(top) and top > 0.5:
        checklist.append(f"concentration: {name} supplies {top*100:.0f}% of the exponent")

    return dict(
        available=True,
        fpr_target=fpr_target,
        top_n=top_n,
        max_usable_bin_um=usable,
        table=table,
        sentences=lines,
        reporting_checklist=checklist,
        identifiability=_safe_ident(f),
    )


def _safe_ident(f: Fit) -> dict:
    try:
        from .core import identifiability
        return identifiability(f)
    except Exception:
        return dict(defined=False, reason="not computed")


# --------------------------------------------------------------------------
# quasi-likelihood: from "how bad is it" to "what do I do"
# --------------------------------------------------------------------------
#
# The functions above answer what the Poisson assumption costs.  A reader who
# accepts the cost still has to act on it, and the action is one line of
# arithmetic: divide the Poisson z by sqrt(phi).  That is the quasi-likelihood
# score test of Wedderburn (1974), and it is the same correction the manuscript
# reports as returning the false-positive rate to 4.9-5.1%.
#
# Stated plainly so it cannot be misread: this is a variance correction, not a
# model.  It does not create independent replicates.  When the two groups being
# compared are whole regions rather than randomly assigned bins, phi estimated
# from the same bins is anti-conservative and the test stays liberal -- the
# manuscript measures 17-36% rejection of true nulls in that design, and no
# per-bin correction repairs it.


def ql_phi(mu: float, alpha: float) -> float:
    """Variance ratio of a negative-binomial count to its Poisson value.

    ``phi = 1 + alpha * mu`` for counts with ``Var(Y) = mu + alpha * mu**2``.
    This is the quantity the report calls the variance inflation.
    """
    mu = float(mu)
    alpha = float(alpha)
    if mu < 0.0:
        raise ValueError("mu must be non-negative")
    if alpha < 0.0:
        raise ValueError("alpha must be non-negative")
    return 1.0 + alpha * mu


def ql_z(z_poisson: float, phi: float) -> float:
    """Quasi-likelihood z for a Poisson z whose variance is inflated by ``phi``.

    The Poisson statistic divides the residual by ``sqrt(mu)`` where it should
    divide by ``sqrt(phi * mu)``, so the corrected statistic is
    ``z_poisson / sqrt(phi)``.  ``phi = 1`` returns the input unchanged.
    """
    phi = float(phi)
    if phi <= 0.0:
        raise ValueError("phi must be positive")
    return float(z_poisson) / math.sqrt(phi)


def ql_p(z_poisson: float, phi: float, two_sided: bool = True) -> float:
    """P-value of the quasi-likelihood z.  Two-sided by default."""
    z = ql_z(z_poisson, phi)
    p = 2.0 * (1.0 - _phi_cdf(abs(z))) if two_sided else 1.0 - _phi_cdf(z)
    return float(min(1.0, max(0.0, p)))


def ql_test(observed: float, expected: float, mu: float, alpha: float,
            two_sided: bool = True) -> dict:
    """Quasi-likelihood Wald test for one summed count, at one bin size.

    Parameters
    ----------
    observed
        Total counts in the region being tested (the sum over its bins).
    expected
        Total counts expected under the null, i.e. ``n_bins * mu``.
    mu
        Mean count per bin at this bin size.
    alpha
        Overdispersion at this bin size, from the fitted ladder.
    two_sided
        Whether to double the upper tail.

    Returns
    -------
    dict with ``z_poisson``, ``phi``, ``z`` and ``p``.

    Notes
    -----
    ``phi`` is taken as known.  When it is estimated from the same bins, the
    correction is anti-conservative for small numbers of bins, and it does not
    repair a design whose groups are whole regions rather than randomised bins.
    """
    observed = float(observed)
    expected = float(expected)
    if expected <= 0.0:
        raise ValueError("expected must be positive")
    if observed < 0.0:
        raise ValueError("observed must be non-negative")
    phi = ql_phi(mu, alpha)
    z_poisson = (observed - expected) / math.sqrt(expected)
    z = ql_z(z_poisson, phi)
    p = ql_p(z_poisson, phi, two_sided=two_sided)
    return dict(z_poisson=float(z_poisson), phi=float(phi), z=float(z), p=p)

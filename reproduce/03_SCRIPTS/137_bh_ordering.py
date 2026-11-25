"""137_bh_ordering.py -- what BH does to a binned comparison.

The raw false-positive rate is not the number a reader sees; a pipeline reports
a BH q value and a top-100 list.  Two things can happen when the variance is
understated:

  (a) BH holds the *proportion* near q, because it is a proportion, and a
      reader concludes the correction saved the analysis;
  (b) the ranking is contaminated, because a null gene with a large alpha has
      its z inflated more than a real gene with a small alpha, so the top of
      the list fills with nulls and BH happily reports them.

This script measures both, at three scales, for Poisson and for quasi-Poisson
with the alpha read from the eight-scale moments.
"""
from __future__ import annotations

import json
from math import erf, sqrt
from pathlib import Path

import numpy as np

DS = np.array([2, 4, 8, 16, 32, 64, 128, 256])
OUT = Path(__file__).resolve().parent.parent / "out"
DATA = Path(__file__).resolve().parent.parent / "data"
Z = 1.959963984540054
SEED = 20261125
G = 400
N_TRUE = 40
REPS = 100


def cdf(x):
    x = np.asarray(x, float)
    return 0.5 * (1.0 + np.vectorize(lambda t: erf(t / sqrt(2.0)))(x))


def p_two(z):
    return 2.0 * (1.0 - cdf(np.abs(z)))


def bh(p, q=0.05):
    """Benjamini-Hochberg; returns a boolean rejection mask."""
    m = p.size
    order = np.argsort(p)
    ranked = p[order]
    thresh = q * np.arange(1, m + 1) / m
    below = np.where(ranked <= thresh)[0]
    out = np.zeros(m, bool)
    if below.size:
        k = below.max()
        out[order[: k + 1]] = True
    return out


def nb_draw(rng, mu, alpha):
    npar = 1.0 / np.maximum(alpha, 1e-12)
    p = 1.0 / (1.0 + alpha * mu)
    return rng.negative_binomial(npar, p).astype(float)


def load(path):
    d = np.load(path)
    mu = np.stack([d[f"mu_{x}"] for x in DS], 1).astype(float)
    var = np.stack([d[f"var_{x}"] for x in DS], 1).astype(float)
    return mu, var, d["n_cells"].astype(float)


def run_scale(mu, var, n, j, reps=REPS):
    m_all = mu[:, j]
    with np.errstate(invalid="ignore", divide="ignore"):
        a_all = (var[:, j] - m_all) / m_all ** 2
    order = np.argsort(-m_all)
    order = order[m_all[order] > 0][:G]
    keep = order[np.isfinite(a_all[order]) & (a_all[order] > 0)]
    m, a = m_all[keep], a_all[keep]
    genes = keep.size
    if genes < 60:
        return None

    nbin = int(min(n[j], 800))
    rng = np.random.default_rng(SEED + j)

    is_true = np.zeros(genes, bool)
    # true effects are spread over the abundance range, not only the top
    idx = rng.choice(genes, size=min(N_TRUE, genes // 4), replace=False)
    is_true[idx] = True
    delta = np.where(is_true, 1.0, 0.0)
    phi = 1.0 + a * m

    acc = {k: [] for k in ("fdr_p", "fdr_q", "rec_p", "rec_q",
                           "topk_true_p", "topk_true_q",
                           "nrej_p", "nrej_q")}
    for _ in range(reps):
        m1 = m
        m2 = m * np.exp(delta)
        y1 = nb_draw(rng, np.broadcast_to(m1, (nbin, genes)).ravel(),
                     np.broadcast_to(a, (nbin, genes)).ravel()).reshape(nbin, genes)
        y2 = nb_draw(rng, np.broadcast_to(m2, (nbin, genes)).ravel(),
                     np.broadcast_to(a, (nbin, genes)).ravel()).reshape(nbin, genes)
        m1h = np.maximum(y1.mean(0), 1e-9)
        m2h = np.maximum(y2.mean(0), 1e-9)
        se = np.sqrt(1.0 / (nbin * m1h) + 1.0 / (nbin * m2h))
        z = np.log(m2h / m1h) / se
        p_p = p_two(z)
        p_q = p_two(z / np.sqrt(phi))

        rp, rq = bh(p_p), bh(p_q)
        n_t = is_true.sum()
        for tag, rej in (("p", rp), ("q", rq)):
            nrej = int(rej.sum())
            fp = int((rej & ~is_true).sum())
            tp = int((rej & is_true).sum())
            acc[f"fdr_{tag}"].append(fp / nrej if nrej else 0.0)
            acc[f"rec_{tag}"].append(tp / n_t if n_t else 0.0)
            acc[f"nrej_{tag}"].append(nrej)
            # precision at the number of planted signals: a clean ranking puts
            # every true gene in the first n_true positions
            k = int(n_t)
            top = np.argsort(p_p if tag == "p" else p_q)[:k]
            acc[f"topk_true_{tag}"].append(int(is_true[top].sum()) / k)

    return {
        "bin_um": float(DS[j]),
        "n_bin": nbin,
        "n_genes": int(genes),
        "n_true": int(is_true.sum()),
        "phi_median": float(np.median(phi)),
        **{k: float(np.mean(v)) for k, v in acc.items()},
    }


def main():
    out = {"per_section": {}, "config": {"G": G, "N_TRUE": N_TRUE, "REPS": REPS,
                                         "delta_true": 1.0, "q": 0.05}}
    for pre, fn in (("P1", "colon_P1_moments.npz"), ("P2", "colon_P2_moments.npz")):
        mu, var, n = load(DATA / fn)
        rows = []
        for j in (2, 4, 5, 6):     # 8, 16, 32, 128 um
            r = run_scale(mu, var, n, j)
            if r:
                rows.append(r)
        out["per_section"][pre] = rows

    (OUT / "bh_ordering.json").write_text(json.dumps(out, indent=2), encoding="utf-8")

    for pre in out["per_section"]:
        print(f"\n== {pre} ==  (G={G}, {N_TRUE} true, logFC=1, q=0.05)")
        print("   bin   phi   FDR_P  FDR_QL  rec_P  rec_QL  nrej_P nrej_QL  top40-true_P  top40-true_QL")
        for r in out["per_section"][pre]:
            print(f"  {r['bin_um']:4.0f} {r['phi_median']:6.2f} "
                  f"{r['fdr_p']*100:6.1f}% {r['fdr_q']*100:6.1f}% "
                  f"{r['rec_p']*100:6.1f}% {r['rec_q']*100:6.1f}% "
                  f"{r['nrej_p']:7.0f} {r['nrej_q']:7.0f} "
                  f"{r['topk_true_p']*100:13.1f}% {r['topk_true_q']*100:13.1f}%")
    print("\nwrote", OUT / "bh_ordering.json")


if __name__ == "__main__":
    main()

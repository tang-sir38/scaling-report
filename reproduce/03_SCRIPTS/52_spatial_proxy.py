# -*- coding: utf-8 -*-
r"""
52_spatial_proxy.py -- is the per-gene exponent a proxy for spatial clustering?

The obvious alternative explanation for a monotone gamma_g ~ abundance curve is
that gamma_g measures how clustered a transcript is: a housekeeping gene fills the
section evenly, a cell-type marker sits in islands, and coarsening averages the
islands away. Under that reading the exponent is spatial structure in disguise,
and the abundance gradient is incidental.

Two per-gene spatial statistics are computed on the 32 um grid, where counts per
bin are large enough for the statistics to be estimable at all:

  xi_g   exponential range from the semivariogram, fitted on seven lags
  I_g    Moran's I on the eight-neighbour adjacency

and both are related to gamma_g before and after abundance is held fixed
(partial correlation on log10 mean, plus within-decile Spearman).

A range that sits on the lowest admissible value carries no information; the
fraction of genes pinned at that bound is reported so the reader can see whether
xi is estimated or merely bounded. Moran's I does not have that failure mode.

Output: 02_OUT/spatial_proxy/spatial_proxy_<tissue>.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "02_OUT" / "spatial_proxy"
LGD = np.log(np.array([2, 4, 8, 16, 32, 64, 128, 256], float))
STRIDE = 4096

OFFSETS = [(0, 1), (1, 1), (0, 2), (2, 2), (0, 4), (4, 4), (0, 8)]
LAGS = np.array([32, 45.25, 64, 90.51, 128, 181.02, 256.0])
XI_FLOOR = 40.0

SETS = {
    "colon": dict(
        bins=ROOT / "01_WORK" / "bins_032um.npz",
        curves=ROOT / "02_OUT" / "gene_curves" / "gene_curves_colon.npz",
    ),
    "lung": dict(
        bins=ROOT / "D3_LungCancer" / "01_WORK" / "bins_032um.npz",
        curves=ROOT / "02_OUT" / "gene_curves" / "gene_curves_lung.npz",
    ),
    "brain": dict(
        bins=ROOT / "D2_MouseBrain" / "01_WORK" / "bins_032um.npz",
        curves=ROOT / "02_OUT" / "gene_curves" / "gene_curves_brain.npz",
    ),
}


def load_bins(path):
    z = np.load(path, allow_pickle=False)
    shape = tuple(int(x) for x in z["shape"])
    X = csr_matrix((z["data"].astype(np.float32), z["indices"], z["indptr"]), shape=shape)
    return X, z["coords"].astype(np.int64), np.array([str(g) for g in z["genes"]])


def pairs_for(coords, dr, dc):
    key = coords[:, 0] * STRIDE + coords[:, 1]
    order = np.argsort(key)
    skey = key[order]
    tgt = (coords[:, 0] + dr) * STRIDE + (coords[:, 1] + dc)
    pos = np.clip(np.searchsorted(skey, tgt), 0, len(skey) - 1)
    hit = skey[pos] == tgt
    i = np.nonzero(hit)[0]
    j = order[pos[hit]]
    return i, j


def semivariograms(X, coords, block=200):
    prs = [pairs_for(coords, dr, dc) for dr, dc in OFFSETS]
    out = np.full((X.shape[0], len(prs)), np.nan, np.float32)
    for s in range(0, X.shape[0], block):
        e = min(s + block, X.shape[0])
        D = X[s:e].toarray()
        for k, (i, j) in enumerate(prs):
            if i.size < 100:
                continue
            d = D[:, i] - D[:, j]
            out[s:e, k] = 0.5 * (d * d).mean(axis=1)
        del D
    return out, np.array([i.size for i, _ in prs])


def fit_range(SV):
    """Exponential model sv(h) = c0 + c1 (1 - exp(-h/xi)), linear in (c0, c1)."""
    grid = np.geomspace(XI_FLOOR, 20000.0, 60)
    G = SV.shape[0]
    best_sse = np.full(G, np.inf)
    best_xi = np.full(G, np.nan)
    best_beta = np.full((2, G), np.nan)
    y = SV.astype(np.float64)
    ytot = np.nansum(y ** 2, axis=1)
    for xi in grid:
        A = np.column_stack([np.ones_like(LAGS), 1.0 - np.exp(-LAGS / xi)])
        AtAinv = np.linalg.inv(A.T @ A)
        beta = AtAinv @ (A.T @ y.T)
        sse = np.nansum((y.T - A @ beta) ** 2, axis=0)
        better = sse < best_sse
        best_sse[better] = sse[better]
        best_xi[better] = xi
        best_beta[:, better] = beta[:, better]
    return best_xi, best_beta.T, best_sse / np.maximum(ytot, 1e-30)


def morans_i(X, coords, block=200):
    """Moran's I on the symmetric union of the 7 lag-1..8 offsets."""
    prs = [pairs_for(coords, dr, dc) for dr, dc in OFFSETS]
    n = X.shape[1]
    rows, cols = [], []
    for i, j in prs:
        rows.append(i)
        cols.append(j)
        rows.append(j)
        cols.append(i)
    r = np.concatenate(rows)
    c = np.concatenate(cols)
    keep = r != c
    W = csr_matrix((np.ones(keep.sum(), np.float32), (r[keep], c[keep])), shape=(n, n))
    wsum = float(W.sum())
    out = np.full(X.shape[0], np.nan)
    for s in range(0, X.shape[0], block):
        e = min(s + block, X.shape[0])
        D = X[s:e].toarray().astype(np.float64)
        mean = D.mean(axis=1, keepdims=True)
        Z = D - mean
        den = (Z ** 2).sum(axis=1)
        num = np.einsum("ij,ij->i", (W @ Z.T).T, Z)
        with np.errstate(invalid="ignore", divide="ignore"):
            out[s:e] = (n / wsum) * num / den
        del D, Z
    return out, wsum, W.nnz


def resid(a, b):
    A = np.column_stack([np.ones_like(b), b])
    coef, *_ = np.linalg.lstsq(A, a, rcond=None)
    return a - A @ coef


def pcorr(a, b, ctrl):
    return float(np.corrcoef(resid(a, ctrl), resid(b, ctrl))[0, 1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue", default="colon", choices=list(SETS))
    args = ap.parse_args()
    tis = args.tissue
    OUT.mkdir(parents=True, exist_ok=True)

    X, coords, genes = load_bins(SETS[tis]["bins"])
    print(f"{tis}: {X.shape[0]} genes x {X.shape[1]} bins", flush=True)

    cur = np.load(SETS[tis]["curves"], allow_pickle=True)
    cg = np.array([str(g) for g in cur["gene"]])
    cmap = {g: i for i, g in enumerate(cg)}
    idx = np.array([cmap.get(g, -1) for g in genes])
    gamma = np.full(len(genes), np.nan)
    mu32 = np.full(len(genes), np.nan)
    ok_c = idx >= 0
    gamma[ok_c] = cur["gamma"][idx[ok_c]].astype(np.float64)
    mu32[ok_c] = cur["mu32"][idx[ok_c]].astype(np.float64)

    sub = np.nonzero(np.isfinite(gamma) & (gamma > 0))[0]
    print(f"  {len(sub)} genes with a computable exponent", flush=True)
    Xs = X[sub]

    SV, npair = semivariograms(Xs, coords)
    print("  pair counts:", npair, flush=True)
    xi_s, beta_s, rel_s = fit_range(SV)

    I_s, wsum, wnnz = morans_i(Xs, coords)
    print(f"  Moran adjacency: {wnnz:,} directed entries, sum w = {wsum:,.0f}", flush=True)

    g = gamma[sub]
    m = mu32[sub]
    good = np.isfinite(m) & (m > 0) & np.isfinite(xi_s) & np.isfinite(I_s)
    lm = np.log10(m[good])
    gg, xx, ii, rr = g[good], np.log10(np.maximum(xi_s[good], 1)), I_s[good], rel_s[good]
    structured = (beta_s[good, 1] > 0) & (rr < 0.5)

    def sp(a, b):
        r = spearmanr(a, b)
        return float(r.statistic), float(r.pvalue)

    r_xi, p_xi = sp(gg, xx)
    r_I, p_I = sp(gg, ii)
    r_mu, _ = sp(gg, lm)

    q = np.quantile(lm, np.linspace(0, 1, 11))
    dec = []
    for k in range(10):
        sel = (lm >= q[k]) & ((lm <= q[k + 1]) if k == 9 else (lm < q[k + 1]))
        if sel.sum() > 200:
            dec.append(dict(decile=k + 1, n=int(sel.sum()),
                            rho_xi=sp(gg[sel], xx[sel])[0],
                            rho_moran=sp(gg[sel], ii[sel])[0]))

    # ---- stratified and semiparametric control for abundance -------------
    # A log-linear partial correlation is not adequate here: both statistics are
    # non-linear in log abundance. Two assumption-light alternatives are reported.
    def band_rho(k):
        e = np.quantile(lm, np.linspace(0, 1, k + 1))
        vals = []
        for i in range(k):
            sel = (lm >= e[i]) & ((lm <= e[i + 1]) if i == k - 1 else (lm < e[i + 1]))
            if sel.sum() > 60:
                vals.append(sp(gg[sel], ii[sel])[0])
        v = np.array(vals)
        return dict(k=k, n_bands=int(v.size), median=float(np.median(v)),
                    lo=float(v.min()), hi=float(v.max()))

    def smooth_resid(y, x, k=40):
        e = np.quantile(x, np.linspace(0, 1, k + 1))
        c = np.zeros_like(y, dtype=float)
        for i in range(k):
            sel = (x >= e[i]) & ((x <= e[i + 1]) if i == k - 1 else (x < e[i + 1]))
            if sel.any():
                c[sel] = y[sel].mean()
        return y - c

    bands = [band_rho(k) for k in (10, 20, 40, 80)]
    rg, ri = smooth_resid(gg, lm), smooth_resid(ii, lm)
    semipar = dict(pearson=float(np.corrcoef(rg, ri)[0, 1]), spearman=sp(rg, ri)[0], k=40)

    res = dict(
        tissue=tis,
        n_ok=int(good.sum()),
        n_structured=int(structured.sum()),
        xi_median_um=float(np.median(xi_s[good][structured])) if structured.any() else None,
        xi_iqr_um=[float(np.percentile(xi_s[good][structured], 25)),
                   float(np.percentile(xi_s[good][structured], 75))] if structured.any() else None,
        xi_frac_at_floor=float(np.mean(xi_s[good] <= XI_FLOOR * 1.10)),
        xi_floor_um=XI_FLOOR,
        moran_median=float(np.median(ii)),
        moran_iqr=[float(np.percentile(ii, 25)), float(np.percentile(ii, 75))],
        rho_gamma_mu=r_mu,
        rho_gamma_xi=r_xi, p_gamma_xi=p_xi,
        rho_gamma_moran=r_I, p_gamma_moran=p_I,
        partial_rho_xi=pcorr(gg, xx, lm),
        partial_rho_moran=pcorr(gg, ii, lm),
        band_rho=bands,
        semiparam=semipar,
        by_decile=dec,
        lags_um=LAGS.tolist(),
        pair_counts=npair.tolist(),
    )

    np.savez_compressed(OUT / f"per_gene_{tis}.npz",
                        gene=genes[sub], gamma=g, mu32=m,
                        xi=xi_s, moran=I_s, sv=SV, rel=rel_s)

    print("\n=== result ===")
    print(f"  rho(gamma, log mu)     {r_mu:+.3f}")
    print(f"  rho(gamma, xi)         {r_xi:+.3f}   partial {res['partial_rho_xi']:+.3f}")
    print(f"  rho(gamma, Moran I)    {r_I:+.3f}   partial {res['partial_rho_moran']:+.3f}")
    print(f"  xi median {res['xi_median_um']} um, at floor {res['xi_frac_at_floor']*100:.0f}%")
    print("  within-decile rho: xi " + " ".join(f"{d['rho_xi']:+.2f}" for d in dec))
    print("                    I  " + " ".join(f"{d['rho_moran']:+.2f}" for d in dec))
    for b in bands:
        print(f"  rho(gamma, I) within {b['k']:>2} abundance bands: median {b['median']:+.3f}"
              f"  ({b['lo']:+.2f} to {b['hi']:+.2f})")
    print(f"  semiparametric partial rho: Pearson {semipar['pearson']:+.3f}"
          f"  Spearman {semipar['spearman']:+.3f}")

    p = OUT / f"spatial_proxy_{tis}.json"
    p.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print("\nwritten", p)


if __name__ == "__main__":
    main()

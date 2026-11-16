# -*- coding: utf-8 -*-
r"""
51_robust_scalar.py -- which summary scalar survives a change of region?

The manuscript argues that a single pooled exponent does not travel. That is a
criticism; this script supplies the constructive half: of the scalars one could
report instead, which is stable across regions of one section?

Five candidates are computed inside each of the four column bands and the two
halves of the colon section, using the same cut rule as 37_splithalf.py
(multiples of 128 array units, so no 256 um block is divided) and the same gene
filter as the manuscript:

  pooled      alpha_pool(d) = sum(var-mu)/sum(mu^2), slope over 8 scales
  median      median of per-gene gamma
  mean        mean of per-gene gamma
  trimmed     mean after dropping the outer 5% on each side
  wmedian     mu^2-weighted median of per-gene gamma

Stability is the SD, range and CV across the four bands. A scalar that is to be
compared between samples has to be stable within one.

Output: 02_OUT/robust_scalar/robust_scalar.json
"""
from __future__ import annotations

import gc
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "02_OUT" / "robust_scalar"
NPZ = ROOT / "01_WORK" / "bins_002um.npz"
DS = [2, 4, 8, 16, 32, 64, 128, 256]
STRIDE = 4096
LGD = np.log(np.array(DS, float))
SCALARS = ("pooled", "median", "mean", "trimmed", "wmedian")


def load(p):
    z = np.load(p, allow_pickle=False)
    shape = tuple(int(x) for x in z["shape"])
    X = csr_matrix((z["data"].astype(np.float32), z["indices"], z["indptr"]), shape=shape)
    return X, z["coords"].astype(np.int64), np.array([str(g) for g in z["genes"]]), z


def moments(X, n):
    s1 = np.asarray(X.sum(axis=1)).ravel().astype(np.float64)
    s2 = np.asarray(X.multiply(X).sum(axis=1)).ravel().astype(np.float64)
    mu = s1 / n
    var = (s2 - s1 ** 2 / n) / (n - 1) if n > 1 else np.zeros_like(mu)
    return mu, var


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    X, coords, genes, z = load(NPZ)
    gs = np.array([str(g) for g in genes])
    keep = ~np.array([s.upper().startswith("MT-") for s in gs])
    X = X[keep].tocsr()
    gs = gs[keep]
    R = coords[:, 0].astype(np.int64)
    C = coords[:, 1].astype(np.int64)
    dom_any = z["dom_any"]
    Rd = (dom_any // STRIDE).astype(np.int64)
    Cd = (dom_any % STRIDE).astype(np.int64)
    print(f"genes={X.shape[0]}  bins={X.shape[1]}  columns {C.min()}..{C.max()}", flush=True)

    AGG, DOMKEY = {}, {}
    for d in DS:
        f = d // 2
        key = (R // f) * 100000 + (C // f)
        uk, inv = np.unique(key, return_inverse=True)
        AGG[d] = csr_matrix((np.ones(key.size, np.float32),
                             (np.arange(key.size), inv.astype(np.int32))),
                            shape=(key.size, uk.size))
        DOMKEY[d] = (Rd // f) * 100000 + (Cd // f)
        del uk, inv
        gc.collect()
    print("aggregation operators ready", flush=True)

    def run(mb, md):
        Xb = X[:, mb]
        G = X.shape[0]
        AL = np.full((G, len(DS)), np.nan)
        MU = np.full((G, len(DS)), np.nan)
        AP = np.empty(len(DS))
        for j, d in enumerate(DS):
            Xa = Xb if d == 2 else Xb @ AGG[d][mb]
            n = int(np.unique(DOMKEY[d][md]).size)
            mu, var = moments(Xa, n)
            num = var - mu
            mu2 = mu ** 2
            with np.errstate(invalid="ignore", divide="ignore"):
                AL[:, j] = num / mu2
            MU[:, j] = mu
            AP[j] = num.sum() / mu2.sum()
            if d != 2:
                del Xa
        del Xb
        gc.collect()

        ok = np.all(np.isfinite(AL) & (AL > 0), axis=1) & (MU[:, -1] >= 0.5)
        y = np.log(AL[ok])
        xm = LGD - LGD.mean()
        gammas = -np.sum(xm[None, :] * (y - y.mean(1, keepdims=True)), axis=1) / np.sum(xm ** 2)
        gp = -np.polyfit(LGD, np.log(AP), 1)[0]

        w = MU[ok, 4] ** 2                       # mu^2 at 32 um, the abundance weight
        o = np.argsort(gammas)
        cw = np.cumsum(w[o])
        cw /= cw[-1]
        wmed = float(gammas[o][np.searchsorted(cw, 0.5)])
        lo, hi = np.percentile(gammas, [5, 95])
        return dict(
            pooled=float(gp),
            median=float(np.median(gammas)),
            mean=float(gammas.mean()),
            trimmed=float(gammas[(gammas >= lo) & (gammas <= hi)].mean()),
            wmedian=wmed,
            n_fit=int(ok.sum()),
            alpha=[float(a) for a in AP],
        )

    full = run(np.ones(X.shape[1], bool), np.ones(dom_any.size, bool))
    print("whole section:", {k: round(full[k], 4) for k in SCALARS}, flush=True)

    c_lo, c_hi = int(C.min()), int(C.max())
    cut = int(np.floor(((c_lo + c_hi) / 2) / 128) * 128)
    edges = np.floor(np.linspace(c_lo, c_hi + 1, 5) / 128) * 128
    # Keep both outer edges on the 128 grid, exactly as 37_splithalf.py does: the
    # residual strip beyond the last multiple is dropped so that no 256 um block is
    # divided and the two scripts describe the same four bands. (Overriding the last
    # edge to c_hi + 1 silently gave a different fourth band and a different SD.)
    print("half cut:", cut, " band edges:", edges.tolist(), flush=True)

    parts, regions = [], []
    for nm, mb, md in (("halfA", C < cut, Cd < cut), ("halfB", C >= cut, Cd >= cut)):
        r = run(mb, md)
        r.update(name=nm, n_bin=int(mb.sum()))
        parts.append(r)
        regions.append(("half", nm, r))
        print(f"  {nm}: " + "  ".join(f"{k}={r[k]:.4f}" for k in SCALARS), flush=True)

    for k in range(4):
        lo, hi = int(edges[k]), int(edges[k + 1])
        mb = (C >= lo) & (C < hi)
        md = (Cd >= lo) & (Cd < hi)
        if mb.sum() < 1000 or np.unique(DOMKEY[256][md]).size < 20:
            print(f"  band {k}: too small, skipped", flush=True)
            continue
        r = run(mb, md)
        r.update(name=f"q{k + 1}", n_bin=int(mb.sum()))
        parts.append(r)
        print(f"  q{k + 1}: " + "  ".join(f"{k2}={r[k2]:.4f}" for k2 in SCALARS), flush=True)

    bands = [p for p in parts if p["name"].startswith("q")]
    summary = {}
    for key in SCALARS:
        v = np.array([b[key] for b in bands])
        summary[key] = dict(
            values=v.tolist(),
            sd=float(v.std(ddof=1)),
            range=float(v.max() - v.min()),
            cv=float(v.std(ddof=1) / v.mean()),
        )

    best = min(SCALARS, key=lambda k: summary[k]["sd"])
    worst = max(SCALARS, key=lambda k: summary[k]["sd"])
    res = dict(
        section="colon",
        scalars=list(SCALARS),
        full={k: full[k] for k in SCALARS},
        n_fit_full=full["n_fit"],
        bands=parts,
        half_cut=cut,
        band_edges=[int(e) for e in edges],
        summary=summary,
        best=dict(name=best, sd=summary[best]["sd"], cv=summary[best]["cv"]),
        worst=dict(name=worst, sd=summary[worst]["sd"], cv=summary[worst]["cv"]),
        ratio_worst_over_best=float(summary[worst]["sd"] / summary[best]["sd"]),
        ratio_pooled_over_best=float(summary["pooled"]["sd"] / summary[best]["sd"]),
    )

    print("\n=== stability across the four bands ===")
    for k in SCALARS:
        s = summary[k]
        print(f"  {k:<9} SD={s['sd']:.4f}  range={s['range']:.4f}  CV={s['cv'] * 100:.1f}%")
    print(f"  best={best}  worst={worst}  worst/best={res['ratio_worst_over_best']:.2f}")

    p = OUT / "robust_scalar.json"
    p.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print("\nwritten", p)


if __name__ == "__main__":
    main()

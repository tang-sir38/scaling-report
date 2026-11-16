"""13_gene_moments.py — 导出每个基因、每个分辨率的 (mu, var)，供后续任意分层复用。

只需跑一次，之后任何"按表达量/按基因类别"的重新分层都不必再读大矩阵。
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix


def moments(x, n_cells):
    s1 = np.asarray(x.sum(axis=1)).ravel().astype(np.float64)
    s2 = np.asarray(x.multiply(x).sum(axis=1)).ravel().astype(np.float64)
    mu = s1 / n_cells
    var = (s2 - s1 ** 2 / n_cells) / (n_cells - 1)
    return mu, var


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--ds", nargs="+", type=int, default=[2, 4, 8, 16, 32, 64, 128, 256])
    ap.add_argument("--domain", default="any")
    ap.add_argument("--exclude-mt", action="store_true")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()

    wd, od = Path(a.workdir), Path(a.outdir)
    od.mkdir(parents=True, exist_ok=True)
    pack = {"ds": a.ds, "mu": {}, "var": {}, "n_cells": {}, "genes": None}

    for d in a.ds:
        p = wd / f"bins_{d:03d}um.npz"
        if not p.exists():
            print(f"  缺 {p.name}", flush=True)
            continue
        z = np.load(p, allow_pickle=True)
        X = csr_matrix((z["data"], z["indices"], z["indptr"]),
                       shape=tuple(int(v) for v in z["shape"]))
        key = f"n_domain_{a.domain}"
        n = int(z[key]) if key in z.files else int(X.shape[1])
        genes = z["genes"]
        keep = np.ones(len(genes), bool)
        if a.exclude_mt:
            gs = [g.decode() if isinstance(g, bytes) else str(g) for g in genes]
            keep = ~np.array([s.upper().startswith("MT-") for s in gs])
        mu, var = moments(X, n)
        pack["mu"][d] = mu[keep]
        pack["var"][d] = var[keep]
        pack["n_cells"][d] = n
        if pack["genes"] is None:
            pack["genes"] = np.array([g.decode() if isinstance(g, bytes) else str(g)
                                      for g in genes])[keep]
        print(f"  d={d:3d}  n_cells={n:>9,}  基因={int(keep.sum())}  "
              f"μ 中位={np.median(mu[keep]):.4g}", flush=True)
        del X, z

    suffix = f"_{a.tag}" if a.tag else ""
    fn = od / f"gene_moments{suffix}.npz"
    np.savez_compressed(fn,
                        ds=np.array(a.ds),
                        genes=pack["genes"],
                        n_cells=np.array([pack["n_cells"].get(d, 0) for d in a.ds]),
                        **{f"mu_{d}": pack["mu"][d] for d in pack["mu"]},
                        **{f"var_{d}": pack["var"][d] for d in pack["var"]})
    print(f"已写出 {fn}  ({fn.stat().st_size/1e6:.1f} MB)")


if __name__ == "__main__":
    main()

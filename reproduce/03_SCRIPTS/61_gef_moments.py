# -*- coding: utf-8 -*-
r"""
61_gef_moments.py -- per-gene moments for a Stereo-seq GEF, on the Visium HD ladder.

Produces the same gene_moments npz that 13_gene_moments.py produces for Visium HD,
so every downstream script (56_section_metrics, 58_pool_sections) runs unchanged.

Contract (must match 13_gene_moments.py exactly):
    mu_g(d)  = S1_g(d) / n(d)
    var_g(d) = (S2_g(d) - S1_g(d)^2 / n(d)) / (n(d) - 1)
with n(d) the number of domain bins at scale d.

Domain: a Stereo-seq `.tissue.gef` holds tissue-filtered records only, so the
domain is the set of bins that contain at least one record -- the direct analogue
of the `in_tissue` "any" convention used for Visium HD.

Scale mapping: bin1 = 1 DNB = 0.5 um, so a d um bin is 2d DNBs on a side and the
ladder 2,4,...,256 um is exact integer aggregation, identical to Visium HD.

Usage:
  python 61_gef_moments.py --bin 00_RAW/stereo/C04042E2_bin1.bin \
      --genes 00_RAW/stereo/C04042E2_genes.npz --out 02_OUT/stereo/gene_moments_stereo.npz \
      --ds 2 4 8 16 32 64 128 256
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

DT = np.dtype([("x", "<i4"), ("y", "<i4"), ("count", "u1")])


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_genes(path: Path) -> dict:
    z = np.load(path, allow_pickle=True)
    return dict(name=z["name"], gid=z["gid"], offset=z["offset"].astype(np.int64),
                count=z["count"].astype(np.int64))


def coord_range(rec, chunk=50_000_000):
    x0 = y0 = 1 << 30
    x1 = y1 = -1
    for a in range(0, len(rec), chunk):
        b = min(a + chunk, len(rec))
        xs, ys = rec["x"][a:b], rec["y"][a:b]
        x0, x1 = min(x0, int(xs.min())), max(x1, int(xs.max()))
        y0, y1 = min(y0, int(ys.min())), max(y1, int(ys.max()))
    return x0, x1, y0, y1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ds", nargs="+", type=int, default=[2, 4, 8, 16, 32, 64, 128, 256])
    ap.add_argument("--exclude-mt", action="store_true")
    args = ap.parse_args()

    rec = np.memmap(args.bin, dtype=DT, mode="r")
    g = load_genes(Path(args.genes))
    assert g["count"].sum() == len(rec), (g["count"].sum(), len(rec))
    log(f"records {len(rec):,}  genes {len(g['name']):,}")

    x0, x1, y0, y1 = coord_range(rec)
    log(f"x {x0}..{x1}  y {y0}..{y1}  (DNB units, 0.5 um)")

    keep = np.ones(len(g["name"]), bool)
    if args.exclude_mt:
        up = np.char.upper(g["name"].astype(str))
        keep = ~np.char.startswith(up, "MT-")

    out = {"ds": np.array(args.ds), "genes": g["name"][keep]}
    n_cells = []

    for d in args.ds:
        t0 = time.time()
        k = 2 * d                                   # DNBs per side
        W = (x1 - x0) // k + 1
        H = (y1 - y0) // k + 1
        nbin = W * H
        mu = np.full(len(g["name"]), np.nan)
        var = np.full(len(g["name"]), np.nan)
        occ = np.zeros(nbin, bool)

        # occupancy: 250M records -> chunk it, never materialise a full-array copy
        for ca in range(0, len(rec), 20_000_000):
            cb = min(ca + 20_000_000, len(rec))
            cx = (rec["x"][ca:cb].astype(np.int64) - x0) // k
            cy = (rec["y"][ca:cb].astype(np.int64) - y0) // k
            occ[cx * H + cy] = True
            del cx, cy
        n = int(occ.sum())
        del occ
        for gi in range(len(g["name"])):
            a = g["offset"][gi]
            b = a + g["count"][gi]
            if b <= a:
                continue
            xs = (rec["x"][a:b] - x0) // k
            ys = (rec["y"][a:b] - y0) // k
            key = xs.astype(np.int64) * H + ys
            c = rec["count"][a:b].astype(np.float64)
            uniq, inv = np.unique(key, return_inverse=True)
            cb = np.bincount(inv, weights=c, minlength=len(uniq))
            s1 = cb.sum()
            s2 = (cb * cb).sum()
            mu[gi] = s1 / n
            var[gi] = (s2 - s1 * s1 / n) / (n - 1)

        out[f"mu_{d}"] = mu[keep]
        out[f"var_{d}"] = var[keep]
        n_cells.append(n)
        log(f"d={d:3d} um  n_domain={n:>9,}  median mu={np.median(mu[keep]):.4g}  "
            f"({time.time()-t0:.0f}s)")

    out["n_cells"] = np.array(n_cells)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.out, **out)
    log(f"wrote {args.out}  ({Path(args.out).stat().st_size/1e6:.1f} MB)")


if __name__ == "__main__":
    main()

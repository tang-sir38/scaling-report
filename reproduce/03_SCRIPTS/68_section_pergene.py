# -*- coding: utf-8 -*-
r"""
68_section_pergene.py -- why does the per-gene median move between two sections
of the same tissue when the pooled exponent does not?

P0 (Visium_HD_Human_Colon_Cancer) and P1 (…_P1) give gamma_pool 0.5547 and
0.5485 -- a difference of 0.006 -- but per-gene medians of 0.4285 and 0.6184, a
difference of 0.190. Two candidate explanations, both testable with the files
already on disk:

  (a) depth. P1 has roughly a third of the UMI of P0. At low depth the 2 um
      overdispersion of a gene is dominated by sampling noise, which inflates
      its slope, and the effect is strongest for low-expression genes. If this
      is the mechanism, the shift must be concentrated in the low-abundance
      deciles and vanish at the top.

  (b) composition. The two sections are different tissue, so the shift could be
      uniform across abundance, i.e. a real change in the per-gene distribution.

The script reports the decile-wise median exponent in each section, the shift,
and a bootstrap CI for the shift in each decile. It does not decide which
explanation is right by assertion: the shape of the shift does.

Usage:
  python 68_section_pergene.py
Output:
  02_OUT/multisection/pergene_shift.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
O = ROOT / "02_OUT"
P1DIR = ROOT / "MS_Visium_HD_Human_Colon_Cancer_P1" / "02_OUT"
DS = np.array([2, 4, 8, 16, 32, 64, 128, 256], float)
LGD = np.log(DS)

FILES = {"colon_P0": O / "gene_moments.npz", "colon_P1": P1DIR / "gene_moments_noMT.npz"}
N_DEC = 10
N_BOOT = 400
RNG = np.random.default_rng(20261115)


def slope(y):
    return -float(np.polyfit(LGD, np.log(y), 1)[0])


def per_gene(path):
    z = np.load(path, allow_pickle=True)
    mu = np.column_stack([z[f"mu_{int(d)}"] for d in DS]).astype(np.float64)
    var = np.column_stack([z[f"var_{int(d)}"] for d in DS]).astype(np.float64)
    genes = np.array([str(g) for g in z["genes"]])
    keep = mu.sum(axis=1) > 0
    mu, var, genes = mu[keep], var[keep], genes[keep]
    with np.errstate(invalid="ignore", divide="ignore"):
        alpha_g = (var - mu) / mu ** 2
    ok = np.all(np.isfinite(alpha_g) & (alpha_g > 0), axis=1) & (mu[:, -1] >= 0.5)
    g = np.full(mu.shape[0], np.nan)
    for i in np.where(ok)[0]:
        g[i] = slope(alpha_g[i])
    return genes, mu[:, -1], g, ok


def boot_median_ci(x, B=N_BOOT):
    if len(x) < 3:
        return None, None
    idx = RNG.integers(0, len(x), size=(B, len(x)))
    m = np.median(x[idx], axis=1)
    return float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


res = {}
tables = {}
for tag, p in FILES.items():
    genes, mu256, g, ok = per_gene(p)
    tables[tag] = dict(genes=genes, mu256=mu256, gamma=g, ok=ok)
    v = g[ok]
    lo, hi = boot_median_ci(v)
    res[tag] = dict(
        n_genes=int(len(genes)), n_fit=int(ok.sum()),
        median=float(np.median(v)), q25=float(np.percentile(v, 25)),
        q75=float(np.percentile(v, 75)),
        median_ci=[lo, hi],
        umi_total=float(mu256.sum()),
        mu256_median=float(np.median(mu256[ok])),
    )

# abundance deciles on the P0 section, then the same edges applied to both
_, mu0, g0, ok0 = (tables["colon_P0"][k] for k in ("genes", "mu256", "gamma", "ok"))
_, mu1, g1, ok1 = (tables["colon_P1"][k] for k in ("genes", "mu256", "gamma", "ok"))
edges = np.quantile(mu0[ok0], np.linspace(0, 1, N_DEC + 1))
edges[0], edges[-1] = 0.0, np.inf   # counts are non-negative

dec = []
for i in range(N_DEC):
    m0 = ok0 & (mu0 >= edges[i]) & (mu0 < edges[i + 1])
    m1 = ok1 & (mu1 >= edges[i]) & (mu1 < edges[i + 1])
    if m0.sum() < 20 or m1.sum() < 20:
        continue
    d = float(np.median(g1[m1]) - np.median(g0[m0]))
    # bootstrap the difference by resampling genes within each decile
    a = g0[m0]; b = g1[m1]
    da = np.median(a[RNG.integers(0, len(a), (N_BOOT, len(a)))], axis=1)
    db = np.median(b[RNG.integers(0, len(b), (N_BOOT, len(b)))], axis=1)
    dd = db - da
    dec.append(dict(
        decile=i + 1,
        mu_lo=float(edges[i]) if np.isfinite(edges[i]) else None,
        mu_hi=float(edges[i + 1]) if np.isfinite(edges[i + 1]) else None,
        n0=int(m0.sum()), n1=int(m1.sum()),
        med0=float(np.median(g0[m0])), med1=float(np.median(g1[m1])),
        shift=d,
        shift_lo=float(np.percentile(dd, 2.5)),
        shift_hi=float(np.percentile(dd, 97.5)),
    ))

shifts = [d["shift"] for d in dec]
lo_dec = [d["shift_lo"] for d in dec]
hi_dec = [d["shift_hi"] for d in dec]

out = dict(
    sections=res,
    deciles=dec,
    shift_overall=float(res["colon_P1"]["median"] - res["colon_P0"]["median"]),
    shift_top_decile=dec[-1]["shift"] if dec else None,
    shift_bottom_decile=dec[0]["shift"] if dec else None,
    n_deciles_ci_excludes_zero=int(sum(1 for l in lo_dec if l > 0)),
    monotone_decreasing=bool(all(shifts[i] >= shifts[i + 1] for i in range(len(shifts) - 1))),
    verdict=("depth-like: shift concentrated at low abundance"
             if dec and dec[0]["shift"] > 2 * abs(dec[-1]["shift"]) else
             "uniform: shift not explained by depth alone"),
    source="68_section_pergene.py",
)

(O / "multisection" / "pergene_shift.json").write_text(
    json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")

print(f"{'dec':>3} {'n0':>6} {'n1':>6} {'P0 med':>7} {'P1 med':>7} {'shift':>7}  95% CI")
for d in dec:
    print(f"{d['decile']:>3} {d['n0']:>6,} {d['n1']:>6,} {d['med0']:>7.3f} "
          f"{d['med1']:>7.3f} {d['shift']:>+7.3f}  [{d['shift_lo']:+.3f}, {d['shift_hi']:+.3f}]")
print("overall shift", round(out["shift_overall"], 4))
print("verdict:", out["verdict"])

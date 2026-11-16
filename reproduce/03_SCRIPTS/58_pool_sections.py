# -*- coding: utf-8 -*-
r"""
58_pool_sections.py -- pool the per-section metrics into the cross-section test.

The manuscript's weakest joint is that every tissue is represented by one section,
so a colon-brain difference cannot be separated from a section-to-section difference.
This script answers that with the sections available:

  1. one row per section: pooled gamma, per-gene median gamma, effective-gene counts,
     the head transcript and its shares, and the within-section spread over four bands;
  2. variance components inside a tissue (one-way random effects over sections,
     bands as replicates): MSW, MSB, sigma_between^2, ICC;
  3. the quantity the claim needs: the colon-brain gap in units of the
     between-section SD;
  4. the cross-section regression of pooled gamma on the head transcript's share of
     the exponent, with a permutation p-value (n is small, so no asymptotics);
  5. minimum detectable difference for n sections per tissue,
     MDE = exp(2.802 * sigma_b * sqrt(2/n)), which is the design statement the
     manuscript can make once sigma_b is measured.

Usage:
  python 58_pool_sections.py
  python 58_pool_sections.py --metrics 02_OUT/section_metrics/*.json MS_*/02_OUT/metrics_*.json
"""
from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "02_OUT" / "multisection"

# label -> tissue. New 10x colon sections arrive as P1/P2/P5; the atlas samples
# keep the names used in the manuscript.
TISSUE = {"colon": "colon_cancer", "lung": "lung_cancer", "brain": "mouse_brain"}


def infer_tissue(label: str) -> str:
    low = label.lower()
    if low in TISSUE:
        return TISSUE[low]
    if low.startswith("p") and low[1:].isdigit():
        return "colon_cancer"
    if low.startswith("colon"):
        return "colon_cancer"
    # Visium_HD_Human_Colon_Cancer_P1 / _P5 ... are additional colon sections
    if "colon" in low and "cancer" in low:
        return "colon_cancer"
    if "lung" in low and "cancer" in low:
        return "lung_cancer"
    if "mouse_brain" in low or "brain" in low:
        return "mouse_brain"
    return low


SH = {"colon": ROOT / "02_OUT" / "splithalf" / "splithalf.json",
      "lung": ROOT / "D3_LungCancer" / "02_OUT" / "splithalf" / "splithalf.json",
      "brain": ROOT / "D2_MouseBrain" / "02_OUT" / "splithalf" / "splithalf.json"}


def find_quarters(label: str) -> float | None:
    """Within-section SD of pooled gamma over four equal-count bands."""
    p = SH.get(label.lower())
    if p is not None and p.exists():
        return float(json.loads(p.read_text(encoding="utf-8"))["quarter_sd"])
    for q in sorted(ROOT.glob(f"MS_*{label}*/02_OUT/splithalf/splithalf.json")):
        return float(json.loads(q.read_text(encoding="utf-8"))["quarter_sd"])
    return None


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    x = x - x.mean()
    y = y - y.mean()
    d = np.sqrt((x ** 2).sum() * (y ** 2).sum())
    return float((x * y).sum() / d) if d > 0 else float("nan")


def perm_p(x: np.ndarray, y: np.ndarray, n_perm: int = 20000, seed: int = 20261114) -> float:
    """Two-sided permutation p for Pearson r: exact when n! is small, else Monte Carlo."""
    import itertools
    if len(x) <= 8:
        obs = abs(pearson(x, y))
        tot = hit = 0
        for perm in itertools.permutations(range(len(x))):
            tot += 1
            if abs(pearson(x, y[list(perm)])) >= obs - 1e-12:
                hit += 1
        return hit / tot
    rng = np.random.default_rng(seed)
    obs = abs(pearson(x, y))
    cnt = 0
    for _ in range(n_perm):
        if abs(pearson(x, rng.permutation(y))) >= obs - 1e-15:
            cnt += 1
    return (cnt + 1) / (n_perm + 1)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--metrics", nargs="*", default=None)
    ap.add_argument("--out", default=str(OUT / "pool.json"))
    a = ap.parse_args()

    paths = a.metrics or (sorted(glob.glob(str(ROOT / "02_OUT" / "section_metrics" / "*.json")))
                          + sorted(glob.glob(str(ROOT / "MS_*" / "02_OUT" / "metrics_*.json"))))
    rows = []
    for p in paths:
        d = json.loads(Path(p).read_text(encoding="utf-8"))
        d["tissue"] = infer_tissue(d["label"])
        d["quarter_sd"] = find_quarters(d["label"])
        rows.append(d)
    if not rows:
        raise SystemExit("no section metrics found")

    rows.sort(key=lambda d: (d["tissue"], d["label"]))
    print(f"{'section':<10}{'tissue':<14}{'gamma_pool':>11}{'gamma_med':>11}"
          f"{'n_eff(alpha)':>13}{'top1':>10}{'share':>8}{'bandSD':>8}")
    for d in rows:
        q = d["quarter_sd"]
        print(f"{d['label']:<10}{d['tissue']:<14}{d['gamma_pool']:>11.4f}"
              f"{d['gamma_median']:>11.4f}{d['n_eff_alpha']:>13.2f}"
              f"{d['top1_gene'][:9]:>10}{d['top1_alpha_share']*100:>7.1f}%"
              f"{(('%8.4f' % q) if q is not None else '       -')}")

    # ---------------------------------------------------------- variance components
    comp = {}
    for tis in sorted({d["tissue"] for d in rows}):
        grp = [d for d in rows if d["tissue"] == tis]
        k = len(grp)
        g = np.array([d["gamma_pool"] for d in grp], float)
        if k >= 2 and all(d["quarter_sd"] is not None for d in grp):
            n_rep = 4
            msw = float(np.mean([d["quarter_sd"] ** 2 for d in grp]))
            msb = float(n_rep * np.var(g, ddof=1))
            sb2 = max(0.0, (msb - msw) / n_rep)
            comp[tis] = dict(k=k, mean=float(g.mean()), sd_between=float(np.std(g, ddof=1)),
                             msw=msw, msb=msb, sigma2_between=sb2,
                             sigma_between=float(np.sqrt(sb2)),
                             within_sd=float(np.sqrt(msw)),
                             icc=float(sb2 / (sb2 + msw)) if (sb2 + msw) > 0 else None)
        else:
            qs = [d["quarter_sd"] for d in grp if d["quarter_sd"] is not None]
            comp[tis] = dict(k=k, mean=float(g.mean()), sd_between=None, msw=None,
                             msb=None, sigma2_between=None, sigma_between=None,
                             within_sd=float(np.sqrt(np.mean(np.square(qs)))) if qs else None,
                             icc=None)

    print("\nvariance components (bands as replicates, one-way random effects)")
    for tis, c in comp.items():
        if c["k"] < 2:
            print(f"  {tis:<14} k=1  mean {c['mean']:.4f}  within {c['within_sd']:.4f}"
                  f"  (no between-section estimate)")
            continue
        if c["sigma_between"] is None:
            print(f"  {tis:<14} k={c['k']}  sd_between {c['sd_between']:.4f}  "
                  f"(band SD missing, ICC not computed)")
            continue
        print(f"  {tis:<14} k={c['k']}  mean {c['mean']:.4f}  sd_between {c['sd_between']:.4f}  "
              f"within {c['within_sd']:.4f}  sigma_b {c['sigma_between']:.4f}  ICC {c['icc']:.2f}")

    # ---------------------------------------------- the gap, in between-section SDs
    gap, pairs = None, {}
    tis_list = sorted(comp)
    for i in range(len(tis_list)):
        for j in range(i + 1, len(tis_list)):
            t1, t2 = tis_list[i], tis_list[j]
            ca, cb = comp[t1], comp[t2]
            d = ca["mean"] - cb["mean"]
            # within-section SD as a lower bound on the between-section SD
            ws = [c["within_sd"] for c in (ca, cb) if c["within_sd"]]
            sw = float(np.sqrt(np.mean(np.square(ws)))) if ws else None
            pairs[f"{t1}-{t2}"] = dict(delta=float(d), within_sd=sw,
                                     ratio=float(abs(d) / sw) if sw else None)
    if "colon_cancer-mouse_brain" in pairs:
        gap = pairs["colon_cancer-mouse_brain"]
        if gap["within_sd"]:
            print(f"\ncolon - brain = {gap['delta']:+.4f}; within-section SD (lower bound) = "
                  f"{gap['within_sd']:.4f}; ratio {gap['ratio']:.2f}")
    for k, v in pairs.items():
        if k == "colon_cancer-mouse_brain":
            continue
        if v["within_sd"]:
            print(f"{k:<28} delta {v['delta']:+.4f}  ratio {v['ratio']:.2f}")

    # ---------------------------------------------------- cross-section regression
    x = np.array([d["top1_alpha_share"] for d in rows], float)
    y = np.array([d["gamma_pool"] for d in rows], float)
    reg = None
    if len(x) >= 3 and np.std(x) > 0:
        r = pearson(x, y)
        b = float(np.polyfit(x, y, 1)[0])
        p = perm_p(x, y)
        reg = dict(n=len(x), r=r, slope=b, perm_p=p)
        print(f"\npooled gamma vs head-transcript share of the exponent: "
              f"n={len(x)}  r={r:+.3f}  slope={b:+.3f}  permutation p={p:.4f}")

    # ------------------------------------------------------------------- MDE table
    mde = {}
    for tis, c in comp.items():
        if not c["sigma_between"]:
            continue
        mde[tis] = {str(n): float(np.exp(2.802 * c["sigma_between"] * np.sqrt(2 / n)))
                    for n in (1, 2, 3, 4, 6, 8)}
    if mde:
        print("\nMDE fold-change between two tissues, by sections per tissue")
        print("  tissue          " + "".join(f"{n:>9}" for n in ("1", "2", "3", "4", "6", "8")))
        for tis, m in mde.items():
            print(f"  {tis:<14}" + "".join(f"{m[n]:>9.2f}" for n in ("1", "2", "3", "4", "6", "8")))

    res = dict(n_sections=len(rows), rows=rows, components=comp,
               colon_brain_gap=gap, pairs=pairs, regression=reg, mde=mde)
    p = Path(a.out)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print("\nwritten", p)


if __name__ == "__main__":
    main()

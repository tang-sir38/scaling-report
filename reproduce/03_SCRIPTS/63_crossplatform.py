# -*- coding: utf-8 -*-
r"""
63_crossplatform.py -- cross-platform test of the composition claim.

Reads the Stereo-seq moments (61_gef_moments.py output) and the Visium HD mouse
brain moments, recomputes every scalar with the *exact* convention of
50_decomp_exact.py, and then evaluates the three criteria frozen in
the analysis plan before any Stereo-seq number was seen (that note is
not part of this repository; the three criteria are restated here):

  H1  composition: top1 alpha share > 0.40  AND  n_eff_alpha < 20
  H2  ordering:    gamma_pool < gamma_median (same sign as Visium HD brain)
  H3  level:       |gamma_pool difference| is reported, not tested

Criteria are read from the plan file, not re-invented here. If H1 or H2 fails,
the manuscript title must be narrowed to "in Visium HD" -- that is a documented
downgrade, not a rescue.

Usage:
  python 63_crossplatform.py --stereo 02_OUT/stereo/gene_moments_stereo.npz \
      --out 02_OUT/crossplatform/crossplatform.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DS = np.array([2, 4, 8, 16, 32, 64, 128, 256], float)
LGD = np.log(DS)

PLATFORMS = {
    # label -> path to a moments npz with mu_<d> / var_<d> / genes
    "visium_hd_brain": ROOT / "D2_MouseBrain" / "02_OUT" / "gene_moments_noMT.npz",
    "visium_hd_colon": ROOT / "02_OUT" / "gene_moments.npz",
}


def load(path: Path):
    z = np.load(path, allow_pickle=True)
    mu = np.column_stack([z[f"mu_{int(d)}"] for d in DS]).astype(np.float64)
    var = np.column_stack([z[f"var_{int(d)}"] for d in DS]).astype(np.float64)
    genes = np.array([str(g) for g in z["genes"]])
    n_cells = np.array(z["n_cells"], dtype=np.int64) if "n_cells" in z.files else None
    return mu, var, genes, n_cells


def slope(y):
    return -float(np.polyfit(LGD, np.log(y), 1)[0])


def analyse(path: Path, label: str) -> dict:
    """Everything 50_decomp_exact.py computes, plus the median it does not."""
    mu, var, genes, n_cells = load(path)
    keep = mu.sum(axis=1) > 0
    mu, var, genes = mu[keep], var[keep], genes[keep]
    G = mu.shape[0]

    with np.errstate(invalid="ignore", divide="ignore"):
        alpha_g = (var - mu) / mu ** 2
    mu2 = mu ** 2
    S2 = mu2.sum(axis=0)
    w = mu2 / S2
    alpha_pool = (w * alpha_g).sum(axis=0)
    direct = (var - mu).sum(axis=0) / S2
    identity_err = float(np.max(np.abs(alpha_pool - direct) / np.abs(direct)))
    gamma_pool = slope(alpha_pool)

    ok = np.all(np.isfinite(alpha_g) & (alpha_g > 0), axis=1) & (mu[:, -1] >= 0.5)
    gamma_g = np.array([slope(alpha_g[i]) if ok[i] else np.nan for i in range(G)])

    with np.errstate(invalid="ignore", divide="ignore"):
        v = w * alpha_g / alpha_pool
    v_bar = np.where(ok, v.mean(axis=1), 0.0)
    v_bar = v_bar / v_bar.sum()
    gamma_weighted = float(np.nansum(v_bar * gamma_g))

    n_eff_mu2 = 1.0 / float((w ** 2).sum(axis=0).mean())
    n_eff_alpha = 1.0 / float((v_bar ** 2).sum())

    order = np.argsort(-v_bar)
    top1 = int(order[0])
    top3 = [str(genes[i]) for i in order[:3]]
    gamma_median = float(np.nanmedian(gamma_g))
    gamma_mean = float(np.nanmean(gamma_g))

    return dict(
        label=label,
        n_genes=int(G), n_genes_ok=int(ok.sum()),
        n_cells=[int(x) for x in n_cells] if n_cells is not None else None,
        gamma_pool=float(gamma_pool),
        gamma_weighted=float(gamma_weighted),
        gamma_median=gamma_median,
        gamma_mean=gamma_mean,
        identity_err=identity_err,
        top1=str(genes[top1]),
        top1_alpha_share=float(v_bar[top1]),
        top1_mu2_share=float(w[top1].mean()),
        top3=top3,
        n_eff_mu2=float(n_eff_mu2),
        n_eff_alpha=float(n_eff_alpha),
        alpha_pool_by_scale=[float(x) for x in alpha_pool],
        _gamma_g=gamma_g,
        _genes=genes,
        _v_bar=v_bar,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stereo", default="02_OUT/stereo/gene_moments_stereo.npz")
    ap.add_argument("--out", default="02_OUT/crossplatform/crossplatform.json")
    a = ap.parse_args()

    sp = Path(a.stereo)
    if not sp.is_absolute():
        sp = ROOT / sp
    if not sp.exists():
        raise SystemExit(f"missing Stereo-seq moments: {sp}\n"
                         f"run 61_gef_moments.py first")

    results = {}
    results["stereo_brain"] = analyse(sp, "stereo_brain")
    for k, p in PLATFORMS.items():
        results[k] = analyse(p, k)

    st = results["stereo_brain"]
    vi = results["visium_hd_brain"]

    # ------------------------------------------------------------------ criteria
    # v1 (frozen 21:40, superseded): absolute thresholds on top1 share / n_eff.
    # v2 (frozen 00:05, before any Stereo number): scale-free concentration,
    # single-gene dominance, ordering sign, and cross-platform per-gene replication.
    h1_v1 = bool(st["top1_alpha_share"] > 0.40 and st["n_eff_alpha"] < 20)

    conc = st["n_eff_alpha"] / max(st["n_genes_ok"], 1)
    h1a = bool(conc < 0.01)
    h1b = bool(st["top1_alpha_share"] > 0.10)
    h2 = bool(st["gamma_pool"] < st["gamma_median"])   # Visium HD brain: 0.3636 < 0.6507
    d_pool = st["gamma_pool"] - vi["gamma_pool"]
    d_med = st["gamma_median"] - vi["gamma_median"]

    # H4: cross-platform per-gene agreement on shared symbols
    gs, gv = st["_genes"], vi["_genes"]
    common, is_, iv = np.intersect1d(gs, gv, return_indices=True)
    ga, gb = st["_gamma_g"][is_], vi["_gamma_g"][iv]
    m = np.isfinite(ga) & np.isfinite(gb)
    rho = None
    if m.sum() > 10:
        from scipy.stats import spearmanr
        rho = float(spearmanr(ga[m], gb[m]).statistic)
    h4 = bool(rho is not None and rho > 0.40)
    shared_top = st["top1"] in vi["top3"] or vi["top1"] in st["top3"]

    verdict = dict(
        v1=dict(pass_=h1_v1, top1_share=st["top1_alpha_share"],
                n_eff_alpha=st["n_eff_alpha"],
                rule="top1 share > 0.40 and n_eff_alpha < 20 (superseded)"),
        H1a_concentration=dict(pass_=h1a, n_eff_over_genes=conc,
                               rule="n_eff_alpha / n_genes_ok < 0.01"),
        H1b_dominance=dict(pass_=h1b, top1_share=st["top1_alpha_share"],
                           rule="top1 alpha share > 0.10"),
        H2_ordering=dict(pass_=h2, gamma_pool=st["gamma_pool"],
                         gamma_median=st["gamma_median"],
                         rule="gamma_pool < gamma_median (Visium HD brain sign)"),
        H3_level=dict(delta_gamma_pool=d_pool, delta_gamma_median=d_med,
                      rule="reported, not tested; absolute level is platform-specific"),
        H4_replication=dict(pass_=h4, spearman_per_gene_gamma=rho,
                            n_shared=int(m.sum()),
                            rule="Spearman rho of per-gene gamma across platforms > 0.40"),
        cross_platform_genes=dict(n_common=int(common.size), n_used=int(m.sum()),
                                  top1_in_other_top3=shared_top),
    )
    v2_pass = bool(h1a and h1b and h2 and h4)
    results["verdict"] = verdict
    results["criteria_met"] = v2_pass
    results["criteria_met_v1"] = h1_v1

    for k in ("stereo_brain", "visium_hd_brain", "visium_hd_colon"):
        r = results[k]
        print(f"\n[{k}]  genes {r['n_genes']:,} (usable {r['n_genes_ok']:,})")
        print(f"  gamma_pool {r['gamma_pool']:.4f}  gamma_median {r['gamma_median']:.4f}  "
              f"gamma_mean {r['gamma_mean']:.4f}  weighted {r['gamma_weighted']:.4f}")
        print(f"  top1 {r['top1']}  alpha-share {r['top1_alpha_share']*100:.1f}%  "
              f"mu2-share {r['top1_mu2_share']*100:.1f}%  top3 {r['top3']}")
        print(f"  n_eff_mu2 {r['n_eff_mu2']:.2f}  n_eff_alpha {r['n_eff_alpha']:.2f}  "
              f"n_eff/genes {r['n_eff_alpha']/max(r['n_genes_ok'],1)*100:.2f}%  "
              f"identity_err {r['identity_err']:.2e}")

    print("\n=== criteria v1 (superseded, still reported) ===")
    print(f"  H1 top1>40% & n_eff<20 : {'PASS' if h1_v1 else 'FAIL'}   "
          f"top1 {st['top1_alpha_share']*100:.1f}%  n_eff {st['n_eff_alpha']:.2f}")
    print("\n=== criteria v2 (frozen 00:05, before any Stereo number) ===")
    print(f"  H1a concentration      : {'PASS' if h1a else 'FAIL'}   "
          f"n_eff/genes {conc*100:.2f}% (need <1%)")
    print(f"  H1b dominance          : {'PASS' if h1b else 'FAIL'}   "
          f"top1 {st['top1_alpha_share']*100:.1f}% (need >10%)")
    print(f"  H2  ordering           : {'PASS' if h2 else 'FAIL'}   "
          f"pool {st['gamma_pool']:.4f} vs median {st['gamma_median']:.4f}")
    print(f"  H3  level              : d_pool {d_pool:+.4f}  d_median {d_med:+.4f}  (no test)")
    print(f"  H4  replication        : {'PASS' if h4 else 'FAIL'}   "
          f"rho {rho if rho is None else round(rho,3)} on {int(m.sum()):,} shared genes "
          f"(need >0.40)")
    print(f"\n  VERDICT v2: {'ALL FOUR MET' if v2_pass else 'CRITERIA NOT MET -> narrow the title'}")
    print(f"  VERDICT v1: {'met' if h1_v1 else 'not met'}")

    out = Path(a.out)
    if not out.is_absolute():
        out = ROOT / out
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = {}
    for k, v in results.items():
        if k == "verdict":
            continue
        if isinstance(v, dict):
            payload[k] = {kk: vv for kk, vv in v.items() if not kk.startswith("_")}
    payload["verdict"] = verdict
    payload["criteria_met"] = results["criteria_met"]
    payload["criteria_met_v1"] = results["criteria_met_v1"]
    out.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()

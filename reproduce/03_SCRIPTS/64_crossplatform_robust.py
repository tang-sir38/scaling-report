# -*- coding: utf-8 -*-
r"""
64_crossplatform_robust.py -- try to break the cross-platform verdict.

Two attacks, both written before the numbers were seen:

  A. contamination: the Stereo-seq top-1 gene Gm42418 is a known rRNA-derived
     locus that dominates mouse capture data and is routinely discarded. If the
     H1 verdict depends on it, the verdict is about a technical artefact.

  B. shared abundance: per-gene gamma correlates with abundance on every
     platform, so a cross-platform rho can be manufactured entirely by both
     platforms ranking genes by abundance the same way. Test the rank
     correlation of gamma *after* removing log abundance, by rank residual.

Neither attack may change the reported verdict silently: both outcomes printed.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def load_module():
    spec = importlib.util.spec_from_file_location("x63", ROOT / "03_SCRIPTS" / "63_crossplatform.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def rank_resid(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Residual of rank(y) on rank(x): removes any monotone abundance effect."""
    from scipy.stats import rankdata
    ry, rx = rankdata(y), rankdata(x)
    A = np.column_stack([np.ones_like(rx), rx])
    beta, *_ = np.linalg.lstsq(A, ry, rcond=None)
    return ry - A @ beta


def main():
    m = load_module()
    DS = m.DS

    def analyse_filtered(path, drop=()):
        mu, var, genes, _ = m.load(path)
        keep = mu.sum(axis=1) > 0
        if drop:
            keep &= ~np.isin(genes, np.array(list(drop)))
        mu, var, genes = mu[keep], var[keep], genes[keep]
        with np.errstate(invalid="ignore", divide="ignore"):
            alpha_g = (var - mu) / mu ** 2
        w = mu ** 2 / (mu ** 2).sum(axis=0)
        alpha_pool = (w * alpha_g).sum(axis=0)
        ok = np.all(np.isfinite(alpha_g) & (alpha_g > 0), axis=1) & (mu[:, -1] >= 0.5)
        gamma_g = np.array([m.slope(alpha_g[i]) if ok[i] else np.nan for i in range(len(genes))])
        with np.errstate(invalid="ignore", divide="ignore"):
            v = w * alpha_g / alpha_pool
        v_bar = np.where(ok, v.mean(axis=1), 0.0)
        v_bar = v_bar / v_bar.sum()
        order = np.argsort(-v_bar)
        return dict(
            gamma_pool=m.slope(alpha_pool),
            gamma_median=float(np.nanmedian(gamma_g)),
            top1=str(genes[order[0]]), top1_share=float(v_bar[order[0]]),
            n_eff_alpha=1.0 / float((v_bar ** 2).sum()),
            n_genes_ok=int(ok.sum()),
            gamma_g=gamma_g, genes=genes, v_bar=v_bar, ok=ok,
            abundance=mu[:, int(np.where(DS == 32)[0][0])],
        )

    st_path = ROOT / "02_OUT" / "stereo" / "gene_moments_stereo.npz"
    vi_path = ROOT / "D2_MouseBrain" / "02_OUT" / "gene_moments_noMT.npz"

    RRNA = ("Gm42418", "Rn45s", "Rn18s", "Rn28s", "Rn4.5s", "Gm26917", "Gm23935")

    out = {}
    print("=== A. contamination attack: drop rRNA-derived loci ===")
    for tag, drop in (("with", ()), ("without", RRNA)):
        r = analyse_filtered(st_path, drop)
        h1a = r["n_eff_alpha"] / max(r["n_genes_ok"], 1) < 0.01
        h1b = r["top1_share"] > 0.10
        print(f"  {tag:8s} genes {r['n_genes_ok']:>6,}  top1 {r['top1']:<10s} "
              f"share {r['top1_share']*100:5.1f}%  n_eff {r['n_eff_alpha']:6.2f}  "
              f"H1a {'PASS' if h1a else 'FAIL'}  H1b {'PASS' if h1b else 'FAIL'}")
        out[f"stereo_{tag}_rrna"] = dict(genes=r["n_genes_ok"], top1=r["top1"],
                                         top1_share=r["top1_share"],
                                         n_eff_alpha=r["n_eff_alpha"],
                                         h1a=bool(h1a), h1b=bool(h1b))
    dropped = [g for g in RRNA if g in set(analyse_filtered(st_path)["genes"].tolist())]
    print(f"  actually dropped: {dropped}")

    print("\n=== B. shared-abundance attack: rho after removing log abundance ===")
    a = analyse_filtered(st_path, RRNA)
    b = analyse_filtered(vi_path)
    common, ia, ib = np.intersect1d(a["genes"], b["genes"], return_indices=True)
    ga, gb = a["gamma_g"][ia], b["gamma_g"][ib]
    ab_a, ab_b = a["abundance"][ia], b["abundance"][ib]
    keep = np.isfinite(ga) & np.isfinite(gb) & (ab_a > 0) & (ab_b > 0)
    ga, gb, ab_a, ab_b = ga[keep], gb[keep], ab_a[keep], ab_b[keep]
    from scipy.stats import spearmanr, pearsonr
    raw = spearmanr(ga, gb).statistic
    ra = rank_resid(ga, np.log(ab_a))
    rb = rank_resid(gb, np.log(ab_b))
    part = pearsonr(ra, rb)[0]
    # abundance itself, for reference
    ab_r = spearmanr(np.log(ab_a), np.log(ab_b)).statistic
    print(f"  n shared genes        {keep.sum():,}")
    print(f"  rho(gamma) raw        {raw:+.3f}")
    print(f"  rho(gamma) | abundance {part:+.3f}   (rank residual on log mu)")
    print(f"  rho(abundance)        {ab_r:+.3f}   <- how much of raw rho could be abundance")
    print(f"  variance explained by abundance alone: {ab_r**2*100:.0f}%")
    out["replication"] = dict(n=int(keep.sum()), rho_raw=float(raw),
                              rho_partial_abundance=float(part),
                              rho_abundance=float(ab_r),
                              h4_raw=bool(raw > 0.40),
                              h4_partial=bool(part > 0.40))

    # concentration of the exponent over the primary (rRNA-included) Stereo
    # analysis, and the same statistic on Visium HD, computed the same way
    st_all = analyse_filtered(st_path)
    def top_share(df, k=100):
        v = np.sort(df["v_bar"])[::-1]
        return float(v[:k].sum() / v.sum())
    out["replication"]["top100_stereo"] = top_share(st_all)
    out["replication"]["top100_visium"] = top_share(b)
    out["replication"]["top100_stereo_without_rrna"] = top_share(a)
    print(f"  top-100 share of alpha: Stereo {out['replication']['top100_stereo']*100:.1f}%  "
          f"Visium {out['replication']['top100_visium']*100:.1f}%  "
          f"(Stereo without rRNA {out['replication']['top100_stereo_without_rrna']*100:.1f}%)")

    print("\n=== verdict under both attacks ===")
    surv = (out["stereo_without_rrna"]["h1a"] and out["stereo_without_rrna"]["h1b"]
            and out["replication"]["h4_partial"])
    print(f"  H1a/H1b after dropping rRNA loci : "
          f"{'PASS' if out['stereo_without_rrna']['h1a'] and out['stereo_without_rrna']['h1b'] else 'FAIL'}")
    print(f"  H4 after removing abundance      : "
          f"{'PASS' if out['replication']['h4_partial'] else 'FAIL'}")
    print(f"  SURVIVES BOTH: {'YES' if surv else 'NO -- report the weaker conclusion'}")

    p = ROOT / "02_OUT" / "crossplatform" / "robustness.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(out, indent=1), encoding="utf-8")

    # per-gene tables, so that Fig. 7c is drawn from the same arrays that
    # produced the correlation rather than from a second pass over the data
    np.savez_compressed(
        p.parent / "pergene_shared.npz",
        genes=common[keep].astype(str),
        gamma_stereo=ga, gamma_visium=gb,
        logmu_stereo=np.log(ab_a), logmu_visium=np.log(ab_b),
        resid_stereo=ra, resid_visium=rb,
    )
    np.savez_compressed(
        p.parent / "vbar_stereo.npz",
        genes=a["genes"].astype(str), v_bar=a["v_bar"], gamma_g=a["gamma_g"],
        ok=a["ok"],
    )
    np.savez_compressed(
        p.parent / "vbar_stereo_primary.npz",
        genes=st_all["genes"].astype(str), v_bar=st_all["v_bar"],
        gamma_g=st_all["gamma_g"], ok=st_all["ok"],
    )
    np.savez_compressed(
        p.parent / "vbar_visium.npz",
        genes=b["genes"].astype(str), v_bar=b["v_bar"], gamma_g=b["gamma_g"],
        ok=b["ok"],
    )
    print(f"  wrote {p}")


if __name__ == "__main__":
    main()

r"""Turn a fit into (a) a machine-readable summary and (b) the reporting rule.

The point of the tool is not the plot. It is that the same four sentences come
out of anybody's data, so two papers can be compared.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from .core import Moments, Fit, fit, identifiability, leave_one_out, warnings_for

Z_ALPHA = 1.959964
Z_POWER = 0.841621


def summarise(m: Moments, f: Fit, *, name: str | None = None,
              loo_top: int = 3, loo_extra: tuple[str, ...] = ()) -> dict:
    """Everything a reader needs, and nothing that needs the raw matrix."""
    ident = identifiability(f)
    top = f.top_table(loo_top)
    loo = []
    for rank, (i, share, gamma) in enumerate(top):
        dropped = leave_one_out(m, f, i)
        loo.append(dict(
            rank=rank + 1, gene=str(m.genes[i]),
            share_of_exponent=share, gamma_gene=gamma,
            gamma_pool_without=float(dropped),
            delta=float(dropped - f.gamma_pool),
        ))
    for g in loo_extra:
        idx = np.where(m.genes == g)[0]
        if idx.size:
            i = int(idx[0])
            dropped = leave_one_out(m, f, i)
            loo.append(dict(rank=None, gene=g,
                            share_of_exponent=float(f.v_bar[i]),
                            gamma_gene=float(f.gamma_g[i]) if f.ok[i] else None,
                            gamma_pool_without=float(dropped),
                            delta=float(dropped - f.gamma_pool)))

    rep = dict(
        schema="scaling-report/1",
        section=name,
        ladder_um=[float(d) for d in m.ds],
        n_genes_in_table=int(m.n_genes),
        n_genes_fitted=int(f.n_fitted),
        gene_filter=dict(mu_at_largest_bin_at_least=f.mu_min),
        gamma_pool=float(f.gamma_pool),
        alpha_pool=[float(a) for a in f.alpha_pool],
        gamma_pool_sensitivity=dict(
            max_abs_delta_over_top=float(max((abs(x["delta"]) for x in loo),
                                             default=0.0)),
            leave_one_out=loo,
        ),
        per_gene=dict(
            median=f.median, q25=f.q25, q75=f.q75,
            iqr=float(f.q75 - f.q25),
            mean=float(np.nanmean(f.gamma_g[f.ok])) if f.n_fitted else None,
            min=float(np.nanmin(f.gamma_g[f.ok])) if f.n_fitted else None,
            max=float(np.nanmax(f.gamma_g[f.ok])) if f.n_fitted else None,
            n=int(f.n_fitted),
        ),
        concentration=dict(
            n_eff_variance=f.n_eff_w,
            n_eff_exponent=f.n_eff_v,
            ratio=f.n_eff_w / f.n_eff_v if f.n_eff_v > 0 else None,
            top1_gene=str(m.genes[f.top1]),
            top1_share=f.top1_share,
            top10_share=f.share_top(10),
            top100_share=f.share_top(100),
        ),
        identifiability=ident,
        mde_template=dict(
            z_alpha=Z_ALPHA, z_power=Z_POWER,
            formula="MDE = (1.96 + 0.84) * sqrt(2) * SD_between / sqrt(n per group)",
            note="SD_between cannot be estimated from one section. Use at least two.",
        ),
        warnings=warnings_for(m, f),
    )
    rep["ordering"] = _ordering(rep)
    return rep


def _ordering(rep: dict) -> dict:
    """Does the pooled value sit above or below the typical transcript?"""
    g, med = rep["gamma_pool"], rep["per_gene"]["median"]
    return dict(
        pooled_minus_median=float(g - med),
        pooled_below_median=bool(g < med),
        reversed=bool(g < med),
        interpretation=(
            "the pooled value is lower than the typical transcript, because the "
            "weight sits on transcripts that decorrelate quickly"
            if g < med else
            "the pooled value is close to or above the typical transcript"),
    )


def rule_text(rep: dict) -> list[str]:
    """The four sentences the paper argues for, filled with this run's numbers."""
    c = rep["concentration"]
    pg = rep["per_gene"]
    ident = rep["identifiability"]
    ladder = "–".join([f"{rep['ladder_um'][0]:g}", f"{rep['ladder_um'][-1]:g}"])
    out = [
        f"Report the pooled exponent together with the ladder ({ladder} µm, "
        f"{len(rep['ladder_um'])} steps) and the gene filter "
        f"(mean ≥ {rep['gene_filter']['mu_at_largest_bin_at_least']:g} counts at the "
        f"largest bin), because the absolute value is a property of both.",

        f"Next to it report the concentration: {c['n_eff_exponent']:.1f} effective "
        f"genes under the exponent weight, and {c['top1_gene']} carrying "
        f"{c['top1_share']*100:.0f}% of it.",

        f"Report the per-gene distribution as a separate sentence "
        f"(median {pg['median']:.3f}, IQR {pg['q25']:.3f}–{pg['q75']:.3f}, "
        f"n = {pg['n']:,}), and never as a substitute for the pooled value: with the "
        f"pooled value held fixed the median can sit anywhere in "
        f"{ident['median_at_fixed_pool_min']:.2f}–{ident['median_at_fixed_pool_max']:.2f}."
        if ident["defined"] else
        f"Report the per-gene distribution as a separate sentence "
        f"(median {pg['median']:.3f}, IQR {pg['q25']:.3f}–{pg['q75']:.3f}, "
        f"n = {pg['n']:,}).",

        "Compare sections on the pooled exponent only, and only with at least two "
        "sections per group; a per-gene median does not reproduce between sections "
        "of one tissue, so it cannot compare two tissues either.",
    ]
    return out


def to_markdown(rep: dict) -> str:
    c, pg, ident = rep["concentration"], rep["per_gene"], rep["identifiability"]
    L = []
    title = rep["section"] or "section"
    L.append(f"# scaling_report · {title}")
    L.append("")
    L.append(f"Ladder: {rep['ladder_um'][0]:g}–{rep['ladder_um'][-1]:g} µm, "
             f"{len(rep['ladder_um'])} steps.  Genes in table: "
             f"{rep['n_genes_in_table']:,}; fitted: {pg['n']:,} "
             f"(filter: mean ≥ {rep['gene_filter']['mu_at_largest_bin_at_least']:g} "
             f"counts at the largest bin).")
    L.append("")
    L.append("| Quantity | Value |")
    L.append("|---|---|")
    L.append(f"| Pooled exponent γ_pool | **{rep['gamma_pool']:.4f}** |")
    L.append(f"| Per-gene median (IQR) | {pg['median']:.3f} "
             f"({pg['q25']:.3f}–{pg['q75']:.3f}) |")
    L.append(f"| Pooled − median | {rep['ordering']['pooled_minus_median']:+.3f} |")
    L.append(f"| Effective genes, exponent weight | {c['n_eff_exponent']:.1f} |")
    L.append(f"| Effective genes, variance weight | {c['n_eff_variance']:.1f} |")
    L.append(f"| Top transcript | {c['top1_gene']}, "
             f"{c['top1_share']*100:.1f}% of the exponent |")
    L.append(f"| Top 10 / 100 share | {c['top10_share']*100:.1f}% / "
             f"{c['top100_share']*100:.1f}% |")
    L.append(f"| Largest single-gene move | "
             f"{rep['gamma_pool_sensitivity']['max_abs_delta_over_top']:.3f} |")
    if ident["defined"]:
        L.append(f"| Median admitted by the same γ_pool | "
                 f"{ident['median_at_fixed_pool_min']:.2f}–"
                 f"{ident['median_at_fixed_pool_max']:.2f} |")
    L.append("")
    L.append("## What this section does not identify")
    L.append("")
    if ident["defined"]:
        L.append(f"Holding γ_pool = {rep['gamma_pool']:.4f} fixed, the per-gene median "
                 f"can be placed anywhere between "
                 f"{ident['median_at_fixed_pool_min']:.2f} and "
                 f"{ident['median_at_fixed_pool_max']:.2f} by moving the heaviest "
                 f"transcript and setting the rest to one constant. "
                 f"d(median)/d(anchor) = {ident['d_median_d_anchor']:+.2f}. "
                 f"The pooled value therefore does not rank tissues by their "
                 f"transcript-level scaling; only the composition part is identified.")
    else:
        L.append("Not computable: the top transcript carries the whole exponent.")
    L.append("")
    L.append("## Reporting rule")
    L.append("")
    for i, s in enumerate(rule_text(rep), 1):
        L.append(f"{i}. {s}")
    L.append("")
    if rep["warnings"]:
        L.append("## Warnings")
        L.append("")
        for w in rep["warnings"]:
            L.append(f"- {w}")
        L.append("")
    L.append("_Generated by scaling_report. Numbers are reproducible from the "
             "moments table alone; no raw matrix is needed after this point._")
    return "\n".join(L) + "\n"


def compare(a: dict, b: dict) -> dict:
    """Difference between two sections, with the uncertainty deliberately withheld."""
    d = b["gamma_pool"] - a["gamma_pool"]
    return dict(
        a=a["section"], b=b["section"],
        gamma_pool_a=a["gamma_pool"], gamma_pool_b=b["gamma_pool"],
        delta=float(d), relative=float(d / a["gamma_pool"]),
        median_shift=float(b["per_gene"]["median"] - a["per_gene"]["median"]),
        note=("This is a difference of two numbers. Its standard error needs the raw "
              "binned matrix, because bins are spatially correlated: resampling bins "
              "independently understates the error by roughly the ratio of the "
              "analysis scale to the correlation length. Resample spatial blocks "
              "(e.g. 512 um superblocks) of each section, or use at least two "
              "sections per group and the MDE in mde_template."),
    )


def write(rep: dict, outdir: str | Path, *, stem: str = "scaling_report") -> dict:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    jp = out / f"{stem}.json"
    mp = out / f"{stem}.md"
    jp.write_text(json.dumps(rep, indent=2, ensure_ascii=False), encoding="utf-8")
    mp.write_text(to_markdown(rep), encoding="utf-8")
    return {"json": str(jp), "markdown": str(mp)}


def mde(sd_between: float, n_per_group: int) -> float:
    """Minimum detectable difference at 80% power, two-sided 5%."""
    if n_per_group < 1:
        raise ValueError("n_per_group must be >= 1")
    return (Z_ALPHA + Z_POWER) * math.sqrt(2) * sd_between / math.sqrt(n_per_group)

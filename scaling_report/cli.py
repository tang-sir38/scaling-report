"""Command line entry point: python -m scaling_report FILE [options]."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .core import Moments, MomentsError, fit
from .report import compare, mde, summarise, to_markdown, write


def _build(m: Moments, args) -> dict:
    f = fit(m, mu_min=args.min_mu)
    rep = summarise(m, f, name=args.name)
    return rep


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="scaling-report",
        description="Report the scaling exponent of binned spatial transcriptomics "
                    "in a way that can be compared between sections.",
        epilog="Exit codes: 0 ok, 2 unreadable input, 3 the pooled exponent is not "
               "positive at some scale.")
    ap.add_argument("moments", type=Path,
                    help=".npz (ds, genes, mu_<d>, var_<d>) or a tidy "
                         ".tsv/.csv with gene, scale, mu, var")
    ap.add_argument("--out", type=Path, default=Path("scaling_report_out"))
    ap.add_argument("--name", default=None, help="label for this section")
    ap.add_argument("--min-mu", type=float, default=0.5,
                    help="gene filter: minimum mean counts at the largest bin "
                         "(default 0.5)")
    ap.add_argument("--compare", type=Path, default=None,
                    help="second moments table; reports the difference and refuses "
                         "to invent a confidence interval")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--json-only", action="store_true")
    a = ap.parse_args(argv)

    try:
        m = Moments.load(a.moments)
    except (OSError, MomentsError) as e:
        print(f"scaling-report: cannot read {a.moments}: {e}", file=sys.stderr)
        return 2

    if a.name is None:
        a.name = a.moments.stem

    try:
        rep = _build(m, a)
    except MomentsError as e:
        print(f"scaling-report: {e}", file=sys.stderr)
        return 3

    paths = write(rep, a.out)

    if a.compare:
        try:
            m2 = Moments.load(a.compare)
        except (OSError, MomentsError) as e:
            print(f"scaling-report: cannot read {a.compare}: {e}", file=sys.stderr)
            return 2
        rep2 = _build(m2, argparse.Namespace(min_mu=a.min_mu, name=a.compare.stem))
        cmp = compare(rep, rep2)
        (a.out / "compare.json").write_text(
            json.dumps(cmp, indent=2), encoding="utf-8")
        with (a.out / "compare.md").open("w", encoding="utf-8") as fh:
            fh.write(f"# {cmp['a']} vs {cmp['b']}\n\n"
                     f"γ_pool {cmp['gamma_pool_a']:.4f} vs "
                     f"{cmp['gamma_pool_b']:.4f}, Δ = {cmp['delta']:+.4f} "
                     f"({cmp['relative']*100:+.1f}%).\n"
                     f"Per-gene median shift: {cmp['median_shift']:+.4f}.\n\n"
                     f"> {cmp['note']}\n")
        if not a.quiet:
            print(f"compare: Δγ_pool = {cmp['delta']:+.4f}, "
                  f"median shift {cmp['median_shift']:+.4f}")

    if not a.quiet:
        c = rep["concentration"]
        print(f"{rep['section']}: γ_pool = {rep['gamma_pool']:.4f} "
              f"| per-gene median {rep['per_gene']['median']:.3f} "
              f"(n={rep['per_gene']['n']:,}) "
              f"| {c['top1_gene']} {c['top1_share']*100:.0f}% "
              f"| n_eff {c['n_eff_exponent']:.1f}")
        for w in rep["warnings"]:
            print(f"  warning: {w}")
        print(f"wrote {paths['json']}")
        if not a.json_only:
            print(f"wrote {paths['markdown']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

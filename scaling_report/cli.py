"""Command line entry point.

Two ways in, one way out.

    scaling-report moments.npz                     # a table you already have
    scaling-report --from-10x run/binned_outputs   # a Space Ranger run
    scaling-report --from-mtx DIR --base-um 2      # one Matrix Market directory
    scaling-report --demo DIR                      # write a synthetic section
    scaling-report --selftest                      # check this machine

Every path writes the same four files: JSON, Markdown, HTML, and the advice.
"""
from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path

from . import __version__
from .core import Moments, MomentsError, fit
from .report import compare, rule_text, summarise, write


def _build(m: Moments, args) -> tuple[dict, dict]:
    from .advice import verdict

    f = fit(m, mu_min=args.min_mu)
    rep = summarise(m, f, name=args.name)
    rep["reporting_rule"] = rule_text(rep)
    adv = verdict(m, f, fpr_target=args.fpr_target)
    return rep, adv


def _emit(rep: dict, adv: dict, out: Path, *, emit_html: bool, quiet: bool,
          meta: dict | None = None) -> dict:
    paths = write(rep, out)
    (out / "advice.json").write_text(
        json.dumps(adv, indent=2, ensure_ascii=False), encoding="utf-8")
    paths["advice"] = str(out / "advice.json")
    if emit_html:
        from .html import to_html, write_html
        doc = to_html(rep, adv, meta=meta or {}, version=__version__)
        p = write_html(doc, out / "scaling_report.html")
        paths["html"] = str(p)
    if not quiet:
        print(f"{rep['section']}: gamma_pool = {rep['gamma_pool']:.4f} | "
              f"per-gene median {rep['per_gene']['median']:.3f} "
              f"(n={rep['per_gene']['n']:,}) | gamma "
              f"{rep['concentration']['n_eff_exponent']:.1f} effective genes")
        for s in adv.get("sentences", []):
            print(f"  {s}")
        for k in ("json", "markdown", "html", "advice"):
            if k in paths:
                print(f"wrote {paths[k]}")
    return paths


def _selftest() -> int:
    """End-to-end check that needs nothing but this package."""
    from .advice import fpr_under_poisson
    from .demo import make_demo
    from .tenx import ladder_from_run

    ok = 0
    fail = 0

    def check(name: str, cond: bool, detail: str = "") -> None:
        nonlocal ok, fail
        if cond:
            ok += 1
            print(f"  PASS  {name}")
        else:
            fail += 1
            print(f"  FAIL  {name}  {detail}")

    print("scaling-report selftest")
    print(f"  version {__version__}")

    # 1. the arithmetic of the false-positive rate, against the paper's value
    check("phi = 1 gives 5%", abs(fpr_under_poisson(1.0) - 0.05) < 1e-9,
          f"{fpr_under_poisson(1.0):.4f}")
    check("phi = 30 gives ~72%", 0.70 < fpr_under_poisson(30.0) < 0.74,
          f"{fpr_under_poisson(30.0):.4f}")

    with tempfile.TemporaryDirectory(prefix="sr_selftest_") as td:
        # 2. a synthetic section with gamma known by construction
        make_demo(Path(td), n=256, n_genes=120, gamma=0.55, seed=7, verbose=False)
        m, meta = ladder_from_run(Path(td), ds=(2, 4, 8, 16, 32, 64))
        check("six scales recovered", m.ds.size == 6, str(m.ds))
        check("gene count preserved", m.n_genes == 120, str(m.n_genes))
        check("no negative moments",
              bool(((m.mu >= 0) & (m.var >= 0)).all()))
        f = fit(m, mu_min=0.05)
        check("gamma recovers the constructed value",
              abs(f.gamma_pool - 0.55) < 0.25, f"{f.gamma_pool:.3f} vs 0.55")
        check("gamma is positive", f.gamma_pool > 0, f"{f.gamma_pool:.3f}")

    print(f"\n{ok} passed, {fail} failed")
    return 0 if fail == 0 else 1


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="scaling-report",
        description="Report the scaling exponent of binned spatial transcriptomics, "
                    "and what a Poisson test would cost at each bin size.",
        epilog="Exit codes: 0 ok, 2 unreadable input, 3 the pooled exponent is not "
               "positive at some scale.")
    ap.add_argument("moments", nargs="?", type=Path,
                    help=".npz (ds, genes, mu_<d>, var_<d>) or a tidy "
                         ".tsv/.csv with gene, scale, mu, var")
    ap.add_argument("--from-10x", type=Path, metavar="DIR",
                    help="a Space Ranger run directory, or its binned_outputs/")
    ap.add_argument("--from-mtx", type=Path, metavar="DIR",
                    help="one Matrix Market directory (matrix.mtx[.gz], "
                         "barcodes.tsv[.gz], features.tsv[.gz])")
    ap.add_argument("--base-um", type=int, default=2,
                    help="bin edge length of the --from-mtx matrix (default 2)")
    ap.add_argument("--demo", type=Path, metavar="DIR",
                    help="write a synthetic section with a known exponent, then stop")
    ap.add_argument("--selftest", action="store_true",
                    help="check this machine end to end, then stop")
    ap.add_argument("--ds", default="2,4,8,16,32,64,128,256",
                    help="bin ladder, comma separated (default 2,...,256)")
    ap.add_argument("--out", type=Path, default=Path("scaling_report_out"))
    ap.add_argument("--name", default=None, help="label for this section")
    ap.add_argument("--min-mu", type=float, default=0.5,
                    help="gene filter: minimum mean counts at the largest bin "
                         "(default 0.5)")
    ap.add_argument("--fpr-target", type=float, default=0.10,
                    help="acceptable false-positive rate for the advice "
                         "(default 0.10)")
    ap.add_argument("--no-html", action="store_true")
    ap.add_argument("--compare", type=Path, default=None,
                    help="second moments table; reports the difference and refuses "
                         "to invent a confidence interval")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)

    if a.selftest:
        return _selftest()

    if a.demo:
        from .demo import make_demo
        make_demo(a.demo, verbose=not a.quiet)
        return 0

    sources = [a.moments is not None, a.from_10x is not None,
               a.from_mtx is not None]
    if sum(sources) == 0:
        ap.error("give a moments table, --from-10x, --from-mtx, --demo or --selftest")
    if sum(sources) > 1:
        ap.error("choose one source")

    ds = tuple(int(x) for x in str(a.ds).split(",") if x.strip())
    meta: dict = {}
    try:
        if a.from_10x:
            from .tenx import ladder_from_run
            m, meta = ladder_from_run(a.from_10x, ds=ds)
            a.name = a.name or Path(a.from_10x).resolve().name
        elif a.from_mtx:
            from .tenx import _resolve_matrix, ladder_from_one_matrix
            matrix = _resolve_matrix(a.from_mtx)
            m = ladder_from_one_matrix(matrix, a.base_um, ds=ds)
            a.name = a.name or Path(a.from_mtx).resolve().name
        else:
            m = Moments.load(a.moments)
            a.name = a.name or a.moments.stem
    except (OSError, MomentsError) as e:
        print(f"scaling-report: cannot read the input: {e}", file=sys.stderr)
        return 2

    try:
        rep, adv = _build(m, a)
    except MomentsError as e:
        print(f"scaling-report: {e}", file=sys.stderr)
        return 3

    _emit(rep, adv, a.out, emit_html=not a.no_html, quiet=a.quiet,
          meta={k: v for k, v in meta.items() if k != "source"})

    if a.compare:
        try:
            m2 = Moments.load(a.compare)
        except (OSError, MomentsError) as e:
            print(f"scaling-report: cannot read {a.compare}: {e}", file=sys.stderr)
            return 2
        rep2, _ = _build(m2, argparse.Namespace(
            min_mu=a.min_mu, name=a.compare.stem, fpr_target=a.fpr_target))
        cmp = compare(rep, rep2)
        (a.out / "compare.json").write_text(
            json.dumps(cmp, indent=2), encoding="utf-8")
        with (a.out / "compare.md").open("w", encoding="utf-8") as fh:
            fh.write(f"# {cmp['a']} vs {cmp['b']}\n\n"
                     f"gamma_pool {cmp['gamma_pool_a']:.4f} vs "
                     f"{cmp['gamma_pool_b']:.4f}, delta = {cmp['delta']:+.4f} "
                     f"({cmp['relative']*100:+.1f}%).\n"
                     f"Per-gene median shift: {cmp['median_shift']:+.4f}.\n\n"
                     f"> {cmp['note']}\n")
        if not a.quiet:
            print(f"compare: delta gamma_pool = {cmp['delta']:+.4f}, "
                  f"median shift {cmp['median_shift']:+.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

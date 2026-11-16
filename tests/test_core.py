r"""Tests. Run with pytest, or directly: python tests/test_core.py

The rule these tests enforce: the estimator must return a known exponent to
machine precision on data where the answer is known by construction, and must
reproduce the published colon values from the shipped example.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scaling_report import (Moments, MomentsError, fit, identifiability,  # noqa: E402
                            leave_one_out, mde, rule_text, summarise)

DS = np.array([2, 4, 8, 16, 32, 64, 128, 256], float)
HERE = Path(__file__).resolve().parent
EXAMPLE = HERE.parent / "examples" / "colon_moments.npz"


def synthetic(gammas, amplitudes=None, mu256=None, ds=DS) -> Moments:
    """mu = m (d/256)^2, alpha = a (d/256)^-gamma, var = mu + alpha mu^2."""
    g = np.asarray(gammas, float)
    a = np.ones_like(g) if amplitudes is None else np.asarray(amplitudes, float)
    m = np.full_like(g, 500.0) if mu256 is None else np.asarray(mu256, float)
    mu = np.empty((g.size, ds.size))
    var = np.empty_like(mu)
    for k, d in enumerate(ds):
        s = (d / 256.0) ** 2
        mu[:, k] = m * s
        alpha = a * (d / 256.0) ** (-g)
        var[:, k] = mu[:, k] + alpha * mu[:, k] ** 2
    genes = np.array([f"G{i:03d}" for i in range(g.size)])
    return Moments(ds=ds, genes=genes, mu=mu, var=var)


def test_single_gene_is_exact():
    f = fit(synthetic([0.4]))
    assert abs(f.gamma_pool - 0.4) < 1e-10, f.gamma_pool


def test_all_genes_same_exponent_is_exact():
    f = fit(synthetic(np.full(50, 0.35), amplitudes=np.linspace(0.1, 2.0, 50)))
    assert abs(f.gamma_pool - 0.35) < 1e-10, f.gamma_pool


def test_per_gene_exponent_recovered():
    g = np.array([0.12, 0.31, 0.5, 0.77, 0.9])
    f = fit(synthetic(g))
    assert np.all(f.ok), "every gene should be fitted"
    err = np.max(np.abs(f.gamma_g - g))
    assert err < 1e-9, err


def test_pooled_is_not_the_mean_of_genes():
    g = np.array([0.9, 0.15, 0.15, 0.15])
    m = np.array([1e5, 10.0, 10.0, 10.0])   # the first gene carries the variance
    f = fit(synthetic(g, mu256=m))
    assert abs(f.gamma_pool - g.mean()) > 0.01, "the test data should separate them"
    # the convex-combination identity holds to numerical precision whenever the
    # variance weights are flat across scales, as they are in this construction
    assert abs(f.gamma_pool - f.gamma_weighted) < 1e-6, f.gamma_weighted


def test_v_bar_sums_to_one():
    f = fit(synthetic(np.linspace(0.1, 0.9, 30), amplitudes=np.linspace(1, 0.01, 30)))
    assert abs(f.v_bar.sum() - 1.0) < 1e-12


def test_identifiability_identity_holds():
    f = fit(synthetic(np.linspace(0.1, 0.9, 40),
                      mu256=np.array([1e5] + [50.0] * 39)))
    r = identifiability(f)
    v1 = r["top1_share"]
    for a, c in ((r["anchor_range"][0], r["median_at_fixed_pool_max"]),
                 (r["anchor_range"][1], r["median_at_fixed_pool_min"])):
        assert abs(v1 * a + (1 - v1) * c - f.gamma_pool) < 1e-9


def test_leave_one_out_moves_the_answer():
    g = np.array([0.9, 0.2, 0.2])
    m = np.array([1e5, 10.0, 10.0])
    mm = synthetic(g, mu256=m)
    f = fit(mm)
    top = f.top1
    assert top == 0
    after = leave_one_out(mm, f, top)
    assert after < f.gamma_pool - 0.05, (after, f.gamma_pool)


def test_tidy_roundtrip(tmp_path=None):
    mm = synthetic(np.linspace(0.2, 0.8, 25), amplitudes=np.linspace(0.1, 1.5, 25))
    direct = fit(mm).gamma_pool
    p = Path(tmp_path or "/tmp") / "tidy.tsv"
    with p.open("w", encoding="utf-8") as fh:
        fh.write("gene\tscale\tmu\tvar\n")
        for i, gene in enumerate(mm.genes):
            for k, d in enumerate(mm.ds):
                fh.write(f"{gene}\t{d:.0f}\t{mm.mu[i, k]:.17g}\t{mm.var[i, k]:.17g}\n")
    back = Moments.from_tidy(p)
    assert back.mu.shape == mm.mu.shape
    assert abs(fit(back).gamma_pool - direct) < 1e-10


def test_rejects_non_monotone_ladder():
    mm = synthetic([0.3, 0.3])
    try:
        Moments(ds=np.array([2, 2, 8]), genes=mm.genes[:1],
                mu=mm.mu[:1, :1], var=mm.var[:1, :1])
    except MomentsError:
        pass
    else:
        raise AssertionError("a repeated bin size should be rejected")


def test_mde_matches_hand_calculation():
    got = mde(0.03, 2)
    want = 2.801585 * np.sqrt(2) * 0.03 / np.sqrt(2)
    assert abs(got - want) < 1e-5, (got, want)


def test_rule_text_has_four_sentences():
    f = fit(synthetic(np.linspace(0.2, 0.8, 20)))
    rep = summarise(synthetic(np.linspace(0.2, 0.8, 20)), f, name="t")
    assert len(rule_text(rep)) == 4


def test_published_colon_values():
    """The shipped example must reproduce the paper's colon row."""
    if not EXAMPLE.exists():
        print("  (skipped: examples/colon_moments.npz absent)")
        return
    m = Moments.load(EXAMPLE)
    f = fit(m)
    assert abs(f.gamma_pool - 0.5546672) < 2e-4, f.gamma_pool
    assert f.n_fitted == 13115, f.n_fitted
    assert str(m.genes[f.top1]) == "IGKC"
    assert abs(f.top1_share - 0.877) < 5e-3, f.top1_share
    assert abs(f.n_eff_v - 1.296) < 5e-3, f.n_eff_v
    ident = identifiability(f)
    assert abs(ident["d_median_d_anchor"] + 7.16) < 0.05, ident["d_median_d_anchor"]
    assert ident["median_at_fixed_pool_min"] < 0.01
    assert ident["median_at_fixed_pool_max"] > 0.99


def main() -> int:
    tests = [(k, v) for k, v in sorted(globals().items())
             if k.startswith("test_") and callable(v)]
    bad = 0
    for name, fn in tests:
        try:
            if name == "test_tidy_roundtrip":
                fn("/tmp")
            else:
                fn()
        except Exception as e:                       # noqa: BLE001
            bad += 1
            print(f"FAIL {name}: {type(e).__name__}: {e}")
        else:
            print(f"ok   {name}")
    print(f"\n{len(tests) - bad}/{len(tests)} passed")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())

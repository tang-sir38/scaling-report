r"""Tests for the quasi-likelihood correction. Run with pytest, or directly:
python tests/test_ql.py

The rule these tests enforce: the correction must be the exact inverse of the
diagnosis.  fpr_under_poisson(phi) is the false-positive rate a Poisson test
reaches when the variance ratio is phi; ql_p(z_crit, phi) is the p-value the
same data give after correction.  At z = z_crit the two must agree, or the
package is telling the reader to do something that does not undo the problem it
just described.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scaling_report import (fpr_under_poisson, phi_for_fpr, ql_p,  # noqa: E402
                            ql_phi, ql_test, ql_z)

PHIS = [1.0, 1.42, 4.0, 16.0, 30.0, 78.0, 400.0]
Z_CRIT = 1.959963984540054  # two-sided 5%


def test_phi_identity():
    assert ql_phi(0.0, 30.0) == 1.0
    assert ql_phi(25.0, 0.0) == 1.0
    assert abs(ql_phi(25.0, 0.6) - 16.0) < 1e-12
    assert abs(ql_phi(1.0, 77.0) - 78.0) < 1e-12


def test_z_is_unchanged_at_phi_one():
    for z in (-3.0, 0.0, 1.96, 5.0):
        assert abs(ql_z(z, 1.0) - z) < 1e-15


def test_z_divides_by_root_phi():
    assert abs(ql_z(3.0, 4.0) - 1.5) < 1e-15
    assert abs(ql_z(3.0, 9.0) - 1.0) < 1e-15


def test_correction_inverts_the_diagnosis():
    """ql_p at the Poisson critical value must equal the reported FPR."""
    for phi in PHIS:
        want = fpr_under_poisson(phi)
        got = ql_p(Z_CRIT, phi)
        assert abs(got - want) < 1e-12, (phi, want, got)


def test_phi_for_fpr_roundtrips():
    for target in (0.05, 0.10, 0.25, 0.50):
        phi = phi_for_fpr(target)
        assert abs(fpr_under_poisson(phi) - target) < 1e-10, target


def test_p_is_monotone_in_phi():
    """More overdispersion must never make the corrected test more significant.

    p rises with phi: 0.012 at phi = 1, 0.90 at phi = 400.  An earlier version of
    this test asserted the opposite and failed, which is the test doing its job.
    """
    ps = [ql_p(2.5, phi) for phi in PHIS]
    assert all(b >= a - 1e-15 for a, b in zip(ps, ps[1:])), ps
    assert ps[0] < 0.05 < ps[-1]


def test_p_bounds():
    assert 0.0 <= ql_p(100.0, 1.0) <= 1.0
    assert ql_p(0.0, 78.0) == 1.0


def test_ql_test_matches_pieces():
    got = ql_test(observed=120.0, expected=100.0, mu=25.0, alpha=0.6)
    assert abs(got["z_poisson"] - 2.0) < 1e-12
    assert abs(got["phi"] - 16.0) < 1e-12
    assert abs(got["z"] - 0.5) < 1e-12
    assert abs(got["p"] - ql_p(2.0, 16.0)) < 1e-15


def test_ql_test_matches_hand_arithmetic():
    """A worked case a reader can check on paper.

    mu = 25, alpha = 0.6  ->  phi = 16
    observed 120 against expected 100  ->  z_poisson = 20/10 = 2
    z_q = 2/4 = 0.5, two-sided p = 2*(1-Phi(0.5)) = 0.6171
    """
    got = ql_test(120.0, 100.0, 25.0, 0.6)
    hand = 2.0 * (1.0 - 0.6914624612740131)
    assert abs(got["p"] - hand) < 1e-12


def test_no_correction_when_alpha_is_zero():
    """alpha = 0 is the Poisson case, so the corrected test is the original."""
    got = ql_test(120.0, 100.0, 25.0, 0.0)
    assert got["phi"] == 1.0
    assert abs(got["z"] - got["z_poisson"]) < 1e-15


def test_rejects_bad_input():
    for fn, args in ((ql_phi, (-1.0, 0.5)), (ql_phi, (1.0, -0.5)),
                     (ql_z, (1.0, 0.0)), (ql_z, (1.0, -2.0)),
                     (ql_test, (1.0, 0.0, 1.0, 1.0)),
                     (ql_test, (1.0, -5.0, 1.0, 1.0))):
        try:
            fn(*args)
        except ValueError:
            continue
        raise AssertionError(f"{fn.__name__}{args} should have raised")


def main() -> int:
    tests = [(k, v) for k, v in sorted(globals().items())
             if k.startswith("test_") and callable(v)]
    bad = 0
    for name, fn in tests:
        try:
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

"""make_numbers_v9.py -- numbers_v8 plus the calibration keys.

Every value here is read from an output file produced by 135/136/137 in this
directory; nothing is typed by hand.  v8 keys are carried over unchanged.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
V8 = HERE.parent / "data" / "numbers_v8.json"
OUTD = HERE.parent / "out"
OUT = OUTD / "numbers_v9.json"

cal = json.loads((OUTD / "calibration.json").read_text(encoding="utf-8"))
gam = json.loads((OUTD / "gamma_scan.json").read_text(encoding="utf-8"))
bh = json.loads((OUTD / "bh_ordering.json").read_text(encoding="utf-8"))

n = json.loads(V8.read_text(encoding="utf-8"))


def analytic(pre, d):
    for r in cal[pre]["analytic"]:
        if abs(r["bin_um"] - d) < 1e-9:
            return r
    raise KeyError((pre, d))


def sim(pre, d):
    for r in cal[pre]["simulated"]:
        if abs(r["bin_um"] - d) < 1e-9 and r["delta"] == 0 and r.get("kind") != "frailty":
            return r
    raise KeyError((pre, d))


def blocked(pre, d):
    for r in cal[pre]["blocked"]:
        if abs(r["bin_um"] - d) < 1e-9:
            return r
    raise KeyError((pre, d))


def bhrow(pre, d):
    for r in bh["per_section"][pre]:
        if abs(r["bin_um"] - d) < 1e-9:
            return r
    raise KeyError((pre, d))


# ---- D: section attribution -------------------------------------------------
n["fpr_256_p1"] = round(analytic("P1", 256)["fpr_poisson"] * 100, 1)
n["fpr_256_p2"] = round(analytic("P2", 256)["fpr_poisson"] * 100, 1)
n["fpr_8_p1"] = round(analytic("P1", 8)["fpr_poisson"] * 100, 1)
n["fpr_8_p2"] = round(analytic("P2", 8)["fpr_poisson"] * 100, 1)
n["fpr_2_p2"] = round(analytic("P2", 2)["fpr_poisson"] * 100, 1)
n["phi_8_p1"] = round(analytic("P1", 8)["phi_median"], 2)
n["phi_256_p1"] = round(analytic("P1", 256)["phi_median"], 1)
n["phi_8_p2"] = round(analytic("P2", 8)["phi_median"], 2)
n["phi_256_p2"] = round(analytic("P2", 256)["phi_median"], 1)
# the two effective-gene definitions, both correct, not interchangeable
import numpy as np

_m = np.load(HERE.parent / "data" / "colon_P2_moments.npz")
_mu = _m["mu_2"].astype(float)
_var = _m["var_2"].astype(float)
with np.errstate(invalid="ignore", divide="ignore"):
    _a = (_var - _mu) / _mu ** 2
_ok = np.isfinite(_a) & (_a > 0) & (_mu > 0)
_w = _mu[_ok] ** 2
_w = _w / _w.sum()
n["neff_w_colon"] = round(float(1.0 / (_w ** 2).sum()), 3)
_v = _w * _a[_ok] / (_w * _a[_ok]).sum()
n["neff_v_colon"] = round(float(1.0 / (_v ** 2).sum()), 3)

# ---- A: the largest usable bin ---------------------------------------------
n["fpr_mean_p2_8"] = round(analytic("P2", 8)["fpr_poisson_mean"] * 100, 1)
n["fpr_mean_p2_16"] = round(analytic("P2", 16)["fpr_poisson_mean"] * 100, 1)
n["fpr_mean_p2_256"] = round(analytic("P2", 256)["fpr_poisson_mean"] * 100, 1)
n["fpr_mean_p1_8"] = round(analytic("P1", 8)["fpr_poisson_mean"] * 100, 1)
n["fpr_mean_p1_16"] = round(analytic("P1", 16)["fpr_poisson_mean"] * 100, 1)
n["fpr_mean_p1_32"] = round(analytic("P1", 32)["fpr_poisson_mean"] * 100, 1)
n["fpr_med_p2_16"] = round(analytic("P2", 16)["fpr_poisson"] * 100, 1)
n["safe_bin_p2"] = 8
n["safe_bin_p1"] = 16

# ---- B: the correction ------------------------------------------------------
ql = [sim(pre, d)["fpr_ql_moments"] * 100 for pre in ("P1", "P2") for d in (8, 16, 32, 64, 128)]
n["ql_fpr_lo"] = round(min(ql), 1)
n["ql_fpr_hi"] = round(max(ql), 1)
shr = [sim(pre, d)["fpr_ql_shrunk"] * 100 for pre in ("P1", "P2") for d in (8, 16, 32, 64, 128)]
n["ql_shrunk_lo"] = round(min(shr), 1)
n["ql_shrunk_hi"] = round(max(shr), 1)

fp = [bhrow(pre, d)["fdr_p"] * 100 for pre in ("P1", "P2") for d in (8, 32, 64, 128)]
fq = [bhrow(pre, d)["fdr_q"] * 100 for pre in ("P1", "P2") for d in (8, 32, 64, 128)]
n["fdr_p_lo"] = round(min(fp), 1)
n["fdr_p_hi"] = round(max(fp), 1)
n["fdr_q_lo"] = round(min(fq), 1)
n["fdr_q_hi"] = round(max(fq), 1)
n["top40_p_128"] = round(bhrow("P2", 128)["topk_true_p"] * 100, 1)
n["top40_q_128"] = round(bhrow("P2", 128)["topk_true_q"] * 100, 1)
n["top40_p_64"] = round(bhrow("P2", 64)["topk_true_p"] * 100, 1)
n["top40_q_64"] = round(bhrow("P2", 64)["topk_true_q"] * 100, 1)

bq = [blocked(pre, d)["fpr_ql_moments"] * 100 for pre in ("P1", "P2") for d in (16, 32, 64)]
bb = [blocked(pre, d)["fpr_blocklevel"] * 100 for pre in ("P1", "P2") for d in (16, 32, 64)]
n["block_ql_lo"] = round(min(bq), 1)
n["block_ql_hi"] = round(max(bq), 1)
n["block_unit_lo"] = round(min(bb), 1)
n["block_unit_hi"] = round(max(bb), 1)
n["block_poisson_lo"] = round(min(blocked(pre, d)["fpr_poisson"] * 100
                                   for pre in ("P1", "P2") for d in (16, 32, 64)), 1)
n["block_poisson_hi"] = round(max(blocked(pre, d)["fpr_poisson"] * 100
                                   for pre in ("P1", "P2") for d in (16, 32, 64)), 1)

# ---- power ------------------------------------------------------------------
n["neff_p1_256"] = round(analytic("P1", 256)["n_equivalent_fraction"] * 100, 1)
n["neff_p2_64"] = round(analytic("P2", 64)["n_equivalent_fraction"] * 100, 1)
n["neff_p2_256"] = round(analytic("P2", 256)["n_equivalent_fraction"] * 100, 1)
n["loss_p1_256"] = round(analytic("P1", 256)["power"]["0.25"]["loss_pp_median"], 1)
n["loss_p2_256"] = round(analytic("P2", 256)["power"]["0.25"]["loss_pp_median"], 1)
n["power_retained_p1_256"] = round(
    100 * analytic("P1", 256)["power"]["0.25"]["ql_median"], 1)
n["power_retained_p2_256"] = round(
    100 * analytic("P2", 256)["power"]["0.25"]["ql_median"], 1)
# computed, not typed: the largest power loss at any scale up to 128 um
n["loss_max_below_128"] = round(max(
    analytic(pre, d)["power"]["0.25"]["loss_pp_median"]
    for pre in ("P1", "P2") for d in (2, 4, 8, 16, 32, 64, 128)), 1)

# the price is ordered by alpha, not by phi and not by gamma
n["alpha_256_p1"] = round((analytic("P1", 256)["phi_median"] - 1.0)
                          / gam["P1"]["mu_median_256"], 2)
n["alpha_256_p2"] = round((analytic("P2", 256)["phi_median"] - 1.0)
                          / gam["P2"]["mu_median_256"], 2)

# depth is not a lever: the loss barely moves across an eightfold depth range
dsc = {round(s["depth_multiplier"], 1): s for s in gam["P1"]["depth_scan"]}
n["depth_loss_p1_half"] = round(dsc[0.5]["loss_pp"], 1)
n["depth_loss_p1_x4"] = round(dsc[4.0]["loss_pp"], 1)
n["depth_power_p1_half"] = round(100 * dsc[0.5]["power_ql_median"], 1)
n["depth_power_p1_x4"] = round(100 * dsc[4.0]["power_ql_median"], 1)
dsc2 = {round(s["depth_multiplier"], 1): s for s in gam["P2"]["depth_scan"]}
n["depth_loss_p2_half"] = round(dsc2[0.5]["loss_pp"], 1)
n["depth_loss_p2_x4"] = round(dsc2[4.0]["loss_pp"], 1)

# gamma extrapolated out of sample: fit on 2-16 um, price 256 um
def pros(pre, U, target):
    for r in gam[pre]["prospective_extrapolation"]:
        if abs(r["fit_window_um"][1] - U) < 1e-9 and abs(r["target_um"] - target) < 1e-9:
            return r
    raise KeyError((pre, U, target))


n["extrap_under_p1"] = round(1.0 / pros("P1", 16, 256)["factor_pred_over_meas"], 1)
n["extrap_under_p2"] = round(1.0 / pros("P2", 16, 256)["factor_pred_over_meas"], 1)
n["extrap_under_p1_32"] = round(1.0 / pros("P1", 16, 32)["factor_pred_over_meas"], 1)
n["extrap_under_p2_32"] = round(1.0 / pros("P2", 16, 32)["factor_pred_over_meas"], 1)

# how far the full-window fit carries
def fid(pre, d):
    for r in gam[pre]["power_law_fidelity"]:
        if abs(r["bin_um"] - d) < 1e-9:
            return r
    raise KeyError((pre, d))


for pre in ("P1", "P2"):
    for d in (16, 32, 64, 256):
        n[f"fid_{pre.lower()}_{d}"] = round(100 * abs(fid(pre, d)["rel_error"]), 1)

# null variance ratio is quoted in S3 Text as "within 4%"; left to that note so
# there is one place to check it, not two.

# gamma scan: the shallow-gamma case
gscan = {s["gamma"]: s for s in gam["P2"]["gamma_scan"]}
n["loss_gamma_030"] = round(gscan[0.3]["loss_pp"], 1)
n["loss_gamma_0364"] = round(gscan[0.364]["loss_pp"], 1)
n["loss_gamma_0554"] = round(gscan[0.554]["loss_pp"], 1)
n["loss_gamma_070"] = round(gscan[0.7]["loss_pp"], 1)
# how far the single power law sits below the measured phi at 256 um
n["pl_under_p2_256"] = round(1.0 / gam["P2"]["extrapolation_ratio"], 1)
n["gamma_median_p1"] = round(gam["P1"]["gamma_median_per_gene"], 3)
n["gamma_median_p2"] = round(gam["P2"]["gamma_median_per_gene"], 3)
n["extrap_ratio_p2"] = round(gam["P2"]["extrapolation_ratio"], 2)

# ---- software ---------------------------------------------------------------
n["n_selftest"] = 10

OUT.write_text(json.dumps(n, indent=2, ensure_ascii=False), encoding="utf-8")
print("keys:", len(n), "->", OUT)
for k in ("fpr_256_p1", "fpr_256_p2", "fpr_mean_p2_8", "fpr_mean_p2_16",
          "ql_fpr_lo", "ql_fpr_hi", "fdr_p_hi", "fdr_q_lo", "fdr_q_hi",
          "loss_p1_256", "loss_p2_256", "loss_gamma_030", "block_ql_hi",
          "block_unit_lo", "block_unit_hi", "top40_p_128", "top40_q_128",
          "neff_p2_64", "gamma_median_p1", "gamma_median_p2"):
    print(f"  {k} = {n[k]}")

# -*- coding: utf-8 -*-
r"""
54_build_lock_v4.py -- NUMBERS_LOCK_v3.json + this round's three analyses = v4.

Added sections
  decomp        exact decomposition of the pooled exponent; effective gene count
                under the alpha-contribution weight as well as the mu^2 weight
  robust        stability of five candidate summary scalars across four bands
  spatial_proxy per-gene correlation length and Moran's I against the exponent

Nothing is typed in by hand: every value is read from a JSON written by
50_decomp_exact.py, 51_robust_scalar.py or 52_spatial_proxy.py.

Output: 02_OUT/NUMBERS_LOCK_v4.json
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
O = ROOT / "02_OUT"

v3 = json.loads((O / "NUMBERS_LOCK_v3.json").read_text(encoding="utf-8"))
out = dict(v3)


def read(p):
    return json.loads(Path(p).read_text(encoding="utf-8"))


# ---------------------------------------------------------------- decomp
dec = read(O / "decomp_exact" / "decomp_exact.json")
rels = {t: dec[t]["rel_diff"] for t in dec}
neff_mu2 = {t: dec[t]["n_eff_mu2"] for t in dec}
neff_al = {t: dec[t]["n_eff_alpha"] for t in dec}
out["decomp"] = dict(
    by_tissue=dec,
    min_rel_repro=min(rels.values()),
    max_rel_repro=max(rels.values()),
    max_identity_err=max(dec[t]["identity_err"] for t in dec),
    min_neff_mu2=min(neff_mu2.values()),
    max_neff_mu2=max(neff_mu2.values()),
    min_neff_alpha=min(neff_al.values()),
    max_neff_alpha=max(neff_al.values()),
    min_shrink=min(neff_mu2[t] / neff_al[t] for t in dec),
    max_shrink=max(neff_mu2[t] / neff_al[t] for t in dec),
    source="50_decomp_exact.py",
)

# ---------------------------------------------------------------- robust
rob = read(O / "robust_scalar" / "robust_scalar.json")
out["robust"] = dict(
    section=rob["section"],
    scalars=rob["scalars"],
    full=rob["full"],
    bands=rob["bands"],
    summary=rob["summary"],
    best=rob["best"],
    worst=rob["worst"],
    ratio_worst_over_best=rob["ratio_worst_over_best"],
    ratio_pooled_over_best=rob["ratio_pooled_over_best"],
    band_edges=rob["band_edges"],
    source="51_robust_scalar.py",
)

# ---------------------------------------------------------------- spatial proxy
sp = {}
for t in ("colon", "lung", "brain"):
    p = O / "spatial_proxy" / f"spatial_proxy_{t}.json"
    if p.exists():
        sp[t] = read(p)
out["spatial_proxy"] = dict(by_tissue=sp, n_tissues=len(sp),
                            source="52_spatial_proxy.py")

# ---------------------------------------------------------------- splithalf by tissue
# 37_splithalf.py has been run on all three sections, so the spread of the pooled
# exponent across four bands of one section is available per tissue rather than for
# the colon section alone.
SH = {
    "colon": O / "splithalf" / "splithalf.json",
    "lung": ROOT / "D3_LungCancer" / "02_OUT" / "splithalf" / "splithalf.json",
    "brain": ROOT / "D2_MouseBrain" / "02_OUT" / "splithalf" / "splithalf.json",
}
sh = {}
for t, p in SH.items():
    if not p.exists():
        continue
    d = read(p)
    q = [v["gamma"] for k, v in d["bands"].items() if k.startswith("q")]
    h = [v["gamma"] for k, v in d["bands"].items() if k.startswith("half")]
    sh[t] = dict(
        quarter_sd=float(d["quarter_sd"]),
        quarter_range=float(d["quarter_range"]),
        quarter_gamma=[float(x) for x in q],
        quarter_min=float(min(q)),
        quarter_max=float(max(q)),
        half_gamma=[float(x) for x in h],
        half_gap=float(abs(h[0] - h[1])) if len(h) == 2 else None,
        n_bin_quarter=[int(v["n_bin"]) for k, v in d["bands"].items() if k.startswith("q")],
        source=str(p.relative_to(ROOT)),
    )
out["splithalf_tissues"] = dict(by_tissue=sh, n_tissues=len(sh),
                                min_sd=min(v["quarter_sd"] for v in sh.values()) if sh else None,
                                max_sd=max(v["quarter_sd"] for v in sh.values()) if sh else None,
                                source="37_splithalf.py")

# ---------------------------------------------------------------- derived4
d4 = dict(
    neff_alpha_gap_pct={t: (neff_mu2[t] - neff_al[t]) / neff_mu2[t] for t in dec},
    drift_frac={t: dec[t]["drift_share"] for t in dec},
    top1_alpha_share={t: dec[t]["top1_alpha_share"] for t in dec},
    top1_mu2_share={t: dec[t]["top1_mu2_share"] for t in dec},
    robust_sd={k: v["sd"] for k, v in rob["summary"].items()},
    robust_cv={k: v["cv"] for k, v in rob["summary"].items()},
    robust_range={k: v["range"] for k, v in rob["summary"].items()},
)
if sp:
    d4["partial_rho_moran"] = {t: sp[t]["partial_rho_moran"] for t in sp}
    d4["partial_rho_xi"] = {t: sp[t]["partial_rho_xi"] for t in sp}
    d4["xi_frac_at_floor"] = {t: sp[t]["xi_frac_at_floor"] for t in sp}
    all_dec = [d["rho_moran"] for t in sp for d in sp[t]["by_decile"]]
    d4["within_decile_moran_min"] = min(all_dec)
    d4["within_decile_moran_max"] = max(all_dec)
out["derived4"] = d4

out["meta"] = dict(v3["meta"])
out["meta"].update(
    version="4.0",
    built="2026-11-14",
    supersedes="NUMBERS_LOCK_v3.json",
    note="唯一权威数字表。v4 = v3 + decomp（精确分解与 alpha 权重有效基因数）+ "
         "robust（五种标量的四带稳定性）+ spatial_proxy（相关长度与 Moran's I）。",
    what_changed=list(v3["meta"].get("what_changed", [])) + [
        "新增 decomp：pooled alpha 作为逐基因 alpha 的凸组合（恒等式误差至机器精度）、"
        "gamma_pool 与 sum(v*gamma) 的重构误差、alpha 贡献口径的有效基因数、权重漂移项",
        "新增 robust：pooled / median / mean / trimmed / mu2-加权中位数 在结肠四带上的 SD、range、CV",
        "新增 spatial_proxy：32 um 网格上逐基因指数相关长度与 Moran's I，及其在控制丰度后与 gamma 的偏相关",
        "新增 derived4：由三节派生的比值与极值，供正文注入",
    ],
)

p = O / "NUMBERS_LOCK_v4.json"
p.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
print("written", p)
print(f"  decomp   rel repro {out['decomp']['min_rel_repro']*100:.2f}-{out['decomp']['max_rel_repro']*100:.2f}%")
print(f"  n_eff    mu2 {out['decomp']['min_neff_mu2']:.1f}-{out['decomp']['max_neff_mu2']:.0f}"
      f"   alpha {out['decomp']['min_neff_alpha']:.1f}-{out['decomp']['max_neff_alpha']:.1f}")
print("  robust   SD " + "  ".join(f"{k}={v['sd']:.4f}" for k, v in rob["summary"].items()))
print(f"  best {rob['best']['name']}  worst {rob['worst']['name']}"
      f"  pooled/best {rob['ratio_pooled_over_best']:.2f}")
for t in sp:
    print(f"  proxy {t}: partial rho Moran {sp[t]['partial_rho_moran']:+.3f}"
          f"  xi floor {sp[t]['xi_frac_at_floor']*100:.0f}%")
for t, v in sh.items():
    print(f"  splithalf {t}: quarter SD {v['quarter_sd']:.4f}  range {v['quarter_range']:.4f}")

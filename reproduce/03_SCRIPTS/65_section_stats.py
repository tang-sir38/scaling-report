# -*- coding: utf-8 -*-
r"""
65_section_stats.py -- section-level inference for the replication analysis.

Reads
  02_OUT/multisection/pool.json              point estimates for four sections
  02_OUT/block_boot/boot_{colon,lung,brain}.json
  MS_..._P1/02_OUT/block_boot/boot_colonP1.json

Writes
  02_OUT/multisection/section_stats.json

Why a separate script: the point estimate of gamma_pool and the mean of the
superblock bootstrap replicates are not the same number (0.5547 vs 0.5559 for
colon). The bootstrap exists to give a standard error, not to replace the
estimate. Everything below therefore takes the estimate from pool.json and the
SE from the 512 um superblock bootstrap, and never mixes the two within one
contrast.
"""
from __future__ import annotations

import json
from math import sqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
O = ROOT / "02_OUT"
P1DIR = ROOT / "MS_Visium_HD_Human_Colon_Cancer_P1" / "02_OUT"

BOOT = {
    "Visium_HD_Human_Colon_Cancer_P1": P1DIR / "block_boot" / "boot_colonP1.json",
    "colon": O / "block_boot" / "boot_colon.json",
    "lung": O / "block_boot" / "boot_lung.json",
    "brain": O / "block_boot" / "boot_brain.json",
}
BLOCK = "512"          # primary block size, frozen before this analysis
N_BOOT = 200
Z_ALPHA = 1.959964     # two-sided 5%
Z_POWER = 0.841621     # 80%


def read(p):
    return json.loads(Path(p).read_text(encoding="utf-8"))


pool = read(O / "multisection" / "pool.json")
rows = {r["label"]: r for r in pool["rows"]}

sec = {}
for label, p in BOOT.items():
    b = read(p)["results"][BLOCK]
    r = rows[label]
    sec[label] = dict(
        tissue=r["tissue"],
        gamma_pool=float(r["gamma_pool"]),
        gamma_median=float(r["gamma_median"]),
        gamma_q25=float(r["gamma_q25"]),
        gamma_q75=float(r["gamma_q75"]),
        r2=float(r["r2"]),
        n_genes=int(r["n_genes"]),
        n_fit=int(r["n_fit"]),
        n_bins=None,          # filled below when the section metrics file is found
        umi_total=float(r["umi_total"]),
        top1_gene=r["top1_gene"],
        top1_alpha_share=float(r["top1_alpha_share"]),
        n_eff_alpha=float(r["n_eff_alpha"]),
        n_eff_mu2=float(r["n_eff_mu2"]),
        quarter_sd=float(r["quarter_sd"]),
        boot_mean=float(b["mean"]),
        boot_sd=float(b["sd"]),
        boot_lo=float(b["q025"]),
        boot_hi=float(b["q975"]),
        boot_nblk=int(b["nblk"]),
    )

# in-tissue bin counts at 2 um: three sections carry them in the v4 lock, the
# repeat section does not, so it is left None rather than guessed.
_lock4 = read(O / "NUMBERS_LOCK_v4.json")
for label, key in (("colon", "colon"), ("lung", "lung"), ("brain", "brain")):
    if key in _lock4.get("datasets", {}):
        sec[label]["n_bins"] = int(_lock4["datasets"][key]["n_bins"])

# ------------------------------------------------------------------ contrasts
def diff(a, b):
    """difference of two independent section estimates; SE from the bootstrap"""
    da = sec[a]["gamma_pool"] - sec[b]["gamma_pool"]
    se = sqrt(sec[a]["boot_sd"] ** 2 + sec[b]["boot_sd"] ** 2)
    return dict(
        a=a, b=b,
        delta=float(da), se_diff=float(se), z=float(da / se),
        ci_lo=float(da - Z_ALPHA * se), ci_hi=float(da + Z_ALPHA * se),
    )


contrasts = dict(
    same_tissue_colon=diff("Visium_HD_Human_Colon_Cancer_P1", "colon"),
    colon_vs_brain=diff("colon", "brain"),
    colon_vs_lung=diff("colon", "lung"),
    lung_vs_brain=diff("lung", "brain"),
)


def diff_regional(a, b):
    """the same contrast under the more conservative error model: the SD of the
    four regional bands inside each section, i.e. the spread the manuscript
    already uses for the design calculation, not the bootstrap SE"""
    da = sec[a]["gamma_pool"] - sec[b]["gamma_pool"]
    se = sqrt(sec[a]["quarter_sd"] ** 2 + sec[b]["quarter_sd"] ** 2)
    return dict(a=a, b=b, delta=float(da), se_diff=float(se), z=float(da / se))


contrasts_regional = {k: diff_regional(v["a"], v["b"]) for k, v in contrasts.items()}

# minimum detectable difference for a two-section design, using the typical
# single-section SE of the two colon sections (the tissue where k = 2)
se_typ = (sec["colon"]["boot_sd"] + sec["Visium_HD_Human_Colon_Cancer_P1"]["boot_sd"]) / 2
mde = (Z_ALPHA + Z_POWER) * sqrt(2) * se_typ

# spread of the four section estimates versus the spread of four bands inside
# one section -- the comparison that says what n = 1 per tissue actually costs
gamma_pool_all = [sec[k]["gamma_pool"] for k in ("colon", "lung", "brain")]
gamma_pool_all.sort()

out = dict(
    block_um=int(BLOCK), n_boot=N_BOOT,
    sections=sec,
    contrasts=contrasts,
    contrasts_regional=contrasts_regional,
    mde=dict(
        se_typ=float(se_typ), value=float(mde),
        z_alpha=Z_ALPHA, z_power=Z_POWER,
        z_used=Z_ALPHA + Z_POWER,
        note="two-section design, 80% power, two-sided 5%",
    ),
    spread=dict(
        colon_P1_minus_colon=float(sec["Visium_HD_Human_Colon_Cancer_P1"]["gamma_pool"]
                                   - sec["colon"]["gamma_pool"]),
        brain_to_colon=float(max(gamma_pool_all) - min(gamma_pool_all)),
        within_section_quarter_sd=[sec[k]["quarter_sd"] for k in
                                   ("Visium_HD_Human_Colon_Cancer_P1", "colon", "lung", "brain")],
        within_section_quarter_sd_mean=float(
            sum(sec[k]["quarter_sd"] for k in
                ("Visium_HD_Human_Colon_Cancer_P1", "colon", "lung", "brain")) / 4),
    ),
    regression=pool["regression"],
    icc_note=("ICC over four equal-count bands is reported by 58_pool_sections.py but is not used: "
              "the bands are non-overlapping spatial subregions, not exchangeable replicates, so the "
              "random-effects model behind ICC does not hold. The superblock bootstrap replaces it."),
    source="65_section_stats.py",
)

(O / "multisection" / "section_stats.json").write_text(
    json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")

print(json.dumps(dict(
    mde=out["mde"]["value"],
    same=out["contrasts"]["same_tissue_colon"],
    cb=out["contrasts"]["colon_vs_brain"],
    typ_se=se_typ,
), indent=1, ensure_ascii=False))

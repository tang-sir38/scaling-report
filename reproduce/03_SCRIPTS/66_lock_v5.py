# -*- coding: utf-8 -*-
r"""
66_lock_v5.py -- NUMBERS_LOCK_v4.json + cross-platform + section replication = v5.

Two blocks are added.

  cross   Stereo-seq mouse brain against the Visium HD mouse brain section:
          the frozen criteria (H1a/H1b/H2/H4), the two falsification attacks,
          and the per-gene replication correlation. Read from
          02_OUT/crossplatform/{crossplatform,robustness}.json.
  multi   four sections, the superblock bootstrap, and the contrasts between
          them. Read from 02_OUT/multisection/{pool,section_stats}.json.
          The contrasts are recomputed by 65_section_stats.py; nothing is
          taken from an earlier summary.

Nothing is typed in by hand.

Output: 02_OUT/NUMBERS_LOCK_v5.json
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
O = ROOT / "02_OUT"

v4 = json.loads((O / "NUMBERS_LOCK_v4.json").read_text(encoding="utf-8"))
out = dict(v4)


def read(p):
    return json.loads(Path(p).read_text(encoding="utf-8"))


# ------------------------------------------------------------------ cross
cp = read(O / "crossplatform" / "crossplatform.json")
rb = read(O / "crossplatform" / "robustness.json")
vd = cp["verdict"]

stereo = cp["stereo_brain"]
visium = cp["visium_hd_brain"]

out["cross"] = dict(
    stereo=dict(
        label="Stereo-seq C04042E2, mouse whole brain, bin1 (0.5 um DNB)",
        n_genes=int(stereo["n_genes"]),
        n_genes_ok=int(stereo["n_genes_ok"]),
        n_cells=[int(x) for x in stereo["n_cells"]],
        gamma_pool=float(stereo["gamma_pool"]),
        gamma_weighted=float(stereo["gamma_weighted"]),
        gamma_median=float(stereo["gamma_median"]),
        gamma_mean=float(stereo["gamma_mean"]),
        identity_err=float(stereo["identity_err"]),
        top1=stereo["top1"],
        top1_alpha_share=float(stereo["top1_alpha_share"]),
        top1_mu2_share=float(stereo["top1_mu2_share"]),
        top3=[str(x) for x in stereo["top3"]],
        n_eff_mu2=float(stereo["n_eff_mu2"]),
        n_eff_alpha=float(stereo["n_eff_alpha"]),
        alpha=[float(x) for x in stereo["alpha_pool_by_scale"]],
    ),
    visium=dict(
        label="Visium HD mouse brain, 2 um bins",
        n_genes=int(visium["n_genes"]),
        n_genes_ok=int(visium["n_genes_ok"]),
        gamma_pool=float(visium["gamma_pool"]),
        gamma_median=float(visium["gamma_median"]),
        top1=visium["top1"],
        top1_alpha_share=float(visium["top1_alpha_share"]),
        n_eff_alpha=float(visium["n_eff_alpha"]),
        alpha=[float(x) for x in visium["alpha_pool_by_scale"]],
    ),
    delta_gamma_pool=float(vd["H3_level"]["delta_gamma_pool"]),
    delta_gamma_median=float(vd["H3_level"]["delta_gamma_median"]),
    concentration=dict(
        n_eff_over_genes=float(vd["H1a_concentration"]["n_eff_over_genes"]),
        h1a_rule=vd["H1a_concentration"]["rule"],
        h1a_pass=bool(vd["H1a_concentration"]["pass_"]),
        h1b_rule=vd["H1b_dominance"]["rule"],
        h1b_pass=bool(vd["H1b_dominance"]["pass_"]),
        h2_rule=vd["H2_ordering"]["rule"],
        h2_pass=bool(vd["H2_ordering"]["pass_"]),
    ),
    without_rrna=dict(
        genes=int(rb["stereo_without_rrna"]["genes"]),
        top1=rb["stereo_without_rrna"]["top1"],
        top1_share=float(rb["stereo_without_rrna"]["top1_share"]),
        n_eff_alpha=float(rb["stereo_without_rrna"]["n_eff_alpha"]),
        h1a_pass=bool(rb["stereo_without_rrna"]["h1a"]),
        h1b_pass=bool(rb["stereo_without_rrna"]["h1b"]),
    ),
    replication=dict(
        n=int(rb["replication"]["n"]),
        rho_raw=float(rb["replication"]["rho_raw"]),
        rho_partial_abundance=float(rb["replication"]["rho_partial_abundance"]),
        rho_abundance=float(rb["replication"]["rho_abundance"]),
        top100_stereo=float(rb["replication"]["top100_stereo"]),
        top100_visium=float(rb["replication"]["top100_visium"]),
        top100_stereo_without_rrna=float(rb["replication"]["top100_stereo_without_rrna"]),
        h4_rule=vd["H4_replication"]["rule"],
        h4_threshold=0.40,   # frozen criterion, not a measurement
        h4_raw_pass=bool(rb["replication"]["h4_raw"]),
        h4_partial_pass=bool(rb["replication"]["h4_partial"]),
    ),
    top100_share_stereo=float(rb["replication"]["top100_stereo"]),
    top100_share_visium=float(rb["replication"]["top100_visium"]),
    n_shared_cells_max=int(max(stereo["n_cells"])),
    source="63_crossplatform.py, 64_crossplatform_robust.py",
)

# ------------------------------------------------------------------ multi
pool = read(O / "multisection" / "pool.json")
ss = read(O / "multisection" / "section_stats.json")
pg = read(O / "multisection" / "pergene_shift.json")

out["multi"] = dict(
    n_sections=int(pool["n_sections"]),
    sections=ss["sections"],
    contrasts=ss["contrasts"],
    contrasts_regional=ss["contrasts_regional"],
    mde=ss["mde"],
    spread=ss["spread"],
    regression=ss["regression"],
    block_um=ss["block_um"],
    n_boot=ss["n_boot"],
    icc=pool["components"],
    icc_note=ss["icc_note"],
    pergene=dict(
        deciles=pg["deciles"],
        shift_overall=float(pg["shift_overall"]),
        shift_top=float(pg["shift_top_decile"]),
        shift_bottom=float(pg["shift_bottom_decile"]),
        n_deciles_ci_excludes_zero=int(pg["n_deciles_ci_excludes_zero"]),
        monotone_decreasing=bool(pg["monotone_decreasing"]),
        verdict=pg["verdict"],
        median_p0=float(pg["sections"]["colon_P0"]["median"]),
        median_p1=float(pg["sections"]["colon_P1"]["median"]),
        median_ci_p0=pg["sections"]["colon_P0"]["median_ci"],
        median_ci_p1=pg["sections"]["colon_P1"]["median_ci"],
        umi_p0=float(pg["sections"]["colon_P0"]["umi_total"]),
        umi_p1=float(pg["sections"]["colon_P1"]["umi_total"]),
        n_fit_p0=int(pg["sections"]["colon_P0"]["n_fit"]),
        n_fit_p1=int(pg["sections"]["colon_P1"]["n_fit"]),
        source="68_section_pergene.py",
    ),
    gamma_pool_min=min(v["gamma_pool"] for v in ss["sections"].values()),
    gamma_pool_max=max(v["gamma_pool"] for v in ss["sections"].values()),
    se_min=min(v["boot_sd"] for v in ss["sections"].values()),
    se_max=max(v["boot_sd"] for v in ss["sections"].values()),
    n_pairs=len(ss["contrasts"]),
    source="65_section_stats.py, 58_pool_sections.py, 68_section_pergene.py",
)

(O / "NUMBERS_LOCK_v5.json").write_text(
    json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")

print("blocks:", len(out), "keys; file chars:",
      len(json.dumps(out, ensure_ascii=False)))
print("cross h1a/h1b/h2/h4raw/h4partial:",
      out["cross"]["concentration"]["h1a_pass"],
      out["cross"]["concentration"]["h1b_pass"],
      out["cross"]["concentration"]["h2_pass"],
      out["cross"]["replication"]["h4_raw_pass"],
      out["cross"]["replication"]["h4_partial_pass"])
print("multi gamma_pool:",
      {k: round(v["gamma_pool"], 4) for k, v in out["multi"]["sections"].items()})
print("mde:", round(out["multi"]["mde"]["value"], 4))

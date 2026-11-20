# Reproducing the analysis

This directory is the working tree in which every number in the manuscript was produced.
It keeps the original layout, so the scripts run without editing: every path inside them is
resolved as `Path(__file__).parents[1]`, which is this directory.

```
reproduce/
├── 01_DATA/        where the downloaded data go (see 01_DATA/README.md)
├── 02_OUT/         derived tables, the numbers lock, the guard output, and
│                   REFS_v4.json, the keyed bibliography the manuscript is built from
├── 03_SCRIPTS/     138 scripts, numbered in the order they were run
├── 07_SUBMISSION/  MANUSCRIPT_v4.md, the template that 98_build_v4.py fills with numbers
├── D2_MouseBrain/                      derived tables for the mouse brain section
├── D3_LungCancer/                      derived tables for the lung section
└── MS_Visium_HD_Human_Colon_Cancer_P1/ derived tables for the second colon section
```

`97_lock_v7.py` reads all four sections, so the three per-section directories have to stay
in place even though the colon section — the one the paper leads with — lives in `02_OUT/`.

Comments and docstrings in the scripts are in Chinese, because that is how they were
written. The first line of each docstring states the question the script answers, in the
form of a claim that can be wrong.

## Two tiers of input

| tier | what | where |
|---|---|---|
| **A** | per-gene moments, per-gene curves, covariogram exponents, split-half, cross-platform and section tables, the numbers lock | shipped here, in `02_OUT/`, 19 MB |
| **B** | the binned count matrices `01_WORK/bins_002um.npz` … `bins_256um.npz` | rebuilt from public data, 41 GB, see `01_DATA/README.md` |

Tier A is enough to re-derive **every number printed in the paper**. Tier B is needed only
to re-derive Tier A itself.

## Re-run

```bash
pip install -r ../requirements.txt
cd 03_SCRIPTS

python 97_lock_v7.py       # -> 02_OUT/NUMBERS_LOCK_v7.json   (the single source of numbers)
python 98_build_v4.py      # -> 07_SUBMISSION/_build_v4/docs/01_Manuscript.docx
python 101_supp_v4.py      # -> 07_SUBMISSION/_build_v4/docs/02_Supplementary_Information.docx
python 100_guard_v4.py     # -> checks, then READY or BLOCKED
```

`100_guard_v4.py` is the last line of defence. It recomputes every derived number from the
lock instead of trusting the ledger, cross-checks the abstract against the body, checks that
the panel letters cited in the body exist in the figure legends, and refuses to print
`READY` if anything disagrees. A run that ends `BLOCKED` names the failing check and the
offending value. The guard is meant to fail: during preparation it caught a mislabelled
figure panel, a table window description that did not match the data, and three bugs in its
own comparison code.

## Which script produces what

**Data and moments**

| script | output |
|---|---|
| `02b_rebin.py` | exact integer re-aggregation 2 → 256 µm (tissue-aware), writes `01_WORK/bins_*.npz` |
| `13_gene_moments.py` | per-gene (μ, var) at every scale, `02_OUT/gene_moments.npz` |
| `20_table1_stats.py` | Table 1 |
| `41_gene_curves.py` | per-gene α_g(d) and γ_g, `gene_curves/` |

**The exponent itself**

| script | question |
|---|---|
| `07_pooled_alpha.py` | pooled α(d) without any gene filter — Fig. 1a |
| `09_pooled_local.py` | the same pooled quantity inside local windows |
| `36_jackknife.py` | leave-one-gene-out, cumulative and random top-k |
| `50_decomp_exact.py` | is γ_pool a convex combination of per-gene γ, arithmetically? two effective-gene counts |
| `51_robust_scalar.py` | five scalars × four column bands, SD and CV |
| `52_spatial_proxy.py` | is γ_g a proxy for clustering? semivariogram range and Moran's I |

**Controls**

| script | question |
|---|---|
| `32_synth_expression.py` | does γ move with abundance alone, on synthetic counts? |
| `33_depth_control.py` | binomial thinning to 46 %, 25 %, 10 % — does γ move? |
| `34_layer_variance.py` | two-way decomposition, expression layer vs tissue |
| `38_fill_sensitivity.py` | all-in-tissue vs full-block domain convention |
| `39_block_boot.py` | superblock bootstrap at 512 / 1024 / 2048 µm |

**Replication and identifiability**

| script | question |
|---|---|
| `37_splithalf.py`, `44_splithalf_genes.py`, `92_alpha_splithalf.py` | within-section split-half, per gene |
| `56_section_metrics.py`, `58_pool_sections.py`, `65_section_stats.py`, `68_section_pergene.py` | second section, section-level inference, why the per-gene median moves |
| `59_gef_probe.py`, `60_gef_pull.py`, `61_gef_moments.py` | Stereo-seq section, HTTP range pull of the expression table |
| `62_identifiability.py` | what γ_pool does **not** identify |
| `63_crossplatform.py`, `64_crossplatform_robust.py` | cross-platform test of the composition claim |

**Mechanism layer (added in this version)**

| script | output |
|---|---|
| `86_theory_alpha_cv.py` | numerical validation of α(A) = CV²₀ · ρ̄(A): synthetic fields, white noise, thinning — Fig. 3a–c |
| `87_variogram_shape.py` | per-gene covariogram exponent β_g, and the local slope of ρ(h) — Figs. 3d, 4a–b |
| `88_region_excise.py` | region excision, per-gene |
| `89_multigene_permute.py` | dose–response: permute the top N transcripts, N = 0…1000 — Fig. 5b, Table S4 |
| `90_window_consistency.py`, `91_window_scan.py`, `94_window_reliability.py` | the 6 × 6 window matrix and the reliability-only choice of window — Fig. 4c, Table S2 |
| `93_fig9_theory.py` | theory figure |
| `95_curve_shape.py`, `96_gamma_from_rhobar.py` | **two failed attempts** to predict the whole α(d) curve — Supplementary Note 1 |

**Assembly**

| script | output |
|---|---|
| `97_lock_v7.py` | `NUMBERS_LOCK_v7.json` — every printed number, recomputed from source, each block carrying the script that produces it |
| `43_figures_v2.py` | Figs. 1, 2 and 6, drawn from the lock |
| `99_fig45_v4.py` | Figs. 3, 4 and 5 |
| `98_build_v4.py` | the manuscript: `{{lock.path:format}}` injected, `[[refkey]]` numbered by first appearance |
| `101_supp_v4.py` | Supplementary Information, tables also taken from the lock |
| `100_guard_v4.py` | the guard |
| `106_downstream_v4.py` | the cost of getting the variance wrong: implied false-positive rate, differential expression and marker mapping across the bin ladder — last Results subsection, Fig. S2, Tables S7–S8 |

Earlier lock builders (`21_build_numbers.py`, `42_build_lock_v3.py`, `54_build_lock_v4.py`,
`66_lock_v5.py`), earlier guards (`70`–`72`) and the earlier manuscript builder
(`78_build_docs_v31.py`) are kept because each lock inherits the previous one verbatim and
the chain has to close. `97_lock_v7.py` is the only one that should be run for the final
numbers.

## Why two scripts in here produce no result

`95_curve_shape.py` and `96_gamma_from_rhobar.py` ask whether the measured covariogram
predicts the whole α(d) curve, not just its slope. It does not, and the reason is a property
of the data rather than of the theory: the covariogram is measured on the 8 µm grid, so its
first lag is already a block average of the 2 µm structure, and block-averaging a second
time overshoots the decay by a factor of 2 to 8. This is stated as a limitation in the main
text and documented in Supplementary Note 1. Both scripts run to completion and write their
outputs; they simply do not support the prediction.

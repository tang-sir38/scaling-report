# scaling_report

**This repository accompanies the manuscript** *Overdispersion in binned spatial
transcriptomics falls 2- to 15-fold, not 16,384-fold, because it is set by tissue
covariance* (Tang et al., under review). It contains two things.

- `scaling_report/` — the reporting procedure as an installable package (numpy only), so
  that a reader can compute the same number for their own section.
- `reproduce/` — the analysis pipeline that produced every number, table and figure in the
  paper, together with the derived tables and the numbers lock file from which the
  manuscript text is filled. Start at [`reproduce/README.md`](reproduce/README.md).

The rest of this file documents the package.

## The package

Report the scaling exponent of **binned spatial transcriptomics** in a form that can be
compared between sections — and see, before writing the sentence, how much of it the
number actually identifies.

Input: mean and variance of per-bin counts for each gene, at a ladder of bin sizes.
Output: one JSON, one Markdown report, and a four-sentence reporting rule filled with
your numbers.

Nothing else is required. No spatial framework, no raw matrix, no AnnData.

```
python -m scaling_report examples/colon_moments.npz --out report/
```

```
Visium HD colon, section P0: γ_pool = 0.5547 | per-gene median 0.428 (n=13,115)
| IGKC 88% | n_eff 1.3
```

## The two things people get wrong

**1. The pooled exponent is a composition statistic, not a tissue property.**
Writing α_pool(d) = Σ_g w_g(d) α_g(d) with w_g = μ_g² / Σ μ² makes it a weighted mean of
per-gene exponents, so it inherits whatever transcript happens to dominate the variance.
In the colon section above, 1.3 genes decide the reported number and one transcript
(IGKC) supplies 88% of it. Report that share, or the number cannot be compared with a
tissue whose immunoglobulin fraction is different.

**2. There are two effective-gene counts, and only one of them is about the exponent.**
`w` is the weight of the *variance*; `v_g = w_g α_g / α_pool` is the weight of the
*exponent*. In the colon section they give 2.5 and 1.3. The paper's argument is that
1/Σw² answers the wrong question.

## What the pooled value does not identify

Hold γ_pool fixed and move the heaviest transcript; the per-gene median moves with it:

```
d(median) / d(anchor) = -v1 / (1 - v1)
```

For the colon section v1 = 0.877, so the slope is −7.16: the median can be pushed
anywhere in 0.00–1.00 while the reported pooled value stays at 0.5547. This is an
existence proof, not an optimisation — put the heavy gene at `a` and every other gene at
`c = (γ_pool − v1·a)/(1 − v1)`, and `c` *is* the median because it is the value of the
vast majority.

Consequence: γ_pool can compare sections; the per-gene median cannot. The tool says this
in the report, with your numbers.

## Install

```
git clone <this repo>
cd scaling-report
pip install -e .
scaling-report --help
```

Python ≥ 3.9, numpy ≥ 1.22. Nothing else.

## Input

**A. Native `.npz`** — the layout the analysis scripts write:

```python
np.savez_compressed("moments.npz",
                    ds=np.array([2, 4, 8, 16, 32, 64, 128, 256]),
                    genes=genes,            # (G,) str
                    mu_2=..., var_2=...,    # (G,) float, one pair per bin size
                    mu_4=..., var_4=..., ...)
```

**B. Tidy `.tsv` / `.csv`** — one row per gene per scale:

```
gene    scale   mu      var
IGKC    2       3.104   18.77
IGKC    4       4.915   31.02
```

Column names are matched case-insensitively; `bin`, `bin_size` and `d_um` are accepted
for `scale`, `mean` and `variance` for `mu` and `var`.

Means and variances are of **per-bin counts** at each bin size, computed over all in-tissue
bins. Aggregating 2 µm bins to 16 µm means summing counts inside each 8×8 block and taking
mean and variance across blocks — not averaging the 2 µm moments.

## Options

| flag | meaning |
|---|---|
| `--out DIR` | where the report goes (default `scaling_report_out`) |
| `--name LABEL` | how the section is named in the report |
| `--min-mu X` | gene filter: minimum mean counts at the largest bin (default `0.5`) |
| `--compare OTHER` | report the difference to a second table, and refuse to give a CI |
| `--quiet` | no stdout summary |

Exit codes: `0` ok, `2` unreadable input, `3` the pooled overdispersion is not positive.

## What it deliberately will not do

It will not give you a confidence interval. Bins are spatially correlated, so resampling
bins independently understates the error by roughly the ratio of the analysis scale to the
correlation length — at 512 µm tiles against a ~30–60 µm correlation length that is a
factor of ten or more. A CI needs either a spatial block bootstrap on the raw matrix (with
block size above the correlation length) or **at least two sections per group**, in which
case

```
MDE = (1.96 + 0.84) · √2 · SD_between / √n
```

is the smallest difference you could have detected. `--compare` prints the difference and
this warning instead of a number that would look like one.

It also will not tell you the absolute exponent of a tissue. The value depends on the bin
ladder, the gene filter and the platform; the tool prints all three next to every number so
that a reader can tell what is being compared.

## Reproducing the published example

`examples/colon_moments.npz` is the Visium HD human colon cancer section (10x Genomics,
spatial-exp 3.0.0) re-aggregated from native 2 µm bins to 2–256 µm. The tool returns:

| quantity | value |
|---|---|
| γ_pool | 0.5547 |
| genes fitted (μ₂₅₆ ≥ 0.5) | 13,115 of 18,074 |
| per-gene median (IQR) | 0.428 (0.338–0.585) |
| top transcript | IGKC, 87.7 % of the exponent |
| effective genes, exponent weight | 1.30 |
| effective genes, variance weight | 2.53 |
| d(median)/d(anchor) | −7.16 |
| median admitted at fixed γ_pool | 0.000–1.000 |

`tests/test_core.py` asserts these values, plus exact recovery on synthetic data where the
per-gene exponent is known by construction (`examples/make_synthetic.py`).

```
python tests/test_core.py        # or: pytest tests/
```

## Reproducing the paper

`reproduce/` mirrors the working tree in which the analysis was run, and ships the derived
tables (19 MB) so that every printed number can be re-derived without downloading the 41 GB
of raw data. [`reproduce/README.md`](reproduce/README.md) gives the run order, the mapping
from script to figure, and where the raw data come from.

One reading note: `reproduce/07_SUBMISSION/` holds the **working template** of the manuscript
(v4 wording), kept because that is the file the numbers are injected into. The text was
rewritten before submission, so the prose there differs from the submitted version; the
numbers lock, the derived tables and the analysis scripts are the ones behind it.

The two analyses that were attempted and failed — predicting the whole α(d) curve from the
measured covariogram, by two independent routes — are kept in that directory instead of
being deleted, and are written up in Supplementary Note 1 of the paper. A negative result
that a reader would otherwise try to reproduce is worth more stated than omitted.

## Citation

If you use this, cite the paper (see `CITATION.cff`) and the version DOI printed on the
Zenodo record.

## Licence

MIT.

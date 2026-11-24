# scaling_report

**One command, from a Space Ranger folder to a decision.**

```
python -m scaling_report --from-10x path/to/binned_outputs --out report/
```

Then open `report/scaling_report.html`. It says, in one sentence, what a Poisson
test would cost you at each bin size, and which bin size is the largest you can
use without inflating your false-positive rate past a level you choose.

That is the whole product. Nothing else needs to be installed: numpy and the
Python standard library, no scipy, no h5py, no anndata, no spatial framework.

---

## What it answers

Three questions, in the order people actually ask them.

| Question | Where the answer is |
|---|---|
| *Can I use a Poisson test at 256 µm?* | the headline of the HTML report |
| *What is the largest bin I can use and still trust a p-value?* | the "largest safe bin" card |
| *How much sample size did I lose by binning?* | the variance-ratio column |

The exponent γ is still reported, with the same decomposition as before, because
it is what the decision is computed from. But it is no longer the first thing a
reader sees.

## Three ways in

**A. A Space Ranger run.** Point at the run directory, its `outs/` folder, or its
`binned_outputs/`. Native `square_002um`, `square_008um` and `square_016um`
matrices are read as they are — no aggregation, no extra memory — and the ladder
is extended upward from the coarsest native scale, where blocks are few.

```
python -m scaling_report --from-10x ~/visium_hd/outs --out report/
```

**B. One Matrix Market directory.** For a platform that deposits only the finest
binning, or a Space Ranger run whose `.h5` is all you kept (re-run Space Ranger
with `--create-binned-outputs`, or convert the `.h5` to `.mtx` once).

```
python -m scaling_report --from-mtx path/to/filtered_feature_bc_matrix --base-um 2
```

**C. A moments table you already computed.** The v1.0 interface, unchanged.

```
python -m scaling_report examples/colon_moments.npz --out report/
```

## What comes out

| File | For |
|---|---|
| `scaling_report.html` | **open this.** Two figures, a decision, a checklist. Prints to A4 |
| `scaling_report.md` | the same, as text, for pasting into a notebook |
| `scaling_report.json` | every number, for a pipeline |
| `advice.json` | the per-scale table and the sentences, for scripting |

The HTML is one file with hand-written SVG: no JavaScript, no external font, no
network. It opens from a USB stick on a machine with nothing installed.

## Before you trust it: run the check

```
python -m scaling_report --selftest
```

Seven checks, a few seconds. It builds a synthetic section whose exponent is known
by construction (a power-law intensity field with Poisson counts on top), reads it
back through the same code path as a real run, and confirms the exponent comes out
where it was put. It also confirms the false-positive arithmetic against the
paper's value (a variance ratio of 30 must give 72%).

If this fails, nothing else the tool says is worth reading.

## Try it without any data

```
python -m scaling_report --demo demo_section
python -m scaling_report --from-10x demo_section --out report/ --ds 2,4,8,16,32,64
```

Writes a 512 × 512 grid of 2 µm bins with 200 genes, laid out exactly like a
Space Ranger output, in about 45 seconds.

## Options

| flag | meaning |
|---|---|
| `--ds 2,4,...,256` | the bin ladder (default `2,4,8,16,32,64,128,256`) |
| `--min-mu X` | gene filter: minimum mean count at the largest bin (default `0.5`) |
| `--fpr-target 0.10` | the false-positive rate you will tolerate; sets "largest safe bin" |
| `--out DIR` | where the report goes (default `scaling_report_out`) |
| `--name LABEL` | how the section is named |
| `--no-html` | skip the HTML |
| `--compare OTHER` | difference to a second table, and refuses to invent a CI |
| `--quiet` | no stdout summary |

Exit codes: `0` ok, `2` unreadable input, `3` the pooled overdispersion is not
positive.

## The arithmetic behind the advice

A negative-binomial count with mean μ and overdispersion α has variance
μ + αμ², so its variance-to-mean ratio is

```
φ = 1 + α μ
```

A test that assumes Poisson divides by √μ where it should divide by √(φμ), so its
z statistic is inflated by √φ, and a two-sided test at level a rejects at

```
FPR(φ) = 2 (1 − Φ(z_crit / √φ)),    z_crit = Φ⁻¹(1 − a/2)
```

which is 5% at φ = 1 and 72% at φ = 30.

**Which genes.** φ is reported as the median over the **2,000 most abundant genes
at each scale**, which is the definition in the paper's S2 Fig. It matters: the
median over every fitted gene is about five times smaller at the coarsest bin,
because the fitted set is dominated by low-count genes whose Poisson term still
hides the overdispersion. Using the wrong set makes this tool disagree with the
paper it accompanies — `advice.json` carries both numbers so the gap is visible.

**Verified against the paper.** On colon section P1 the tool returns φ = 1.083 at
8 µm and 29.5 at 256 µm, i.e. 6.0% and 71.8% false positives; the manuscript
prints 1.08, 30, 6.0% and 72%. On the index section P2 the same code returns 78%
at 256 µm, which is a property of that section, not a discrepancy.

## What the exponent does and does not identify

`gamma_pool` is a composition statistic. Writing α_pool(d) = Σ w_g α_g with
w_g = μ_g² / Σμ² makes it a weighted mean of per-gene exponents, so it inherits
whatever transcript dominates the variance. In the colon section 2.5 genes carry
the variance behind the reported number and 1.3 carry the exponent itself, and
one transcript (IGKC) supplies 88% of it. The report prints that share.

Two effective-gene counts are carried, because they answer different questions:
`w` is the weight of the *variance*, `v_g = w_g α_g / α_pool` is the weight of
the *exponent*. The paper reports 2.5; 1.3 is the sharper count for the exponent.

Holding γ_pool fixed, the per-gene median can be moved anywhere in 0.00–1.00 by
relocating the heaviest transcript — an existence proof, not an optimisation. The
pooled value can compare sections; the per-gene median cannot.

## What it deliberately will not do

- **No confidence interval.** Bins are spatially correlated, so resampling bins
  independently understates the error by roughly the ratio of the analysis scale
  to the correlation length. A CI needs a spatial block bootstrap on the raw
  matrix, or at least two sections per group, in which case
  `MDE = (1.96 + 0.84)·√2·SD_between/√n` is the smallest difference you could
  have detected. `--compare` prints the difference and this warning.
- **No absolute exponent for a tissue.** The value depends on the ladder, the gene
  filter and the platform; all three are printed next to it.
- **No `.h5`.** Only Matrix Market text. If that is all you have, convert once.

## Install

```
pip install git+https://github.com/tang-sir38/scaling-report.git
scaling-report --help
```

Or from a clone:

```
git clone https://github.com/tang-sir38/scaling-report.git
cd scaling-report
pip install -e .
```

Python ≥ 3.9, numpy ≥ 1.22. Nothing else.

## Reproducing the published example

`examples/colon_moments.npz` is the Visium HD human colon cancer index section
(10x Genomics, spatial-exp 3.0.0) re-aggregated from native 2 µm bins.

| quantity | value |
|---|---|
| γ_pool | 0.5547 |
| genes fitted (μ₂₅₆ ≥ 0.5) | 13,115 of 18,074 |
| per-gene median (IQR) | 0.428 (0.338–0.585) |
| top transcript | IGKC, 87.7 % of the exponent |
| effective genes, exponent weight | 1.30 |
| false positives at 256 µm | 78 % |
| largest bin under 10 % | 8 µm |

`tests/test_core.py` asserts the v1.0 values; `tests/test_tenx.py` checks the
aggregation against arithmetic small enough to do by hand (a 4 × 4 grid summed to
4 µm must give 14, 22, 46, 54).

```
python tests/test_core.py
python tests/test_tenx.py
```

## Reproducing the paper

`reproduce/` mirrors the working tree in which the analysis was run, and ships the
derived tables so that every printed number can be re-derived without downloading
the 41 GB of raw data. [`reproduce/README.md`](reproduce/README.md) gives the run
order, the mapping from script to figure, and where the raw data come from.

`reproduce/07_SUBMISSION/` holds the working template of the manuscript (v4
wording), kept because that is the file the numbers are injected into; the prose
was rewritten before submission.

The two analyses that were attempted and failed — predicting the whole α(d) curve
from the measured covariogram, by two independent routes — are kept rather than
deleted, and are written up in Supplementary Note 1.

## Citation

Cite the paper (see `CITATION.cff`) and the version DOI on the Zenodo record.

## Licence

MIT.

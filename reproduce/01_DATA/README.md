# Data

Nothing in this directory is redistributed. All five sections are public.

| section | platform | source |
|---|---|---|
| Visium HD Human Colon Cancer | Visium HD, spatial-exp 3.0.0 | 10x Genomics public dataset of that name |
| Visium HD Human Colon Cancer P1 | Visium HD, second section, same block | 10x Genomics |
| Visium HD Human Lung Cancer | Visium HD, spatial-exp 3.0.0 | 10x Genomics |
| Visium HD Mouse Brain | Visium HD, spatial-exp 3.0.0 | 10x Genomics |
| Stereo-seq Mouse Whole Brain | Stereo-seq | STOmicsDB, dataset `C04042E2` |

The 10x datasets are downloaded from <https://www.10xgenomics.com/datasets> by searching the
dataset name; each provides a `binned_outputs/` directory containing
`square_002um/filtered_feature_bc_matrix.h5`. The Stereo-seq section is distributed through
the STOmicsDB portal (<https://db.cngb.org/stomics/>) under `C04042E2`; only its expression
table is needed, and `59_gef_probe.py` / `60_gef_pull.py` fetch it by HTTP range request
rather than by downloading the whole archive.

## Expected layout

Each section gets its own working directory, one level above `03_SCRIPTS`:

```
reproduce/
├── 01_DATA/Visium_HD_Human_Colon_Cancer/binned_outputs/square_002um/...
├── 01_WORK/          written by 02b_rebin.py: bins_002um.npz ... bins_256um.npz
├── 02_OUT/           written by 13_gene_moments.py and everything downstream
└── 03_SCRIPTS/
```

For the colon section the working directory is `reproduce/` itself, matching how the
analysis was run: `01_WORK/bins_002um.npz` is 2.0 GB, and the eight re-aggregations together
are 8.5 GB. The mouse brain and lung sections use the same layout with their own prefix
(`D2_MouseBrain/`, `D3_LungCancer/`); the scripts take the section as a command-line
argument.

## What Tier B costs

41 GB of binned matrices, 8.5 GB per human section. The re-aggregation in `02b_rebin.py` is
exact — every 2 µm bin is assigned to exactly one coarser bin by integer division of its
array coordinates, with no interpolation and no read discarded — so this step changes the
resolution and nothing else. That property is what makes the exponent comparable across the
ladder, and it is why nothing is resampled in this pipeline.

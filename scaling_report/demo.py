"""Synthetic Visium HD directory with a known exponent.

Two uses, both practical.

A reader who has not yet run the tool can run it on this and see that the
pipeline works on their machine before pointing it at 40 GB of real counts.

A reviewer who wants to know whether the aggregation is correct can check it
against a case where the answer is known by construction: the intensity field is
drawn from a power-law spectrum S(k) ~ k^-beta, block averaging then gives

    Var(block mean) ~ L^(beta - 2)   =>   gamma = 2 - beta

so asking for gamma = 0.55 means beta = 1.45, and a pipeline that reports 0.55
back has neither lost the correlation structure nor invented one.

The counts themselves are Poisson draws, so the Poisson term is present in the
data exactly as it is in a real section. It is the block averaging that makes it
vanish, which is the point being checked.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np


def _power_law_field(n: int, beta: float, rng: np.random.Generator) -> np.ndarray:
    """Unit-variance field with spectrum k^-beta on an n x n grid."""
    k = np.fft.fftfreq(n)
    kx, ky = np.meshgrid(k, k, indexing="ij")
    kmag = np.sqrt(kx * kx + ky * ky)
    kmag[0, 0] = np.inf
    amp = kmag ** (-beta / 2.0)
    amp[0, 0] = 0.0
    noise = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    field = np.real(np.fft.ifft2(amp * noise))
    sd = field.std()
    if sd == 0:
        raise RuntimeError("degenerate field; raise the grid size")
    return (field - field.mean()) / sd


def make_demo(out: Path, *, n: int = 512, n_genes: int = 200,
              gamma: float = 0.55, mean_umi: float = 1.0,
              seed: int = 20261123, verbose: bool = True) -> Path:
    """Write a Space Ranger-shaped directory with a known gamma.

    ``n`` is the number of 2 um bins on each side, so 512 means a 1 mm square.
    Genes are given a spread of expression levels and of correlation length, as
    in a real section; the *pooled* exponent is set by the spectrum and is what
    the demo checks against.
    """
    out = Path(out)
    root = out / "binned_outputs" / "square_002um" / "filtered_feature_bc_matrix"
    root.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(seed)
    beta = 2.0 - gamma

    genes = [f"GENE{i:04d}" for i in range(n_genes)]
    # expression: a few dominant transcripts, the rest low — the concentration
    # the paper warns about is present here too
    weights = rng.pareto(1.2, size=n_genes) + 0.05
    weights = weights / weights.sum()
    gene_scale = mean_umi * weights * n_genes * 0.15

    rows_all, cols_all, gene_all, val_all = [], [], [], []
    for gi in range(n_genes):
        # per-gene correlation length varies: high-expression genes decorrelate
        # faster, which is the empirical pattern the paper measures
        b = beta + 0.35 * (gi / max(n_genes - 1, 1))
        z = _power_law_field(n, b, rng)
        lam = gene_scale[gi] * (1.0 + 0.45 * z)
        lam = np.clip(lam, 0.0, None)
        counts = rng.poisson(lam)
        r, c = np.nonzero(counts)
        if r.size == 0:
            continue
        rows_all.append(r.astype(np.int32))
        cols_all.append(c.astype(np.int32))
        gene_all.append(np.full(r.size, gi, dtype=np.int32))
        val_all.append(counts[r, c].astype(np.int64))

    rows = np.concatenate(rows_all)
    cols = np.concatenate(cols_all)
    gene = np.concatenate(gene_all)
    val = np.concatenate(val_all)

    # Space Ranger writes one column per barcode, in row-major order of bins that
    # carry at least one UMI. Reproduce that ordering so the reader is exercised
    # on the same layout it will meet in the wild.
    flat = rows.astype(np.int64) * n + cols
    uniq_bins, bin_inverse = np.unique(flat, return_inverse=True)
    brow = (uniq_bins // n).astype(np.int32)
    bcol = (uniq_bins % n).astype(np.int32)

    order = np.lexsort((gene, bin_inverse))
    gene_o = gene[order] + 1          # Matrix Market is 1-based
    bin_o = bin_inverse[order] + 1
    val_o = val[order]

    with gzip.open(root / "matrix.mtx.gz", "wt", encoding="utf-8", newline="\n") as fh:
        fh.write("%%MatrixMarket matrix coordinate integer general\n")
        fh.write("%metadata_json: {\"software_version\": \"scaling_report demo\"}\n")
        fh.write(f"{n_genes} {uniq_bins.size} {val_o.size}\n")
        block = 1_000_000
        for s in range(0, val_o.size, block):
            e = min(s + block, val_o.size)
            for g, b, v in zip(gene_o[s:e], bin_o[s:e], val_o[s:e]):
                fh.write(f"{g} {b} {v}\n")

    with gzip.open(root / "barcodes.tsv.gz", "wt", encoding="utf-8", newline="\n") as fh:
        for r, c in zip(brow, bcol):
            fh.write(f"s_{r:05d}_{c:05d}-1\n")

    with gzip.open(root / "features.tsv.gz", "wt", encoding="utf-8", newline="\n") as fh:
        for i, g in enumerate(genes):
            fh.write(f"ENSG{i:011d}\t{g}\tGene Expression\n")

    if verbose:
        print(f"demo written to {out}")
        print(f"  grid {n}x{n} 2 um bins, {n_genes} genes, "
              f"{val_o.size:,} non-zero entries")
        print(f"  built with gamma = {gamma:.2f} (spectrum exponent "
              f"beta = {beta:.2f})")
        print(f"  read it with:  scaling-report --from-10x {out}")
    return out

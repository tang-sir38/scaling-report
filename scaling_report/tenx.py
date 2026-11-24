"""Read Space Ranger output directly, without a spatial framework.

The package in v1.0 asked for a table of per-bin means and variances. That table
is the *result* of the hardest step, so a reader with a Visium HD experiment had
to write the aggregation themselves before the tool would start. This module
removes that step.

Two inputs are accepted.

**A. A Space Ranger run directory** (``--from-10x``).  Space Ranger writes
``binned_outputs/square_002um``, ``square_008um`` and ``square_016um``; each
holds a count matrix over bins of that edge length.  Every one of those is a
valid scale of the ladder, so the three scales cost nothing to obtain: read the
matrix, take the mean and variance of each gene across bins, done.  The ladder is
then extended upwards (16 -> 32 -> ... -> 256 um) by summing the 16 um bins into
blocks, which is cheap because there are only ~80,000 of them.

**B. A single Matrix Market directory** (``--from-mtx``).  One matrix, one bin
size; the ladder is built upwards from there.  This is the path for platforms
that deposit only the finest binning, and it is slower because the finest matrix
is the largest one.

Nothing outside numpy and the standard library is imported. No h5py, no scipy,
no anndata: the count matrix is read as Matrix Market text and aggregated with
``bincount``.  ``.h5`` files are not read; if that is all you have, run
``spaceranger`` once more with ``--create-binned-outputs`` and use the ``.mtx``
that comes with it.

Barcodes.  A Visium HD square-bin barcode looks like ``s_00012_00345-1``, where
the two numbers are the zero-based row and column of the bin on the finest grid.
The grid geometry is read from the barcodes themselves, so no
``tissue_positions.parquet`` is needed.  If your barcodes do not follow that
pattern the reader says so instead of guessing.
"""
from __future__ import annotations

import gzip
import io
import json
import re
from pathlib import Path

import numpy as np

from .core import DEFAULT_DS, Moments, MomentsError

# --------------------------------------------------------------------------
# locating files
# --------------------------------------------------------------------------

_BARCODE_RE = re.compile(r"^[a-zA-Z]*_?(\d+)_(\d+)[-_]\d+$")
_MATRIX_NAMES = (
    "filtered_feature_bc_matrix/matrix.mtx.gz",
    "filtered_feature_bc_matrix/matrix.mtx",
    "raw_feature_bc_matrix/matrix.mtx.gz",
    "raw_feature_bc_matrix/matrix.mtx",
    "matrix.mtx.gz",
    "matrix.mtx",
)


def _open_text(path: Path):
    if path.suffix == ".gz":
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _resolve_matrix(root: Path) -> Path:
    """Find the count matrix under a Space Ranger binned-output directory."""
    if root.is_file():
        return root
    for rel in _MATRIX_NAMES:
        p = root / rel
        if p.exists():
            return p
    # one level down, the usual layout: <run>/binned_outputs/square_016um
    hits = sorted(root.glob("**/matrix.mtx*"))
    hits = [h for h in hits if h.suffix in (".mtx", ".gz")]
    if len(hits) == 1:
        return hits[0]
    if len(hits) > 1:
        listed = "\n  ".join(str(h) for h in hits[:8])
        raise MomentsError(
            f"{root} holds {len(hits)} count matrices; point at one of them:\n  {listed}")
    raise MomentsError(f"no matrix.mtx or matrix.mtx.gz under {root}")


def _sibling(matrix: Path, stem: str) -> Path:
    """barcodes.tsv / features.tsv next to a matrix, gzipped or not."""
    for cand in (matrix.parent / f"{stem}.tsv.gz", matrix.parent / f"{stem}.tsv"):
        if cand.exists():
            return cand
    raise MomentsError(f"{stem}.tsv[.gz] not found next to {matrix}")


# --------------------------------------------------------------------------
# reading one matrix
# --------------------------------------------------------------------------


def _read_lines(path: Path) -> list[str]:
    with _open_text(path) as fh:
        return [ln.rstrip("\n") for ln in fh if ln.strip()]


def read_barcodes(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Return (row, col) of every barcode, in file order.

    Raises if the barcodes do not encode a grid position, because every later
    step depends on that geometry being right.
    """
    lines = _read_lines(path)
    if not lines:
        raise MomentsError(f"{path} is empty")
    rows = np.empty(len(lines), dtype=np.int32)
    cols = np.empty(len(lines), dtype=np.int32)
    for i, ln in enumerate(lines):
        m = _BARCODE_RE.match(ln)
        if m is None:
            raise MomentsError(
                f"barcode {ln!r} in {path.name} does not encode a grid position. "
                "Visium HD square bins look like s_00012_00345-1; if your platform "
                "names bins differently, export the moments table instead "
                "(see README, 'Input').")
        rows[i] = int(m.group(1))
        cols[i] = int(m.group(2))
    return rows, cols


def read_features(path: Path) -> np.ndarray:
    names = []
    for ln in _read_lines(path):
        parts = ln.split("\t")
        names.append(parts[1] if len(parts) > 1 else parts[0])
    if not names:
        raise MomentsError(f"{path} is empty")
    return np.array(names, dtype=object)


def read_matrix(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, int, int]:
    """Read a Matrix Market coordinate file.

    Returns (gene_index, bin_index, count, n_genes, n_bins) with zero-based
    indices.  The file is 1-based and gene-major (column of bins per gene row),
    but the storage order does not matter here because every consumer sorts.
    """
    with _open_text(path) as fh:
        header = fh.readline().strip()
        if not header.startswith("%%MatrixMarket"):
            raise MomentsError(f"{path} is not a Matrix Market file")
        if "coordinate" not in header.lower():
            raise MomentsError(f"{path} is not in coordinate format")
        line = ""
        while True:
            line = fh.readline()
            if not line:
                raise MomentsError(f"{path} ends before the size line")
            if not line.startswith("%"):
                break
        n_rows, n_cols, nnz = (int(x) for x in line.split()[:3])
        body = np.fromstring(fh.read(), dtype=np.int64, sep=" ")

    if body.size != 3 * nnz:
        raise MomentsError(
            f"{path}: header says {nnz} entries, found {body.size // 3}")
    body = body.reshape(nnz, 3)
    # Matrix Market is 1-based; Space Ranger writes genes x bins
    gene = body[:, 0].astype(np.int32) - 1
    binc = body[:, 1].astype(np.int32) - 1
    val = body[:, 2].astype(np.float64)
    return gene, binc, val, n_rows, n_cols


# --------------------------------------------------------------------------
# moments of one matrix
# --------------------------------------------------------------------------


def moments_of(gene: np.ndarray, val: np.ndarray, n_genes: int,
               n_bins: int, genes: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Mean and variance of each gene's count across all bins.

    Computed from the first two moments, so the bin x gene matrix is never
    materialised:

        mean_g = sum_b x_bg / B
        var_g  = sum_b x_bg^2 / B - mean_g^2

    ``B`` counts every bin in the matrix, including bins with no UMI for that
    gene, which is what makes the variance a property of the tissue rather than
    of the expressed subset.
    """
    if n_bins <= 0:
        raise MomentsError("a count matrix with no bins")
    s1 = np.bincount(gene, weights=val, minlength=n_genes)
    s2 = np.bincount(gene, weights=val * val, minlength=n_genes)
    mu = s1 / n_bins
    var = s2 / n_bins - mu * mu
    # rounding can push a variance a hair below zero when a gene is constant
    var = np.maximum(var, 0.0)
    return mu, var


# --------------------------------------------------------------------------
# building the ladder
# --------------------------------------------------------------------------


def _block_sum(bin_r: np.ndarray, bin_c: np.ndarray, binc: np.ndarray,
               gene: np.ndarray, val: np.ndarray, n_genes: int, ratio: int
               ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Sum ``ratio x ratio`` blocks of bins.

    Two levels are carried separately and must not be confused:

    * per bin — ``bin_r``/``bin_c`` hold the grid position of each bin, so their
      length is the number of bins;
    * per non-zero entry — ``binc``/``gene``/``val``, whose length is the number
      of stored entries and which index into the bin arrays.

    Returns (new_r, new_c, new_binc, new_gene, new_val): the first two give the
    grid position of each coarse bin; the rest is the aggregated matrix, with the
    same two-level convention.
    """
    if ratio == 1:
        return bin_r, bin_c, binc, gene, val

    r = bin_r.astype(np.int64) // ratio
    c = bin_c.astype(np.int64) // ratio
    width = int(c.max()) + 1
    key = r * width + c

    uniq, bin_map = np.unique(key, return_inverse=True)
    new_r = (uniq // width).astype(np.int32)
    new_c = (uniq % width).astype(np.int32)
    new_binc = bin_map[binc].astype(np.int32)

    k2 = new_binc.astype(np.int64) * n_genes + gene.astype(np.int64)
    u2, inv2 = np.unique(k2, return_inverse=True)
    tot = np.bincount(inv2, weights=val)
    g2 = (u2 % n_genes).astype(np.int32)
    b2 = (u2 // n_genes).astype(np.int32)
    return new_r, new_c, b2, g2, tot


def ladder_from_one_matrix(matrix: Path, base_um: int,
                           ds: tuple[int, ...] = DEFAULT_DS) -> Moments:
    """Build the whole ladder from a single finest-resolution matrix."""
    bc = _sibling(matrix, "barcodes")
    ft = _sibling(matrix, "features")
    rows, cols = read_barcodes(bc)
    genes = read_features(ft)
    gene, binc, val, n_genes, n_bins = read_matrix(matrix)

    if n_bins != rows.size:
        raise MomentsError(
            f"{matrix}: {n_bins} columns but {rows.size} barcodes")
    if n_genes != genes.size:
        raise MomentsError(
            f"{matrix}: {n_genes} matrix rows but {genes.size} feature names")

    ds = tuple(int(d) for d in ds if d >= base_um)
    if not ds:
        raise MomentsError("no requested bin size is >= the matrix resolution")

    mu_cols, var_cols = [], []
    r, c, b, g, v = rows, cols, binc, gene, val
    cur = base_um
    for target in ds:
        if target != cur:
            ratio = target // cur
            if ratio < 1 or target % cur:
                raise MomentsError(
                    f"cannot reach {target} um from {cur} um by integer blocks")
            r, c, b, g, v = _block_sum(r, c, b, g, v, n_genes, ratio)
            cur = target
        mu, var = moments_of(g, v, n_genes, int(r.size), genes)
        mu_cols.append(mu)
        var_cols.append(var)
    return Moments(ds=np.array(ds, float), genes=genes,
                   mu=np.column_stack(mu_cols), var=np.column_stack(var_cols))


def find_binned_outputs(root: Path) -> dict[int, Path]:
    """Map bin edge length -> directory, for a Space Ranger run directory.

    Accepts either the run root (looks under ``binned_outputs/`` and
    ``outs/binned_outputs/``) or the ``binned_outputs`` directory itself.
    """
    if root.is_file():
        raise MomentsError("--from-10x expects a directory, not a file")
    candidates = [root / "binned_outputs", root / "outs" / "binned_outputs", root]
    found: dict[int, Path] = {}
    for base in candidates:
        if not base.is_dir():
            continue
        for d in sorted(base.iterdir()):
            if not d.is_dir():
                continue
            m = re.search(r"square_(\d+)um", d.name)
            if m is None:
                continue
            um = int(m.group(1))
            if any((d / rel).exists() for rel in _MATRIX_NAMES):
                found.setdefault(um, d)
    if not found:
        raise MomentsError(
            f"no square_<N>um binned output with a count matrix under {root}. "
            "Expected binned_outputs/square_002um, square_008um, square_016um.")
    return found


def ladder_from_run(root: Path, ds: tuple[int, ...] = DEFAULT_DS) -> tuple[Moments, dict]:
    """Build the ladder from every binned output a Space Ranger run provides.

    The native scales are free: no aggregation, no extra memory. The ladder is
    then extended upwards from the coarsest native scale, where blocks are few.
    """
    found = find_binned_outputs(root)
    native = sorted(found)
    ds = tuple(int(d) for d in ds)

    mu_cols, var_cols, used = [], [], []
    for um in native:
        if um not in ds:
            continue
        matrix = _resolve_matrix(found[um])
        bc = _sibling(matrix, "barcodes")
        ft = _sibling(matrix, "features")
        rows, cols = read_barcodes(bc)
        genes = read_features(ft)
        gene, binc, val, n_genes, n_bins = read_matrix(matrix)
        if n_bins != rows.size:
            raise MomentsError(f"{matrix}: {n_bins} columns, {rows.size} barcodes")
        mu, var = moments_of(gene, val, n_genes, n_bins, genes)
        mu_cols.append(mu)
        var_cols.append(var)
        used.append(um)

    if not used:
        raise MomentsError(
            f"none of the native bin sizes {native} is on the requested ladder {list(ds)}")

    # extend above the coarsest native scale
    base = max(used)
    upper = [int(d) for d in ds if d > base]
    if upper:
        matrix = _resolve_matrix(found[base])
        rows, cols = read_barcodes(_sibling(matrix, "barcodes"))
        genes = read_features(_sibling(matrix, "features"))
        gene, binc, val, n_genes, n_bins = read_matrix(matrix)
        if n_bins != rows.size:
            raise MomentsError(f"{matrix}: {n_bins} columns, {rows.size} barcodes")
        r, c, b, g, v = rows, cols, binc, gene, val
        cur = base
        for target in upper:
            ratio = target // cur
            if ratio < 1 or target % cur:
                raise MomentsError(f"cannot reach {target} um from {cur} um")
            r, c, b, g, v = _block_sum(r, c, b, g, v, n_genes, ratio)
            cur = target
            mu, var = moments_of(g, v, n_genes, int(r.size), genes)
            mu_cols.append(mu)
            var_cols.append(var)
            used.append(target)

    order = np.argsort(used)
    ds_arr = np.array([used[i] for i in order], float)
    mu = np.column_stack([mu_cols[i] for i in order])
    var = np.column_stack([var_cols[i] for i in order])
    meta = dict(native=native, used=sorted(used), source=str(root))
    return Moments(ds=ds_arr, genes=genes, mu=mu, var=var), meta

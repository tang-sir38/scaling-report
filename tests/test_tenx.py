"""Tests for the Space Ranger reader.

The aggregation is the part a reader cannot check by eye, so it is checked here
against a matrix small enough to sum by hand.
"""
from __future__ import annotations

import gzip
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scaling_report.tenx import (_block_sum, ladder_from_run, moments_of,
                                 read_barcodes, read_matrix)


def _write_mtx(path: Path, genes: int, bins: int, rows, cols, vals) -> None:
    """rows/cols are zero-based (gene, bin); the file is 1-based."""
    order = np.lexsort((np.asarray(cols), np.asarray(rows)))
    r = np.asarray(rows)[order]
    c = np.asarray(cols)[order]
    v = np.asarray(vals)[order]
    with gzip.open(path, "wt", encoding="utf-8", newline="\n") as fh:
        fh.write("%%MatrixMarket matrix coordinate integer general\n")
        fh.write(f"{genes} {bins} {v.size}\n")
        for a, b, x in zip(r, c, v):
            fh.write(f"{a+1} {b+1} {x}\n")


def test_block_sum_matches_hand_arithmetic():
    """A 4x4 grid of 2 um bins summed to 4 um must equal the 2x2 block sums."""
    # 4x4 grid, 2 genes, one entry per (gene, bin) for simplicity
    n = 4
    rows = np.array([0, 1, 2, 3] * 4, dtype=np.int32)      # grid row of each bin
    cols = np.array([0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3],
                    dtype=np.int32)                          # grid col
    # one non-zero entry per bin: gene 0 everywhere, count = row*4+col+1
    binc = np.arange(16, dtype=np.int32)
    gene = np.zeros(16, dtype=np.int32)
    val = (rows * 4 + cols + 1).astype(np.float64)

    nr, nc, nb, ng, nv = _block_sum(rows, cols, binc, gene, val, 2, ratio=2)

    # coarse grid is 2x2 -> four bins, each holding a 2x2 block
    assert nr.size == 4, nr
    assert nc.size == 4, nc
    # block (0,0) holds rows 0-1, cols 0-1 -> counts 1+2+5+6 = 14
    lookup = {(int(a), int(b)): 0.0 for a, b in zip(nr, nc)}
    for a, b, x in zip(nb, ng, nv):
        lookup[(int(nr[a]), int(nc[a]))] += float(x)
    assert lookup[(0, 0)] == 14.0, lookup
    assert lookup[(0, 1)] == 22.0, lookup   # 3+4+7+8
    assert lookup[(1, 0)] == 46.0, lookup   # 9+10+13+14
    assert lookup[(1, 1)] == 54.0, lookup   # 11+12+15+16


def test_block_sum_merges_the_same_gene():
    """Two entries of one gene in one block become one entry, not two."""
    rows = np.array([0, 0], dtype=np.int32)
    cols = np.array([0, 1], dtype=np.int32)
    binc = np.array([0, 1], dtype=np.int32)
    gene = np.array([3, 3], dtype=np.int32)
    val = np.array([5.0, 7.0])
    nr, nc, nb, ng, nv = _block_sum(rows, cols, binc, gene, val, 10, ratio=2)
    assert nb.size == 1 and ng[0] == 3 and nv[0] == 12.0, (nb, ng, nv)


def test_moments_match_direct_computation():
    """moments_of must equal numpy's mean and var over the dense column."""
    rng = np.random.default_rng(0)
    n_genes, n_bins = 3, 40
    gene, binc, val = [], [], []
    dense = np.zeros((n_bins, n_genes))
    for b in range(n_bins):
        for g in range(n_genes):
            if rng.random() < 0.6:
                x = int(rng.integers(1, 6))
                dense[b, g] = x
                gene.append(g)
                binc.append(b)
                val.append(x)
    mu, var = moments_of(np.array(gene), np.array(val, float), n_genes, n_bins,
                         np.array(["a", "b", "c"]))
    assert np.allclose(mu, dense.mean(axis=0))
    # variance on a finite sample (ddof=0), as the paper defines it
    assert np.allclose(var, dense.var(axis=0))


def test_barcode_round_trip():
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "barcodes.tsv.gz"
        with gzip.open(p, "wt", encoding="utf-8", newline="\n") as fh:
            fh.write("s_00000_00000-1\ns_00000_00003-1\ns_00007_00002-1\n")
        r, c = read_barcodes(p)
        assert list(r) == [0, 0, 7]
        assert list(c) == [0, 3, 2]


def test_reader_rejects_a_barcode_it_cannot_place():
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "barcodes.tsv.gz"
        with gzip.open(p, "wt", encoding="utf-8", newline="\n") as fh:
            fh.write("AAACCTGAGCACCCAC-1\n")
        try:
            read_barcodes(p)
        except Exception as e:
            assert "grid position" in str(e)
        else:
            raise AssertionError("a non-spatial barcode must not be accepted")


def test_ladder_end_to_end_on_a_tiny_run():
    """A full 2 -> 8 um ladder on a run directory laid out like Space Ranger."""
    with tempfile.TemporaryDirectory() as td:
        root = Path(td) / "binned_outputs" / "square_002um" / "filtered_feature_bc_matrix"
        root.mkdir(parents=True)
        n = 4
        rows, cols, gene, val = [], [], [], []
        binc = 0
        for r in range(n):
            for c in range(n):
                rows.append(r)
                cols.append(c)
                gene.append(0)
                val.append(1 + r * n + c)
                binc += 1
        _write_mtx(root / "matrix.mtx.gz", 1, 16,
                   np.array(gene), np.arange(16), np.array(val))
        with gzip.open(root / "barcodes.tsv.gz", "wt", encoding="utf-8",
                       newline="\n") as fh:
            for r, c in zip(rows, cols):
                fh.write(f"s_{r:05d}_{c:05d}-1\n")
        with gzip.open(root / "features.tsv.gz", "wt", encoding="utf-8",
                       newline="\n") as fh:
            fh.write("ENSG00000000001\tG1\tGene Expression\n")

        m, meta = ladder_from_run(Path(td), ds=(2, 4, 8))
        assert list(m.ds) == [2.0, 4.0, 8.0], m.ds
        assert m.n_genes == 1
        # at 8 um the single bin holds 1+..+16 = 136
        assert abs(m.mu[0, -1] - 136.0) < 1e-9, m.mu
        # one bin -> variance across bins is zero by definition
        assert abs(m.var[0, -1]) < 1e-9, m.var


if __name__ == "__main__":
    fails = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"PASS  {name}")
            except Exception as e:
                fails += 1
                print(f"FAIL  {name}: {type(e).__name__}: {e}")
    print(f"\n{fails} failed")
    raise SystemExit(1 if fails else 0)

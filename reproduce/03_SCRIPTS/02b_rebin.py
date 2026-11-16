r"""
02b_rebin.py — Visium HD 2um bin 精确重聚合（组织感知版）。

与 02_rebin.py 的区别（这是关键）:
  02_rebin 只按 barcode 分块求和，分母（格子数）留给下游用"外接矩形"推。
  真实切片不是矩形，组织只占外接矩形的 ~78%。若把组织外的空白格当格子，
  粗分辨率下空格子更多，方差被系统性抬高，alpha 虚高，斜率会被伪造成偏离 -2。
  本脚本改用 tissue_positions.parquet 的 in_tissue 标记定义分析域：
    n_domain_any  = 至少含 1 个 in_tissue 2um 格的粗格数（主口径）
    n_domain_full = 全部 16 个 2um 格都在组织内的粗格数（严格口径，做敏感性对照）
  两个口径都存，下游先跑主口径，再用严格口径复核斜率是否稳定。

重聚合本身仍是精确整数分块求和：d um 的一个格 = (d/2)^2 个 2um 格原地相加，
不插值、不降采样、不丢 reads。

内存: 逐基因 bincount，峰值约 O(nnz_2um) ≈ 3 GB，13.9 GB 内存可跑。

用法:
  python 02b_rebin.py --selftest
  python 02b_rebin.py --h5 <h5> --parquet <parquet> --outdir <dir> --d 2 4 8 16 32 64
  python 02b_rebin.py ... --genes-limit 500        # 只跑前 500 基因，用于测速
"""
import argparse
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, csc_matrix

ROOT = Path(__file__).resolve().parents[1]

KEY_STRIDE = 4096          # 列索引上限 3350，留余量
ROW_SL, COL_SL = slice(8, 13), slice(14, 19)


def parse_barcodes_fast(bb):
    """barcode 定宽切片，向量化。格式 s_002um_RRRRR_CCCCC-1（共 21 字符）。"""
    b = np.asarray(bb)
    if b.dtype.kind == "O":
        b = b.astype("S21")
    if b.dtype.kind == "U":
        b = b.astype("S")
    L = b.dtype.itemsize
    if L != 21:
        raise ValueError(f"barcode 定宽异常：itemsize={L}，示例 {b[0]!r}")
    a = np.ascontiguousarray(b).view("S1").reshape(-1, L)
    rows = a[:, ROW_SL].copy().view("S5").ravel().astype(np.int32)
    cols = a[:, COL_SL].copy().view("S5").ravel().astype(np.int32)
    return rows, cols


def domain_keys(prq_path, factor):
    """组织内 2um 格 -> 粗格 key（去重）。返回 (any_keys, full_keys, info)。"""
    df = pd.read_parquet(prq_path, columns=["in_tissue", "array_row", "array_col"])
    r = df["array_row"].to_numpy(np.int32)
    c = df["array_col"].to_numpy(np.int32)
    t = df["in_tissue"].to_numpy(np.int8) == 1
    r, c = r[t], c[t]

    r0 = (r // factor).astype(np.int32)
    c0 = (c // factor).astype(np.int32)
    key = (r0.astype(np.int32) * np.int32(KEY_STRIDE) + c0).astype(np.int32)
    any_keys = np.unique(key)

    if factor == 1:
        full_keys = any_keys
    else:
        # 一个粗格要"完整在组织内"：它包含的 factor^2 个 2um 格全在组织内
        cnt = pd.Series(key).value_counts()
        full_keys = np.sort(cnt.index.to_numpy()[cnt.to_numpy() == factor * factor])
    return any_keys, full_keys, dict(n_tissue_2um=int(t.sum()))


def rebin_one(X, indptr_i32, indices_i32, rows_h5, cols_h5, factor, verbose=True):
    """逐基因把 2um 计数聚到 factor 网格。返回 (csr, coords_obs, stats)。"""
    rows_h5 = rows_h5
    cols_h5 = cols_h5
    n_genes = X.shape[0]

    if factor == 1:
        coords = np.column_stack([rows_h5, cols_h5]).astype(np.int32)
        tot = np.asarray(X.sum(axis=0)).ravel()
        coords = coords[tot > 0]
        keep = tot > 0
        Xk = X[:, keep].tocsr()
        return Xk, coords, dict(nnz=int(Xk.nnz))

    key_h5 = ((rows_h5 // factor).astype(np.int32) * np.int32(KEY_STRIDE)
              + (cols_h5 // factor).astype(np.int32)).astype(np.int32)

    out_idx, out_dat, out_ptr = [], [], np.zeros(n_genes + 1, dtype=np.int64)
    pos = 0
    data_all = X.data if hasattr(X, "data") else None
    for g in range(n_genes):
        a, b = int(indptr_i32[g]), int(indptr_i32[g + 1])
        if b > a:
            ck = key_h5[indices_i32[a:b]]
            v = data_all[a:b].astype(np.float32)
            uq, inv = np.unique(ck, return_inverse=True)
            sm = np.bincount(inv, weights=v, minlength=uq.size).astype(np.float32)
            pos += uq.size
            out_ptr[g + 1] = pos
            out_idx.append(uq)
            out_dat.append(sm)
        else:
            out_ptr[g + 1] = pos
        if verbose and (g + 1) % 4000 == 0:
            print(f"      {g+1}/{n_genes} 基因, 累计 nnz={pos:,}", flush=True)

    idx_all = np.concatenate(out_idx) if out_idx else np.zeros(0, np.int32)
    dat_all = np.concatenate(out_dat) if out_dat else np.zeros(0, np.float32)

    # 全局 key -> 紧凑列号
    ukeys, inv2 = np.unique(idx_all, return_inverse=True)
    Xc = csr_matrix((dat_all, inv2.astype(np.int32), out_ptr), shape=(n_genes, ukeys.size))
    rr = (ukeys // KEY_STRIDE).astype(np.int32)
    cc = (ukeys % KEY_STRIDE).astype(np.int32)
    coords = np.column_stack([rr, cc])
    return Xc, coords, dict(nnz=int(Xc.nnz))


def selftest():
    """合成带空洞的组织掩膜，核对重聚合与暴力求和逐元素相等。"""
    rng = np.random.default_rng(11)
    G, R, C = 7, 30, 26
    rows = np.repeat(np.arange(R), C).astype(np.int32)
    cols = np.tile(np.arange(C), R).astype(np.int32)
    keep = ~((rows > 20) & (cols < 6))            # 挖掉一块，模拟组织外
    keep &= ~((rows % 7 == 0) & (cols % 5 == 0))  # 撒一些孔洞
    rows, cols = rows[keep], cols[keep]

    n = rows.size
    dense = np.zeros((G, n), dtype=np.float32)
    for j in range(n):
        k = rng.integers(0, 4)
        if k:
            gs = rng.choice(G, size=k, replace=False)
            dense[gs, j] = rng.integers(1, 9, size=k)
    X = csr_matrix(dense)

    ok = True
    for factor in (2, 3, 4):
        keys = (rows // factor).astype(np.int64) * KEY_STRIDE + (cols // factor)
        uq, inv = np.unique(keys, return_inverse=True)
        brute = np.zeros((G, uq.size), dtype=np.float64)
        for gi, k in enumerate(uq):
            brute[:, gi] = dense[:, keys == k].sum(axis=1)
        # 全零的粗格不进入稀疏矩阵（其为零由 n_domain 补），对照时同样剔除
        nz_col = brute.sum(axis=0) > 0
        brute, uq = brute[:, nz_col], uq[nz_col]

        # 与 rebin_one 同逻辑
        out_idx, out_dat, out_ptr = [], [], np.zeros(G + 1, np.int64)
        pos = 0
        for g in range(G):
            idx = np.nonzero(dense[g])[0]
            if idx.size:
                uq2, inv2 = np.unique(keys[idx], return_inverse=True)
                sm = np.bincount(inv2, weights=dense[g, idx], minlength=uq2.size)
                pos += uq2.size
                out_ptr[g + 1] = pos
                out_idx.append(uq2)
                out_dat.append(sm)
            else:
                out_ptr[g + 1] = pos
        ia = np.concatenate(out_idx)
        dd = np.concatenate(out_dat)
        u2, inv3 = np.unique(ia, return_inverse=True)
        got = csr_matrix((dd, inv3, out_ptr), shape=(G, u2.size)).toarray()

        # 对齐列顺序再比
        order_got = np.argsort(u2)
        order_ref = np.argsort(uq)
        same = np.array_equal(got[:, order_got], brute[:, order_ref])
        cons = abs(dense.sum() - got.sum()) < 1e-4
        print(f"  factor={factor}: 形状 {got.shape} vs {brute.shape}  逐元素相等={same}  总数守恒={cons}")
        ok = ok and same and cons
    print("SELFTEST:", "PASS" if ok else "FAIL")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--h5")
    ap.add_argument("--parquet")
    ap.add_argument("--outdir", default=str(ROOT / "01_WORK"))
    ap.add_argument("--d", nargs="+", type=int, default=[2, 4, 8, 16, 32, 64])
    ap.add_argument("--genes-limit", type=int, default=0)
    a = ap.parse_args()

    if a.selftest:
        raise SystemExit(0 if selftest() else 1)

    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with h5py.File(a.h5, "r") as f:
        m = f["matrix"]
        shape = tuple(int(x) for x in m["shape"][:])
        data = m["data"][:]
        indices = m["indices"][:].astype(np.int32)
        indptr = m["indptr"][:].astype(np.int64)
        bb = m["barcodes"][:]
        genes = [n.decode() if isinstance(n, bytes) else str(n) for n in m["features/name"][:]]
        in_tissue_bc = m["barcodes"][:]

    if a.genes_limit:
        n_keep = min(a.genes_limit, shape[0])
    else:
        n_keep = shape[0]

    rows, cols = parse_barcodes_fast(bb)
    # 10x h5 是 CSC：shape=(genes, barcodes)，indptr 长度 = barcodes+1，
    # indices 存的是基因号。必须先按 CSC 构建再转 CSR，否则整张矩阵会转置。
    X = csc_matrix((data.astype(np.float32), indices, indptr), shape=shape).tocsr()
    del data, indices, indptr
    assert X.shape[0] == shape[0] and X.indptr.size == shape[0] + 1
    if n_keep < shape[0]:
        X = X[:n_keep].tocsr()
        genes = genes[:n_keep]
    print(f"读入: {X.shape[0]} genes x {X.shape[1]} bins, nnz={X.nnz:,}", flush=True)
    print(f"网格: row {rows.min()}-{rows.max()}, col {cols.min()}-{cols.max()}  "
          f"= {(rows.max()+1)*2} x {(cols.max()+1)*2} um", flush=True)
    print(f"总 UMI = {X.sum():,.0f}；每 2um 格均值 = {X.sum()/shape[1]:.1f}", flush=True)

    for d in a.d:
        f_ = d // 2
        anyk, fullk, info = domain_keys(a.parquet, f_)
        print(f"\n== d={d}um (factor={f_}) ==", flush=True)
        print(f"  组织内 2um 格 {info['n_tissue_2um']:,}；"
              f"n_domain_any={anyk.size:,}；n_domain_full={fullk.size:,}", flush=True)
        Xc, coords, st = rebin_one(X, np.asarray(X.indptr), np.asarray(X.indices),
                                   rows, cols, f_)
        tag = f"{d:03d}um"
        np.savez(
            outdir / f"bins_{tag}.npz",
            data=Xc.data, indices=Xc.indices, indptr=Xc.indptr,
            shape=np.array(Xc.shape), coords=coords, genes=np.array(genes),
            n_domain_any=np.int64(anyk.size), n_domain_full=np.int64(fullk.size),
            dom_any=anyk, dom_full=fullk,
            total_umi_all=np.float64(Xc.sum()),
        )
        tot = np.asarray(Xc.sum(axis=0)).ravel()
        print(f"  写出 bins_{tag}.npz  观测 bin={Xc.shape[1]:,}  nnz={st['nnz']:,}  "
              f"UMI={tot.sum():,.0f}  每 bin 均值={tot.mean():.1f}", flush=True)


if __name__ == "__main__":
    main()

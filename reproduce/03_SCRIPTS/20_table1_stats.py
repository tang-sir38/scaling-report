r"""
20_table1_stats.py — Table 1 的基本统计，全部从 2 µm 稀疏矩阵现场算。

输出 <outdir>/table1_stats.json，字段：
  n_genes_raw      原始基因数（矩阵行数）
  n_genes_noMT     剔除线粒体（MT- / mt- 前缀）后的基因数
  n_bins_domain    组织内 2 µm 格数（domain="any" 口径，与 07/09 号脚本一致）
  zero_frac        组织内零计数格占比（空 bin 计 0，不计缺失）
  umi_total        组织内总 UMI
  umi_per_bin      每 2 µm 格平均 UMI
  umi_median       每 2 µm 格 UMI 中位数

用法：
  python 20_table1_stats.py --npz 01_WORK/bins_002um.npz --out 02_OUT/table1_stats.json
"""
import argparse
import json
from pathlib import Path

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    z = np.load(args.npz, allow_pickle=False)
    n_genes, n_bins = int(z["shape"][0]), int(z["shape"][1])
    genes = np.asarray(z["genes"])
    mt = np.array([str(g).upper().startswith("MT-") for g in genes])

    n_domain = int(z["n_domain_any"])
    dom = np.asarray(z["dom_any"]).astype(bool) if "dom_any" in z.files else np.ones(n_bins, bool)
    # 矩阵列 = 组织内「至少 1 个计数」的格；零计数格不进稀疏矩阵，但在分母里计 0。
    cols_are_nonempty = dom.size != n_bins

    # 分块累加列和，避免一次性构造 8.7e6 x nnz 的大数组
    colsum = np.zeros(n_bins, dtype=np.float64)
    indptr = z["indptr"]
    data = z["data"]
    indices = z["indices"]
    step = 2_000_000  # nnz 分块
    for g0 in range(0, n_genes, 200):
        g1 = min(g0 + 200, n_genes)
        a, b = int(indptr[g0]), int(indptr[g1])
        if b <= a:
            continue
        idx = indices[a:b]
        val = data[a:b]
        for s in range(0, idx.size, step):
            e = min(s + step, idx.size)
            colsum += np.bincount(idx[s:e], weights=val[s:e], minlength=n_bins)

    cs = colsum if cols_are_nonempty else colsum[dom]
    n_zero = int((cs == 0).sum())
    if cols_are_nonempty and n_zero:
        raise SystemExit("假定列即非空格，但出现零列，需改映射")
    zero_count = n_zero + (n_domain - n_bins if cols_are_nonempty else 0)
    out = {
        "npz": str(Path(args.npz).name),
        "n_genes_raw": n_genes,
        "n_genes_noMT": int((~mt).sum()),
        "n_bins_domain": n_domain,
        "n_bins_matrix": n_bins,
        "n_zero_bins": zero_count,
        "zero_frac": float(zero_count / n_domain),
        "umi_total": float(cs.sum()),
        "umi_per_bin": float(cs.sum() / n_domain),
        "umi_median": float(np.median(cs)),
    }
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

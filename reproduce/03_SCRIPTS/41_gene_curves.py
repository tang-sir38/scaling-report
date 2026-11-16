r"""
41_gene_curves.py — 把逐基因曲线落盘，供出图与复核。

31_gene_decomp.py 只存了汇总（中位、四分位），图需要点本身。
本脚本从已落盘的逐基因矩 (mu_d, var_d) 直接复算 alpha_g(d) 与 gamma_g，
筛选口径与 31 完全一致（8 档 alpha 全正、mu_256 >= 0.5），
并把汇总与 gene_decomp.json 对齐打印，作为口径自检。

输出 <outdir>/gene_curves_<tag>.npz
"""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

SETS = {
    "colon": str(ROOT / "02_OUT" / "gene_moments.npz"),
    "lung": str(ROOT / "D3_LungCancer" / "02_OUT" / "gene_moments_noMT.npz"),
    "brain": str(ROOT / "D2_MouseBrain" / "02_OUT" / "gene_moments_noMT.npz"),
}
OUT = Path(str(ROOT / "02_OUT" / "gene_curves"))
LOCK = json.load(open(str(ROOT / "02_OUT" / "gene_decomp" / "gene_decomp.json"), encoding="utf-8"))
MIN_MU256 = 0.5


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for tag, p in SETS.items():
        z = np.load(p, allow_pickle=True)
        ds = np.array([int(x) for x in z["ds"]], float)
        genes = np.array([str(g) for g in z["genes"]])
        MU = np.column_stack([z[f"mu_{int(d)}"] for d in ds])
        VA = np.column_stack([z[f"var_{int(d)}"] for d in ds])
        AL = (VA - MU) / np.where(MU > 0, MU ** 2, np.nan)

        ok = np.all(np.isfinite(AL) & (AL > 0), axis=1) & (MU[:, -1] >= MIN_MU256)
        idx = np.where(ok)[0]
        y = np.log(AL[idx])
        xd = np.log(ds)
        xm = xd - xd.mean()
        b = (y * xm).sum(1) / (xm ** 2).sum()
        gamma = np.full(genes.size, np.nan)
        gamma[idx] = -b

        g_ok = gamma[idx]
        rec = dict(n_fit=int(idx.size), n_genes=int(genes.size),
                   median=float(np.median(g_ok)),
                   q25=float(np.percentile(g_ok, 25)), q75=float(np.percentile(g_ok, 75)),
                   p95=float(np.percentile(g_ok, 95)), gmax=float(g_ok.max()))
        ref = LOCK[tag]["gamma_gene"]
        drift = {k: rec[k] - ref[k] for k in ("median", "q25", "q75", "p95", "gmax")}
        print(f"[{tag}] n_fit={rec['n_fit']} (lock {LOCK[tag]['n_fit']})  "
              f"med={rec['median']:.4f} (lock {ref['median']:.4f})  "
              f"max drift {max(abs(v) for v in drift.values()):.2e}")
        if rec["n_fit"] != LOCK[tag]["n_fit"]:
            raise SystemExit(f"{tag} 筛选口径与 lock 不一致，停止")

        np.savez_compressed(
            OUT / f"gene_curves_{tag}.npz",
            gene=genes.astype("U16"), ds=ds,
            mu32=MU[:, int(np.where(ds == 32)[0][0])].astype(np.float32),
            mu256=MU[:, -1].astype(np.float32),
            gamma=gamma.astype(np.float32),
            fitted=ok,
        )
        print(f"  写出 gene_curves_{tag}.npz")


if __name__ == "__main__":
    main()

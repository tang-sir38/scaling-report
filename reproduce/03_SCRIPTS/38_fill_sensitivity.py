r"""
38_fill_sensitivity.py — 填充率敏感性：边界半满粗格会不会把 alpha 抬高？

问题
----
主口径 n_domain_any 把"至少含 1 个 in_tissue 2um 格"的粗格都算作分析单元。
在 256 um，一个粗格含 128^2=16384 个 2um 格；切片边界上的粗格只有一部分在组织内，
计数系统性偏低，于是粗格之间的方差里混进了"几何填充率"这一项，不是生物学。
审稿人一定会问："16,384 倍是名义面积，你真的有那么多格？"

做法
----
严格口径 n_domain_full：只保留"全部 2um 格都在组织内"的粗格，分母用其个数。
两种口径都用同一批 in-tissue 基因、同一套矩估计，只换列与分母。
若两口径的 gamma 之差 < 0.02，则主结论不依赖填充率假设。

同时输出每尺度的
  n_domain_any / n_domain_full / 填充率 / 观测非零列数 / 名义面积倍率 vs 实际格数倍率

用法
----
python 38_fill_sensitivity.py --npzdir <dir> --tag colon --outdir <dir>
"""
import argparse
import gc
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

ROOT = Path(__file__).resolve().parents[1]

STRIDE = 4096
DS = [2, 4, 8, 16, 32, 64, 128, 256]


def load(p):
    z = np.load(p, allow_pickle=False)
    shape = tuple(int(x) for x in z["shape"])
    X = csr_matrix((z["data"].astype(np.float32), z["indices"], z["indptr"]), shape=shape)
    return X, z


def moments_sum(X, vec, chunk=1500):
    """per-gene sum 与 sum-of-squares，只对 vec==1 的列求和（不复制列）。"""
    g = X.shape[0]
    s1 = np.empty(g, np.float64)
    s2 = np.empty(g, np.float64)
    for a in range(0, g, chunk):
        b = min(a + chunk, g)
        Xc = X[a:b]
        s1[a:b] = np.asarray(Xc @ vec).ravel()
        s2[a:b] = np.asarray(Xc.multiply(Xc) @ vec).ravel()
        del Xc
    return s1, s2


def pooled_alpha(s1, s2, n):
    mu = s1 / n
    var = (s2 - s1 ** 2 / n) / (n - 1)
    num = float(np.sum(var - mu))
    den = float(np.sum(mu ** 2))
    return num / den, mu, var


def fit_gamma(al):
    ds = np.array([d for d in DS if al.get(d) and al[d] > 0], float)
    aa = np.array([al[d] for d in DS if al.get(d) and al[d] > 0], float)
    c = np.polyfit(np.log(ds), np.log(aa), 1)
    return float(-c[0]), float(1 - np.var(np.log(aa) - np.polyval(c, np.log(ds))) / np.var(np.log(aa)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npzdir", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--outdir", default=str(ROOT / "02_OUT" / "fill_sensitivity"))
    a = ap.parse_args()
    od = Path(a.outdir)
    od.mkdir(parents=True, exist_ok=True)

    out = dict(tag=a.tag, ds=DS, by_scale={})
    al_any, al_full = {}, {}

    for d in DS:
        p = Path(a.npzdir) / f"bins_{d:03d}um.npz"
        X, z = load(p)
        # 与主分析同口径：剔除线粒体基因
        genes = np.array([str(s) for s in z["genes"]])
        keep = ~np.array([s.upper().startswith("MT-") for s in genes])
        if (~keep).sum():
            X = X[keep].tocsr()
        # coords 已经是 d 尺度网格上的 (row, col)（02b 里 key = (row2//f)*STRIDE + col2//f，
        # 再拆回 r//f 与 c//f），不可再除一次 f。
        coords = z["coords"].astype(np.int64)
        colkey = coords[:, 0] * STRIDE + coords[:, 1]

        n_any = int(z["n_domain_any"])
        n_full = int(z["n_domain_full"])
        full_keys = z["dom_full"] if "dom_full" in z.files else z["dom_any"]
        if d == 2:
            # 2um 尺度：dom_full 与 dom_any 相同（无"半格"概念）
            in_full = np.ones(colkey.size, bool)
            n_full = int(np.unique(np.concatenate([colkey, full_keys])).size)
        else:
            in_full = np.isin(colkey, full_keys)

        v_any = np.ones(colkey.size, np.float32)
        v_full = in_full.astype(np.float32)

        s1a, s2a = moments_sum(X, v_any)
        aa, mua, vara = pooled_alpha(s1a, s2a, n_any)
        s1f, s2f = moments_sum(X, v_full)
        af, muf, varf = pooled_alpha(s1f, s2f, n_full)

        al_any[d] = aa
        al_full[d] = af
        out["by_scale"][str(d)] = dict(
            n_domain_any=n_any, n_domain_full=n_full,
            fill_frac=n_full / n_any,
            n_cols_observed=int(colkey.size),
            n_cols_full=int(in_full.sum()),
            alpha_any=float(aa), alpha_full=float(af),
            ratio_full_over_any=float(af / aa),
            umi=s1a.sum(),
        )
        print(f"  [{a.tag}] d={d:3d}  n_any={n_any:,}  n_full={n_full:,}  "
              f"fill={n_full/n_any:.3f}  alpha_any={aa:.4f}  alpha_full={af:.4f}  "
              f"ratio={af/aa:.4f}", flush=True)
        del X, z
        gc.collect()

    g_any, r2_any = fit_gamma(al_any)
    g_full, r2_full = fit_gamma(al_full)
    out["gamma_any"] = g_any
    out["r2_any"] = r2_any
    out["gamma_full"] = g_full
    out["r2_full"] = r2_full
    out["gamma_delta"] = g_full - g_any
    out["alpha_any"] = {str(k): v for k, v in al_any.items()}
    out["alpha_full"] = {str(k): v for k, v in al_full.items()}
    out["nominal_area_ratio"] = 16384.0
    out["actual_bin_ratio"] = out["by_scale"]["2"]["n_domain_any"] / out["by_scale"]["256"]["n_domain_any"]
    out["actual_bin_ratio_full"] = out["by_scale"]["2"]["n_domain_full"] / out["by_scale"]["256"]["n_domain_full"]

    json.dump(out, open(od / f"fill_{a.tag}.json", "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print(f"\n[{a.tag}] gamma_any={g_any:.4f}  gamma_full={g_full:.4f}  "
          f"delta={g_full-g_any:+.4f}")
    print(f"  名义面积倍率 {16384}  实际格数倍率 {out['actual_bin_ratio']:.0f}  "
          f"严格口径 {out['actual_bin_ratio_full']:.0f}")
    print(f"写出 {od/'fill_'+a.tag+'.json'}")


if __name__ == "__main__":
    main()

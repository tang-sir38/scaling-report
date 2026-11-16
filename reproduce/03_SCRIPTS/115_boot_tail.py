#!/usr/bin/env python3
r"""
115_boot_tail.py — 尾段分辨率的 bin 级 bootstrap，用来判断下降倍数是否稳。

为什么只做尾段
--------------
新组织（小鼠小肠）算出的下降倍数是 15.59，超出现有上界 14.31 约 9%。这个数字
由两个端点决定：alpha(2um) 与 alpha(256um)。其中 256um 只有 476 个有观测的格子，
是全链条里样本最少的一档，也是 128->256 段斜率最陡（-0.97，其它组织约 -0.6）的
来源。所以在动标题之前，必须先把 256um 那一档的不确定性量出来。

2um 有 547 万个格子、1.34 亿 UMI，抽样误差相对可忽略；而且对它做 bootstrap 需要
在 1.28 亿非零元的矩阵上反复做 X.multiply(X)，单次就要 2 GB 上下，本机跑不动。
所以这里只对 64/128/256 三档做重采样，把"终点那一档到底有多不确定"钉死。

方法
----
对 bin 做有放回重采样（权重 k_j ~ Multinomial(n, 1/n)），则

    s1* = X  @ k          per-gene 总计数
    s2* = X2 @ k          per-gene 平方和
    n*  = sum(k)          口径内格子数（含零格，与主流程口径一致）
    mu* = s1*/n*,  var* = (s2* - s1*^2/n*)/(n*-1)
    alpha* = SUM(var* - mu*) / SUM(mu*^2)

和 07_pooled_alpha.py 的 pooled_alpha 完全同一公式，只是把求和换成加权求和，
所以 bootstrap 分布与主流程的估计量是同一个东西。

用法
----
  python 115_boot_tail.py --tissue mouse_small_intestine --ds 64 128 256 --n 200
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

ROOT_DEFAULT = Path(r"E:\A_spatial_20261116")


def load_bins(p: Path):
    z = np.load(p, allow_pickle=False)
    shape = tuple(int(x) for x in z["shape"])
    X = csr_matrix((z["data"], z["indices"], z["indptr"]), shape=shape)
    n_cells = int(z["n_domain_any"]) if "n_domain_any" in z.files else int(X.shape[1])
    genes = z["genes"]
    return X, n_cells, genes


def alpha_from_moments(s1, s2, n):
    """与 07_pooled_alpha.py 同式：alpha = SUM(var-mu)/SUM(mu^2)。"""
    mu = s1 / n
    var = (s2 - s1 ** 2 / n) / (n - 1)
    num = float(np.sum(var - mu))
    den = float(np.sum(mu ** 2))
    return (num / den if den > 0 else np.nan), num, den


def boot_one(X, X2, rng):
    """一次 bootstrap：返回 alpha*。X 已按 --exclude-mt 口径过滤。"""
    ncol = X.shape[1]
    k = rng.multinomial(ncol, np.full(ncol, 1.0 / ncol)).astype(np.float64)
    s1 = np.asarray(X @ k).ravel()
    s2 = np.asarray(X2 @ k).ravel()
    return alpha_from_moments(s1, s2, float(k.sum()))[0]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue", required=True)
    ap.add_argument("--root", default=str(ROOT_DEFAULT))
    ap.add_argument("--ds", nargs="+", type=int, default=[64, 128, 256])
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--seed", type=int, default=20261116)
    a = ap.parse_args()

    wd = Path(a.root) / "WORK" / a.tissue
    out = Path(a.root) / "OUT" / a.tissue
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    res = {"tissue": a.tissue, "n_boot": a.n, "seed": a.seed, "per_d": {}}
    alpha_hat = {}
    for d in a.ds:
        p = wd / f"bins_{d:03d}um.npz"
        if not p.exists():
            print(f"缺 {p.name}，跳过", flush=True)
            continue
        X, n_cells, genes = load_bins(p)
        gs = [g.decode() if isinstance(g, bytes) else str(g) for g in genes]
        keep = ~np.array([s.upper().startswith("MT-") for s in gs])
        X = X[keep]                       # 与主流程 --exclude-mt 同口径
        X2 = X.multiply(X).tocsr()
        s1 = np.asarray(X.sum(axis=1)).ravel().astype(np.float64)
        s2 = np.asarray(X2.sum(axis=1)).ravel().astype(np.float64)
        a_hat, num, den = alpha_from_moments(s1, s2, n_cells)
        alpha_hat[d] = a_hat

        bs = np.empty(a.n)
        for i in range(a.n):
            bs[i] = boot_one(X, X2, rng)
        lo, hi = np.percentile(bs, [2.5, 97.5])
        res["per_d"][str(d)] = dict(
            n_cells=n_cells, n_obs=int(X.shape[1]), alpha_hat=a_hat,
            boot_mean=float(bs.mean()), boot_sd=float(bs.std(ddof=1)),
            ci95=[float(lo), float(hi)], cv=float(bs.std(ddof=1) / a_hat),
        )
        print(f"d={d:>3}um  bins={X.shape[1]:>6,}  alpha={a_hat:.4f}  "
              f"boot_sd={bs.std(ddof=1):.4f}  95%CI=[{lo:.4f}, {hi:.4f}]  "
              f"CV={bs.std(ddof=1)/a_hat*100:.1f}%", flush=True)
        del X, X2

    # 端点联合：2um 视为固定（抽样误差可忽略），只把 256um 的不确定性传过去
    if "2" in alpha_hat or True:
        a2 = None
        p2 = wd / "bins_002um.npz"
        if p2.exists():
            # 只取 2um 的 alpha 点估计；不重采样（矩阵太大）
            X2m, n2, genes2 = load_bins(p2)
            gs2 = [g.decode() if isinstance(g, bytes) else str(g) for g in genes2]
            k2 = ~np.array([s.upper().startswith("MT-") for s in gs2])
            X2m = X2m[k2]
            s1 = np.asarray(X2m.sum(axis=1)).ravel().astype(np.float64)
            s2 = np.asarray(X2m.multiply(X2m).sum(axis=1)).ravel().astype(np.float64)
            a2, _, _ = alpha_from_moments(s1, s2, n2)
            del X2m
        if a2 is not None and "256" in res["per_d"]:
            ci = res["per_d"]["256"]["ci95"]
            res["drop_fold"] = dict(
                point=a2 / alpha_hat[256],
                ci_from_alpha256=[a2 / ci[1], a2 / ci[0]],
                alpha_2um=a2,
            )
            print(f"\n下降倍数 点估计 {a2/alpha_hat[256]:.2f}  "
                  f"区间 [{a2/ci[1]:.2f}, {a2/ci[0]:.2f}]（只传 256um 误差）")

    (out / "boot_tail.json").write_text(
        json.dumps(res, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"\n写入 {out / 'boot_tail.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

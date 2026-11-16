r"""
39_block_boot.py — 块自助：gamma_pool 的抽样不确定度随"块多大才算独立"如何变化。

为什么做
--------
主稿用四带极差（0.131）与区域 SD（0.056）作为不确定度，并明确说"这是下界"。
审稿人会问：为什么不直接自助？因为 512 um 瓦片小于表达的相关长度，重采样它们等于
假装有上百个独立重复。本脚本把这句话做成一条曲线：块边长 512 / 1024 / 2048 um，
各 B=200 次，看 gamma_pool 的 SD 怎么随块变大而变化。

预期读法
--------
若 SD 随块边长单调上升并在 2 mm（约切片三分之一）附近接近区域 SD，则"一张切片
内部不存在可作为独立重复的尺度"这句话被量化，四带极差不再是一个随手取的数。

做法
----
对每个尺度 d、每个块边长 bs（bs 必须是 256 的整数倍，保证所有 d<=256 的粗格不被切断）：
  1. 用 dom_any（in_tissue 粗格 key）统计每个超块含多少个 d 尺度域格 -> n_b
  2. 用 X @ A 得到每超块每基因的计数和 S、平方和 Q（A 是 ncol x nblk 指示矩阵）
重采样：mult ~ Multinomial(nblk, 1/nblk)，
  S* = S @ mult, Q* = Q @ mult, n* = n_b @ mult
  alpha*(d) = sum_g[(Q*-S*^2/n*)/(n*-1) - S*/n*] / sum_g (S*/n*)^2
  gamma* = -slope(log alpha* ~ log d)

用法
----
python 39_block_boot.py --npzdir <dir> --tag colon --outdir <dir>
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


def block_stats(X, sid_obs, sid_all, chunk=1500):
    """返回 S (genes x nblk), Q (genes x nblk)。"""
    nblk = sid_all.size
    idx = np.searchsorted(sid_all, sid_obs)
    A = csr_matrix((np.ones(sid_obs.size, np.float32),
                    (np.arange(sid_obs.size, dtype=np.int32), idx.astype(np.int32))),
                   shape=(sid_obs.size, nblk))
    g = X.shape[0]
    S = np.zeros((g, nblk), np.float64)
    Q = np.zeros((g, nblk), np.float64)
    for a in range(0, g, chunk):
        b = min(a + chunk, g)
        Xc = X[a:b]
        S[a:b] = (Xc @ A).toarray()
        Q[a:b] = (Xc.multiply(Xc) @ A).toarray()
        del Xc
    return S, Q


def alpha_from(S, Q, n):
    mu = S / n
    var = (Q - S ** 2 / n) / (n - 1)
    num = float(np.sum(var - mu))
    den = float(np.sum(mu ** 2))
    return num / den


def fit_gamma(ds, aa):
    ds = np.asarray(ds, float)
    aa = np.asarray(aa, float)
    ok = np.isfinite(aa) & (aa > 0)
    c = np.polyfit(np.log(ds[ok]), np.log(aa[ok]), 1)
    return float(-c[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npzdir", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--outdir", default=str(ROOT / "02_OUT" / "block_boot"))
    ap.add_argument("--bs", nargs="+", type=int, default=[512, 1024, 2048])
    ap.add_argument("--B", type=int, default=200)
    ap.add_argument("--seed", type=int, default=20261114)
    a = ap.parse_args()
    od = Path(a.outdir)
    od.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    # ---- 一次读入所有尺度的块统计，按块边长分开存
    store = {bs: dict(S=[], Q=[], n=[], Sfull=None) for bs in a.bs}
    print(f"[{a.tag}] 各尺度域格统计 ...", flush=True)
    for d in DS:
        X, z = load(Path(a.npzdir) / f"bins_{d:03d}um.npz")
        genes = np.array([str(s) for s in z["genes"]])
        keep = ~np.array([s.upper().startswith("MT-") for s in genes])
        if (~keep).sum():
            X = X[keep].tocsr()
        coords = z["coords"].astype(np.int64)
        # coords 已经是 d 尺度网格索引；dom_any 的关键字同样是 d 尺度 (r,c) 打包
        r0, c0 = coords[:, 0], coords[:, 1]
        domk = z["dom_any"].astype(np.int64)
        dr, dc = domk // STRIDE, domk % STRIDE

        for bs in a.bs:
            fb = bs // d
            if fb < 1:
                raise ValueError(f"block {bs} 小于尺度 {d}")
            sid_obs = (r0 // fb) * STRIDE + (c0 // fb)
            sid_dom = (dr // fb) * STRIDE + (dc // fb)
            sid_all = np.unique(np.concatenate([sid_obs, sid_dom]))
            n_b = np.bincount(np.searchsorted(sid_all, sid_dom),
                              minlength=sid_all.size).astype(np.float64)
            S, Q = block_stats(X, sid_obs, sid_all)
            store[bs]["S"].append(S)
            store[bs]["Q"].append(Q)
            store[bs]["n"].append(n_b)
            print(f"  d={d:3d} bs={bs:4d}: 超块 {sid_all.size:4d}  "
                  f"域格/块 中位 {np.median(n_b):8.0f}", flush=True)
        del X, z
        gc.collect()

    out = dict(tag=a.tag, B=a.B, block_sizes=a.bs, ds=DS, results={})
    for bs in a.bs:
        st = store[bs]
        n_b = np.vstack(st["n"])          # (8, nblk)
        nblk = n_b.shape[1]
        p = 1.0 / nblk
        gams = []
        for _ in range(a.B):
            mult = rng.multinomial(nblk, np.full(nblk, p)) if nblk < 40 else \
                rng.multinomial(nblk, np.full(nblk, p))
            nn = n_b @ mult
            if np.any(nn <= 1):
                continue
            al = []
            for j in range(len(DS)):
                al.append(alpha_from(st["S"][j] @ mult, st["Q"][j] @ mult, nn[j]))
            gams.append(fit_gamma(DS, al))
        g = np.array(gams)
        out["results"][str(bs)] = dict(
            nblk=int(nblk), n_valid=int(g.size),
            mean=float(g.mean()), sd=float(g.std(ddof=1)),
            q025=float(np.quantile(g, 0.025)), q975=float(np.quantile(g, 0.975)),
            ci_width=float(np.quantile(g, 0.975) - np.quantile(g, 0.025)),
            gamma=g.tolist())
        print(f"  bs={bs:4d} nblk={nblk:4d}  SD={g.std(ddof=1):.4f}  "
              f"95%CI=[{np.quantile(g,0.025):.3f},{np.quantile(g,0.975):.3f}]", flush=True)

    json.dump(out, open(od / f"boot_{a.tag}.json", "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("写出", od / ("boot_" + a.tag + ".json"))


if __name__ == "__main__":
    main()

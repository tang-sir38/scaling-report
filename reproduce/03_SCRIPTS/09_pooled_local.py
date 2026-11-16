r"""
09_pooled_local.py — 局部窗口内的**池化**过散度（补齐口径不一致）。

问题
----
主文用的估计量是池化的：alpha_pool = SUM_g(Var_g - mu_g) / SUM_g mu_g^2，
不筛选基因。但 04_local_alpha.py 仍是逐基因口径（要求每个基因在每档都有足够计数），
两套数不能直接写在正文里互相引用。

本脚本把 04 的"256 um 方块内方差"接到 07 的"池化"上：
  var_g^local = SUM_b (S2_gb - S1_gb^2/n_b) / SUM_b (n_b - 1)
  mu_g        = TOT_g / N_domain
  alpha_pool^local = SUM_g (var_g^local - mu_g) / SUM_g mu_g^2

与全局版的唯一差别是 Var 换成了"块内方差"，即已扣除块间（域）差异。

用法
----
python 09_pooled_local.py --workdir <root>\01_WORK --outdir <root>\02_OUT --ds 8 16 32 64 --block 256
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

ROOT = Path(__file__).resolve().parents[1]

KEY_STRIDE = 4096


def load_npz(p):
    z = np.load(p, allow_pickle=True)
    X = csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(int(x) for x in z["shape"]))
    return X, z["coords"], z


def block_sums(X, coords, dom_rc, block_um, d, chunk=2000):
    """返回每基因的 (num, tot, den, n_tot)。

    num_g = SUM_b (S2_gb - S1_gb^2/n_b)；den = SUM_b (n_b - 1)
    """
    reps = max(1, block_um // d)
    k_obs = (coords[:, 0].astype(np.int64) // reps) * KEY_STRIDE + (coords[:, 1].astype(np.int64) // reps)
    k_dom = (dom_rc[:, 0].astype(np.int64) // reps) * KEY_STRIDE + (dom_rc[:, 1].astype(np.int64) // reps)

    ublk, n_b = np.unique(k_dom, return_counts=True)
    nb_map = dict(zip(ublk.tolist(), n_b.tolist()))
    pos = np.searchsorted(ublk, k_obs)
    pos_c = np.clip(pos, 0, ublk.size - 1)
    ok = ublk[pos_c] == k_obs
    n_blk = ublk.size

    B = csr_matrix(
        (np.ones(int(ok.sum()), dtype=np.float64), (pos_c[ok], np.nonzero(ok)[0])),
        shape=(n_blk, coords.shape[0]),
    )
    nb_vec = np.array([nb_map.get(int(k), 0) for k in ublk], dtype=np.float64)
    keep = nb_vec > 1
    den = float(nb_vec[keep].sum() - keep.sum())
    w = nb_vec[keep]

    G = X.shape[0]
    num = np.zeros(G, dtype=np.float64)
    tot = np.zeros(G, dtype=np.float64)
    for a in range(0, G, chunk):
        b = min(a + chunk, G)
        Xc = X[a:b].astype(np.float64)
        S1 = np.asarray((Xc @ B.T).todense())[:, keep]
        X2 = Xc.copy()
        X2.data **= 2
        S2 = np.asarray((X2 @ B.T).todense())[:, keep]
        num[a:b] = (S2 - S1 ** 2 / w).sum(axis=1)
        tot[a:b] = np.asarray(Xc.sum(axis=1)).ravel()
        del Xc, X2, S1, S2
    n_tot = float(nb_vec.sum())
    return num, tot, den, n_tot, int(keep.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--ds", nargs="+", type=int, default=[8, 16, 32, 64])
    ap.add_argument("--block", type=int, default=256, help="固定方块边长（um）")
    ap.add_argument("--block-mult", type=int, default=0,
                    help="若 >0，方块边长 = mult * d（尺度不变口径，推荐）")
    ap.add_argument("--domain", default="any", choices=["any", "full"])
    ap.add_argument("--exclude-mt", action="store_true")
    a = ap.parse_args()

    wd, od = Path(a.workdir), Path(a.outdir)
    od.mkdir(parents=True, exist_ok=True)

    rows, keep = {}, None
    for d in a.ds:
        p = wd / f"bins_{d:03d}um.npz"
        if not p.exists():
            print(f"  缺 {p.name}，跳过")
            continue
        X, coords, z = load_npz(p)
        dom = z[f"dom_{a.domain}"].astype(np.int64)
        dom_rc = np.column_stack([dom // KEY_STRIDE, dom % KEY_STRIDE]).astype(np.int32)
        blk = a.block_mult * d if a.block_mult > 0 else a.block
        num, tot, den, n_tot, nblk = block_sums(X, coords, dom_rc, blk, d)

        if a.exclude_mt:
            if keep is None:
                gs = [g.decode() if isinstance(g, bytes) else str(g) for g in z["genes"]]
                keep = ~np.array([s.upper().startswith("MT-") for s in gs])
            num, tot = num[keep], tot[keep]

        var = num / den
        mu = tot / n_tot
        alpha_pool = float(np.sum(var - mu) / np.sum(mu ** 2))
        rows[d] = dict(alpha_pool=alpha_pool, n_blocks=nblk, n_domain=n_tot,
                       total_umi=float(tot.sum()))
        print(f"  d={d:3d}um  方块={blk:5d}um ({nblk:4d} 个)  域内格={int(n_tot):>8d}  "
              f"alpha_pool_local={alpha_pool:.5f}", flush=True)
        del X, num, tot

    ds = sorted(rows)
    if len(ds) < 3:
        print("分辨率不足 3 档")
        return
    x = np.log(np.array(ds, float))
    y = np.array([np.log(rows[d]["alpha_pool"]) for d in ds])
    slope = float(np.polyfit(x, y, 1)[0])
    seg = []
    for i in range(len(ds) - 1):
        r = rows[ds[i + 1]]["alpha_pool"] / rows[ds[i]]["alpha_pool"]
        seg.append(dict(frm=ds[i], to=ds[i + 1], ratio=float(r),
                        slope=float(np.log(r) / np.log(ds[i + 1] / ds[i]))))
        print(f"  {ds[i]:>3}->{ds[i+1]:>3}um  比值 {r:.3f}  斜率 {seg[-1]['slope']:+.3f}")
    print(f"\n  == 池化局部斜率（{'方块=' + str(a.block_mult) + 'x d' if a.block_mult else '方块=' + str(a.block) + 'um'}，"
          f"d={ds[0]}-{ds[-1]}）= {slope:+.3f} ==")

    out = dict(block_um=(None if a.block_mult else a.block), block_mult=a.block_mult,
               domain=a.domain, ds=ds, per_d=rows,
               slope=slope, segments=seg, exclude_mt=bool(a.exclude_mt))
    suffix = "_noMT" if a.exclude_mt else ""
    tag = f"{a.block_mult}xd" if a.block_mult else f"{a.block}um"
    fn = f"pooled_local_{tag}{suffix}.json"
    json.dump(out, open(od / fn, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(f"已写出 {od/fn}")


if __name__ == "__main__":
    main()

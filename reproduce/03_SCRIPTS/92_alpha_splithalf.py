r"""92_alpha_splithalf.py — 把两半各自的逐基因 alpha 曲线落盘

问题
----
91_window_scan 显示：gamma 的拟合窗口会显著改变与 beta 的相关
（2-256 时 r=0.70-0.78，16-256 时 r=0.91-0.94）。
若直接选后者，就是"用结果挑窗口"，审稿人必打。
需要一个**与 beta 完全无关**的窗口判据。

解法
----
用两半可靠性。44_splithalf_genes.py 已经证明：把切片按列劈成两半，
各自独立算 gamma，两半之间的相关就是"测量可靠性"，
它只反映估计噪声，与 beta 毫无关系。
但 44 只存了全窗口的 gamma_a/gamma_b，没存中间量，无法算别的窗口。

本脚本复用 44 的精确划分与聚合口径（STRIDE、DOMKEY、列中点切、moments 全同），
只多存一步：两半各自的 alpha 曲线 (G,8)。之后任意窗口的可靠性都可复算。

口径自检
--------
对 colon，用存下的 AL_a/AL_b 按 2-256 重算 gamma 并核对 44 的 r=0.7925 与 n=11805；
不一致即停（防止悄悄换了口径）。

用法
----
  python 92_alpha_splithalf.py --tag colon
  python 92_alpha_splithalf.py --tag lung
  python 92_alpha_splithalf.py --tag brain
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[1]
DS = [2, 4, 8, 16, 32, 64, 128, 256]
STRIDE = 4096

SETS = {
    "colon": (ROOT / "01_WORK/bins_002um.npz", ROOT / "02_OUT"),
    "lung": (ROOT / "D3_LungCancer/01_WORK/bins_002um.npz", ROOT / "D3_LungCancer/02_OUT"),
    "brain": (ROOT / "D2_MouseBrain/01_WORK/bins_002um.npz", ROOT / "D2_MouseBrain/02_OUT"),
}


def moments(Xa, n):
    s1 = np.asarray(Xa.sum(axis=1)).ravel().astype(np.float64)
    s2 = np.asarray(Xa.multiply(Xa).sum(axis=1)).ravel().astype(np.float64)
    mu = s1 / n
    var = (s2 - s1 ** 2 / n) / (n - 1)
    return mu, var


def fit_gamma_win(al, ds_win):
    """与 41/44 同式的 OLS：log alpha ~ log d，需窗口内全正有限。"""
    j = np.array([DS.index(int(d)) for d in ds_win])
    sub = al[:, j]
    ok = np.all(np.isfinite(sub) & (sub > 0), axis=1)
    x = np.log(np.array(ds_win, float))
    xm = x - x.mean()
    y = np.log(np.where(ok[:, None], sub, 1.0))
    ym = y - y.mean(axis=1, keepdims=True)
    g = -np.sum(xm[None, :] * ym, axis=1) / np.sum(xm ** 2)
    out = np.full(al.shape[0], np.nan)
    out[ok] = g[ok]
    return out, ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="colon", choices=list(SETS))
    ap.add_argument("--outdir", default=None)
    a = ap.parse_args()

    npz_p, default_out = SETS[a.tag]
    od = Path(a.outdir) if a.outdir else (ROOT / "02_OUT/variogram")
    od.mkdir(parents=True, exist_ok=True)

    z = np.load(npz_p, allow_pickle=False)
    shape = tuple(int(x) for x in z["shape"])
    X = csr_matrix((z["data"].astype(np.float32), z["indices"], z["indptr"]), shape=shape)
    coords = z["coords"].astype(np.int64)
    genes = np.array([str(s) for s in z["genes"]])
    keep = ~np.array([s.upper().startswith("MT-") for s in genes])
    X = X[keep].tocsr()
    gs = genes[keep]
    domk = z["dom_any"].astype(np.int64)
    R, C = coords[:, 0], coords[:, 1]
    Rd, Cd = domk // STRIDE, domk % STRIDE
    del z

    AGG, DOMKEY = {}, {}
    for d in DS:
        f = d // 2
        key = (R // f) * 100000 + (C // f)
        uk, inv = np.unique(key, return_inverse=True)
        AGG[d] = csr_matrix((np.ones(key.size, np.float32),
                             (np.arange(key.size), inv.astype(np.int32))),
                            shape=(key.size, uk.size))
        DOMKEY[d] = (Rd // f) * 100000 + (Cd // f)
        del uk, inv
    print(f"[{a.tag}] {X.shape[0]} genes x {X.shape[1]} bins", flush=True)

    c_lo, c_hi = int(C.min()), int(C.max())
    cut = int(np.floor(((c_lo + c_hi) / 2) / 128) * 128)

    def run(mb, md):
        Xb = X[:, mb]
        AL = np.empty((X.shape[0], len(DS)))
        MU = np.empty((X.shape[0], len(DS)))
        nb = np.zeros(len(DS), np.int64)
        for j, d in enumerate(DS):
            Xa = Xb if d == 2 else Xb @ AGG[d][mb]
            n = int(np.unique(DOMKEY[d][md]).size)
            nb[j] = n
            mu, var = moments(Xa, n)
            MU[:, j] = mu
            with np.errstate(invalid="ignore", divide="ignore"):
                AL[:, j] = (var - mu) / mu ** 2
            if d != 2:
                del Xa
        return AL, MU, nb

    print(f"  切点 col={cut}；左 {int((C < cut).sum()):,} / 右 {int((C >= cut).sum()):,} bin",
          flush=True)
    AL_a, MU_a, nb_a = run(C < cut, Cd < cut)
    print("  左半完成", flush=True)
    AL_b, MU_b, nb_b = run(C >= cut, Cd >= cut)
    print("  右半完成", flush=True)

    # ---- 口径自检：全窗口可靠性必须与 44 完全一致
    chk = None
    if a.tag == "colon":
        ga, oka = fit_gamma_win(AL_a, DS)
        gb, okb = fit_gamma_win(AL_b, DS)
        both = oka & okb
        r = float(pearsonr(ga[both], gb[both]).statistic)
        ref = json.loads((ROOT / "02_OUT/splithalf/splithalf.json")
                         .read_text(encoding="utf-8"))["reliability"]
        chk = dict(n=int(both.sum()), pearson=r, ref_n=int(ref["n"]),
                   ref_pearson=float(ref["pearson"]))
        print(f"  自检 n={chk['n']} (ref {chk['ref_n']})  r={r:.6f} (ref {chk['ref_pearson']:.6f})")
        if chk["n"] != chk["ref_n"] or abs(r - chk["ref_pearson"]) > 1e-6:
            raise SystemExit("与 44 口径不一致，停止")

    np.savez_compressed(
        od / f"splithalf_alpha_{a.tag}.npz",
        gene=gs.astype("U16"), ds=np.array(DS),
        al_a=AL_a.astype(np.float32), al_b=AL_b.astype(np.float32),
        mu_a=MU_a.astype(np.float32), mu_b=MU_b.astype(np.float32),
        nb_a=nb_a, nb_b=nb_b,
    )
    meta = dict(tag=a.tag, cut=int(cut), n_genes=int(gs.size),
                n_bins_left=int((C < cut).sum()), n_bins_right=int((C >= cut).sum()),
                n_bins_per_d_left=nb_a.tolist(), n_bins_per_d_right=nb_b.tolist(),
                selfcheck=chk)
    (od / f"splithalf_alpha_{a.tag}.json").write_text(
        json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"写出 {od}/splithalf_alpha_{a.tag}.npz", flush=True)


if __name__ == "__main__":
    main()

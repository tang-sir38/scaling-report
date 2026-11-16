r"""
44_splithalf_genes.py — 两半之间逐基因 gamma 的可靠性，落盘点本身，并按表达层分层。

37_splithalf.py 只存了汇总（r = 0.79，n = 11,805），图需要点，
审稿意见 ⑤（"低表达基因的 gamma 就是噪声"）需要分层答案。
本脚本重跑同一套两半划分，但把逐基因 gamma 与 mu32 一起存下来，
并核对汇总与 37 一致（否则说明口径变了，直接停）。

输出 <outdir>/splithalf_genes.npz 与 reliability_layers.json
"""
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[1]
NPZ = ROOT / "01_WORK/bins_002um.npz"
OUT = ROOT / "02_OUT/splithalf_genes"
DS = [2, 4, 8, 16, 32, 64, 128, 256]
STRIDE = 4096
EDGES = np.array([0.01, 0.0316, 0.1, 0.316, 1.0, 3.16, 10.0, 31.6, 100.0])


def moments(Xa, n):
    s1 = np.asarray(Xa.sum(axis=1)).ravel().astype(np.float64)
    s2 = np.asarray(Xa.multiply(Xa).sum(axis=1)).ravel().astype(np.float64)
    mu = s1 / n
    var = (s2 - s1 ** 2 / n) / (n - 1)
    return mu, var


def fit_gamma(al):
    ok = np.all(np.isfinite(al) & (al > 0), axis=1)
    x = np.log(np.array(DS, float))
    xm = x - x.mean()
    y = np.log(np.where(ok[:, None], al, 1.0))
    ym = y - y.mean(axis=1, keepdims=True)
    g = -np.sum(xm[None, :] * ym, axis=1) / np.sum(xm ** 2)
    out = np.full(al.shape[0], np.nan)
    out[ok] = g[ok]
    return out, ok


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    z = np.load(NPZ, allow_pickle=False)
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

    AGG, DOMKEY = {}, {}
    for d in DS:
        f = d // 2
        key = R * 0 + 0  # placeholder to keep linters quiet
        key = (R // f) * 100000 + (C // f)
        uk, inv = np.unique(key, return_inverse=True)
        AGG[d] = csr_matrix((np.ones(key.size, np.float32),
                             (np.arange(key.size), inv.astype(np.int32))),
                            shape=(key.size, uk.size))
        DOMKEY[d] = (Rd // f) * 100000 + (Cd // f)
        del uk, inv
    print(f"读入 {X.shape[0]} genes x {X.shape[1]} bins", flush=True)

    c_lo, c_hi = int(C.min()), int(C.max())
    cut = int(np.floor(((c_lo + c_hi) / 2) / 128) * 128)

    def run(mb, md):
        Xb = X[:, mb]
        AL = np.empty((X.shape[0], len(DS)))
        for j, d in enumerate(DS):
            Xa = Xb if d == 2 else Xb @ AGG[d][mb]
            n = int(np.unique(DOMKEY[d][md]).size)
            mu, var = moments(Xa, n)
            with np.errstate(invalid="ignore", divide="ignore"):
                AL[:, j] = (var - mu) / mu ** 2
            if d != 2:
                del Xa
        return AL

    print(f"切点 col={cut}；左半 {int((C < cut).sum()):,} bin，右半 {int((C >= cut).sum()):,}", flush=True)
    AL_a = run(C < cut, Cd < cut)
    AL_b = run(C >= cut, Cd >= cut)
    ga, oka = fit_gamma(AL_a)
    gb, okb = fit_gamma(AL_b)
    both = oka & okb

    r = float(pearsonr(ga[both], gb[both]).statistic)
    rho = float(spearmanr(ga[both], gb[both]).statistic)
    mad = float(np.median(np.abs(ga[both] - gb[both])))
    print(f"逐基因可靠性 n={int(both.sum()):,} r={r:.4f} rho={rho:.4f} MAD={mad:.4f}")
    ref = json.loads((ROOT / "02_OUT/splithalf/splithalf.json").read_text(encoding="utf-8"))["reliability"]
    print(f"  lock: n={ref['n']} r={ref['pearson']:.4f} MAD={ref['mad']:.4f}")
    if int(both.sum()) != int(ref["n"]) or abs(r - ref["pearson"]) > 1e-6:
        raise SystemExit("与 37_splithalf.py 的口径不一致，停止")

    # 32 um 的 mu（与逐基因曲线同一把尺）
    X32 = X @ AGG[32]
    n32 = int(np.unique(DOMKEY[32][np.ones(Cd.size, bool)]).size)
    mu32 = np.asarray(X32.sum(axis=1)).ravel() / n32

    lay = np.full(ga.size, -1, np.int8)
    for i in range(EDGES.size - 1):
        m = (mu32 >= EDGES[i]) & (mu32 < EDGES[i + 1])
        lay[m] = i

    out = dict(n=int(both.sum()), pearson=r, spearman=rho, mad=mad,
               layers=[], note="层边界同 Table 3；r 与 SD 均在 log gamma 上算？否：线性")
    for i in range(EDGES.size - 1):
        m = both & (lay == i)
        if m.sum() < 30:
            out["layers"].append(dict(i=i, lo=float(EDGES[i]), hi=float(EDGES[i + 1]),
                                      n=int(m.sum()), pearson=None, mad=None))
            continue
        rr = float(pearsonr(ga[m], gb[m]).statistic)
        out["layers"].append(dict(i=i, lo=float(EDGES[i]), hi=float(EDGES[i + 1]),
                                  n=int(m.sum()), pearson=rr,
                                  mad=float(np.median(np.abs(ga[m] - gb[m]))),
                                  med_a=float(np.median(ga[m])), med_b=float(np.median(gb[m]))))
        print(f"  层 {EDGES[i]:.3g}-{EDGES[i+1]:.3g}: n={int(m.sum()):5d}  r={rr:.3f}  "
              f"MAD={out['layers'][-1]['mad']:.3f}")

    np.savez_compressed(OUT / "splithalf_genes.npz",
                        gene=gs.astype("U16"), mu32=mu32.astype(np.float32),
                        gamma_a=ga.astype(np.float32), gamma_b=gb.astype(np.float32),
                        both=both, layer=lay)
    json.dump(out, open(OUT / "reliability_layers.json", "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print(f"写出 {OUT/'splithalf_genes.npz'} 与 reliability_layers.json")


if __name__ == "__main__":
    main()

r"""
37_splithalf.py — 同一张切片内，把组织按列切成若干带，各自重算 gamma。

为什么这是最硬的一击
--------------------
原稿的核心比较是"每个组织一个样本，比较 gamma"。若 gamma 是组织性质，同一张切片
内部（细胞组成相近、只差区域）的 gamma 应当基本一致；若 gamma 主要由"哪些转录本
占据 sum(mu^2)"决定，则同一张切片的不同区域就可能给出不同的 gamma——而这些区域的
组织结构并没有本质差别。

做法
----
按 256 um 网格线把域切成竖带（256 um = 128 个 2 um bin，切点取 128 的整数倍，
保证所有尺度 d<=256（f=d/2<=128）的粗格都不被切断）。
每条带单独算逐基因矩 -> pooled alpha(8 档) -> gamma。
附带：两半之间逐基因 gamma 的可靠性（Pearson/Spearman/MAD）。

输出 <outdir>/splithalf.json
"""
import argparse, json, gc
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
from scipy.stats import spearmanr, pearsonr

ROOT = Path(__file__).resolve().parents[1]


def load_npz(p):
    z = np.load(p, allow_pickle=False)
    shape = tuple(int(x) for x in z["shape"])
    X = csr_matrix((z["data"].astype(np.float32), z["indices"], z["indptr"]), shape=shape)
    dom_any = z["dom_any"] if "dom_any" in z.files else None
    return X, z["coords"], z["genes"], z["n_domain_any"], dom_any


def per_gene_moments(Xa, n_cells):
    s1 = np.asarray(Xa.sum(axis=1)).ravel().astype(np.float64)
    s2 = np.asarray(Xa.multiply(Xa).sum(axis=1)).ravel().astype(np.float64)
    mu = s1 / n_cells
    var = (s2 - s1 ** 2 / n_cells) / (n_cells - 1) if n_cells > 1 else np.zeros_like(mu)
    return mu, var


def fit_gamma(ds, al):
    ok = np.all(np.isfinite(al) & (al > 0), axis=1)
    x = np.log(ds)
    xm = x - x.mean()
    y = np.log(np.where(ok[:, None], al, 1.0))
    ym = y - y.mean(axis=1, keepdims=True)
    g = -np.sum(xm[None, :] * ym, axis=1) / np.sum(xm ** 2)
    out = np.full(al.shape[0], np.nan)
    out[ok] = g[ok]
    return out, ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", default=str(ROOT / "01_WORK" / "bins_002um.npz"))
    ap.add_argument("--outdir", default=str(ROOT / "02_OUT" / "splithalf"))
    ap.add_argument("--ds", nargs="+", type=int, default=[2, 4, 8, 16, 32, 64, 128, 256])
    ap.add_argument("--stride", type=int, default=4096)
    a = ap.parse_args()
    od = Path(a.outdir)
    od.mkdir(parents=True, exist_ok=True)

    X, coords, genes, n_dom2, dom_any = load_npz(a.npz)
    gs = np.array([str(g) for g in genes])
    keep = ~np.array([s.upper().startswith("MT-") for s in gs])
    X = X[keep].tocsr()
    gs = gs[keep]
    # coords 是二维 (row, col)；dom_any 是扁平 key = row*W + col，W=4096（与 02b 一致）
    R = coords[:, 0].astype(np.int64)
    C = coords[:, 1].astype(np.int64)
    Rd = (dom_any // a.stride).astype(np.int64)
    Cd = (dom_any % a.stride).astype(np.int64)
    print(f"读入 {X.shape[0]} genes x {X.shape[1]} bins；列 {C.min()}..{C.max()}", flush=True)

    ds = a.ds
    AGG, DOMKEY = {}, {}
    for d in ds:
        f = d // 2
        key = (R // f) * 100000 + (C // f)
        uk, inv = np.unique(key, return_inverse=True)
        AGG[d] = csr_matrix((np.ones(key.size, np.float32),
                             (np.arange(key.size), inv.astype(np.int32))),
                            shape=(key.size, uk.size))
        DOMKEY[d] = (Rd // f) * 100000 + (Cd // f)
        del uk, inv
        print(f"  d={d:3d} 观测粗格 {AGG[d].shape[1]:8d}  域格 {np.unique(DOMKEY[d]).size:8d}", flush=True)

    def run_band(mb, md):
        """mb: 2um bin 的布尔；md: 域格 key 的布尔。返回 pooled gamma、alpha、逐基因 alpha。"""
        Xb = X[:, mb]
        EX = np.empty((X.shape[0], len(ds)))
        MU2 = np.empty_like(EX)
        for j, d in enumerate(ds):
            if d == 2:
                Xa = Xb
            else:
                Xa = Xb @ AGG[d][mb]
            n = int(np.unique(DOMKEY[d][md]).size)
            mu, var = per_gene_moments(Xa, n)
            EX[:, j] = var - mu
            MU2[:, j] = mu ** 2
            if d != 2:
                del Xa
                gc.collect()
        ap_ = EX.sum(0) / MU2.sum(0)
        c = np.polyfit(np.log(ds), np.log(ap_), 1)
        with np.errstate(invalid="ignore", divide="ignore"):
            AL = EX / MU2
        del Xb
        gc.collect()
        return float(-c[0]), ap_.tolist(), AL

    out = dict(source=str(a.npz), ds=ds, genes=int(X.shape[0]), bands={})

    c_lo, c_hi = int(C.min()), int(C.max())
    cut = int(np.floor(((c_lo + c_hi) / 2) / 128) * 128)
    if cut <= c_lo:
        cut = c_lo + 128
    print(f"两半切点 col={cut}", flush=True)

    ALA = {}
    for nm, mb, md in (("halfA", C < cut, Cd < cut), ("halfB", C >= cut, Cd >= cut)):
        g, ap_, AL = run_band(mb, md)
        ALA[nm] = AL
        out["bands"][nm] = dict(gamma=g, alpha=ap_, n_bin=int(mb.sum()),
                                n_dom=[int(np.unique(DOMKEY[d][md]).size) for d in ds])
        print(f"  {nm}: n_bin={int(mb.sum()):,}  gamma={g:.4f}", flush=True)

    ga, oka = fit_gamma(np.array(ds, float), ALA["halfA"])
    gb, okb = fit_gamma(np.array(ds, float), ALA["halfB"])
    both = oka & okb
    out["reliability"] = dict(
        n=int(both.sum()),
        pearson=float(pearsonr(ga[both], gb[both]).statistic),
        spearman=float(spearmanr(ga[both], gb[both]).statistic),
        mad=float(np.median(np.abs(ga[both] - gb[both]))),
        med_a=float(np.median(ga[both])), med_b=float(np.median(gb[both])))
    r = out["reliability"]
    print(f"  逐基因可靠性 n={r['n']} pearson={r['pearson']:.3f} spearman={r['spearman']:.3f}"
          f" MAD={r['mad']:.3f}", flush=True)
    del ALA
    gc.collect()

    edges = np.floor(np.linspace(c_lo, c_hi + 1, 5) / 128) * 128
    for i in range(4):
        lo, hi = int(edges[i]), int(edges[i + 1])
        mb = (C >= lo) & (C < hi)
        md = (Cd >= lo) & (Cd < hi)
        if mb.sum() < 1000:
            continue
        g, ap_, _ = run_band(mb, md)
        out["bands"][f"q{i+1}"] = dict(gamma=g, alpha=ap_, n_bin=int(mb.sum()),
                                       n_dom=[int(np.unique(DOMKEY[d][md]).size) for d in ds])
        print(f"  q{i+1} [{lo},{hi}): n_bin={int(mb.sum()):,}  gamma={g:.4f}", flush=True)

    gq = [v["gamma"] for k, v in out["bands"].items() if k.startswith("q")]
    if len(gq) > 1:
        out["quarter_range"] = float(max(gq) - min(gq))
        out["quarter_sd"] = float(np.std(gq, ddof=1))
        print(f"  四等分极差 {out['quarter_range']:.4f}  SD {out['quarter_sd']:.4f}", flush=True)

    json.dump(out, open(od / "splithalf.json", "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print(f"\n写出 {od/'splithalf.json'}")


if __name__ == "__main__":
    main()

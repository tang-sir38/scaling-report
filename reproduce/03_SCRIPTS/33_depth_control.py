r"""
33_depth_control.py — 深度对照：把 colon 稀释到 brain 的测序深度，gamma 会变吗？

为什么需要
----------
审稿人一定会问：脑的指数更浅（0.36 vs 0.55），是不是因为它的 2 um bin 只有 14 UMI
（结肠 30.5、肺 24.5）？低计数会不会把 alpha 压低、把曲线压平？

理论预期
--------
若计数是二项稀释（每个 UMI 以概率 p 保留），则 mu -> p*mu，Var -> p*mu + p^2*Var_bio，
于是 Var - mu -> p^2*(Var_bio)，mu^2 -> p^2*mu^2，alpha = (Var-mu)/mu^2 **不变**。
所以 gamma 应在稀释下严格不变；任何变化都来自估计偏差或基因筛除。

做法
----
对 2 um 矩阵做 binomial thinning（p = 目标深度/现有深度），重新聚合到 8 个尺度，
用与主分析完全相同的 alpha_pool 口径重算 gamma。跑 3 档：p=0.46（脑的深度比）、0.25、0.10。

输出 <outdir>/depth_control.json
"""
import argparse, json
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix

ROOT = Path(__file__).resolve().parents[1]

def load_npz(p):
    z = np.load(p, allow_pickle=False)
    shape = tuple(int(x) for x in z["shape"])
    X = csr_matrix((z["data"].astype(np.float32), z["indices"], z["indptr"]), shape=shape)
    dom_any = z["dom_any"] if "dom_any" in z.files else None
    return X, z["coords"], z["genes"], z["n_domain_any"], dom_any

def domain_blocks(dom_any, factor, stride=4096):
    """由 in_tissue 的 2um 粗格 key 推 d 尺度的域格数（含全零粗格）。"""
    r = (dom_any // stride).astype(np.int64)
    c = (dom_any % stride).astype(np.int64)
    k = (r // factor) * 100000 + (c // factor)
    return int(np.unique(k).size)

def agg_matrix(coords, factor, stride=8192):
    """返回 (n_bin, n_block) 的聚合矩阵：X @ A -> X_agg。"""
    key = (coords[:, 0] // factor).astype(np.int64) * stride + (coords[:, 1] // factor)
    uk, inv = np.unique(key, return_inverse=True)
    n_bin = coords.shape[0]
    A = csr_matrix((np.ones(n_bin, np.float32), (np.arange(n_bin, dtype=np.int32), inv.astype(np.int32))),
                   shape=(n_bin, uk.size))
    return A, uk.size

def pooled(X, n_cells):
    s1 = np.asarray(X.sum(axis=1)).ravel().astype(np.float64)
    s2 = np.asarray(X.multiply(X).sum(axis=1)).ravel().astype(np.float64)
    mu = s1 / n_cells
    var = (s2 - s1 ** 2 / n_cells) / (n_cells - 1)
    num = float(np.sum(var - mu)); den = float(np.sum(mu ** 2))
    return num / den if den > 0 else np.nan, float(s1.sum())

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", default=str(ROOT / "01_WORK" / "bins_002um.npz"))
    ap.add_argument("--outdir", default=str(ROOT / "02_OUT" / "depth_control"))
    ap.add_argument("--ps", nargs="+", type=float, default=[0.46, 0.25, 0.10])
    ap.add_argument("--ds", nargs="+", type=int, default=[2, 4, 8, 16, 32, 64, 128, 256])
    ap.add_argument("--seed", type=int, default=99)
    a = ap.parse_args()
    od = Path(a.outdir); od.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    X0, coords, genes, n_dom2, dom_any = load_npz(a.npz)
    # 剔除线粒体基因，与主分析同口径
    gs = np.array([str(g) for g in genes])
    keep = ~np.array([s.upper().startswith("MT-") for s in gs])
    if (~keep).sum():
        X0 = X0[keep].tocsr()
        print(f"  剔除 MT 基因 {int((~keep).sum())} 个，保留 {int(keep.sum())} 个", flush=True)
    print(f"读入 {X0.shape[0]} genes x {X0.shape[1]} bins, nnz={X0.nnz:,}, 2um 域格 {n_dom2:,}", flush=True)
    # 各尺度的域格数（必须用 in_tissue 定义，含全零粗格；否则低尺度 alpha 会偏）
    dom = {}
    A = {}
    for d in a.ds:
        f = d // 2
        Ai, nb_obs = agg_matrix(coords, f)
        A[d] = Ai
        if dom_any is not None:
            dom[d] = domain_blocks(dom_any, f)
        else:
            dom[d] = nb_obs
            print(f"  警告: {a.npz} 无 dom_any，d={d} 用观测粗格数 {nb_obs:,}（可能偏小）", flush=True)
    print("各尺度域格数:", {d: int(dom[d]) for d in a.ds}, flush=True)

    out = dict(source=str(a.npz), n_domain={str(d): int(dom[d]) for d in a.ds}, runs={})

    def run(tag, p):
        if p >= 1.0:
            X = X0.copy()
        else:
            X = X0.copy()
            ip = X.indptr; ix = X.indices; dt = X.data
            new = rng.binomial(dt.astype(np.int64), p).astype(np.float32)
            X.data = new
            X.eliminate_zeros()
        al = {}
        for d in a.ds:
            Xa = X if d == 2 else (X @ A[d])
            av, tot = pooled(Xa, dom[d])
            al[d] = float(av)
            print(f"    [{tag}] d={d:3d}  alpha={av:.4f}  UMI={tot:.3e}", flush=True)
            if d != 2:
                del Xa
        dd = np.array([d for d in a.ds if al[d] > 0], float)
        aa = np.array([al[d] for d in a.ds if al[d] > 0], float)
        ly = np.log(aa)
        cff = np.polyfit(np.log(dd), ly, 1)
        g = float(-cff[0])
        r2 = float(1 - np.var(ly - np.polyval(cff, np.log(dd))) / np.var(ly))
        out["runs"][tag] = dict(p=p, alpha=al, gamma=g, r2=r2)
        print(f"  [{tag}] p={p:.2f}  gamma={g:+.4f}  R2={r2:.4f}", flush=True)

    run("full", 1.0)
    for p in a.ps:
        run(f"p{p:g}", p)

    g0 = out["runs"]["full"]["gamma"]
    out["summary"] = {k: dict(p=v["p"], gamma=v["gamma"], delta=float(v["gamma"] - g0))
                      for k, v in out["runs"].items()}
    json.dump(out, open(od / "depth_control.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(f"\n写出 {od/'depth_control.json'}")
    print("判据：若各 p 的 gamma 与 full 之差 < 0.03，则深度不能解释组织间差异。")

if __name__ == "__main__":
    main()

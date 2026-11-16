r"""
32_synth_expression.py — 关键对照：gamma 随表达量变化，是真实结构还是低计数偏差？

问题
----
真实数据里，低表达层的 gamma 高（~1.0-1.3），高表达层低（~0.3-0.4）。
两种解释：
  H1 真实：低表达基因的空间结构更细碎（细胞特异），高表达基因更平滑（管家）。
  H2 偏差：低计数下 alpha 的矩估计有偏，gamma 被系统性抬高。
若 H2 成立，则"gamma 依赖表达量"是伪影，pooled gamma 才是真信号。

设计
----
生成一个**已知空间结构**（幂律相关场，相关指数 beta 固定）的合成组织，
让全部基因共享同一个 beta（即同一个真实 gamma），只让表达量跨 4 个数量级变化。
若估计无偏，各表达层的 gamma_hat 应相同；若随表达量变化，则 H2 成立。

场的生成
--------
二维功率谱 S(k) ∝ k^(beta-2)（对应 C(h) ∝ h^-beta），
低频截断在 k_min = 2*pi/L（L = 视窗边长），等效最大相关尺度。
用 FFT 合成，再 exp 转成强度 lambda(x)（对数正态，避免负值）。

计数
----
NB: Var = mu + alpha*mu^2，alpha 固定（模拟泊松之上的固定过散）。
聚合与主分析完全一致（矩形分块求和，零格保留）。
"""
import argparse, json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

def make_field(n, beta, rng, kmin_frac=1.0):
    """二维幂律相关场，标准化到 mean=0, var=1（对数域）。"""
    ky = np.fft.fftfreq(n)[:, None]
    kx = np.fft.fftfreq(n)[None, :]
    k = np.sqrt(kx ** 2 + ky ** 2)
    k[0, 0] = np.inf
    S = k ** ((beta - 2) / 2.0)                  # 幅度谱
    S[k < kmin_frac / n] = 0.0                   # 低频截断 = 最大相关尺度
    S[0, 0] = 0.0
    w = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    f = np.fft.ifft2(w * S).real
    f = (f - f.mean()) / (f.std() + 1e-12)
    return f

def nb_counts(lam, alpha, rng):
    """alpha=0 时纯 Poisson；否则 NB（Var=lam+alpha*lam^2），用 Gamma-Poisson 混合。
    注意：alpha>0 的过散是**每个 bin 独立**的，聚合后按 d^-2 衰减，
    会掩盖强度场的空间结构。检验空间标度时必须用 alpha=0。"""
    if alpha <= 0:
        return rng.poisson(lam).astype(np.float32)
    shape = 1.0 / alpha
    g = rng.gamma(shape, scale=lam * alpha, size=lam.shape)   # mean = lam
    return rng.poisson(g).astype(np.float32)

def pooled_alpha_stack(counts, ds):
    """counts: (n,n) float。返回 {d: alpha_pool}，口径与主分析一致（全域、零保留）。"""
    out = {}
    n = counts.shape[0]
    for d in ds:
        f = d // 2
        if n % f:
            continue
        b = counts.reshape(n // f, f, n // f, f).sum(axis=(1, 3))
        m = b.mean()
        v = b.var(ddof=1)
        out[d] = float((v - m) / (m ** 2)) if m > 0 else np.nan
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=str(ROOT / "02_OUT" / "synth_expr"))
    ap.add_argument("--n", type=int, default=1024)
    ap.add_argument("--ngenes", type=int, default=240)
    ap.add_argument("--beta", type=float, default=0.55)
    ap.add_argument("--alpha", type=float, default=1.0)
    ap.add_argument("--cv", type=float, default=0.3, help="强度场的变异系数（越小越接近线性，gamma 越接近 beta）")
    ap.add_argument("--transform", default="linear", choices=["linear", "log"])
    ap.add_argument("--seed", type=int, default=2026)
    a = ap.parse_args()
    od = Path(a.outdir); od.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    ds = [2, 4, 8, 16, 32, 64, 128, 256]

    n = a.n
    # 基因基线强度：log-normal 跨 4 个数量级（2um bin 均值 0.01 -> 100 UMI）
    logmu = rng.uniform(np.log(0.01), np.log(100.0), a.ngenes)
    base = np.exp(logmu)

    layers = [(0.01, 0.03), (0.03, 0.1), (0.1, 0.3), (0.3, 1.0), (1.0, 3.0), (3.0, 10.0), (10.0, 100.0)]
    per_layer = {f"{lo:g}-{hi:g}": [] for lo, hi in layers}
    all_gamma = []
    detail = []

    for gi in range(a.ngenes):
        f = make_field(n, a.beta, rng)
        if a.transform == "log":
            lam = base[gi] * np.exp(0.8 * f)
        else:
            lam = base[gi] * np.maximum(1.0 + a.cv * f, 0.05)
        lam = lam * (base[gi] / lam.mean())       # 精确校准均值
        c = nb_counts(lam, a.alpha, rng)
        al = pooled_alpha_stack(c, ds)
        dd = np.array([d for d in ds if d in al], float)
        aa = np.array([al[d] for d in ds if d in al], float)
        if np.any(aa <= 0) or not np.all(np.isfinite(aa)):
            continue
        g = float(-np.polyfit(np.log(dd), np.log(aa), 1)[0])
        mu2 = float(c.mean())
        all_gamma.append(g)
        detail.append(dict(i=gi, mu2=mu2, gamma=g))
        for lo, hi in layers:
            if lo <= mu2 < hi:
                per_layer[f"{lo:g}-{hi:g}"].append(g)
                break

    res = dict(beta_true=a.beta, alpha_true=a.alpha, n_grid=n, n_genes=a.ngenes,
               gamma_all=dict(med=float(np.median(all_gamma)),
                              q25=float(np.percentile(all_gamma, 25)),
                              q75=float(np.percentile(all_gamma, 75)),
                              sd=float(np.std(all_gamma))),
               by_layer={k: dict(n=len(v),
                                 med=float(np.median(v)) if v else None,
                                 q25=float(np.percentile(v, 25)) if v else None,
                                 q75=float(np.percentile(v, 75)) if v else None)
                         for k, v in per_layer.items()},
               detail=detail)
    json.dump(res, open(od / "synth_expression.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)

    print(f"真值 beta={a.beta} => 期望 gamma≈{a.beta:.2f}（C(h)∝h^-beta、mu∝d^2 时 alpha∝d^-beta）")
    print(f"全部 {len(all_gamma)} 个基因：gamma 中位 {res['gamma_all']['med']:.3f} "
          f"IQR [{res['gamma_all']['q25']:.3f},{res['gamma_all']['q75']:.3f}] sd={res['gamma_all']['sd']:.3f}")
    print("按 2um 每 bin 均值分层：")
    for k, v in res["by_layer"].items():
        if v["n"]:
            print(f"  mu∈{k:>10s}  n={v['n']:3d}  gamma 中位 {v['med']:.3f}  "
                  f"IQR [{v['q25']:.3f},{v['q75']:.3f}]")
    # 关键判据：最高层与最低层的中位之差
    ks = [k for k, v in res["by_layer"].items() if v["n"]]
    if len(ks) >= 2:
        lo, hi = res["by_layer"][ks[0]]["med"], res["by_layer"][ks[-1]]["med"]
        print(f"\n最低层 {ks[0]} gamma={lo:.3f}  最高层 {ks[-1]} gamma={hi:.3f}  差={lo-hi:+.3f}")
        print("若差接近 0，说明无偏，真实数据的'表达量依赖'是结构差异；" \
              "若差明显为正，说明低计数把 gamma 抬高（偏差）。")

if __name__ == "__main__":
    main()

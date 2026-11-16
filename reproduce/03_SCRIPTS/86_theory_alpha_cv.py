"""86_theory_alpha_cv.py — T2b 数值验证：alpha(A) = CV^2 * rhobar(A)

目的
----
把施工总纲 §3 的解析框架放到已知真值的合成场上检验。这是全文新理论的地基。

理论（连续，平稳，Poisson 复合）
--------------------------------
设表达强度场 X(s) 平稳，E[X]=m，Var[X]=sigma^2，相关函数 rho(h)。
Bin B 面积 A，条件于 X 的计数 N(B) ~ Poisson(int_B X)。

    E[N]   = m*A
    Var[N] = m*A + sigma^2 * A^2 * rhobar(A)
    rhobar(A) = (1/A^2) * int_B int_B rho(|s-t|) ds dt

    alpha(A) = (Var[N] - E[N]) / E[N]^2 = (sigma^2/m^2) * rhobar(A) = CV^2 * rhobar(A)

两份性质：
  1. Poisson 项在 (Var - E) 里精确约掉 => alpha 与测序深度无关（代数必然）。
  2. 归一化 alpha 曲线 alpha(d)/alpha(2um) 就是 bin 内平均相关函数。

标度（2D）
----------
  rho(h) ~ h^-beta  =>  rhobar(A) ~ A^(-beta/2) = d^(-beta)  =>  gamma = beta
  白噪声（点状）    =>  rhobar(A) ~ A^-1 = d^-2            =>  gamma = 2
  指数（有 xi）      =>  小 d 平台，大 d 转 d^-2，无单一 gamma

强度场用对数正态（严格正、无截断）
----------------------------------
  X = m * exp(sigma*F - sigma^2/2)，F 零均值单位方差高斯
  => E[X]=m, CV^2 = exp(sigma^2)-1, 且 log X 的相关结构 = F 的相关结构

在离散网格上精确算 rhobar
-------------------------
  bin = k x k 格。lag (a,b) 在 bin 内相距该 lag 的点对数 = (k-|a|)(k-|b|)
  rhobar_k = sum_{|a|,|b|<k} (k-|a|)(k-|b|) C(a,b) / k^4

用法
----
  python 86_theory_alpha_cv.py --outdir 02_OUT/theory --n 1024 --n-real 8
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

BIN_SIZES = [1, 2, 4, 8, 16, 32, 64]
FIT_KMIN, FIT_KMAX = 2, 64


# ---------------------------------------------------------------- 场生成
def _freq_grid(n: int):
    f = np.fft.fftfreq(n)
    KY, KX = np.meshgrid(f, f, indexing="ij")
    return np.sqrt(KY ** 2 + KX ** 2)


def make_gauss_field(n: int, kind: str, param: float, seed: int, n_real: int = 1):
    """零均值单位方差高斯场（周期边界）。"""
    rng = np.random.default_rng(seed)
    K = _freq_grid(n)
    K[0, 0] = np.inf
    if kind == "powerlaw":
        beta = float(param)
        if not (0 < beta < 2):
            raise ValueError("powerlaw 要求 0 < beta < 2")
        amp = K ** ((beta - 2.0) / 2.0)
    elif kind == "exponential":
        xi = float(param)
        amp = (1.0 + (2.0 * np.pi * xi * K) ** 2) ** (-0.75)
    elif kind == "white":
        amp = np.ones_like(K)
    else:
        raise ValueError(kind)
    amp[0, 0] = 0.0

    out = np.empty((n_real, n, n), dtype=np.float64)
    for i in range(n_real):
        W = np.fft.fft2(rng.normal(size=(n, n)))
        f = np.real(np.fft.ifft2(W * amp))
        f -= f.mean()
        s = f.std()
        if s > 0:
            f /= s
        out[i] = f
    return out


def to_lognormal(F, sigma):
    """零均值单位方差高斯 -> 均值1、CV^2=exp(sigma^2)-1 的对数正态场。"""
    return np.exp(sigma * F - 0.5 * sigma ** 2)


# ---------------------------------------------------------------- 协方差
def acov_periodic(fields):
    """周期边界下的自协方差（跨 realization 平均，已除以 n^2）。

    必须先中心化：否则 lag-0 给的是 E[X^2] 而不是 Var[X]，均值项会把 rho(h)
    整体撑平（rho(h) -> Cov(h)/(Var+mean^2)），rhobar 随之失真。
    """
    X = fields - fields.mean(axis=(-2, -1), keepdims=True)
    Fk = np.fft.fft2(X, axes=(-2, -1))
    return np.real(np.fft.ifft2(np.abs(Fk) ** 2, axes=(-2, -1))).mean(axis=0) / fields.shape[-1] ** 2


def rhobar_from_cov(acov_norm, k):
    """用精确 lag 权重从归一化协方差算 bin 内平均相关 rhobar_k。"""
    n = acov_norm.shape[0]
    tot = 0.0
    for a in range(-k + 1, k):
        wa = k - abs(a)
        for b in range(-k + 1, k):
            tot += wa * (k - abs(b)) * acov_norm[a % n, b % n]
    return tot / (k ** 4)


def rho_ring(acov_norm, max_lag):
    """按欧氏距离环平均得到 rho(h)，h=0..max_lag（周期边界）。"""
    n = acov_norm.shape[0]
    rho = np.zeros(max_lag + 1)
    cnt = np.zeros(max_lag + 1)
    for a in range(0, max_lag + 1):
        for b in range(-max_lag, max_lag + 1):
            h = int(round(np.hypot(a, b)))
            if h > max_lag:
                continue
            rho[h] += acov_norm[a % n, b % n]
            cnt[h] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(cnt > 0, rho / np.maximum(cnt, 1), np.nan)


# ---------------------------------------------------------------- 工具
def alpha_curve(field, ds, lam_scale, rng):
    """把强度场聚合成 k x k bin、抽 Poisson，返回 mu/var/alpha。"""
    n = field.shape[0]
    out = {}
    for k in ds:
        nb = n // k
        if nb < 8:
            continue
        B = field[: nb * k, : nb * k].reshape(nb, k, nb, k).sum(axis=(1, 3))
        cnt = rng.poisson(lam_scale * B).astype(np.float64)
        mu = cnt.mean()
        var = cnt.var(ddof=1)
        if mu > 0:
            out[k] = dict(mu=float(mu), var=float(var),
                          alpha=float((var - mu) / mu ** 2))
    return out


def fit_slope(xs, ys, kmin=None, kmax=None):
    xs = np.asarray(xs, float)
    ys = np.asarray(ys, float)
    m = np.isfinite(ys) & (ys > 0)
    if kmin is not None:
        m &= xs >= kmin
    if kmax is not None:
        m &= xs <= kmax
    if m.sum() < 3:
        return np.nan, np.nan, np.nan
    x, y = np.log(xs[m]), np.log(ys[m])
    A = np.vstack([x, np.ones_like(x)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ coef
    ss_res = ((y - pred) ** 2).sum()
    ss_tot = ((y - y.mean()) ** 2).sum()
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    se = np.sqrt(ss_res / max(len(x) - 2, 1) / ((x - x.mean()) ** 2).sum())
    return float(coef[0]), float(se), float(r2)


def run_case(kind, param, sigma, n, n_real, lam_scale, seed):
    """跑一个相关结构案例，返回 alpha 曲线 + 协方差诊断。"""
    F = make_gauss_field(n, kind, param, seed, n_real=n_real)
    X = np.stack([to_lognormal(F[i], sigma) for i in range(n_real)])
    cv2_theory = float(np.exp(sigma ** 2) - 1)
    ac = acov_periodic(X)
    # 样本 CV^2：恒等式里的 CV^2 必须与 rho 同源于这一批场，否则重尾分布的
    # 样本方差收敛慢会带来常数比例偏差（sigma=2 时约 -16%）
    cv2_sample = float(ac[0, 0] / (X.mean() ** 2))
    ac_n = ac / ac[0, 0]
    acc = {}
    for i in range(n_real):
        c = alpha_curve(X[i], BIN_SIZES, lam_scale, np.random.default_rng(seed + 991 + i))
        for k, v in c.items():
            acc.setdefault(k, []).append(v["alpha"])
    ks = sorted(acc)
    alphas = [float(np.nanmean(acc[k])) for k in ks]
    rows = []
    for k, al in zip(ks, alphas):
        rb = rhobar_from_cov(ac_n, k)
        pred = cv2_sample * rb
        rows.append(dict(k=k, alpha=al, rhobar=float(rb), predicted=float(pred),
                         rel_err=float((al - pred) / pred)))
    rho = rho_ring(ac_n, max_lag=min(64, n // 4))
    beta_emp, beta_se, beta_r2 = fit_slope(np.arange(1, 33), rho[1:33])
    return dict(rows=rows, cv2=cv2_sample, cv2_theory=cv2_theory,
                rho_h=rho[:65].tolist(),
                beta_emp=-beta_emp if np.isfinite(beta_emp) else np.nan,
                beta_se=beta_se, beta_r2=beta_r2, ks=ks, alphas=alphas)


# ---------------------------------------------------------------- 主流程
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="02_OUT/theory")
    ap.add_argument("--n", type=int, default=1024)
    ap.add_argument("--n-real", type=int, default=8)
    ap.add_argument("--sigma", type=float, default=0.8,
                    help="对数正态 sigma；CV^2 = exp(sigma^2)-1")
    ap.add_argument("--mean-count", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=20261115)
    a = ap.parse_args()

    od = Path(a.outdir)
    od.mkdir(parents=True, exist_ok=True)
    res = {"config": vars(a), "experiments": {}}

    # ---------- E1 核心恒等式 ----------
    print("=" * 78)
    print(f"E1 核心恒等式  alpha(A) = CV^2 * rhobar(A)   [lognormal, sigma={a.sigma}, "
          f"CV²={np.exp(a.sigma**2)-1:.3f}]")
    print("=" * 78)
    print(f"{'beta':>5} {'gamma_fit':>11} {'±se':>7} {'r2':>7} "
          f"{'beta_emp':>9} {'max|RE|':>9}")
    e1 = {}
    for beta in [0.3, 0.45, 0.55, 0.7, 0.9, 1.2, 1.5]:
        r = run_case("powerlaw", beta, a.sigma, a.n, a.n_real, a.mean_count,
                     a.seed + int(beta * 1000))
        g, gse, gr2 = fit_slope(r["ks"], r["alphas"], FIT_KMIN, FIT_KMAX)
        max_re = max(abs(x["rel_err"]) for x in r["rows"])
        e1[f"beta_{beta}"] = dict(true_beta=beta, gamma_fit=g, gamma_se=gse,
                                  gamma_r2=gr2, beta_emp=r["beta_emp"],
                                  beta_emp_r2=r["beta_r2"], rows=r["rows"],
                                  max_abs_rel_err=float(max_re))
        print(f"{beta:>5.2f} {g:>11.3f} {gse:>7.3f} {gr2:>7.4f} "
              f"{r['beta_emp']:>9.3f} {max_re:>9.2%}")
    res["experiments"]["E1_identity"] = e1

    # ---------- E2 三种相关结构 ----------
    print()
    print("=" * 78)
    print("E2 三种相关结构：标度指数")
    print("=" * 78)
    print(f"{'case':<18} {'gamma(2-64)':>12} {'r2':>7} {'gamma(16-64)':>13} {'r2':>7}")
    e2 = {}
    cases = [("white", 0.0), ("exponential", 2.0), ("exponential", 8.0),
             ("exponential", 32.0), ("exponential", 128.0),
             ("powerlaw", 0.55), ("powerlaw", 1.2)]
    for kind, param in cases:
        r = run_case(kind, param, a.sigma, a.n, a.n_real, a.mean_count,
                     a.seed + 7 + int(param * 10))
        g1, _, r21 = fit_slope(r["ks"], r["alphas"], 2, 64)
        g2, _, r22 = fit_slope(r["ks"], r["alphas"], 16, 64)
        tag = f"{kind}_{param:g}"
        e2[tag] = dict(ks=r["ks"], alphas=r["alphas"], gamma_2_64=g1, r2_2_64=r21,
                       gamma_16_64=g2, r2_16_64=r22, beta_emp=r["beta_emp"])
        print(f"{tag:<18} {g1:>12.3f} {r21:>7.4f} {g2:>13.3f} {r22:>7.4f}")
    res["experiments"]["E2_families"] = e2

    # ---------- E3 深度不变性 ----------
    print()
    print("=" * 78)
    print("E3 深度不变性（thinning 1000 倍）")
    print("=" * 78)
    e3 = {}
    F = make_gauss_field(a.n, "powerlaw", 0.55, a.seed + 55, n_real=a.n_real)
    X = np.stack([to_lognormal(F[i], a.sigma) for i in range(a.n_real)])
    for scale in [30.0, 10.0, 3.0, 1.0, 0.3, 0.1, 0.03]:
        acc = {}
        for i in range(a.n_real):
            c = alpha_curve(X[i], [2, 8, 32, 64], scale,
                            np.random.default_rng(a.seed + 200 + i))
            for k, v in c.items():
                acc.setdefault(k, []).append(v["alpha"])
        ks = sorted(acc)
        al = [float(np.nanmean(acc[k])) for k in ks]
        g, se, r2 = fit_slope(ks, al)
        e3[f"{scale:g}"] = dict(ks=ks, alphas=al, gamma=g, gamma_se=se, r2=r2,
                                mean_count_at_2um=scale * float(X[:, ::2, ::2].mean()))
        print(f"  ×{scale:<6g} γ={g:+.4f}±{se:.4f}  α(2)={al[0]:.5f}  α(64)={al[-1]:.5f}")
    gs = [v["gamma"] for v in e3.values()]
    res["experiments"]["E3_depth"] = dict(rows=e3, gamma_range=float(max(gs) - min(gs)),
                                          gamma_sd=float(np.std(gs)))
    print(f"  => γ 跨 1000 倍深度变幅 = {max(gs)-min(gs):.2e}")

    # ---------- E4 归一化塌缩 ----------
    print()
    print("=" * 78)
    print("E4 归一化 alpha 曲线 = 平均相关函数")
    print("=" * 78)
    e4 = {}
    for beta in [0.3, 0.55, 0.9]:
        r = run_case("powerlaw", beta, a.sigma, a.n, a.n_real, a.mean_count,
                     a.seed + 900 + int(beta * 10))
        al = np.array(r["alphas"])
        norm = (al / al[0]).tolist()
        e4[f"beta_{beta}"] = dict(ks=r["ks"], alpha_norm=norm)
        print(f"  beta={beta:.2f}  " + "  ".join(f"{k}:{v:.3f}" for k, v in zip(r["ks"], norm)))
    res["experiments"]["E4_collapse"] = e4

    # ---------- E5 CV 扫描 ----------
    print()
    print("=" * 78)
    print("E5 CV 扫描：恒等式是否随 CV 失效")
    print("=" * 78)
    e5 = {}
    for sigma in [0.3, 0.6, 0.8, 1.2, 1.6, 2.0]:
        r = run_case("powerlaw", 0.55, sigma, a.n, max(4, a.n_real // 2),
                     a.mean_count, a.seed + 400 + int(sigma * 10))
        max_re = max(abs(x["rel_err"]) for x in r["rows"])
        g, se, r2 = fit_slope(r["ks"], r["alphas"], FIT_KMIN, FIT_KMAX)
        e5[f"sigma_{sigma:g}"] = dict(cv2=r["cv2"], max_abs_rel_err=float(max_re),
                                       gamma=g, gamma_se=se, r2=r2, rows=r["rows"])
        print(f"  sigma={sigma:<4g} CV²={r['cv2']:9.3f}  max|RE|={max_re:6.2%}  γ={g:+.3f}")
    res["experiments"]["E5_cv_sweep"] = e5

    # ---------- 判定 ----------
    all_re = [v["max_abs_rel_err"] for v in e1.values()] + \
             [v["max_abs_rel_err"] for v in e5.values()]
    max_re = float(max(all_re))
    # 核心一致性检验：alpha 的拟合指数 gamma 是否等于场本身协方差的指数 beta_emp
    g = np.array([abs(v["gamma_fit"]) for v in e1.values()])
    b = np.array([v["beta_emp"] for v in e1.values()])
    diff = g - b
    slope, _, r2 = fit_slope(b, g)
    verdict = dict(
        max_abs_rel_err_overall=max_re,
        median_abs_rel_err=float(np.median(all_re)),
        threshold=0.05,
        pass_identity=bool(max_re < 0.05),
        depth_gamma_range=res["experiments"]["E3_depth"]["gamma_range"],
        gamma_vs_beta_emp_slope=slope,
        gamma_vs_beta_emp_r2=r2,
        gamma_minus_beta_emp_max_abs=float(np.max(np.abs(diff))),
        gamma_minus_beta_emp_mean=float(np.mean(diff)),
        gamma_vs_nominal_beta_max_dev=float(max(
            abs(abs(v["gamma_fit"]) - v["true_beta"]) for v in e1.values())),
        note=("gamma 对 k(线性 bin 尺寸)拟合，A=k^2，故斜率直接等于 beta。"
              "对照基准应是用同一批场实测的协方差指数 beta_emp，"
              "而不是名义参数 beta（离散化与 Nyquist 截断使 beta_emp 本身略高于名义值）。"),
    )
    res["verdict"] = verdict
    print()
    print("=" * 78)
    print(f"恒等式最大相对误差 = {max_re:.2%}  (阈值 5%)  "
          f"=> {'PASS' if verdict['pass_identity'] else 'FAIL'}")
    print(f"gamma vs beta_emp: 斜率={slope:.3f} (r²={r2:.5f})，"
          f"|差| 最大={verdict['gamma_minus_beta_emp_max_abs']:.3f} "
          f"平均={verdict['gamma_minus_beta_emp_mean']:+.3f}")
    print(f"gamma vs 名义 beta: 最大偏离={verdict['gamma_vs_nominal_beta_max_dev']:.3f} "
          f"（离散化所致，不代表理论失效）")
    print(f"深度不变性：γ 变幅 = {verdict['depth_gamma_range']:.2e}")
    print("=" * 78)

    (od / "theory_check.json").write_text(
        json.dumps(res, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"写出 {od/'theory_check.json'}")


if __name__ == "__main__":
    main()

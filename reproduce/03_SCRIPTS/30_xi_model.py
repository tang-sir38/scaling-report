r"""
30_xi_model.py — 有限相关长度模型：把 gamma 变成可估的组织参数。

原理
----
对基因 g，bin（边长 d 的方块 A）内计数 Y 的过散度
    alpha_g(d) = Var(Lambda_A) / mu_A^2 ,  Lambda_A = ∫_A lambda_g(x) dx
若 lambda_g 二阶平稳，协方差 sigma_g^2 * rho(h/xi)，则
    Var(Lambda_A) = sigma_g^2 * ∫∫_A rho(|x-y|/xi) dx dy = sigma_g^2 * d^4 * H(d/xi)
    mu_A = lambda_bar_g * d^2
于是
    alpha_g(d) = (sigma_g^2/lambda_bar_g^2) * H(d/xi) = c_g * H(d/xi)
H 是无量纲形状函数：
    H(w) = ∫_{-1}^{1}∫_{-1}^{1} (1-|s|)(1-|t|) rho(w*sqrt(s^2+t^2)) ds dt
极限：w→0  H→1（全相关，gamma=0）；w→∞  H→2*pi/w^2（独立，gamma=2）。
pooled alpha 因为权重 mu^2 与 d 无关，形状仍是 H 的混合：
    alpha_pool(d) = sum_g w_g c_g H(d/xi_g) / sum_g w_g ,  w_g = lambda_bar_g^2
若 xi_g 与 g 无关，alpha_pool(d) ∝ H(d/xi)。

模型
----
M1 纯幂律    alpha = A d^-gamma                     (2 参数)
M2 指数相关  alpha = A H_exp(d/xi)                  (2 参数)
M3 高斯相关  alpha = A H_gauss(d/xi)                (2 参数)
M4 指数+幂律尾（Matérn 型，nu 固定 1.5）             (3 参数)
用 AIC 比较；xi 用 profile likelihood 给区间（F 检验阈值）。

输出
----
  <outdir>/xi_model.json
  <outdir>/xi_model.png
"""
import argparse, json, sys
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.interpolate import CubicSpline

ROOT = Path(__file__).resolve().parents[1]

# ---------- 形状函数 H(w) ----------
def _gl_nodes(n=240):
    x, w = np.polynomial.legendre.leggauss(n)
    return x, w

def build_H(kind="exp", n=240, wmax=400.0, nw=200):
    """在 log 网格上预计算 H(w)，返回可插值的 CubicSpline（外推用渐近式）。"""
    x, wq = _gl_nodes(n)
    s = x[:, None]                      # (n,1)
    t = x[None, :]                      # (1,n)
    wt = (wq[:, None] * wq[None, :])    # (n,n)
    taper = (1.0 - np.abs(s)) * (1.0 - np.abs(t))
    r = np.sqrt(s ** 2 + t ** 2)
    wgrid = np.concatenate([[0.0], np.geomspace(1e-3, wmax, nw)])
    H = np.empty_like(wgrid)
    for i, w in enumerate(wgrid):
        if kind == "exp":
            k = np.exp(-w * r)
        elif kind == "gauss":
            k = np.exp(-0.5 * (w * r) ** 2)
        else:
            raise ValueError(kind)
        H[i] = float(np.sum(wt * taper * k))
    # 归一化：H(0) 应精确为 1（∫∫(1-|s|)(1-|t|) = 1）
    H = H / H[0]
    spl = CubicSpline(np.log(wgrid[1:]), np.log(H[1:]), extrapolate=True)
    return spl, wgrid, H

def H_eval(spl, w):
    w = np.asarray(w, float)
    out = np.empty_like(w)
    small = w < 1e-3
    out[small] = 0.0                      # log H → 0
    big = ~small
    if np.any(big):
        out[big] = spl(np.log(w[big]))
    return np.exp(out)

def dlogH_dlogw(spl, w, h=1e-4):
    """d ln H / d ln w，数值中心差分。"""
    return (spl(np.log(w) + h) - spl(np.log(w) - h)) / (2 * h)

# ---------- 拟合 ----------
def fit_quad(d, a):
    """log alpha = a0 + a1 log d + a2 (log d)^2。a2>0 表示斜率随尺度上升（向上弯）。"""
    from scipy.stats import f as fdist
    x, y = np.log(d), np.log(a)
    c = np.polyfit(x, y, 2)
    pred = np.polyval(c, x)
    rss = float(np.sum((y - pred) ** 2))
    lin = np.polyfit(x, y, 1)
    rss_lin = float(np.sum((y - np.polyval(lin, x)) ** 2))
    n = len(d)
    # 嵌套 F 检验：二次项是否显著
    F = ((rss_lin - rss) / 1) / (rss / (n - 3)) if rss > 0 else np.inf
    p = float(fdist.sf(F, 1, n - 3))
    return dict(model="quad", k=3,
                params=dict(a0=float(c[2]), a1=float(c[1]), a2=float(c[0])),
                rss=rss, rss_lin=rss_lin, F=float(F), p_F=p,
                p_lt_0p05=bool(p < 0.05),
                vertex_logd=float(-c[1] / (2 * c[0])) if c[0] != 0 else None,
                vertex_d=float(np.exp(-c[1] / (2 * c[0]))) if c[0] != 0 else None)


def fit_powerlaw(d, a):
    x, y = np.log(d), np.log(a)
    p = np.polyfit(x, y, 1)
    pred = np.polyval(p, x)
    rss = float(np.sum((y - pred) ** 2))
    return dict(model="power", k=2, params=dict(A=float(np.exp(p[1])), gamma=float(-p[0])),
                rss=rss, resid=float(np.max(np.abs(y - pred))))

def fit_xi(d, a, spl, kind):
    y = np.log(a)
    def resid(p):
        logA, logxi = p
        return y - (logA + np.log(H_eval(spl, d / np.exp(logxi))))
    best = None
    for lx0 in (np.log(20), np.log(200), np.log(2000)):
        r = least_squares(resid, [np.log(a[0]), lx0], method="lm", max_nfev=4000)
        if best is None or r.cost < best.cost:
            best = r
    rss = float(np.sum(best.fun ** 2))
    return dict(model=kind, k=2,
                params=dict(A=float(np.exp(best.x[0])), xi=float(np.exp(best.x[1]))),
                rss=rss, resid=float(np.max(np.abs(best.fun))))

def profile_xi(d, a, spl, kind, n=400, lo=1.0, hi=1e6):
    """固定 xi，解析最优 A = exp(mean(y - logH))；返回 (xi grid, rss, best)。"""
    y = np.log(a)
    xs = np.geomspace(lo, hi, n)
    rss = np.empty(n)
    for i, xi in enumerate(xs):
        lh = np.log(H_eval(spl, d / xi))
        A = np.mean(y - lh)
        rss[i] = float(np.sum((y - lh - A) ** 2))
    j = int(np.argmin(rss))
    return xs, rss, float(xs[j]), float(rss[j])

def aic(rss, n, k):
    return float(n * np.log(rss / n) + 2 * k)

# ---------- 主流程 ----------
def analyse(tag, pj, outdir, nboot=400, seed=7, max_d=None):
    j = json.load(open(pj))
    ds = [d for d in j["ds"] if (max_d is None or d <= max_d)]
    alpha = np.array([j["per_d"][str(d)]["alpha_pool"] for d in ds], float)
    ncell = np.array([j["per_d"][str(d)]["n_cells"] for d in ds], float)
    d = np.array(ds, float)
    # pooled alpha 的近似 SE：由总计数场的过散度与有效 bin 数给出
    # 采用与 quadrant 一致的相对误差量级（保守 3%），仅用于加权展示，不用于 CI
    res = dict(tag=tag, ds=ds, alpha=alpha.tolist(), n_cells=ncell.tolist(), models={})

    spl_e, wg, Hg = build_H("exp")
    spl_g, _, _ = build_H("gauss")
    # 自检：H(0)=1；H(w) 在 w 大时应趋近 2*pi/w^2
    chk = [(w, float(Hg[np.argmin(np.abs(wg - w))]), 2 * np.pi / w ** 2) for w in (10, 30, 100)]
    res["H_selfcheck"] = [dict(w=w, H=float(h), asymptote=float(a_)) for w, h, a_ in chk]

    res["models"]["power"] = fit_powerlaw(d, alpha)
    res["models"]["quad"] = fit_quad(d, alpha)
    res["models"]["exp"] = fit_xi(d, alpha, spl_e, "exp")
    res["models"]["gauss"] = fit_xi(d, alpha, spl_g, "gauss")

    n = len(d)
    for k, m in res["models"].items():
        m["aic"] = aic(m["rss"], n, m["k"])
        m.pop("resid", None)
    best_aic = min(m["aic"] for m in res["models"].values())
    for k, m in res["models"].items():
        m["dAIC"] = float(m["aic"] - best_aic)

    # profile likelihood for xi (指数模型)
    xs, rss, xib, rssb = profile_xi(d, alpha, spl_e, "exp")
    thr = rssb * (1.0 + 5.99 / (n - 2))          # F(1, n-2, 0.95) = 5.99
    ok = xs[rss <= thr]
    res["xi_profile"] = dict(
        xi_grid=xs.tolist(), rss=rss.tolist(), xi_hat=xib,
        ci_lo=float(ok.min()) if ok.size else None,
        ci_hi=float(ok.max()) if ok.size else None,
        threshold=float(thr), rss_min=rssb,
        note=("ci_hi=hi 表示上界未闭合（数据只给出下界）" if ok.size and ok.max() >= xs[-1] * 0.99 else ""))

    # 局部斜率：实测 vs 二次模型（ξ 模型的预测单独存，便于显示它为何被拒）
    xi = res["models"]["exp"]["params"]["xi"]
    w = d / xi
    pred_local_xi = -dlogH_dlogw(spl_e, w)
    obs_local = np.array([j["segments"][f"{ds[i]}-{ds[i+1]}"] if isinstance(j["segments"], dict)
                          else j["segments"][i]["slope"] for i in range(len(ds) - 1)], float)
    seg_mid = np.sqrt(d[:-1] * d[1:])
    pred_seg_xi = -dlogH_dlogw(spl_e, seg_mid / xi)
    c = res["models"]["quad"]["params"]
    pred_seg_quad = -(c["a1"] + 2 * c["a2"] * np.log(seg_mid))
    res["local_slope"] = dict(
        seg_mid=seg_mid.tolist(), observed=(-obs_local).tolist(),
        predicted_xi=pred_seg_xi.tolist(), predicted_quad=pred_seg_quad.tolist(),
        rmse_xi=float(np.sqrt(np.mean(((-obs_local) - pred_seg_xi) ** 2))),
        rmse_quad=float(np.sqrt(np.mean(((-obs_local) - pred_seg_quad) ** 2))))

    # 外推：达到 gamma=1.5 所需的 bin
    w_target = None
    ww = np.geomspace(1e-2, 1e4, 20000)
    ss = -dlogH_dlogw(spl_e, ww)
    idx = np.where(ss >= 1.5)[0]
    if idx.size:
        w_target = float(ww[idx[0]])
    res["extrapolation"] = dict(
        w_at_slope_1p5=w_target,
        d_at_slope_1p5=(w_target * xi) if w_target else None,
        note="模型外推，超出观测窗口 256 um，仅作方向性预测")

    # bootstrap（对 8 个点做残差重抽样，给 gamma 与 xi 的经验区间）
    rng = np.random.default_rng(seed)
    y = np.log(alpha)
    fit_pred = np.log(res["models"]["exp"]["params"]["A"]) + np.log(
        H_eval(spl_e, d / xi))
    resid = y - fit_pred
    gs, xis = [], []
    for _ in range(nboot):
        yy = fit_pred + rng.choice(resid, size=n, replace=True)
        try:
            m = fit_xi(d, np.exp(yy), spl_e, "exp")
            xis.append(m["params"]["xi"])
            gs.append(-np.polyfit(np.log(d), yy, 1)[0])
        except Exception:
            pass
    if gs:
        res["boot"] = dict(
            gamma_med=float(np.median(gs)), gamma_lo=float(np.percentile(gs, 2.5)),
            gamma_hi=float(np.percentile(gs, 97.5)),
            xi_med=float(np.median(xis)), xi_lo=float(np.percentile(xis, 2.5)),
            xi_hi=float(np.percentile(xis, 97.5)), n=int(len(gs)))
    return res

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=str(ROOT / "02_OUT" / "xi_model"))
    ap.add_argument("--max-d", type=int, default=None,
                    help="只用 d <= max_d 的尺度（2-256 um 口径传 256）")
    ap.add_argument("--outname", default="xi_model.json")
    a = ap.parse_args()
    od = Path(a.outdir); od.mkdir(parents=True, exist_ok=True)
    sets = {
        "colon": str(ROOT / "02_OUT" / "pooled_alpha_noMT.json"),
        "lung": str(ROOT / "D3_LungCancer" / "02_OUT" / "pooled_alpha_noMT.json"),
        "brain": str(ROOT / "D2_MouseBrain" / "02_OUT" / "pooled_alpha_noMT.json"),
    }
    out = {}
    for tag, pj in sets.items():
        if not Path(pj).exists():
            print(f"缺 {pj}"); continue
        r = analyse(tag, pj, od, max_d=a.max_d)
        out[tag] = r
        print(f"\n== {tag} ==")
        for k, m in r["models"].items():
            extra = f" xi={m['params'].get('xi', float('nan')):.1f}" if "xi" in m["params"] else ""
            print(f"  {k:6s} rss={m['rss']:.5f} AIC={m['aic']:.2f} dAIC={m['dAIC']:.2f}{extra}")
        xp = r["xi_profile"]
        print(f"  xi_hat={xp['xi_hat']:.1f}  CI=[{xp['ci_lo']:.1f}, {xp['ci_hi']:.1f}]  {xp['note']}")
        ls = r["local_slope"]
        print(f"  局部斜率实测  {np.round(ls['observed'],3)}")
        print(f"  二次模型预测  {np.round(ls['predicted_quad'],3)}  RMSE={ls['rmse_quad']:.3f}")
        print(f"  xi模型预测    {np.round(ls['predicted_xi'],3)}  RMSE={ls['rmse_xi']:.3f}")
        if "boot" in r:
            b = r["boot"]
            print(f"  boot gamma [{b['gamma_lo']:.3f},{b['gamma_hi']:.3f}]  xi [{b['xi_lo']:.1f},{b['xi_hi']:.1f}]")
    json.dump(out, open(od / a.outname, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(f"\n写出 {od/'xi_model.json'}")

if __name__ == "__main__":
    main()

r"""87_variogram_shape.py — T3：用独立路径测空间协方差指数 beta_g，与 gamma_g 对照

为什么这是"独立"的
------------------
gamma_g（现稿已有）来自**聚合**：把 2 um 网格重采样到 d um，看 alpha_g(d) 的
衰减斜率。它只用到每个尺度上的 (mu_d, var_d) 两个数。

beta_g（本脚本新算）来自**空间**：在 8 um 网格上直接测同一基因的自协方差
rho_g(h) 随距离 h 的衰减。它用到场的逐点配对结构。

理论预测（86_theory_alpha_cv.py 已在合成场验证）：
    alpha(A) = CV^2 * rhobar(A)，而 rhobar(A) 是 rho(h) 的块内平均
    => alpha 的标度指数 gamma = rho 的标度指数 beta
所以 gamma_g 与 beta_g 应当一致。若一致，"alpha 指数 = 场的几何性质"成立，
gamma 就有了解释与可预测性；若不一致，说明聚合路径有额外偏差，需如实报告。

同时做两个配套检验
------------------
P1 塌缩：alpha_g(d)/alpha_g(8um) 是否在所有基因上塌缩为同一函数（= 相关函数）
P3 occupancy：gamma_g 与"检测到该基因的 bin 比例"的相关是否强于与丰度的相关

输出
----
02_OUT/variogram/beta_g_<tag>.npz   逐基因 beta、rho(h)、拟合质量
02_OUT/variogram/shape_<tag>.json   汇总、P1/P2/P3 检验结果

用法
----
  python 87_variogram_shape.py --tag colon --bins 01_WORK/bins_008um.npz \
      --moments 02_OUT/gene_moments.npz --curves 02_OUT/gene_curves/gene_curves_colon.npz
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

KEY_STRIDE = 4096
ROOT = Path(__file__).resolve().parents[1]

SETS = {
    "colon": dict(bins="01_WORK/bins_008um.npz",
                  moments="02_OUT/gene_moments.npz",
                  curves="02_OUT/gene_curves/gene_curves_colon.npz",
                  ds=8),
    "lung": dict(bins="D3_LungCancer/01_WORK/bins_008um.npz",
                 moments="D3_LungCancer/02_OUT/gene_moments_noMT.npz",
                 curves="02_OUT/gene_curves/gene_curves_lung.npz",
                 ds=8),
    "brain": dict(bins="D2_MouseBrain/01_WORK/bins_008um.npz",
                  moments="D2_MouseBrain/02_OUT/gene_moments_noMT.npz",
                  curves="02_OUT/gene_curves/gene_curves_brain.npz",
                  ds=8),
}


# ------------------------------------------------------------------ 场与协方差
def load_field_stack(npz_path: Path, gene_idx: np.ndarray):
    """返回 (F, M, genes)。F: (n_gene, H, W) float32；M: (H,W) bool 域掩码。"""
    z = np.load(npz_path, allow_pickle=True)
    shape = tuple(int(v) for v in z["shape"])
    X = csr_matrix((z["data"], z["indices"], z["indptr"]), shape=shape)
    coords = z["coords"].astype(np.int64)
    genes = np.array([str(g) for g in z["genes"]])

    # 网格范围必须取 coords 与 dom_any 的并集：零计数的域内 bin 只在 dom_any 里，
    # 只用 coords 会漏掉它们，导致 H/W 偏小、索引越界（lung 曾因此报错）。
    keys = z["dom_any"].astype(np.int64) if "dom_any" in z.files else None
    if keys is not None:
        kr, kc = keys // KEY_STRIDE, keys % KEY_STRIDE
        r0 = int(min(coords[:, 0].min(), kr.min()))
        c0 = int(min(coords[:, 1].min(), kc.min()))
        r1 = int(max(coords[:, 0].max(), kr.max()))
        c1 = int(max(coords[:, 1].max(), kc.max()))
    else:
        r0, c0 = int(coords[:, 0].min()), int(coords[:, 1].min())
        r1, c1 = int(coords[:, 0].max()), int(coords[:, 1].max())
    H, W = r1 - r0 + 1, c1 - c0 + 1

    M = np.zeros((H, W), dtype=bool)
    if keys is not None:
        M[kr - r0, kc - c0] = True
    else:
        M[coords[:, 0] - r0, coords[:, 1] - c0] = True

    pos_r = coords[:, 0] - r0
    pos_c = coords[:, 1] - c0

    sub = X[gene_idx, :]                     # 只取需要的基因
    F = np.zeros((len(gene_idx), H, W), dtype=np.float32)
    for j in range(sub.shape[0]):
        row = sub.getrow(j)
        rr = pos_r[row.indices]
        cc = pos_c[row.indices]
        F[j, rr, cc] = row.data
    return F, M, genes[gene_idx], (r0, c0, H, W)


def masked_acov(F, M):
    """掩码自协方差。返回 (cov, w) 两个 (H,W) 数组，cov 为按配对数归一化的协方差。"""
    Mf = M.astype(np.float64)
    A = F.astype(np.float64) * Mf
    FA = np.fft.rfft2(A)
    FM = np.fft.rfft2(Mf)
    num = np.fft.irfft2(FA * np.conj(FA), s=F.shape)
    den = np.fft.irfft2(FM * np.conj(FM), s=F.shape)
    mean = A.sum() / Mf.sum()
    with np.errstate(invalid="ignore", divide="ignore"):
        cov = np.where(den > 1.0, num / den - mean ** 2, np.nan)
    return cov, den


def masked_acov_cached(F, FM_conj, den, mask_sum):
    """同一个掩码下复用 fft(M)，省一次 FFT（M 固定时唯一正确的做法）。"""
    Mf = FM_conj[0]
    A = F.astype(np.float64)
    FA = np.fft.rfft2(A)
    num = np.fft.irfft2(FA * np.conj(FA), s=F.shape)
    mean = A.sum() / mask_sum
    with np.errstate(invalid="ignore", divide="ignore"):
        cov = np.where(den > 1.0, num / den - mean ** 2, np.nan)
    return cov


def make_ring_index(max_lag, W):
    """预计算环平均的扁平索引与 lag 归属，替掉 Python 双重循环。"""
    idx, hid = [], []
    for a in range(0, max_lag + 1):
        for b in range(-max_lag, max_lag + 1):
            h = int(round(np.hypot(a, b)))
            if h > max_lag:
                continue
            idx.append(a * W + (b % W))
            hid.append(h)
    return np.array(idx, dtype=np.int64), np.array(hid, dtype=np.int64)


def ring_average_fast(A, idx, hid, max_lag):
    flat = A.ravel()[idx]
    ok = np.isfinite(flat)
    if not ok.any():
        return np.full(max_lag + 1, np.nan)
    tot = np.bincount(hid[ok], weights=flat[ok], minlength=max_lag + 1)
    cnt = np.bincount(hid[ok], minlength=max_lag + 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(cnt > 0, tot / np.maximum(cnt, 1), np.nan)


def ring_average(A, max_lag):
    """把 2D lag 图按欧氏距离环平均（只用 lag>=0 的半平面）。"""
    out = np.zeros(max_lag + 1)
    cnt = np.zeros(max_lag + 1)
    for a in range(0, max_lag + 1):
        for b in range(-max_lag, max_lag + 1):
            h = int(round(np.hypot(a, b)))
            if h > max_lag:
                continue
            v = A[a, b % A.shape[1]]
            if np.isfinite(v):
                out[h] += v
                cnt[h] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(cnt > 0, out / np.maximum(cnt, 1), np.nan)


def fit_beta(rho, lags, kmin, kmax):
    """log-log 拟合 rho ~ h^-beta，只用 rho>0 的 lag。"""
    h = np.asarray(lags, float)
    r = np.asarray(rho, float)
    m = np.isfinite(r) & (r > 0) & (h >= kmin) & (h <= kmax)
    if m.sum() < 3:
        return np.nan, np.nan, np.nan
    x, y = np.log(h[m]), np.log(r[m])
    A = np.vstack([x, np.ones_like(x)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ coef
    ss_res = ((y - pred) ** 2).sum()
    ss_tot = ((y - y.mean()) ** 2).sum()
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return float(-coef[0]), float(coef[1]), float(r2)


# ------------------------------------------------------------------ 主流程
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="colon", choices=list(SETS))
    ap.add_argument("--outdir", default="02_OUT/variogram")
    ap.add_argument("--min-mu", type=float, default=0.05,
                    help="8um 下平均计数阈值（低于此则场太稀疏，rho 估不准）")
    ap.add_argument("--max-lag", type=int, default=32)
    ap.add_argument("--fit-kmin", type=int, default=2)
    ap.add_argument("--fit-kmax", type=int, default=24)
    ap.add_argument("--max-genes", type=int, default=0, help="0=全部合格基因")
    a = ap.parse_args()

    cfg = SETS[a.tag]
    od = Path(a.outdir)
    od.mkdir(parents=True, exist_ok=True)

    # --- 基因筛选（与 gene_curves 口径对齐 + 稀疏度门槛）---
    mo = np.load(ROOT / cfg["moments"], allow_pickle=True)
    ds_all = [int(x) for x in mo["ds"]]
    genes_mo = np.array([str(g) for g in mo["genes"]])
    dsel = cfg["ds"]
    mu_d = mo[f"mu_{dsel}"]
    gc = np.load(ROOT / cfg["curves"], allow_pickle=True)
    gc_genes = np.array([str(g) for g in gc["gene"]])
    gamma_all = gc["gamma"]
    fitted = gc["fitted"]

    # 对齐两个基因表
    pos = {g: i for i, g in enumerate(genes_mo)}
    order = np.array([pos[g] for g in gc_genes])
    mu8 = mu_d[order]

    keep = fitted & np.isfinite(gamma_all) & (mu8 >= a.min_mu)
    idx_curves = np.where(keep)[0]
    if a.max_genes and idx_curves.size > a.max_genes:
        rng = np.random.default_rng(20261115)
        idx_curves = np.sort(rng.choice(idx_curves, a.max_genes, replace=False))
    gene_names = gc_genes[idx_curves]
    gamma = gamma_all[idx_curves]
    print(f"[{a.tag}] 合格基因 {idx_curves.size} / {gc_genes.size} "
          f"(fitted={int(fitted.sum())}, mu_{dsel}>={a.min_mu})")

    # --- 建场并逐基因算协方差 ---
    idx_in_moments = order[idx_curves]
    F, M, genes_used, geom = load_field_stack(ROOT / cfg["bins"], idx_in_moments)
    r0, c0, H, W = geom
    print(f"  网格 {H}x{W}，域内 {int(M.sum()):,} 格，取 {F.shape[0]} 基因")

    lags = np.arange(0, a.max_lag + 1)
    beta = np.full(F.shape[0], np.nan)
    beta_icpt = np.full(F.shape[0], np.nan)
    beta_r2 = np.full(F.shape[0], np.nan)
    rho_mat = np.full((F.shape[0], a.max_lag + 1), np.nan)
    alpha_norm = np.full((F.shape[0], 8), np.nan)   # P1 用
    cv2 = np.full(F.shape[0], np.nan)
    occupancy = np.full(F.shape[0], np.nan)

    # P1 用的 alpha 曲线（直接从 moments 取，不重算）
    AL = np.full((F.shape[0], 8), np.nan)
    ds_arr = np.array(ds_all, float)
    for j, im in enumerate(idx_in_moments):
        mu_v = np.array([mo[f"mu_{d}"][im] for d in ds_all])
        va_v = np.array([mo[f"var_{d}"][im] for d in ds_all])
        with np.errstate(invalid="ignore", divide="ignore"):
            AL[j] = (va_v - mu_v) / np.where(mu_v > 0, mu_v ** 2, np.nan)

    # 预计算：掩码侧的 FFT 只需一次（M 固定）；环索引也只建一次
    Mf = M.astype(np.float64)
    FM = np.fft.rfft2(Mf)
    den = np.fft.irfft2(FM * np.conj(FM), s=M.shape)
    mask_sum = float(Mf.sum())
    if den.max() <= 0:
        raise SystemExit("掩码自相关全为 0，检查几何")
    ring_idx, ring_hid = make_ring_index(a.max_lag, W)
    print(f"  预计算完成：环条目 {ring_idx.size}，fft 复用已开启")

    for j in range(F.shape[0]):
        f = F[j]
        A = f.astype(np.float64)
        FA = np.fft.rfft2(A)
        num = np.fft.irfft2(FA * np.conj(FA), s=f.shape)
        mean = A.sum() / mask_sum
        with np.errstate(invalid="ignore", divide="ignore"):
            cov = np.where(den > 1.0, num / den - mean ** 2, np.nan)
        var0 = cov[0, 0]
        if not np.isfinite(var0) or var0 <= 0:
            continue
        rho = ring_average_fast(cov, ring_idx, ring_hid, a.max_lag) / var0
        rho_mat[j] = rho
        b, ic, r2 = fit_beta(rho, lags, a.fit_kmin, a.fit_kmax)
        beta[j], beta_icpt[j], beta_r2[j] = b, ic, r2
        fv = f[M]
        cv2[j] = var0 / (fv.mean() ** 2)
        occupancy[j] = float((fv > 0).mean())
        if j % 200 == 0:
            print(f"    {j}/{F.shape[0]}  beta={b:.3f} r2={r2:.3f}", flush=True)

    # P1: 归一化 alpha 曲线
    with np.errstate(invalid="ignore", divide="ignore"):
        alpha_norm = AL / AL[:, [int(np.where(ds_arr == dsel)[0][0])]]

    # --- 汇总 ---
    ok = np.isfinite(beta) & np.isfinite(gamma) & (beta_r2 > 0.5)
    print(f"\n  有效拟合 {int(ok.sum())} / {F.shape[0]}")

    def corr(x, y, m):
        x, y = x[m], y[m]
        if x.size < 5:
            return np.nan
        return float(np.corrcoef(x, y)[0, 1])

    from scipy.stats import spearmanr
    pear = corr(beta, gamma, ok)
    spear = float(spearmanr(beta[ok], gamma[ok]).statistic) if ok.sum() >= 5 else np.nan

    # 同质斜率：gamma 与 beta 的差
    diff = gamma[ok] - beta[ok]

    # P3: occupancy vs abundance
    logmu = np.log10(np.maximum(mu8[idx_curves], 1e-12))
    r_occ = corr(gamma, occupancy, ok)
    r_mu = corr(gamma, logmu, ok)
    # 偏相关 gamma ~ occupancy | logmu
    def partial(x, y, z, m):
        xs, ys, zs = x[m], y[m], z[m]
        if xs.size < 10:
            return np.nan
        bxy = np.polyfit(zs, xs, 1); byy = np.polyfit(zs, ys, 1)
        rx = xs - np.polyval(bxy, zs); ry = ys - np.polyval(byy, zs)
        return float(np.corrcoef(rx, ry)[0, 1])
    p_occ = partial(gamma, occupancy, logmu, ok)
    p_mu = partial(gamma, logmu, occupancy, ok)

    summary = dict(
        tag=a.tag, n_genes=int(F.shape[0]), n_valid=int(ok.sum()),
        grid=[H, W], domain_cells=int(M.sum()),
        min_mu=a.min_mu, fit_lag=[a.fit_kmin, a.fit_kmax],
        beta=dict(median=float(np.nanmedian(beta[ok])),
                  q25=float(np.percentile(beta[ok], 25)),
                  q75=float(np.percentile(beta[ok], 75)),
                  mean=float(np.nanmean(beta[ok]))),
        gamma=dict(median=float(np.nanmedian(gamma[ok])),
                   q25=float(np.percentile(gamma[ok], 25)),
                   q75=float(np.percentile(gamma[ok], 75))),
        gamma_minus_beta=dict(median=float(np.median(diff)),
                              mean=float(np.mean(diff)),
                              sd=float(np.std(diff))),
        P2_correlation=dict(pearson=pear, spearman=spear,
                            note="gamma(聚合路径) vs beta(空间协方差路径)"),
        P3_occupancy=dict(r_gamma_occupancy=r_occ, r_gamma_logmu=r_mu,
                          partial_gamma_occupancy_given_mu=p_occ,
                          partial_gamma_mu_given_occupancy=p_mu),
        P1_collapse=dict(
            n_genes=int(ok.sum()),
            alpha_norm_median=np.nanmedian(alpha_norm[ok], axis=0).tolist(),
            alpha_norm_q25=np.nanpercentile(alpha_norm[ok], 25, axis=0).tolist(),
            alpha_norm_q75=np.nanpercentile(alpha_norm[ok], 75, axis=0).tolist(),
            ds=ds_all,
        ),
    )

    print("\n" + "=" * 72)
    print(f"  beta(空间)  中位 {summary['beta']['median']:.3f} "
          f"[{summary['beta']['q25']:.3f}, {summary['beta']['q75']:.3f}]")
    print(f"  gamma(聚合) 中位 {summary['gamma']['median']:.3f} "
          f"[{summary['gamma']['q25']:.3f}, {summary['gamma']['q75']:.3f}]")
    print(f"  gamma - beta 中位 {summary['gamma_minus_beta']['median']:+.3f} "
          f"(SD {summary['gamma_minus_beta']['sd']:.3f})")
    print(f"  P2 相关: Pearson {pear:+.3f}  Spearman {spear:+.3f}")
    print(f"  P3: gamma~occupancy {r_occ:+.3f} | gamma~logmu {r_mu:+.3f} "
          f"| 偏(occ|mu) {p_occ:+.3f} | 偏(mu|occ) {p_mu:+.3f}")
    print("=" * 72)

    np.savez_compressed(od / f"beta_g_{a.tag}.npz",
                        gene=gene_names.astype("U16"), beta=beta, beta_r2=beta_r2,
                        beta_icpt=beta_icpt, rho=rho_mat, lags=lags,
                        gamma=gamma, mu8=mu8[idx_curves], occupancy=occupancy,
                        cv2=cv2, alpha_norm=alpha_norm, ds=ds_arr,
                        valid=ok)
    (od / f"shape_{a.tag}.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"写出 {od}/beta_g_{a.tag}.npz 与 shape_{a.tag}.json")


if __name__ == "__main__":
    main()

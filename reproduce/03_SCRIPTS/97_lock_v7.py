r"""97_lock_v7.py — 数字锁 v7：把机制层（理论/协方差/操控/窗口）并入唯一权威数字表

原则（沿用 v2-v6）
------------------
1. 只增不改：v6 的每个键原样搬运，新键追加。
2. 现场重算：所有数字从源 npz/json 现算，不手抄。
3. 每个块带 source，指到能重跑出它的脚本。

v7 新增
-------
  theory      合成场验证（恒等式误差、白噪声、深度不变性、gamma vs beta_emp）
  rho_shape   rho(h) 的局部斜率与 lag 漂移（脑为何没有单一 beta）
  gamma_beta  逐基因 gamma vs beta 的对应（含 log-log Pearson、Deming 斜率）
  window_scan gamma/beta 双侧窗口矩阵
  intervene   单基因置换 + 多基因剂量-反应

用法
----
  python 97_lock_v7.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "02_OUT"
V6 = OUT / "NUMBERS_LOCK_v6.json"
V7 = OUT / "NUMBERS_LOCK_v7.json"
DOWNSTREAM = OUT / "downstream" / "downstream_colon.json"

SETS = {
    "colon": ("02_OUT/gene_moments.npz", "02_OUT/variogram/beta_g_colon.npz",
              "02_OUT/variogram/shape_colon.json"),
    "lung": ("D3_LungCancer/02_OUT/gene_moments_noMT.npz", "02_OUT/variogram/beta_g_lung.npz",
             "02_OUT/variogram/shape_lung.json"),
    "brain": ("D2_MouseBrain/02_OUT/gene_moments_noMT.npz", "02_OUT/variogram/beta_g_brain.npz",
              "02_OUT/variogram/shape_brain.json"),
}


def r(x):
    return None if x is None else (float(x) if not isinstance(x, (int, str, bool)) else x)


def fit_slope(x, y):
    A = np.vstack([x, np.ones_like(x)]).T
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ c
    ss = ((y - pred) ** 2).sum()
    tot = ((y - y.mean()) ** 2).sum()
    return float(-c[0]), float(1 - ss / tot if tot > 0 else np.nan)


def deming(x, y):
    """errors-in-both-variables 斜率：真实斜率应接近 1 才说明两个估计同尺度。"""
    sx2, sy2 = np.var(x, ddof=1), np.var(y, ddof=1)
    sxy = np.cov(x, y, ddof=1)[0, 1]
    d = sy2 - sx2
    return float((d + np.sqrt(d ** 2 + 4 * sxy ** 2)) / (2 * sxy))


def block_theory():
    p = OUT / "theory/theory_check.json"
    d = json.loads(p.read_text(encoding="utf-8"))
    e1 = d["experiments"]["E1_identity"]
    e2 = d["experiments"]["E2_families"]
    e3 = d["experiments"]["E3_depth"]
    v = d["verdict"]
    fams = {}
    if isinstance(e2, dict):
        for k, val in e2.items():
            if isinstance(val, dict):
                fams[k] = {kk: r(vv) for kk, vv in val.items()
                           if isinstance(vv, (int, float))}
    return dict(
        source="86_theory_alpha_cv.py",
        identity=dict(max_abs_rel_err=float(v["max_abs_rel_err_overall"]),
                      median_abs_rel_err=float(v["median_abs_rel_err"]),
                      threshold=float(v["threshold"]),
                      passed=bool(v["pass_identity"]),
                      note="alpha = CV^2 * rhobar 在合成对数正态场上的最大相对误差"),
        gamma_vs_beta_emp=dict(slope=float(v["gamma_vs_beta_emp_slope"]),
                               r2=float(v["gamma_vs_beta_emp_r2"]),
                               mean_diff=float(v["gamma_minus_beta_emp_mean"]),
                               max_abs_diff=float(v["gamma_minus_beta_emp_max_abs"])),
        gamma_vs_nominal_beta_max_dev=float(v["gamma_vs_nominal_beta_max_dev"]),
        depth_gamma_range=float(v["depth_gamma_range"]),
        depth_note="1000 倍稀释下 gamma 的变幅（应接近 0：alpha 与深度无关）",
        families=fams,
        e1_raw={kk: r(vv) for kk, vv in e1.items() if isinstance(vv, (int, float))} if isinstance(e1, dict) else {},
        e3_raw={kk: r(vv) for kk, vv in e3.items() if isinstance(vv, (int, float))} if isinstance(e3, dict) else {},
    )


def block_rho_shape():
    out = {}
    for tag, (_, bt_p, _) in SETS.items():
        bt = np.load(ROOT / bt_p, allow_pickle=True)
        rho = bt["rho"].astype(float)
        lags = bt["lags"].astype(float)
        v = bt["valid"].astype(bool)
        med = np.nanmedian(rho[v], axis=0)
        rows = []
        for i in range(1, len(lags) - 1):
            if med[i] > 0 and med[i + 1] > 0:
                ls = -(np.log(med[i + 1]) - np.log(med[i])) / np.log(lags[i + 1] / lags[i])
            else:
                ls = np.nan
            rows.append(dict(lag=int(lags[i]), h_um=int(lags[i] * 8),
                             med_rho=float(med[i]),
                             local_slope=None if not np.isfinite(ls) else float(ls)))
        # 局部斜率对 lag 的 Spearman（显著为正 = 不是幂律）
        sl = [(x["h_um"], x["local_slope"]) for x in rows if x["local_slope"] is not None]
        hh = np.array([a for a, _ in sl], float)
        ss = np.array([b for _, b in sl], float)
        trend = float(spearmanr(hh, ss).statistic)
        # 分段中位斜率
        seg = {}
        for lo, hi in [(8, 16), (24, 64), (72, 128), (136, 248)]:
            m = (hh >= lo) & (hh <= hi)
            if m.sum():
                seg[f"{lo}-{hi}um"] = float(np.median(ss[m]))
        b, b_r2 = fit_slope(np.log(lags[2:]), np.log(np.where(med[2:] > 0, med[2:], np.nan)))
        out[tag] = dict(n_valid=int(v.sum()), by_lag=rows,
                        slope_vs_lag_spearman=trend,
                        slope_by_range=seg,
                        powerlaw_full_range=bool(abs(trend) < 0.5),
                        median_rho_at_lag1=float(med[1]),
                        note=("局部斜率随 lag 单调上升 = rho(h) 不是幂律，不存在单一 beta；"
                              "slope_by_range 给出各段的中位局部斜率"))
    return dict(source="97_lock_v7.py 从 beta_g_*.npz 现场算", by_tissue=out)


def block_gamma_beta():
    out = {}
    for tag, (mo_p, bt_p, sh_p) in SETS.items():
        mo = np.load(ROOT / mo_p, allow_pickle=True)
        ds = np.array([int(x) for x in mo["ds"]], float)
        genes_mo = np.array([str(g) for g in mo["genes"]])
        MU = np.column_stack([mo[f"mu_{int(d)}"] for d in ds])
        VA = np.column_stack([mo[f"var_{int(d)}"] for d in ds])
        with np.errstate(invalid="ignore", divide="ignore"):
            AL = (VA - MU) / np.where(MU > 0, MU ** 2, np.nan)

        bt = np.load(ROOT / bt_p, allow_pickle=True)
        bg = np.array([str(g) for g in bt["gene"]])
        beta = bt["beta"].astype(float)
        valid = bt["valid"].astype(bool)
        pos = {g: i for i, g in enumerate(genes_mo)}
        ii = np.array([pos[g] for g in bg])
        A = AL[ii]

        ok = np.all(np.isfinite(A) & (A > 0), axis=1)
        x = np.log(ds)
        xm = x - x.mean()
        y = np.log(np.where(ok[:, None], A, 1.0))
        ym = y - y.mean(axis=1, keepdims=True)
        g = np.full(bg.size, np.nan)
        g[ok] = -np.sum(xm[None, :] * ym, axis=1)[ok] / np.sum(xm ** 2)

        m = valid & np.isfinite(g) & np.isfinite(beta) & (g > 0) & (beta > 0)
        xb, yg = np.log(beta[m]), np.log(g[m])
        sh = json.loads((ROOT / sh_p).read_text(encoding="utf-8"))
        out[tag] = dict(
            n=int(m.sum()),
            gamma_median=float(np.median(g[m])),
            beta_median=float(np.median(beta[m])),
            diff_median=float(np.median(g[m] - beta[m])),
            diff_mean=float(np.mean(g[m] - beta[m])),
            diff_sd=float(np.std(g[m] - beta[m])),
            pearson_raw=float(pearsonr(g[m], beta[m]).statistic),
            spearman=float(spearmanr(g[m], beta[m]).statistic),
            pearson_loglog=float(pearsonr(yg, xb).statistic),
            ols_slope_loglog=float(np.polyfit(xb, yg, 1)[0]),
            deming_slope_loglog=deming(xb, yg),
            beta_window_lag=[2, 24],
            beta_window_note="87_variogram_shape.py 默认 fit-kmin 2 / fit-kmax 24（8um 网格 → 16-192um）",
            source="87_variogram_shape.py (beta) + 41_gene_curves.py 口径 (gamma)，本脚本合并",
            shape_json=dict(
                P2_pearson=sh["P2_correlation"]["pearson"],
                P2_spearman=sh["P2_correlation"]["spearman"],
                P3=sh["P3_occupancy"],
            ),
        )
    return dict(source="97_lock_v7.py 现场算", by_tissue=out,
                note=("deming_slope_loglog 是 errors-in-both-variables 斜率："
                      "两个估计都有噪声时，OLS 斜率会被系统性压低；Deming 无偏。"
                      "理论预测 log gamma = log beta，即斜率 1、截距 0"))


def block_window_scan():
    out = {}
    for tag in SETS:
        p = OUT / f"variogram/windowscan_{tag}.json"
        if not p.exists():
            continue
        d = json.loads(p.read_text(encoding="utf-8"))
        out[tag] = dict(n=d["n"], grid=d["grid"],
                        best_by_pearson=d["best_by_pearson"],
                        best_by_absdiff=d["best_by_absdiff"],
                        beta_drift_with_lag=d["beta_drift_with_lag"])
    return dict(source="91_window_scan.py", by_tissue=out,
                note=("行=gamma 拟合窗口，列=beta lag 窗口。"
                      "注意：窗口不是按 r 挑的，现稿已预先钉死 2-256（见正文）；"
                      "本表只作敏感性"))


def block_downstream():
    """下流分析代价：方差误设（解析）、差异表达（实证）、映射代理（实证）。

    原始数据由 106_downstream_v4.py 产出；这里只做结构化与派生，不重算统计量。
    """
    if not DOWNSTREAM.exists():
        raise SystemExit(f"缺少 {DOWNSTREAM}，先跑 106_downstream_v4.py")
    d = json.loads(DOWNSTREAM.read_text(encoding="utf-8"))

    od = {}
    for s, v in d["part_a_overdispersion"].items():
        od[str(int(s))] = dict(alpha=v["alpha_median"], mu=v["mu_median"],
                               phi=v["phi_median"], fpr=v["fpr_median"],
                               neff=v["n_eff_ratio"])
    de = {str(int(r["scale"])): {k: r[k] for k in
          ("n_a", "n_b", "top100_recall", "jaccard_sig", "n_sig_ref", "n_sig_s",
           "logfc_slope", "logfc_r", "logfc_slope_raw", "logfc_raw_r",
           "stat_spearman")}
          for r in d["part_b_de"]}
    mp = d["part_c_mapping"]

    ks = sorted(od, key=int)
    lo, hi = ks[0], ks[-1]

    def mvals(scale, field, ver):
        return [mp[scale][pn][ver][field] for pn in mp[scale]
                if ver in mp[scale][pn]]

    raw_slopes = [v for s in mp for pn in mp[s] for v in mvals(s, "slope", "raw")]
    raw_rs = [v for s in mp for pn in mp[s] for v in mvals(s, "pearson", "raw")]
    log_slopes = [v for s in mp for pn in mp[s] for v in mvals(s, "slope", "log1p")]
    # 浆细胞集含 IGKC（占 α 权重 66%），空间聚集最强，单列出来避免用宽范围掩盖结构
    not_plasma = [mp[s][pn]["raw"]["slope"] for s in mp for pn in mp[s]
                  if pn != "plasma_B" and "raw" in mp[s][pn]]

    return dict(
        source="106_downstream_v4.py", section=d["section"],
        ref_scale_um=8, region_definition="marker top/bottom 30% at 16 um",
        od=od, de=de, map=mp,
        summary=dict(
            phi_lo=od[lo]["phi"], phi_hi=od[hi]["phi"],
            fpr_lo=od[lo]["fpr"], fpr_hi=od[hi]["fpr"],
            neff_hi=od[hi]["neff"], mu_hi=od[hi]["mu"], mu_lo=od[lo]["mu"],
            jacc_32=de["32"]["jaccard_sig"], jacc_16=de["16"]["jaccard_sig"],
            top100_16=de["16"]["top100_recall"], top100_64=de["64"]["top100_recall"],
            sig_ref=de["16"]["n_sig_ref"], sig_128=de["128"]["n_sig_s"],
            slope_log_64=de["64"]["logfc_slope"], slope_raw_64=de["64"]["logfc_slope_raw"],
            r_64=de["64"]["logfc_r"],
            map_raw_min=min(raw_slopes), map_raw_max=max(raw_slopes),
            map_raw_min_x=min(not_plasma), map_raw_max_x=max(not_plasma),
            map_raw_r_min=min(raw_rs), map_log_min=min(log_slopes),
            map_log_max=max(log_slopes),
        ),
    )


def block_intervene():
    out = {}
    p = OUT / "excise/excise_colon.json"
    if p.exists():
        d = json.loads(p.read_text(encoding="utf-8"))
        cases = {}
        for case, obj in d["pooled"]["cases"].items():
            cases[case] = {kk: float(vv) for kk, vv in obj.items()
                           if isinstance(vv, (int, float)) and not isinstance(vv, bool)}
        genes = {}
        for gn, v in d["genes"].items():
            genes[gn] = dict(
                total=float(v["total"]), occupancy=float(v["occupancy"]),
                hot_share_of_counts=float(v["hot_share_of_counts"]),
                baseline_gamma=float(v["baseline"]["gamma"]),
                permute_gamma=float(v["permute"]["gamma_median"]),
                permute_lo=float(v["permute"]["gamma_lo"]),
                permute_hi=float(v["permute"]["gamma_hi"]),
                excise_gamma=float(v["excise"]["gamma"]),
                isolate_gamma=float(v["isolate"]["gamma"]),
                delta_permute=float(v["permute"]["gamma_median"] - v["baseline"]["gamma"]),
            )
        out["single_gene"] = dict(
            source="88_region_excise.py", nperm=int(d["nperm"]),
            pooled=dict(target=d["pooled"]["target"],
                        gamma_baseline=float(d["pooled"]["gamma_baseline"]),
                        delta_permute=float(d["pooled"]["delta_permute"]),
                        delta_excise=float(d["pooled"]["delta_excise"]),
                        delta_isolate=float(d["pooled"]["delta_isolate"]),
                        cases=cases),
            genes=genes,
            summary=dict(
                n_genes=len(genes),
                delta_min=float(min(v["delta_permute"] for v in genes.values())),
                delta_max=float(max(v["delta_permute"] for v in genes.values())),
                n_positive=int(sum(1 for v in genes.values() if v["delta_permute"] > 0)),
                base_min=float(min(v["baseline_gamma"] for v in genes.values())),
                base_max=float(max(v["baseline_gamma"] for v in genes.values())),
                note="六个基因各自独立置换；全部为正 = 方向一致，不是单点证据"),
            design=("每个基因独立在 2um 域内随机置换其非零计数位置，"
                    "边缘分布（总量、非零值集合）完全不变，只打乱空间位置。"
                    "因此 gamma 的变化只能归因于空间结构"))

    p2 = OUT / "excise/multigene_permute.json"
    if p2.exists():
        d = json.loads(p2.read_text(encoding="utf-8"))
        out["dose_response"] = dict(
            source="89_multigene_permute.py", ds=d["ds"],
            curves={k: dict(n_permuted=v["n_permuted"], mu2_share=v["mu2_share"],
                            gamma=v["gamma"], r2=v["r2"])
                    for k, v in d["curves"].items()})
    return out


def main():
    v6 = json.loads(V6.read_text(encoding="utf-8"))
    v7 = dict(v6)  # 只增不改：全部键原样继承
    added = []

    v7["theory"] = block_theory(); added.append("theory")
    v7["rho_shape"] = block_rho_shape(); added.append("rho_shape")
    v7["gamma_beta"] = block_gamma_beta(); added.append("gamma_beta")
    v7["window_scan"] = block_window_scan(); added.append("window_scan")
    v7["intervene"] = block_intervene(); added.append("intervene")
    v7["downstream"] = block_downstream(); added.append("downstream")

    v7["meta"] = dict(v6.get("meta", {}))
    v7["meta"]["version"] = 7.0
    v7["meta"]["built"] = "2026-11-15"
    v7["meta"]["supersedes"] = "NUMBERS_LOCK_v6.json"
    v7["meta"]["note"] = (
        "唯一权威数字表。v7 = v6 + 机制层：theory（合成场验证恒等式与标度）、"
        "rho_shape（rho(h) 的局部斜率漂移）、gamma_beta（逐基因对应，含 Deming 斜率）、"
        "window_scan（双侧窗口敏感性）、intervene（单基因置换 + 剂量-反应）、"
        "downstream（方差误设的假阳性代价、差异表达与映射的跨尺度偏差）。"
        "v6 及以前的每个键原样保留。")
    v7["meta"]["what_changed"] = list(v6.get("meta", {}).get("what_changed", [])) + [
        f"新增 {k}（v7）" for k in added
    ]

    V7.write_text(json.dumps(v7, indent=1, ensure_ascii=False), encoding="utf-8")
    print(f"写出 {V7}")
    print(f"顶层键 {len(v6)} -> {len(v7)}；新增：{', '.join(added)}")

    # ---- 摘要打印，便于核对
    t = v7["theory"]
    print(f"\n[theory] 恒等式最大相对误差 {t['identity']['max_abs_rel_err']:.4f} "
          f"(阈值 {t['identity']['threshold']})  gamma-vs-beta_emp r2={t['gamma_vs_beta_emp']['r2']:.4f} "
          f"slope={t['gamma_vs_beta_emp']['slope']:.3f}")
    for tag, d in v7["rho_shape"]["by_tissue"].items():
        print(f"[rho_shape] {tag:>6} 斜率-vs-lag Spearman={d['slope_vs_lag_spearman']:+.3f} "
              f"幂律={d['powerlaw_full_range']} 分段={d['slope_by_range']}")
    for tag, d in v7["gamma_beta"]["by_tissue"].items():
        print(f"[gamma_beta] {tag:>6} n={d['n']:>5} g={d['gamma_median']:.3f} "
              f"b={d['beta_median']:.3f} diff={d['diff_median']:+.3f} "
              f"r={d['pearson_raw']:.3f} rho={d['spearman']:.3f} "
              f"r_log={d['pearson_loglog']:.3f} deming={d['deming_slope_loglog']:.3f}")
    s = v7["downstream"]["summary"]
    print(f"[downstream] phi {s['phi_lo']:.2f}->{s['phi_hi']:.1f}  "
          f"FPR {s['fpr_lo']:.3f}->{s['fpr_hi']:.3f}  n_eff {s['neff_hi']:.3f}  "
          f"Jacc@32={s['jacc_32']:.2f}  top100@64={s['top100_64']:.2f}  "
          f"slope log/raw @64={s['slope_log_64']:.2f}/{s['slope_raw_64']:.2f}")


if __name__ == "__main__":
    main()

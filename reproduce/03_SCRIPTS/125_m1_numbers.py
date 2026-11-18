r"""125_m1_numbers.py — M1 的四个数，现场从原始产物算出，不抄任何已有表格。

产出 02_OUT/m1_v551_numbers.json，供图与稿共用。四项：

1. pooled 口径修正：`excise`／`isolate` 两个操作会改变目标基因的 μ（进而改变它在
   pooled α 里的 μ² 权重）。旧脚本（88_region_excise.py）换曲线却没换权重，
   报出的 0.502／0.818 是"权重冻结"的估计量，不是合并场的 α。
   这里按定义 α_pool = Σ μ²α / Σ μ² 重算，权重取操作后的 μ。
   单基因时 α_pool ≡ 该基因的 α（Σ 只有一项），故 isolate 的 pooled 值必然等于
   per-gene isolate γ，即 Table S3 的 0.839——这一条不需要算术，是恒等式。

2. P1 复现：fold、Δγ、块自助 z。

3. 七片 fold 与端点值（给 Fig 1b）。

4. 0.818/0.839 的口径核验：把三种可能的读法都算出来，写清楚哪个对应哪句话。
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "02_OUT" / "m1_v551_numbers.json"


def pooled_alpha(mu, var):
    """α_pool(d) = Σ_g μ_g² α_g / Σ_g μ_g² — 等价于对合并场直接算 (Var−E)/E²。"""
    with np.errstate(invalid="ignore", divide="ignore"):
        alpha = (var - mu) / np.where(mu > 0, mu ** 2, np.nan)
    w = np.where(np.isfinite(alpha), mu ** 2, 0.0)
    num = np.nansum(np.where(np.isfinite(alpha), w * alpha, 0.0), axis=0)
    den = np.nansum(np.where(np.isfinite(alpha), w, 0.0), axis=0)
    return num / den, alpha, w


def fit_gamma(ds, alpha):
    ds = np.asarray(ds, float)
    al = np.asarray(alpha, float)
    m = np.isfinite(al) & (al > 0)
    x, y = np.log(ds[m]), np.log(al[m])
    A = np.vstack([x, np.ones_like(x)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ coef
    r2 = 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    return float(-coef[0]), float(r2)


def main() -> int:
    z = np.load(ROOT / "02_OUT" / "gene_moments.npz", allow_pickle=True)
    ds = [int(x) for x in z["ds"]]
    genes = np.array([str(g) for g in z["genes"]])
    MU = np.column_stack([z[f"mu_{d}"] for d in ds])
    VA = np.column_stack([z[f"var_{d}"] for d in ds])
    ap_base, AL, W = pooled_alpha(MU, VA)
    g_base, r2_base = fit_gamma(ds, ap_base)

    ex = json.loads((ROOT / "02_OUT/excise/excise_colon.json").read_text(encoding="utf-8"))
    ti = int(np.where(genes == "IGKC")[0][0])
    ig = ex["genes"]["IGKC"]

    def curve_arrays(op):
        c = ig[op]["curve"]
        mu = np.array([c[str(d)]["mu"] for d in ds])
        al = np.array([c[str(d)]["alpha"] for d in ds])
        return mu, al

    # 其余基因的贡献：S_other = Σ_{g≠t} μ²α，W_other = Σ_{g≠t} μ²
    keep = np.isfinite(AL[ti]) & (AL[ti] > 0)
    S_all = np.where(np.isfinite(AL), W * np.where(np.isfinite(AL), AL, 0), 0).sum(axis=0)
    W_all = np.where(np.isfinite(AL), W, 0).sum(axis=0)
    S_other = S_all - np.where(keep, W[ti] * np.where(keep, AL[ti], 0), 0)
    W_other = W_all - np.where(keep, W[ti], 0)

    cases = {}
    for op in ("baseline", "permute", "excise", "isolate"):
        if op == "permute":
            mu = MU[ti].copy()
            al = np.full(len(ds), np.nan)
            # permute 不改变 μ，逐次曲线未存中位 α(d)；用 20 次里存的 gamma 中位数与
            # 已存的 pooled case（权重本来就对）核对，不重算。
            cases[op] = dict(note="μ 不变，旧权重正确；沿用 excise_colon.json 的 pooled case",
                             gamma_frozen=ex["pooled"]["cases"]["permute"]["gamma"])
            continue
        mu, al = curve_arrays(op)
        num = S_other + mu ** 2 * al
        den = W_other + mu ** 2
        ap = num / den
        g, r2 = fit_gamma(ds, ap)
        # 冻结权重的旧口径（复现旧脚本的行为），留作对照
        num_f = S_other + W[ti] * al
        den_f = W_other + W[ti]
        gf, _ = fit_gamma(ds, num_f / den_f)
        cases[op] = dict(gamma=g, r2=r2, gamma_frozen_weight=gf,
                         mu=mu.tolist(), alpha=al.tolist(), alpha_pool=ap.tolist(),
                         weight_change_2um=float((mu[0] / MU[ti][0]) ** 2))

    # 单基因恒等式核验：isolate 的 pooled 必须等于该基因自己的 γ
    iso_gene = ig["isolate"]["gamma"]
    ident_gap = abs(cases["isolate"]["gamma"] - iso_gene)

    rec = {
        "source": "125_m1_numbers.py",
        "ds": ds,
        "pooled_baseline": dict(gamma=g_base, r2=r2_base, alpha_pool=ap_base.tolist()),
        "pooled_cases": cases,
        "identity_isolate_minus_pergene": ident_gap,
        "note": "gamma 为 α ∝ d^-gamma 的指数（正数）；lock 里 slope 存成负数",
    }

    # ---------- P1 复现 ----------
    p1 = json.loads(Path(
        r"E:\SPATIAL_2026\MS_Visium_HD_Human_Colon_Cancer_P1\02_OUT"
        r"\metrics_Visium_HD_Human_Colon_Cancer_P1.json").read_text(encoding="utf-8"))
    b1 = json.loads(Path(
        r"E:\SPATIAL_2026\MS_Visium_HD_Human_Colon_Cancer_P1\02_OUT"
        r"\block_boot\boot_colonP1.json").read_text(encoding="utf-8"))
    b0 = json.loads((ROOT / "02_OUT/block_boot/boot_colon.json").read_text(encoding="utf-8"))
    a = p1["alpha"]
    fold_p1 = a[0] / a[-1]
    sec = {s["key"]: s for s in json.loads(
        (ROOT / "02_OUT/extra_tissues_20261116.json").read_text(encoding="utf-8"))["sections"]}
    colon = sec["colon"]
    r1, r0 = b1["results"]["512"], b0["results"]["512"]
    d_gamma = r0["mean"] - r1["mean"]
    z = d_gamma / np.sqrt(r0["sd"] ** 2 + r1["sd"] ** 2)
    rec["p1_replicate"] = dict(
        label=p1["label"],
        alpha_2um=a[0], alpha_256um=a[-1], fold=fold_p1,
        gamma=p1["gamma_pool"], gamma_median=p1["gamma_median"],
        colon_fold=colon["drop_fold"], colon_gamma=-colon["gamma"],
        fold_gap_pct=100 * (colon["drop_fold"] - fold_p1) / colon["drop_fold"],
        delta_gamma=colon["gamma"] + p1["gamma_pool"],
        boot_p1=[r1["mean"], r1["sd"]], boot_colon=[r0["mean"], r0["sd"]],
        boot_delta=d_gamma, boot_z=float(z),
        n_genes=p1["n_genes"], n_fit=p1["n_fit"],
        source="E:/SPATIAL_2026/MS_Visium_HD_Human_Colon_Cancer_P1/02_OUT",
    )

    # ---------- 七片端点（Fig 1b）----------
    rec["seven"] = [
        dict(key=s["key"], label=s["label"], fold=s["drop_fold"],
             gamma=-s["gamma"], alpha_by_d={k: v for k, v in s["alpha_by_d"].items()})
        for s in sec.values()]

    OUT.write_text(json.dumps(rec, ensure_ascii=False, indent=1), encoding="utf-8")

    print("== pooled（权重随操作更新）==")
    print(f"  baseline            γ={g_base:.4f}  r2={r2_base:.4f}")
    for op in ("permute", "excise", "isolate"):
        c = cases[op]
        if "gamma" in c:
            print(f"  {op:<18s} γ={c['gamma']:.4f}  (冻结权重旧值 {c['gamma_frozen_weight']:.4f})"
                  f"  μ²权重×{c['weight_change_2um']:.3f}")
        else:
            print(f"  {op:<18s} γ={c['gamma_frozen']:.4f}  (μ 不变，无修正)")
    print(f"  单基因恒等式 |pooled_isolate − per-gene isolate| = {ident_gap:.2e}")
    print("\n== P1 复现 ==")
    for k, v in rec["p1_replicate"].items():
        print(f"  {k:20s} {v}")
    print("\n== 七片 ==")
    for s in rec["seven"]:
        print(f"  {s['key']:24s} fold={s['fold']:6.2f} γ={s['gamma']:.3f}")
    print(f"\nwritten {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

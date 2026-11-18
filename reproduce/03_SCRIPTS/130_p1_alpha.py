r"""130_p1_alpha.py — 从 P1 的逐基因矩重建池化 α(d) 曲线（v5.8 用）。

背景
----
v5.7 把 P5 画进 Fig 1，却把 P1（13.1 倍）留在图外。两者的差别只有"P5 最大"，
审稿人会问为什么挑最大的画。v5.8 改为九片全画，因此需要 P1 的八档池化 α。

P1 没有留下 pooled_alpha.json（当年只存了 2 µm 与 256 µm 两个锚点和 γ），
但 gene_moments_noMT.npz 里 mu_d / var_d 齐全，池化 α 可直接重建：

    alpha_pool(d) = SUM_g ( var_g(d) - mu_g(d) ) / SUM_g mu_g(d)^2

口径与 07_pooled_alpha.py 完全一致（无基因筛选，全部 18,074 基因参与）。

三道回验（不过就停，不出数）
--------------------------
  alpha(2 µm)   = 29.161157756622114   (m1_v551_numbers.json :: p1_replicate)
  alpha(256 µm) =  2.226584791156037
  fold          = 13.096809909260953
  gamma         =  0.5485115671359143

输出
----
  02_OUT/p1_alpha_curve_v58.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

NPZ = Path(r"E:\SPATIAL_2026\MS_Visium_HD_Human_Colon_Cancer_P1\02_OUT\gene_moments_noMT.npz")
OUT = Path(__file__).resolve().parents[1] / "02_OUT" / "p1_alpha_curve_v58.json"

# 锚点（来自 m1_v551_numbers.json :: p1_replicate，锁死）
A2_TRUE = 29.161157756622114
A256_TRUE = 2.226584791156037
FOLD_TRUE = 13.096809909260953
GAMMA_TRUE = 0.5485115671359143

TOL = 5e-6      # 相对误差容限：矩重建与当年逐档累加只应有浮点级差异


def main():
    z = np.load(NPZ, allow_pickle=False)
    ds = [int(x) for x in z["ds"]]
    if ds != [2, 4, 8, 16, 32, 64, 128, 256]:
        raise SystemExit(f"d 序列不是八档：{ds}")

    alpha, n_cells = {}, {}
    for d in ds:
        mu = z[f"mu_{d}"].astype(np.float64)
        var = z[f"var_{d}"].astype(np.float64)
        num = float(np.sum(var - mu))
        den = float(np.sum(mu ** 2))
        alpha[d] = num / den
        n_cells[d] = int(z["n_cells"][ds.index(d)])
        print(f"  d={d:>3} µm  bins={n_cells[d]:>9,}  alpha_pool={alpha[d]:.6f}", flush=True)

    fold = alpha[2] / alpha[256]
    x = np.log(np.array(ds, float))
    y = np.log(np.array([alpha[d] for d in ds]))
    slope = float(np.polyfit(x, y, 1)[0])
    gamma = -slope

    # ---- 三道回验 ----
    checks = [
        ("alpha(2 µm)", alpha[2], A2_TRUE),
        ("alpha(256 µm)", alpha[256], A256_TRUE),
        ("fold 2→256", fold, FOLD_TRUE),
        ("gamma", gamma, GAMMA_TRUE),
    ]
    bad = []
    for name, got, want in checks:
        rel = abs(got - want) / abs(want)
        ok = rel <= TOL
        print(f"  {name:14s} 重建 {got:.6f}  锁定 {want:.6f}  相对误差 {rel:.2e}  {'OK' if ok else 'FAIL'}")
        if not ok:
            bad.append(name)
    if bad:
        raise SystemExit(f"回验失败：{bad}，不写输出")

    umi = float(np.sum(z[f"mu_{2}"]) * n_cells[2])
    payload = {
        "key": "colon_p1",
        "label": "Visium HD Human Colon Cancer P1",
        "origin": "130_p1_alpha.py :: gene_moments_noMT.npz（三锚点回验通过）",
        "alpha_2um": alpha[2],
        "alpha_256um": alpha[256],
        "drop_fold": fold,
        "gamma": -gamma,          # 与 extra_tissues_20261116.json 同号约定（拟合斜率为负）
        "n_bins_2um": n_cells[2],
        "umi_total": umi,
        "umi_per_bin": umi / n_cells[2],
        "n_genes_noMT": int(z["genes"].shape[0]),
        "alpha_by_d": {str(d): alpha[d] for d in ds},
        "checks": {name: dict(rebuild=got, locked=want) for name, got, want in checks},
    }
    OUT.write_text(json.dumps(payload, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"\n  写出 {OUT}")
    print(f"  P1: {fold:.2f}-fold, gamma={gamma:.4f}, {umi / n_cells[2]:.1f} transcripts/bin")


if __name__ == "__main__":
    main()

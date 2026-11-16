r"""
42_build_lock_v3.py — 汇总为唯一数字表 NUMBERS_LOCK_v3.json。

v3 = v2（七个新分析）+ v1 的 synth / spread / window（v2 未搬过来的三节）
     + 本轮新增：fill（填充率敏感性）、boot（块自助）、layer_weights（层权重）
     + derived3：所有派生量（F 检验 p、四象限 SD、层覆盖、样本量换算）

规则：正文任何数字必须能在 v3 里找到键；v3 的每个数都由脚本算出，不手填。
"""
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import f as fdist

ROOT = Path(__file__).resolve().parents[1]
O = ROOT / "02_OUT"

v2 = json.loads((O / "NUMBERS_LOCK_v2.json").read_text(encoding="utf-8"))
v1 = json.loads((O / "NUMBERS_LOCK_v1.json").read_text(encoding="utf-8"))

out = dict(v2)
out["meta"] = dict(
    version="3.0",
    built="2026-11-14",
    supersedes="NUMBERS_LOCK_v2.json",
    note="唯一权威数字表；正文任何数字须在此有键。v3 在 v2 之上并入 v1 的 synth/spread/window，"
         "并新增 fill / boot / layer_weights 三节与 derived3。",
    what_changed=[
        "并入 v1 的 synth（独立细胞合成对照）、spread（四象限/自助离散度）、window（拟合窗口敏感性）",
        "新增 fill：严格填充率口径（只保留完整在组织内的粗格）下的 alpha(d) 与 gamma",
        "新增 boot：512/1024/2048 um 超块自助 B=200 的 gamma 分布",
        "新增 layer_weights：6 个表达层各占 Sigma(mu^2) 的份额（32 um 口径）",
        "新增 derived3：F 检验 p、四象限 SD、层覆盖、样本量换算、各操控对 gamma 的效应量",
    ],
)
out["legacy"] = dict(synth=v1["synth"], spread=v1["spread"], window=v1["window"],
                     source="NUMBERS_LOCK_v1.json")

# ---------------------------------------------------------------- 统一窗口到 2-256 um
# v2 的 models 节用 9 个尺度（含 512 um，colon 仅 155 个 bin）拟合，与正文声明的
# 2-256 um 窗口不一致；这里换成 8 点重跑结果，旧的留档为 models_9pt。
p8 = O / "xi_model_8pt" / "xi_model_8pt.json"
if p8.exists():
    m8 = json.loads(p8.read_text(encoding="utf-8"))
    out["models_9pt"] = dict(by_tissue=v2["models"],
                             note="9 点拟合（含 512 um），仅留档；正文不使用")
    out["models"] = {t: dict(
        aic={k: v["aic"] for k, v in m8[t]["models"].items()},
        daic={k: v["dAIC"] for k, v in m8[t]["models"].items()},
        rss={k: v["rss"] for k, v in m8[t]["models"].items()},
        quad=m8[t]["models"]["quad"]["params"],
        quad_F=m8[t]["models"]["quad"]["F"],
        quad_p_value=m8[t]["models"]["quad"]["p_F"],
        xi_hat=m8[t]["xi_profile"]["xi_hat"],
        xi_ci=[m8[t]["xi_profile"]["ci_lo"], m8[t]["xi_profile"]["ci_hi"]],
        local_slope=dict(seg_mid=m8[t]["local_slope"]["seg_mid"],
                         observed=m8[t]["local_slope"]["observed"],
                         predicted_xi=m8[t]["local_slope"]["predicted_xi"],
                         predicted_quad=m8[t]["local_slope"]["predicted_quad"],
                         rmse_xi=m8[t]["local_slope"]["rmse_xi"],
                         rmse_quad=m8[t]["local_slope"]["rmse_quad"]),
        xi_profile=m8[t]["xi_profile"],
        source="30_xi_model.py --max-d 256 (xi_model_8pt)",
    ) for t in ("colon", "lung", "brain")}
    out["meta"]["what_changed"].append(
        "models 节改为 8 点（2-256 um）拟合，与正文窗口一致；9 点版本留档为 models_9pt")

# ---------------------------------------------------------------- fill
fill = {}
for t in ("colon", "lung", "brain"):
    p = O / "fill_sensitivity" / f"fill_{t}.json"
    if p.exists():
        fill[t] = json.loads(p.read_text(encoding="utf-8"))
out["fill"] = dict(by_tissue=fill, source="38_fill_sensitivity.py")

# ---------------------------------------------------------------- boot
boot = {}
for t in ("colon", "lung", "brain"):
    p = O / "block_boot" / f"boot_{t}.json"
    if p.exists():
        boot[t] = json.loads(p.read_text(encoding="utf-8"))
out["boot"] = dict(by_tissue=boot, source="39_block_boot.py")

# ---------------------------------------------------------------- layer weights
p = O / "layer_maps" / "layer_weights.json"
if p.exists():
    out["layer_weights"] = dict(by_tissue=json.loads(p.read_text(encoding="utf-8")),
                                source="40_layer_maps.py")

# ---------------------------------------------------------------- per-gene reliability by layer
p = O / "splithalf_genes" / "reliability_layers.json"
if p.exists():
    out["reliability_layers"] = json.loads(p.read_text(encoding="utf-8"))

# ---------------------------------------------------------------- layers / variance（由 34 本次输出覆盖）
plv = O / "layer_variance" / "layer_variance.json"
if plv.exists():
    lv = json.loads(plv.read_text(encoding="utf-8"))
    out["layers"] = dict(
        edges=lv["edges"], ref_d=lv["ref_d"], by_tissue=lv["by_tissue"], below=lv["below"],
        source="34_layer_variance.py --lo -2 --hi 2 --step 0.5 --ref-d 32")
    out["variance"] = dict(pooled=lv["variance"], median=lv["variance_median"])
    out["meta"]["what_changed"].append(
        "layers 节改由 34 本次输出覆盖：每层新增 share（占 Sigma(mu^2)）与 n_gamma，"
        "并补一行 >=10 counts/bin 的合并层，使表中权重合计到 100%")

# ---------------------------------------------------------------- 主导基因自身指数 + CV 隐含指数
for t in ("colon", "lung", "brain"):
    pg = O / "gene_curves" / f"gene_curves_{t}.npz"
    if pg.exists():
        z = np.load(pg, allow_pickle=True)
        names = np.array([str(x) for x in z["gene"]])
        top = out["jackknife"][t]["top1_gene"]
        i = int(np.where(names == top)[0][0])
        out["jackknife"][t]["top1_gamma"] = float(z["gamma"][i])
        out["jackknife"][t]["top1_mu32"] = float(z["mu32"][i])

# ---------------------------------------------------------------- derived3
d3 = {}
for t in ("colon", "lung", "brain"):
    m = out["models"][t]
    d3[f"{t}_quad_F"] = float(m["quad_F"])
    d3[f"{t}_quad_F_p"] = float(m["quad_p_value"])
    d3[f"{t}_power_dAIC"] = float(m["daic"]["power"])
    d3[f"{t}_xi_dAIC_min"] = float(min(m["daic"]["exp"], m["daic"]["gauss"]))
    d3[f"{t}_xi_dAIC_max"] = float(max(m["daic"]["exp"], m["daic"]["gauss"]))
    q = [abs(v) for v in v2["datasets"][t]["quadrants"].values()]
    d3[f"{t}_quad_sd"] = float(np.std(q, ddof=1))
    d3[f"{t}_quad_min"] = float(min(q))
    d3[f"{t}_quad_max"] = float(max(q))
    d3[f"{t}_quad_mean"] = float(np.mean(q))
    d3[f"{t}_gamma_from_window32"] = float(v1["window"][t]["gamma_from_32"])
    d3[f"{t}_window_delta"] = float(v1["window"][t]["delta_2_to_32"])
    d3[f"{t}_window_range"] = float(v1["window"][t]["range"])
    cvv = v2["datasets"][t]["cv"]
    d3[f"{t}_cv_gamma"] = float(2.0 * math.log(cvv["2"] / cvv["256"]) / math.log(128.0))
    d3[f"{t}_pooled_minus_median"] = float(
        out["jackknife"][t]["gamma_full"] - out["gene_level"][t]["gamma_gene"]["median"])
    if "layer_weights" in out:
        lw = out["layer_weights"]["by_tissue"][t]
        d3[f"{t}_coverage6"] = float(lw["coverage_6layers"])
        d3[f"{t}_coverage8"] = float(lw["coverage_all8"])

# 半宽与样本量：以四带/四象限 SD 作为"区域内 SD"，它是样本间 SD 的下界
sd_col = v2["splithalf"]["quarter_sd"]
d3["colon_band_sd"] = float(sd_col)
d3["colon_band_range"] = float(v2["splithalf"]["quarter_range"])
d3["tissue_gap_colon_brain"] = float(v2["derived"]["colon_brain_gamma_gap"])
d3["tissue_gap_colon_lung"] = float(v2["derived"]["colon_lung_gamma_gap"])
d3["halfwidth95_regional_2sec"] = float(1.96 * math.sqrt(2) * sd_col)
d3["gap_over_regional_halfwidth"] = float(d3["tissue_gap_colon_brain"] / d3["halfwidth95_regional_2sec"])

for t in ("colon", "lung", "brain"):
    rms = v2["datasets"][t]["rms_gamma_cut"]
    d3[f"{t}_rms_cut"] = float(rms)
    for dlt in (0.19, 0.10, 0.05, 0.02):
        key = f"{t}_n_for_delta{str(dlt).replace('.', '')}"
        d3[key] = float((2.802 * rms / dlt) ** 2)
    d3[f"{t}_mde1_from_rms"] = float(2.802 * rms)
    d3[f"{t}_mde1_from_band"] = float(2.802 * sd_col)

d3["colon_n_for_delta_gap"] = float(
    (2.802 * v2["datasets"]["colon"]["rms_gamma_cut"] / d3["tissue_gap_colon_lung"]) ** 2)

d3["depth_delta_max"] = float(v2["derived"]["depth_gamma_delta_max"])
d3["synth_bias"] = float(v2["synthetic"]["expression"]["gamma_all"]["med"]
                         - v2["synthetic"]["expression"]["beta_true"])
d3["synth_range_across_layers"] = float(v2["derived"]["synth_gamma_range_across_layers"])
d3["layer_range_pooled"] = float(v2["variance"]["pooled"]["range_layer"])
d3["layer_range_median"] = float(v2["variance"]["median"]["range_layer"])
d3["tissue_range_pooled"] = float(v2["variance"]["pooled"]["range_tissue"])
d3["var_ratio_pooled"] = float(v2["variance"]["pooled"]["range_layer"] / v2["variance"]["pooled"]["range_tissue"])
d3["var_ratio_median"] = float(v2["variance"]["median"]["range_layer"] / v2["variance"]["median"]["range_tissue"])
d3["sigma_mu2_note"] = "coverage6 = 6 个表达层（mu32<=10）覆盖的 Sigma(mu^2) 份额"

out["derived3"] = d3

p = O / "NUMBERS_LOCK_v3.json"
p.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
print(f"写出 {p}  ({p.stat().st_size/1024:.0f} KB, {len(out)} 节)")
print("fill:", {k: round(v["gamma_delta"], 4) for k, v in fill.items()})
print("boot:", {k: {b: round(v["results"][b]["sd"], 4) for b in v["results"]} for k, v in boot.items()})
print("F p:", {t: round(d3[f'{t}_quad_F_p'], 5) for t in ('colon', 'lung', 'brain')})
print("quad sd:", {t: round(d3[f'{t}_quad_sd'], 4) for t in ('colon', 'lung', 'brain')})
print("coverage6:", {t: round(d3.get(f'{t}_coverage6', -1), 4) for t in ('colon', 'lung', 'brain')})
print("halfwidth:", round(d3["halfwidth95_regional_2sec"], 4),
      "gap/hw:", round(d3["gap_over_regional_halfwidth"], 3))

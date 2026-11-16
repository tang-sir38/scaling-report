r"""
21_build_numbers.py — 把全部分析 JSON 汇总成唯一数字源 NUMBERS_LOCK_v1.json。

正文/摘要/图注/表格/投稿信里的每个数字都由 30_render.py 从这里取值注入，
不允许在文档脚本里出现字面量数字（守卫 27_guard_numbers.py 会反查）。

用法：python 21_build_numbers.py --out 02_OUT/NUMBERS_LOCK_v1.json
"""
import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SETS = {
    "colon": {"alpha": ROOT / "02_OUT/pooled_alpha_noMT.json",
              "local": ROOT / "02_OUT/pooled_local_8xd_noMT.json",
              "boot": ROOT / "02_OUT/bootstrap_512um_noMT.json",
              "power": ROOT / "02_OUT/power_table_colon.json",
              "t1": ROOT / "02_OUT/table1_stats.json",
              "alpha_mt": ROOT / "02_OUT/pooled_alpha.json"},
    "lung": {"alpha": ROOT / "D3_LungCancer/02_OUT/pooled_alpha_noMT.json",
             "local": ROOT / "D3_LungCancer/02_OUT/pooled_local_8xd_noMT.json",
             "boot": ROOT / "D3_LungCancer/02_OUT/bootstrap_512um_noMT.json",
             "power": ROOT / "D3_LungCancer/02_OUT/power_table_lung.json",
             "t1": ROOT / "D3_LungCancer/02_OUT/table1_stats.json"},
    "brain": {"alpha": ROOT / "D2_MouseBrain/02_OUT/pooled_alpha_noMT.json",
              "local": ROOT / "D2_MouseBrain/02_OUT/pooled_local_8xd_noMT.json",
              "boot": ROOT / "D2_MouseBrain/02_OUT/bootstrap_512um_noMT.json",
              "power": ROOT / "D2_MouseBrain/02_OUT/power_table_brain.json",
              "t1": ROOT / "D2_MouseBrain/02_OUT/table1_stats.json",
              "alpha_mt": ROOT / "D2_MouseBrain/02_OUT/pooled_alpha.json"},
}


def load(p):
    return json.loads(Path(p).read_text(encoding="utf-8"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "02_OUT/NUMBERS_LOCK_v1.json"))
    args = ap.parse_args()

    out = {"datasets": {}, "synth": {}, "layers": {}, "meta": {}}

    for tag, files in SETS.items():
        a = load(files["alpha"])
        loc = load(files["local"])
        boot = load(files["boot"])
        pw = load(files["power"])
        t1 = load(files["t1"])
        ds = sorted(int(k) for k in a["per_d"])
        rec = {
            "n_bins": t1["n_bins_domain"],
            "n_bins_matrix": t1["n_bins_matrix"],
            "n_genes_raw": t1["n_genes_raw"],
            "n_genes_noMT": t1["n_genes_noMT"],
            "zero_frac": t1["zero_frac"],
            "umi_total": t1["umi_total"],
            "umi_per_bin": t1["umi_per_bin"],
            "umi_median": t1["umi_median"],
            "ds": ds,
            "alpha": {str(d): a["per_d"][str(d)]["alpha_pool"] for d in ds},
            "alpha_local": {str(d): loc["per_d"][str(d)]["alpha_pool"] for d in ds if str(d) in loc["per_d"]},
            "gamma": a["slope_all"],
            "gamma_local": loc["slope"],
            "segments": {f"{s['frm']}-{s['to']}": s["slope"] for s in a["segments"]},
            "boot_ci": boot["slope_ci95"],
            "boot_sd": boot["slope_boot_sd"],
            "quadrants": {k: v["slope"] for k, v in boot["quadrants"].items()},
            "quadrant_alpha": {k: v["alpha"] for k, v in boot["quadrants"].items()},
            "power": {},
            "rms_gamma_random": pw["rms_gamma_random"],
            "rms_gamma_cut": pw["rms_gamma_cut"],
            "deff_gamma": pw["deff_gamma"],
            "mde_gamma": pw["mde_gamma"],
            "gamma_power_recheck": pw["slope_point"],
        }
        for d in ds:
            p = pw["per_d"].get(str(d))
            if not p:
                continue
            rec["power"][str(d)] = {
                "rel_err": p["median_rel_err"],
                "spearman": p["spearman"],
                "n_genes": p["n_genes_reliability"],
                "frac_ge_0p1": p.get("frac_ge_0p1"),
                "mde_1v1": p["mde_pooled_1v1"],
                "mde_2v2": p["mde_pooled_2v2"],
                "alpha_point": p["alpha_pool_point"],
                "n_bins": p["n_domain"],
            }
        # 象限极差
        qs = list(rec["quadrants"].values())
        rec["quad_min"], rec["quad_max"] = min(qs), max(qs)
        rec["quad_range"] = max(qs) - min(qs)
        # 含 MT 的敏感性
        if "alpha_mt" in files:
            rec["gamma_with_mt"] = load(files["alpha_mt"])["slope_all"]
        out["datasets"][tag] = rec

    cv = load(ROOT / "02_OUT/cv_table_3sets.json")
    for tag in out["datasets"]:
        out["datasets"][tag]["cv"] = {d: cv[tag][d]["cv"] for d in cv[tag]}
        out["datasets"][tag]["cv_poisson"] = {d: cv[tag][d]["cv_poisson"] for d in cv[tag]}

    syn = load(ROOT / "02_OUT/_synth_power/power_table_synth_indep.json")
    halves = [v for k, v in syn["slope_halves"].items()]
    vals = [x for pair in halves for x in pair.values()]
    out["synth"]["indep"] = {
        "gamma": syn["slope_point"],
        "rms_cut": syn["rms_gamma_cut"],
        "rms_random": syn["rms_gamma_random"],
        "deff": syn["deff_gamma"],
        "min": min(vals),
        "max": max(vals),
        "mde": syn["mde_gamma"]["n1"],
        "n_halves": len(vals),
    }
    for tag, path in [("tissue", "_synth_tissue_check/gate_results.json"),
                      ("sparse", "_synth_sparse_check/gate_results.json")]:
        p = ROOT / "02_OUT" / path
        d = load(p)
        z = d.get("z_slope_all") or d.get("z_slope")
        out["synth"][tag] = {"n_genes": z["n_genes"], "median_slope": z["median_slope"],
                             "q25": z["q25"], "q75": z["q75"]}
    # 低表达层合成：逐基因 CI 来自 15 号脚本的 bootstrap
    le = load(ROOT / "02_OUT/_synth_lowexp_check/layer_ci_lowexp.json")["layers"]
    le = [r for r in le if r["hi"] == 0.1][0]
    out["synth"]["lowexp"] = {"n_genes": le["n"], "median_slope": le["slope"],
                              "lo": le["ci"][0], "hi_ci": le["ci"][1]}

    lc = load(ROOT / "02_OUT/layer_compare.json")
    out["layers"] = {
        "edges": lc["edges"],
        "ref_d": lc["ref_d"],
        "shared_genes": {k: {"n": v["n"], "slope": v["slope"]} for k, v in lc["shared_genes"].items()},
        "all_genes": {k: {"slope": v["slope"]} for k, v in lc["all_genes"].items()},
        "by_expression": {k: [{"lo": r["lo"], "hi": r["hi"], "slope": r["slope"], "n": r["n"]} for r in v]
                          for k, v in lc["by_expression"].items()},
    }
    out["meta"]["n_datasets"] = 3
    out["meta"]["ds_grid"] = [2, 4, 8, 16, 32, 64, 128, 256]

    Path(args.out).write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("written", args.out)
    for tag, r in out["datasets"].items():
        print(f"{tag:6s} gamma={r['gamma']:.4f} local={r['gamma_local']:.4f} "
              f"ci=[{r['boot_ci'][0]:.4f},{r['boot_ci'][1]:.4f}] quad=[{r['quad_min']:.3f},{r['quad_max']:.3f}]")


if __name__ == "__main__":
    main()

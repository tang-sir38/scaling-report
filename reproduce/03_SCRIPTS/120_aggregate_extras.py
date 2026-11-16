#!/usr/bin/env python3
r"""
120_aggregate_extras.py — 把四个补充组织的实测结果汇总成一份可注入的锁块。

为什么要有这一步
----------------
正文里的每一个倍数都必须来自 `NUMBERS_LOCK_v7.json`，而 97_lock_v7.py 只重建
原来那三片（colon / lung / brain）＋Stereo-seq。新组织的结果散在
`E:\A_spatial_20261116\OUT\<组织>\` 里，谁也不知道最终区间是哪几个数取的。

这个脚本做三件事，全部现场算，不抄任何文档里的数：
  1. 从每个组织的 filtered_feature_bc_matrix.h5 直接读 2 µm 格子数、总 UMI、基因数
     （基因数扣掉线粒体基因，口径与正文 "we removed mitochondrial genes" 一致）；
  2. 从 pooled_alpha_noMT.json 读 α(2 µm) 与 α(256 µm)，相除得下降倍数；
     与 boot_tail.json 的点估计交叉核对，差 >1e-6 就报错；
  3. 与原三片（取自锁）合起来算区间与 floor 取整值，写成
     `02_OUT/extra_tissues_20261116.json`，并用 --write-lock 注入锁的
     `extra_tissues` 块。

用法
----
  python 120_aggregate_extras.py                 # 只算，写独立汇总
  python 120_aggregate_extras.py --write-lock    # 同时写进 NUMBERS_LOCK_v7.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
REPRO = HERE.parent
LOCK = REPRO / "02_OUT" / "NUMBERS_LOCK_v7.json"
SUMMARY = REPRO / "02_OUT" / "extra_tissues_20261116.json"
WS = Path(r"E:\A_spatial_20261116")

# 四个补充组织。label / tissue / condition 会进 Table S1，必须与数据集本身相符。
EXTRAS = {
    "mouse_small_intestine": dict(
        label="Visium HD Mouse Small Intestine", species="Mouse",
        tissue="Small intestine", condition="normal",
        dataset="Visium_HD_Mouse_Small_Intestine", source="10x Genomics (public)",
    ),
    "mouse_kidney": dict(
        label="Visium HD Mouse Kidney", species="Mouse",
        tissue="Kidney", condition="normal",
        dataset="Visium_HD_Mouse_Kidney", source="10x Genomics (public)",
    ),
    "human_pancreas": dict(
        label="Visium HD Human Pancreas", species="Human",
        tissue="Pancreas", condition="normal",
        dataset="Visium_HD_Human_Pancreas", source="10x Genomics (public)",
    ),
    "mouse_embryo": dict(
        label="Visium HD Mouse Embryo", species="Mouse",
        tissue="Embryo", condition="developing",
        dataset="Visium_HD_Mouse_Embryo", source="10x Genomics (public)",
    ),
}

BASE_SECTIONS = ("colon", "lung", "brain")


def find_h5(ext: Path) -> Path:
    h5s = sorted(ext.rglob("square_002um/**/filtered_feature_bc_matrix.h5"))
    if not h5s:
        raise SystemExit(f"找不到 square_002um 的 h5: {ext}")
    return h5s[0]


def counts_from_h5(h5: Path) -> dict:
    """2 µm 格子数、总 UMI、扣 MT 后的基因数。全部从矩阵本身读。"""
    with h5py.File(h5, "r") as f:
        m = f["matrix"]
        shape = tuple(int(x) for x in m["shape"][:])
        total = float(np.asarray(m["data"][:]).sum())
        names = [x.decode("utf-8", "replace") for x in m["features/name"][:]]
    mt = [g for g in names if g.lower().startswith("mt-")]
    return dict(n_bins_2um=shape[1], umi_total=total,
                n_genes_raw=shape[0], n_genes_noMT=shape[0] - len(mt),
                n_mt=len(mt))


def read_out(key: str) -> tuple[dict, dict | None]:
    out = WS / "OUT" / key
    pooled = out / "pooled_alpha_noMT.json"
    if not pooled.exists():
        raise SystemExit(f"缺 {pooled}，先跑 114_pilot_tissue.py --tissue {key}")
    d = json.loads(pooled.read_text(encoding="utf-8"))
    boot = out / "boot_tail.json"
    b = json.loads(boot.read_text(encoding="utf-8")) if boot.exists() else None
    return d, b


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-lock", action="store_true")
    a = ap.parse_args()

    lock_raw = LOCK.read_text(encoding="utf-8")
    lock = json.loads(lock_raw)
    ds = lock["datasets"]

    sections: list[dict] = []

    # 原三片：倍数从锁里算，不读文档
    for k in BASE_SECTIONS:
        al = ds[k]["alpha"]
        a2, a256 = float(al["2"]), float(al["256"])
        sections.append(dict(
            key=k, label={"colon": "Visium HD Human Colon Cancer",
                          "lung": "Visium HD Human Lung Cancer",
                          "brain": "Visium HD Mouse Brain"}[k],
            origin="NUMBERS_LOCK_v7.json :: datasets",
            alpha_2um=a2, alpha_256um=a256, drop_fold=a2 / a256,
            gamma=float(ds[k]["gamma"]),
            n_bins_2um=int(ds[k]["n_bins"]), umi_total=float(ds[k]["umi_total"]),
            umi_per_bin=float(ds[k]["umi_per_bin"]),
            n_genes_noMT=int(ds[k]["n_genes_noMT"]),
            alpha_by_d={str(d): float(v) for d, v in al.items()
                        if int(d) in (2, 4, 8, 16, 32, 64, 128, 256)},
            boot=None,
        ))

    # 四个新组织：全部现场读
    for key, meta in EXTRAS.items():
        h5 = find_h5(WS / "RAW" / key / "extracted")
        c = counts_from_h5(h5)
        pooled, boot = read_out(key)
        a2 = float(pooled["per_d"]["2"]["alpha_pool"])
        a256 = float(pooled["per_d"]["256"]["alpha_pool"])
        drop = a2 / a256
        rec = dict(key=key, label=meta["label"], species=meta["species"],
                   tissue=meta["tissue"], condition=meta["condition"],
                   dataset=meta["dataset"], source=meta["source"],
                   origin=f"E:/A_spatial_20261116/OUT/{key}",
                   alpha_2um=a2, alpha_256um=a256, drop_fold=drop,
                   gamma=float(pooled["slope_all"]),
                   n_bins_2um=c["n_bins_2um"], umi_total=c["umi_total"],
                   umi_per_bin=c["umi_total"] / c["n_bins_2um"],
                   n_genes_raw=c["n_genes_raw"], n_genes_noMT=c["n_genes_noMT"],
                   n_mt=c["n_mt"],
                   h5=str(h5),
                   alpha_by_d={str(d): float(pooled["per_d"][str(d)]["alpha_pool"])
                               for d in (2, 4, 8, 16, 32, 64, 128, 256)},
                   boot=None)
        if boot is not None:
            bd = boot["drop_fold"]
            assert abs(bd["point"] - drop) < 1e-6, (
                f"{key}: 与 bootstrap 点估计不一致 {bd['point']:.6f} vs {drop:.6f}")
            rec["boot"] = dict(n=int(boot["n_boot"]), seed=int(boot["seed"]),
                               point=bd["point"], ci95=bd["ci_from_alpha256"],
                               alpha_256_ci95=boot["per_d"]["256"]["ci95"],
                               n_cells_256=int(boot["per_d"]["256"]["n_cells"]))
        sections.append(rec)

    drops = [s["drop_fold"] for s in sections]
    lo_t, hi_t = min(drops), max(drops)
    n_sec = len(sections)

    # floor 取整（标题惯例：5.82→5、14.31→14，实测已验证）
    floor_lo, floor_hi = math.floor(lo_t), math.floor(hi_t)

    summary = dict(
        generated_by="120_aggregate_extras.py",
        lock_sha256_before=hashlib.sha256(lock_raw.encode("utf-8")).hexdigest(),
        n_sections=n_sec,
        sections=sections,
        range=dict(min=lo_t, max=hi_t, min_key=min(sections, key=lambda s: s["drop_fold"])["key"],
                   max_key=max(sections, key=lambda s: s["drop_fold"])["key"],
                   floor_min=floor_lo, floor_max=floor_hi,
                   title_slice=f"{floor_lo}- to {floor_hi}-fold"),
        note=("drop_fold = alpha(2um)/alpha(256um)，逐片现场算；"
              "floor 取整是标题惯例，不是四舍五入"),
    )

    if a.write_lock:
        lock["extra_tissues"] = summary
        LOCK.write_text(json.dumps(lock, ensure_ascii=False, indent=1), encoding="utf-8")
        print(f"写入 {LOCK}")

    SUMMARY.write_text(json.dumps(summary, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"写入 {SUMMARY}\n")

    print(f"{'section':26s} {'a(2um)':>8s} {'a(256)':>8s} {'drop':>7s} {'floor':>5s}  CI")
    for s in sorted(sections, key=lambda x: x["drop_fold"]):
        ci = ""
        if s.get("boot"):
            ci = f"[{s['boot']['ci95'][0]:.2f}, {s['boot']['ci95'][1]:.2f}]"
        print(f"{s['label'][11:]:26s} {s['alpha_2um']:8.3f} {s['alpha_256um']:8.3f} "
              f"{s['drop_fold']:7.2f} {math.floor(s['drop_fold']):5d}  {ci}")
    print(f"\n区间 {lo_t:.2f}-{hi_t:.2f}，floor → {floor_lo}- to {floor_hi}-fold")
    print(f"片数 {n_sec}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

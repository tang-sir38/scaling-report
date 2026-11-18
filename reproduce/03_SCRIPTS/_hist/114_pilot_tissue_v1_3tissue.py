#!/usr/bin/env python3
r"""
114_pilot_tissue.py — 新组织试点：体检 → 重聚合 → 池化指数 → 判读。

为什么先试点一片
----------------
现有三个组织的 alpha 从 2 um 到 256 um 分别下降 14.31（结肠）、13.76（肺）、
5.82（脑）倍，标题的 "5- to 14-fold" 就是这三个数取的区间。加新组织之前必须
先知道它的下降倍数落不落在 5.82-14.31 里 —— 落在里面，标题不用动；落在外面，
标题、摘要、正文至少三处要改。所以先跑一片，拿到数再决定扩不扩。

这个脚本不做任何新分析，只是把 02b_rebin.py 和 07_pooled_alpha.py 串起来，
中间加一道内存体检和一道判读，免得跑完不知道怎么读。

用法
----
  python 114_pilot_tissue.py --tissue mouse_small_intestine --check-only
  python 114_pilot_tissue.py --tissue mouse_small_intestine
  python 114_pilot_tissue.py --tissue mouse_small_intestine --skip-rebin   # 只重算指数

前置
----
  113_fetch_tissue.py 已下完并解包，即
  <root>/RAW/<tissue>/extracted/ 下有 binned_outputs/square_002um/
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import h5py

HERE = Path(__file__).resolve().parent
ROOT_DEFAULT = Path(r"E:\A_spatial_20261116")
REDUCE = [2, 4, 8, 16, 32, 64, 128, 256]

# 现有三个组织的下降倍数（2 um -> 256 um），取自 02_OUT/NUMBERS_LOCK_v7.json
KNOWN = {"colon": 14.31, "lung": 13.76, "brain": 5.82}
LO, HI = min(KNOWN.values()), max(KNOWN.values())


def find_inputs(extracted: Path) -> tuple[Path, Path]:
    """在解包目录里找 2 um 的 h5 与 tissue_positions.parquet。"""
    h5s = sorted(extracted.rglob("square_002um/**/filtered_feature_bc_matrix.h5"))
    if not h5s:
        h5s = sorted(extracted.rglob("**/filtered_feature_bc_matrix.h5"))
        h5s = [p for p in h5s if "square_002um" in str(p).replace("\\", "/")]
    if not h5s:
        raise SystemExit(f"在 {extracted} 下找不到 square_002um 的 filtered_feature_bc_matrix.h5")
    h5 = h5s[0]
    prqs = sorted(h5.parent.rglob("tissue_positions.parquet"))
    if not prqs:
        prqs = sorted(extracted.rglob("**/tissue_positions.parquet"))
    if not prqs:
        raise SystemExit(f"在 {extracted} 下找不到 tissue_positions.parquet")
    return h5, prqs[0]


def check(h5: Path) -> dict:
    """只读元数据，不加载矩阵。用来估算峰值内存。"""
    with h5py.File(h5, "r") as f:
        m = f["matrix"]
        shape = tuple(int(x) for x in m["shape"][:])
        nnz = int(m["data"].shape[0])
        dtype = m["data"].dtype
        genes = len(m["features/name"])
    item = dtype.itemsize
    # 02b_rebin 的峰值：原始 data+indices+indptr 一份，CSC->CSR 转换再一份
    raw = nnz * (item + 4)
    peak = raw * 2.6
    return dict(shape=shape, genes=genes, nnz=nnz, dtype=str(dtype),
                raw_gb=raw / 2**30, peak_gb=peak / 2**30)


def free_gb() -> float:
    out = subprocess.run(
        ["powershell", "-NoProfile", "-Command",
         "(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory"],
        capture_output=True, text=True)
    try:
        return int(out.stdout.strip()) / 1024 / 1024
    except Exception:
        return -1.0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue", required=True)
    ap.add_argument("--root", default=str(ROOT_DEFAULT))
    ap.add_argument("--check-only", action="store_true")
    ap.add_argument("--skip-rebin", action="store_true")
    ap.add_argument("--force", action="store_true", help="内存不足也跑")
    a = ap.parse_args()

    root = Path(a.root)
    ext = root / "RAW" / a.tissue / "extracted"
    if not ext.exists():
        raise SystemExit(f"未解包: {ext}\n先跑 113_fetch_tissue.py --tissue {a.tissue} --extract")
    h5, prq = find_inputs(ext)
    print(f"h5      {h5}")
    print(f"parquet {prq}")

    c = check(h5)
    print(f"矩阵    {c['genes']:,} genes x {c['shape'][1]:,} bins  nnz={c['nnz']:,}  {c['dtype']}")
    print(f"预估峰值内存 {c['peak_gb']:.1f} GiB（原始读入 {c['raw_gb']:.1f}）")
    fg = free_gb()
    print(f"当前可用内存 {fg:.1f} GiB")
    verdict = "ok"
    if fg > 0 and fg < c["peak_gb"] * 1.15:
        verdict = "tight"
        print(f"  !! 可用内存不足峰值需求（{fg:.1f} < {c['peak_gb'] * 1.15:.1f}）")
        print("     先关掉浏览器 / Word / 多余 WorkBuddy 窗口，再重跑；")
        print("     或加 --force 硬跑（有 OOM 风险，进程会被系统杀掉）")
    if a.check_only:
        print(json.dumps(dict(inputs=dict(h5=str(h5), parquet=str(prq)),
                              check=c, free_gb=fg, verdict=verdict),
                         ensure_ascii=False, indent=1))
        return 0
    if verdict == "tight" and not a.force:
        raise SystemExit("内存不足，已中止。释放内存后重跑，或加 --force")

    work = root / "WORK" / a.tissue
    out = root / "OUT" / a.tissue
    work.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)

    if not a.skip_rebin:
        cmd = [sys.executable, str(HERE / "02b_rebin.py"),
               "--h5", str(h5), "--parquet", str(prq),
               "--outdir", str(work), "--d"] + [str(d) for d in REDUCE]
        print("\n重聚合：", " ".join(cmd))
        r = subprocess.run(cmd)
        if r.returncode != 0:
            raise SystemExit(f"02b_rebin 失败 rc={r.returncode}")

    cmd = [sys.executable, str(HERE / "07_pooled_alpha.py"),
           "--workdir", str(work), "--outdir", str(out),
           "--ds"] + [str(d) for d in REDUCE] + ["--exclude-mt"]
    print("\n池化指数：", " ".join(cmd))
    r = subprocess.run(cmd)
    if r.returncode != 0:
        raise SystemExit(f"07_pooled_alpha 失败 rc={r.returncode}")

    res = json.loads((out / "pooled_alpha_noMT.json").read_text(encoding="utf-8"))
    per_d = res["per_d"]
    missing = [d for d in ("2", "256") if d not in per_d]
    if missing:
        raise SystemExit(f"缺分辨率档 {missing}；实际有 {sorted(per_d, key=int)}")
    a2 = float(per_d["2"]["alpha_pool"])
    a256 = float(per_d["256"]["alpha_pool"])
    drop = a2 / a256
    gamma = res.get("slope_all")

    inside = LO <= drop <= HI
    print("\n" + "=" * 66)
    print(f"{a.tissue}: alpha(2um)={a2:.3f}  alpha(256um)={a256:.3f}  下降 {drop:.2f} 倍")
    print(f"gamma = {gamma}")
    print(f"现有区间 {LO:.2f}-{HI:.2f}（脑 5.82 / 肺 13.76 / 结肠 14.31）")
    print("判读：" + ("落在区间内 —— 标题 \"5- to 14-fold\" 不用改"
                     if inside else "落在区间外 —— 标题的倍数区间必须改"))
    print("=" * 66)

    (out / "pilot_verdict.json").write_text(json.dumps(dict(
        tissue=a.tissue, alpha_2um=a2, alpha_256um=a256, drop_fold=drop,
        gamma=gamma, existing_range=[LO, HI], inside_range=inside,
        known=KNOWN,
    ), ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"写入 {out / 'pilot_verdict.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

# -*- coding: utf-8 -*-
"""301_p5_extract.py — 从 Visium HD P5 的 binned_outputs.tar.gz 里只取需要的两个成员。

tar.gz 不能随机访问，必须顺序读完整个流；但只把 square_002um 的
filtered_feature_bc_matrix.h5 与 spatial/tissue_positions.parquet 落到磁盘，
其余（raw_probe、图像、8/16 µm）全部丢弃，省空间。

用法：
  python 301_p5_extract.py            # 解压
  python 301_p5_extract.py --dry-run  # 只列成员清单
"""
from __future__ import annotations

import argparse
import os
import sys
import tarfile
import time
from pathlib import Path

TAR = Path(r"E:\visium_p5\Visium_HD_Human_Colon_Cancer_P5_binned_outputs.tar.gz")
TOTAL = 16568446864
DEST = Path(r"E:\A_spatial_20261116\RAW\human_colon_cancer_p5\extracted")

WANT_SUFFIX = (
    "square_002um/filtered_feature_bc_matrix.h5",
    "square_002um/spatial/tissue_positions.parquet",
)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    if not TAR.exists():
        raise SystemExit(f"还没下完: {TAR}")
    size = TAR.stat().st_size
    print(f"tar {size:,} B  ({size/1e9:.2f} GB)  expect {TOTAL:,}"
          f"  {'OK' if size == TOTAL else 'SIZE MISMATCH'}", flush=True)

    DEST.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    got, seen = [], 0
    with tarfile.open(TAR, "r|gz") as tf:
        for m in tf:
            seen += 1
            name = m.name.replace("\\", "/")
            if a.dry_run:
                if seen <= 12 or any(name.endswith(s) for s in WANT_SUFFIX):
                    print(f"  {m.size:>13,}  {name}")
                continue
            if not any(name.endswith(s) for s in WANT_SUFFIX):
                continue
            rel = name.split("binned_outputs/", 1)[-1]
            out = DEST / "binned_outputs" / rel
            out.parent.mkdir(parents=True, exist_ok=True)
            f = tf.extractfile(m)
            if f is None:
                continue
            with open(out, "wb") as fo:
                while True:
                    c = f.read(1 << 22)
                    if not c:
                        break
                    fo.write(c)
            got.append((str(out.relative_to(DEST)), out.stat().st_size))
            print(f"  extracted {out.name}  {out.stat().st_size:,} B  "
                  f"{time.time()-t0:.0f}s", flush=True)

    print(f"\n{seen} members scanned in {time.time()-t0:.0f}s")
    for p, n in got:
        print(f"  {p}  {n:,} B")
    if not a.dry_run and len(got) != len(WANT_SUFFIX):
        raise SystemExit(f"只拿到 {len(got)}/{len(WANT_SUFFIX)} 个成员")
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""Stage the v5.5.1 package directory from the delivered v5.5 package.

v5.5 is read only. v5.5.1 differs in exactly six figure files (Fig1, Fig5,
FigS3; pdf + png each), two docx (manuscript, supplementary) and their two
exported PDFs, plus the audit trail. Everything else is carried forward byte
for byte — re-deriving figures nobody asked about would move numbers silently.

Usage: python 117_stage_v551.py [--force]
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

SRC = Path(r"C:\Users\唐朝华\Desktop\23空间分辨率\A题_空间分辨率_投稿包_v5.5_20261116")
DST = Path(r"C:\Users\唐朝华\Desktop\23空间分辨率\A题_空间分辨率_投稿包_v5.5.1_20261117")
BUILD = Path(r"C:\Users\唐朝华\Desktop\scaling-report\reproduce\07_SUBMISSION\_build_v551")
FIGS = BUILD / "figs"


def main() -> int:
    if not SRC.exists():
        raise SystemExit(f"v5.5 包不存在: {SRC}")
    if DST.exists():
        if "--force" in sys.argv:
            shutil.rmtree(DST)
        else:
            raise SystemExit(f"目标已存在，加 --force 重建: {DST}")
    DST.mkdir(parents=True, exist_ok=True)

    n = 0
    for p in sorted(SRC.rglob("*")):
        if not p.is_file():
            continue
        target = DST / p.relative_to(SRC)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, target)
        n += 1
    print(f"carried forward {n} files  {SRC.name} -> {DST.name}")

    # 六张图换成 v5.5.1 的新版；只要有一张缺就停，不做部分替换
    need = ["Fig1.pdf", "Fig1.png", "Fig5.pdf", "Fig5.png",
            "FigS3.pdf", "FigS3.png"]
    missing = [f for f in need if not (FIGS / f).exists()]
    if missing:
        raise SystemExit(f"新图未生成: {missing}")
    for f in need:
        shutil.copy2(FIGS / f, DST / "04_图件" / f)
    print(f"replaced {len(need)} figure files from {FIGS}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

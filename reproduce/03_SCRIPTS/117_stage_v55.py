"""Stage the v5.5 package directory from the delivered v5.4 package.

v5.5 differs from v5.4 in three docx, the supplementary table (four new rows),
one new figure, and the numbers lock (which gains an `extra_tissues` block).
Everything else — the six main figures, the older audit trail, the derived
tables — is carried forward byte for byte, because re-deriving them was not part
of this change and a rebuild would silently move numbers nobody asked about.

The v5.4 directory is only read. Nothing in it is renamed, moved or deleted.
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

SRC = Path(r"C:\Users\唐朝华\Desktop\23空间分辨率\A题_空间分辨率_投稿包_v5.4_20261116")
DST = Path(r"C:\Users\唐朝华\Desktop\23空间分辨率\A题_空间分辨率_投稿包_v5.5_20261116")
LOCK = Path(r"C:\Users\唐朝华\Desktop\scaling-report\reproduce\02_OUT\NUMBERS_LOCK_v7.json")


def main() -> int:
    if not SRC.exists():
        raise SystemExit(f"v5.4 包不存在: {SRC}")
    DST.mkdir(parents=True, exist_ok=True)

    n = 0
    for p in sorted(SRC.rglob("*")):
        if not p.is_file():
            continue
        rel = p.relative_to(SRC)
        target = DST / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, target)
        n += 1
    print(f"carried forward {n} files  {SRC.name} -> {DST.name}")

    lock = json.loads(LOCK.read_text(encoding="utf-8"))
    if "extra_tissues" not in lock:
        raise SystemExit("锁里还没有 extra_tissues，先跑 120_aggregate_extras.py --write-lock")
    (DST / "05_Audit").mkdir(parents=True, exist_ok=True)
    (DST / "05_Audit" / "NUMBERS_LOCK_v7.json").write_bytes(LOCK.read_bytes())
    print(f"lock copied (top-level keys {len(lock)}), extra_tissues present: "
          f"{len(lock['extra_tissues']['sections'])} sections")
    return 0


if __name__ == "__main__":
    sys.exit(main())

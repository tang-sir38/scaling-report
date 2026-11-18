"""Rebuild MANIFEST.json for the v5.5.1 package so that it lists every file with
its byte count and sha256.

The manifest has to be rebuilt after the guard has written its own output into
05_Audit, because the guard checks that the manifest covers every file present
— including the guard's own report. The guard's JSON carries no timestamp, so
building the manifest and re-running the guard converges instead of looping.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

PKG = Path(r"C:\Users\唐朝华\Desktop\23空间分辨率\A题_空间分辨率_投稿包_v5.5.1_20261117")
NAME = PKG.name


def main() -> int:
    files = []
    for p in sorted(PKG.rglob("*")):
        if not p.is_file():
            continue
        if p.name == "MANIFEST.json":
            continue
        if "__pycache__" in p.parts or p.suffix in {".pyc", ".pyo"}:
            continue
        files.append(dict(
            file=p.relative_to(PKG).as_posix(),
            bytes=p.stat().st_size,
            sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
        ))
    rec = dict(package=NAME, files=files)
    (PKG / "MANIFEST.json").write_text(
        json.dumps(rec, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"written {PKG / 'MANIFEST.json'}  {len(files)} files")
    return 0


if __name__ == "__main__":
    sys.exit(main())

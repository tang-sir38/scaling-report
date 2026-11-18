"""112_pack_v56.py — 重建 MANIFEST 并把 v5.6 打成交付 zip。

规则与 v5.5.1 一致：一个顶层目录、每个 hash 现场重算、.sha256 sidecar、
项目目录与桌面两份逐字节比对。v5.4 / v5.5 / v5.5.1 的包一律不动。
"""
from __future__ import annotations

import hashlib
import json
import sys
import zipfile
from pathlib import Path

SRC = Path(r"C:\Users\唐朝华\Desktop\23空间分辨率\A题_空间分辨率_投稿包_v5.6_20261117")
NAME = "A题_空间分辨率_投稿包_v5.6_20261117"
ZIP_NAME = "－－23－－空间分辨率_投稿包_v5.6_20261117.zip"
DEST_DIRS = [
    Path(r"C:\Users\唐朝华\Desktop\23空间分辨率"),
    Path(r"C:\Users\唐朝华\Desktop"),
]
BUILD = Path(r"C:\Users\唐朝华\Desktop\scaling-report\reproduce\07_SUBMISSION"
             r"\_build_v56\audit")


def collect() -> list[Path]:
    out = []
    for p in sorted(SRC.rglob("*")):
        if not p.is_file():
            continue
        if any(part in {"__pycache__", ".pytest_cache"} for part in p.parts):
            continue
        if p.suffix in {".pyc", ".pyo"}:
            continue
        if p.name == "MANIFEST.json":
            continue
        out.append(p)
    return out


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        while True:
            c = f.read(1 << 20)
            if not c:
                break
            h.update(c)
    return h.hexdigest()


def write_manifest(files: list[Path]) -> dict:
    rec = dict(package=NAME,
               files=[dict(file=p.relative_to(SRC).as_posix(),
                           bytes=p.stat().st_size, sha256=sha(p))
                      for p in files])
    (SRC / "MANIFEST.json").write_text(
        json.dumps(rec, ensure_ascii=False, indent=1), encoding="utf-8")
    return rec


def write_zip(target: Path, files: list[Path]) -> str:
    if target.exists():
        target.unlink()
    with zipfile.ZipFile(target, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for f in files:
            arc = f"{NAME}/{f.relative_to(SRC).as_posix()}"
            zi = zipfile.ZipInfo(arc, date_time=(2026, 11, 17, 0, 0, 0))
            zi.compress_type = zipfile.ZIP_DEFLATED
            zi.external_attr = (0o644 & 0xFFFF) << 16
            z.writestr(zi, f.read_bytes())
    return hashlib.sha256(target.read_bytes()).hexdigest()


def main() -> int:
    files = collect()
    man = write_manifest(files)
    all_files = files + [SRC / "MANIFEST.json"]

    print(f"source {SRC}")
    print(f"  {len(all_files)} files, {sum(p.stat().st_size for p in all_files):,} bytes")

    shas = {}
    for d in DEST_DIRS:
        target = d / ZIP_NAME
        s = write_zip(target, all_files)
        shas[str(target)] = s
        (d / (ZIP_NAME + ".sha256")).write_text(f"{s}  {ZIP_NAME}\n", encoding="utf-8")
        print(f"  wrote {target}  {target.stat().st_size:,} bytes  sha256 {s[:16]}")

    uniq = set(shas.values())
    assert len(uniq) == 1, f"两个 zip 不同: {shas}"

    with zipfile.ZipFile(DEST_DIRS[0] / ZIP_NAME) as z:
        names = [n for n in z.namelist() if not n.endswith("/")]
        tops = {n.split("/")[0] for n in names}
        assert tops == {NAME}, f"顶层: {tops}"
        assert len(names) == len(all_files), f"{len(names)} vs {len(all_files)}"
        bad = []
        for e in man["files"]:
            arc = f"{NAME}/{e['file']}"
            if arc not in z.namelist():
                bad.append(("missing", e["file"]))
                continue
            data = z.read(arc)
            if len(data) != e["bytes"] or \
                    hashlib.sha256(data).hexdigest() != e["sha256"]:
                bad.append(("mismatch", e["file"]))
        assert not bad, f"与 manifest 不符: {bad[:5]}"
        print(f"  回读: {len(names)} 项, 单顶层, {len(man['files'])} 个 hash 全对")

    rec = dict(zip=ZIP_NAME, sha256=list(uniq)[0],
               bytes=(DEST_DIRS[0] / ZIP_NAME).stat().st_size,
               entries=len(all_files), top_level=NAME,
               copies=[str(d / ZIP_NAME) for d in DEST_DIRS],
               manifest_entries=len(man["files"]),
               verification="zip 回读；每项对 MANIFEST.json 核 hash；两份逐字节相同",
               previous="v5.5.1 原包保留：－－23－－空间分辨率_投稿包_v5.5.1_20261117.zip")
    BUILD.mkdir(parents=True, exist_ok=True)
    (BUILD / "pack_v5_6.json").write_text(
        json.dumps(rec, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"  两份一致: {list(uniq)[0]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""Pack v5.4 into the delivery zip.

Follows the same rules as 107_pack_v5.py / 109_pack_v52.py: a single top-level
directory inside the archive, no build residue, every hash recomputed here, and
a .sha256 sidecar written next to the archive. The zip is written twice — once
to the project directory and once to the Desktop — and the two are compared
byte for byte, because a delivery that differs between the place it is built
and the place it is picked up is the one failure mode that no downstream check
would catch.
"""
from __future__ import annotations

import hashlib
import json
import zipfile
from pathlib import Path

SRC = Path(r"C:\Users\唐朝华\Desktop\23空间分辨率\A题_空间分辨率_投稿包_v5.4_20261116")
NAME = "A题_空间分辨率_投稿包_v5.4_20261116"
ZIP_NAME = "－－23－－空间分辨率_投稿包_v5.4_20261116.zip"
DEST_DIRS = [
    Path(r"C:\Users\唐朝华\Desktop\23空间分辨率"),
    Path(r"C:\Users\唐朝华\Desktop"),
]
SKIP_DIRS = {"__pycache__", ".pytest_cache"}
SKIP_SUFFIX = {".pyc", ".pyo"}


def collect() -> list[Path]:
    out = []
    for p in sorted(SRC.rglob("*")):
        if not p.is_file():
            continue
        if any(part in SKIP_DIRS for part in p.parts):
            continue
        if p.suffix in SKIP_SUFFIX:
            continue
        out.append(p)
    return out


def write_zip(target: Path, files: list[Path]) -> str:
    if target.exists():
        target.unlink()
    with zipfile.ZipFile(target, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for f in files:
            arc = f"{NAME}/{f.relative_to(SRC).as_posix()}"
            zi = zipfile.ZipInfo(arc, date_time=(2026, 11, 16, 0, 0, 0))
            zi.compress_type = zipfile.ZIP_DEFLATED
            zi.external_attr = (0o644 & 0xFFFF) << 16
            z.writestr(zi, f.read_bytes())
    return hashlib.sha256(target.read_bytes()).hexdigest()


def main() -> None:
    files = collect()
    man = json.loads((SRC / "MANIFEST.json").read_text(encoding="utf-8"))
    listed = {e["file"] for e in man["files"]}
    actual = {p.relative_to(SRC).as_posix() for p in files} - {"MANIFEST.json"}
    assert listed == actual, f"manifest mismatch: {sorted(listed ^ actual)[:5]}"

    print(f"source {SRC}")
    print(f"  {len(files)} files, {sum(p.stat().st_size for p in files):,} bytes")

    shas = {}
    for d in DEST_DIRS:
        target = d / ZIP_NAME
        sha = write_zip(target, files)
        shas[str(target)] = sha
        size = target.stat().st_size
        print(f"  wrote {target}  {size:,} bytes  sha256 {sha[:16]}")
        (d / (ZIP_NAME + ".sha256")).write_text(
            f"{sha}  {ZIP_NAME}\n", encoding="utf-8"
        )

    uniq = set(shas.values())
    assert len(uniq) == 1, f"the two archives differ: {shas}"
    print(f"  two archives identical: {list(uniq)[0]}")

    # read the archive back and verify every entry against the manifest
    with zipfile.ZipFile(DEST_DIRS[0] / ZIP_NAME) as z:
        names = [n for n in z.namelist() if not n.endswith("/")]
        tops = {n.split("/")[0] for n in names}
        assert tops == {NAME}, f"top-level entries: {tops}"
        assert len(names) == len(files), f"{len(names)} entries vs {len(files)} files"
        bad = []
        for e in man["files"]:
            arc = f"{NAME}/{e['file']}"
            if arc not in z.namelist():
                bad.append(("missing", e["file"]))
                continue
            data = z.read(arc)
            if len(data) != e["bytes"] or hashlib.sha256(data).hexdigest() != e["sha256"]:
                bad.append(("mismatch", e["file"]))
        assert not bad, f"archive does not match manifest: {bad[:5]}"
        print(f"  archive re-read: {len(names)} entries, single top level, "
              f"all {len(man['files'])} hashes match")

    rec = dict(
        zip=ZIP_NAME,
        sha256=list(uniq)[0],
        bytes=(DEST_DIRS[0] / ZIP_NAME).stat().st_size,
        entries=len(files),
        top_level=NAME,
        copies=[str(d / ZIP_NAME) for d in DEST_DIRS],
        manifest_entries=len(man["files"]),
        verification="archive re-read; every entry hashed against MANIFEST.json; two copies byte-identical",
    )
    out = Path(r"C:\Users\唐朝华\Desktop\scaling-report\reproduce\07_SUBMISSION\_build_v54\audit")
    out.mkdir(parents=True, exist_ok=True)
    (out / "pack_v5_4.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()

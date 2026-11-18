r"""133_figembed_v58.py — 把主稿/SI 里内嵌的 Fig1、Fig S3 换成九片版。

规则（v5.7 踩过的坑）
--------------------
- 按**字节**认旧图，不按顺序认：顺序一变就换错图。
- 换 blob 的同时按新图的像素比同步 extent 高度，宽度保持不变。
- 只换内嵌副本不算完，还必须重新导 PDF（132_export_pdf_v58.py）。
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import docx
from docx.shared import Emu
from PIL import Image

NEW = Path(r"C:\Users\唐朝华\Desktop\23空间分辨率\A题_空间分辨率_投稿包_v5.8_20261118")
OLD = Path(r"C:\Users\唐朝华\Desktop\23空间分辨率\A题_空间分辨率_投稿包_v5.7_20261118")

PKG = NEW
MAN = PKG / "02_主稿" / "01_Manuscript.docx"
SI = PKG / "03_补充材料" / "02_Supplementary_Information.docx"


def swap(doc_path: Path, old_png: Path, new_png: Path) -> dict:
    d = docx.Document(str(doc_path))
    old_sha = hashlib.sha256(old_png.read_bytes()).hexdigest()
    blob = new_png.read_bytes()
    target = None
    for part in d.part.package.iter_parts():
        try:
            b = part.blob
        except Exception:
            continue
        if part.content_type.startswith("image/") and hashlib.sha256(b).hexdigest() == old_sha:
            target = part
            break
    if target is None:
        raise SystemExit(f"{doc_path.name}: 找不到 {old_png.name} 的旧字节")
    shape = None
    for sh in d.inline_shapes:
        rid = sh._inline.graphic.graphicData.pic.blipFill.blip.embed
        if d.part.related_parts[rid] is target:
            shape = sh
            break
    if shape is None:
        raise SystemExit(f"{doc_path.name}: 找不到引用该图的 inline shape")

    w, h = Image.open(new_png).size
    width = shape.width
    target._blob = blob
    shape.width = Emu(width)
    shape.height = Emu(int(round(width * h / w)))
    d.save(str(doc_path))
    return dict(png=[w, h], width_mm=round(width / 36000, 1),
                height_mm=round(shape.height / 36000, 1),
                new_sha=hashlib.sha256(blob).hexdigest()[:16])


def main() -> int:
    r1 = swap(MAN, OLD / "04_图件" / "Fig1.png", NEW / "04_图件" / "Fig1.png")
    print("  MAN Fig1 ", r1)
    r2 = swap(SI, OLD / "04_图件" / "FigS3.png", NEW / "04_图件" / "FigS3.png")
    print("  SI  FigS3", r2)

    # 回读：包内每个 image part 的 sha 必须与 04_图件 的对应文件一致
    for doc_path, name in ((MAN, "Fig1"), (SI, "FigS3")):
        want = hashlib.sha256((NEW / "04_图件" / f"{name}.png").read_bytes()).hexdigest()
        d = docx.Document(str(doc_path))
        got = {hashlib.sha256(p.blob).hexdigest() for p in d.part.package.iter_parts()
               if p.content_type.startswith("image/")}
        print(f"  {'OK ' if want in got else 'FAIL'} {doc_path.name} 内嵌 {name} 与 04_图件 同字节")
        assert want in got
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

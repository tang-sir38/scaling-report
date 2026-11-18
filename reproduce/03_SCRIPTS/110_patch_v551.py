r"""110_patch_v551.py — v5.5 → v5.5.1 的全部文字改动，一处不动 docx 结构。

五个洞，逐条对应：

1. Fig 1b 七条曲线 + 端点标倍数 —— 图由 123_fig1_v551.py 画；这里只改图注，
   让图注说清"数在末端"这件事（正文标题的 2- to 15-fold 自此有图上落点）。
2. Fig 1a 几何 —— 同上，图注改写成 8×8 个 2 µm 原生格、16 µm 整块、
   16,384 个原生格构成一个 256 µm bin；三个互相矛盾的标尺不再出现。
3. Fig 5a 的"95% range" —— 图由 124_fig5_v551.py 画（range 画出端帽）；
   图注删掉"horizontal bar spans the 95% range"，改成连线是什么、range 有多宽。
4. 0.818 vs Table S3 0.839 —— 统一到 0.839。理由：0.818 来自 88_region_excise.py
   把 isolate 曲线代进 pooled 公式时**沿用了操作前的 μ² 权重**；isolate 会改变
   IGKC 的 μ（2 µm 处 0.728→0.264），权重必须跟着改。改成正确权重后 pooled
   isolate = 0.723，而"矩阵里只剩 IGKC 一个基因"的 pooled 值按定义恒等于该基因
   自己的 γ = 0.839，即 Table S3 的 "γ, isolated"。正文改成后者，与表一致。
   同一 bug 使 excise 的 0.502 应为 0.494 —— 一并改（125_m1_numbers.py 复算，
   两个旧值都在该 json 的 gamma_frozen_weight 里留档）。
5. P1 复现数 —— Methods 早就写了 "The second colon section is a replicate"，
   但 Results 从来没报它的 fold。补一句（13.1 vs 14.3、Δγ=0.006、z=0.20），
   并从别处删等量词。

另外：Fig S3 因 Fig 1b 已画全七片而成了重复图，改画绝对值的 α(d)
（126_fig_s3_v551.py），图注随之改写；SI 内嵌的那张 image3.png 一并替换。

用法：
  python 110_patch_v551.py --dry-run     # 只打印 diff 与词数，不写
  python 110_patch_v551.py               # 写入
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import docx
from docx.shared import Emu

PKG = Path(r"C:\Users\唐朝华\Desktop\23空间分辨率\A题_空间分辨率_投稿包_v5.5.1_20261117")
BUILD = Path(r"C:\Users\唐朝华\Desktop\scaling-report\reproduce\07_SUBMISSION\_build_v551")
OUT = BUILD / "audit" / "patch_v5_5_1.json"

MAN = PKG / "02_主稿" / "01_Manuscript.docx"
SI = PKG / "03_补充材料" / "02_Supplementary_Information.docx"

# ---------------------------------------------------------------- 正文改动
# (段号, 原文片段, 新文片段, 说明)
BODY = [
    (41,
     "the pooled exponent falls to 0.502, and taking IGKC as the only gene gives 0.818.",
     "the pooled exponent falls to 0.494, and keeping only those regions gives IGKC itself 0.839.",
     "0.502→0.494（权重随操作更新）；0.818→0.839（=Table S3 的 γ, isolated，即只剩该基因时的 pooled 值）"),
    (46,
     "A second colon section, captured separately, gives 0.549 ± 0.026 against 0.554;",
     "A second section of the same block, captured separately, falls 13.1-fold against "
     "14.3-fold, with γ = 0.549 ± 0.026 against 0.554, a difference of 0.006 "
     "(block bootstrap z = 0.20);",
     "把 P1 的 fold 写进 Results，兑现 Methods 的 replicate"),
    (56,
     "What this design still cannot do is test the prediction where it would bite, "
     "below 8 µm, which is a measurement we cannot make.",
     "What this design still cannot do is test the prediction where it would bite, below 8 µm.",
     "腾词：删掉与末段重复的 'which is a measurement we cannot make'"),
    (63,
     "changes the pooled exponent without changing the tissue, so the new value belongs in the result.",
     "changes the pooled exponent without changing the tissue.",
     "腾词：'the new value belongs in the result' 已由上一段（binned result should carry…）承担"),
    (57,
     "geostatistics has treated the variance of a block as an integral of the covariance for "
     "decades, ecology has measured variance against grain for as long,",
     "geostatistics has treated block variance as a covariance integral for decades, ecology "
     "has measured variance against grain as long,",
     "腾词：同一句压缩 6 词，论点不变"),
    (60,
     "so it establishes that arrangement matters, not what a particular arrangement would do.",
     "so it establishes that arrangement matters, not which arrangement does.",
     "腾词：3 词"),
]

# ---------------------------------------------------------------- 图注改动
CAPTIONS = [
    (154, 1,
     " (a) Exact integer re-aggregation: native 2 µm bins are summed into 4, 8, 16, 32, "
     "64, 128 and 256 µm blocks on the same grid, without interpolation and without "
     "discarding reads. (b) Pooled overdispersion relative to its value at 2 µm, against "
     "bin size, for colon (dark blue), lung (mid blue) and mouse brain (ochre), three of "
     "the seven sections; all seven are in Fig. S3. Dashed grey, α ∝ d",
     " (a) Exact integer re-aggregation. One square is one native 2 µm bin and the "
     "outlined block is 8 × 8 of them, 16 µm on a side; a 256 µm bin holds 128 × 128 "
     "native bins, 16,384, summed as integers with no interpolation and no read "
     "discarded. (b) Pooled overdispersion relative to its value at 2 µm, against bin "
     "size, for all seven sections; the number at the right end of each curve is the "
     "fall between 2 and 256 µm, from 15.6-fold (small intestine) to 2.1-fold "
     "(pancreas). The same curves in absolute units are in Fig. S3. Dashed grey, α ∝ d",
     "Fig 1 图注：七条曲线 + 端点倍数；panel a 几何改写"),
    (154, 3,
     ", the counting-statistics limit that merging is meant to approach; dotted grey, "
     "α ∝ d⁰. If the working assumption of the field held, the three curves would lie "
     "on the steep grey line.",
     ", the counting-statistics limit that merging is meant to approach; dotted grey, "
     "α ∝ d⁰. If the working assumption of the field held, the seven curves would lie "
     "on the steep grey line.",
     "三→七"),
    (162, 1,
     "; the horizontal bar spans the 95% range over permutations, and the genes are "
     "ordered by their observed value.",
     "; the thin bar through the filled circle is the 95% range over 20 permutations, "
     "at most 0.076 wide and mostly hidden by the symbol, and the genes are ordered by "
     "their observed value.",
     "Fig 5a：连线与 range 分开说，range 给出实际宽度"),
]

SI_PARAS = [
    (36, 0, "Figure S3. The same decay holds on all seven sections.",
     "Figure S3. The same curves on an absolute scale."),
    (38, 1,
     " α(d) normalised to its value at 2 µm, on the same axes and with the same grey "
     "guides as Fig. 1b, for all seven sections. Exact integer re-aggregation to the same "
     "eight scales was applied to each, so the curves differ only in tissue and in depth, "
     "and the four added sections were fixed before their exponents were computed. Every "
     "curve stays above the counting limit at every scale.",
     " α(d) in absolute units for all seven sections, on the same log–log axes and the "
     "same eight scales. Exact integer re-aggregation was applied to each, so the curves "
     "differ only in tissue and in depth, and the four added sections were fixed before "
     "their exponents were computed. At 2 µm the seven span 21-fold; at 256 µm they span "
     "3.8-fold, because the sections whose exponent is shallow start low and stay low "
     "while the steep ones fall to meet them. The level at 2 µm depends on composition "
     "and depth and is not comparable between sections, which is why Fig. 1b normalises; "
     "the fall is."),
]

FIG_DIR = PKG / "04_图件"
FIG_S3 = FIG_DIR / "FigS3.png"


def body_words(paras: list[str]) -> int:
    return len(" ".join(paras[15:65]).split())


def replace_once(text: str, old: str, new: str, where: str) -> str:
    n = text.count(old)
    if n != 1:
        raise SystemExit(f"{where}: 片段出现 {n} 次（要求 1 次）\n  {old[:80]!r}")
    return text.replace(old, new)


def patch_manuscript(dry: bool) -> dict:
    d = docx.Document(str(MAN))
    paras = [p.text for p in d.paragraphs]
    before_body = body_words(paras)
    before_abs = len(paras[13].split())

    for idx, old, new, why in BODY:
        p = d.paragraphs[idx]
        if len(p.runs) != 1:
            raise SystemExit(f"段 {idx} 有 {len(p.runs)} 个 run，不能整段替换")
        p.runs[0].text = replace_once(p.runs[0].text, old, new, f"正文段 {idx}")

    for idx, ridx, old, new, why in CAPTIONS:
        r = d.paragraphs[idx].runs[ridx]
        r.text = replace_once(r.text, old, new, f"图注段 {idx} run{ridx}")

    paras = [p.text for p in d.paragraphs]
    after_body = body_words(paras)
    after_abs = len(paras[13].split())

    if after_body > 3500:
        raise SystemExit(f"正文超限: {after_body} > 3500")
    if after_abs > 150:
        raise SystemExit(f"摘要超限: {after_abs} > 150")
    if len(paras) != 165:
        raise SystemExit(f"段落数变了: {len(paras)}")

    if not dry:
        d.save(str(MAN))
        pdf = MAN.with_suffix(".pdf")
        if not pdf.exists():
            print(f"  [warn] 还没有导出 PDF: {pdf}")

    return dict(body_before=before_body, body_after=after_body,
                abstract_before=before_abs, abstract_after=after_abs,
                paragraphs=len(paras))


def patch_supplementary(dry: bool) -> dict:
    d = docx.Document(str(SI))
    for idx, ridx, old, new in SI_PARAS:
        r = d.paragraphs[idx].runs[ridx]
        r.text = replace_once(r.text, old, new, f"SI 段 {idx} run{ridx}")

    # 内嵌的 Figure S3 换成绝对值版
    old = Path(r"C://Users//唐朝华//Desktop//23空间分辨率"
               r"\A题_空间分辨率_投稿包_v5.5_20261116\04_图件\FigS3.png")
    rec = replace_embedded(d, old, FIG_S3, width_emu=5903999, dry=dry)
    if not dry:
        d.save(str(SI))
    return rec


def replace_embedded(doc, old_png: Path, new_png: Path, width_emu: int,
                     dry: bool) -> dict:
    """把 docx 里那张与 old_png 逐字节相同的图换成 new_png，并同步 extent。

    按字节认图，不按顺序认——顺序一变就换错图，而字节比对是可验证的。
    """
    from PIL import Image
    import hashlib

    old_sha = hashlib.sha256(old_png.read_bytes()).hexdigest()
    blob = new_png.read_bytes()
    new_sha = hashlib.sha256(blob).hexdigest()
    target = None
    for part in doc.part.package.iter_parts():
        try:
            b = part.blob
        except Exception:
            continue
        if part.content_type.startswith("image/") and \
                hashlib.sha256(b).hexdigest() == old_sha:
            target = part
            break
    if target is None:
        raise SystemExit(f"docx 里找不到 {old_png.name} 的旧字节（sha {old_sha[:12]}）")
    shape = None
    for sh in doc.inline_shapes:
        rid = sh._inline.graphic.graphicData.pic.blipFill.blip.embed
        if doc.part.related_parts[rid] is target:
            shape = sh
            break
    if shape is None:
        raise SystemExit("找不到引用该图件的 inline shape")

    w, h = Image.open(new_png).size
    new_w, new_h = Emu(width_emu), Emu(int(round(width_emu * h / w)))
    if not dry:
        target._blob = blob
        shape.width = new_w
        shape.height = new_h
    return dict(png=[w, h], width_mm=round(new_w / 36000, 1),
                height_mm=round(new_h / 36000, 1),
                new_sha256=new_sha[:16], old_sha256=old_sha[:16],
                part=str(target.partname))


def patch_figures_in_manuscript(dry: bool) -> list[dict]:
    """主稿正文里的 Fig1 / Fig5 也各有内嵌副本，必须一起换。

    这两张图在 04_图件 里是独立文件，但 Word 文档里存的是内嵌副本；
    只换 04_图件 不换内嵌副本，PDF 里还是旧图（v5.5.1 第一轮就踩了）。
    """
    legacy = Path(r"C://Users//唐朝华//Desktop//23空间分辨率"
                  r"\A题_空间分辨率_投稿包_v5.5_20261116\04_图件")
    d = docx.Document(str(MAN))
    out = []
    for name in ("Fig1", "Fig5"):
        rec = replace_embedded(d, legacy / f"{name}.png", FIG_DIR / f"{name}.png",
                               width_emu=5943600, dry=dry)
        rec["figure"] = name
        out.append(rec)
    if not dry:
        d.save(str(MAN))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    dry = a.dry_run
    OUT.parent.mkdir(parents=True, exist_ok=True)

    print(f"package {PKG.name}  dry_run={dry}")
    m = patch_manuscript(dry)
    print(f"  body {m['body_before']} -> {m['body_after']} / 3500"
          f"    abstract {m['abstract_before']} -> {m['abstract_after']} / 150")
    figs = patch_figures_in_manuscript(dry)
    for f in figs:
        print(f"  manuscript {f['figure']}: {f['png']} px -> {f['width_mm']} x "
              f"{f['height_mm']} mm  (was {f['old_sha256']})")
    s = patch_supplementary(dry)
    print(f"  SI image -> {s['png']} px, {s['width_mm']} x {s['height_mm']} mm, "
          f"sha {s['new_sha256']}")

    rec = dict(change="v5.5 -> v5.5.1: M1 five holes",
               manuscript=m, supplementary=s, manuscript_figures=figs,
               body_edits=[dict(para=i, why=w) for i, _, _, w in BODY],
               caption_edits=[dict(para=i, run=r, why=w) for i, r, _, _, w in CAPTIONS],
               numbers_source="02_OUT/m1_v551_numbers.json (125_m1_numbers.py)",
               figures=dict(Fig1="123_fig1_v551.py", Fig5="124_fig5_v551.py",
                            FigS3="126_fig_s3_v551.py"))
    if not dry:
        (OUT).write_text(json.dumps(rec, ensure_ascii=False, indent=1), encoding="utf-8")
        print(f"  written {OUT}")
    else:
        print("  dry run: nothing written")
    return 0


if __name__ == "__main__":
    sys.exit(main())

# -*- coding: utf-8 -*-
"""110_patch_v56.py — v5.5.1 → v5.6：把 M2 的单分子检查写进 SI Note 1。

改动范围（最小版，不动 Results / Discussion / 摘要 / 图）：
  MAN 第 59 段：limitation 归因从"单一原因"放宽为"部分原因"，Table 引用加 S9
  MAN 第 64 段：点尺度测量"已经做了、并没有解决"，指向 Note 1
  SI  第 20 段：'a hard limit' → 'one limit, and not the only one'
  SI  第 21 段：加一句"限制不随网格消失"
  SI  新增一节 + Table S9：单分子坐标上的第三平台检查

原则：v5.5.1 目录一字不动，只写 v5.6。
用法：
  python 110_patch_v56.py --dry-run
  python 110_patch_v56.py
"""
from __future__ import annotations

import argparse
import copy
import json
import sys
from pathlib import Path

import docx
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt

OLD = Path(r"C:\Users\唐朝华\Desktop\23空间分辨率\A题_空间分辨率_投稿包_v5.5.1_20261117")
PKG = Path(r"C:\Users\唐朝华\Desktop\23空间分辨率\A题_空间分辨率_投稿包_v5.6_20261117")
BUILD = Path(r"C:\Users\唐朝华\Desktop\scaling-report\reproduce\07_SUBMISSION\_build_v56")
OUT = BUILD / "audit" / "patch_v5_6.json"

MAN = PKG / "02_主稿" / "01_Manuscript.docx"
SI = PKG / "03_补充材料" / "02_Supplementary_Information.docx"
PRED = Path(r"E:\xen_p2\out\xen_predict_closing10.json")

# ---------------------------------------------------------------- 主稿
BODY = [
    (59,
     "fail for that reason, and are reported in full (Table S5, S6; Fig. S1)",
     "fail, partly for that reason, and are reported in full (Table S5, S6, S9; Fig. S1)",
     "归因从单一原因放宽为部分原因；Table S9 = 单分子重跑"),
    (64,
     "rather than counting them in a 2 µm square; imaging-based methods already "
     "reach that scale [4–6].",
     "rather than counting them; imaging-based methods reach that scale [4–6], "
     "and one such section analysed here reproduces the failure (Note 1).",
     "点尺度测量已经做过，没有解决；删 'in a 2 µm square' 与 'already' 腾词"),
]

# ---------------------------------------------------------------- SI 既有段落
SI_PARAS = [
    (20,
     "but a hard limit of this data set. Measuring",
     "but one limit of this data set, and not the only one (below). Measuring",
     "'硬限制'→'其中一个限制'"),
    (21,
     "cannot be verified. What can be compared",
     "cannot be verified, and the single-molecule check below shows that the limit "
     "does not disappear when the grid does. What can be compared",
     "限制不随网格消失"),
]

# ---------------------------------------------------------------- SI 新增一节
# (文本, bold, 字号pt, 说明)
NEW = [
    ("A third platform: single-molecule coordinates", True, 10.5, "小标题"),
    ("The cross-chemistry comparison in the main text uses Stereo-seq and Visium HD, "
     "both of which count transcripts inside a grid square. A platform that localises "
     "individual transcripts allows the same two quantities to be measured without a "
     "grid at all, and therefore asks whether the split reported there belongs to "
     "counting chemistries or to the binned statistic itself.", False, 10, "引言"),
    ("We used Xenium 2.0 human colon cancer (section P2), 66,248,117 transcripts in a "
     "422-gene panel, positions reported at 0.5 µm. After filtering to qv ≥ 20 and to "
     "transcripts inside the tissue mask, 50,594,830 remain. Exponents were computed "
     "exactly as for Visium HD, by re-aggregating the same transcripts into 2, 4, 8, "
     "16, 32, 64, 128 and 256 µm blocks; the pooled exponent is the µ²-weighted sum "
     "used throughout.", False, 10, "数据与口径"),
    ("The pooled exponent does not transfer. On the 408 genes shared between the panel "
     "and the colon section, the exponent over 2–64 µm is 0.42 on Visium HD and 0.10 "
     "on Xenium, each weighted by its own genes — a gap of 0.32, or 0.34 with both "
     "weighted by the Visium HD genes — "
     "against 0.20 between Stereo-seq and Visium HD on matched "
     "whole-transcriptome data. The gap is not a detection-efficiency artefact: "
     "thinning the Visium HD counts gene by gene to the Xenium per-gene totals moves "
     "the exponent by 0.001. It is partly a panel effect — restricting the same "
     "counts to those 408 genes moves the exponent by 0.089 — which is why the "
     "comparison is made on the shared set.", False, 10, "pooled 不转移"),
    ("The per-gene ordering does transfer. Over the 403 genes that pass the fit on "
     "both platforms, Spearman ρ = 0.871 (p < 10⁻¹²⁰), against 0.769 between the two "
     "counting chemistries, and the median per-gene exponent is 0.54 on Visium HD "
     "and 0.41 on Xenium. The two platforms therefore disagree on the pooled number "
     "and agree on the ordering — the same split the main text reports.", False, 10,
     "per-gene 转移"),
    ("The sub-8 µm range is not usable. There α(d) is not monotone: it peaks at 4 µm "
     "(α = 3.11) and turns negative at 0.5 µm (α = −0.45), the signature of the "
     "localisation kernel smoothing the point density rather than of tissue "
     "covariance. The local slope is +0.56 over 1–8 µm and −0.19 over 8–64 µm, so "
     "the two ranges cannot be fitted as one power law.", False, 10, "8 µm 以下不可用"),
    ("For the same reason the integral prediction cannot be tested there either. "
     "Integrating the measured g(r) over the disc of radius d and comparing with the "
     "observed α(d) of the pooled transcript field, normalised at 2 µm, gives a "
     "predicted slope of −0.343 against an observed −0.024 over 1–256 µm; the "
     "largest deviation is 0.83 log10 units at 256 µm, a factor of 6.8 (Table S9). "
     "Removing the 8 µm grid therefore does not rescue the integral. The residual "
     "failure is consistent with the stationarity assumption that the displacement "
     "average requires: the intensity field of this tumour carries gradients over "
     "more than 100 µm, which a locally normalised g(r) averages away while the "
     "block variance retains them.", False, 10, "积分预测仍失败"),
    ("We do not compare exponents across these two platforms on all genes, because "
     "the panel does not contain IGKC, IGHG1 or JCHAIN and because two-point "
     "functions have already been compared across these platforms (Ilin et al., "
     "bioRxiv 2026, doi:10.64898/2026.10.09.757915); what is asked here is what "
     "their integral means for the binned variance, which is a different quantity "
     "and the one every study reports.", False, 10, "划界"),
    ("Table S9. The integral prediction repeated on single-molecule coordinates. "
     "Xenium 2.0 human colon cancer (P2), 422-gene panel, qv ≥ 20, 50,594,830 "
     "transcripts inside the tissue mask. Column 2: α predicted from the integral of "
     "g(r) over the disc of radius d on a 0.5 µm raster. Column 3: α observed after "
     "exact re-aggregation of the same transcripts to d. Column 4: log10 of the "
     "ratio after normalising both to 2 µm, so zero at that scale by construction. "
     "Slopes over 1–256 µm: predicted −0.343, observed −0.024.", False, 9, "表注"),
]


def body_words(paras: list[str]) -> int:
    return len(" ".join(paras[15:65]).split())


def replace_once(text: str, old: str, new: str, where: str) -> str:
    n = text.count(old)
    if n != 1:
        raise SystemExit(f"{where}: 片段出现 {n} 次（要求 1 次）\n  {old[:90]!r}")
    return text.replace(old, new)


def make_run(template_p, text: str, bold: bool, pt: float):
    """复制模板段落的 pPr，造一个新段落。"""
    p = copy.deepcopy(template_p._p)
    for child in list(p):
        if child.tag != qn('w:pPr'):
            p.remove(child)
    r = OxmlElement('w:r')
    rPr = OxmlElement('w:rPr')
    if bold:
        b = OxmlElement('w:b')
        rPr.append(b)
    sz = OxmlElement('w:sz')
    sz.set(qn('w:val'), str(int(round(pt * 2))))
    rPr.append(sz)
    rf = OxmlElement('w:rFonts')
    rf.set(qn('w:ascii'), 'Times New Roman')
    rf.set(qn('w:hAnsi'), 'Times New Roman')
    rPr.append(rf)
    r.append(rPr)
    t = OxmlElement('w:t')
    t.set(qn('xml:space'), 'preserve')
    t.text = text
    r.append(t)
    p.append(r)
    return p


def build_table9(doc) -> object:
    """按 Table S5（tables[4]）的版式造一张 4 列 × 10 行的表。"""
    pred = json.loads(PRED.read_text(encoding="utf-8"))
    rows = pred["rows"]
    # 2 µm 归一
    ref = [r for r in rows if abs(r["d"] - 2.0) < 1e-9][0]
    scale = ref["alpha_true"] / ref["alpha_pred"]

    src = doc.tables[4]._tbl
    new_tbl = copy.deepcopy(src)
    # 删第 5 列
    for tr in new_tbl.findall(qn('w:tr')):
        tcs = tr.findall(qn('w:tc'))
        if len(tcs) > 4:
            tr.remove(tcs[-1])
    grid = new_tbl.find(qn('w:tblGrid'))
    gcs = grid.findall(qn('w:gridCol'))
    while len(gcs) > 4:
        grid.remove(gcs[-1])
        gcs = grid.findall(qn('w:gridCol'))
    # 行数 → 10
    trs = new_tbl.findall(qn('w:tr'))
    while len(trs) < 10:
        new_tbl.append(copy.deepcopy(trs[1]))
        trs = new_tbl.findall(qn('w:tr'))
    while len(trs) > 10:
        new_tbl.remove(trs[-1])
        trs = new_tbl.findall(qn('w:tr'))

    from docx.table import Table
    tbl = Table(new_tbl, doc.part)
    head = ["d (µm)", "α predicted", "α observed", "log10 deviation"]
    for c, v in enumerate(head):
        tbl.cell(0, c).text = v
    for i, r in enumerate(rows, start=1):
        pn = r["alpha_pred"] * scale
        ov = r["alpha_true"]
        dev = 0.0 if abs(r["d"] - 2.0) < 1e-9 else \
            float(__import__("math").log10(abs(pn / ov)))
        ds = "0.000" if abs(dev) < 5e-4 else f"{dev:+.3f}".replace("-", "\u2212")
        vals = [f"{r['d']:g}", f"{r['alpha_pred']:.3f}", f"{ov:.3f}", ds]
        for c, v in enumerate(vals):
            tbl.cell(i, c).text = v
    for row in tbl.rows:
        for cell in row.cells:
            for p in cell.paragraphs:
                for run in p.runs:
                    run.font.size = Pt(9)
    return new_tbl


def patch_manuscript(dry: bool) -> dict:
    d = docx.Document(str(MAN))
    paras = [p.text for p in d.paragraphs]
    before_body, before_abs = body_words(paras), len(paras[13].split())

    for idx, old, new, why in BODY:
        p = d.paragraphs[idx]
        if len(p.runs) != 1:
            raise SystemExit(f"段 {idx} 有 {len(p.runs)} 个 run")
        p.runs[0].text = replace_once(p.runs[0].text, old, new, f"正文段 {idx}")

    paras = [p.text for p in d.paragraphs]
    after_body, after_abs = body_words(paras), len(paras[13].split())

    if after_body > 3500:
        raise SystemExit(f"正文超限: {after_body} > 3500")
    if after_abs > 150:
        raise SystemExit(f"摘要超限: {after_abs} > 150")
    if len(paras) != 165:
        raise SystemExit(f"主稿段落数变了: {len(paras)}")
    if not dry:
        d.save(str(MAN))
    return dict(body_before=before_body, body_after=after_body,
                abstract_before=before_abs, abstract_after=after_abs,
                paragraphs=len(paras))


def patch_supplementary(dry: bool) -> dict:
    d = docx.Document(str(SI))
    n_par_before = len(d.paragraphs)

    for idx, old, new, why in SI_PARAS:
        p = d.paragraphs[idx]
        done = False
        for r in p.runs:
            if old in r.text:
                r.text = replace_once(r.text, old, new, f"SI 段 {idx}")
                done = True
                break
        if not done:
            raise SystemExit(f"SI 段 {idx} 找不到片段: {old[:70]!r}")

    anchor = d.paragraphs[21]._p
    for text, bold, pt, why in NEW:
        el = make_run(d.paragraphs[21], text, bold, pt)
        anchor.addnext(el)
        anchor = el
    tbl = build_table9(d)
    # 表放在表注之后（表注是最后一节，anchor 指向它）
    anchor.addnext(tbl)

    n_par_after = len(d.paragraphs)
    if not dry:
        d.save(str(SI))
    return dict(paras_before=n_par_before, paras_after=n_par_after,
                tables_before=8, tables_after=len(d.tables),
                new_section_words=len(" ".join(t for t, b, p, w in NEW).split()))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    dry = a.dry_run
    OUT.parent.mkdir(parents=True, exist_ok=True)

    if not PKG.exists():
        raise SystemExit(f"v5.6 目录不存在: {PKG}（先 cp -r v5.5.1 v5.6）")

    print(f"old {OLD.name}  (untouched)")
    print(f"new {PKG.name}  dry_run={dry}")
    m = patch_manuscript(dry)
    print(f"  MAN body {m['body_before']} -> {m['body_after']} / 3500"
          f"    abstract {m['abstract_before']} -> {m['abstract_after']} / 150"
          f"    paras {m['paragraphs']}")
    s = patch_supplementary(dry)
    print(f"  SI paras {s['paras_before']} -> {s['paras_after']}"
          f"    tables 8 -> {s['tables_after']}"
          f"    new section {s['new_section_words']} words")

    rec = dict(change="v5.5.1 -> v5.6: single-molecule check enters SI Note 1",
               manuscript=m, supplementary=s,
               body_edits=[dict(para=i, why=w) for i, _, _, w in BODY],
               si_edits=[dict(para=i, why=w) for i, _, _, w in SI_PARAS],
               new_section="A third platform: single-molecule coordinates + Table S9",
               numbers_source="E:/xen_p2/out/xen_predict_closing10.json, "
                              "xen_alpha_closing10.json, gamma_compare.json, "
                              "xen_thin06b.json, M2_SUMMARY.json")
    if not dry:
        (OUT).write_text(json.dumps(rec, ensure_ascii=False, indent=1),
                         encoding="utf-8")
        print(f"  written {OUT}")
    else:
        print("  dry run: nothing written")
    return 0


if __name__ == "__main__":
    sys.exit(main())

# -*- coding: utf-8 -*-
"""110_patch_v57.py — v5.6 → v5.7：第三片结肠癌（P5）入稿 + P1 身份更正。

三件事，逐条可核：
  1. 数据集身份：v5.6 把 P1 写成 "a second section of the same tumour block"，是错的。
     P1 是另一个患者（Y 染色体基因 per-bin 计数全为 0 → 女性；P2 是男性）。
     改为 "two further colon cancers from different patients"。
  2. 第三片入稿：P5 的 pooled fold = 19.86（八片最大），标题区间 2–15 → 2–20。
     正文/摘要/图注/表格/图件同步，一处都不能漏。
  3. 逐基因与池化的**双向解耦**：P1 对 P2 是 pooled 差 0.006、per-gene 中位差 0.19；
     P5 对 P2 是 pooled 差 0.075、per-gene 中位差 0.012。两个方向都解耦，
     这是"pooled 描述组成不描述组织"最直接的一对证据，写进正文 [46] 与 SI Note 2。

原则：v5.6 目录一字不动，只写 v5.7。字数上限 3500 / 150，超了直接退出。
用法：
  python 110_patch_v57.py --dry-run
  python 110_patch_v57.py
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import sys
from pathlib import Path

import docx
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt

OLD = Path(r"C:\Users\唐朝华\Desktop\23空间分辨率\A题_空间分辨率_投稿包_v5.6_20261117")
PKG = Path(r"C:\Users\唐朝华\Desktop\23空间分辨率\A题_空间分辨率_投稿包_v5.7_20261118")
BUILD = Path(r"C:\Users\唐朝华\Desktop\scaling-report\reproduce\07_SUBMISSION\_build_v57")
OUT = BUILD / "audit" / "patch_v5_7.json"

MAN = PKG / "02_主稿" / "01_Manuscript.docx"
SI = PKG / "03_补充材料" / "02_Supplementary_Information.docx"
FIGDIR = BUILD / "04_图件"

P5A = Path(r"E:\A_spatial_20261116\OUT\human_colon_cancer_p5\pooled_alpha_noMT.json")
P5B = Path(r"E:\A_spatial_20261116\OUT\human_colon_cancer_p5\boot_tail.json")
P5C = Path(r"E:\A_spatial_20261116\OUT\human_colon_cancer_p5\pergene_three.json")
P5D = Path(r"E:\A_spatial_20261116\OUT\human_colon_cancer_p5\block_boot\boot_colonp5.json")


# ---------------------------------------------------------------- 主稿文本
def man_body(p5):
    """(段号, 旧, 新, 说明)。所有数字由 p5 字典注入。"""
    fold = p5["fold"]
    lo, hi = p5["fold_ci"]
    g, pg = p5["gamma"], p5["pergene"]
    d_pool_12 = abs(g["P1"] - g["P2"])
    d_pg_12 = abs(pg["P1"] - pg["P2"])
    d_pool_52 = abs(g["P5"] - g["P2"])
    d_pg_52 = abs(pg["P5"] - pg["P2"])
    return [
        (0,
         "falls 2- to 15-fold, not 16,384-fold",
         "falls 2- to 20-fold, not 16,384-fold",
         "标题区间上界 15.6 → 19.9"),
        (13,
         "the native 2 \u00b5m bins of seven Visium HD sections",
         "the native 2 \u00b5m bins of eight Visium HD sections",
         "摘要：七片 → 八片"),
        (13,
         "falls of 2.1- to 15.6-fold",
         "falls of 2.1- to 19.9-fold",
         "摘要：区间同步"),
        (20,
         "\u22120.554",
         f"\u2212{g['P2']:.3f}",
         "四舍五入：colon 0.554 → 0.555（上标 run 单独一格）"),
        (20,
         "plus four further sections (normal pancreas, small intestine, kidney, embryo)",
         "plus five further sections (normal pancreas, small intestine, kidney, "
         "embryo, and a third colon tumour)",
         "引言：4 → 5 片"),
        (25,
         "exponent \u03b3 = 0.554 (Fig. 1b)",
         f"exponent \u03b3 = {g['P2']:.3f} (Fig. 1b)",
         "四舍五入同 [20]"),
        (25,
         "Four further sections fall 15.6-fold (small intestine), 3.5-fold (kidney), "
         "2.1-fold (pancreas) and 2.3-fold (embryo). The largest has a bootstrap 95% "
         "confidence interval of 14.3 to 18.1, which overlaps the colon value.",
         f"Five further sections fall 15.6- (small intestine), 3.5- (kidney), "
         f"2.1- (pancreas), 2.3- (embryo) and {fold:.1f}-fold (the third colon "
         f"tumour). Its bootstrap interval, {lo:.1f} to {hi:.1f}, does not overlap "
         f"the 14.3 above.",
         "Results：加 P5；最大值的 CI 换人"),
        (46,
         "Reproducibility follows the same split. A second section of the same block, "
         "captured separately, falls 13.1-fold against 14.3-fold, with \u03b3 = 0.549 "
         "\u00b1 0.026 against 0.554, a difference of 0.006 (block bootstrap z = 0.20); "
         "its per-gene median is 0.618 against 0.428, a gap of 0.19 that no abundance "
         "correction removes.",
         f"Reproducibility splits in both directions. Two further colon tumours from "
         f"different patients fall 13.1- and {fold:.1f}-fold against 14.3-fold. "
         f"Against the 14.3-fold section, the 13.1-fold one moves the pooled exponent "
         f"by {d_pool_12:.3f} and the per-gene median by {d_pg_12:.2f}, while the "
         f"{fold:.1f}-fold one moves the pooled exponent by {d_pool_52:.3f} and the "
         f"per-gene median by {d_pg_52:.3f}; neither summary tracks the other "
         f"(Note 2).",
         "重复性：一片动 pooled 不动 per-gene，另一片反过来"),
        (46,
         "Within a single section, cutting the colon into four column bands moves the "
         "pooled exponent over 0.131 and the per-gene median over 0.098, with "
         "between-band standard deviations of 0.056 and 0.042.",
         "Within one section, four column bands move the pooled exponent over 0.131 "
         "and the per-gene median over 0.098 (band SDs 0.056, 0.042).",
         "腾词：四带那句压缩 9 词"),
        (60,
         "and there is one section per tissue",
         "and only colon is replicated",
         "limitation：不再是每组织一片"),
        (67,
         "We analysed eight Visium HD sections and one Stereo-seq section.",
         "We analysed nine Visium HD sections and one Stereo-seq section.",
         "Methods：8 → 9 片"),
        (67,
         "a second section of the same tumour block (Visium_HD_Human_Colon_Cancer_P1)",
         "two further colon cancers from different patients "
         "(Visium_HD_Human_Colon_Cancer_P1 and _P5)",
         "P1 身份更正：不是同一 block 的第二片"),
        (67,
         "The second colon section is a replicate; the other seven enter the scaling "
         "analysis, one per tissue.",
         "The two further colon sections are replicates; the other seven enter the "
         "scaling analysis, one per tissue.",
         "replicate 单数 → 复数"),
        (67,
         "The four sections added here were fixed before their exponents were computed",
         "The four further tissues were fixed before their exponents were computed",
         "预注册的是四个组织，P1/P5 是后来的重复片"),
        (88,
         "The two colon sections were processed by the identical pipeline.",
         "The three colon sections were processed by the identical pipeline.",
         "Methods：两片 → 三片"),
        (102,
         "The eight Visium HD sections and the Stereo-seq section analysed here",
         "The nine Visium HD sections and the Stereo-seq section analysed here",
         "Data availability：8 → 9 片（P1 也算 analysed）"),
        (154,
         "against bin size, for all seven sections",
         "against bin size, for all eight sections",
         "Fig 1 注：七条 → 八条"),
        (154,
         "from 15.6-fold (small intestine) to 2.1-fold (pancreas)",
         f"from {fold:.1f}-fold (colon P5) to 2.1-fold (pancreas)",
         "Fig 1 注：区间两端换人"),
        (154,
         "If the working assumption of the field held, the seven curves would lie",
         "If the working assumption of the field held, the eight curves would lie",
         "Fig 1 注：curves 数"),
    ]


# ---------------------------------------------------------------- SI 既有段落
def si_paras(p5):
    return [
        (1,
         "falls 2- to 15-fold, not 16,384-fold",
         "falls 2- to 20-fold, not 16,384-fold",
         "SI 标题同步"),
        (4,
         "The colon P1 section is a second section of the same tumour block and is "
         "used only as a replicate for the pooled exponent; it is not part of the "
         "eight-scale analysis.",
         "The two further colon sections, P1 and P5, come from different patients and "
         "are used as replicates: P1 for the pooled and per-gene comparison, P5 also "
         "as the largest fall in Fig. 1b.",
         "Table S1 注：P1 身份 + P5 用途"),
        (4,
         "All four Visium HD sections were obtained as Space Ranger 2 \u00b5m output",
         "All Visium HD sections were obtained as Space Ranger 2 \u00b5m output",
         "Table S1 注：删过期的 'four'"),
        (47,
         "for all seven sections, on the same log\u2013log axes",
         "for all eight sections, on the same log\u2013log axes",
         "Fig S3 注：七条 → 八条"),
        (47,
         "At 2 \u00b5m the seven span 21-fold; at 256 \u00b5m they span 3.8-fold",
         "At 2 \u00b5m the eight span 27-fold; at 256 \u00b5m they span 3.8-fold",
         "Fig S3 注：跨度的数现场算过"),
        (47,
         "and the four added sections were fixed before their exponents were computed",
         "and the four sections of the first extension were fixed before their "
         "exponents were computed, the two colon replicates being added afterwards "
         "by the same pipeline",
         "Fig S3 注：预注册的是四片，重复片是后来的"),
    ]


# ---------------------------------------------------------------- SI 新增 Note 2
def _minus(x):
    """带符号的三位小数，负号用 U+2212（Word 排版用减号，不是连字符）。"""
    return f"\u2212{abs(x):.3f}" if x < 0 else f"+{x:.3f}"


def _plus(x):
    return _minus(x)


def note2(p5):
    g = p5["gamma"]
    pg = p5["pergene"]
    rho = p5["rho"]
    return [
        ("Supplementary Note 2. Replication within one tumour type.", True, 10.5, "标题"),
        ("The three colon sections are the only place in this data set where the same "
         "tumour type was captured more than once on the same platform, by the same "
         "pipeline, from three different patients. They therefore separate two things "
         "that a single section confounds: how much the pooled exponent moves between "
         "patients, and how much a typical gene's exponent moves over the same "
         "comparison.", False, 10, "为什么做"),
        (f"Pooled, the three fall 13.1-, 14.3- and {p5['fold']:.1f}-fold between 2 and "
         f"256 \u00b5m, with exponents {g['P1']:.3f} \u00b1 {p5['se']['P1']:.3f}, "
         f"{g['P2']:.3f} and {g['P5']:.3f} \u00b1 {p5['se']['P5']:.3f} "
         f"(block bootstrap over 512 \u00b5m superblocks, 200 resamples). Per gene, the "
         f"medians over the genes that pass the fit in each section are "
         f"{pg['P1']:.3f}, {pg['P2']:.3f} and {pg['P5']:.3f} "
         f"({p5['nfit']['P1']:,}, {p5['nfit']['P2']:,} and {p5['nfit']['P5']:,} genes). "
         f"The two summaries disagree in both directions. P1 against P2 moves the "
         f"pooled exponent by {abs(g['P1'] - g['P2']):.3f} but the per-gene median by "
         f"{abs(pg['P1'] - pg['P2']):.3f}; P5 against P2 moves the pooled exponent by "
         f"{abs(g['P5'] - g['P2']):.3f} but the per-gene median by "
         f"{abs(pg['P5'] - pg['P2']):.3f}.", False, 10, "全部数字"),
        (f"Of the {p5['rho_n_all']:,} genes that pass the fit in all three, the "
         f"pairwise correlations use the {rho['n']:,} with a finite value in every "
         f"pairing; over those, the per-gene exponents correlate at "
         f"\u03c1 = {rho['P1_P2']['spearman']:.3f} (P1\u2013P2), "
         f"{rho['P1_P5']['spearman']:.3f} (P1\u2013P5) and "
         f"{rho['P2_P5']['spearman']:.3f} (P2\u2013P5), and the median differences are "
         f"{_minus(rho['P1_P2']['gap'])}, {_minus(rho['P1_P5']['gap'])} and "
         f"{_plus(rho['P2_P5']['gap'])}. The ordering therefore replicates in all "
         f"three pairings, most tightly between the two sections whose pooled "
         f"exponents differ most. That is the direction the main text reports, now "
         f"with two independent replicates rather than one.", False, 10, "排序一致性"),
        ("One caveat bounds the P1 comparison. Its sequencing depth is 13.0 transcripts "
         "per 2 \u00b5m bin against 30.5 in P2 and 27.6 in P5, so the 0.19 gap in the "
         "per-gene median could in principle be a depth effect. Binomial thinning of "
         "the P2 field to 0.1 transcripts per bin moves the pooled exponent by 0.001 "
         "(Fig. 3d), and the per-gene median gap survives matching on mean abundance, "
         "but a depth-matched re-analysis of a deeper P1 library was not available and "
         "the gap is reported as measured, not corrected.", False, 10, "深度混杂"),
    ]


def body_words(paras):
    return len(" ".join(paras[15:65]).split())


def replace_once(text, old, new, where):
    n = text.count(old)
    if n != 1:
        raise SystemExit(f"{where}: 片段出现 {n} 次（要求 1 次）\n  {old[:100]!r}")
    return text.replace(old, new)


def make_run(template_p, text, bold, pt):
    p = copy.deepcopy(template_p._p)
    for child in list(p):
        if child.tag != qn('w:pPr'):
            p.remove(child)
    r = OxmlElement('w:r')
    rPr = OxmlElement('w:rPr')
    if bold:
        rPr.append(OxmlElement('w:b'))
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


def set_runs(par, text, where):
    """把一个段落的全部文字塞进第一个 run，其余 run 清空（保住段落级格式）。"""
    if not par.runs:
        raise SystemExit(f"{where}: 段落没有 run")
    par.runs[0].text = text
    for r in par.runs[1:]:
        r.text = ""


def patch_table_s1(doc, p5, dry):
    t = doc.tables[0]
    if len(t.rows) != 10:
        raise SystemExit(f"Table S1 行数变了: {len(t.rows)}")
    # 第 9 行 = P1（0 是表头）
    p1 = t.rows[9]
    vals = ["Visium HD Human Colon Cancer P1", "Visium HD", "Human",
            "Colon cancer, second patient", "6.5 \u00d7 6.5 mm",
            f"{p5['p1_bins']:,}", f"{p5['p1_perbin']:.1f}", "18,074",
            "10x Genomics (public)"]
    for c, v in enumerate(vals):
        if c < len(p1.cells):
            set_runs(p1.cells[c].paragraphs[0], v, f"Table S1 行9列{c}")
    # 复制该行造 P5 行
    new_tr = copy.deepcopy(p1._tr)
    t._tbl.append(new_tr)
    from docx.table import _Row
    p5row = _Row(new_tr, t)
    vals5 = ["Visium HD Human Colon Cancer P5", "Visium HD", "Human",
             "Colon cancer, third patient", "6.5 \u00d7 6.5 mm",
             f"{p5['p5_bins']:,}", f"{p5['p5_perbin']:.1f}", "18,074",
             "10x Genomics (public)"]
    for c, v in enumerate(vals5):
        if c < len(p5row.cells):
            set_runs(p5row.cells[c].paragraphs[0], v, f"Table S1 行10列{c}")
    return dict(rows_before=10, rows_after=len(t.rows),
                p1_row=[c.text for c in t.rows[9].cells],
                p5_row=[c.text for c in t.rows[10].cells])


def replace_embedded(doc, old_png, new_png, width_emu, dry):
    from PIL import Image
    from docx.shared import Emu
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
    keep_w = int(shape.width) if width_emu is None else int(width_emu)
    new_w, new_h = Emu(keep_w), Emu(int(round(keep_w * h / w)))
    old_mm = (round(shape.width / 36000, 1), round(shape.height / 36000, 1))
    if not dry:
        target._blob = blob
        shape.width = new_w
        shape.height = new_h
    return dict(png=[w, h],
                old_mm=old_mm,
                new_mm=(round(new_w / 36000, 1), round(new_h / 36000, 1)),
                new_sha256=new_sha[:16], old_sha256=old_sha[:16],
                part=str(target.partname))


def patch_manuscript(p5, dry):
    d = docx.Document(str(MAN))
    paras = [p.text for p in d.paragraphs]
    before_body, before_abs = body_words(paras), len(paras[13].split())
    if len(paras) != 165:
        raise SystemExit(f"主稿段落数不是 165: {len(paras)}")

    for idx, old, new, why in man_body(p5):
        p = d.paragraphs[idx]
        if len(p.runs) == 1:
            p.runs[0].text = replace_once(p.runs[0].text, old, new, f"正文段 {idx}")
        else:
            done = False
            for r in p.runs:
                if old in r.text:
                    r.text = replace_once(r.text, old, new, f"正文段 {idx}")
                    done = True
                    break
            if not done:
                raise SystemExit(f"正文段 {idx} 有 {len(p.runs)} 个 run 且找不到片段"
                                 f"\n  {old[:90]!r}")

    paras = [p.text for p in d.paragraphs]
    after_body, after_abs = body_words(paras), len(paras[13].split())
    if after_body > 3500:
        raise SystemExit(f"正文超限: {after_body} > 3500（多 {after_body - 3500} 词）")
    if after_abs > 150:
        raise SystemExit(f"摘要超限: {after_abs} > 150")

    figs = []
    for name, w in (("Fig1", 5943600),):
        figs.append(replace_embedded(d, OLD / "04_图件" / f"{name}.png",
                                     FIGDIR / f"{name}.png", None, dry))
    if not dry:
        d.save(str(MAN))
    return dict(body_before=before_body, body_after=after_body,
                abstract_before=before_abs, abstract_after=after_abs,
                paragraphs=len(paras), figures=figs)


def patch_supplementary(p5, dry):
    d = docx.Document(str(SI))
    n0 = len(d.paragraphs)

    for idx, old, new, why in si_paras(p5):
        p = d.paragraphs[idx]
        done = False
        for r in p.runs:
            if old in r.text:
                r.text = replace_once(r.text, old, new, f"SI 段 {idx}")
                done = True
                break
        if not done:
            raise SystemExit(f"SI 段 {idx} 找不到片段: {old[:80]!r}")

    tbl = patch_table_s1(d, p5, dry)

    anchor = d.paragraphs[-1]._p
    for text, bold, pt, why in note2(p5):
        el = make_run(d.paragraphs[-1], text, bold, pt)
        anchor.addnext(el)
        anchor = el

    figs = [replace_embedded(d, OLD / "04_图件" / "FigS3.png",
                             FIGDIR / "FigS3.png", None, dry)]
    if not dry:
        d.save(str(SI))
    return dict(paras_before=n0, paras_after=len(d.paragraphs),
                table=tbl, figures=figs,
                note2_words=len(" ".join(t for t, b, p, w in note2(p5)).split()))


def load_numbers():
    pa = json.loads(P5A.read_text(encoding="utf-8"))
    pb = json.loads(P5B.read_text(encoding="utf-8"))
    pc = json.loads(P5C.read_text(encoding="utf-8"))
    pd = json.loads(P5D.read_text(encoding="utf-8"))
    lock = json.loads((Path(r"C:\Users\唐朝华\Desktop\scaling-report\reproduce"
                            r"\02_OUT\NUMBERS_LOCK_v7.json")).read_text(encoding="utf-8"))
    ps = lock["panel_subset"]
    a2 = float(pa["per_d"]["2"]["alpha_pool"])
    a256 = float(pa["per_d"]["256"]["alpha_pool"])
    fold = a2 / a256
    ci = pb["per_d"]["256"]["ci95"]          # boot_tail 给的是 alpha256 的区间
    return dict(
        fold=fold,
        fold_ci=(a2 / ci[1], a2 / ci[0]),
        alpha2=a2, alpha256=a256,
        gamma={"P1": float(ps["colonP1"]["gamma_full"]),
               "P2": float(ps["colon"]["gamma_full"]),
               "P5": -float(pa["slope_all"])},
        se={"P1": 0.026, "P5": float(pd["results"]["512"]["sd"])},
        pergene={k: float(v["median"]) for k, v in pc["per_section"].items()},
        nfit={k: int(v["n_fit"]) for k, v in pc["per_section"].items()},
        rho=dict(n=int(pc["pairs"]["P1_P2"]["n"]), **pc["pairs"]),
        rho_n_all=int(pc["common_genes"]),
        p5_boot=dict(sd=float(pd["results"]["512"]["sd"]),
                     mean=float(pd["results"]["512"]["mean"]),
                     q025=float(pd["results"]["512"]["q025"]),
                     q975=float(pd["results"]["512"]["q975"]),
                     nblk=int(pd["results"]["512"]["nblk"]), B=int(pd["B"])),
        p1_bins=8115564, p1_perbin=12.971991595408527,
        p5_bins=8665632, p5_perbin=239246304 / 8665632,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    dry = a.dry_run
    OUT.parent.mkdir(parents=True, exist_ok=True)
    if not PKG.exists():
        raise SystemExit(f"v5.7 目录不存在: {PKG}")

    p5 = load_numbers()
    print(f"P5: fold={p5['fold']:.3f}  CI={p5['fold_ci'][0]:.2f}-{p5['fold_ci'][1]:.2f}"
          f"  gamma={p5['gamma']['P5']:.4f}")
    print(f"per-gene median: {p5['pergene']}")

    m = patch_manuscript(p5, dry)
    print(f"  MAN body {m['body_before']} -> {m['body_after']} /3500"
          f"    abstract {m['abstract_before']} -> {m['abstract_after']} /150")
    for f in m["figures"]:
        print(f"    Fig1 {f['old_mm']} -> {f['new_mm']} mm  sha {f['new_sha256']}")

    s = patch_supplementary(p5, dry)
    print(f"  SI paras {s['paras_before']} -> {s['paras_after']}"
          f"    Table S1 rows {s['table']['rows_before']} -> {s['table']['rows_after']}"
          f"    Note 2 {s['note2_words']} words")
    for f in s["figures"]:
        print(f"    FigS3 {f['old_mm']} -> {f['new_mm']} mm  sha {f['new_sha256']}")

    rec = dict(change="v5.6 -> v5.7: third colon section (P5) enters the paper; "
                      "P1 identity corrected from 'same block' to 'second patient'",
               manuscript=m, supplementary=s,
               body_edits=[dict(para=i, why=w) for i, _, _, w in man_body(p5)],
               si_edits=[dict(para=i, why=w) for i, _, _, w in si_paras(p5)],
               numbers=p5,
               numbers_source=[str(P5A), str(P5B), str(P5C), str(P5D)])
    if not dry:
        OUT.write_text(json.dumps(rec, ensure_ascii=False, indent=1), encoding="utf-8")
        print(f"  written {OUT}")
    else:
        print("  dry run: nothing written")
    return 0


if __name__ == "__main__":
    sys.exit(main())

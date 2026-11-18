# -*- coding: utf-8 -*-
"""111_guard_v57.py — v5.7 守卫。只对本次构建出的这几份文件有效。

改一个字就要重跑；换一份文件也要重跑（v5.6 的 59/59 只对 v5.6 那两份文件成立）。

不引用 110 的任何中间值：每个数字都从 E:/A_spatial_20261116/OUT 与
reproduce/02_OUT 的 JSON 现场重算，再与 docx 的 OOXML 文本、PDF 的渲染文本、
Table S1 的单元格、投稿信正文逐处比对。

用法：python 111_guard_v57.py
退出码 0 = 全过；1 = 有失败项。
"""
from __future__ import annotations

import hashlib
import json
import math
import re
import sys
from pathlib import Path

import docx
import pymupdf

OLD = Path(r"C:\Users\唐朝华\Desktop\23空间分辨率\A题_空间分辨率_投稿包_v5.6_20261117")
PKG = Path(r"C:\Users\唐朝华\Desktop\23空间分辨率\A题_空间分辨率_投稿包_v5.7_20261118")
REPRO = Path(r"C:\Users\唐朝华\Desktop\scaling-report\reproduce")
P5OUT = Path(r"E:\A_spatial_20261116\OUT\human_colon_cancer_p5")

MAN = PKG / "02_主稿" / "01_Manuscript.docx"
SI = PKG / "03_补充材料" / "02_Supplementary_Information.docx"
COVER = PKG / "01_投稿信" / "00_Cover_Letter.docx"
MAN_PDF = PKG / "02_主稿" / "01_Manuscript.pdf"
SI_PDF = PKG / "03_补充材料" / "02_Supplementary_Information.pdf"
COVER_PDF = PKG / "01_投稿信" / "00_Cover_Letter.pdf"
FIGD = PKG / "04_图件"

CHECKS: list[tuple[str, bool, str]] = []


def ck(name, ok, detail=""):
    CHECKS.append((name, bool(ok), str(detail)))


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def docx_text(p: Path) -> str:
    d = docx.Document(str(p))
    parts = [x.text for x in d.paragraphs]
    for t in d.tables:
        for r in t.rows:
            parts.extend(c.text for c in r.cells)
    return "\n".join(parts)


def pdf_text(p: Path) -> str:
    with pymupdf.open(str(p)) as d:
        return re.sub(r"\s+", " ", "".join(pg.get_text() for pg in d))


def flat(s: str) -> str:
    return re.sub(r"\s+", " ", s)


def main() -> int:
    # ---------------------------------------------------------- 0 原包未动
    old_man = json.loads((OLD / "MANIFEST.json").read_text(encoding="utf-8"))
    bad = [e["file"] for e in old_man["files"]
           if not (OLD / e["file"]).exists()
           or (OLD / e["file"]).stat().st_size != e["bytes"]
           or sha(OLD / e["file"]) != e["sha256"]]
    ck("v5.6 原包逐字节未动", not bad, f"changed={bad[:4]}")

    # ---------------------------------------------------------- 1 源头数字
    pa = json.loads((P5OUT / "pooled_alpha_noMT.json").read_text("utf-8"))
    pb = json.loads((P5OUT / "boot_tail.json").read_text("utf-8"))
    pc = json.loads((P5OUT / "pergene_three.json").read_text("utf-8"))
    pd = json.loads((P5OUT / "block_boot" / "boot_colonp5.json").read_text("utf-8"))
    lock = json.loads((REPRO / "02_OUT" / "NUMBERS_LOCK_v7.json").read_text("utf-8"))
    pss = lock["panel_subset"]

    a2 = float(pa["per_d"]["2"]["alpha_pool"])
    a256 = float(pa["per_d"]["256"]["alpha_pool"])
    fold5 = a2 / a256
    ci = pb["per_d"]["256"]["ci95"]
    lo, hi = a2 / ci[1], a2 / ci[0]
    ck("P5 α(2)=40.335", abs(a2 - 40.3352) < 5e-3, f"{a2:.4f}")
    ck("P5 α(256)=2.031", abs(a256 - 2.03071) < 5e-4, f"{a256:.5f}")
    ck("P5 fold=19.86 → 文本 19.9", f"{fold5:.1f}" == "19.9", f"{fold5:.4f}")
    ck("P5 fold CI 16.2–23.4", f"{lo:.1f}==16.2" and f"{hi:.1f}==23.4",
       f"{lo:.2f}-{hi:.2f}")

    g1 = float(pss["colonP1"]["gamma_full"])
    g2 = float(pss["colon"]["gamma_full"])
    g5 = -float(pa["slope_all"])
    ck("γ P1/P2/P5 = 0.549/0.555/0.629",
       abs(g1 - 0.5485) < 5e-4 and abs(g2 - 0.5547) < 5e-4 and abs(g5 - 0.6290) < 5e-4,
       f"{g1:.4f}/{g2:.4f}/{g5:.4f}")
    ck("正文 0.554 → 0.555 已统一", f"{g2:.3f}" == "0.555", f"{g2:.4f}")

    pg = {k: float(v["median"]) for k, v in pc["per_section"].items()}
    nf = {k: int(v["n_fit"]) for k, v in pc["per_section"].items()}
    ck("per-gene 中位 0.618/0.428/0.440",
       f"{pg['P1']:.3f}" == "0.618" and f"{pg['P2']:.3f}" == "0.428"
       and f"{pg['P5']:.3f}" == "0.440",
       f"{pg['P1']:.3f}/{pg['P2']:.3f}/{pg['P5']:.3f}")
    ck("n_fit 12,650 / 13,115 / 13,387",
       (nf["P1"], nf["P2"], nf["P5"]) == (12650, 13115, 13387), str(nf))
    pr = pc["pairs"]
    ck("ρ 0.658/0.686/0.760",
       f"{pr['P1_P2']['spearman']:.3f}" == "0.658"
       and f"{pr['P1_P5']['spearman']:.3f}" == "0.686"
       and f"{pr['P2_P5']['spearman']:.3f}" == "0.760",
       f"{pr['P1_P2']['spearman']:.3f}/{pr['P1_P5']['spearman']:.3f}/"
       f"{pr['P2_P5']['spearman']:.3f}")
    ck("共同基因 12,149（配对用 12,148）",
       pc["common_genes"] == 12149 and pr["P1_P2"]["n"] == 12148,
       f"{pc['common_genes']}/{pr['P1_P2']['n']}")
    sd5 = float(pd["results"]["512"]["sd"])
    ck("P5 γ 的 512 µm 块自助 SD=0.031", f"{sd5:.3f}" == "0.031", f"{sd5:.4f}")

    d_pool_12 = abs(g1 - g2)
    d_pg_12 = abs(pg["P1"] - pg["P2"])
    d_pool_52 = abs(g5 - g2)
    d_pg_52 = abs(pg["P5"] - pg["P2"])
    ck("双向解耦：pooled 0.006/0.074，per-gene 0.19/0.012",
       f"{d_pool_12:.3f}" == "0.006" and f"{d_pool_52:.3f}" == "0.074"
       and f"{d_pg_12:.2f}" == "0.19" and f"{d_pg_52:.3f}" == "0.012",
       f"{d_pool_12:.3f}/{d_pool_52:.3f} | {d_pg_12:.3f}/{d_pg_52:.3f}")

    # 三片各自的 fold（从 lock 的 alpha 数组重算，不用 extra_tissues 的 drop_fold）
    p1a = [float(v) for v in pss["colonP1"]["alpha"]]
    p2a = [float(v) for v in pss["colon"]["alpha"]]
    fold1, fold2 = p1a[0] / p1a[-1], p2a[0] / p2a[-1]
    ck("P1 fold 13.1", f"{fold1:.1f}" == "13.1", f"{fold1:.3f}")
    ck("P2 fold 14.3", f"{fold2:.1f}" == "14.3", f"{fold2:.3f}")
    ck("三片 fold 极差 6.77 且比值 1.52×",
       abs((fold5 - fold1) - 6.766) < 5e-3 and abs(fold5 / fold1 - 1.5166) < 5e-3,
       f"{fold5 - fold1:.3f} / {fold5 / fold1:.3f}")

    # ---------------------------------------------------------- 2 主稿结构
    dm = docx.Document(str(MAN))
    paras = [p.text for p in dm.paragraphs]
    body = len(" ".join(paras[15:65]).split())
    abstract = len(paras[13].split())
    ck("主稿段落数 165", len(paras) == 165, len(paras))
    ck("正文词数 ≤3500", body <= 3500, f"{body}")
    ck("摘要词数 ≤150", abstract <= 150, f"{abstract}")

    mt = flat(docx_text(MAN))
    for s, n, why in [
        ("falls 2- to 20-fold, not 16,384-fold", 1, "标题区间"),
        ("eight Visium HD sections from human and mouse", 1, "摘要片数"),
        ("falls of 2.1- to 19.9-fold", 1, "摘要区间"),
        ("plus five further sections", 1, "[20] 引言"),
        ("and a third colon tumour", 1, "[20]"),
        ("Five further sections fall 15.6-", 1, "[25]"),
        ("3.5- (kidney), 2.1- (pancreas), 2.3- (embryo)", 1, "[25] 另三片"),
        ("fall 13.1- and 19.9-fold against 14.3-fold", 1, "[46] 三片 fold"),
        ("the pooled exponent by 0.006 and the per-gene median by 0.19", 1,
         "[46] P1 对 P2"),
        ("the pooled exponent by 0.074 and the per-gene median by 0.012", 1,
         "[46] P5 对 P2"),
        ("19.9-fold (the third colon tumour)", 1, "[25]"),
        ("Its bootstrap interval, 16.2 to 23.4", 1, "[25]"),
        ("Reproducibility splits in both directions", 1, "[46]"),
        ("Two further colon tumours from different patients", 1, "[46]"),
        ("neither summary tracks the other (Note 2)", 1, "[46] 指向 SI Note 2"),
        ("only colon is replicated", 1, "[60]"),
        ("We analysed nine Visium HD sections", 1, "[67]"),
        ("(Visium_HD_Human_Colon_Cancer_P1 and _P5)", 1, "[67]"),
        ("The two further colon sections are replicates", 1, "[67]"),
        ("The four further tissues were fixed", 1, "[67]"),
        ("The three colon sections were processed", 1, "[88]"),
        ("for all eight sections", 1, "Fig 1 注"),
        ("from 19.9-fold (colon P5) to 2.1-fold (pancreas)", 1, "Fig 1 注"),
        ("the eight curves would lie", 1, "Fig 1 注"),
        ("The nine Visium HD sections and the Stereo-seq section", 1,
         "Data availability"),
        # 旧串必须清零
        ("2- to 15-fold", 0, "旧标题"),
        ("2.1- to 15.6-fold", 0, "旧摘要区间"),
        ("seven Visium HD sections", 0, "旧片数"),
        ("a second section of the same tumour block", 0, "旧 P1 身份"),
        ("The second colon section is a replicate", 0, "旧单数"),
        ("all seven sections", 0, "旧图注"),
        ("the seven curves would lie", 0, "旧图注"),
    ]:
        c = mt.count(s)
        ck(f"主稿 '{s[:42]}' ×{n}", c == n, f"{c}  [{why}]")

    # ---------------------------------------------------------- 3 SI 结构
    ds = docx.Document(str(SI))
    sp = [p.text for p in ds.paragraphs]
    st = flat(docx_text(SI))
    ck("SI 段落 53", len(sp) == 53, len(sp))
    ck("SI 表格 9", len(ds.tables) == 9, len(ds.tables))
    for s, n, why in [
        ("falls 2- to 20-fold, not 16,384-fold", 1, "SI 标题"),
        ("The two further colon sections, P1 and P5", 1, "Table S1 注"),
        ("P5 also as the largest fall in Fig. 1b", 1, "Table S1 注"),
        ("All Visium HD sections were obtained", 1, "Table S1 注"),
        ("for all eight sections", 1, "Fig S3 注"),
        ("At 2 µm the eight span 27-fold", 1, "Fig S3 注"),
        ("the two colon replicates being added afterwards", 1, "Fig S3 注"),
        ("Supplementary Note 2. Replication within one tumour type.", 1, "Note 2"),
        ("13.1-, 14.3- and 19.9-fold", 1, "Note 2 pooled"),
        ("0.549 ± 0.026, 0.555 and 0.629 ± 0.031", 1, "Note 2 γ"),
        ("0.618, 0.428 and 0.440", 1, "Note 2 per-gene"),
        ("moves the pooled exponent by 0.006 but the per-gene median by 0.190", 1,
         "Note 2 P1 对 P2"),
        ("moves the pooled exponent by 0.074 but the per-gene median by 0.012", 1,
         "Note 2 P5 对 P2"),
        ("12,650, 13,115 and 13,387 genes", 1, "Note 2 n_fit"),
        ("the 12,148 with a finite value", 1, "Note 2 ρ n"),
        ("ρ = 0.658 (P1–P2), 0.686 (P1–P5) and 0.760 (P2–P5)", 1, "Note 2 ρ"),
        ("−0.194, −0.189 and +0.005", 1, "Note 2 Δ"),
        ("13.0 transcripts per 2 µm bin against 30.5 in P2 and 27.6 in P5", 1,
         "Note 2 深度"),
        # 旧串
        ("a second section of the same tumour block", 0, "旧 P1 身份"),
        ("All four Visium HD sections", 0, "旧 'four'"),
        ("seven span 21-fold", 0, "旧跨度"),
        ("the seven curves", 0, "旧曲线数"),
    ]:
        c = st.count(s)
        ck(f"SI '{s[:42]}' ×{n}", c == n, f"{c}  [{why}]")

    # ---------------------------------------------------------- 4 Table S1
    t1 = ds.tables[0]
    ck("Table S1 11 行 9 列", len(t1.rows) == 11 and len(t1.columns) == 9,
       f"{len(t1.rows)}x{len(t1.columns)}")
    r9 = [c.text.strip() for c in t1.rows[9].cells]
    r10 = [c.text.strip() for c in t1.rows[10].cells]
    ck("Table S1 P1 行（第二患者）",
       r9[0].endswith("P1") and r9[3] == "Colon cancer, second patient"
       and r9[5] == "8,115,564" and r9[6] == "13.0" and r9[7] == "18,074", str(r9))
    ck("Table S1 P5 行（第三患者）",
       r10[0].endswith("P5") and r10[3] == "Colon cancer, third patient"
       and r10[5] == "8,665,632" and r10[6] == "27.6" and r10[7] == "18,074", str(r10))
    ck("Table S1 无 'same block'", all("same block" not in c.text for r in t1.rows
                                       for c in r.cells), "")

    # ---------------------------------------------------------- 5 载体覆盖
    ck("主稿 PDF 页数 28", pymupdf.open(str(MAN_PDF)).page_count == 28,
       pymupdf.open(str(MAN_PDF)).page_count)
    ck("SI PDF 页数 9", pymupdf.open(str(SI_PDF)).page_count == 9,
       pymupdf.open(str(SI_PDF)).page_count)
    ck("投稿信 PDF 页数 1", pymupdf.open(str(COVER_PDF)).page_count == 1,
       pymupdf.open(str(COVER_PDF)).page_count)

    mpt = pdf_text(MAN_PDF)
    for s in ["2- to 20-fold", "2.1- to 19.9-fold", "eight Visium HD sections",
              "Five further sections fall 15.6-", "19.9-fold (the third colon tumour)",
              "16.2 to 23.4", "Reproducibility splits in both directions",
              "only colon is replicated", "nine Visium HD sections",
              "19.9-fold (colon P5)", "eight curves would lie"]:
        ck(f"主稿 PDF 含 '{s[:36]}'", mpt.count(s) >= 1, mpt.count(s))
    for s in ["2- to 15-fold", "2.1- to 15.6-fold", "seven Visium HD sections",
              "same tumour block", "15.6-fold (small intestine) to 2.1-fold"]:
        ck(f"主稿 PDF 无旧串 '{s[:30]}'", mpt.count(s) == 0, mpt.count(s))

    spt = pdf_text(SI_PDF)
    for s in ["Supplementary Note 2", "Replication within one tumour type",
              "13.1-, 14.3- and 19.9-fold", "0.549 ± 0.026, 0.555 and 0.629 ± 0.031",
              "0.618, 0.428 and 0.440", "12,148", "0.658", "0.686", "0.760",
              "27-fold", "2- to 20-fold", "Colon cancer, second patient",
              "Colon cancer, third patient"]:
        ck(f"SI PDF 含 '{s[:36]}'", spt.count(s) >= 1, spt.count(s))
    for s in ["seven span 21-fold", "same tumour block", "All four Visium"]:
        ck(f"SI PDF 无旧串 '{s[:30]}'", spt.count(s) == 0, spt.count(s))

    cpt = pdf_text(COVER_PDF)
    ct = flat(docx_text(COVER))
    for s in ["2- to 20-fold", "eight Visium HD sections", "19.9-fold",
              "2- to 20-fold range", "0.13 to 0.63",
              "Three sections of one tumour separate along two axes"]:
        ck(f"投稿信含 '{s[:36]}'", (s in ct and s in cpt),
           f"docx={ct.count(s)} pdf={cpt.count(s)}")
    for s in ["2- to 15-fold", "seven Visium HD sections", "Two sections of one tumour agree"]:
        ck(f"投稿信无旧串 '{s[:30]}'", ct.count(s) == 0 and cpt.count(s) == 0,
           f"{ct.count(s)}/{cpt.count(s)}")

    # ---------------------------------------------------------- 6 图件与内嵌副本
    fj = json.loads((REPRO / "07_SUBMISSION" / "_build_v57" / "audit" /
                     "fig_v57.json").read_text("utf-8"))
    ck("Fig 1b/S3 曲线数 = 8", fj["n_sections"] == 8, fj["n_sections"])
    ck("Fig1b 端点 fold 含 19.863", abs(float(fj["fold"]["colon_p5"]) - 19.863) < 5e-3,
       fj["fold"]["colon_p5"])
    ck("Fig S3 2 µm 跨 26.9 倍", abs(fj["fig_s3"]["ratio_2um"] - 26.94) < 0.05,
       fj["fig_s3"]["ratio_2um"])
    ck("Fig S3 256 µm 跨 3.8 倍", abs(fj["fig_s3"]["ratio_256um"] - 3.83) < 0.05,
       fj["fig_s3"]["ratio_256um"])

    for doc, name in ((MAN, "Fig1"), (SI, "FigS3")):
        d = docx.Document(str(doc))
        want = sha(FIGD / f"{name}.png")
        hit = any(sha_bytes == want for sha_bytes in
                  (hashlib.sha256(p.blob).hexdigest()
                   for p in d.part.package.iter_parts()
                   if p.content_type.startswith("image/")))
        ck(f"{doc.name} 内嵌 {name} == 04_图件/{name}.png", hit, want[:16])
    ck("04_图件 Fig1.pdf/png 存在", (FIGD / "Fig1.pdf").exists()
       and (FIGD / "Fig1.png").exists(), "")
    ck("04_图件 FigS3.pdf/png 存在", (FIGD / "FigS3.pdf").exists()
       and (FIGD / "FigS3.png").exists(), "")

    # ---------------------------------------------------------- 7 Xenium 一节仍在
    PRED = Path(r"E:\xen_p2\out\xen_predict_closing10.json")
    pred = json.loads(PRED.read_text("utf-8"))
    rows = pred["rows"]
    dsv = [r["d"] for r in rows]
    pv = [r["alpha_pred"] for r in rows]
    ov = [r["alpha_true"] for r in rows]
    x = [math.log(v) for v in dsv]
    xm = sum(x) / len(x)

    def slope(yv):
        ym = sum(yv) / len(yv)
        return (sum((a - xm) * (b - ym) for a, b in zip(x, yv))
                / sum((a - xm) ** 2 for a in x))
    s_pred = slope([math.log(v) for v in pv])
    s_obs = slope([math.log(v) for v in ov])
    ck("Note 1 斜率 predicted −0.343", abs(s_pred + 0.34346) < 5e-4, f"{s_pred:.5f}")
    ck("Note 1 斜率 observed −0.024", abs(s_obs + 0.02447) < 5e-4, f"{s_obs:.5f}")
    i2 = dsv.index(2.0)
    sc = ov[i2] / pv[i2]
    devs = [0.0 if d == 2.0 else math.log10(abs(p * sc / o))
            for d, p, o in zip(dsv, pv, ov)]
    mx = max(abs(v) for v in devs)
    ck("Note 1 最大 log10 偏差 0.831", abs(mx - 0.83122) < 5e-4, f"{mx:.5f}")

    t9 = ds.tables[6]
    cell_bad = []
    for i, r in enumerate(rows, start=1):
        got = [c.text.strip() for c in t9.rows[i].cells]
        dv = devs[i - 1]
        want = [f"{r['d']:g}", f"{r['alpha_pred']:.3f}", f"{r['alpha_true']:.3f}",
                "0.000" if abs(dv) < 5e-4 else f"{dv:+.3f}".replace("-", "\u2212")]
        if want != got:
            cell_bad.append((want, got))
    ck("Table S9 九行与 JSON 复算一致", not cell_bad, str(cell_bad[:1]))

    # ---------------------------------------------------------- 8 输出
    fails = [c for c in CHECKS if not c[1]]
    w = max(len(c[0]) for c in CHECKS)
    for name, ok, detail in CHECKS:
        print(f"  {'PASS' if ok else 'FAIL'}  {name:<{w}}  {detail}")
    print(f"\n{len(CHECKS) - len(fails)}/{len(CHECKS)} PASS")

    rec = dict(generator="111_guard_v57.py", n_checks=len(CHECKS),
               n_pass=len(CHECKS) - len(fails), pass_all=not fails,
               files={str(p.relative_to(PKG)): sha(p) for p in
                      (MAN, SI, COVER, MAN_PDF, SI_PDF, COVER_PDF)},
               checks=[dict(name=n, ok=o, detail=d) for n, o, d in CHECKS])
    out = PKG / "05_Audit" / "guard_v5_7.json"
    out.write_text(json.dumps(rec, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"wrote {out}")
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())

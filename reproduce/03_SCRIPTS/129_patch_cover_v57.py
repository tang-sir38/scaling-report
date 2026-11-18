# -*- coding: utf-8 -*-
"""129_patch_cover_v57.py — 投稿信与投稿指南同步到 v5.7。

投稿信只有三处：标题区间、七片→八片那句、以及"两片一致"改成"两轴分离"。
投稿指南是内部手册（不上传），改状态表、门户 Title、checklist、攻击表，
并补一节 "v5.7 相对 v5.6 改了什么"。
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import docx

PKG = Path(r"C:\Users\唐朝华\Desktop\23空间分辨率\A题_空间分辨率_投稿包_v5.7_20261118")
COVER = PKG / "01_投稿信" / "00_Cover_Letter.docx"
GUIDE = PKG / "00_投稿指南.md"

EDITS = [
    ("Overdispersion in binned spatial transcriptomics falls 2- to 15-fold, "
     "not 16,384-fold",
     "Overdispersion in binned spatial transcriptomics falls 2- to 20-fold, "
     "not 16,384-fold",
     "标题"),
    ("Across seven Visium HD sections, one per tissue, and one Stereo-seq section, "
     "pooled overdispersion between 2 and 256 \u00b5m falls 2.1-fold in human pancreas "
     "and 15.6-fold in mouse small intestine \u2014 a 2- to 15-fold range \u2014 against "
     "the 16,384-fold the counting argument predicts; the exponents run 0.13 to 0.55.",
     "Across eight Visium HD sections \u2014 seven tissues, colon captured three times "
     "\u2014 and one Stereo-seq section, pooled overdispersion between 2 and 256 \u00b5m "
     "falls 2.1-fold in human pancreas and 19.9-fold in a third colon tumour \u2014 a "
     "2- to 20-fold range \u2014 against the 16,384-fold the counting argument "
     "predicts; the exponents run 0.13 to 0.63.",
     "正文段落：七片 → 八片、区间、指数范围"),
    ("Two sections of one tumour agree; one tissue measured on two chemistries does "
     "not; and the same statistic ranks three tissues in opposite orders when pooled "
     "and when taken per gene.",
     "Three sections of one tumour separate along two axes: one moves the pooled "
     "exponent by 0.006 and the per-gene median by 0.19, another the reverse; one "
     "tissue on two chemistries does not agree; and the same statistic ranks three "
     "tissues in opposite orders when pooled and when taken per gene.",
     "\u201c两片一致\u201d → \u201c两轴分离\u201d"),
]

GUIDE_EDITS = [
    ("# 投稿指南 · A 题（空间分辨率 / 分箱标度律）v5.6",
     "# 投稿指南 · A 题（空间分辨率 / 分箱标度律）v5.7", "标题"),
    ("版本 v5.6 ｜ 2026-11-17 ｜ 目标刊", "版本 v5.7 ｜ 2026-11-18 ｜ 目标刊", "日期"),
    ("| 正文字数（Introduction→Discussion，含小节标题） | 3,496 | 3,500 | 余 4 |",
     "| 正文字数（Introduction→Discussion，含小节标题） | 3,499 | 3,500 | 余 1 |",
     "字数"),
    ("| 补充表 | S1–S9（**S9 为 v5.6 新增**） | 不限 | — |",
     "| 补充表 | S1–S9（S9 为 v5.6 新增；S1 为 11 行，含 P1/P5） | 不限 | — |", "补充表"),
    ("| 分析切片 | **7 片进入标度分析**（包内共 8 个 Visium HD 切片 + 1 个 Stereo-seq） | — | — |",
     "| 分析切片 | **8 片入图**（7 个组织 + 第三片结肠；包内共 9 个 Visium HD 切片 + 1 个 Stereo-seq，P1 只作重复不入图） | — | — |",
     "切片数"),
    ("| 包内文件 | 37 | — | — |", "| 包内文件 | 40（以 MANIFEST.json 为准） | — | — |", "文件数"),
    ("守卫状态：**59/59 PASS，0 失败**，见 `05_Audit/guard_v5_6.json`",
     "守卫状态：**见 `05_Audit/guard_v5_7.json`**（v5.6 为 59/59，见 `guard_v5_6.json`）",
     "守卫"),
    ("| Title | Overdispersion in binned spatial transcriptomics falls 2- to 15-fold, "
     "not 16,384-fold, because it is set by tissue covariance | 与稿件逐字一致 |",
     "| Title | Overdispersion in binned spatial transcriptomics falls 2- to 20-fold, "
     "not 16,384-fold, because it is set by tissue covariance | 与稿件逐字一致 |",
     "门户 Title"),
    ("| \"只加了几个组织，广度还是不够\" | 7 片、7 种组织、三类状态（两种癌、四种正常、一种发育中）。命题的单位不是组织而是相关结构；每片都有 Table S1 的来源与格数 |",
     "| \"只加了几个组织，广度还是不够\" | 8 片入图、7 种组织、三类状态（三种癌、四种正常、一种发育中），且结肠有三片不同患者可作组内对照。命题的单位不是组织而是相关结构；每片都有 Table S1 的来源与格数 |",
     "攻击表-广度"),
    ("| \"组织是挑出来的，正好凑出 2–15\" | 四片在算指数前定下（Methods、Fig. S3 图注各写一次）；且最大一片的 CI 与结肠重合，本文没有用它去主张任何排序 |",
     "| \"组织是挑出来的，正好凑出 2–20\" | 四片在算指数前定下（Methods、Fig. S3 图注各写一次）；P1/P5 是同一数据集里的重复片，加入后最大一片（P5，19.9）的 CI 与 P2（14.3）不重叠，本文主动写出这一分歧，没有用它去主张任何排序 |",
     "攻击表-挑组织"),
    ("- [ ] 系统里 Title 是 `2- to 15-fold`（**不是**旧版的 `5- to 14-fold`）",
     "- [ ] 系统里 Title 是 `2- to 20-fold`（**不是**旧版的 `2- to 15-fold`）", "checklist"),
    ("- [ ] 补充材料上传后翻到最后：**Figure S3 在**，S1–S9 与 Note 1 都在",
     "- [ ] 补充材料上传后翻到最后：**Figure S3 在**，S1–S9 与 Note 1、Note 2 都在（Note 2 是 v5.7 新增，第 8–9 页）",
     "checklist-Note2"),
    ("- [ ] Table S1 是 **10 行**（少了就是拿成 v5.4 的包）",
     "- [ ] Table S1 是 **11 行**（含 P1、P5 两行重复片；少了就是拿成 v5.6 的包）",
     "checklist-表"),
]


def replace_once(text: str, old: str, new: str, where: str) -> str:
    n = text.count(old)
    if n != 1:
        raise SystemExit(f"{where}: 出现 {n} 次（要求 1 次）\n  {old[:90]!r}")
    return text.replace(old, new)


def patch_cover():
    d = docx.Document(str(COVER))
    for old, new, why in EDITS:
        hit = False
        for p in d.paragraphs:
            if old in p.text:
                if len(p.runs) == 1:
                    p.runs[0].text = replace_once(p.runs[0].text, old, new, why)
                    hit = True
                    break
                # 跨 run：把整段文字压到第一个 run（投稿信全是同一格式）
                full = p.text
                if full.count(old) != 1:
                    raise SystemExit(f"{why}: 段内出现 {full.count(old)} 次")
                p.runs[0].text = full.replace(old, new)
                for r in p.runs[1:]:
                    r.text = ""
                hit = True
                break
        if not hit:
            raise SystemExit(f"{why}: 找不到\n  {old[:80]!r}")
    total = sum(len(p.text.split()) for p in d.paragraphs)
    d.save(str(COVER))
    return dict(words=total, edits=[w for _, _, w in EDITS])


def patch_guide():
    s = GUIDE.read_text(encoding="utf-8")
    for old, new, why in GUIDE_EDITS:
        s = replace_once(s, old, new, why)
    s = s.replace(
        "---\n\n## 零之一、v5.6 相对 v5.5.1 改了什么（只有 SI 一处 + 主稿两句）",
        """---

## 零之一、v5.7 相对 v5.6 改了什么（第三片结肠 + 一处身份更正）

1. **P1 身份更正（事实错误）**。v5.6 写 "a second section of the same tumour block"，
   是错的。P1 是另一个患者：其 2 µm 矩阵里 DDX3Y/KDM5D/EIF1AY/UTY/USP9Y/NLGN4Y
   六个 Y 基因的 per-bin 计数全部为 0，而 P2 是男性。10x 的样本命名
   （P1/P2/P3/P5）也是患者编号。已在正文 [67]、[46] 与 Table S1 改正。
2. **第三片结肠癌 P5 入稿**。pooled 2→256 µm 下降 **19.86 倍**（八片最大），
   八片区间 2.12–19.86 → 标题 `2- to 20-fold`，摘要 `2.1- to 19.9-fold`，
   Fig 1b 与 Fig S3 由七条曲线改成八条（同组织两片用同色、实心/空心区分）。
3. **双向解耦（新的实证）**。三片结肠、三个患者、同平台同流程：
   - P1 vs P2：pooled γ 差 0.006，per-gene 中位差 0.190
   - P5 vs P2：pooled γ 差 0.074，per-gene 中位差 0.012
   两个方向都解耦，比原来"两片一致"更有力。正文 [46] 重写，SI 新增 Note 2。
4. **四舍五入统一**：colon 的 γ 0.554 → 0.555（[20]、[25]；与 [45]、SI 一致）。

v5.6 相对 v5.5.1 的改动见下一节（保留作历史）。

---

## 零之二、v5.6 相对 v5.5.1 改了什么（只有 SI 一处 + 主稿两句）""", 1)
    GUIDE.write_text(s, encoding="utf-8")
    return dict(chars=len(s))


def main() -> int:
    c = patch_cover()
    g = patch_guide()
    print(f"cover: {c['words']} words")
    print(f"guide: {g['chars']} chars")
    return 0


if __name__ == "__main__":
    sys.exit(main())

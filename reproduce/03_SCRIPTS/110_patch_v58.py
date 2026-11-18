r"""110_patch_v58.py — 在 v5.7 交付包上做 v5.8 增量修订（不重建、只改差异）。

七处硬伤（v5.7 遗留）
--------------------
1 [23] 小标题 "at most fifteen-fold"         → twenty-fold（与标题 20 倍冲突）
2 [55] 讨论 "two- to fifteen-fold"            → two- to twenty-fold
3 [13] 摘要 "α ∝ d−0.13 to d−0.55"            → d−0.63（P5 拟合 0.629 已是上界）
4 [67] 方法 "其余七片进 scaling"              → 与 Fig 1 九条曲线一致
5 SI   Note 2 "0.1 transcripts per bin / 0.001" → 到 10% counts / 0.0004
6 [63] "a halving of depth by 0.0004"        → 0.00004（0.0004 是 10% 那次）
7 投稿信 "eight sections — seven tissues, colon ×3" → nine（7 组织含结肠×3 = 9 片）

连带同步（八片 → 九片）
----------------------
摘要、引言、Results、Methods、Fig 1 图注、Fig S3 图注、投稿信、投稿指南。

数字一律现场注入或来自审计文件，不手抄。
"""
from __future__ import annotations

import json
from pathlib import Path

import docx

PKG = Path(r"C:\Users\唐朝华\Desktop\23空间分辨率\A题_空间分辨率_投稿包_v5.8_20261118")
AUD = Path(r"C:\Users\唐朝华\Desktop\scaling-report\reproduce\07_SUBMISSION\_build_v58\audit\fig_v58.json")

MAN = PKG / "02_主稿" / "01_Manuscript.docx"
SI = PKG / "03_补充材料" / "02_Supplementary_Information.docx"
COVER = PKG / "01_投稿信" / "00_Cover_Letter.docx"
GUIDE = PKG / "00_投稿指南.md"


def replace_in_para(p, old: str, new: str) -> bool:
    """跨 run 的精确替换：只重写受影响的那几个 run，保留其余格式。"""
    runs = p.runs
    text = "".join(r.text for r in runs)
    idx = text.find(old)
    if idx < 0:
        return False
    end = idx + len(old)
    pos = 0
    spans = []
    for r in runs:
        spans.append((pos, pos + len(r.text)))
        pos += len(r.text)

    first = next(i for i, (a, b) in enumerate(spans) if b > idx)
    last = next(i for i, (a, b) in enumerate(spans) if b >= end)
    a0, _ = spans[first]
    _, b1 = spans[last]
    if first == last:
        r = runs[first]
        r.text = r.text[: idx - a0] + new + r.text[end - a0:]
    else:
        r0 = runs[first]
        r0.text = r0.text[: idx - a0] + new
        for i in range(first + 1, last):
            runs[i].text = ""
        rl = runs[last]
        rl.text = rl.text[end - spans[last][0]:]
    return True


def apply(doc, edits, tag):
    for label, old, new in edits:
        hits = [p for p in doc.paragraphs if old in "".join(r.text for r in p.runs)]
        if len(hits) != 1:
            raise SystemExit(f"[{tag}] {label}: 命中 {len(hits)} 段，应为 1\n  old={old[:80]!r}")
        if not replace_in_para(hits[0], old, new):
            raise SystemExit(f"[{tag}] {label}: 替换失败")
        print(f"  [{tag}] {label}")
    for tb in doc.tables:
        for row in tb.rows:
            for c in row.cells:
                for p in c.paragraphs:
                    for label, old, new in edits:
                        if old in "".join(r.text for r in p.runs):
                            raise SystemExit(f"[{tag}] {label} 也命中了表格，需人工确认")


def man_edits(aud):
    return [
        ("摘要 八→九",
         "the native 2 µm bins of eight Visium HD sections",
         "the native 2 µm bins of nine Visium HD sections"),
        ("摘要 指数上界 0.55→0.63",
         "with α ∝ d−0.13 to d−0.55",
         "with α ∝ d−0.13 to d−0.63"),
        ("引言 五→六片",
         "plus five further sections (normal pancreas, small intestine, kidney, "
         "embryo, and a third colon tumour)",
         "plus six further sections (normal pancreas, small intestine, kidney, "
         "embryo, and two further colon tumours)"),
        ("Results 小标题 fifteen→twenty",
         "Overdispersion falls at most fifteen-fold, not sixteen thousand-fold",
         "Overdispersion falls at most twenty-fold, not sixteen thousand-fold"),
        ("Results 五→六片并补 13.1",
         "Five further sections fall 15.6- (small intestine), 3.5- (kidney), "
         "2.1- (pancreas), 2.3- (embryo) and 19.9-fold (the third colon tumour). "
         "Its bootstrap interval, 16.2 to 23.4, does not overlap the 14.3 above.",
         "Six further sections fall 19.9- and 13.1- (two further colon tumours), "
         "15.6- (small intestine), 3.5- (kidney), 2.3- (embryo) and 2.1-fold "
         "(pancreas). The largest has a bootstrap interval of 16.2 to 23.4, "
         "which does not overlap the 14.3 above."),
        ("讨论 two- to fifteen→twenty",
         "overdispersion falls two- to fifteen-fold rather than 16,384-fold",
         "overdispersion falls two- to twenty-fold rather than 16,384-fold"),
        ("讨论 halving 的数值",
         "a halving of depth by 0.0004",
         "a halving of depth by 0.00004"),
        ("方法 七片进 scaling→九片全进",
         "The two further colon sections are replicates; the other seven enter "
         "the scaling analysis, one per tissue.",
         "Three are colon adenocarcinomas from different patients; the other six "
         "tissues contribute one section each, and all nine enter the scaling "
         "analysis."),
        ("Fig1 注 八→九 (a)",
         "for all eight sections; the number at the right end",
         "for all nine sections; the number at the right end"),
        ("Fig1 注 八→九 (b)",
         "the eight curves would lie on the steep grey line",
         "the nine curves would lie on the steep grey line"),
    ]


def si_edits(aud):
    return [
        ("FigS3 注 八→九",
         "α(d) in absolute units for all eight sections, on the same log–log axes "
         "and the same eight scales.",
         "α(d) in absolute units for all nine sections, on the same log–log axes "
         "and the same eight bin sizes."),
        ("FigS3 注 跨度主语",
         "At 2 µm the eight span 27-fold",
         "At 2 µm the nine span 27-fold"),
        ("Note 2 稀释口径与数值",
         "Binomial thinning of the P2 field to 0.1 transcripts per bin moves the "
         "pooled exponent by 0.001 (Fig. 3d)",
         "Binomial thinning of the P2 field to a tenth of its counts moves the "
         "pooled exponent by 0.0004 (Methods)"),
    ]


def cover_edits(aud):
    return [
        ("投稿信 八→九",
         "Across eight Visium HD sections — seven tissues, colon captured three times",
         "Across nine Visium HD sections — seven tissues, colon captured three times"),
    ]


def main() -> int:
    aud = json.loads(AUD.read_text(encoding="utf-8"))
    assert aud["n_sections"] == 9, aud["n_sections"]

    print("manuscript")
    d = docx.Document(str(MAN))
    before = len(" ".join(p.text for p in d.paragraphs).split())
    apply(d, man_edits(aud), "MAN")
    d.save(str(MAN))
    d2 = docx.Document(str(MAN))
    after = len(" ".join(p.text for p in d2.paragraphs).split())
    print(f"  词数 {before} → {after}（正文口径 15:65 另行统计）")

    print("supplementary")
    s = docx.Document(str(SI))
    apply(s, si_edits(aud), "SI")
    s.save(str(SI))

    print("cover letter")
    c = docx.Document(str(COVER))
    apply(c, cover_edits(aud), "COVER")
    c.save(str(COVER))

    print("guide")
    g = GUIDE.read_text(encoding="utf-8")
    for label, old, new in [
        ("指南 八片最大→九片", "（八片最大）", "（九片最大）"),
        ("指南 八片区间→九片", "八片区间 2.12–19.86", "九片区间 2.12–19.86"),
        ("指南 8 片入图→9 片", "8 片入图、7 种组织", "9 片入图、7 种组织"),
    ]:
        n = g.count(old)
        if n != 1:
            raise SystemExit(f"[GUIDE] {label}: 命中 {n} 次，应为 1")
        g = g.replace(old, new)
        print(f"  [GUIDE] {label}")
    GUIDE.write_text(g, encoding="utf-8")

    print("\n回读校验")
    blob_m = " ".join(p.text for p in docx.Document(str(MAN)).paragraphs)
    blob_s = " ".join(p.text for p in docx.Document(str(SI)).paragraphs)
    blob_c = " ".join(p.text for p in docx.Document(str(COVER)).paragraphs)
    must = [
        (blob_m, "nine Visium HD sections", 1),
        (blob_m, "d−0.13 to d−0.63", 1),
        (blob_m, "at most twenty-fold", 1),
        (blob_m, "two- to twenty-fold", 1),
        (blob_m, "Six further sections", 1),
        (blob_m, "all nine enter the scaling analysis", 1),
        (blob_m, "the nine curves would lie", 1),
        (blob_m, "a halving of depth by 0.00004", 1),
        (blob_s, "for all nine sections", 1),
        (blob_s, "the nine span 27-fold", 1),
        (blob_s, "a tenth of its counts", 1),
        (blob_c, "Across nine Visium HD sections", 1),
    ]
    for blob, s_, k in must:
        c_ = blob.count(s_)
        print(f"  {'OK ' if c_ == k else 'FAIL'} {s_[:46]:48s} {c_}")
    forbidden = ["fifteen-fold", "2- to 15-fold", "five further sections",
                 "eight Visium HD sections", "the eight curves", "the eight span",
                 "0.1 transcripts per bin", "by 0.001 (Fig. 3d)"]
    for f in forbidden:
        cm, cs = blob_m.count(f), blob_s.count(f)
        print(f"  {'OK ' if cm == 0 and cs == 0 else 'FAIL'} 已清除 {f!r}  (MAN {cm} / SI {cs})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

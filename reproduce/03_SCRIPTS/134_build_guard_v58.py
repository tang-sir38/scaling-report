r"""134_build_guard_v58.py — 由 v5.7 守卫生成 v5.8 守卫（保留 135 项，改差异项）。

为什么要"生成"而不是手抄
------------------------
v5.7 守卫 135 项全部通过，却漏掉了 "fifteen-fold"（换了个词形）。这次除了把
受影响的断言改成新值，还加了一类**黑名单**检查：任何旧词形只要出现就 FAIL，
不再依赖"我记得它可能怎么写"。

所有替换都断言命中次数，命中数不对就停，不生成半成品。
"""
from __future__ import annotations

from pathlib import Path

SRC = Path(__file__).resolve().parent / "111_guard_v57.py"
DST = Path(__file__).resolve().parent / "111_guard_v58.py"

REPL = [
    # ---- 路径 / 版本 ----
    ("A题_空间分辨率_投稿包_v5.7_20261118", "A题_空间分辨率_投稿包_v5.8_20261118", 1),
    ('generator="111_guard_v57.py"', 'generator="111_guard_v58.py"', 1),
    ('guard_v5_7.json', 'guard_v5_8.json', 1),
    ('_build_v57" / "audit" /\n                     "fig_v57.json"',
     '_build_v58" / "audit" /\n                     "fig_v58.json"', 1),

    # ---- 主稿断言：八→九，五→六 ----
    ('"eight Visium HD sections from human and mouse", 1, "摘要片数"',
     '"nine Visium HD sections from human and mouse", 1, "摘要片数"', 1),
    ('"plus five further sections", 1, "[20] 引言"',
     '"plus six further sections", 1, "[20] 引言"', 1),
    ('"and a third colon tumour", 1, "[20]"',
     '"and two further colon tumours", 1, "[20]"', 1),
    ('"Five further sections fall 15.6-", 1, "[25]"',
     '"Six further sections fall 15.6-", 1, "[25]"', 1),
    ('"3.5- (kidney), 2.1- (pancreas), 2.3- (embryo)", 1, "[25] 另三片"',
     '"3.5- (kidney), 2.3- (embryo), 2.1- (pancreas)", 1, "[25] 另三片"', 1),
    ('"19.9-fold (the third colon tumour)", 1, "[25]"',
     '"13.1- and 19.9-fold (two further colon tumours)", 1, "[25]"', 1),
    ('"Its bootstrap interval, 16.2 to 23.4", 1, "[25]"',
     '"The 19.9-fold section\'s interval, 16.2 to 23.4", 1, "[25]"', 1),
    ('"Two further colon tumours from different patients", 1, "[46]"',
     '"Two further colon tumours fall 13.1- and 19.9-fold", 1, "[46]"', 1),
    ('"The two further colon sections are replicates", 1, "[67]"',
     '"Three are colon adenocarcinomas from different patients", 1, "[67]"', 1),
    ('"for all eight sections", 1, "Fig 1 注"',
     '"for all nine sections", 1, "Fig 1 注"', 1),
    ('"the eight curves would lie", 1, "Fig 1 注"',
     '"the nine curves would lie", 1, "Fig 1 注"', 1),
    # ---- SI 断言 ----
    ('"for all eight sections", 1, "Fig S3 注"',
     '"for all nine sections", 1, "Fig S3 注"', 1),
    ('"At 2 µm the eight span 27-fold", 1, "Fig S3 注"',
     '"At 2 µm the nine span 27-fold", 1, "Fig S3 注"', 1),
    # ---- 图件审计 ----
    ('ck("Fig 1b/S3 曲线数 = 8", fj["n_sections"] == 8, fj["n_sections"])',
     'ck("Fig 1b/S3 曲线数 = 9", fj["n_sections"] == 9, fj["n_sections"])', 1),
    ('abs(fj["fig_s3"]["ratio_256um"] - 3.83) < 0.05',
     'abs(fj["fig_s3"]["ratio_256um"] - 3.848) < 5e-3', 1),
    # ---- PDF 必含串 ----
    ('"2- to 20-fold", "2.1- to 19.9-fold", "eight Visium HD sections",\n'
     '              "Five further sections fall 15.6-", "19.9-fold (the third colon tumour)",',
     '"2- to 20-fold", "2.1- to 19.9-fold", "nine Visium HD sections",\n'
     '              "Six further sections fall 15.6-",\n'
     '              "13.1- and 19.9-fold (two further colon tumours)",', 1),
    ('"19.9-fold (colon P5)", "eight curves would lie"',
     '"19.9-fold (colon P5)", "nine curves would lie"', 1),
    # 手稿带行号，长句会被行号切断（"16.2 79 to 23.4"），PDF 级只查短串
    ('"16.2 to 23.4", "Reproducibility splits in both directions"',
     '"19.9-fold section", "Reproducibility splits in both directions"', 1),
    # ---- 投稿信 ----
    ('"2- to 20-fold", "eight Visium HD sections", "19.9-fold",',
     '"2- to 20-fold", "nine Visium HD sections", "19.9-fold",', 1),
]

# 新增检查：插在 "主稿 PDF 含旧串" 那段之后
EXTRA_PDF_BLACK = '''    for s in ["fifteen-fold", "the eight curves", "five further sections",
              "0.1 transcripts per bin", "a bootstrap interval of 16.2"]:
        ck(f"主稿 PDF 无旧词形 '{s[:30]}'", mpt.count(s) == 0, mpt.count(s))
'''

# docx 级黑名单：v5.7 就是只查了 PDF 与精确串，才让 "fifteen-fold" 活下来
MAN_BLACK = '''    for s in ["fifteen-fold", "the eight curves", "the eight span",
              "five further sections", "eight Visium HD sections",
              "0.1 transcripts per bin", "2- to 15-fold"]:
        ck(f"主稿 docx 无旧词形 '{s[:28]}'", mt.count(s) == 0, mt.count(s))
'''

SI_BLACK = '''    for s in ["fifteen-fold", "the eight span", "0.1 transcripts per bin",
              "a bootstrap interval of 16.2"]:
        ck(f"SI docx 无旧词形 '{s[:28]}'", st.count(s) == 0, st.count(s))
'''

EXTRA_FIG = '''    ck("Fig1b 端点 fold 含 13.097（P1）",
       abs(float(fj["fold"]["colon_p1"]) - 13.097) < 5e-3, fj["fold"]["colon_p1"])
'''

EXTRA_NOTE2 = '''        ("a tenth of its counts moves the pooled exponent by 0.0004", 1,
         "Note 2 稀释口径"),
'''


def main() -> int:
    s = SRC.read_text(encoding="utf-8")
    for old, new, n in REPL:
        c = s.count(old)
        if c != n:
            raise SystemExit(f"命中 {c} 次（应 {n}）:\n  {old[:90]!r}")
        s = s.replace(old, new)

    anchor = '    for s in ["2- to 15-fold", "2.1- to 15.6-fold", "seven Visium HD sections",'
    if s.count(anchor) != 1:
        raise SystemExit("找不到 PDF 旧串段锚点")
    s = s.replace(anchor, EXTRA_PDF_BLACK + anchor)

    # v5.8 新要求：小标题必须是 twenty-fold（v5.7 这处没有 docx 检查）
    a23 = '        ("falls 2- to 20-fold, not 16,384-fold", 1, "标题区间"),'
    if s.count(a23) != 1:
        raise SystemExit("找不到标题断言锚点")
    s = s.replace(a23, a23 + '\n        ("at most twenty-fold", 1, "[23] 小标题"),')

    # 两处 docx 黑名单：插在各自 "c = ..." 循环之后
    a_mt = ('        c = mt.count(s)\n'
            '        ck(f"主稿 \'{s[:42]}\' ×{n}", c == n, f"{c}  [{why}]")\n')
    if s.count(a_mt) != 1:
        raise SystemExit("找不到主稿循环锚点")
    s = s.replace(a_mt, a_mt + MAN_BLACK)

    a_st = ('        c = st.count(s)\n'
            '        ck(f"SI \'{s[:42]}\' ×{n}", c == n, f"{c}  [{why}]")\n')
    if s.count(a_st) != 1:
        raise SystemExit("找不到 SI 循环锚点")
    s = s.replace(a_st, a_st + SI_BLACK)

    anchor2 = '    ck("Fig S3 2 µm 跨 26.9 倍"'
    if s.count(anchor2) != 1:
        raise SystemExit("找不到 Fig S3 锚点")
    s = s.replace(anchor2, EXTRA_FIG + anchor2)

    anchor3 = '        # 旧串\n        ("a second section of the same tumour block", 0, "旧 P1 身份"),'
    if s.count(anchor3) != 1:
        raise SystemExit("找不到 SI 旧串锚点")
    s = s.replace(anchor3, EXTRA_NOTE2 + anchor3)

    DST.write_text(s, encoding="utf-8")
    print(f"written {DST.name}  {len(s.splitlines())} lines")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

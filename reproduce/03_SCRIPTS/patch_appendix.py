# -*- coding: utf-8 -*-
"""S1_Appendix v7.3 -> v8.1：脚本化补丁，逐条断言唯一命中。

只改段落文本，不动表格。改动原因见 07_审核报告：
  A. 附录标题仍是旧稿标题（v7.3 标题已被 v8 替换）
  B. Contents 里 S1 Fig / S3 Fig 与正文图注互换（正文引用的是图注，故 Contents 错）
  C. P2 的 γ 在附录写 0.555，主稿写 0.554：两条流水线基因集不同（18,074 vs 18,061），
     统一到主分析值并注明重复流水线的返回值
"""
import sys, docx, json

SRC = r'C:\Users\唐朝华\WorkBuddy\A73_rewrite\_src\03_补充材料\S1_Appendix.docx'
DST = r'C:\Users\唐朝华\WorkBuddy\A73_rewrite\v8\03_S1_Appendix_v8.docx'
N = json.load(open(r'C:\Users\唐朝华\WorkBuddy\A73_rewrite\v8\numbers_v8.json', encoding='utf-8'))

TITLE_OLD = ('Overdispersion in binned spatial transcriptomics falls 2- to 20-fold, '
             'not 16,384-fold, because the exponent is the decay rate of tissue covariance')
TITLE_NEW = ('Binning in spatial transcriptomics removes two- to twenty-fold of the excess '
             'variance, not the sixteen-thousand-fold that counting statistics predict')

REPS = [
    # A 标题
    (TITLE_OLD, TITLE_NEW),
    # B Contents 两行互换（先改成占位，避免互相覆盖）
    ('S1 Fig. The two failed curve predictions.', '@@FIG_A@@'),
    ('S3 Fig. The same curves on an absolute scale.', 'S3 Fig. The two failed curve predictions.'),
    ('@@FIG_A@@', 'S1 Fig. The same curves on an absolute scale.'),
    # C P2 的 γ
    ('with exponents 0.549 ± 0.026, 0.555 and 0.629 ± 0.031',
     'with exponents 0.549 ± 0.026, 0.554 and 0.629 ± 0.031'),
    # D 链接复核日期：本轮 10/10 全部 200（2026-11-24），据实更新
    ('All ten URLs returned HTTP 200 on 20 November 2026.',
     'All ten URLs returned HTTP 200 on 24 November 2026.'),
    # E 数字源文件名：主稿写的是 numbers_v8.json，附录不能还写 v7 的锁名
    ('Every number in the main text is read from the numbers lock file (NUMBERS_LOCK_v7.json) that ships with the code.',
     'Every number in the main text is read from the numbers file (numbers_v8.json) that ships with the code.'),
]

TAIL_ANCHOR = ('Pooled, the three fall 13.1-, 14.3- and 19.9-fold between 2 and 256 µm')
TAIL_ADD = (' P2 is quoted at 0.554 from the main analysis; recomputed through the replication '
            'pipeline, which applies its own gene filter, it returns 0.555.')

d = docx.Document(SRC)
paras = d.paragraphs
hit = {k: 0 for k, _ in REPS}
tail = 0
for p in paras:
    t = p.text
    if not t.strip():
        continue
    for old, new in REPS:
        if old in t:
            hit[old] += 1
            if old.startswith('@@') or new.startswith('@@'):
                p.text = t.replace(old, new)
                t = p.text
            else:
                p.text = t.replace(old, new)
                t = p.text
    if TAIL_ANCHOR in t and '0.554 from the main analysis' not in t:
        p.text = t + TAIL_ADD
        tail += 1

for old, _ in REPS:
    if hit[old] != 1:
        sys.exit(f'FATAL: 断言失败 {old[:40]!r} 命中 {hit[old]} 次')
if tail != 1:
    sys.exit(f'FATAL: 尾注锚点命中 {tail} 次')

d.save(DST)
print('saved', DST)
for p in d.paragraphs:
    if p.text.strip().startswith(('Binning in spatial', 'S1 Fig.', 'S3 Fig.')):
        print(' >', p.text[:110])

"""build_v8.py — 从 MS_v8.md 生成主稿 docx。
- 引用 [[CITE:key]] 按首现序自动编号（禁手写编号）
- 数字 {{key}} / {{key:.3f}} 从 numbers_v8.json 注入，缺键或引用未定义即中止
- 版式：Letter / 2.54cm / TNR 12 / 双倍行距 / 连续行号 / 页脚页码
"""
import json, re, sys, os
from docx import Document
from docx.shared import Pt, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.enum.section import WD_SECTION
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

BASE = os.path.dirname(os.path.abspath(__file__))
MD = os.path.join(BASE, 'MS_v8.md')
N = json.load(open(os.path.join(BASE, 'numbers_v8.json'), encoding='utf-8'))
REFS = json.load(open(os.path.join(BASE, 'refs_v8.json'), encoding='utf-8'))

src = open(MD, encoding='utf-8').read()

# ---------- 1. 引用编号：按首现序 ----------
order = []
for m in re.finditer(r'\[\[CITE:([a-zA-Z0-9_]+)\]\]', src):
    k = m.group(1)
    if k not in REFS:
        sys.exit(f'FATAL: undefined reference key {k}')
    if k not in order:
        order.append(k)
if len(order) != len(REFS):
    missing = [k for k in REFS if k not in order]
    sys.exit(f'FATAL: {len(missing)} reference(s) never cited: {missing}')
num = {k: i + 1 for i, k in enumerate(order)}

# 连续三个以上引用折叠为 [a–c]
def cite_repl(m):
    return f'[{num[m.group(1)]}]'
src = re.sub(r'\[\[CITE:([a-zA-Z0-9_]+)\]\]', cite_repl, src)
# 折叠相邻引用组：解析成数字集合，再输出紧凑形式
def _nums(inner):
    out = []
    for part in inner.split(','):
        part = part.strip()
        if '–' in part or '-' in part:
            a, b = re.split(r'[–-]', part)
            out += list(range(int(a), int(b) + 1))
        elif part.isdigit():
            out.append(int(part))
    return out

def _compact(ns):
    ns = sorted(set(ns)); groups = []; cur = [ns[0]]
    for x in ns[1:]:
        if x == cur[-1] + 1: cur.append(x)
        else: groups.append(cur); cur = [x]
    groups.append(cur)
    parts = []
    for g in groups:
        if len(g) >= 3: parts.append(f'{g[0]}–{g[-1]}')
        else: parts += [str(x) for x in g]
    return '[' + ','.join(parts) + ']'

CITE_RUN = re.compile(r'(\[[0-9,\u2013\-]+\])((?:\s*\[[0-9,\u2013\-]+\])+)')
def _merge(m):
    all_nums = _nums(m.group(1)[1:-1])
    for g in re.findall(r'\[([0-9,\u2013\-]+)\]', m.group(0)):
        all_nums += _nums(g)
    return _compact(all_nums)
for _ in range(4):
    src = CITE_RUN.sub(_merge, src)
src = re.sub(r'\[([0-9,\u2013\-]+)\]', lambda m: _compact(_nums(m.group(1))), src)

# ---------- 2. 数字注入 ----------
def num_repl(m):
    body = m.group(1)
    if ':' in body:
        key, spec = body.split(':', 1)
    else:
        key, spec = body, ''
    if key not in N:
        sys.exit(f'FATAL: missing number key {key}')
    v = N[key]
    if spec:
        return format(v, spec)
    if isinstance(v, bool):
        return str(v)
    if isinstance(v, int):
        return f'{v:,}'
    if isinstance(v, float):
        s = f'{v:.10g}'
        return s
    return str(v)

src = re.sub(r'\{\{([^}]+)\}\}', num_repl, src)
if '{{' in src or '}}' in src:
    sys.exit('FATAL: unresolved placeholder remains')

# ---------- 3. 块解析 ----------
lines = src.split('\n')
blocks = []
buf = []
for ln in lines:
    if ln.startswith('# '):
        if buf: blocks.append(('p', '\n'.join(buf))); buf = []
        blocks.append(('title', ln[2:].strip()))
    elif ln.startswith('## '):
        if buf: blocks.append(('p', '\n'.join(buf))); buf = []
        blocks.append(('h1', ln[3:].strip()))
    elif ln.startswith('### '):
        if buf: blocks.append(('p', '\n'.join(buf))); buf = []
        blocks.append(('h2', ln[4:].strip()))
    elif ln.strip() == '':
        if buf: blocks.append(('p', '\n'.join(buf))); buf = []
    else:
        buf.append(ln)
if buf: blocks.append(('p', '\n'.join(buf)))

# ---------- 4. docx ----------
doc = Document()
sec = doc.sections[0]
sec.page_width, sec.page_height = Cm(21.59), Cm(27.94)
for a in ('left_margin', 'right_margin', 'top_margin', 'bottom_margin'):
    setattr(sec, a, Cm(2.54))

st = doc.styles['Normal']
st.font.name = 'Times New Roman'
st.font.size = Pt(12)
st.element.rPr.rFonts.set(qn('w:eastAsia'), 'Times New Roman')
pf = st.paragraph_format
pf.line_spacing_rule = WD_LINE_SPACING.DOUBLE
pf.space_before = Pt(0)
pf.space_after = Pt(0)

# 连续行号
sectPr = sec._sectPr
ln = OxmlElement('w:lnNumType'); ln.set(qn('w:countBy'), '1'); ln.set(qn('w:restart'), 'continuous')
sectPr.append(ln)

# 页脚页码
footer = sec.footer.paragraphs[0]
footer.alignment = WD_ALIGN_PARAGRAPH.CENTER
run = footer.add_run()
fld = OxmlElement('w:fldSimple'); fld.set(qn('w:instr'), 'PAGE')
run._r.addnext(fld)

def add_par(text, bold=False, size=12, align=None, space_after=0, line=None, italic=False):
    p = doc.add_paragraph()
    if align is not None: p.alignment = align
    p.paragraph_format.space_after = Pt(space_after)
    if line:
        p.paragraph_format.line_spacing_rule = line
    r = p.add_run(text)
    r.bold = bold; r.italic = italic; r.font.size = Pt(size)
    r.font.name = 'Times New Roman'
    return p

def add_rich(text, size=12, bold_prefix=None):
    """处理 **bold** 片段与上标标记 ^x^"""
    p = doc.add_paragraph()
    p.paragraph_format.space_after = Pt(0)
    tokens = re.split(r'(\*\*[^*]+\*\*|\^[^^]+\^)', text)
    for tk in tokens:
        if not tk: continue
        if tk.startswith('**') and tk.endswith('**'):
            r = p.add_run(tk[2:-2]); r.bold = True
        elif tk.startswith('^') and tk.endswith('^') and len(tk) > 2:
            r = p.add_run(tk[1:-1]); r.font.superscript = True
        else:
            r = p.add_run(tk)
        r.font.size = Pt(size); r.font.name = 'Times New Roman'
    return p

first_h1 = True
for kind, text in blocks:
    if kind == 'title':
        add_par(text, bold=True, size=16, space_after=12, line=WD_LINE_SPACING.ONE_POINT_FIVE)
    elif kind == 'h1':
        if text in ('Abstract', 'Author summary'):
            add_par(text, bold=True, size=12, space_after=6)
        elif text == 'References':
            add_par(text, bold=True, size=12, space_after=6)
            for k in order:
                p = doc.add_paragraph()
                p.paragraph_format.space_after = Pt(0)
                p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE
                r = p.add_run(f'{num[k]}. {REFS[k]}'); r.font.size = Pt(12); r.font.name = 'Times New Roman'
            REF_DONE = True
        else:
            add_par(text, bold=True, size=12, space_after=6)
    elif kind == 'h2':
        add_par(text, bold=True, size=12, space_after=6)
    else:
        t = text.replace('\n', ' ')
        if t.strip() == '[[REFS]]':
            continue
        # 图注：以 **Fig 开头
        if t.startswith('**Fig'):
            add_rich(t, size=11)
        elif t.startswith('**S1 Appendix'):
            add_rich(t, size=12)
        elif 'These authors contributed equally' in t or 'Corresponding author' in t:
            add_rich(t.replace('\\*', '*'), size=10)
        elif t.startswith('Department of Gynecology'):
            add_rich(t, size=10)
        elif t.startswith('Correspondence:') or t.startswith('Short title:'):
            add_par(t, size=10, space_after=2)
        elif t.startswith('Chaohua Tang'):
            add_rich(t.replace('\\*', '*'), size=12)
        else:
            add_rich(t.replace('\\*', '*'), size=12)

out = os.path.join(BASE, '01_Manuscript_v8.docx')
doc.save(out)
print('saved', out)
print('citations:', len(order), 'words(md):', len(re.findall(r"[A-Za-z0-9µαβγ'’\-]+", src)))

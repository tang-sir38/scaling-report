r"""98_build_v4.py — 由模板 + 数字锁 v7 生成 v4.0 手稿（NCS Analysis 规格）

设计
----
文本与数字分离：
  07_SUBMISSION/_v4/MANUSCRIPT_v4.md   散文 + {{锁路径:格式}} + [[refkey]]
  02_OUT/NUMBERS_LOCK_v7.json          唯一权威数字
  02_OUT/REFS_v4.json                  键控文献表

本脚本做四件事：
  1. 由锁现场构造扁平 digest（键名面向散文，路径留在 ledger 里）
  2. 解析占位符并按 Python 格式串注入；缺键即报错，不静默
  3. 文献按**正文首现顺序**编号，输出参考文献表
  4. 渲染 docx（A4 / 12pt / 双倍行距 / 行号），并写出 ledger 供守卫复核

产物
----
  07_SUBMISSION/_build_v4/docs/01_Manuscript.docx
  07_SUBMISSION/_build_v4/numbers_ledger_v4.json
  07_SUBMISSION/_build_v4/stats_v4.json     （词数 / 图位数 / 文献数 / 自检）

用法
----
  python 98_build_v4.py
"""
from __future__ import annotations

import json
import re
from pathlib import Path

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt

ROOT = Path(__file__).resolve().parents[1]
LOCK = ROOT / "02_OUT/NUMBERS_LOCK_v7.json"
REFS = ROOT / "02_OUT/REFS_v4.json"
TPL = ROOT / "07_SUBMISSION/_v4/MANUSCRIPT_v4.md"
BUILD = ROOT / "07_SUBMISSION/_build_v4"
DOCS = BUILD / "docs"
FIGS = BUILD / "figs"

MAIN_TEXT_LIMIT = 3500
MAX_ITEMS = 6
MAX_REFS = 50


# --------------------------------------------------------------- 1. digest
def get(d, path):
    cur = d
    for part in path.split("."):
        if isinstance(cur, dict) and part in cur:
            cur = cur[part]
        elif isinstance(cur, list) and part.isdigit():
            cur = cur[int(part)]
        else:
            raise KeyError(f"锁里没有路径 {path}（断在 {part}）")
    return cur


TR: dict[str, str] = {}      # 键 -> 施于锁值的变换名（守卫据此重放，禁靠人工告知答案）


def build_digest(L):
    D, P = {}, {}

    def put(key, path, tf=None, scale=1.0):
        v = get(L, path)
        if tf is not None:
            v = tf(v)
            TR[key] = getattr(tf, "__name__", "tf")
        D[key] = float(v) * scale
        P[key] = path + (f"  ×{scale}" if scale != 1.0 else "")

    for t, k in (("colon", "c"), ("lung", "l"), ("brain", "b")):
        a2 = get(L, f"datasets.{t}.alpha.2")
        a256 = get(L, f"datasets.{t}.alpha.256")
        put(f"{k}.a2", f"datasets.{t}.alpha.2")
        put(f"{k}.a256", f"datasets.{t}.alpha.256")
        D[f"{k}.fall"] = a2 / a256
        P[f"{k}.fall"] = f"datasets.{t}.alpha.2 / alpha.256"
        put(f"{k}.gamma", f"datasets.{t}.gamma", abs)
        put(f"{k}.r2", f"datasets.{t}.r2")
        put(f"{k}.nbins", f"datasets.{t}.n_bins")
        put(f"{k}.umi", f"datasets.{t}.umi_per_bin")
        put(f"{k}.pois", f"derived.{t}_obs_over_poisson_256")

    # 第二张结肠切片
    p1 = "multi.sections.Visium_HD_Human_Colon_Cancer_P1."
    put("c2.gamma", p1 + "gamma_pool")
    put("c2.sd", p1 + "boot_sd")
    put("c2.med", p1 + "gamma_median")
    D["c2.meddiff"] = get(L, p1 + "gamma_median") - get(L, "gene_level.colon.gamma_gene.median")
    P["c2.meddiff"] = p1 + "gamma_median − gene_level.colon.gamma_gene.median"

    # 下流分析代价（106 -> 97.downstream）
    ds = "downstream.summary."
    for k, src in (("ds.phi_lo", "phi_lo"), ("ds.phi_hi", "phi_hi"),
                   ("ds.fpr_lo", "fpr_lo"), ("ds.fpr_hi", "fpr_hi"),
                   ("ds.neff", "neff_hi"), ("ds.jacc32", "jacc_32"),
                   ("ds.top16", "top100_16"), ("ds.top64", "top100_64"),
                   ("ds.sigref", "sig_ref"), ("ds.sig128", "sig_128"),
                   ("ds.sl64", "slope_log_64"), ("ds.sr64", "slope_raw_64"),
                   ("ds.r64", "r_64"), ("ds.rawmin", "map_raw_min"),
                   ("ds.rawmax", "map_raw_max"), ("ds.rawmaxx", "map_raw_max_x"),
                   ("ds.rawrmin", "map_raw_r_min"),
                   ("ds.logmin", "map_log_min"), ("ds.logmax", "map_log_max")):
        put(k, ds + src)

    # Stereo-seq
    put("st.gamma", "cross.stereo.gamma_pool")
    put("st.med", "cross.stereo.gamma_median")
    put("st.delta", "cross.delta_gamma_pool")
    put("st.n", "cross.replication.n")
    put("st.rho", "cross.replication.rho_raw")
    put("st.rhopart", "cross.replication.rho_partial_abundance")

    # 逐基因
    for t, k in (("colon", "c"), ("lung", "l"), ("brain", "b")):
        g = f"gene_level.{t}.gamma_gene."
        c = f"gene_level.{t}.composition."
        put(f"g.{k}.pooled", g + "pooled")
        put(f"g.{k}.med", g + "median")
        put(f"g.{k}.neff", c + "n_eff")
        put(f"g.{k}.top100", c + "top100_share")
        put(f"g.{k}.share1", f"derived4.top1_alpha_share.{t}")
        put(f"g.{k}.gt1", g + "frac_gt_1")
        put(f"g.{k}.gt15", g + "frac_gt_1p5")
    put("g.c.top1share", "decomp.by_tissue.colon.top1_mu2_share")

    # 理论
    D["th.idmaxpct"] = get(L, "theory.identity.max_abs_rel_err") * 100
    P["th.idmaxpct"] = "theory.identity.max_abs_rel_err ×100"
    D["th.idmedpct"] = get(L, "theory.identity.median_abs_rel_err") * 100
    P["th.idmedpct"] = "theory.identity.median_abs_rel_err ×100"
    put("th.wn", "theory.families.white_0.gamma_2_64", abs)
    put("th.r2", "theory.gamma_vs_beta_emp.r2")
    put("th.depth", "theory.depth_gamma_range")
    put("th.ex2a", "theory.families.exponential_2.gamma_2_64", abs)
    put("th.ex2b", "theory.families.exponential_2.gamma_16_64", abs)
    put("th.ex8a", "theory.families.exponential_8.gamma_2_64", abs)
    put("th.ex8b", "theory.families.exponential_8.gamma_16_64", abs)

    # gamma vs beta
    tot = 0
    for t, k in (("colon", "c"), ("lung", "l"), ("brain", "b")):
        gb = f"gamma_beta.by_tissue.{t}."
        put(f"gb.{k}.r", gb + "pearson_raw")
        put(f"gb.{k}.rho", gb + "spearman")
        put(f"gb.{k}.diff", gb + "diff_median")
        put(f"gb.{k}.g", gb + "gamma_median")
        put(f"gb.{k}.b", gb + "beta_median")
        tot += int(get(L, gb + "n"))
    D["gb.total"] = tot
    P["gb.total"] = "gamma_beta.by_tissue.*.n 之和"

    # rho 形状
    for t, k in (("colon", "c"), ("brain", "b")):
        put(f"rho.{k}.trend", f"rho_shape.by_tissue.{t}.slope_vs_lag_spearman")
    sl_c = [r["local_slope"] for r in get(L, "rho_shape.by_tissue.colon.by_lag")
            if r["local_slope"] is not None]
    D["rho.c.lo"], D["rho.c.hi"] = float(min(sl_c)), float(max(sl_c))
    P["rho.c.lo"] = "rho_shape.by_tissue.colon.by_lag.local_slope min"
    P["rho.c.hi"] = "rho_shape.by_tissue.colon.by_lag.local_slope max"
    put("rho.b.lo", "rho_shape.by_tissue.brain.slope_by_range.8-16um")
    put("rho.b.hi", "rho_shape.by_tissue.brain.slope_by_range.136-248um")

    # 干预
    put("iv.igkc.base", "intervene.single_gene.genes.IGKC.baseline_gamma")
    put("iv.igkc.perm", "intervene.single_gene.genes.IGKC.permute_gamma")
    put("iv.igkc.lo", "intervene.single_gene.genes.IGKC.permute_lo")
    put("iv.igkc.hi", "intervene.single_gene.genes.IGKC.permute_hi")
    put("iv.pool.base", "intervene.single_gene.pooled.gamma_baseline")
    put("iv.pool.perm", "intervene.single_gene.pooled.cases.permute.gamma")
    put("iv.pool.excise", "intervene.single_gene.pooled.cases.excise.gamma")
    put("iv.pool.isolate", "intervene.single_gene.pooled.cases.isolate.gamma")
    put("iv.pool.delta", "intervene.single_gene.pooled.delta_permute")
    for i, n in enumerate([0, 10, 30, 100, 300, 1000]):
        put(f"iv.dn{i}", f"intervene.dose_response.curves.{n}.gamma")
    put("iv.sh1", "intervene.dose_response.curves.10.mu2_share")
    put("iv.dmin", "intervene.single_gene.summary.delta_min")
    put("iv.dmax", "intervene.single_gene.summary.delta_max")

    # 稳健性
    put("rob.pooled.sd", "robust.summary.pooled.sd")
    put("rob.pooled.range", "robust.summary.pooled.range")
    put("rob.median.sd", "robust.summary.median.sd")
    put("rob.median.range", "robust.summary.median.range")

    # 设计
    for t, k in (("colon", "c"), ("lung", "l"), ("brain", "b")):
        put(f"des.mde.{k}", f"datasets.{t}.mde_gamma.n1")
    put("des.depth", "derived3.depth_delta_max")
    put("des.drop", "jackknife.colon.loo.max_abs_delta")
    put("des.gmax", "derived.lung_gamma_max")
    D["des.gap_cl"] = abs(D["c.gamma"] - D["l.gamma"])
    P["des.gap_cl"] = "|c.gamma − l.gamma|"
    D["des.gap_cb"] = abs(D["c.gamma"] - D["b.gamma"])
    P["des.gap_cb"] = "|c.gamma − b.gamma|"

    return D, P


# --------------------------------------------------------------- 2. 注入
PH = re.compile(r"\{\{([^}:]+?)(?::([^}]+))?\}\}")


def resolve(text, D, P, ledger):
    def sub(m):
        key, fmt = m.group(1).strip(), (m.group(2) or ".3f").strip()
        if key not in D:
            raise KeyError(f"模板要 {key}，digest 里没有")
        v = D[key]
        try:
            s = format(v, fmt) if fmt != "," else format(int(round(v)), ",")
        except ValueError as e:
            raise ValueError(f"格式串 {fmt!r} 不能用于 {key}={v}: {e}")
        ledger[key] = dict(value=v, formatted=s, fmt=fmt, lock_path=P[key],
                           tf=TR.get(key))
        return s
    return PH.sub(sub, text)


CITE = re.compile(r"\[\[([^\]]+)\]\]")


def number_citations(text, REFS, order):
    def sub(m):
        keys = [k.strip() for k in m.group(1).split(",")]
        nums = []
        for k in keys:
            if k not in REFS:
                raise KeyError(f"未定义的文献键 {k}")
            if k not in order:
                order.append(k)
            nums.append(order.index(k) + 1)
        nums.sort()
        # 连续 3 个以上压成区间
        out, i = [], 0
        while i < len(nums):
            j = i
            while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
                j += 1
            if j - i >= 2:
                out.append(f"{nums[i]}–{nums[j]}")
            else:
                out.extend(str(x) for x in nums[i:j + 1])
            i = j + 1
        return "[" + ",".join(out) + "]"
    return CITE.sub(sub, text)


# --------------------------------------------------------------- 3. docx
def set_base(doc, size=12, spacing=2.0):
    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(size)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), "SimSun")
    pf = st.paragraph_format
    pf.line_spacing = spacing
    pf.space_after = Pt(0)
    for s in doc.sections:
        s.top_margin = s.bottom_margin = Cm(2.54)
        s.left_margin = s.right_margin = Cm(2.54)


def add_line_numbers(doc):
    for s in doc.sections:
        sectPr = s._sectPr
        ln = OxmlElement("w:lnNumType")
        ln.set(qn("w:countBy"), "1")
        ln.set(qn("w:restart"), "continuous")
        sectPr.append(ln)


def add_page_numbers(doc):
    for s in doc.sections:
        p = s.footer.paragraphs[0]
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        r = p.add_run()
        for tag, txt in (("begin", None), (None, "PAGE"), ("end", None)):
            if tag:
                fc = OxmlElement("w:fldChar"); fc.set(qn("w:fldCharType"), tag)
                r._r.append(fc)
            else:
                it = OxmlElement("w:instrText"); it.set(qn("xml:space"), "preserve")
                it.text = txt
                r._r.append(it)


SUP = re.compile(r"\^(\([^)]*\)|[−\-+]?[0-9]+(?:\.[0-9]+)?|[−\-+]?[α-ωΑ-Ω])")
SUB = re.compile(r"_([A-Za-z0-9]+)")


def add_runs(p, text, bold=False, italic=False, size=None):
    """行内标记：**bold**；^2 / ^−0.554 走上标；N_A / μ_gd 走下标。

    投稿件里不该出现 d^-2、N_A 这种字面写法，Word 能排的就必须排。
    """
    for chunk in re.split(r"(\*\*[^*]+\*\*)", text):
        if not chunk:
            continue
        b = chunk.startswith("**") and chunk.endswith("**")
        body = chunk[2:-2] if b else chunk
        for j, piece in enumerate(SUP.split(body)):
            if piece == "":
                continue
            if j % 2 == 1:                       # 指数
                r = p.add_run(piece[1:-1] if piece.startswith("(") else piece)
                r.font.superscript = True
                r.bold = bold or b
                r.italic = italic
                if size:
                    r.font.size = Pt(size)
                continue
            for k, sub in enumerate(SUB.split(piece)):
                if sub == "":
                    continue
                r = p.add_run(sub)
                if k % 2 == 1:                   # 下标
                    r.font.subscript = True
                r.bold = bold or b
                r.italic = italic
                if size:
                    r.font.size = Pt(size)


def para(doc, text, bold=False, italic=False, size=None, align="justify",
         space_after=0, spacing=2.0, indent=True):
    p = doc.add_paragraph()
    p.alignment = {"justify": WD_ALIGN_PARAGRAPH.JUSTIFY,
                   "left": WD_ALIGN_PARAGRAPH.LEFT,
                   "center": WD_ALIGN_PARAGRAPH.CENTER}[align]
    p.paragraph_format.line_spacing = spacing
    p.paragraph_format.space_after = Pt(space_after)
    if indent and align == "justify":
        p.paragraph_format.first_line_indent = Cm(0.75)
    add_runs(p, text, bold=bold, italic=italic, size=size)
    return p


def heading(doc, text, size=13):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(12)
    p.paragraph_format.space_after = Pt(6)
    p.paragraph_format.line_spacing = 1.5
    p.paragraph_format.keep_with_next = True   # 标题不许孤零零留在页脚
    r = p.add_run(text)
    r.bold = True
    r.font.size = Pt(size)
    return p


def build_docx(blocks, refs_list, out):
    doc = Document()
    set_base(doc)
    add_line_numbers(doc)
    add_page_numbers(doc)

    def emit(sec_name, lines):
        if sec_name == "TITLE":
            for t in lines:
                para(doc, t, bold=True, size=15, align="center", spacing=1.5)
            return
        if sec_name == "AUTHORS":
            for t in lines:
                para(doc, t, align="center", spacing=1.5)
            return
        if sec_name in ("RUNNING_TITLE", "AFFILIATION", "CORRESPONDENCE", "KEYWORDS"):
            heading(doc, sec_name.replace("_", " ").title(), size=11)
            for t in lines:
                para(doc, t, size=11, spacing=1.5)
            return
        if sec_name == "FIGURE_LEGENDS":
            heading(doc, "Figure legends", size=14)
            for t in lines:
                m = re.match(r"\*\*Fig\.\s*(\d+)\.", t)
                if m:                      # 图注前插入该图，图与注不分离
                    num = int(m.group(1))
                    png = FIGS / f"Fig{num}.png"
                    if png.exists():
                        pic = doc.add_paragraph()
                        pic.alignment = WD_ALIGN_PARAGRAPH.CENTER
                        pic.paragraph_format.space_before = Pt(14)
                        pic.paragraph_format.space_after = Pt(4)
                        pic.paragraph_format.keep_with_next = True
                        pic.add_run().add_picture(str(png), width=Cm(16.5))
                para(doc, t, size=11, spacing=1.5)
            return
        if sec_name == "REFERENCES":
            heading(doc, "References")
            for i, t in enumerate(refs_list, 1):
                p = doc.add_paragraph()
                p.paragraph_format.line_spacing = 1.5
                p.paragraph_format.space_after = Pt(3)
                p.paragraph_format.left_indent = Cm(0.9)
                p.paragraph_format.first_line_indent = Cm(-0.9)
                p.add_run(f"{i}. {t}")   # 文献串不做上下标解析，避免误伤 DOI/基因名
            return
        heading(doc, sec_name.replace("_", " ").title())
        for t in lines:
            para(doc, t)

    for name, lines in blocks:
        if name in ("BODY", "METHODS"):
            heading(doc, name.title(), size=14)
            for t in lines:
                if t.startswith("### "):
                    heading(doc, t[4:], size=12)
                elif t.startswith("## "):
                    heading(doc, t[3:], size=12)
                elif t.startswith("  "):
                    p = doc.add_paragraph()
                    p.paragraph_format.line_spacing = 1.5
                    p.paragraph_format.left_indent = Cm(1.2)
                    add_runs(p, t.strip(), italic=True)
                else:
                    para(doc, t)
        else:
            emit(name, lines)

    out.parent.mkdir(parents=True, exist_ok=True)
    props = doc.core_properties
    props.author = "Chaohua Tang"
    props.title = "Binning is set by spatial covariance and by which transcripts dominate"
    doc.save(str(out))


# --------------------------------------------------------------- 4. main
def main():
    docx_ = json.loads(LOCK.read_text(encoding="utf-8"))
    refs = json.loads(REFS.read_text(encoding="utf-8"))["refs"]
    raw = TPL.read_text(encoding="utf-8")

    D, P = build_digest(docx_)
    ledger: dict[str, dict] = {}

    # 切节
    blocks, cur, buf = [], None, []
    for line in raw.splitlines():
        if line.startswith("# "):
            if cur:
                blocks.append((cur, [x for x in buf if x.strip()]))
            cur, buf = line[2:].strip(), []
        else:
            buf.append(line)
    if cur:
        blocks.append((cur, [x for x in buf if x.strip()]))

    names = [b[0] for b in blocks]
    for req in ("TITLE", "ABSTRACT", "BODY", "METHODS", "FIGURE_LEGENDS", "DECLARATIONS", "REFERENCES"):
        if req not in names:
            raise SystemExit(f"模板缺节 {req}（现有 {names}）")

    order: list[str] = []
    out_blocks, main_words, abs_words = [], 0, 0
    for name, lines in blocks:
        new = []
        for t in lines:
            t = resolve(t, D, P, ledger)
            t = number_citations(t, refs, order)
            new.append(t)
        if name == "BODY":
            txt = " ".join(new)
            txt = re.sub(r"^#{2,3} ", "", txt, flags=re.M)
            main_words = len([w for w in txt.split() if w.strip()])
        if name == "ABSTRACT":
            abs_words = len(" ".join(new).split())
        out_blocks.append((name, new))

    # 参考文献表按首现顺序
    refs_list = [refs[k] for k in order]
    unused = [k for k in refs if k not in order]

    DOCS.mkdir(parents=True, exist_ok=True)
    build_docx(out_blocks, refs_list, DOCS / "01_Manuscript.docx")

    n_items = len([1 for nm, ls in out_blocks if nm == "FIGURE_LEGENDS"
                   for t in ls if t.startswith("**Fig.")])
    stats = dict(main_text_words=main_words, abstract_words=abs_words,
                 n_display_items=n_items, n_refs=len(refs_list),
                 n_cited_keys=len(order), unused_ref_keys=unused,
                 n_placeholders_used=len(ledger),
                 limits=dict(main_text=MAIN_TEXT_LIMIT, items=MAX_ITEMS, refs=MAX_REFS))
    (BUILD / "stats_v4.json").write_text(json.dumps(stats, indent=2, ensure_ascii=False),
                                         encoding="utf-8")
    (BUILD / "numbers_ledger_v4.json").write_text(
        json.dumps(ledger, indent=1, ensure_ascii=False), encoding="utf-8")

    print(f"写出 {DOCS/'01_Manuscript.docx'}")
    print(f"  正文词数 {main_words}（上限 {MAIN_TEXT_LIMIT}）")
    print(f"  摘要词数 {abs_words}（目标 100-150）")
    print(f"  图位 {n_items}（上限 {MAX_ITEMS}）")
    print(f"  文献 {len(refs_list)}（上限 {MAX_REFS}）；未引用 {len(unused)} 个")
    print(f"  注入占位符 {len(ledger)} 个，全部有锁路径")
    bad = [k for k, v in stats.items() if False]
    if main_words > MAIN_TEXT_LIMIT:
        print(f"  !! 正文超限 {main_words - MAIN_TEXT_LIMIT} 词")
    if n_items > MAX_ITEMS:
        print(f"  !! 图位超限")
    if len(refs_list) > MAX_REFS:
        print(f"  !! 文献超限")


if __name__ == "__main__":
    main()

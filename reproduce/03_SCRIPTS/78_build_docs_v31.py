# -*- coding: utf-8 -*-
r"""
78_build_docs_v31.py — A 题 v3.1 正文 / 投稿信 / 补充材料 / 投稿指南，数字全部由 NUMBERS_LOCK_v5.json 注入。

产物（写到 --outdir）：
  01_Manuscript.docx       Title / Abstract / Background / Results / Discussion / Methods /
                           Figure legends / Declarations / References / Table 1-3
  02_Cover_Letter.docx     1 页
  04_Supplementary.docx    Supplementary methods, Table S1-S5, Fig S1-S5
  00_投稿指南.md            中文操作手册
numbers_ledger_v2.json     每个格式化数字的出处（供守卫复核）

用法：python 45_build_docs_v2.py --outdir ...
"""
import argparse
import json
import math
import re
from pathlib import Path

from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "07_SUBMISSION/_build_v3.1"

# ----------------------------------------------------------------- 数字注入
class NS(dict):
    """dict 包装：允许数字下标（2 -> "2"）、区间下标（2:4 -> "2-4"），
    以及非数字键字典的位置下标（用于按层顺序取表）。"""
    __getattr__ = dict.__getitem__

    def __missing__(self, k):
        if isinstance(k, int):
            if dict.__contains__(self, str(k)):
                return dict.__getitem__(self, str(k))
            vals = list(dict.values(self))
            if 0 <= k < len(vals):
                return vals[k]
        if isinstance(k, slice):
            alt = "%s-%s" % (k.start, k.stop)
            if dict.__contains__(self, alt):
                return dict.__getitem__(self, alt)
        raise KeyError(k)


def wrap(o):
    if isinstance(o, dict):
        return NS({k: wrap(v) for k, v in o.items()})
    if isinstance(o, list):
        return [wrap(v) for v in o]
    return o


RAW = json.loads((ROOT / "02_OUT/NUMBERS_LOCK_v5.json").read_text(encoding="utf-8"))
N = wrap(RAW)
X = RAW["cross"]      # plain dict: the tables in the supplement index it directly
MU = RAW["multi"]
LEDGER = []
LOADED = False


class _Lazy(dict):
    """N 在文献加载后才可用；先占位，避免 import 顺序问题。"""
    def __getattr__(self, k):
        return N[k]


def PP(x):
    """depth.runs 的键是 "p0.46" 这种形式，用 PP(0.46) 取。"""
    return "p%g" % x


def num(expr, nd=2, kind="f"):
    v = eval(expr, {"__builtins__": {}},
             {"N": N, "d": N.datasets, "D3": N.derived3, "J": N.jackknife,
              "F": N.fill.by_tissue, "B": N.boot.by_tissue, "LW": N.layer_weights.by_tissue,
              "L": N.layers, "S": N.synthetic, "V": N.variance, "M": N.models,
              "RL": N.reliability_layers, "n1": "n1", "n3": "n3", "PP": PP,
              "X": N.cross, "MU": N.multi,
              "abs": abs, "math": math, "np": math, "min": min, "max": max,
              "sum": sum, "len": len, "round": round,
              "log": math.log, "exp": math.exp, "sqrt": math.sqrt})
    if kind == "f":
        s = f"{v:.{nd}f}"
    elif kind == "pct":
        s = f"{v * 100:.{nd}f}%"
    elif kind == "int":
        s = f"{int(round(v)):,}"
    elif kind == "g":
        s = f"{v:g}"
    else:
        raise ValueError(kind)
    LEDGER.append({"expr": expr, "nd": nd, "kind": kind, "text": s})
    return s


# ----------------------------------------------------------------- 参考文献
def load_refs():
    return json.loads((BUILD / "refs_lock_v2.json").read_text(encoding="utf-8"))["ok"]


REFS = load_refs()
CITE_ORDER = []


def cite(m):
    keys = [k.strip() for k in m.group(1).split(",")]
    nums = []
    for k in keys:
        if k not in REFS:
            raise KeyError(f"unknown ref {k}")
        if k not in CITE_ORDER:
            CITE_ORDER.append(k)
        nums.append(CITE_ORDER.index(k) + 1)
    nums = sorted(set(nums))
    out, i = [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        out.append(f"{nums[i]}\u2013{nums[j]}" if j - i >= 2 else
                   ("%d,%d" % (nums[i], nums[j]) if j > i else str(nums[i])))
        i = j + 1
    return "[" + ",".join(out) + "]"


CITE_RE = re.compile(r"\[\[([a-zA-Z0-9_,\s]+)\]\]")


def R(s):
    return CITE_RE.sub(cite, s)


# ----------------------------------------------------------------- 文本
TITLE = ("The pooled noise\u2013resolution exponent of binned spatial transcriptomics "
         "is a composition statistic")
RUNNING = "The pooled exponent is a composition statistic"

_R = ("N.gene_level.%s.gamma_gene.p95/N.gene_level.%s.gamma_gene.p05")
RATIO_MIN = "min(" + ",".join(_R % (t, t) for t in ("colon", "lung", "brain")) + ")"
RATIO_MAX = "max(" + ",".join(_R % (t, t) for t in ("colon", "lung", "brain")) + ")"
_Q = "N.gene_level.%s.gamma_gene.q75/N.gene_level.%s.gamma_gene.q25"
QR_MIN = "min(" + ",".join(_Q % (t, t) for t in ("colon", "lung", "brain")) + ")"
QR_MAX = "max(" + ",".join(_Q % (t, t) for t in ("colon", "lung", "brain")) + ")"

AUTHORS = ("Chaohua Tang\u00b9,*,**, Fangfang Wu\u00b9,*, Yingying Zheng\u00b9,*, Lingling Zhan\u00b9, "
           "Hong Tan\u00b9, Lijing Lin\u00b9, Meimei Wen\u00b9")
AFFIL = ("\u00b9 Department of Gynecology, The 95th Hospital of Putian, No. 485 Dongyan Road, "
         "Licheng District, Putian City, Fujian Province 351100, China")
CORR = ("Chaohua Tang, Department of Gynecology, The 95th Hospital of Putian, "
        "No. 485 Dongyan Road, Licheng District, Putian City, Fujian Province 351100, China. "
        "Tel: +86 13860904395. Email: tong@fjmu.edu.cn. ORCID: 0009-0004-7080-8993.")
KEYWORDS = ("spatial transcriptomics; Visium HD; overdispersion; scaling; "
            "transcript abundance; experimental design")

ABSTRACT = {
    "Background": (
        "Binned spatial transcriptomics buys counts with resolution, on the expectation that "
        "overdispersion falls as bin area grows: \u03b1 \u221d d^\u22122. That expectation is counting "
        "statistics, yet the fitted exponent is reported as a property of the tissue. It has been measured "
        "only at one or two resolutions."),
    "Results": (
        f"We re-aggregated native 2 \u00b5m bins of three Visium HD sections \u2014 human colon cancer, "
        f"human lung cancer and mouse brain \u2014 to 256 \u00b5m in eight exact steps. Across a "
        f"{num('(256/2)**2', 0, 'int')}-fold increase in area, \u03b1 fell "
        f"{num('N.derived.colon_alpha_ratio_2_to_256', 1)}-, {num('N.derived.lung_alpha_ratio_2_to_256', 1)}- and "
        f"{num('N.derived.brain_alpha_ratio_2_to_256', 1)}-fold, not {num('(256/2)**2', 0, 'int')}-fold. No single "
        f"correlation length survives: a finite-\u03be model trails the best fit by "
        f"{num('min(D3.colon_xi_dAIC_min,D3.lung_xi_dAIC_min,D3.brain_xi_dAIC_min)', 0, 'int')}\u2013"
        f"{num('max(D3.colon_xi_dAIC_max,D3.lung_xi_dAIC_max,D3.brain_xi_dAIC_max)', 0, 'int')} AIC units. "
        f"The exponent belongs to transcripts, not tissues: per-gene "
        f"\u03b3 varies {num(RATIO_MIN, 1)}\u2013{num(RATIO_MAX, 1)}-fold between the 5th and 95th percentile and "
        f"falls monotonically with abundance in all three sections "
        f"(\u03c1 = \u2212{num('abs(J.colon.gamma_vs_mu.spearman)', 2)}, "
        f"\u2212{num('abs(J.lung.gamma_vs_mu.spearman)', 2)}, "
        f"\u2212{num('abs(J.brain.gamma_vs_mu.spearman)', 2)}), reproduces between halves of one section "
        f"(r = {num('N.splithalf.reliability.pearson', 2)}). A second section of "
        f"the same tumour reproduces the pooled value (difference "
        f"{num('abs(MU.contrasts.same_tissue_colon.delta)', 3)}) but not the per-gene median "
        f"({num('MU.pergene.shift_overall', 3)}); a brain section on a different chemistry reproduces the "
        f"per-gene structure (\u03c1 = {num('X.replication.rho_raw', 2)}) at a different absolute level. "
        f"It is a \u03bc\u00b2-weighted mean carried by "
        f"{num('min(N.gene_level.colon.composition.n_eff,N.gene_level.lung.composition.n_eff,N.gene_level.brain.composition.n_eff)', 1)}\u2013"
        f"{num('max(N.gene_level.colon.composition.n_eff,N.gene_level.lung.composition.n_eff,N.gene_level.brain.composition.n_eff)', 1)} "
        f"effective genes "
        f"({num('N.decomp.min_neff_alpha', 1)}\u2013{num('N.decomp.max_neff_alpha', 1)} when weighted by each "
        f"transcript's share of the exponent)."),
    "Conclusions": (
        "Report a per-gene distribution, an effective gene count and regional replicates, not \u201cthe "
        "exponent of a tissue\u201d. The pooled value is the only summary that travels between sections and "
        "between platforms; a per-gene summary belongs in a distribution, not in a single number. "
        "Because every transcript has \u03b3 < 2, merging bins cannot deliver the "
        "precision that counting statistics promise, and depth buys almost nothing; the budget belongs in "
        "sections."),
}

BACKGROUND = [
    ("p", "Spatial transcriptomics now resolves expression at the scale of a cell or below. Visium HD tiles a "
          "section with a contiguous 2 \u00b5m grid; Stereo-seq uses DNA nanoball arrays; Slide-seq, seqFISH, "
          "MERFISH, CosMx and Xenium reach comparable or finer scales "
          "[[stahl2016,stickels2021,chen2022stereo,eng2019seqfish,chen2015merfish,he2022cosmx,janesick2023]]. "
          "None of these assays counts a transcript in every bin. A 2 \u00b5m bin holds tens of transcripts, so "
          "each bin is a sparse sample of the expression field beneath it [[moses2022,rao2021]]."),
    ("p", "The field handles that sparsity with one move: merge neighbouring bins. The merged unit is treated as "
          "a more precise measurement of the same quantity, because a larger area should contain more "
          "transcripts. If transcripts fell independently, the variance of an expression estimate would fall as "
          "1/area, and the overdispersion \u03b1 = (Var \u2212 \u03bc)/\u03bc\u00b2 would fall as d^\u22122. The "
          "move is everywhere: Visium HD workflows merge 2 \u00b5m bins into 8 or 16 \u00b5m before clustering "
          "[[sthd2025,bin2cell2024]], Stereo-seq analyses merge the 220 nm grid into 25\u201350 \u00b5m bins "
          "[[chen2022stereo]], benchmarking studies recommend matching bin size to the question [[li2022bench]], "
          "and best-practice reviews treat resolution as a knob to trade against depth "
          "[[moses2022,luecken2019,heumos2023]]. Count models put the variance of a feature in inverse "
          "proportion to its counts [[anders2010,love2014,law2014,grun2014,hafemeister2019,ahlmanneltze2023]]. "
          "The zeros once read as dropout are themselves better explained by an overdispersed count "
          "distribution than by zero inflation [[svensson2020]]."),
    ("p", "Merging is not the only response to that sparsity. Sub-spot and single-cell profiles can be recovered "
          "by deconvolution or by alignment to a reference [[zhao2021bayesspace,klesh2022cell2loc,cable2022,"
          "biancalani2021]], counts can be denoised or smoothed [[eraslan2019,dries2021,palla2022]], and "
          "spatially variable genes can be called without ever summarizing a whole section "
          "[[svensson2018,edsgard2018,sun2020spark,weber2023nnsvg]]. Each of these methods takes the binned or "
          "per-spot count that the assay delivers and borrows strength across units. None of them changes how that "
          "count is generated, so none of them changes the exponent that governs how its overdispersion decays "
          "with area. What they change is the precision of a cell-type or gene-level estimate at a fixed bin size, "
          "which is a different quantity from the one compared between tissues."),
    ("p", "Written as a scaling law, \u03b1(d) \u221d d^\u2212\u03b3, the exponent \u03b3 is a design "
          "parameter. At \u03b3 = 2, coarsening from 2 to 256 \u00b5m cuts noise 16,384-fold and resolution is "
          "nearly free; at \u03b3 \u2248 0.5 the same move costs resolution and buys one order of magnitude. The "
          "same number is also used descriptively, as the noise\u2013resolution exponent of a tissue, and "
          "compared between tissues and conditions. Those two uses are compatible only if \u03b3 is a property "
          "of the tissue \u2014 and that is the assumption we test."),
    ("p", "The ingredients for testing it exist outside genomics. Ecology has measured how variance depends on "
          "grain size since the 1980s, and its variance\u2013grain curves are usually shallower than d^\u22122 "
          "[[wiens1989]]; geography names the same dependence the modifiable areal unit problem [[jelinski1996]]; "
          "spatial autocorrelation is known to shrink the effective number of independent samples "
          "[[legendre1993]]. Inside genomics, the mean\u2013variance dependency of expression is a recognised "
          "confounder for any statistic that compares variability [[eling2019,brennecke2013,grun2014,sarkar2021]], "
          "and single-cell studies learned to test at the donor level rather than the cell level "
          "[[squair2021,zimmerman2021,murphy2022,svensson2017,ess2014]]. Reviews of the field's open problems "
          "begin from counts already in hand and treat their grain as settled [[lahnemann2020]]. What has not "
          "been done is to vary grain "
          "continuously across two orders of magnitude inside one section, and to decompose the resulting "
          "exponent into the parts that belong to the tissue and the parts that belong to the transcripts."),
    ("p", "Here we do that in three Visium HD sections \u2014 human colon cancer, human lung cancer and mouse "
          "brain \u2014 by exact integer re-aggregation of native 2 \u00b5m bins to 256 \u00b5m in eight steps. "
          "We ask four questions in order. How fast does noise actually fall? Is there a single correlation "
          "length behind the curve? Whose exponent is it? And what follows for design? The answers change the "
          "interpretation of a number the field already reports, and leave its design arithmetic intact."),
]

RESULTS = [
    ("h2", "Noise falls by one order of magnitude, not four, across the usable range of bin size"),
    ("p", f"Re-aggregating the native 2 \u00b5m bins by exact integer summation (Fig. 1a), from 2 to "
          f"256 \u00b5m \u2014 a {num('256/2', 0, 'int')}-fold change in d and a "
          f"{num('(256/2)**2', 0, 'int')}-fold change in bin area \u2014 pooled overdispersion fell from "
          f"{num('d.colon.alpha[2]', 1)} to {num('d.colon.alpha[256]', 2)} in colon, from "
          f"{num('d.lung.alpha[2]', 1)} to {num('d.lung.alpha[256]', 2)} in lung and from "
          f"{num('d.brain.alpha[2]', 2)} to {num('d.brain.alpha[256]', 2)} in brain. That is "
          f"{num('N.derived.colon_alpha_ratio_2_to_256', 1)}-, {num('N.derived.lung_alpha_ratio_2_to_256', 1)}- and "
          f"{num('N.derived.brain_alpha_ratio_2_to_256', 1)}-fold, against {num('(256/2)**2', 0, 'int')}-fold under "
          f"d^\u22122 (Fig. 1b). A power law describes the eight points with R\u00b2 = "
          f"{num('d.colon.r2', 3)}, {num('d.lung.r2', 3)} and {num('d.brain.r2', 3)}, giving \u03b3 = "
          f"{num('J.colon.gamma_full', 2)}, {num('J.lung.gamma_full', 2)} and "
          f"{num('J.brain.gamma_full', 2)} \u2014 {num('2/J.colon.gamma_full', 1)}-, "
          f"{num('2/J.lung.gamma_full', 1)}- and {num('2/J.brain.gamma_full', 1)}-fold shallower than counting "
          f"statistics (Fig. 1b). Removing between-domain variance by within-block estimation gives "
          f"{num('abs(d.colon.gamma_local)', 2)}, {num('abs(d.lung.gamma_local)', 2)} and "
          f"{num('abs(d.brain.gamma_local)', 2)}, so the ordering survives but the level does not (Fig. S4)."),
    ("p", f"Two properties of the measurement bound what can be claimed. First, the coarse end is thin: at "
          f"256 \u00b5m the colon section yields {num('d.colon.power[256].n_bins', 0, 'int')} bins, and only "
          f"{num('F.colon.by_scale[256].fill_frac', 1, 'pct')} of them lie fully inside the tissue mask; the "
          f"lung and brain sections are similar ({num('F.lung.by_scale[256].fill_frac', 1, 'pct')} and "
          f"{num('F.brain.by_scale[256].fill_frac', 1, 'pct')}). Second, that boundary convention is not "
          f"neutral. Recomputing with only fully-in-tissue blocks raises \u03b3 by "
          f"{num('F.colon.gamma_delta', 2)} in colon, {num('F.lung.gamma_delta', 2)} in lung and "
          f"{num('F.brain.gamma_delta', 2)} in brain (Table S5). We use the all-in-tissue convention "
          f"throughout because it uses every bin; the strict convention is reported alongside every claim that "
          f"depends on the level of \u03b3."),
    ("p", f"The exponent also depends on the fit window, and we fix that window once: 2\u2013256 \u00b5m, all "
          f"eight scales. Dropping the four finest scales and fitting from 32 \u00b5m upward would give "
          f"{num('D3.colon_gamma_from_window32', 2)}, {num('D3.lung_gamma_from_window32', 2)} and "
          f"{num('D3.brain_gamma_from_window32', 2)} instead (Fig. S1). Every comparison below uses the same "
          f"window, and \u03b3 is read as an effective exponent over 2\u2013256 \u00b5m rather than as an "
          f"asymptotic slope."),

    ("h2", "A single correlation length is rejected"),
    ("p", f"If one correlation length \u03be set the curve, \u03b1 would follow a two-parameter shape "
          f"H(d/\u03be) that flattens below \u03be and steepens above it. It does not. A quadratic in log d beats "
          f"a pure power law by \u0394AIC = {num('M.colon.daic.power', 1)}, {num('M.lung.daic.power', 1)} and "
          f"{num('M.brain.daic.power', 1)} in the three sections, with F = {num('M.colon.quad_F', 1)} "
          f"(p = {num('M.colon.quad_p_value', 4)}), {num('M.lung.quad_F', 1)} "
          f"(p = {num('M.lung.quad_p_value', 4)}) and {num('M.brain.quad_F', 1)} "
          f"(p = {num('M.brain.quad_p_value', 4)}) on one and five degrees of freedom (Fig. 2b; Table S2). Both "
          f"finite-\u03be shapes trail the quadratic by at least "
          f"{num('min(D3.colon_xi_dAIC_min,D3.lung_xi_dAIC_min,D3.brain_xi_dAIC_min)', 0, 'int')} AIC units in "
          f"every section. Forced to fit, the exponential model returns "
          f"\u03be = {num('M.colon.xi_hat', 0, 'int')} \u00b5m "
          f"[{num('M.colon.xi_ci[0]', 0, 'int')}, {num('M.colon.xi_ci[1]', 0, 'int')}], "
          f"{num('M.lung.xi_hat', 0, 'int')} \u00b5m [{num('M.lung.xi_ci[0]', 0, 'int')}, "
          f"{num('M.lung.xi_ci[1]', 0, 'int')}] and {num('M.brain.xi_hat', 0, 'int')} \u00b5m "
          f"[{num('M.brain.xi_ci[0]', 0, 'int')}, {num('M.brain.xi_ci[1]', 0, 'int')}]: three intervals that "
          f"overlap, so the three sections cannot be ordered by correlation length (Fig. 2c)."),
    ("p", f"The local slope makes the same point in a form that does not depend on any model. Observed slopes "
          f"rise across the range in both tumours \u2014 colon from "
          f"{num('abs(d.colon.segments[2:4])', 2)} at 2\u20134 \u00b5m to "
          f"{num('abs(d.colon.segments[32:64])', 2)} at 32\u201364 \u00b5m, lung from "
          f"{num('abs(d.lung.segments[2:4])', 2)} to {num('abs(d.lung.segments[64:128])', 2)} at "
          f"64\u2013128 \u00b5m \u2014 while brain stays between "
          f"{num('abs(d.brain.segments[2:4])', 2)} and {num('abs(d.brain.segments[128:256])', 2)} at "
          f"every scale (Fig. 2a). A finite-\u03be model predicts that the slope should climb steeply \u2014 to "
          f"{num('max(M.lung.local_slope.predicted_xi)', 2)} at the coarsest segment in lung \u2014 while the "
          f"observed slope never exceeds {num('max(M.lung.local_slope.observed)', 2)}. At the level of single "
          f"genes the rejection is "
          f"harder still: only {num('N.gene_level.colon.xi_fit.frac_pref_xi', 1, 'pct')}, "
          f"{num('N.gene_level.lung.xi_fit.frac_pref_xi', 1, 'pct')} and "
          f"{num('N.gene_level.brain.xi_fit.frac_pref_xi', 1, 'pct')} of genes prefer a finite \u03be to a pure "
          f"power law by more than 2 AIC units."),

    ("h2", "The exponent is a property of transcripts, and it is indexed by abundance"),
    ("p", f"Fitting each gene separately over the same eight scales gives a distribution, not a number. Among "
          f"genes that pass the filter \u2014 positive \u03b1 at all eight scales and at least 0.5 counts per "
          f"256 \u00b5m bin ({num('N.gene_level.colon.n_fit', 0, 'int')}, "
          f"{num('N.gene_level.lung.n_fit', 0, 'int')} and {num('N.gene_level.brain.n_fit', 0, 'int')} genes) "
          f"\u2014 the median \u03b3 is {num('N.gene_level.colon.gamma_gene.median', 2)}, "
          f"{num('N.gene_level.lung.gamma_gene.median', 2)} and "
          f"{num('N.gene_level.brain.gamma_gene.median', 2)}, with interquartile ranges of "
          f"{num('N.gene_level.colon.gamma_gene.q25', 2)}\u2013{num('N.gene_level.colon.gamma_gene.q75', 2)}, "
          f"{num('N.gene_level.lung.gamma_gene.q25', 2)}\u2013{num('N.gene_level.lung.gamma_gene.q75', 2)} and "
          f"{num('N.gene_level.brain.gamma_gene.q25', 2)}\u2013{num('N.gene_level.brain.gamma_gene.q75', 2)} "
          f"(Fig. 3a,b). The 95th percentile is {num('N.gene_level.colon.gamma_gene.p95', 2)}, "
          f"{num('N.gene_level.lung.gamma_gene.p95', 2)} and {num('N.gene_level.brain.gamma_gene.p95', 2)}; the "
          f"largest single value anywhere is {num('max(N.gene_level.colon.gamma_gene.gmax,N.gene_level.lung.gamma_gene.gmax,N.gene_level.brain.gamma_gene.gmax)', 2)}. "
          f"No transcript in any section reaches the counting limit. The spread is not a tail of noisy genes: "
          f"{num('N.gene_level.colon.gamma_gene.frac_gt_1', 1, 'pct')}, "
          f"{num('N.gene_level.lung.gamma_gene.frac_gt_1', 1, 'pct')} and "
          f"{num('N.gene_level.brain.gamma_gene.frac_gt_1', 1, 'pct')} of genes exceed \u03b3 = 1, and the "
          f"ratio of the 75th to the 25th percentile is {num(QR_MIN, 1)} to {num(QR_MAX, 1)} in every section."),
    ("p", f"\u03b3 falls monotonically with abundance. Spearman correlation between \u03b3_g and the mean count "
          f"per 32 \u00b5m bin is \u2212{num('abs(J.colon.gamma_vs_mu.spearman)', 2)}, "
          f"\u2212{num('abs(J.lung.gamma_vs_mu.spearman)', 2)} and "
          f"\u2212{num('abs(J.brain.gamma_vs_mu.spearman)', 2)}. In layers of one decade the median travels from "
          f"\u03b3 = {num('J.colon.gamma_vs_mu.by_decade[0].med_gamma', 2)} for genes below 0.01 counts per bin "
          f"to {num('J.colon.gamma_vs_mu.by_decade[4].med_gamma', 2)} for genes between 10 and 100 in colon; in "
          f"lung from {num('J.lung.gamma_vs_mu.by_decade[0].med_gamma', 2)} to "
          f"{num('J.lung.gamma_vs_mu.by_decade[4].med_gamma', 2)}, and in brain from "
          f"{num('J.brain.gamma_vs_mu.by_decade[0].med_gamma', 2)} to "
          f"{num('J.brain.gamma_vs_mu.by_decade[4].med_gamma', 2)} (Fig. 3d; Table 3). The three curves are offset but "
          f"parallel; the ordering of the sections is the same at every layer \u2014 brain first, colon second, "
          f"lung third, in all six (Table 3). A two-way decomposition of those eighteen values puts "
          f"{num('V.pooled.frac_layer', 1, 'pct')} of their variance on the layer, "
          f"{num('V.pooled.frac_tissue', 1, 'pct')} on the section and "
          f"{num('V.pooled.frac_resid', 1, 'pct')} on their interaction."),
    ("p", f"Three controls establish that the gradient is a property of the data and not of the estimator, and a "
          f"fourth \u2014 leaving one transcript out \u2014 is reported with the pooled value below. Depth: "
          f"binomial thinning of the colon section to {num('N.depth.runs[PP(0.46)].p', 0, 'pct')} of its UMIs "
          f"changes \u03b3 by {num('N.depth.runs[PP(0.46)].gamma', 4)} against "
          f"{num('N.depth.runs.full.gamma', 4)}, a shift of "
          f"{num('abs(N.depth.runs[PP(0.46)].gamma-N.depth.runs.full.gamma)', 5)}; at "
          f"{num('N.depth.runs[PP(0.1)].p', 0, 'pct')} of depth the shift is "
          f"{num('abs(N.depth.runs[PP(0.1)].gamma-N.depth.runs.full.gamma)', 5)} (Fig. 5c). This is algebraic: "
          f"binomial thinning multiplies numerator and denominator by the same p\u00b2, so \u03b1 is invariant. "
          f"Synthetic abundance: on a correlated field with a known exponent of "
          f"{num('S.expression.beta_true', 2)} and pure Poisson counting over "
          f"{num('S.expression.n_genes', 0, 'int')} genes spanning four orders of magnitude, the recovered "
          f"\u03b3 is flat at {num('S.expression.gamma_all.med', 2)} "
          f"across layers (range {num('D3.synth_range_across_layers', 3)}; Fig. 3d, dashed; Table S3), so the "
          f"estimator does not manufacture an abundance dependence. Reliability: between the two halves of one "
          f"section the per-gene exponents correlate at r = {num('N.splithalf.reliability.pearson', 3)} "
          f"(n = {num('N.splithalf.reliability.n', 0, 'int')}, MAD = "
          f"{num('N.splithalf.reliability.mad', 3)}; Fig. 3c), and the reliability holds inside every expression "
          f"layer (Table S6)."),
    ("p", f"The synthetic control also calibrates the level. It recovers "
          f"{num('S.expression.gamma_all.med', 2)} when the input is {num('S.expression.beta_true', 2)} "
          f"(Fig. 3d, dashed): a "
          f"bias of {num('D3.synth_bias', 2)} caused by fitting a single power law to a curve that is still "
          f"curving inside the window. The absolute value of \u03b3 therefore carries only relative meaning, "
          f"which is one more reason not to compare \u201cthe exponent\u201d of one study with another. "
          f"Comparisons within a window, between genes or between regions, are unaffected."),
    ("p", f"A last control asks whether the exponent is a restatement of spatial clustering. It is not, and the "
          f"three sections disagree about how much texture contributes. In colon and brain the per-gene exponent "
          f"falls as Moran's I rises, and the association survives when abundance is held fixed: inside twenty "
          f"equal-count bands the median correlation is "
          f"\u2212{num('abs(N.spatial_proxy.by_tissue.colon.band_rho[1].median)', 2)} in colon and "
          f"\u2212{num('abs(N.spatial_proxy.by_tissue.brain.band_rho[1].median)', 2)} in brain, and a "
          f"semiparametric partial correlation that removes a binned mean of log abundance gives "
          f"\u2212{num('abs(N.spatial_proxy.by_tissue.colon.semiparam.pearson)', 2)} and "
          f"\u2212{num('abs(N.spatial_proxy.by_tissue.brain.semiparam.pearson)', 2)} (Table S9). In lung the "
          f"abundance gradient is present but texture adds nothing: the within-band median is "
          f"\u2212{num('abs(N.spatial_proxy.by_tissue.lung.band_rho[1].median)', 2)} and the partial "
          f"\u2212{num('abs(N.spatial_proxy.by_tissue.lung.semiparam.pearson)', 2)}. Where the association "
          f"appears its sign is the one the geometry requires: coarse structure survives binning, so its "
          f"overdispersion decays slowly and its \u03b3 is small, while structure at the finest scale is averaged "
          f"away and its \u03b3 is large. A gradient that is present in all three sections, while the texture "
          f"association is present in two, cannot be reduced to the amount of clustering measured at one lag. A "
          f"fitted correlation length cannot arbitrate either, because "
          f"{num('N.spatial_proxy.by_tissue.colon.xi_frac_at_floor', 1, 'pct')} of colon genes return a range at "
          f"the lowest admissible value \u2014 the variogram is bounded rather than estimated \u2014 and the range "
          f"carries no information about \u03b3 once abundance is controlled (partial \u03c1 = "
          f"{num('N.spatial_proxy.by_tissue.colon.partial_rho_xi', 2)}). Abundance and spatial texture are "
          f"separate axes, and only the first is shared by the three sections."),

    ("h2", "The pooled exponent is carried by a handful of genes, and it moves between regions of one section"),
    ("p", f"Pooled \u03b1 weights every gene by \u03bc\u00b2, so the pooled exponent is a weighted average in "
          f"which a few transcripts dominate. The effective number of genes, 1/\u03a3w\u00b2, is "
          f"{num('N.gene_level.colon.composition.n_eff', 1)} in colon, {num('N.gene_level.lung.composition.n_eff', 1)} in lung and "
          f"{num('N.gene_level.brain.composition.n_eff', 1)} in brain. The top 100 genes carry "
          f"{num('N.gene_level.colon.composition.top100_share', 1, 'pct')}, "
          f"{num('N.gene_level.lung.composition.top100_share', 1, 'pct')} and "
          f"{num('N.gene_level.brain.composition.top100_share', 1, 'pct')} of \u03a3\u03bc\u00b2; half of the weight sits on "
          f"{num('N.gene_level.colon.composition.n_for_50pct', 0, 'int')}, {num('N.gene_level.lung.composition.n_for_50pct', 0, 'int')} "
          f"and {num('N.gene_level.brain.composition.n_for_50pct', 0, 'int')} genes and nine tenths on "
          f"{num('N.gene_level.colon.composition.n_for_90pct', 0, 'int')}, "
          f"{num('N.gene_level.lung.composition.n_for_90pct', 0, 'int')} and "
          f"{num('N.gene_level.brain.composition.n_for_90pct', 0, 'int')} (Fig. 4a)."),
    ("p", f"At every scale the pooled quantity is a convex combination of the per-gene ones. Because "
          f"\u03b1_pool(d) = \u03a3_g w_g(d)\u03b1_g(d) with w_g = \u03bc_g\u00b2/\u03a3\u03bc\u00b2, recomputing "
          f"the pooled value gene by gene returns the direct estimator to machine precision. Fitting the eight "
          f"scales of that sum gives {num('N.decomp.by_tissue.colon.gamma_weighted', 4)}, "
          f"{num('N.decomp.by_tissue.lung.gamma_weighted', 4)} and "
          f"{num('N.decomp.by_tissue.brain.gamma_weighted', 4)} against the reported "
          f"{num('J.colon.gamma_full', 3)}, {num('J.lung.gamma_full', 3)} and {num('J.brain.gamma_full', 3)}, "
          f"residuals of {num('N.decomp.min_rel_repro', 1, 'pct')} to "
          f"{num('N.decomp.max_rel_repro', 1, 'pct')} (Fig. S6a; Table S7). The residual has two sources, both bookkeeping: "
          f"genes without a computable exponent are absent from the mean, and the weights drift as d grows, which "
          f"accounts for {num('N.decomp.by_tissue.colon.drift_share', 0, 'pct')} of the distance between the "
          f"weighted and the unweighted summary in colon and "
          f"{num('N.decomp.by_tissue.lung.drift_share', 0, 'pct')} in lung. The pooled exponent is a weighted "
          f"mean of the per-gene exponents, and which transcripts carry the weights is what fixes its value."),
    ("p", f"The weights that decide the exponent are not the ones that decide the variance. Gene g's share of "
          f"the pooled numerator is v_g(d) = w_g(d)\u03b1_g(d)/\u03b1_pool(d), and under those weights the "
          f"effective number of genes is {num('N.decomp.by_tissue.colon.n_eff_alpha', 1)}, "
          f"{num('N.decomp.by_tissue.lung.n_eff_alpha', 1)} and "
          f"{num('N.decomp.by_tissue.brain.n_eff_alpha', 1)} \u2014 below the "
          f"{num('N.gene_level.colon.composition.n_eff', 1)}, {num('N.gene_level.lung.composition.n_eff', 1)} "
          f"and {num('N.gene_level.brain.composition.n_eff', 1)} implied by \u03a3\u03bc\u00b2 by a factor of "
          f"{num('N.decomp.min_shrink', 1)} to {num('N.decomp.max_shrink', 1)}. The reason is the gradient "
          f"itself: \u03b3 rises as abundance falls, so a scarce transcript contributes more to the exponent "
          f"than its share of signal suggests. In colon IGKC holds {num('J.colon.top1_weight', 1, 'pct')} of "
          f"\u03a3\u03bc\u00b2 but {num('N.decomp.by_tissue.colon.top1_alpha_share', 1, 'pct')} of the exponent, "
          f"and the three heaviest contributors hold "
          f"{num('N.decomp.by_tissue.colon.top3_alpha_share', 1, 'pct')} (Table S7); in lung the same transcript "
          f"moves from {num('J.lung.top1_weight', 1, 'pct')} to "
          f"{num('N.decomp.by_tissue.lung.top1_alpha_share', 1, 'pct')}. In brain the two orderings disagree "
          f"outright: the heaviest transcript by \u03bc\u00b2 is Camk2n1, while the largest single "
          f"contributor to the exponent is Ttr, at "
          f"{num('N.decomp.by_tissue.brain.top1_alpha_share', 1, 'pct')}."),
    ("p", f"The consequence is visible when genes are added back one at a time in order of weight. In colon, one "
          f"gene gives \u03b3 = {num('J.colon.cumulative.gamma[0]', 3)} and all "
          f"{num('J.colon.n_genes', 0, 'int')} give {num('J.colon.cumulative.gamma[-1]', 3)}: adding "
          f"the remaining genes moves the answer by "
          f"{num('abs(J.colon.cumulative.gamma[-1]-J.colon.cumulative.gamma[0])', 3)}. In brain the same single "
          f"gene gives {num('J.brain.cumulative.gamma[0]', 3)} and the full set "
          f"{num('J.brain.cumulative.gamma[-1]', 3)}; it takes "
          f"{num('J.brain.cumulative.k[10]', 0, 'int')} genes to reach "
          f"{num('J.brain.cumulative.gamma[10]', 3)} (Fig. 4b). Leaving one gene out confirms the point from the "
          f"other side: dropping IGKC lowers colon by {num('J.colon.loo.max_abs_delta', 3)} and lung by "
          f"{num('J.lung.loo.max_abs_delta', 3)}, dropping Camk2n1 raises brain by "
          f"{num('J.brain.loo.max_abs_delta', 3)}, and in colon and lung exactly "
          f"{num('J.colon.loo.n_delta_gt_0p05', 0, 'int')} of "
          f"{num('J.colon.n_genes', 0, 'int')} genes moves the exponent by more than 0.05 (Fig. 4c). The same "
          f"reported number therefore means different things: in colon \u201c\u03b3 = 0.55\u201d is very nearly "
          f"the exponent of one immunoglobulin transcript, and in brain the same number is an average over "
          f"{num('N.decomp.by_tissue.brain.n_eff_alpha', 0, 'int')} effective genes."),
    ("p", f"The rank order of the sections reverses depending on which of the two statistics is quoted. Pooled, they "
          f"order colon ({num('J.colon.gamma_full', 3)}), lung ({num('J.lung.gamma_full', 3)}), brain "
          f"({num('J.brain.gamma_full', 3)}). Per gene, they order brain "
          f"({num('N.gene_level.brain.gamma_gene.median', 2)}), colon "
          f"({num('N.gene_level.colon.gamma_gene.median', 2)}), lung "
          f"({num('N.gene_level.lung.gamma_gene.median', 2)}), and brain is the steepest in all six abundance "
          f"layers (Table 3). The section with the flattest pooled exponent therefore has the steepest genes, "
          f"and the section with the second-highest pooled exponent has the flattest genes. The sign of the "
          f"difference between the two statistics follows the heaviest transcript: the pooled value sits "
          f"{num('D3.colon_pooled_minus_median', 2)} above the per-gene median in colon and "
          f"{num('D3.lung_pooled_minus_median', 2)} above it in lung, where IGKC carries "
          f"{num('J.colon.top1_weight', 1, 'pct')} and {num('J.lung.top1_weight', 1, 'pct')} of \u03a3\u03bc\u00b2, "
          f"and {num('abs(D3.brain_pooled_minus_median)', 2)} below it in brain, where Camk2n1 carries "
          f"{num('J.brain.top1_weight', 1, 'pct')}. The two transcripts differ in the way the arithmetic "
          f"requires: IGKC has an exponent of {num('J.colon.top1_gamma', 2)} in colon and Camk2n1 one of "
          f"{num('J.brain.top1_gamma', 2)} in brain. An ordering of sections by pooled \u03b3 is therefore not "
          f"an ordering by the steepness of their transcripts."),
    ("p", f"Pooled \u03b3 is also not stable within a section. Cutting the colon section into four column bands "
          f"on the 256 \u00b5m grid \u2014 so that no coarse bin is divided \u2014 and refitting each band gives "
          f"{num('N.splithalf.bands.q1.gamma', 3)}, {num('N.splithalf.bands.q2.gamma', 3)}, "
          f"{num('N.splithalf.bands.q3.gamma', 3)} and {num('N.splithalf.bands.q4.gamma', 3)}, a range of "
          f"{num('N.splithalf.quarter_range', 3)} and a standard deviation of "
          f"{num('N.splithalf.quarter_sd', 3)} (Fig. 5a; Table S4). The same split applied to the lung and brain "
          f"sections gives {num('N.splithalf_tissues.by_tissue.lung.quarter_sd', 3)} and "
          f"{num('N.splithalf_tissues.by_tissue.brain.quarter_sd', 3)}, so the colon value is not an outlier "
          f"among sections. Two halves cut through the centroid give "
          f"{num('N.splithalf.bands.halfA.gamma', 3)} and {num('N.splithalf.bands.halfB.gamma', 3)}. Resampling "
          f"superblocks of 512, 1024 and 2048 \u00b5m gives an SD of "
          f"{num('B.colon.results[512].sd', 3)}, {num('B.colon.results[1024].sd', 3)} and "
          f"{num('B.colon.results[2048].sd', 3)} (Fig. 5d); that resampling measures local heterogeneity only "
          f"and stays below the spread across bands, so a bootstrap over tiles would understate the uncertainty "
          f"of one section. For "
          f"comparison, the difference between the colon and brain sections is "
          f"{num('D3.tissue_gap_colon_brain', 3)}. Treating the regional SD as a lower bound on the "
          f"between-section SD \u2014 one section cannot supply more \u2014 the 95% half-width for a two-section "
          f"comparison is {num('D3.halfwidth95_regional_2sec', 3)}, so the tissue gap clears the threshold by a "
          f"factor of {num('D3.gap_over_regional_halfwidth', 2)} (Fig. 5b). A second partition into four "
          f"2\u00d72 quadrants, which does cut coarse blocks at the boundary, gives a wider range "
          f"({num('D3.colon_quad_min', 2)}\u2013{num('D3.colon_quad_max', 2)}, SD "
          f"{num('D3.colon_quad_sd', 3)}; Table 2) and a half-width of "
          f"{num('1.96*sqrt(2)*D3.colon_quad_sd', 3)}, against which the gap clears by only "
          f"{num('D3.tissue_gap_colon_brain/(1.96*sqrt(2)*D3.colon_quad_sd)', 2)}. Using the segmentation-based rms "
          f"instead ({num('d.colon.rms_gamma_cut', 3)}) the half-width is "
          f"{num('1.96*sqrt(2)*d.colon.rms_gamma_cut', 3)}, and the gap does not clear it at all. Under the "
          f"strict fill convention the gap narrows to "
          f"{num('F.colon.gamma_full-F.brain.gamma_full', 3)}. Whichever estimate is used, a single section per "
          f"tissue cannot establish a difference in exponent between tissues."),
    ("p", f"The two statements fit together, and the decomposition above says how. Gene-level exponents are "
          f"reproducible (r = {num('N.splithalf.reliability.pearson', 2)}) while the sample-level summary is not "
          f"(\u03c3 = {num('N.splithalf.quarter_sd', 3)} between bands of one section; Fig. 3c and Fig. 5a). A "
          f"weighted mean inherits the instability of its weights, and the weights \u2014 which transcripts "
          f"happen to be most abundant \u2014 differ between regions and between patients [[wu2021]]. That is arithmetic, "
          f"not a difference in the steepness of the transcripts themselves, and it is the reason the reported "
          f"number does not travel."),

    ("h2", "A scalar that travels between regions"),
    ("p", f"Cutting the colon section into the same four bands and recomputing five candidate summaries inside "
          f"each ranks them clearly (Table S8; Fig. S6b). Pooled \u03b3 has a standard deviation of "
          f"{num('N.robust.summary.pooled.sd', 3)} across bands and a coefficient of variation of "
          f"{num('N.robust.summary.pooled.cv', 1, 'pct')}; the per-gene mean has "
          f"{num('N.robust.summary.mean.sd', 3)} ({num('N.robust.summary.mean.cv', 1, 'pct')}), the trimmed "
          f"mean {num('N.robust.summary.trimmed.sd', 3)} ({num('N.robust.summary.trimmed.cv', 1, 'pct')}) and "
          f"the per-gene median {num('N.robust.summary.median.sd', 3)} "
          f"({num('N.robust.summary.median.cv', 1, 'pct')}). A \u03bc\u00b2-weighted median is the worst of the "
          f"five at {num('N.robust.summary.wmedian.sd', 3)} "
          f"({num('N.robust.summary.wmedian.cv', 1, 'pct')}), because it inherits the instability of the "
          f"weights without averaging over genes. The per-gene mean and the trimmed mean are indistinguishable "
          f"at this precision ({num('N.derived4.robust_sd.mean', 4)} against "
          f"{num('N.derived4.robust_sd.trimmed', 4)}); between regions of one section the per-gene mean is the "
          f"better summary, because it discards no genes and needs no tuning constant. It is "
          f"{num('N.robust.ratio_pooled_over_best', 1)} times more stable than the pooled value across the four "
          f"bands and remains a property of the transcripts. Stability between regions is not the same "
          f"property as stability between sections, and the next section tests the second."),

    ("h2", "The pooled exponent replicates between sections; the per-gene summary does not"),
    ("p", f"That test needs a second section of the same tissue carried through the identical pipeline. "
          f"Visium_HD_Human_Colon_Cancer_P1 is such a section: a different block face of the same "
          f"tumour type, released separately. We processed it with the same scripts and the same filters and "
          f"added it to the three original sections. Precision is now estimated by a spatial block bootstrap: "
          f"the section is cut into {num('MU.block_um', 0, 'int')} \u00b5m superblocks, "
          f"{num('MU.n_boot', 0, 'int')} resamples are drawn with replacement, and \u03b3 is refitted inside "
          f"each. Superblocks keep neighbouring bins together, so the estimate carries the spatial "
          f"correlation that a bin-level bootstrap would discard."),
    ("p", f"Across the four sections pooled \u03b3 is "
          f"{num('MU.sections.colon.gamma_pool', 3)} \u00b1 {num('MU.sections.colon.boot_sd', 3)} (colon), "
          f"{num('MU.sections.Visium_HD_Human_Colon_Cancer_P1.gamma_pool', 3)} \u00b1 "
          f"{num('MU.sections.Visium_HD_Human_Colon_Cancer_P1.boot_sd', 3)} (the same tissue, second section), "
          f"{num('MU.sections.lung.gamma_pool', 3)} \u00b1 {num('MU.sections.lung.boot_sd', 3)} (lung) and "
          f"{num('MU.sections.brain.gamma_pool', 3)} \u00b1 {num('MU.sections.brain.boot_sd', 3)} (brain); "
          f"the intervals are in Fig. 8a, the eight points behind each in Fig. 8b and Table S11. The two colon "
          f"sections differ by "
          f"{num('abs(MU.contrasts.same_tissue_colon.delta)', 3)}, which is "
          f"{num('abs(MU.contrasts.same_tissue_colon.z)', 1)} standard errors, while colon and brain differ by "
          f"{num('MU.contrasts.colon_vs_brain.delta', 3)}, or "
          f"{num('MU.contrasts.colon_vs_brain.z', 1)} standard errors. The smallest difference a design with "
          f"two sections per group could detect at 80% power is "
          f"{num('MU.mde.value', 3)}, so the tissue gap is "
          f"{num('MU.contrasts.colon_vs_brain.delta/MU.mde.value', 1)} times that minimum and the section gap "
          f"is a fraction of it."),
    ("p", f"The per-gene summary behaves in the opposite way. The median exponent is "
          f"{num('MU.pergene.median_p0', 3)} in the first colon section and "
          f"{num('MU.pergene.median_p1', 3)} in the second, a difference of "
          f"{num('MU.pergene.shift_overall', 3)} \u2014 "
          f"{num('MU.pergene.shift_overall/abs(MU.contrasts.same_tissue_colon.delta)', 0, 'int')} times the "
          f"difference in the pooled value between the same two sections (Fig. 8c, Table S12). The gap is not "
          f"a level offset that one calibration could absorb: the two per-gene distributions cross. Genes in "
          f"the lowest abundance decile have a median exponent "
          f"{num('MU.pergene.shift_bottom', 3)} higher in the second section, genes in the highest decile "
          f"{num('MU.pergene.shift_top', 3)} lower, and the interval for the shift excludes zero in "
          f"{num('MU.pergene.n_deciles_ci_excludes_zero', 0, 'int')} of the ten deciles. The second section "
          f"also carries {num('MU.pergene.umi_p0/MU.pergene.umi_p1', 1)} times fewer counts "
          f"(Table S11), and a shallower section inflates the exponent of poorly detected genes, whose "
          f"small-scale overdispersion is sampling noise. We therefore read the per-gene difference as a "
          f"property of the measurement rather than of the tissue, and do not correct it."),
    ("p", f"Two consequences follow. First, one section is enough to estimate the pooled exponent of that "
          f"section, and the estimate is reproducible: the between-section difference is "
          f"{num('abs(MU.contrasts.same_tissue_colon.delta)', 3)} where the between-tissue difference is "
          f"{num('MU.contrasts.colon_vs_brain.delta', 3)}. Second, the recommendation to report a per-gene "
          f"summary instead of a pooled value does not survive this test. The per-gene mean is the more stable "
          f"summary inside a section, where the gene set is held fixed, and the less stable one between "
          f"sections, where it is not. Neither scalar is portable on its own: the per-gene distribution has to "
          f"be reported, and Table 2 is arranged for that."),

    ("h2", "The same composition under a different capture chemistry"),
    ("p", f"A second platform tests whether the composition belongs to the tissue or to the assay. "
          f"Stereo-seq C04042E2 is a mouse whole-brain section captured in situ on a 0.5 \u00b5m DNA "
          f"nanoball array, with a chemistry that shares no probe design with Visium HD. We aggregated its "
          f"native 0.5 \u00b5m DNBs into square blocks of side d = 2, 4, \u2026, 256 \u00b5m, the same "
          f"eight-scale ladder as the Visium sections, and applied the identical filter and estimator to the "
          f"resulting moments."),
    ("p", f"The Stereo-seq section gives \u03b3 = {num('X.stereo.gamma_pool', 3)} against "
          f"{num('X.visium.gamma_pool', 3)} for the Visium HD mouse brain section, and a per-gene median of "
          f"{num('X.stereo.gamma_median', 3)} against {num('X.visium.gamma_median', 3)}. The absolute level "
          f"therefore does not travel between platforms, as expected if the exponent is a composition "
          f"statistic whose weights depend on capture efficiency (Fig. 7a). What travels is the structure. The "
          f"ordering "
          f"that defines the statistic holds: the pooled value lies far below the per-gene median "
          f"({num('X.stereo.gamma_pool', 3)} against {num('X.stereo.gamma_median', 3)}), as it does on "
          f"Visium HD brain. Concentration holds: {num('X.stereo.n_eff_alpha', 2)} effective genes out of "
          f"{num('X.stereo.n_genes_ok', 0, 'int')} tested \u2014 "
          f"{num('X.concentration.n_eff_over_genes', 2, 'pct')} \u2014 with "
          f"{num('X.stereo.top1_alpha_share', 1, 'pct')} of the pooled exponent carried by one transcript. On "
          f"both platforms the top 100 genes carry "
          f"{num('X.top100_share_stereo', 0, 'pct')} and {num('X.top100_share_visium', 0, 'pct')} of it "
          f"(Fig. 7b)."),
    ("p", f"The heaviest transcript is Gm42418, an rRNA-derived locus routinely discarded as contamination, so "
          f"the concentration result could be an artefact of one junk gene. It is not. Dropping the "
          f"rRNA-derived loci (Gm42418, Gm26917) leaves the top gene as Ttr with "
          f"{num('X.without_rrna.top1_share', 1, 'pct')} of the exponent and "
          f"{num('X.without_rrna.n_eff_alpha', 2)} effective genes, and both concentration criteria still "
          f"pass. Ttr is also the heaviest transcript in the Visium HD brain section, so after removing the "
          f"contaminant the two platforms agree on which transcript dominates \u2014 a stronger result than "
          f"the criterion written in advance."),
    ("p", f"The per-gene exponents themselves replicate. Over "
          f"{num('X.replication.n', 0, 'int')} genes measured on both platforms, \u03b3_g has a rank "
          f"correlation of {num('X.replication.rho_raw', 3)}. Both platforms also rank genes by abundance in "
          f"the same order (\u03c1 = {num('X.replication.rho_abundance', 3)}), and abundance alone could "
          f"manufacture a correlation of that size. Removing it does not: after rank-residualising each "
          f"platform's exponent on its own log abundance, the correlation is still "
          f"{num('X.replication.rho_partial_abundance', 3)} (Fig. 7c). The agreement is in the per-gene "
          f"structure, not in a shared ranking by expression level."),
    ("p", f"Read together, the two platforms separate what the field conflates. The number \u03b3 is "
          f"platform-specific \u2014 {num('X.delta_gamma_pool', 3)} apart for the same tissue on two "
          f"chemistries \u2014 while the composition that produces it is not. A title that speaks of binned "
          f"spatial transcriptomics rather than of Visium HD is, on this evidence, the accurate one."),


    ("h2", "What follows for design"),
    ("p", f"The design consequences do not depend on the decomposition above; they need only the fact that "
          f"every observed \u03b3 is at most {num('max(N.gene_level.colon.gamma_gene.gmax,N.gene_level.lung.gamma_gene.gmax,N.gene_level.brain.gamma_gene.gmax)', 2)}. "
          f"The coefficient of variation of a bin's total count needs no fit and falls in colon from "
          f"{num('d.colon.cv[2]', 0, 'pct')} at 2 \u00b5m to "
          f"{num('d.colon.cv[256]', 1, 'pct')} at 256 \u00b5m, against a Poisson prediction of "
          f"{num('d.colon.cv_poisson[2]', 1, 'pct')} to {num('d.colon.cv_poisson[256]', 3, 'pct')} "
          f"(Fig. 6a). At the coarsest scale the observed CV is "
          f"{num('N.derived.colon_obs_over_poisson_256', 0, 'int')}-fold the Poisson value in colon, "
          f"{num('N.derived.lung_obs_over_poisson_256', 0, 'int')}-fold in lung and "
          f"{num('N.derived.brain_obs_over_poisson_256', 0, 'int')}-fold in brain. It is a different "
          f"statistic, not a second estimate of the exponent: it averages per-gene overdispersion with weights "
          f"equal to the squared share of the total count, whereas pooled \u03b1 weights by the squared share of "
          f"\u03a3\u03bc\u00b2, which in colon puts {num('J.colon.top1_weight', 1, 'pct')} of the weight on one "
          f"transcript. The exponent implied by its own decline is {num('D3.colon_cv_gamma', 2)}, "
          f"{num('D3.lung_cv_gamma', 2)} and {num('D3.brain_cv_gamma', 2)}, below the \u03b1-based value in "
          f"every section, and we use it only as a design quantity. Sequencing deeper cannot move this "
          f"curve, because depth enters only through the Poisson term, which at 2 \u00b5m is "
          f"{num('d.colon.cv_poisson[2]/d.colon.cv[2]', 2)} of the observed CV in colon."),
    ("p", f"Per-gene estimates peak in the middle of the range, never at the native resolution. The median "
          f"relative error of a gene's \u03b1 is "
          f"{num('d.colon.power[2].rel_err', 0, 'pct')} at 2 \u00b5m in colon, "
          f"{num('d.lung.power[2].rel_err', 0, 'pct')} in lung and "
          f"{num('d.brain.power[2].rel_err', 0, 'pct')} in brain, falling to "
          f"{num('d.colon.power[32].rel_err', 0, 'pct')}, "
          f"{num('d.lung.power[128].rel_err', 0, 'pct')} and "
          f"{num('d.brain.power[32].rel_err', 0, 'pct')} at 32, 128 and 32 \u00b5m respectively, and rising "
          f"again past that point in colon and brain (Table S1)."),
    ("p", f"For comparisons, the currency is sections. With one section per group the minimum detectable "
          f"difference in \u03b3 at 80% power is {num('d.colon.mde_gamma[n1]', 3)}, "
          f"{num('d.lung.mde_gamma[n1]', 3)} and {num('d.brain.mde_gamma[n1]', 3)}; with three it is "
          f"{num('d.colon.mde_gamma[n3]', 3)}, {num('d.lung.mde_gamma[n3]', 3)} and "
          f"{num('d.brain.mde_gamma[n3]', 3)} (Fig. 6b). The "
          f"{num('D3.tissue_gap_colon_lung', 3)} difference between the two human tumours "
          f"would need {num('D3.colon_n_for_delta_gap', 0, 'int')} sections per group in colon; a difference of "
          f"0.05 needs {num('D3.colon_n_for_delta005', 0, 'int')} and one of 0.02 needs "
          f"{num('D3.colon_n_for_delta002', 0, 'int')}. Bin size does not substitute: at 2 \u00b5m the minimum "
          f"detectable fold-change in pooled \u03b1 between two sections is "
          f"{num('d.colon.power[2].mde_1v1', 2)} in colon and at 256 \u00b5m it is still "
          f"{num('d.colon.power[256].mde_1v1', 2)}."),
    ("p", "Sample-size calculations of this form already exist for other spatial assays [[geomex2023]]. The "
          "difference here is the variance that enters them, plotted against sections per group in Fig. 6b. It is "
          "not the Poisson variance of a count, which "
          "shrinks with bin area, but the spread of \u03b3 between regions of one section, which does not."),
    ("p", f"Ranking the manipulations by how far they move \u03b3 puts the design levers in order: depth "
          f"{num('D3.depth_delta_max', 4)}; dropping the heaviest gene {num('J.colon.loo.max_abs_delta', 3)}; "
          f"changing region within one section {num('N.splithalf.quarter_range', 3)}; changing the fit window "
          f"{num('D3.lung_window_delta', 3)} (lung, the largest of the three); changing expression layer "
          f"{num('V.pooled.range_layer', 3)} (Fig. 6c). Three of those five are analyst choices, not biology. "
          f"None of the five approaches the factor of "
          f"{num('2/min(J.colon.gamma_full,J.lung.gamma_full,J.brain.gamma_full)', 2)} that separates the "
          f"measured exponent from the counting limit \u2014 which is why the design conclusion holds "
          f"regardless of how the attribution question is settled."),
]

DISCUSSION = [
    ("p", f"Binned spatial transcriptomics does not behave the way its defaults assume. Over a "
          f"{num('256/2', 0, 'int')}-fold range of bin size the noise falls as d^\u2212"
          f"{num('J.colon.gamma_full', 2)} in a human colon tumour, d^\u2212{num('J.lung.gamma_full', 2)} in a "
          f"human lung tumour and d^\u2212{num('J.brain.gamma_full', 2)} in mouse brain, not as d^\u22122. A "
          f"{num('(256/2)**2', 0, 'int')}-fold increase in bin area therefore buys a "
          f"{num('N.derived.brain_alpha_ratio_2_to_256', 1)}- to {num('N.derived.colon_alpha_ratio_2_to_256', 1)}-fold "
          f"reduction in overdispersion rather than a four-order-of-magnitude one. But the more consequential "
          f"finding is not the size of the shortfall. It is that the shortfall is not a constant of the tissue, "
          f"and that most of what a section reports is set by which transcripts dominate it."),
    ("p", f"Three independent views say the same thing. Per gene, the exponent varies "
          f"{num(RATIO_MIN, 1)}\u2013{num(RATIO_MAX, 1)}-fold between the 5th and 95th percentile within a "
          f"single section and decreases monotonically with abundance. The pooled value is a "
          f"\u03bc\u00b2-weighted average of those per-gene values, with an effective gene count between "
          f"{num('N.decomp.min_neff_alpha', 1)} and {num('N.decomp.max_neff_alpha', 0, 'int')} once the "
          f"weights are those the exponent actually carries; "
          f"and the average moves by {num('N.splithalf.quarter_range', 3)} between four regions of one section. "
          f"A quantity whose value is set by which transcripts happen to dominate \u03a3\u03bc\u00b2 is a "
          f"composition statistic. Two studies that report \u201cthe exponent of tissue X\u201d are comparing "
          f"their cell-type mixtures and library composition, and only in part their architecture. In our own "
          f"data the "
          f"two human tumours differ by {num('D3.tissue_gap_colon_lung', 3)}, which is smaller than the "
          f"change produced by dropping one gene in either of them. The gradient has a spatial counterpart in "
          f"two of the three sections and not in the third (Table S9), so abundance is not simply a proxy for "
          f"how coarse a transcript's pattern is."),
    ("p", "That reading also disposes of a model the field has been implicitly using. If a tissue had a "
          "characteristic correlation length, a two-parameter shape H(d/\u03be) would fit better than a "
          "flexible curve in log d and the fitted \u03be would separate tissues. Neither happens: the quadratic "
          "wins by 15 to 18 AIC units over a pure power law and the finite-\u03be shapes lose by 20 to 64, "
          "while the three \u03be intervals overlap over most of their range. What the curve shows is a "
          "continuum of correlation, consistent with the shallow variance\u2013grain slopes that ecology has "
          "measured for decades [[wiens1989]] and with the scale dependence that the modifiable areal unit "
          "problem describes for any spatially aggregated statistic [[jelinski1996]]. Spatial autocorrelation "
          "reduces the effective sample size [[legendre1993]], and in a binned section the reduction is severe: "
          "the design effect of half-section differences relative to random halves is "
          f"{num('d.colon.deff_gamma', 1)} in colon and {num('d.brain.deff_gamma', 1)} in brain."),
    ("p", f"For an analyst the reporting rule is short, and it can be applied to existing tables. First, do not "
          f"report a single exponent for a sample; report the per-gene distribution \u2014 median, interquartile "
          f"range and the fraction of genes above 1 \u2014 which is reproducible between halves of a section. Quote "
          f"its mean or its trimmed mean rather than the pooled value, which varies "
          f"{num('N.robust.ratio_pooled_over_best', 1)} times as much between regions of one section (Table S8). "
          f"Second, report an effective gene count and say which weights it uses: 1/\u03a3w\u00b2 counts the "
          f"variance weights, 1/\u03a3v\u00b2 the weights the exponent itself carries, and the two differ by up "
          f"to {num('N.decomp.max_shrink', 1)}-fold. Anything below ten means the pooled value is a single-gene "
          f"measurement. Third, report at least two spatially separated regions per section, or state the "
          f"regional SD, because that is the only internal estimate of how far the number travels. A paper that "
          f"gives median \u03b3, its IQR, an effective count and a regional range can be compared with another "
          f"such paper; a paper that gives one exponent cannot."),
    ("p", "The design conclusions survive all of this because they rest on an inequality rather than on an "
          f"attribution. Every transcript measured here has \u03b3 at most "
          f"{num('max(N.gene_level.colon.gamma_gene.gmax,N.gene_level.lung.gamma_gene.gmax,N.gene_level.brain.gamma_gene.gmax)', 2)}, "
          "so merging bins cannot deliver "
          "the variance reduction that counting statistics promise, and the discrepancy grows with the area "
          "merged: the four-order-of-magnitude gain that motivates coarse binning is unavailable at any bin "
          "size. Depth is nearly free to vary \u2014 binomial thinning to a tenth of the UMIs moves \u03b3 by "
          f"{num('abs(N.depth.runs[PP(0.1)].gamma-N.depth.runs.full.gamma)', 4)} \u2014 so samples of different "
          "depth can be compared on \u03b3 without correction, "
          "which is a positive methodological result in its own right. What remains expensive is biological "
          "replication, and the minimum detectable difference falls only as 1/\u221an in the number of "
          "sections. The practical instruction is therefore to spend on sections and on per-gene reporting, "
          "not on depth or on coarser bins."),
    ("p", f"Six limits bound the claims. Three tissues have one section each, and the one tissue with two "
          f"gives a difference of {num('abs(MU.contrasts.same_tissue_colon.delta)', 3)} \u2014 small enough to "
          f"say the pooled value reproduces, not small enough to replace replication. The regional SD remains a "
          f"lower bound on between-patient variation and every design number derived from it stays optimistic. "
          f"The tumours are cancers and the brain sample is normal, so species and tissue state are confounded "
          f"and we do not attribute the between-tissue offset to either. The absolute level of \u03b3 does not "
          f"survive a change of platform \u2014 {num('X.delta_gamma_pool', 3)} between the two brain sections "
          f"below \u2014 and also depends on the fit window, as the synthetic calibration shows "
          f"({num('S.expression.gamma_all.med', 2)} recovered for an input of "
          f"{num('S.expression.beta_true', 2)}). The level of \u03b3 also depends on how boundary bins are "
          f"counted, by {num('F.brain.gamma_delta', 2)} in the worst case here. Fifth, the 32 \u00b5m grid "
          f"resolves the spatial statistics only at scales of tens of micrometres and above: "
          f"{num('N.spatial_proxy.by_tissue.colon.xi_frac_at_floor', 1, 'pct')} of variogram ranges in the colon "
          f"section sit on that bound, so the texture result in Table S9 describes structure that survives "
          f"binning rather than decomposing it, and it holds in two of the three sections. None of these limits touches the two claims we rely "
          f"on \u2014 that per-gene exponents are reproducible and abundance-ordered, and that no transcript "
          f"approaches the counting limit."),
    ("p", f"The replication results change what the paper can claim and what it cannot. Between sections of "
          f"one tumour, the pooled exponent moves by {num('abs(MU.contrasts.same_tissue_colon.delta)', 3)} "
          f"while the gap between tumours and brain is {num('MU.contrasts.colon_vs_brain.delta', 3)}; that is "
          f"a reproducibility statement about a one-section design, and it is the first one this literature has "
          f"for \u03b3. It does not extend to the per-gene summary: its median moves by "
          f"{num('MU.pergene.shift_overall', 3)} between the same two sections, and it moves in opposite "
          f"directions at the two ends of the abundance range. A reader who wants a single portable number "
          f"therefore has one, and it is the pooled value \u2014 with the caveat that what makes it portable is "
          f"also what makes it uninformative about most transcripts."),
    ("p", f"The Stereo-seq section tests the other axis. Its exponent is "
          f"{num('X.delta_gamma_pool', 3)} away from the Visium value for the same organ, so the number is not "
          f"a property of mouse brain. Its composition is: "
          f"{num('X.concentration.n_eff_over_genes', 2, 'pct')} of genes effective, "
          f"{num('X.stereo.top1_alpha_share', 0, 'pct')} of the exponent on one transcript ("
          f"{num('X.without_rrna.top1_share', 0, 'pct')} on Ttr once rRNA-derived loci are removed), and "
          f"per-gene exponents correlated at {num('X.replication.rho_raw', 2)} raw and "
          f"{num('X.replication.rho_partial_abundance', 2)} after abundance is removed. That is what a "
          f"composition statistic should do when the chemistry changes: the weights move and the structure "
          f"holds. It also settles a question the Visium data could not, because both platforms put the same "
          f"transcript at the top once the contaminant is gone."),
    ("p", "\u201cBinning buys counts, not precision\u201d remains true, but not because the tissue imposes a "
          "noise floor. The exponent belongs to transcripts, and a pooled value is a weighted mean of them, so "
          "what a section reports depends on where its heaviest transcripts sit relative to the window. In "
          f"colon the weight sits on IGKC, which is confined to clusters: its own exponent is "
          f"{num('J.colon.top1_gamma', 2)}, the pooled value lands "
          f"{num('D3.colon_pooled_minus_median', 2)} above the per-gene median, and coarsening averages the "
          "clusters away. In brain the weight sits on Camk2n1, which is spread across the section: its exponent "
          f"is {num('J.brain.top1_gamma', 2)}, the pooled value lands "
          f"{num('abs(D3.brain_pooled_minus_median)', 2)} below the median, and coarsening does not touch it. "
          "The same arithmetic is why the two tumours agree on the pooled statistic and disagree on the "
          "per-gene one. The measurement to report is the distribution, and the resource to buy is sections."),
]

METHODS = [
    ("h2", "Datasets"),
    ("p", f"Three public Visium HD datasets generated by 10x Genomics (spatial-exp 3.0.0): "
          f"Visium_HD_Human_Colon_Cancer, Visium_HD_Human_Lung_Cancer and Visium_HD_Mouse_Brain, one section "
          f"each. The colon tissue contributes a fourth section, Visium_HD_Human_Colon_Cancer_P1, released "
          f"separately; it enters only the replication analysis and is described under that heading. A fifth "
          f"section, Stereo-seq C04042E2 mouse brain, supplies the second platform. We used the binned outputs "
          f"at native 2 \u00b5m resolution and the in-tissue mask supplied with them. In-tissue bins number {num('d.colon.n_bins', 0, 'int')} (colon), "
          f"{num('d.lung.n_bins', 0, 'int')} (lung) and {num('d.brain.n_bins', 0, 'int')} (brain); genes number "
          f"{num('d.colon.n_genes_raw', 0, 'int')} before and {num('d.colon.n_genes_noMT', 0, 'int')} after "
          f"removal of mitochondrial genes in the two human samples, and {num('d.brain.n_genes_raw', 0, 'int')} "
          f"and {num('d.brain.n_genes_noMT', 0, 'int')} in mouse brain. Total UMI per section are "
          f"{num('d.colon.umi_total', 0, 'int')}, {num('d.lung.umi_total', 0, 'int')} and "
          f"{num('d.brain.umi_total', 0, 'int')} (Table 1). All analyses exclude mitochondrial genes unless "
          f"stated."),
    ("h2", "Exact re-aggregation and the analysis domain"),
    ("p", "Bins are aggregated on the native grid into non-overlapping square blocks of side d, aligned to the "
          "array origin, for d = 2, 4, 8, 16, 32, 64, 128 and 256 \u00b5m. Counts are summed as integers; "
          "in-tissue bins with no counts are retained as zeros rather than dropped, so the denominator at every "
          "scale is the number of in-tissue blocks, including blocks with no counts. Blocks are keyed on a grid "
          "with a stride of 4096 positions per row, which exceeds the largest array coordinate, so no key "
          "collision is possible. The implementation was checked against a brute-force groupby on the full "
          "matrix: sums are element-wise identical, including in masked regions with holes."),
    ("p", f"Two domain conventions are defined. The primary one counts every block that contains at least one "
          f"in-tissue 2 \u00b5m bin. The strict one counts only blocks whose 2 \u00b5m bins all lie inside the "
          f"tissue. At 256 \u00b5m the strict convention retains "
          f"{num('F.colon.by_scale[256].fill_frac', 1, 'pct')} of the colon grid, "
          f"{num('F.lung.by_scale[256].fill_frac', 1, 'pct')} of lung and "
          f"{num('F.brain.by_scale[256].fill_frac', 1, 'pct')} of brain; the difference between "
          f"conventions is reported in Table S5 and in the first Results subsection."),
    ("h2", "Pooled overdispersion"),
    ("p", "For each gene g at scale d we compute the mean \u03bc_g and variance Var_g across in-tissue blocks, "
          "then \u03b1_pool(d) = \u03a3_g [Var_g(d) \u2212 \u03bc_g(d)] / \u03a3_g \u03bc_g(d)\u00b2. With the "
          "n\u22121 denominator, E[Var_g] = \u03bc_g + \u03b1_g\u03bc_g\u00b2, so the ratio is an unbiased moment "
          "estimator of the \u03bc\u00b2-weighted mean overdispersion. Weighting by \u03bc\u00b2 is "
          "inverse-variance weighting under Poisson noise; no gene filter is applied and every gene in the "
          "matrix contributes. Exponents are ordinary least squares slopes of log \u03b1 on log d over the eight "
          "scales, reported as positive \u03b3 in \u03b1 \u221d d^\u2212\u03b3."),
    ("h2", "Within-block estimation"),
    ("p", "To remove variance that comes from differences between tissue domains, we partition each section into "
          "blocks of side 8d and compute the pooled overdispersion from within-block variance only. Because the "
          "block side scales with d, every block contains an 8\u00d78 grid of d-blocks at every scale, so the "
          "estimator sees the same geometry throughout. A fixed block size shrinks the variance as d grows and "
          "produced a spurious steepening to \u22121.01 in earlier work; the scale-invariant version does not."),
    ("h2", "Per-gene exponents, filtering and reliability"),
    ("p", f"For each gene we compute \u03b1_g(d) = [Var_g(d) \u2212 \u03bc_g(d)]/\u03bc_g(d)\u00b2 at the eight "
          f"scales and fit the same log\u2013log slope. A gene enters the distribution if \u03b1_g is finite and "
          f"positive at all eight scales and \u03bc_g at 256 \u00b5m is at least 0.5 counts per bin; this retains "
          f"{num('N.gene_level.colon.n_fit', 0, 'int')}, {num('N.gene_level.lung.n_fit', 0, 'int')} and "
          f"{num('N.gene_level.brain.n_fit', 0, 'int')} genes. The filter is on \u03bc at the coarsest scale, not "
          f"on the fitted quantity, and the controls in the text show that it does not create the abundance "
          f"gradient. Reliability is measured by cutting each section into two halves at the column nearest the "
          f"centroid that is a multiple of the coarsest block (256 \u00b5m), fitting each half independently and "
          f"correlating the per-gene exponents of genes that pass the filter in both halves."),
    ("h2", "Decomposition of the pooled exponent"),
    ("p", "The pooled estimator is rewritten gene by gene as \u03b1_pool(d) = \u03a3_g w_g(d)\u03b1_g(d) with "
          "w_g = \u03bc_g\u00b2/\u03a3\u03bc\u00b2, which holds by construction; the identity was checked against "
          "the direct estimator at all eight scales in all three sections. Gene g's share of the pooled "
          "numerator is v_g(d) = w_g(d)\u03b1_g(d)/\u03b1_pool(d), averaged over the eight scales and normalised "
          "over the genes that carry a computable exponent; the effective count under those weights is "
          "1/\u03a3v\u00b2, and it is reported alongside 1/\u03a3w\u00b2 rather than instead of it. Five summary "
          "scalars \u2014 pooled \u03b3, the per-gene median, the per-gene mean, the 5% trimmed mean and the "
          "\u03bc\u00b2-weighted median \u2014 were computed inside each of the four column bands of the colon "
          "section, using the same cuts and the same gene filter as the reliability analysis."),
    ("h2", "Spatial statistics"),
    ("p", "Per-gene spatial structure is measured on the 32 \u00b5m grid, where counts per bin are large enough "
          "for the statistics to be estimable at all. The semivariogram is computed at seven lags along the axis "
          "and the diagonal (32 to 256 \u00b5m) and fitted with an exponential model plus a nugget. The range is "
          "found by a profile search on a grid from 40 \u00b5m to 20 mm, so a value at the lower bound is a bound "
          "and not an estimate. Moran's I is computed on the symmetric eight-neighbour adjacency at the same "
          "resolution. Both statistics are related to the per-gene exponent by Spearman correlation, by partial "
          "correlation controlling log\u2081\u2080 mean counts, and within deciles of mean counts."),
    ("h2", "Expression layers and variance decomposition"),
    ("p", f"Genes are binned by mean counts per 32 \u00b5m bin into half-decade layers with common edges for all "
          f"three sections. Pooled \u03b1 is recomputed within each layer, and a two-way decomposition without "
          f"replication (6 layers \u00d7 3 sections) partitions the variance of the layer-level exponents into "
          f"layer, section and residual components. The same decomposition is repeated on the per-gene median. "
          f"The decomposition uses the six layers in which all three sections retain at least twenty genes with "
          f"a computable exponent. Those layers cover {num('LW.colon.coverage_6layers', 1, 'pct')}, "
          f"{num('LW.lung.coverage_6layers', 1, 'pct')} and {num('LW.brain.coverage_6layers', 1, 'pct')} of "
          f"\u03a3\u03bc\u00b2 in colon, lung and brain, so the decomposition describes the per-gene population "
          f"and does not decompose the pooled value. Table 3 adds a seventh row for "
          f"\u03bc\u2083\u2082 \u2265 10 counts per bin, which carries the remaining "
          f"{num('L.by_tissue.colon[8].share', 1, 'pct')}, {num('L.by_tissue.lung[8].share', 1, 'pct')} and "
          f"{num('L.by_tissue.brain[8].share', 1, 'pct')} of the weight but too few genes for a median."),
    ("h2", "Model comparison"),
    ("p", "Four models are fitted to the eight pooled \u03b1 values: a pure power law (2 parameters), a quadratic "
          "in log d (3), an exponential correlation model \u03b1 = A\u00b7H(d/\u03be) (2) and a Gaussian "
          "correlation model (2). In the two correlation models H is the exact double integral of the "
          "correlation function over a square window, normalised so that H(0) = 1 and H(w) \u2192 2\u03c0/w\u00b2. "
          "Models are compared by AIC. The "
          "quadratic term is tested by an F test on one and five degrees of freedom. The correlation length is "
          "given by a profile likelihood over a grid from 1 to 10\u2076 \u00b5m, with a 95% threshold from the F "
          "distribution. The same per-gene comparison uses the two 2-parameter shapes and a grid from 0.5 to "
          "10\u2076 \u00b5m."),
    ("h2", "Synthetic abundance control"),
    ("p", f"Counts were simulated on a Gaussian correlated field with a power-law spectrum of exponent "
          f"{num('S.expression.beta_true', 2)} at "
          f"{num('S.expression.n_genes', 0, 'int')} genes spanning four orders of magnitude in mean expression, "
          f"with pure Poisson counting and the same eight-scale re-aggregation and fitting pipeline. If the "
          f"estimator manufactured an abundance dependence, the recovered exponent would drift across "
          f"expression layers; it does not (Fig. S2). The same simulation calibrates the absolute level: it "
          f"returns {num('S.expression.gamma_all.med', 2)} for an input of "
          f"{num('S.expression.beta_true', 2)}, because a single power law fitted inside a still-curving window "
          f"is steeper than the asymptotic slope."),
    ("h2", "Depth, fill and regional controls"),
    ("p", f"Depth is controlled by binomial thinning of the native 2 \u00b5m counts to "
          f"{num('N.depth.runs[PP(0.46)].p', 0, 'pct')}, {num('N.depth.runs[PP(0.25)].p', 0, 'pct')} and "
          f"{num('N.depth.runs[PP(0.1)].p', 0, 'pct')} of the original UMIs, followed by the full re-aggregation and "
          f"fitting pipeline. Fill is controlled by recomputing \u03b1 from blocks that lie entirely inside the "
          f"tissue. Regional variation is measured by cutting each section into four column bands on the "
          f"coarsest block grid and by cutting it into two halves; cuts are placed on multiples of 256 \u00b5m so "
          f"that no coarse block is divided. A block bootstrap over superblocks of 512, 1024 and 2048 \u00b5m "
          f"(200 resamples each) is reported as a function of block size, because resampling units smaller than "
          f"the correlation length treats correlated blocks as independent."),
    ("h2", "Uncertainty and minimum detectable difference"),
    ("p", f"Two error models are used and both are reported. The first is the spread of four regional bands "
          f"inside a section, which is the conservative one: it treats the differences between regions as "
          f"between-section variability, even though they share donor, block and library. The second is a "
          f"spatial block bootstrap used for the section-replication analysis below: the section is cut into "
          f"contiguous superblocks of 512 \u00b5m, blocks are resampled with replacement, and \u03b3 is refitted "
          f"inside each of the {num('MU.n_boot', 0, 'int')} replicates. The blocks are an order of magnitude "
          f"larger than the fitted correlation length (\u03be = "
          f"{num('M.colon.xi_hat', 2)}\u2013"
          f"{num('max(M.colon.xi_hat,M.lung.xi_hat,M.brain.xi_hat)', 1)} \u00b5m), so a block carries the "
          f"within-block spatial correlation, which an i.i.d. tile bootstrap of the same tile size does not. "
          f"The block bootstrap measures the sampling error of one section; it does not measure the "
          f"difference between sections, which is why the two colon sections are used for that. Minimum "
          f"detectable differences at 80% power and a two-sided 5% level are MDE = exp(k\u00b7SD/\u221an) with "
          f"k = 2.802 and n sections per group."),
    ("h2", "Second section and section-level replication"),
    ("p", f"Visium_HD_Human_Colon_Cancer_P1 was downloaded, filtered and aggregated by exactly the scripts "
          f"used for the first colon section, with no change to the gene filter, the domain convention or the "
          f"estimator. Its per-scale moments were reduced with the same function, and the four sections were "
          f"compared on \u03b3_pool, the per-gene median and the effective gene count. Contrasts are "
          f"differences of section estimates with a standard error of "
          f"\u221a(SE\u00b2_a + SE\u00b2_b); the SE is taken from the block bootstrap in the primary analysis "
          f"and from the regional-band SD in the conservative one, and both are tabulated (Table S11). The "
          f"minimum detectable difference for a two-section design is ({num('MU.mde.z_alpha', 2)} + "
          f"{num('MU.mde.z_power', 2)})\u00b7\u221a2\u00b7SE, with "
          f"SE the mean of the two colon sections. An intraclass correlation over the four bands was computed "
          f"and is reported in Table S11, but is not used: the bands are non-overlapping spatial subregions, "
          f"not exchangeable replicates, so the random-effects model that defines ICC does not hold."),
    ("h2", "Per-gene comparison between sections"),
    ("p", f"Per-gene exponents were recomputed for both colon sections with the filter used throughout "
          f"(finite and positive at all eight scales, mean at 256 \u00b5m at least 0.5). Abundance deciles "
          f"were defined on the first section's mean counts at 256 \u00b5m and applied unchanged to the second, "
          f"so that a decile holds genes of the same abundance in both. The median exponent was computed in "
          f"each decile for each section, and the shift was bootstrapped by resampling genes within the "
          f"decile ({num('MU.n_boot', 0, 'int')} replicates). Deciles are used rather than a single median "
          f"because a median over all genes mixes a change in the exponent with a change in which genes are "
          f"above the filter, and the two have to be separated before either can be interpreted."),
    ("h2", "Stereo-seq section and cross-platform comparison"),
    ("p", f"C04042E2 is a Stereo-seq mouse whole-brain section distributed through the STOmicsDB demo "
          f"portal, captured on a 0.5 \u00b5m DNA nanoball array. The release stores expression as flat records of "
          f"(x int32, y int32, count uint8) per bin size, 9 bytes each, with offsets in a companion gene "
          f"table. We read the eight bin sizes that correspond to 2, 4, \u2026, 256 \u00b5m blocks "
          f"({num('max(X.stereo.n_cells)', 0, 'int')} records at 2 \u00b5m), accumulated first and second moments "
          f"per gene and per scale in one streaming pass, and then applied the Visium filter and estimator "
          f"unchanged. No gene was removed except in the contamination attack, which drops the rRNA-derived "
          f"loci Gm42418 and Gm26917."),
    ("p", f"Criteria were frozen before the Stereo-seq moments existed, and are listed with their outcomes in "
          f"Table S10. They are scale-free on purpose, because an absolute exponent is a property of the "
          f"platform. The concentration criteria are (H1a) the effective number of genes under the "
          f"\u03b1-contribution weight stays below 1% of the genes tested and (H1b) one transcript carries "
          f"more than 10% of the pooled exponent. The ordering criterion is (H2) the pooled value lies below "
          f"the per-gene median, as on Visium HD brain. The replication criterion is (H4) per-gene exponents "
          f"correlate between platforms at \u03c1 > {num('X.replication.h4_threshold', 2)} over genes "
          f"measured on both. H4 is tested twice, "
          f"once raw and once after each platform's exponent is rank-residualised on its own log abundance, "
          f"because both platforms rank genes by abundance and that alone can produce a correlation. An "
          f"earlier version of H1, which required the top transcript to exceed 40% of the exponent, was "
          f"discarded before the data were seen to fail or pass it, and its outcome is reported alongside."),
    ("h2", "Software"),
    ("p", "Python 3.13 with numpy 2.5, scipy 1.18 and matplotlib 3.11. Every number in the manuscript, including "
          "all supplementary tables, is generated from a locked JSON number table by scripts deposited with the "
          "code; the text contains no number that is not an injection from that table."),
]

LEGENDS = [
    ("Fig. 1", "Noise falls one order of magnitude, not four. (a) Exact integer re-aggregation: native 2 \u00b5m "
               "bins are summed into 4, 8, 16, 32, 64, 128 and 256 \u00b5m blocks on the same grid, with no "
               "interpolation and no dropped reads. (b) Pooled overdispersion normalised to its value at 2 \u00b5m, "
               "against bin size, for colon (dark blue), lung (mid blue) and mouse brain (ochre). Grey dashed: "
               "\u03b1 \u221d d^\u22122, the counting-statistics limit that merging bins is meant to approach; "
               "grey dotted: \u03b1 \u221d d\u2070. If the field's working assumption held, the data would lie on "
               "the steep grey line."),
    ("Fig. 2", "No single correlation length. (a) Local slope of log \u03b1 against log d within each segment, "
               "plotted at the geometric mean of the segment endpoints; colon (dark blue), lung (mid blue), mouse "
               "brain (ochre). Solid: observed. Dashed: prediction of the "
               "quadratic model. Dotted: prediction of the exponential correlation model with the fitted \u03be. "
               f"The finite-\u03be prediction climbs to {num('M.lung.local_slope.predicted_xi[-1]', 2)} while the "
               f"observed slope never exceeds "
               f"{num('max(max(M.colon.local_slope.observed),max(M.lung.local_slope.observed),max(M.brain.local_slope.observed))', 2)}. "
               "(b) \u0394AIC of four models against the best model, by section, in the same colours. (c) Profile "
               "\u0394RSS for the "
               "correlation length, with the 95% threshold (dashed) and the profile interval (bar) for each "
               "section; the three intervals overlap."),
    ("Fig. 3", "The exponent belongs to transcripts. (a) Per-gene exponent against mean counts per 32 \u00b5m bin "
               "for the three sections; lines are medians within abundance bins. (b) Distribution of per-gene "
               "exponents; white dot, median; bar, interquartile range; dashed, pooled value for that section. "
               "(c) Per-gene exponent in the two halves of the colon section; r, Pearson correlation. (d) Median "
               "exponent by abundance decade for the three sections (solid) and for the synthetic control with a "
               "known input exponent of 0.55 (dashed); the control is flat, the data are not."),
    ("Fig. 4", "The pooled exponent is carried by a few transcripts. (a) Cumulative share of \u03a3\u03bc\u00b2 "
               "against genes ranked by \u03bc\u00b2, for colon (dark blue), lung (mid blue) and mouse brain "
               "(ochre); dots mark the number of genes needed for half of the "
               "weight. (b) Exponent recomputed from the top-k genes only (solid circles) against the random-k "
               "95% range (shaded), in the same colours, with the full-gene value as a dashed line; the first "
               "point of each curve is the heaviest transcript alone, which sits above the full-set value in "
               "colon and lung and below it in brain. (c) Largest "
               "change in exponent "
               "produced by leaving one gene out; the gene responsible is named below each bar, in the order "
               "colon, lung, brain, and the note at the top gives the number of genes that move the exponent by "
               "more than 0.05."),
    ("Fig. 5", "One section does not pin the exponent down. (a) Exponent recomputed from four column bands and "
               "two halves of the colon section; the dashed line is the whole-section value. (b) Difference in "
               "exponent between sections (bars) against the 95% half-width implied by the regional SD "
               "(whiskers); the dotted line is the regional range itself. (c) Exponent after binomial thinning of "
               "the colon section; the vertical axis is magnified to show the shift. (d) SD of the exponent across "
               "200 block-bootstrap "
               "resamples, as a function of superblock side; the dashed line is the SD across the four bands."),
    ("Fig. 6", "What each design choice can buy. (a) Coefficient of variation of total counts per bin (solid) and "
               "the Poisson prediction for the same counts (dashed), against bin size, for colon (dark blue), lung "
               "(mid blue) and mouse brain (ochre). (b) Minimum detectable "
               "difference in \u03b3 at 80% power against the number of sections per group, computed from the "
               "regional SD of each section in the same colours; dotted lines mark the observed colon\u2013lung "
               "and colon\u2013brain gaps. (c) Change in exponent produced by one manipulation, on a "
               "logarithmic axis. Depth is three orders of magnitude below the others; no manipulation "
               "approaches the gap between the measured exponent and the counting limit."),
    ("Fig. 7", "The composition survives a change of capture chemistry. (a) Pooled overdispersion normalised to "
               "its value at 2 \u00b5m against bin size: colon (dark blue) and mouse brain (ochre) on Visium HD, "
               "mouse brain on Stereo-seq (purple, dashed). The Stereo-seq ladder is built from 0.5 \u00b5m "
               "nanoballs rather than 2 \u00b5m bins, and its absolute level differs from the Visium section of "
               "the same tissue; grey dotted, \u03b1 \u221d d\u207b\u00b2. (b) Cumulative share of the pooled "
               "exponent against genes ranked by their contribution, for mouse brain on both platforms. Both "
               f"curves leave the origin at {num('X.stereo.top1_alpha_share', 0, 'pct')} against "
               f"{num('X.visium.top1_alpha_share', 0, 'pct')} and reach half the exponent within the top few "
               "dozen genes. (c) Per-gene exponent on Stereo-seq against the same gene on Visium HD, over "
               f"{num('X.replication.n', 0, 'int')} genes measured on both; colour is the number of genes per "
               "cell on a logarithmic scale, and the diagonal is shown for reference. \u03c1 is the rank "
               f"correlation, \u03c1|abundance the same correlation after each platform's exponent is "
               "rank-residualised on its own log abundance."),
    ("Fig. 8", "The pooled exponent replicates between sections; the per-gene summary does not. (a) Pooled "
               "exponent of four sections with 95% intervals from a 512 \u00b5m superblock bootstrap "
               f"(B = {num('MU.n_boot', 0, 'int')}); the two colon sections are the same colour and different "
               "markers. The difference between them is annotated with the difference between colon and brain "
               "and the minimum detectable difference of a two-section design. (b) The eight points behind each "
               "estimate, with the fitted value in the legend. (c) Median per-gene exponent in abundance "
               "deciles, first colon section against second, with the deciles defined on the first. The two "
               "distributions cross: the second section is higher at low abundance and lower at high abundance."),

    ("Fig. S1", "The exponent depends on the fit window. \u03b3 fitted from a varying first bin size to a fixed "
                "last bin size of 256 \u00b5m."),
    ("Fig. S2", "The synthetic control does not manufacture an abundance dependence. Recovered per-gene exponent "
                "by expression layer for a simulated field with an input exponent of 0.55 and Poisson counting; "
                "points are medians and bars are interquartile ranges."),
    ("Fig. S3", "Region-to-region variation in \u03b1(d). Thick lines: four column bands of the colon section "
                "fitted separately. Grey lines: the two halves."),
    ("Fig. S4", "Global and within-block estimates of \u03b1(d), by section. Solid: global; dashed: within-block, "
                "which removes between-domain variance. The within-block estimator lowers the exponent in all "
                "three sections and by the most in brain."),
    ("Fig. S5", "The transcript that dominates \u03a3\u03bc\u00b2 is spatially localised in colon and broadly "
                "expressed in brain. Counts per 32 \u00b5m bin for IGKC and for all genes in the colon section, and "
                "for Camk2n1 and all genes in the mouse brain section. Colour scale is clipped at the 99th "
                "percentile."),
    ("Fig. S6", "(a) Effective number of genes behind the pooled exponent, computed with the variance weights "
                "(\u03bc\u00b2, pale) and with the weights the exponent itself carries (\u03b1 contribution, "
                "solid); the figure above each pair is the residual left when the pooled exponent is rebuilt as "
                "the weighted mean of the per-gene exponents. (b) Five candidate summary scalars computed inside "
                "each of the four column bands of Fig. 5a (dots), with the mean across bands (rule) and the "
                "standard deviation across bands (number at right)."),
]

DECLARATIONS = [
    ("Funding", "This work received no external funding."),
    ("Competing interests", "The authors declare no competing interests."),
    ("Ethics approval and consent to participate",
     "Not applicable. This study reanalyses publicly available, de-identified data generated by 10x Genomics; no "
     "human participants or animals were involved and no new data were collected."),
    ("Consent for publication", "Not applicable."),
    ("Data availability",
     "Four sections are public and were downloaded from 10x Genomics (spatial-exp 3.0.0): "
     "Visium_HD_Human_Colon_Cancer and Visium_HD_Human_Colon_Cancer_P1, Visium_HD_Human_Lung_Cancer and "
     "Visium_HD_Mouse_Brain (https://www.10xgenomics.com/datasets). The Stereo-seq mouse brain section C04042E2 "
     "is public through the STOmicsDB demonstration portal. Every number in the manuscript, including all "
     "supplementary tables, is generated from a locked JSON number table that ships with the analysis code in "
     "the submission package."),
    ("Code availability",
     "Analysis code (Python 3.13, numpy 2.5, scipy 1.18, matplotlib 3.11), the locked number table and the "
     "scripts that regenerate every figure, table and number in the text are provided in the submission package "
     "(05_Code). The locked table alone regenerates all figures and documents and re-runs the three consistency "
     "checks; the table itself is rebuilt from the public inputs by the upstream scripts in the same directory. "
     "The reporting procedure itself is released as a standalone package, "
     "05_Code/scaling_report-1.0.0 (MIT licence, numpy only), which takes a table of per-gene means and "
     "variances at a bin ladder and returns the pooled exponent, the per-gene distribution, the number of "
     "genes behind the reported value and the range of per-gene medians that value admits; it reproduces "
     "every entry of Table 3 from the moments tables and requires no raw matrix."),
    ("Authors' contributions",
     "C.T. designed the study, wrote the analysis code, ran all analyses and wrote the manuscript. F.W. and Y.Z. "
     "verified the statistical pipeline and the number table. L.Z., H.T., L.L. and M.W. contributed to "
     "interpretation and critically revised the manuscript. All authors approved the final version."),
    ("Acknowledgements", "We thank the 10x Genomics public datasets programme for making the three sections available."),
    ("Use of AI",
     "Language editing was assisted by ChatGPT (OpenAI). All analyses, numerical results and interpretations were "
     "produced and verified by the authors; no AI system generated data or conclusions."),
]


# ----------------------------------------------------------------- docx 工具
def set_base_style(doc, size=11, spacing=1.5):
    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(size)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), "宋体")
    pf = st.paragraph_format
    pf.line_spacing = spacing
    pf.space_after = Pt(6)
    pf.space_before = Pt(0)


def add_page_numbers(doc):
    for section in doc.sections:
        p = section.footer.paragraphs[0]
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        run = p.add_run()
        for el, attr in [("w:fldChar", {"w:fldCharType": "begin"}),
                         ("w:instrText", None),
                         ("w:fldChar", {"w:fldCharType": "end"})]:
            e = OxmlElement(el)
            if attr:
                for k, v in attr.items():
                    e.set(qn(k), v)
            else:
                e.set(qn("xml:space"), "preserve")
                e.text = " PAGE "
            run._r.append(e)


def para(doc, text, bold=False, italic=False, size=None, align="justify", space_after=6, spacing=1.5):
    p = doc.add_paragraph()
    p.alignment = {"justify": WD_ALIGN_PARAGRAPH.JUSTIFY, "left": WD_ALIGN_PARAGRAPH.LEFT,
                   "center": WD_ALIGN_PARAGRAPH.CENTER}[align]
    p.paragraph_format.space_after = Pt(space_after)
    p.paragraph_format.line_spacing = spacing
    r = p.add_run(text)
    r.bold = bold
    r.italic = italic
    if size:
        r.font.size = Pt(size)
    return p


def heading(doc, text, size=12, space_before=12):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(space_before)
    p.paragraph_format.space_after = Pt(4)
    r = p.add_run(text)
    r.bold = True
    r.font.size = Pt(size)
    return p


def shade(cell, hexcolor):
    tcPr = cell._tc.get_or_add_tcPr()
    sh = OxmlElement("w:shd")
    sh.set(qn("w:val"), "clear")
    sh.set(qn("w:fill"), hexcolor)
    tcPr.append(sh)


def make_table(doc, rows, widths=None, font=7.5, first_col_bold=False, autofit=False):
    t = doc.add_table(rows=len(rows), cols=len(rows[0]))
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = autofit
    for i, row in enumerate(rows):
        for j, val in enumerate(row):
            cell = t.cell(i, j)
            cell.text = ""
            p = cell.paragraphs[0]
            p.paragraph_format.space_after = Pt(0.8)
            p.paragraph_format.space_before = Pt(0.8)
            p.paragraph_format.line_spacing = 1.0
            r = p.add_run(str(val))
            r.font.size = Pt(font)
            r.font.name = "Times New Roman"
            if i == 0:
                r.bold = True
                shade(cell, "E8EEF4")
            if first_col_bold and j == 0 and i > 0:
                r.bold = True
    if widths:
        for j, w in enumerate(widths):
            for row in t.rows:
                row.cells[j].width = Cm(w)
    return t


def table1_rows():
    d = N.datasets
    rows = [
        ["Metric", "Human colon cancer", "Human lung cancer", "Mouse brain"],
        ["10x dataset", "Visium_HD_Human_Colon_Cancer", "Visium_HD_Human_Lung_Cancer", "Visium_HD_Mouse_Brain"],
        ["Species / tissue state", "Human / cancer", "Human / cancer", "Mouse / normal"],
        ["In-tissue 2 \u00b5m bins", num("d.colon.n_bins", 0, "int"), num("d.lung.n_bins", 0, "int"),
         num("d.brain.n_bins", 0, "int")],
        ["Blocks at 256 \u00b5m (inside / full)", 
         f"{num('d.colon.power[256].n_bins', 0, 'int')} / {num('F.colon.by_scale[256].n_domain_full', 0, 'int')}",
         f"{num('d.lung.power[256].n_bins', 0, 'int')} / {num('F.lung.by_scale[256].n_domain_full', 0, 'int')}",
         f"{num('d.brain.power[256].n_bins', 0, 'int')} / {num('F.brain.by_scale[256].n_domain_full', 0, 'int')}"],
        ["Zero-count in-tissue bins", num("d.colon.zero_frac", 1, "pct"), num("d.lung.zero_frac", 1, "pct"),
         num("d.brain.zero_frac", 1, "pct")],
        ["Genes (all / after MT removal)",
         f"{num('d.colon.n_genes_raw', 0, 'int')} / {num('d.colon.n_genes_noMT', 0, 'int')}",
         f"{num('d.lung.n_genes_raw', 0, 'int')} / {num('d.lung.n_genes_noMT', 0, 'int')}",
         f"{num('d.brain.n_genes_raw', 0, 'int')} / {num('d.brain.n_genes_noMT', 0, 'int')}"],
        ["Total UMI (all genes)", num("d.colon.umi_total", 0, "int"), num("d.lung.umi_total", 0, "int"),
         num("d.brain.umi_total", 0, "int")],
        ["UMI per 2 \u00b5m bin, mean (median)",
         f"{num('d.colon.umi_per_bin', 1)} ({num('d.colon.umi_median', 0, 'int')})",
         f"{num('d.lung.umi_per_bin', 1)} ({num('d.lung.umi_median', 0, 'int')})",
         f"{num('d.brain.umi_per_bin', 1)} ({num('d.brain.umi_median', 0, 'int')})"],
        ["\u03b1 at 2 \u00b5m \u2192 \u03b1 at 256 \u00b5m",
         f"{num('d.colon.alpha[2]', 1)} \u2192 {num('d.colon.alpha[256]', 2)}",
         f"{num('d.lung.alpha[2]', 1)} \u2192 {num('d.lung.alpha[256]', 2)}",
         f"{num('d.brain.alpha[2]', 2)} \u2192 {num('d.brain.alpha[256]', 2)}"],
        ["\u03b3, pooled 2\u2013256 \u00b5m", num("J.colon.gamma_full", 3), num("J.lung.gamma_full", 3),
         num("J.brain.gamma_full", 3)],
        ["\u03b3, within block (8d)", num("abs(d.colon.gamma_local)", 3), num("abs(d.lung.gamma_local)", 3),
         num("abs(d.brain.gamma_local)", 3)],
        ["CV of counts, 2 \u2192 256 \u00b5m",
         f"{num('d.colon.cv[2]', 0, 'pct')} \u2192 {num('d.colon.cv[256]', 0, 'pct')}",
         f"{num('d.lung.cv[2]', 0, 'pct')} \u2192 {num('d.lung.cv[256]', 0, 'pct')}",
         f"{num('d.brain.cv[2]', 0, 'pct')} \u2192 {num('d.brain.cv[256]', 0, 'pct')}"],
    ]
    return rows


def table2_rows():
    rows = [["Section", "\u03b3 pooled", "Regional \u03b3 (2\u00d72 quadrants)", "SD across regions",
             "Per-gene \u03b3, median (IQR)", "Effective genes", "Heaviest transcript (share)",
             "\u03b3 without it"]]
    for t, lab in (("colon", "Colon"), ("lung", "Lung"), ("brain", "Brain")):
        gl = N.gene_level[t]["gamma_gene"]
        jk = N.jackknife[t]
        rows.append([
            lab,
            num(f"J.{t}.gamma_full", 3),
            f"{num(f'abs(D3.{t}_quad_min)', 2)}\u2013{num(f'abs(D3.{t}_quad_max)', 2)}",
            num(f"D3.{t}_quad_sd", 3),
            f"{num(f'N.gene_level.{t}.gamma_gene.median', 2)} "
            f"({num(f'N.gene_level.{t}.gamma_gene.q25', 2)}\u2013{num(f'N.gene_level.{t}.gamma_gene.q75', 2)})",
            num(f"N.gene_level.{t}.composition.n_eff", 1),
            f"{jk['top1_gene']} ({num(f'J.{t}.top1_weight', 1, 'pct')})",
            num(f"J.{t}.gamma_full + J.{t}.loo.worst_delta", 3),
        ])
    return rows


def table3_rows():
    rows = [["\u03bc\u2083\u2082 (counts/bin)", "Genes fitted C / L / B", "\u03b3 median colon",
             "\u03b3 median lung", "\u03b3 median brain", "\u03a3\u03bc\u00b2 colon", "\u03a3\u03bc\u00b2 lung",
             "\u03a3\u03bc\u00b2 brain"]]
    for i in range(6):
        rows.append([
            f"{num(f'LW.colon.rows[{i}].lo', 2, 'g')}\u2013{num(f'LW.colon.rows[{i}].hi', 2, 'g')}",
            f"{num(f'L.by_tissue.colon[{i}].n_gamma', 0, 'int')} / "
            f"{num(f'L.by_tissue.lung[{i}].n_gamma', 0, 'int')} / "
            f"{num(f'L.by_tissue.brain[{i}].n_gamma', 0, 'int')}",
            num(f"L.by_tissue.colon[{i}].gamma_median", 2),
            num(f"L.by_tissue.lung[{i}].gamma_median", 2),
            num(f"L.by_tissue.brain[{i}].gamma_median", 2),
            num(f"L.by_tissue.colon[{i}].share", 3, 'pct'),
            num(f"L.by_tissue.lung[{i}].share", 3, 'pct'),
            num(f"L.by_tissue.brain[{i}].share", 3, 'pct'),
        ])
    rows.append([
        "\u226510",
        f"{num('L.by_tissue.colon[8].n_gamma', 0, 'int')} / "
        f"{num('L.by_tissue.lung[8].n_gamma', 0, 'int')} / "
        f"{num('L.by_tissue.brain[8].n_gamma', 0, 'int')}",
        "\u2014", "\u2014", "\u2014",
        num("L.by_tissue.colon[8].share", 1, 'pct'),
        num("L.by_tissue.lung[8].share", 1, 'pct'),
        num("L.by_tissue.brain[8].share", 1, 'pct'),
    ])
    return rows


def _wc(s):
    return len(re.findall(r"[A-Za-z][A-Za-z'\u2019\-]*", s))


def build_manuscript(outdir):
    doc = Document()
    set_base_style(doc)
    sec = doc.sections[0]
    sec.page_height, sec.page_width = Cm(29.7), Cm(21.0)
    sec.top_margin = sec.bottom_margin = Cm(2.5)
    sec.left_margin = sec.right_margin = Cm(2.5)
    add_page_numbers(doc)

    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(60)
    r = p.add_run("Research article")
    r.italic = True
    r.font.size = Pt(11)
    para(doc, TITLE, bold=True, size=16, align="left", space_after=14, spacing=1.2)
    para(doc, f"Running title: {RUNNING}", italic=True, size=10, align="left", space_after=16)
    para(doc, AUTHORS, size=11, align="left", space_after=2)
    para(doc, "* These authors contributed equally to this work.", size=9, align="left", space_after=0)
    para(doc, "** Corresponding author.", size=9, align="left", space_after=10)
    para(doc, "Affiliation", bold=True, size=10, align="left", space_after=2)
    para(doc, AFFIL, size=10, align="left", space_after=10)
    para(doc, "Correspondence", bold=True, size=10, align="left", space_after=2)
    para(doc, CORR, size=10, align="left", space_after=10)
    para(doc, "Keywords", bold=True, size=10, align="left", space_after=2)
    para(doc, KEYWORDS, size=10, align="left", space_after=0)

    doc.add_page_break()
    heading(doc, "Abstract", size=13, space_before=0)
    for k, v in ABSTRACT.items():
        p = doc.add_paragraph()
        p.paragraph_format.space_after = Pt(6)
        p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        rr = p.add_run(f"{k}: ")
        rr.bold = True
        p.add_run(R(v))

    doc.add_page_break()
    heading(doc, "Background", size=13, space_before=0)
    for _, t in BACKGROUND:
        para(doc, R(t))

    heading(doc, "Results", size=13)
    for kind, t in RESULTS:
        if kind == "h2":
            heading(doc, R(t), size=11, space_before=10)
        else:
            para(doc, R(t))

    heading(doc, "Discussion", size=13)
    for _, t in DISCUSSION:
        para(doc, R(t))

    heading(doc, "Methods", size=13)
    for kind, t in METHODS:
        if kind == "h2":
            heading(doc, R(t), size=11, space_before=10)
        else:
            para(doc, R(t))

    doc.add_page_break()
    heading(doc, "Figure legends", size=13, space_before=0)
    for tag, t in LEGENDS:
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        p.paragraph_format.space_after = Pt(8)
        rr = p.add_run(f"{tag}. ")
        rr.bold = True
        p.add_run(R(t))

    heading(doc, "Declarations", size=13)
    for k, v in DECLARATIONS:
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        p.paragraph_format.space_after = Pt(6)
        rr = p.add_run(f"{k}: ")
        rr.bold = True
        p.add_run(R(v))

    doc.add_page_break()
    heading(doc, "References", size=13, space_before=0)
    for i, key in enumerate(CITE_ORDER, 1):
        para(doc, f"{i}. {fmt_ref(REFS[key])}", size=9, align="left", space_after=3, spacing=1.0)

    doc.add_page_break()
    heading(doc, "Table 1. Datasets and basic statistics", size=12, space_before=0)
    para(doc, "All values are computed from the native 2 \u00b5m matrices. \u03b1, overdispersion; \u03b3, "
              "effective exponent in \u03b1 \u221d d^\u2212\u03b3 over 2\u2013256 \u00b5m; CV, coefficient of "
              "variation of total counts per bin; MT, mitochondrial genes. Blocks at 256 \u00b5m are given for "
              "the two domain conventions: containing at least one in-tissue bin, and lying entirely inside the "
              "tissue.",
         size=8, align="left", space_after=6, spacing=1.0)
    make_table(doc, table1_rows(), widths=[4.1, 4.5, 4.5, 4.3], font=7, first_col_bold=True)

    heading(doc, "Table 2. What to report for a binned sample", size=12, space_before=14)
    para(doc, "Pooled \u03b3 is the value the field reports; the remaining columns are the quantities that make "
              "it interpretable. Regional \u03b3 is the range across four spatial quadrants and SD is their "
              "standard deviation; the four column bands used in Fig. 5a, cut on the coarse grid, give a range of "
              f"{num('N.splithalf.quarter_range', 3)} and SD {num('N.splithalf.quarter_sd', 3)} in colon. "
              "Effective genes is 1/\u03a3w\u00b2 with w = \u03bc\u2083\u2082\u00b2/\u03a3"
              "\u03bc\u2083\u2082\u00b2, the count implied by the variance weights; the same quantity computed "
              f"with the weights the exponent itself carries, 1/\u03a3v\u00b2, is "
              f"{num('N.decomp.by_tissue.colon.n_eff_alpha', 1)}, "
              f"{num('N.decomp.by_tissue.lung.n_eff_alpha', 1)} and "
              f"{num('N.decomp.by_tissue.brain.n_eff_alpha', 1)} (Table S7). The last column is the pooled "
              "\u03b3 recomputed with the heaviest transcript removed.",
         size=8, align="left", space_after=6, spacing=1.0)
    make_table(doc, table2_rows(), widths=[1.9, 1.6, 3.4, 2.0, 3.0, 1.8, 2.6, 1.6], font=7, first_col_bold=True)

    heading(doc, "Table 3. The exponent by expression layer", size=12, space_before=14)
    para(doc, "Layers are defined on mean counts per 32 \u00b5m bin with common edges for all three sections. "
              "The \u03b3 columns are within-layer medians over genes with a computable exponent, and the "
              "second column counts those genes (colon / lung / brain); the layer itself can hold more, as in "
              f"the lowest one, where {num('L.by_tissue.colon[0].n_gamma', 0, 'int')} of "
              f"{num('L.by_tissue.colon[0].n', 0, 'int')} colon genes can be fitted. Columns 6\u20138 give each "
              "layer's share of \u03a3\u03bc\u00b2, which is the weight that layer carries in the pooled value. "
              f"The six layers cover {num('D3.colon_coverage6', 1, 'pct')}, {num('D3.lung_coverage6', 1, 'pct')} "
              f"and {num('D3.brain_coverage6', 1, 'pct')} of \u03a3\u03bc\u00b2; the \u226510 row collects the "
              f"remaining {num('L.by_tissue.colon[8].share', 1, 'pct')}, "
              f"{num('L.by_tissue.lung[8].share', 1, 'pct')} and "
              f"{num('L.by_tissue.brain[8].share', 1, 'pct')} in "
              f"{num('L.by_tissue.colon[8].n_gamma', 0, 'int')}, "
              f"{num('L.by_tissue.lung[8].n_gamma', 0, 'int')} and "
              f"{num('L.by_tissue.brain[8].n_gamma', 0, 'int')} genes, too few for a median to mean anything. "
              f"In colon IGKC alone carries {num('J.colon.top1_weight', 1, 'pct')} of \u03a3\u03bc\u00b2, more "
              "than the six layers combined. Genes below 0.01 counts per bin are omitted; they carry a "
              "negligible share.",
         size=8, align="left", space_after=6, spacing=1.0)
    make_table(doc, table3_rows(), widths=[3.0, 3.2, 2.2, 2.2, 2.2, 2.0, 2.0, 2.0], font=7, first_col_bold=True)

    out = Path(outdir) / "01_Manuscript.docx"
    doc.save(out)
    return out


def fmt_ref(r):
    au = r.get("authors") or []
    names = []
    for a in au[:6]:
        fam = a.get("family", "")
        giv = a.get("given", "")
        ini = "".join([p[0] for p in giv.replace("-", " ").split() if p])
        names.append(f"{fam} {ini}".strip())
    if len(au) > 6:
        names.append("et al.")
    auth = ", ".join(names).rstrip(".")
    vol = r.get("volume") or ""
    page = r.get("page") or ""
    title = (r["title"] or "").strip().rstrip(".").replace("&amp;", "&")
    journal = (r["journal"] or "").replace("&amp;", "&")
    sep = "" if title.endswith(("?", "!")) else "."
    bits = [f"{auth}. {title}{sep} {journal}"]
    if r.get("year"):
        bits[-1] += f" {r['year']}"
    if vol:
        bits[-1] += f";{vol}"
    if page:
        bits[-1] += f":{page}"
    bits[-1] += f". doi:{r['doi']}"
    return "".join(bits)


def build_cover(outdir):
    doc = Document()
    set_base_style(doc, size=11, spacing=1.15)
    sec = doc.sections[0]
    sec.page_height, sec.page_width = Cm(29.7), Cm(21.0)
    sec.top_margin = sec.bottom_margin = Cm(2.2)
    sec.left_margin = sec.right_margin = Cm(2.4)

    para(doc, "14 November 2026", size=10, align="left", space_after=10)
    para(doc, "The Editors, GigaScience", size=10, align="left", space_after=10)
    para(doc, "Dear Editors,", size=11, align="left", space_after=8)

    body = [
        f"We submit as a Research article \u201c{TITLE}\u201d. It is a measurement study of four public "
        f"Visium HD sections \u2014 two of one colon tumour, plus lung cancer and mouse brain \u2014 and a "
        f"Stereo-seq mouse brain section, re-aggregated to 256 \u00b5m in eight exact steps on both "
        f"platforms.",

        f"The manuscript corrects how a statistic the field already reports is read. Overdispersion "
        f"is expected to fall as d^\u22122 when bins are merged for counts; the fitted exponent is then quoted as "
        f"\u201cthe exponent of a tissue\u201d and compared between samples. We find it falls as d^"
        f"\u2212{num('J.colon.gamma_full', 2)} in colon and d^\u2212{num('J.brain.gamma_full', 2)} in brain, that "
        f"no single correlation length survives model comparison, and that the exponent is a property of "
        f"transcripts, not tissues: it varies {num(RATIO_MIN, 1)}\u2013{num(RATIO_MAX, 1)}-fold between "
        f"the 5th and 95th percentile of genes, falls monotonically with "
        f"abundance in all three sections, and is reproducible between halves of one section "
        f"(r = {num('N.splithalf.reliability.pearson', 2)}). Being a "
        f"\u03bc\u00b2-weighted mean it is carried by "
        f"{num('min(N.gene_level.colon.composition.n_eff,N.gene_level.lung.composition.n_eff,N.gene_level.brain.composition.n_eff)', 1)} to "
        f"{num('max(N.gene_level.colon.composition.n_eff,N.gene_level.lung.composition.n_eff,N.gene_level.brain.composition.n_eff)', 0, 'int')} "
        f"effective genes, which is why the pooled ranking of the three reverses their per-gene ranking. We also show what was only asserted: the pooled exponent is a "
        f"weighted mean of the per-gene ones, rebuilt to within "
        f"{num('N.decomp.min_rel_repro', 1, 'pct')}\u2013{num('N.decomp.max_rel_repro', 1, 'pct')}, and the per-gene mean is "
        f"{num('N.robust.ratio_pooled_over_best', 1)} times more stable across regions than the pooled value.",

        f"Replication is the other half. A second section of the same tumour gives the same "
        f"pooled value (difference {num('abs(MU.contrasts.same_tissue_colon.delta)', 3)}) but not the same "
        f"per-gene median ({num('MU.pergene.shift_overall', 3)}), so the number that travels between sections "
        f"is the pooled one; and a brain section captured by a different chemistry reproduces the per-gene "
        f"structure (\u03c1 = {num('X.replication.rho_raw', 2)} over "
        f"{num('X.replication.n', 0, 'int')} genes) at a different absolute level. Merging bins still cannot "
        f"deliver the variance reduction counting statistics promise, at any bin size, because every "
        f"transcript has an exponent well below 2; depth is almost irrelevant, and the number of sections is "
        f"what has to be bought. The paper ends with a reporting rule that applies to an existing table.",

        "It suits the journal: a measurement problem on public data whose result changes how a default is "
        "reported. Every number is regenerated from a locked table by the deposited code; the controls are "
        "synthetic fields with a known answer.",

        "The manuscript is not under consideration elsewhere, is unpublished, and is approved by all authors. It "
        "reanalyses public de-identified data, so no ethics approval was required. We declare no competing "
        "interests and no funding. ChatGPT assisted with language editing; the analyses remain ours.",

        "Suggested reviewers, for spatial statistics and Visium HD expertise: Dr Shiquan "
        "Sun (Xi'an Jiaotong University, China), Dr Krzysztof Pola\u0144ski (Wellcome Sanger Institute, UK) and "
        "Dr Lambda Moses (California Institute of Technology, USA). We request no exclusions. The subject sits "
        "between statistical genomics and spatial methods; Prof. Lana Garmire (University of Michigan) or Dr "
        "Longqi Liu (BGI Research, Hangzhou) could handle it.",

        "Yours sincerely,",
    ]
    for b in body:
        para(doc, b, size=10, align="justify", space_after=5)

    para(doc, "Chaohua Tang, MD", bold=True, size=10.5, align="left", space_after=0)
    para(doc, "Department of Gynecology, The 95th Hospital of Putian, Putian, Fujian 351100, China \u00b7 "
              "tong@fjmu.edu.cn \u00b7 ORCID 0009-0004-7080-8993", size=9, align="left", space_after=0)

    out = Path(outdir) / "02_Cover_Letter.docx"
    doc.save(out)
    return out


def build_supplementary(outdir):
    doc = Document()
    set_base_style(doc)
    sec = doc.sections[0]
    sec.page_height, sec.page_width = Cm(29.7), Cm(21.0)
    sec.top_margin = sec.bottom_margin = Cm(2.5)
    sec.left_margin = sec.right_margin = Cm(2.5)
    add_page_numbers(doc)

    heading(doc, "Supplementary information", size=15, space_before=0)
    para(doc, TITLE, bold=True, size=12, align="left", space_after=12)
    para(doc, "Contents: Supplementary methods M1\u2013M6; Tables S1\u2013S9; Figures S1\u2013S6.",
         size=10, align="left")

    heading(doc, "Supplementary methods", size=13)
    sm = [
        ("M1. Why the pooled estimator is unbiased.",
         "For a gene g with mean \u03bc_g and overdispersion \u03b1_g across blocks, Var_g = \u03bc_g + "
         "\u03b1_g\u03bc_g\u00b2. With the n\u22121 denominator the sample variance is unbiased for this quantity, "
         "so Var_g \u2212 \u03bc_g is unbiased for \u03b1_g\u03bc_g\u00b2. Summing over genes and dividing by "
         "\u03a3\u03bc_g\u00b2 therefore estimates the \u03bc\u00b2-weighted average of \u03b1_g exactly, without a "
         "gene filter. The weight is also the inverse of the variance of each gene's contribution under Poisson "
         "noise, so highly expressed genes dominate the estimate and rare genes contribute almost nothing. This "
         "is why the pooled exponent is stable at 2 \u00b5m, where per-gene estimates are not."),
        ("M2. Why the block side scales with d.",
         "A within-block estimator with a fixed block size contains fewer and fewer d-blocks per block as d "
         "grows, so the within-block variance is computed on progressively smaller samples and is biased "
         "downward. On real colon data that artefact produced an apparent steepening of the exponent to "
         "\u22121.01, which disappeared when the block side was set to 8d. At 8d every block always holds 64 "
         "d-blocks, so the sampling geometry is identical at every scale."),
        ("M3. Why bootstrap intervals over blocks are too narrow.",
         "Resampling blocks treats them as independent units. A 512 \u00b5m block is large against the fitted "
         "correlation length of a few tens of micrometres, but small against the section, so resampling cannot "
         "reproduce the regional differences in composition that set the spread: the SD across four quadrants is "
         f"{num('min(D3.colon_quad_sd/B.colon.results[512].sd,D3.lung_quad_sd/B.lung.results[512].sd,D3.brain_quad_sd/B.brain.results[512].sd)', 1)} "
         f"to {num('max(D3.colon_quad_sd/B.colon.results[512].sd,D3.lung_quad_sd/B.lung.results[512].sd,D3.brain_quad_sd/B.brain.results[512].sd)', 1)} "
         "times the bootstrap SD. We therefore report the regional SD and the half-section rms as uncertainty, "
         "and treat the block bootstrap as a lower bound only."),
        ("M4. Block bootstrap over superblocks.",
         "To show how the uncertainty of \u03b3 depends on the size of the resampling unit, each section was "
         "tiled into superblocks of 512, 1024 and 2048 \u00b5m, chosen as multiples of 256 \u00b5m so that no "
         "coarse block is cut. For each superblock and scale we accumulated per-gene counts and sums of squares; "
         "200 resamples were drawn with replacement over superblocks, and \u03b1(d) and \u03b3 were recomputed "
         "from the combined moments. In colon the SD rises with superblock side, from "
         f"{num('B.colon.results[512].sd', 3)} to {num('B.colon.results[2048].sd', 3)}, without reaching the "
         f"{num('N.splithalf.quarter_sd', 3)} observed across the four bands even at 2 mm, which is the quantitative "
         "form of the statement that a section contains no internal replicate."),
        ("M5. Minimum detectable difference.",
         "Two independent sections differ with standard deviation equal to the rms of the half-section difference "
         "(each half has SD rms/\u221a2, so the difference has SD rms). For a two-sided test at 5% with 80% "
         "power, MDE = exp(2.802 \u00d7 rms/\u221an) on the log scale, where n is the number of sections per "
         "group. The same formula applied to pooled log \u03b1 gives the fold-change MDE."),
        ("M6. Synthetic data.",
         "Counts are Poisson draws on a Gaussian correlated field with a power-law spectrum, over "
         f"{num('S.expression.n_genes', 0, 'int')} genes whose means span four orders of magnitude, and with the "
         "same eight-scale re-aggregation and fitting pipeline applied without modification. The input exponent "
         f"is {num('S.expression.beta_true', 2)}. A separate family of negative-binomial simulations with "
         "independent cells (true exponent 2.00) was used in earlier work to confirm that the estimator recovers "
         "the counting limit when the counting limit is true."),
    ]
    for h, t in sm:
        heading(doc, h, size=10.5, space_before=8)
        para(doc, t)

    heading(doc, "Table S1. Per-resolution statistics", size=12, space_before=14)
    para(doc, "For each bin size: number of in-tissue blocks (colon / lung / brain), pooled \u03b1, coefficient "
              "of variation of counts per block, and median relative error of per-gene \u03b1.",
         size=8, align="left", space_after=6, spacing=1.0)
    rows = [["d (\u00b5m)", "Blocks C / L / B", "\u03b1 colon", "\u03b1 lung", "\u03b1 brain",
             "CV colon", "CV lung", "CV brain", "Err colon", "Err lung", "Err brain"]]
    for dd in [2, 4, 8, 16, 32, 64, 128, 256]:
        rows.append([str(dd),
                     f"{num(f'd.colon.power[{dd}].n_bins', 0, 'int')} / "
                     f"{num(f'd.lung.power[{dd}].n_bins', 0, 'int')} / "
                     f"{num(f'd.brain.power[{dd}].n_bins', 0, 'int')}",
                     num(f"d.colon.alpha[{dd}]", 2), num(f"d.lung.alpha[{dd}]", 2),
                     num(f"d.brain.alpha[{dd}]", 2),
                     num(f"d.colon.cv[{dd}]", 0, "pct"), num(f"d.lung.cv[{dd}]", 0, "pct"),
                     num(f"d.brain.cv[{dd}]", 0, "pct"),
                     num(f"d.colon.power[{dd}].rel_err", 0, "pct"),
                     num(f"d.lung.power[{dd}].rel_err", 0, "pct"),
                     num(f"d.brain.power[{dd}].rel_err", 0, "pct")])
    make_table(doc, rows, font=7)

    heading(doc, "Table S2. Model comparison over 2\u2013256 \u00b5m", size=12, space_before=14)
    para(doc, "Four models fitted to the eight pooled \u03b1 values. \u0394AIC is relative to the best model in "
              "that section. The F test compares the quadratic against the pure power law on one and five "
              "degrees of freedom; \u03be is the profile-likelihood estimate with its 95% interval.",
         size=8, align="left", space_after=6, spacing=1.0)
    rows = [["Section", "Model", "k", "RSS", "AIC", "\u0394AIC", "\u03be (\u00b5m)", "95% CI"]]
    for t, lab in (("colon", "Colon"), ("lung", "Lung"), ("brain", "Brain")):
        for m, mlab in (("power", "power law"), ("quad", "quadratic"), ("exp", "exp \u03be"), ("gauss", "Gauss \u03be")):
            xi = num(f"M.{t}.xi_hat", 1) if m == "exp" else (
                num(f"M.{t}.xi_ci[0]", 1) + "\u2013" + num(f"M.{t}.xi_ci[1]", 1) if m == "exp" else "\u2014")
            if m == "exp":
                ci = num(f"M.{t}.xi_ci[0]", 1) + "\u2013" + num(f"M.{t}.xi_ci[1]", 1)
                xiv = num(f"M.{t}.xi_hat", 1)
            else:
                ci, xiv = "\u2014", "\u2014"
            rows.append([lab, mlab, str(2 if m in ("power", "exp", "gauss") else 3),
                         num(f"M.{t}.rss.{m}", 4), num(f"M.{t}.aic.{m}", 2), num(f"M.{t}.daic.{m}", 2),
                         xiv, ci])
    make_table(doc, rows, font=7)
    para(doc, f"F test for the quadratic term: colon F = {num('M.colon.quad_F', 1)}, "
              f"p = {num('M.colon.quad_p_value', 4)}; lung F = {num('M.lung.quad_F', 1)}, "
              f"p = {num('M.lung.quad_p_value', 4)}; brain F = {num('M.brain.quad_F', 1)}, "
              f"p = {num('M.brain.quad_p_value', 4)}.",
         size=8, align="left", space_after=6, spacing=1.0)

    heading(doc, "Table S3. Synthetic abundance control", size=12, space_before=14)
    para(doc, "Poisson counts on a correlated field with an input exponent of 0.55; the same pipeline was applied. "
              "Medians and interquartile ranges are over genes within each expression layer.",
         size=8, align="left", space_after=6, spacing=1.0)
    rows = [["Expression layer (counts/bin)", "Genes", "Median \u03b3", "Q25", "Q75"]]
    for k, lab in zip(N.synthetic.expression.by_layer,
                      ["0.01\u20130.03", "0.03\u20130.1", "0.1\u20130.3", "0.3\u20131", "1\u20133", "3\u201310", "10\u2013100"]):
        rows.append([lab, num(f"S.expression.by_layer[{k!r}].n", 0, "int"),
                     num(f"S.expression.by_layer[{k!r}].med", 3),
                     num(f"S.expression.by_layer[{k!r}].q25", 3),
                     num(f"S.expression.by_layer[{k!r}].q75", 3)])
    rows.append(["All genes", num("S.expression.n_genes", 0, "int"),
                 num("S.expression.gamma_all.med", 3), num("S.expression.gamma_all.q25", 3),
                 num("S.expression.gamma_all.q75", 3)])
    make_table(doc, rows, widths=[4.6, 2.2, 2.6, 2.4, 2.4], font=7.5, first_col_bold=True)

    heading(doc, "Table S4. Column bands of the colon section", size=12, space_before=14)
    para(doc, "Bands are cut on multiples of 256 \u00b5m so that no coarse block is divided. \u03b1(d) is given "
              "for all eight scales, with the number of in-tissue blocks per band.",
         size=8, align="left", space_after=6, spacing=1.0)
    rows = [["Band", "Blocks", "d = 2", "4", "8", "16", "32", "64", "128", "256", "\u03b3"]]
    for k, lab in (("q1", "q1"), ("q2", "q2"), ("q3", "q3"), ("q4", "q4"), ("halfA", "half A"), ("halfB", "half B")):
        b = N.splithalf.bands[k]
        rows.append([lab, num(f"N.splithalf.bands.{k}.n_bin", 0, "int")] +
                    [num(f"N.splithalf.bands.{k}.alpha[{i}]", 2) for i in range(8)] +
                    [num(f"N.splithalf.bands.{k}.gamma", 3)])
    make_table(doc, rows, font=6.8)

    heading(doc, "Table S5. Fill convention and block bootstrap", size=12, space_before=14)
    para(doc, "Top: pooled \u03b1 under the two domain conventions, and the fraction of blocks that lie entirely "
              "inside the tissue. Bottom: distribution of \u03b3 over 200 block-bootstrap resamples, for superblocks "
              "of 512, 1024 and 2048 \u00b5m, in each section.",
         size=8, align="left", space_after=6, spacing=1.0)
    rows = [["Section", "d (\u00b5m)", "Full fraction", "\u03b1, all in-tissue", "\u03b1, full blocks only",
             "Ratio"]]
    for t, lab in (("colon", "Colon"), ("lung", "Lung"), ("brain", "Brain")):
        for dd in [2, 8, 32, 64, 128, 256]:
            rows.append([lab if dd == 2 else "", str(dd),
                         num(f"F.{t}.by_scale[{dd}].fill_frac", 3),
                         num(f"F.{t}.by_scale[{dd}].alpha_any", 2),
                         num(f"F.{t}.by_scale[{dd}].alpha_full", 2),
                         num(f"F.{t}.by_scale[{dd}].ratio_full_over_any", 3)])
    make_table(doc, rows, font=7)
    para(doc, f"Effect of the convention on the fitted exponent: colon {num('F.colon.gamma_delta', 3)}, "
              f"lung {num('F.lung.gamma_delta', 3)}, brain {num('F.brain.gamma_delta', 3)}.",
         size=8, align="left", space_after=8, spacing=1.0)
    if N.get("boot"):
        rows = [["Section", "Superblock (\u00b5m)", "Blocks", "Mean \u03b3", "SD", "2.5%", "97.5%"]]
        for t, lab in (("colon", "Colon"), ("lung", "Lung"), ("brain", "Brain")):
            if t not in N.boot.by_tissue:
                continue
            for bs in N.boot.by_tissue[t].block_sizes:
                rows.append([lab if bs == N.boot.by_tissue[t].block_sizes[0] else "", str(bs),
                             num(f"B.{t}.results[{bs}].nblk", 0, "int"),
                             num(f"B.{t}.results[{bs}].mean", 3),
                             num(f"B.{t}.results[{bs}].sd", 4),
                             num(f"B.{t}.results[{bs}].q025", 3),
                             num(f"B.{t}.results[{bs}].q975", 3)])
        make_table(doc, rows, widths=[2.2, 3.0, 2.2, 2.4, 2.2, 2.2, 2.2], font=7.5, first_col_bold=True)
    para(doc, f"For comparison, the SD across the four bands is {num('N.splithalf.quarter_sd', 3)} and the range "
              f"is {num('N.splithalf.quarter_range', 3)}.", size=8, align="left", space_after=6, spacing=1.0)

    heading(doc, "Table S6. Per-gene reliability by expression layer (colon)", size=12, space_before=14)
    para(doc, "Per-gene exponents fitted independently in the two halves of the colon section, stratified by mean "
              "counts per 32 \u00b5m bin. Reliability is the Pearson correlation between halves.",
         size=8, align="left", space_after=6, spacing=1.0)
    rows = [["Expression layer (counts/bin)", "Genes", "r", "MAD", "Median, half A", "Median, half B"]]
    for lay in N.reliability_layers.layers:
        _i = lay["i"]
        lab = (f"{num(f'RL.layers[{_i}].lo', 2, 'g')}"
               f"\u2013{num(f'RL.layers[{_i}].hi', 2, 'g')}")
        if lay["pearson"] is None:
            rows.append([lab, num(f"RL.layers[{lay['i']}].n", 0, "int"), "\u2014", "\u2014", "\u2014", "\u2014"])
        else:
            rows.append([lab, num(f"RL.layers[{lay['i']}].n", 0, "int"),
                         num(f"RL.layers[{lay['i']}].pearson", 2),
                         num(f"RL.layers[{lay['i']}].mad", 3),
                         num(f"RL.layers[{lay['i']}].med_a", 2),
                         num(f"RL.layers[{lay['i']}].med_b", 2)])
    make_table(doc, rows, widths=[4.6, 2.0, 1.8, 1.8, 3.0, 3.0], font=7.5, first_col_bold=True)

    heading(doc, "Table S7. Decomposition of the pooled exponent", size=12, space_before=14)
    para(doc, "The pooled exponent is rebuilt as the weighted mean of the per-gene exponents, "
              "\u03a3_g v_g\u03b3_g with v_g averaged over the eight scales and normalised over the genes that "
              "carry a computable exponent. Residual is the relative difference from the reported pooled value. "
              "n_eff is 1/\u03a3w\u00b2 (variance weights) or 1/\u03a3v\u00b2 (weights the exponent carries). "
              "The last column is the share of the distance between the weighted and the unweighted summary that "
              "is attributable to the weights themselves drifting with scale, rather than to the weights at a "
              "fixed scale.",
         size=8, align="left", space_after=6, spacing=1.0)
    rows = [["Section", "\u03b3 pooled", "\u03a3 v\u00b7\u03b3", "Residual", "n_eff, \u03bc\u00b2",
             "n_eff, \u03b1", "Leading gene (\u03b1 share)", "Its \u03bc\u00b2 share", "Weight drift"]]
    for t, lab in (("colon", "Colon"), ("lung", "Lung"), ("brain", "Brain")):
        rows.append([lab,
                     num(f"N.decomp.by_tissue.{t}.gamma_pool", 4),
                     num(f"N.decomp.by_tissue.{t}.gamma_weighted", 4),
                     num(f"N.decomp.by_tissue.{t}.rel_diff", 2, "pct"),
                     num(f"N.decomp.by_tissue.{t}.n_eff_mu2", 1),
                     num(f"N.decomp.by_tissue.{t}.n_eff_alpha", 1),
                     f"{N.decomp.by_tissue[t]['top1_gene']} "
                     f"({num(f'N.decomp.by_tissue.{t}.top1_alpha_share', 1, 'pct')})",
                     num(f"N.decomp.by_tissue.{t}.top1_mu2_share", 1, "pct"),
                     num(f"N.decomp.by_tissue.{t}.drift_share", 0, "pct")])
    make_table(doc, rows, font=7, first_col_bold=True)

    heading(doc, "Table S8. Stability of five summary scalars across four bands (colon)", size=12, space_before=14)
    para(doc, "Each scalar is recomputed inside each of the four column bands of Fig. 5a and inside the whole "
              "section. SD and CV are across the four bands. A scalar intended for comparison between samples has "
              "to be stable within one.",
         size=8, align="left", space_after=6, spacing=1.0)
    bidx = [i for i, b in enumerate(N.robust.bands) if b["name"].startswith("q")]
    bnames = ["q1", "q2", "q3", "q4"]
    rows = [["Scalar", "q1", "q2", "q3", "q4", "Whole section", "SD", "CV"]]
    for key, lab in (("pooled", "pooled \u03b3"), ("median", "per-gene median"),
                     ("mean", "per-gene mean"), ("trimmed", "5% trimmed mean"),
                     ("wmedian", "\u03bc\u00b2-weighted median")):
        rows.append([lab] +
                    [num(f"N.robust.bands[{i}].{key}", 3) for i in bidx] +
                    [num(f"N.robust.full.{key}", 3),
                     num(f"N.robust.summary.{key}.sd", 3),
                     num(f"N.robust.summary.{key}.cv", 1, "pct")])
    make_table(doc, rows, font=7, first_col_bold=True)

    if N.spatial_proxy.by_tissue:
        heading(doc, "Table S9. Spatial statistics against the per-gene exponent", size=12, space_before=14)
        para(doc, "All statistics are computed on the 32 \u00b5m grid over genes with a computable exponent. "
                  "\u03be is the exponential range of the semivariogram, with the fraction of genes that return "
                  "the lowest admissible value shown in the last column; I is Moran's I on the eight-neighbour "
                  "adjacency. Within-band is the median Spearman correlation inside 20 equal-count bands of mean "
                  "expression, with the range across bands in brackets; semiparametric is the correlation of the "
                  "residuals after removing a binned mean of log abundance in 40 bands.",
             size=8, align="left", space_after=6, spacing=1.0)
        rows = [["Section", "Genes", "\u03c1(\u03b3, log \u03bc)", "\u03c1(\u03b3, \u03be)", "partial \u03be",
                 "\u03c1(\u03b3, I)", "Within-band \u03c1(\u03b3, I)", "Semiparametric",
                 "median \u03be (\u00b5m)", "At floor"]]
        for t, lab in (("colon", "Colon"), ("lung", "Lung"), ("brain", "Brain")):
            if t not in N.spatial_proxy.by_tissue:
                continue
            rows.append([lab,
                         num(f"N.spatial_proxy.by_tissue.{t}.n_ok", 0, "int"),
                         num(f"N.spatial_proxy.by_tissue.{t}.rho_gamma_mu", 2),
                         num(f"N.spatial_proxy.by_tissue.{t}.rho_gamma_xi", 2),
                         num(f"N.spatial_proxy.by_tissue.{t}.partial_rho_xi", 2),
                         num(f"N.spatial_proxy.by_tissue.{t}.rho_gamma_moran", 2),
                         f"{num(f'N.spatial_proxy.by_tissue.{t}.band_rho[1].median', 2)} "
                         f"({num(f'N.spatial_proxy.by_tissue.{t}.band_rho[1].lo', 2)}, "
                         f"{num(f'N.spatial_proxy.by_tissue.{t}.band_rho[1].hi', 2)})",
                         num(f"N.spatial_proxy.by_tissue.{t}.semiparam.pearson", 2),
                         num(f"N.spatial_proxy.by_tissue.{t}.xi_median_um", 0, "int"),
                         num(f"N.spatial_proxy.by_tissue.{t}.xi_frac_at_floor", 1, "pct")])
        make_table(doc, rows, font=6.8, first_col_bold=True)

    # ---------------- Table S10: cross-platform criteria
    heading(doc, "Table S10. Criteria for the cross-platform comparison, frozen before the Stereo-seq "
                 "moments existed", size=12, space_before=14)
    para(doc, "The criteria are scale-free because an absolute exponent is a property of the platform. Each was "
              "written down before the Stereo-seq section was processed; the version of H1 in the first row of "
              "the frozen list was replaced before any result was seen, and both are reported. The two attacks "
              "were also specified in advance and neither changes the verdict.",
         size=8, align="left", space_after=6, spacing=1.0)
    rows = [["Criterion", "Rule", "Stereo-seq brain", "Visium HD brain", "Verdict"]]
    rows.append(["H1a effective genes",
                 "n_eff(\u03b1 weight) / genes tested < 1%",
                 num("X.concentration.n_eff_over_genes", 2, "pct"),
                 num("X.visium.n_eff_alpha/X.visium.n_genes_ok", 2, "pct"),
                 "PASS" if X["concentration"]["h1a_pass"] else "FAIL"])
    rows.append(["H1b dominance",
                 "heaviest transcript > 10% of the exponent",
                 num("X.stereo.top1_alpha_share", 1, "pct"),
                 num("X.visium.top1_alpha_share", 1, "pct"),
                 "PASS" if X["concentration"]["h1b_pass"] else "FAIL"])
    rows.append(["H2 ordering",
                 "pooled value < per-gene median",
                 f"{num('X.stereo.gamma_pool', 3)} < {num('X.stereo.gamma_median', 3)}",
                 f"{num('X.visium.gamma_pool', 3)} < {num('X.visium.gamma_median', 3)}",
                 "PASS" if X["concentration"]["h2_pass"] else "FAIL"])
    rows.append(["H3 level",
                 "not tested; absolute values may differ",
                 num("X.stereo.gamma_pool", 3),
                 num("X.visium.gamma_pool", 3),
                 f"\u0394 = {num('X.delta_gamma_pool', 3)}"])
    rows.append(["H4 per-gene replication",
                 f"rank \u03c1 of \u03b3_g > {num('X.replication.h4_threshold', 2)} over shared genes",
                 f"{num('X.replication.rho_raw', 3)} (n = {num('X.replication.n', 0, 'int')})",
                 "same genes",
                 "PASS" if X["replication"]["h4_raw_pass"] else "FAIL"])
    rows.append(["H4 after abundance",
                 "same, on the rank residual of log abundance",
                 num("X.replication.rho_partial_abundance", 3),
                 "same genes",
                 "PASS" if X["replication"]["h4_partial_pass"] else "FAIL"])
    rows.append(["Attack A: rRNA",
                 "drop Gm42418, Gm26917",
                 f"{X['without_rrna']['top1']} {num('X.without_rrna.top1_share', 1, 'pct')}, "
                 f"n_eff {num('X.without_rrna.n_eff_alpha', 2)}",
                 f"{X['visium']['top1']} {num('X.visium.top1_alpha_share', 1, 'pct')}",
                 "PASS" if X["without_rrna"]["h1a_pass"] and X["without_rrna"]["h1b_pass"] else "FAIL"])
    rows.append(["Attack B: abundance",
                 "\u03c1 of log abundance between platforms",
                 num("X.replication.rho_abundance", 3),
                 "\u2014",
                 f"{num('X.replication.rho_abundance**2', 0, 'pct')} of raw \u03c1"])
    rows.append(["Top 100 genes",
                 "share of the pooled exponent",
                 num("X.top100_share_stereo", 1, "pct"),
                 num("X.top100_share_visium", 1, "pct"),
                 "reported"])
    make_table(doc, rows, font=6.6, first_col_bold=True)

    # ---------------- Table S11: four sections
    heading(doc, "Table S11. Four sections: estimates, intervals and contrasts", size=12, space_before=14)
    para(doc, f"Intervals are 95% from a spatial block bootstrap over contiguous "
              f"{num('MU.block_um', 0, 'int')} \u00b5m superblocks with B = {num('MU.n_boot', 0, 'int')} "
              f"resamples. The two error models give different answers and both are shown: the block bootstrap "
              f"measures the sampling error of one section, the regional SD measures the spread of four bands "
              f"inside it and is the conservative choice. ICC over four bands is reported for completeness and "
              f"not used, because the bands are non-overlapping subregions rather than exchangeable replicates.",
         size=8, align="left", space_after=6, spacing=1.0)
    rows = [["Section", "Tissue", "Genes tested", "Fitted", "\u03b3 pooled", "SE", "95% CI",
             "Per-gene median", "Effective genes", "UMI"]]
    for k, lab, tis in (("colon", "Colon P0", "colon cancer"),
                        ("Visium_HD_Human_Colon_Cancer_P1", "Colon P1", "colon cancer"),
                        ("lung", "Lung", "lung cancer"),
                        ("brain", "Brain", "mouse, normal")):
        s = MU["sections"][k]
        rows.append([lab, tis,
                     num(f"MU.sections.{k}.n_genes", 0, "int"),
                     num(f"MU.sections.{k}.n_fit", 0, "int"),
                     num(f"MU.sections.{k}.gamma_pool", 3),
                     num(f"MU.sections.{k}.boot_sd", 3),
                     f"{num(f'MU.sections.{k}.boot_lo', 3)}\u2013{num(f'MU.sections.{k}.boot_hi', 3)}",
                     num(f"MU.sections.{k}.gamma_median", 3),
                     num(f"MU.sections.{k}.n_eff_alpha", 2),
                     num(f"MU.sections.{k}.umi_total", 0, "int")])
    make_table(doc, rows, font=6.6, first_col_bold=True)

    rows = [["Contrast", "\u0394", "SE, block bootstrap", "z", "SE, regional SD", "z"]]
    for key, lab in (("same_tissue_colon", "Colon P1 \u2212 P0"),
                     ("colon_vs_lung", "Colon \u2212 lung"),
                     ("colon_vs_brain", "Colon \u2212 brain"),
                     ("lung_vs_brain", "Lung \u2212 brain")):
        a = MU["contrasts"][key]
        b = MU["contrasts_regional"][key]
        rows.append([lab, num(f"MU.contrasts.{key}.delta", 3),
                     num(f"MU.contrasts.{key}.se_diff", 3),
                     num(f"MU.contrasts.{key}.z", 1),
                     num(f"MU.contrasts_regional.{key}.se_diff", 3),
                     num(f"MU.contrasts_regional.{key}.z", 1)])
    rows.append(["MDE, two sections, 80%", num("MU.mde.value", 3),
                 num("MU.mde.se_typ", 3), "\u2014", "\u2014", "\u2014"])
    make_table(doc, rows, font=6.6, first_col_bold=True)

    # ---------------- Table S12: per-gene shift by abundance
    heading(doc, "Table S12. Per-gene exponent by abundance decile, two sections of one tumour",
             size=12, space_before=14)
    para(doc, "Deciles are defined on the mean counts at 256 \u00b5m of the first section and applied unchanged "
              "to the second, so a decile holds genes of comparable abundance in both. The shift is the "
              "difference of medians, second minus first; intervals resample genes within the decile.",
         size=8, align="left", space_after=6, spacing=1.0)
    rows = [["Decile", "Mean counts at 256 \u00b5m", "Genes P0", "Genes P1",
             "Median P0", "Median P1", "Shift", "95% CI"]]
    nd = len(MU["pergene"]["deciles"])
    for i in range(nd):
        lo = num(f"MU.pergene.deciles[{i}].mu_lo", 2)
        hi = ("\u2014" if i == nd - 1 else num(f"MU.pergene.deciles[{i}].mu_hi", 2))
        rows.append([num(f"MU.pergene.deciles[{i}].decile", 0, "g"), f"{lo}\u2013{hi}",
                     num(f"MU.pergene.deciles[{i}].n0", 0, "int"),
                     num(f"MU.pergene.deciles[{i}].n1", 0, "int"),
                     num(f"MU.pergene.deciles[{i}].med0", 3),
                     num(f"MU.pergene.deciles[{i}].med1", 3),
                     num(f"MU.pergene.deciles[{i}].shift", 3),
                     f"{num(f'MU.pergene.deciles[{i}].shift_lo', 3)} to "
                     f"{num(f'MU.pergene.deciles[{i}].shift_hi', 3)}"])
    rows.append(["All genes", "\u2014",
                 num("MU.pergene.n_fit_p0", 0, "int"),
                 num("MU.pergene.n_fit_p1", 0, "int"),
                 num("MU.pergene.median_p0", 3),
                 num("MU.pergene.median_p1", 3),
                 num("MU.pergene.shift_overall", 3), "\u2014"])
    make_table(doc, rows, font=6.6, first_col_bold=True)

    heading(doc, "Supplementary figures", size=13, space_before=14)
    para(doc, "Each figure is reproduced below at 300 dpi and is also supplied as a separate PDF in 03_Figures/. "
              "Legends are listed in the main manuscript.", size=10, align="left")
    for tag, width in (("FigS1", Cm(11.0)), ("FigS2", Cm(11.0)), ("FigS3", Cm(11.0)),
                       ("FigS4", Cm(16.5)), ("FigS5", Cm(16.5)), ("FigS6", Cm(16.5))):
        png = BUILD / "figs" / f"{tag}.png"
        if not png.exists():
            raise FileNotFoundError(png)
        doc.add_picture(str(png), width=width)
        doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER

    out = Path(outdir) / "04_Supplementary.docx"
    doc.save(out)
    return out


def build_guide(outdir):
    txt = f"""# 投稿指南 · A 题（空间转录组噪声标度律）v2.2

版本 v2.2 ｜ 2026-11-14 ｜ 目标刊：**GigaScience**（Research article）
本文件是操作手册，不是稿件的一部分，不上传。

---

## 零、这一版改了什么

v1.0 的核心归因是错的，而且是被自己的数据推翻的：原稿说"这个指数由组织结构决定、组织之间可比"，
实测指数由**哪几个转录本最丰**决定，组织身份只解释 {num('V.pooled.frac_tissue', 1, 'pct')}。v2.0 换主角：

| | v1.0 | v2.2 | v3.1 |
|---|---|---|
| 主张 | 组织设定噪声分辨率指数 | 汇总量是**组成统计量**；逐基因层面由丰度决定；设计结论不变 | 组成在**换平台、换切片**下成立；可移植的只有 pooled 那一个数 |
| 证据 | 一条 α(d) 曲线 + 四象限 | 逐基因 γ、权重分解、留一基因、分带、深度、合成、填充率、块自助、排序反转 | 前项 + 第二平台（Stereo-seq）+ 第二张切片 + 跨片逐基因分解 |
| 图 | 6 主 4 补 | 6 主 6 补 | 8 主 6 补（新增跨平台、切片复现） |
| 表 | 2 主 3 补 | 3 主 9 补 | 3 主 12 补 |
| 数字 | `NUMBERS_LOCK_v1.json` | `NUMBERS_LOCK_v5.json`（v1 的 synth/spread/window 已并入） | 同一份 v5，增补 `cross` / `multi` 两块 |

### v3.1 只加了一件事：报告流程变成能跑的包

1. **`scaling_report` 1.0.0 独立成包**（MIT，只依赖 numpy）。输入是逐基因在 bin 阶梯上的均值与
   方差，输出 pooled 指数、逐基因分布、两个有效基因数、同一 pooled 值允许的中位数区间，以及
   四条填好数字的报告规则。不需要原始矩阵，换平台也能用。
2. **它复现本文的数字，不是"参考实现"**：结肠 0.5547、IGKC 87.7%、n=13,115、中位数 0.428
   （IQR 0.338–0.585）、d(med)/d(anchor) = −7.16、区间 0.000–1.000；两片对照 Δ = −0.0062、
   逐基因中位数 +0.190。README 与测试里逐条钉死。
3. **12 项测试**，其中三项用"答案由构造而来"的合成数据要求机器精度还原（单基因、全体同指数、
   逐基因各自指数），一项用已发表的结肠数字做回归。合成数据也随包发布。
4. **它拒绝给单切片置信区间**，把 v3.0 第 4 条那个错误（i.i.d. 瓦片自助）写进了工具的行为里：
   要么超块自助，要么每组建 ≥2 片，要么只报差值。审稿人拿到工具也就拿到了这条约束。
5. **顺手修掉两处过期全称判断**：Data availability 还写着 "three datasets"（现为四张 10x 切片
   加一张 Stereo-seq）；投稿指南第 3 步的操作命令仍指向 v2 时代的 45/47/48/49/31 号脚本。
   两处都是"文档没跟着数据长"，不是数字错。

### v3.0 在 v2.2 基础上加了五件事（两项新数据、一处结论反转、一处旧说法纠错）

1. **第二平台：Stereo-seq 小鼠全脑**（C04042E2，0.5 µm DNB，2–256 µm 八级同梯）。
   判据在跑数之前写死（H1a/H1b/H2/H4，Table S10）：有效基因占受检基因 **0.02%**、头号转录本
   占指数 **54.0%**、pooled **0.565** < 逐基因中位数 **1.532**、11,230 个共测基因上逐基因 ρ =
   **0.769**。两项证伪攻击都过：剔掉 rRNA 衍生位点后头号基因变成 **Ttr（27.0%）**——正是 Visium
   脑片的头号基因，两平台在去污染后**指向同一个转录本**；把两平台各自的 log 丰度秩残差化后
   ρ 仍为 **0.615**，不是"都按丰度排序"的假象。γ 的绝对水平差 **0.202**，这正是标题该说的是
   "组成"而不是"某个平台"的证据（新增 Fig. 7、Table S10）。
2. **第二张切片：同组织两片**。P1（12.32 GB）走完全相同的脚本与过滤，四片一起做
   512 µm 超块自助法（B=200）：结肠 0.555±0.012、P1 0.549±0.026、肺 0.546±0.015、
   脑 0.364±0.011。同组织差 **0.006（z=0.21）**，结肠−脑差 **0.191（z=11.5）**，
   两片设计的 MDE **0.075**——组织差是它的 2.5 倍。"n=1/组织"这条批评第一次有了正面回答
   （新增 Fig. 8、Table S11）。
3. **一处结论按数据反转，不是修补**：逐基因中位数在两片之间差 **0.190**（0.428 对 0.618），
   是同一对切片 pooled 差值的 31 倍；而且在丰度最低的十分位高 **0.050**、最高的十分位低
   **0.029**——两端反向，所以**不是**一个可以校准掉的平移（Table S12、Fig. 8c）。
   因此 v2.2 那句"推荐报逐基因均值，它比 pooled 更稳"只在**片内**成立（片内基因集固定），
   到**片间**就反过来了。正文改成两维口径：片内用逐基因均值 + IQR，片间用 pooled，
   并明说单一标量都不可移植。
4. **旧 Methods 里有一句和论文自己的结果矛盾，删掉重写**。v2.2 写"512 µm 瓦片小于相关长度，
   所以自助区间偏窄"；但论文自己的 Fig. 2 估出 ξ = 30.8–63.8 µm，512 µm 是它的 8–17 倍。
   区间偏窄的真实原因是 i.i.d. 瓦片把空间异质性当成零，不是瓦片太小。改写后同时给出两种误差
   模型（超块自助 / 区域带 SD），结肠−脑在保守口径下 z=2.83，结论不变但更弱，两个都报。
   ICC 也按原计划算了（0.00），但结构性不成立（四条带是互不重叠的空间子区，不可交换），
   正文与 Table S11 里都写明"算了但不用"。
5. **摘要与局限段同步**：摘要加"第二片复现 pooled、不复现逐基因中位数"与"换化学体系后
   逐基因结构复现、绝对水平不复现"两句，结论段改成"可移植的是 pooled"；局限段原来说
   "all three datasets come from one platform"，现在不成立，改成按平台差异陈述。

### v2.2 在 v2.1 基础上补了八处（前五处把"推理"换成"算过"，后三处补齐承诺与口径）
1. **汇总量现在是算出来的加权平均，不是论证出来的**。α_pool 在每个尺度上都是逐基因 α 的凸组合，
   恒等式误差到机器精度；把八个尺度分别重构再加权，得到 0.5563 / 0.5513 / 0.3732，与报告的
   0.555 / 0.546 / 0.364 相差 0.3%–2.6%（新增 Table S7、Fig. S6a）。残差来源写清楚：入算基因
   集合不同 + 权重随尺度漂移（结肠 4%、肺 48%）。
2. **有效基因数换成两个口径**。论文原来只报 1/Σw²（2.5 / 7.7 / 82.3），但指数真正携带的权重是
   v_g = w_gα_g/α_pool，按它算是 **1.3 / 1.7 / 31.2**。结肠 IGKC 从占 Σμ² 的 62.4% 升到占指数的
   87.7%，前三名合计 92.0%；脑里两个排序甚至指向不同基因（Σμ² 是 Camk2n1，指数是 Ttr 12.7%）。
   方向和主张一致，但更锋利，也堵住"你的权重到底是哪个"这一问。
3. **给了可替代的标量**，批判变建设。结肠四带重算五种标量：pooled γ 的带间 SD 是 0.056
   （CV 9.7%），逐基因均值 0.038（7.6%）、5% 截尾均值 0.038、中位数 0.042；μ² 加权中位数最差
   0.115（22.0%）。**推荐报逐基因均值 + IQR，它比 pooled 稳 1.5 倍**（新增 Table S8、Fig. S6b）。
   均值与截尾均值在这个精度上分不出高下（0.0377 对 0.0376），推荐均值的理由是它不丢基因、不带
   调参常数，不是因为它更稳。
4. **空间纹理这条写成了"两有一无"，没有硬说成机制**。32 µm 上逐基因 Moran's I 与 γ 在控制丰度
   后仍相关：结肠 −0.58、脑 −0.72（半参数偏相关 −0.41 / −0.62），**肺没有（−0.12 / −0.07）**。
   丰度梯度三片都在、纹理关联只有两片，所以写"丰度不能还原为单一滞后上的聚集程度"，而不是
   "γ 就是空间自相关"。变异函数的相关长度在本分辨率估不出来（47.1% 顶在下界），如实写成
   "有界而非估计"（新增 Table S9）。
5. **一处分母不一致的隐患顺手修掉**：原文"in brain it is an average over a thousand genes"与新的
   有效基因数 31.2 冲突，改成"over 31 effective genes"。
6. **三条躺在文献库里、正文却没引的条目补到位**（43 → 46 条，现在 46/46 都有引用位置）。
   不是凑数，各补在论证真正需要它的地方：Svensson 2020（液滴 scRNA-seq 的零不是"零膨胀"，
   是过散计数）接在计数模型那句之后；Lähnemann 2020（单细胞数据科学十一大挑战）用来说明
   "开放问题清单从拿到计数矩阵开始，把颗粒度当既定"；Wu 2021（乳腺癌单细胞 + 空间图谱）用来
   支撑"权重不只在区域之间不同，在患者之间也不同"。此前这三条在本版文献表里悬空，构建日志会报
   `uncited`。
7. **投稿信压到 1 页（469 词）**。v2.1 的投稿信实测是 2 页（492 词），而指南写"1 页"——承诺与
   事实不符，本轮按逐行测量收干净：删的是重复表述，不是内容，4 位建议人（3 审稿人 + 2 编委）、
   三条声明、全部关键数字都还在。
8. **修掉一处自己造出来的口径分裂，代价是结论变弱**。同一句"结肠四带"在两处给出了两个 SD：
   前面 0.056、新节 0.064。根因是 51 号脚本把四带的首尾边界改成了数据真实范围
   （`edges[0], edges[-1] = c_lo, c_hi + 1`），违反它自己写明的"按 128 的整数倍切、不切 256 块"
   规则；只有第四带吃到这条越界，所以 q1–q3 一直对得上、只有 q4 不同。统一回预先写下的规则后：
   pooled SD 0.056（CV 9.7%），逐基因均值比 pooled 稳的倍数从 **2.5 降到 1.5**，"最稳标量"也从
   均值变成均值与截尾均值打平。取更弱的那个口径，是因为规则先写在代码里、而结果后出现。
   顺带把 37 号跑满三套数据：肺 0.051、脑 0.038，结肠不是孤例，这条已写进正文。

### v2.1 在 v2.0 基础上改了五处（都是被自己的数字逼出来的）

1. **标题换了**：`…is a composition statistic`。原来的 "not tissue architecture" 站不住——结肠那个数
   恰恰是浆细胞成团这种结构的产物。
2. **承认排序反转**（Results 新增一段，摘要与投稿信各加一句）：pooled 排序是结肠 0.555 > 肺 0.546 >
   脑 0.364，逐基因中位数排序是脑 0.65 > 结肠 0.43 > 肺 0.35，**恰好相反**。这是全文最硬的一条证据，
   v2.0 却没写。
3. **结论段的机制句改对了**：原来写"主导转录本在数百微米上平滑"，与 Fig. S5 的"结肠局灶、脑内铺开"
   自相矛盾。现在按证据写：IGKC 自身 γ = {num('J.colon.top1_gamma', 2)}（局灶，把 pooled 抬高
   {num('D3.colon_pooled_minus_median', 2)}），Camk2n1 γ = {num('J.brain.top1_gamma', 2)}（铺开，把 pooled
   压低 {num('abs(D3.brain_pooled_minus_median)', 2)}）。
4. **CV 不再是"同一测量的无模型版"**：它给出的隐含指数是 0.25 / 0.21 / 0.32，比 α 口径低一半。现在
   明确写成"另一个统计量"，并把三个数摆出来，避免审稿人自己算出来。
5. **Table 3 补全**：基因数改成**实际入算数**（最低层 862/1,339，原来列 1,337 会让人以为全都用了），
   并补一行 `≥10`，把剩下 {num('L.by_tissue.colon[8].share', 1, 'pct')} / {num('L.by_tissue.lung[8].share', 1, 'pct')} /
   {num('L.by_tissue.brain[8].share', 1, 'pct')} 的权重交代清楚——权重合计到 100%。顺带把悬空的
   "(below)" 换成真实的方差分解 81.5% / 1.6% / 16.9%。

---

## 一、这个包里有什么

```
00_投稿指南.md        本文件（不上传）
01_Manuscript.docx    正文 + 图注 + 声明 + 参考文献 + 表 1–3
02_Cover_Letter.docx  投稿信，1 页（含 3 位建议审稿人 + 2 位建议编委）
03_Figures/           Fig1–Fig8 + FigS1–FigS6，PDF（矢量）+ PNG（300 dpi）
04_Supplementary.docx 补充方法 M1–M9、表 S1–S12、图 S1–S6
05_Code/              全部分析脚本 + NUMBERS_LOCK_v5.json + README（图件与文档可由锁复现）
05_Code/scaling_report-1.0.0/  独立报告工具（MIT，numpy），只需矩量表即可复现 Table 3
```

---

## 二、投稿前必做（约 20 分钟，只有这一件是真必须）

**建代码仓库。** GigaScience 要求数据与代码可复现，代码不在投稿时公开，编辑很可能直接退回。

1. 打开 https://github.com/new ，仓库名建议 `spatial-noise-scaling`，选 **Public**，许可 MIT 或 GPL-3。
2. 上传 `05_Code` 整个目录 + `NUMBERS_LOCK_v5.json`，以及 `42_build_lock_v3.py` 读取的 `02_OUT` 下各
   JSON / npz（gene_moments、gene_curves、splithalf_genes、layer_maps、xi_model_8pt、fill_sensitivity、
   block_boot、layer_variance）。README 已按两级写明：只给锁，能重出全部图件文档并重跑三道守卫；
   要重算锁本身，得先跑上游脚本。
3. 记下仓库地址。**不要手改 docx**——正文里每个数字都由脚本注入，手改会与数字锁脱钩。
   正确做法：在 `78_build_docs_v31.py` 的 `("Code availability", ...)` 段末尾补一句
   `"The repository is at <你的地址>."`，然后依次重跑

   ```
   python 78_build_docs_v31.py      # 重出四件文档
   python 70_guard_numbers_v3.py    # 数字守卫
   python 71_guard_language_v3.py   # 语言守卫
   python 72_guard_carriers_v3.py   # 载体守卫
   python 73_package_v3.py          # 重新打包
   ```

4. 三个守卫都 PASS 之后再上传；任何一条 FAIL 都说明稿子与数字锁不一致。

> 包里的 `05_Code/README.md` 已写清三行命令（`42` → `43` → `45`）与预期输出（Fig1–Fig8、Table 1–3、FigS1–S6、Table S1–S12），
> 直接把整个 `05_Code/` 推到仓库即可。

---

## 三、系统里每一步填什么

投稿入口：`https://academic.oup.com/gigascience` → 右上 **Submit** → 跳到 ScholarOne。
用 **tong@fjmu.edu.cn** 注册（通讯作者必须用这个邮箱）。

| 步骤 | 填什么 |
|---|---|
| Article type | **Research**（不要选 Data Note / Technical Note） |
| Title | 复制 `01_Manuscript.docx` 第一页标题 |
| Running title | The pooled exponent is a composition statistic |
| Abstract | 复制 Abstract 整段（Background / Results / Conclusions 三个小标题保留） |
| Keywords | 6 个，首页已列 |
| Authors | Chaohua Tang¹,\\*,\\*\\*、Fangfang Wu¹,\\*、Yingying Zheng¹,\\*、Lingling Zhan¹、Hong Tan¹、Lijing Lin¹、Meimei Wen¹；前三位勾 equal contribution，唐朝华勾 corresponding |
| Affiliation | Department of Gynecology, The 95th Hospital of Putian, No. 485 Dongyan Road, Licheng District, Putian City, Fujian Province 351100, China（逐字照抄） |
| Funding | None |
| Cover letter | 上传 `02_Cover_Letter.docx` |
| Manuscript | 上传 `01_Manuscript.docx`（单文件，含标题页、正文、图注、声明、文献、三张表） |
| Figures | 逐张上传 `03_Figures/Fig1.pdf` … `Fig6.pdf`，类型 Figure |
| Supplementary | 上传 `04_Supplementary.docx`，类型 Supplementary Material |
| Suggested reviewers | 3 位，见投稿信末段。**邮箱请在 ScholarOne 里按姓名检索后选**，不要手填猜测地址 |
| Opposed reviewers | 留空 |
| Data availability | 选 "Data available in public repository"，填三个 10x 数据集名 |
| Code availability | 填第二步的仓库地址（先按第二步重新生成再上传） |
| 是否之前投过 | 否 |

---

## 四、投稿前 5 分钟自检

- [ ] 通讯邮箱 tong@fjmu.edu.cn，ORCID 0009-0004-7080-8993
- [ ] 7 位作者顺序、单位拼写与稿件完全一致
- [ ] 6 主图 + 5 补充图 PDF 都能打开
- [ ] 稿件里没有中文、没有批注、没有修订痕迹
- [ ] 摘要 ≤ 250 词（首页不写字数，系统里会算）
- [ ] Code availability 里的网址能点开

---

## 五、这一版的证据链（审稿人问什么，看哪张图）

| 审稿人会问 | 答案在哪 |
|---|---|
| 指数到底多浅？ | Fig 1b，γ = 0.55 / 0.55 / 0.36 |
| 有没有单一相关长度？ | Fig 2a–c，ΔAIC ≥ 20，ξ 区间重叠 |
| 逐基因的 γ 是不是噪声？ | Fig 3c（两半 r = 0.79）+ Fig S2（合成无漂移）+ 深度 Thinning |
| 为什么 pooled 值不能比？ | Fig 4a–c（n_eff 2.5–82，去掉 IGKC 动 0.081）+ Fig 5a（同一切片四带差 0.131） |
| 16,384 倍是名义值？ | Table S5 上半（严格填充口径） |
| 一块切片内部有几个独立重复？ | Fig 5d（块自助 SD 随块边长上升） |
| 那到底该怎么设计？ | Fig 6b（MDE vs 切片数）、Fig 6c（各操控的效应量排序） |

---

## 六、投稿之后

- 初审一般 2–6 周。GigaScience 实行开放同行评议，审稿报告会随文章公开。
- 收到意见后把原文发我，按"先锐评、再一次改到位"的流程处理。
- 录用后 APC 约 **$2,638**（≈¥1.9 万），走医院报销，留存 invoice 与付款凭证。

---

## 七、这一版做了什么、没做什么

**做了**：三套数据统一到 2–256 µm 口径重算（模型比较原来是 9 点含 512 µm，已改 8 点）；
全部数字由脚本从 `NUMBERS_LOCK_v5.json` 现场注入，跨文档逐字一致；六张主图 + 六张补充图重画；
三张主表 + 九张补充表重建；新增七项分析（填充率敏感性、块自助、层权重、两半逐基因可靠性分层、pooled 指数精确分解、五种汇总标量的四带稳定性、逐基因空间纹理与 Moran's I）；
参考文献 43 → 46 条（新增 Svensson 2020、Lähnemann 2020、Wu 2021），逐条核 DOI，**46/46 全部在正文有引用位置**；
投稿信重写并给出 3 位审稿人 + 2 位编委，压到 1 页。

**没做（需要你拍板）**：
1. 代码仓库地址（唯一硬缺口，见第二节）。
2. 第 4 套数据、影像平台交叉验证——结论不依赖平台，加了两周也改不了判断。
"""
    out = Path(outdir) / "00_投稿指南.md"
    out.write_text(txt, encoding="utf-8")
    return out


def count_words(doc_path):
    doc = Document(doc_path)
    txt = "\n".join(p.text for p in doc.paragraphs)
    return len(re.findall(r"[A-Za-z][A-Za-z'\-]*", txt))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=str(BUILD / "docs"))
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 先做一遍干跑，确定文献编号（正文顺序），再正式生成
    for _, t in BACKGROUND:
        R(t)
    for k, t in RESULTS:
        R(t)
    for _, t in DISCUSSION:
        R(t)
    for k, t in METHODS:
        R(k if k == "h2" else t)
    for k, v in ABSTRACT.items():
        R(v)
    for tag, t in LEGENDS:
        R(t)

    m = build_manuscript(outdir)
    c = build_cover(outdir)
    s = build_supplementary(outdir)
    g = build_guide(outdir)

    (BUILD / "numbers_ledger_v2.json").write_text(
        json.dumps(LEDGER, ensure_ascii=False, indent=1), encoding="utf-8")
    print("manuscript", m, count_words(m), "words total")
    print("cover", c)
    print("supp", s)
    print("guide", g)
    print("refs cited", len(CITE_ORDER), "of", len(REFS))
    print("uncited:", [k for k in REFS if k not in CITE_ORDER])
    print("ledger entries", len(LEDGER))


if __name__ == "__main__":
    main()

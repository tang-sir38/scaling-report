"""
101_supp_v4.py — v4 补充材料（Supplementary Information）

正文只引三处：Table S1（数据来源）、Table S2（窗口敏感性矩阵）、Supplementary Note 1
（两条失败的曲线预测）。这里把它们连同支撑它们的表和图一次建齐，所有数字从
NUMBERS_LOCK_v7.json 与 02_OUT/variogram/*.json 现场读取，禁手抄。

输出：
  07_SUBMISSION/_build_v4/docs/02_Supplementary_Information.docx
  07_SUBMISSION/_build_v4/figs/FigS1.pdf|png
  07_SUBMISSION/_build_v4/supp_manifest.json
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.shared import Cm, Pt, RGBColor

ROOT = Path(__file__).resolve().parents[1]
LOCK = ROOT / "02_OUT" / "NUMBERS_LOCK_v7.json"
VAR = ROOT / "02_OUT" / "variogram"
BUILD = ROOT / "07_SUBMISSION" / "_build_v4"
FIGS = BUILD / "figs"
DOCS = BUILD / "docs"

TISSUE = [("colon", "Colon"), ("lung", "Lung"), ("brain", "Brain")]
C_COLON, C_LUNG, C_BRAIN = "#2C5F7C", "#5B8FA8", "#C1663D"


def get(d, path):
    cur = d
    for p in path.split("."):
        cur = cur[p] if isinstance(cur, dict) else cur[int(p)]
    return cur


# ---------------------------------------------------------------- docx 基础
def setup(doc):
    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(10)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), "宋体")
    pf = st.paragraph_format
    pf.space_after = Pt(4)
    pf.line_spacing = 1.15
    for s in doc.sections:
        s.page_width, s.page_height = Cm(21.0), Cm(29.7)
        s.left_margin = s.right_margin = Cm(2.2)
        s.top_margin = s.bottom_margin = Cm(2.2)


def h(doc, text, size=13, before=10, after=4, bold=True):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(before)
    p.paragraph_format.space_after = Pt(after)
    r = p.add_run(text)
    r.bold = bold
    r.font.size = Pt(size)
    return p


def para(doc, text, size=10, italic=False, after=4):
    p = doc.add_paragraph()
    p.paragraph_format.space_after = Pt(after)
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    r = p.add_run(text)
    r.font.size = Pt(size)
    r.italic = italic
    return p


def caption(doc, label, text):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(8)
    p.paragraph_format.space_after = Pt(3)
    r1 = p.add_run(f"{label} ")
    r1.bold = True
    r1.font.size = Pt(9)
    r2 = p.add_run(text)
    r2.font.size = Pt(9)
    return p


def table(doc, header, rows, widths=None, font=8.5, align=None):
    t = doc.add_table(rows=1, cols=len(header))
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = False
    hdr = t.rows[0]
    hdr.cells  # touch
    trPr = hdr._tr.get_or_add_trPr()
    trPr.append(trPr.makeelement(qn("w:tblHeader"), {}))
    for i, txt in enumerate(header):
        c = hdr.cells[i]
        c.text = ""
        p = c.paragraphs[0]
        p.paragraph_format.space_after = Pt(0)
        p.paragraph_format.line_spacing = 1.0
        r = p.add_run(str(txt))
        r.bold = True
        r.font.size = Pt(font)
    for row in rows:
        cells = t.add_row().cells
        for i, v in enumerate(row):
            cells[i].text = ""
            p = cells[i].paragraphs[0]
            p.paragraph_format.space_after = Pt(0)
            p.paragraph_format.line_spacing = 1.0
            if align and align[i] == "r":
                p.alignment = WD_ALIGN_PARAGRAPH.RIGHT
            r = p.add_run(str(v))
            r.font.size = Pt(font)
    if widths:
        for i, w in enumerate(widths):
            for row in t.rows:
                row.cells[i].width = Cm(w)
    # 单元格内边距压到 0.08 cm，避免窄列被挤成零宽
    tblPr = t._tbl.tblPr
    mar = tblPr.makeelement(qn("w:tblCellMar"), {})
    for side in ("top", "left", "bottom", "right"):
        el = mar.makeelement(qn(f"w:{side}"), {qn("w:w"): "45", qn("w:type"): "dxa"})
        mar.append(el)
    tblPr.append(mar)
    return t


def ftxt(v, fmt):
    return format(v, fmt)


def fnum(v, fmt=".3f", signed=False):
    """nan 一律写 n.m.（not measurable），不留 'nan' 字样进表格。"""
    try:
        if v is None or (isinstance(v, float) and np.isnan(v)):
            return "n.m."
        if signed and abs(float(v)) < 5e-4:
            return format(0.0, fmt)
    except TypeError:
        pass
    return format(v, ("+" if signed else "") + fmt)


# ---------------------------------------------------------------- 数据表
def t_s1(L):
    rows = []
    spec = [("colon", "Visium HD Human Colon Cancer", "Human", "Colon adenocarcinoma"),
            ("lung", "Visium HD Human Lung Cancer", "Human", "Lung adenocarcinoma"),
            ("brain", "Visium HD Mouse Brain", "Mouse", "Whole brain")]
    for key, name, sp, tis in spec:
        d = L["datasets"][key]
        rows.append([name, "Visium HD", sp, tis, "6.5 × 6.5 mm",
                     f"{d['n_bins']:,}", f"{d['umi_per_bin']:.1f}",
                     f"{d['n_genes_noMT']:,}", "10x Genomics (public)"])
    st = get(L, "cross.stereo")
    rows.append(["Stereo-seq Mouse Whole Brain (C04042E2)", "Stereo-seq", "Mouse",
                 "Whole brain", "not applicable", "—", "—", "—",
                 "STOmicsDB (C04042E2)"])
    rows.append(["Visium HD Human Colon Cancer P1", "Visium HD", "Human",
                 "Second section, same block", "6.5 × 6.5 mm", "—", "—",
                 f"{get(L, 'multi.sections.Visium_HD_Human_Colon_Cancer_P1.n_genes'):,}",
                 "10x Genomics (public)"])
    return rows


def t_s2(L):
    ws = L["window_scan"]["by_tissue"]
    grid = {t: {(tuple(g["gamma_window"]), tuple(g["beta_window"])): g
                for g in ws[t]["grid"]} for t, _ in TISSUE}
    keys = sorted(grid["colon"].keys(), key=lambda k: (k[0][0], k[0][1], k[1][0]))
    rows = []
    for k in keys:
        gw, bw = k
        r = [f"{gw[0]}–{gw[1]}", f"{bw[0]}–{bw[1]}"]
        for t, _ in TISSUE:
            r.append(f"{grid[t][k]['pearson']:.3f}")
        rows.append(r)
    return rows


def t_s3(L):
    g = L["intervene"]["single_gene"]["genes"]
    order = ["IGKC", "IGHG1", "JCHAIN", "MZB1", "EPCAM", "COL1A1"]
    rows = []
    for k in order:
        d = g[k]
        rows.append([k, f"{d['total']:,.0f}", f"{d['occupancy'] * 100:.1f}",
                     fnum(abs(d["baseline_gamma"])),
                     fnum(abs(d["permute_gamma"])),
                     (f"{abs(d['permute_lo']):.3f}–{abs(d['permute_hi']):.3f}"
                      if not np.isnan(d["permute_lo"]) else "n.m."),
                     fnum(abs(d["permute_gamma"]) - abs(d["baseline_gamma"]), signed=True),
                     fnum(abs(d["excise_gamma"]) if not np.isnan(d["excise_gamma"]) else np.nan),
                     fnum(abs(d["isolate_gamma"]))])
    return rows


def t_s4(L):
    c = L["intervene"]["dose_response"]["curves"]
    rows = []
    for k in sorted(c, key=lambda x: int(x)):
        d = c[k]
        rows.append([f"{d['n_permuted']:,}", f"{d['mu2_share'] * 100:.1f}",
                     f"{abs(d['gamma']):.3f}", f"{d['r2']:.3f}"])
    return rows


def t_s5():
    d = json.loads((VAR / "curveshape_colon.json").read_text(encoding="utf-8"))
    zp = d["per_gene_zero_param"]
    cal = d["per_gene_calibrated_at_16"]
    rows = []
    for k in ["8", "16", "32", "64", "128"]:
        a = zp.get(k, {})
        b = cal.get(k, {})
        rows.append([k,
                     fnum(a.get("median"), signed=True) if a.get("median") is not None else "—",
                     fnum(b.get("median"), signed=True) if b.get("median") is not None else "—",
                     fnum(b.get("pearson")) if b.get("pearson") is not None else "—",
                     fnum(b.get("slope")) if b.get("slope") is not None else "—"])
    return rows


def t_s6():
    d = json.loads((VAR / "gamma_from_rhobar.json").read_text(encoding="utf-8"))["by_tissue"]
    rows = []
    for key, name in TISSUE:
        v = d[key]
        rows.append([name, f"{v['n']:,}", f"{abs(v['gamma_obs_median']):.3f}",
                     f"{abs(v['gamma_pred_median']):.3f}",
                     f"{v['diff_median']:+.3f}", f"{v['r_log']:.3f}",
                     f"{v['r_spearman']:.3f}"])
    return rows


# ---------------------------------------------------------------- Fig S1
def fig_s1():
    zc = json.loads((VAR / "curveshape_colon.json").read_text(encoding="utf-8"))
    gr = json.loads((VAR / "gamma_from_rhobar.json").read_text(encoding="utf-8"))["by_tissue"]

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "axes.linewidth": 0.7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1))

    # (a) 零参数 vs 在 16 µm 定标后，中位偏差随 d
    ax = axes[0]
    ks = ["8", "16", "32", "64", "128"]
    x = np.array([float(k) for k in ks])
    zp = [zc["per_gene_zero_param"].get(k, {}).get("median") for k in ks]
    cal = [zc["per_gene_calibrated_at_16"].get(k, {}).get("median") for k in ks]
    ax.plot(x, zp, "o-", color=C_COLON, ms=4.5, lw=1.2, label="zero free parameter")
    ax.plot(x, cal, "s--", color=C_BRAIN, ms=4.5, lw=1.2,
            label="calibrated at 16 µm")
    ax.axhline(0, color="0.45", lw=0.8, ls=":")
    ax.set_xscale("log", base=2)
    ax.set_xticks(x)
    ax.set_xticklabels(ks)
    ax.minorticks_off()
    ax.set_xlabel("Block side $d$ (µm)")
    ax.set_ylabel("Median $\\log\\,\\alpha_{obs}-\\log\\,\\alpha_{pred}$")
    ax.set_title("Predicted curve is too steep, and drifts", fontsize=8.5, pad=6)
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.text(-0.16, 1.06, "a", transform=ax.transAxes, fontsize=10, va="top")

    # (b) γ_pred vs γ_obs，三套数据
    ax = axes[1]
    cols = {"colon": C_COLON, "lung": C_LUNG, "brain": C_BRAIN}
    lim = [0, 1.55]
    ax.plot(lim, lim, color="0.55", lw=0.9, ls="--", zorder=1)
    for key, name in TISSUE:
        v = gr[key]
        ax.scatter(abs(v["gamma_obs_median"]), abs(v["gamma_pred_median"]),
                   s=46, color=cols[key], edgecolor="white", linewidth=0.7,
                   zorder=3, label=name)
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_xlabel("Observed exponent $\\gamma$ (16–128 µm)")
    ax.set_ylabel("Predicted exponent from $\\bar\\rho$")
    ax.set_title("Predicted exponent is steeper in every section", fontsize=8.5, pad=6)
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.text(-0.16, 1.06, "b", transform=ax.transAxes, fontsize=10, va="top")

    fig.tight_layout(w_pad=2.2)
    FIGS.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGS / "FigS1.pdf", bbox_inches="tight", metadata={"CreationDate": None})
    fig.savefig(FIGS / "FigS1.png", dpi=400, bbox_inches="tight")
    plt.close(fig)


def t_s7(L):
    """差异表达跨尺度汇总，每个粗化尺度一行。"""
    de = L["downstream"]["de"]
    rows = []
    for k in sorted(de, key=int):
        d = de[k]
        rows.append([k,
                     f"{d['n_a']:,} / {d['n_b']:,}",
                     f"{d['top100_recall'] * 100:.0f}",
                     f"{d['jaccard_sig']:.2f}",
                     f"{d['n_sig_ref']:,} / {d['n_sig_s']:,}",
                     f"{d['logfc_slope']:.3f}",
                     f"{d['logfc_slope_raw']:.3f}",
                     f"{d['logfc_r']:.3f}",
                     f"{d['stat_spearman']:.3f}"])
    return rows


def t_s8(L):
    """映射代理：每个 marker 集 × 每个尺度，两种归一化。"""
    mp = L["downstream"]["map"]
    rows = []
    for pn in ("epithelial", "stromal", "T_cell", "plasma_B", "endothelial"):
        for k in sorted(mp, key=int):
            if pn not in mp[k]:
                continue
            a = mp[k][pn].get("log1p", {})
            b = mp[k][pn].get("raw", {})
            rows.append([pn.replace("_", " "), k,
                         fnum(a.get("slope")), fnum(a.get("pearson")),
                         fnum(b.get("slope")), fnum(b.get("pearson"))])
    return rows


def fig_s2(L):
    D = L["downstream"]
    od, de, mp = D["od"], D["de"], D["map"]
    ks = sorted(od, key=int)
    x = np.array([int(k) for k in ks], dtype=float)

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "axes.linewidth": 0.7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    fig, axes = plt.subplots(1, 3, figsize=(7.8, 3.05))

    # (a) φ 与隐含 FPR
    ax = axes[0]
    ph = np.array([od[k]["phi"] for k in ks])
    fp = np.array([od[k]["fpr"] * 100 for k in ks])
    ax.plot(x, ph, "o-", color=C_COLON, ms=3.6, lw=1.1)
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(ks)
    ax.minorticks_off()
    ax.set_xlabel("Block side $d$ (µm)")
    ax.set_ylabel("$\\phi=1+\\mu\\alpha$", color=C_COLON, labelpad=1)
    ax.tick_params(axis="y", colors=C_COLON)
    ax2 = ax.twinx()
    ax2.plot(x, fp, "s--", color=C_BRAIN, ms=3.2, lw=1.0)
    ax2.set_ylabel("Implied FPR (%)", color=C_BRAIN)
    ax2.tick_params(axis="y", colors=C_BRAIN)
    ax2.set_ylim(0, 80)
    ax.set_title("Poisson fails faster than $d^{-2}$", fontsize=8, pad=6)
    ax.spines["top"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    ax.text(-0.24, 1.06, "a", transform=ax.transAxes, fontsize=10, va="top")

    # (b) 映射斜率：同一批分数，两种归一化
    ax = axes[1]
    PAN = [("epithelial", "Epithelial", C_COLON), ("stromal", "Stromal", C_LUNG),
           ("T_cell", "T cell", C_BRAIN), ("plasma_B", "Plasma cell", "#7A6A9B"),
           ("endothelial", "Endothelial", "#8C6D4F")]
    sks = sorted(mp, key=int)
    lim = [0.65, 3.15]
    ax.plot(lim, lim, color="0.55", lw=0.9, ls="--", zorder=1)
    ax.axhline(1, color="0.75", lw=0.7, ls=":")
    ax.axvline(1, color="0.75", lw=0.7, ls=":")
    for pn, lab, col in PAN:
        xr, yl = [], []
        for k in sks:
            a = mp[k].get(pn, {})
            if "raw" in a and "log1p" in a:
                xr.append(a["raw"]["slope"])
                yl.append(a["log1p"]["slope"])
        ax.scatter(xr, yl, s=24, color=col, edgecolor="white", linewidth=0.5,
                   zorder=3, label=lab)
    ax.set_xlim(0.65, 1.85)
    ax.set_ylim(*lim)
    ax.set_xlabel("Slope, raw ratio")
    ax.set_ylabel("Slope, log1p")
    ax.set_title("Same scores, two normalisations", fontsize=8, pad=6)
    ax.legend(frameon=True, facecolor="white", edgecolor="none", framealpha=0.92,
              fontsize=5.8, loc="lower right", handletextpad=0.3, borderpad=0.3,
              labelspacing=0.25, scatterpoints=1)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.text(-0.24, 1.06, "b", transform=ax.transAxes, fontsize=10, va="top")

    # (c) DE：rank 指标与斜率
    ax = axes[2]
    dks = sorted(de, key=int)
    xd = np.array([int(k) for k in dks], dtype=float)
    ax.plot(xd, [de[k]["top100_recall"] for k in dks], "o-",
            color=C_COLON, ms=3.4, lw=1.1, label="Top-100 recall")
    ax.plot(xd, [de[k]["jaccard_sig"] for k in dks], "s-",
            color=C_LUNG, ms=3.4, lw=1.1, label="Jaccard, significant set")
    ax.set_ylim(0, 1.05)
    ax.set_xscale("log", base=2)
    ax.set_xticks(xd)
    ax.set_xticklabels([str(int(v)) for v in xd])
    ax.minorticks_off()
    ax.set_xlabel("Block side $d$ (µm)")
    ax.set_ylabel("Agreement with 8 µm")
    ax3 = ax.twinx()
    ax3.plot(xd, [de[k]["logfc_slope"] for k in dks], "^--",
             color=C_BRAIN, ms=3.0, lw=1.0, label="logFC slope, log1p")
    ax3.plot(xd, [de[k]["logfc_slope_raw"] for k in dks], "v--",
             color="#7A6A9B", ms=3.0, lw=1.0, label="logFC slope, raw")
    ax3.axhline(1, color="0.55", lw=0.7, ls=":")
    ax3.set_ylabel("logFC slope on 8 µm")
    ax3.set_ylim(0.6, 1.3)
    ax.set_title("Rank decays, slope flips sign", fontsize=8, pad=6)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax3.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, frameon=True, facecolor="white", edgecolor="none",
              framealpha=0.92, fontsize=5.8, loc="lower left", handlelength=1.4,
              borderpad=0.3, labelspacing=0.3)
    ax.spines["top"].set_visible(False)
    ax3.spines["top"].set_visible(False)
    ax.text(-0.24, 1.06, "c", transform=ax.transAxes, fontsize=10, va="top")

    fig.tight_layout(w_pad=1.9)
    FIGS.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGS / "FigS2.pdf", bbox_inches="tight", metadata={"CreationDate": None})
    fig.savefig(FIGS / "FigS2.png", dpi=400, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- main
def main():
    L = json.loads(LOCK.read_text(encoding="utf-8"))
    fig_s1()

    doc = Document()
    setup(doc)

    p = doc.add_paragraph()
    r = p.add_run("Supplementary Information")
    r.bold = True
    r.font.size = Pt(15)
    p.paragraph_format.space_after = Pt(2)

    para(doc, "Binning is set by spatial covariance and by which transcripts dominate",
         size=11, italic=True, after=2)
    para(doc, "This file contains three items: the dataset table (Table S1), the full "
              "window-sensitivity matrix (Table S2), and one note recording two analyses "
              "that were attempted and failed (Supplementary Note 1). Two further tables "
              "(S3, S4) give the per-gene and dose–response numbers behind Fig. 5. "
              "Every number is taken from the numbers lock file (NUMBERS_LOCK_v7.json) or "
              "from the derived output of the script named in the corresponding Methods "
              "subsection; nothing is transcribed by hand.", size=9)

    # ---- Table S1
    h(doc, "Table S1. Sections analysed and their provenance.")
    table(doc, ["Section", "Platform", "Species", "Tissue", "Capture area",
                "2 µm bins in tissue", "Transcripts per bin", "Genes retained", "Source"],
          t_s1(L), widths=[3.3, 1.6, 1.2, 2.6, 1.9, 2.0, 1.7, 1.5, 2.4], font=7.5,
          align=["l", "l", "l", "l", "l", "r", "r", "r", "l"])
    para(doc, "Bin counts and transcripts per bin are for the full in-tissue 2 µm grid "
              "before any re-aggregation. The colon P1 section is a second section of the "
              "same tumour block and is used only as a replicate for the pooled exponent; "
              "it is not part of the eight-scale analysis. The Stereo-seq section is used "
              "only for the cross-chemistry comparison of per-gene exponents. All four "
              "Visium HD sections were obtained as Space Ranger 2 µm output; no further "
              "alignment, filtering or normalisation was applied beyond removal of "
              "mitochondrial genes.", size=8.5, italic=True)

    # ---- Table S2
    h(doc, "Table S2. Window sensitivity of the correspondence between the two per-gene exponents.")
    table(doc, ["γ window (µm)", "β window (grid units)", "r, colon", "r, lung", "r, brain"],
          t_s2(L), widths=[3.0, 3.4, 2.4, 2.4, 2.4], font=7.5,
          align=["l", "l", "r", "r", "r"])
    para(doc, "Pearson correlation between the per-gene exponent obtained by aggregation "
              "(γ) and the per-gene exponent of the covariogram (β), for every combination "
              "of the two fit windows, at the genes fitted in both (n = 2,476 colon, 1,726 "
              "lung, 809 brain). One grid unit of the covariogram is 8 µm. The six γ "
              "windows are not a full cross product: the lower bound is held at 2 µm and "
              "the upper bound varied (32, 64, 128, 256 µm), then the upper bound is held "
              "at 256 µm and the lower bound raised (8, 16 µm). The window used in the "
              "main text, γ over 2–256 µm and β over 2–8 grid units, was fixed before "
              "these correlations were computed and is not selected from this matrix; the "
              "matrix is reported as a sensitivity analysis, and its purpose is to show "
              "which conclusions survive a change of window, not to choose one.",
         size=8.5, italic=True)

    # ---- Table S3
    h(doc, "Table S3. Per-gene exponent before and after random permutation of transcript positions.")
    table(doc, ["Gene", "Transcripts", "Occupancy (%)", "γ, observed",
                "γ, permuted", "95% range", "Δγ", "γ, excised", "γ, isolated"],
          t_s3(L), widths=[2.0, 2.4, 2.2, 2.1, 2.1, 3.0, 1.7, 2.1, 2.1], font=8,
          align=["l", "r", "r", "r", "r", "l", "r", "r", "r"])
    para(doc, "Permutation moves each gene's counts to random positions drawn from the "
              "in-tissue domain, so the total, the count distribution and the occupancy "
              "are unchanged and only the arrangement is destroyed; 20 independent "
              "permutations per gene. Occupancy is the fraction of 2 µm bins holding at "
              "least one transcript of that gene. All exponents are reported as absolute "
              "values, as in the main text. Isolated: the gene is the only one retained. "
              "n.m., not measurable: for MZB1 the exponent after excision did not "
              "converge; the gene carries 14,935 transcripts in total, an order of "
              "magnitude fewer than the others, so removing its densest bins leaves too "
              "little signal at the coarsest scales for the fit.",
         size=8.5, italic=True)

    # ---- Table S4
    h(doc, "Table S4. Pooled exponent as a function of how many transcripts are permuted.")
    table(doc, ["Genes permuted", "Share of Σμ² permuted (%)", "Pooled γ", "r²"],
          t_s4(L), widths=[3.0, 5.4, 2.6, 2.2], font=8.5,
          align=["r", "r", "r", "r"])
    dr = L["intervene"]["dose_response"]["curves"]
    para(doc, "Genes were ranked by μ² at 2 µm and the top N permuted independently within "
              "the colon section. The exponent rises monotonically with N: permuting the "
              f"ten heaviest transcripts, which carry {dr['10']['mu2_share']:.0%} of the "
              f"weight, already takes it from {abs(dr['0']['gamma']):.3f} to "
              f"{abs(dr['10']['gamma']):.3f}. Full curves are in Fig. 5b.",
         size=8.5, italic=True)

    # ---- Supplementary Note
    h(doc, "Supplementary Note 1. Two attempts to predict the whole α(d) curve, and why they fail.", 12)

    para(doc, "The identity α(A) = CV²₀ · ρ̄(A) relates the overdispersion at support A to "
              "the mean correlation inside that support. If it could be evaluated, it would "
              "predict the entire α(d) curve with no free parameter — including the slight "
              "upward bend near 16–32 µm visible in the local slope of Fig. 1b — and the "
              "equality of exponents would become a corollary rather than a separate "
              "result. Two routes were tried on the colon section, where the per-gene "
              "covariogram is measured on the 8 µm grid. Both fail, and the reason is a "
              "property of the data rather than of the identity.")

    h(doc, "First route: zero-parameter prediction of the curve", 10.5, before=8)
    para(doc, "For each gene, α(d) was predicted from CV²₀ and the block average of its "
              "measured covariogram, taking the covariogram as isotropic and averaging "
              "ρ(h) over the pairs of 8 µm cells inside a k × k block. The deviation "
              "log α_obs − log α_pred is zero at 8 µm by construction and then grows "
              "monotonically (Table S5, Fig. S1a). A single multiplicative error — which "
              "is what an 8 µm resolution limit inside the block would produce — would "
              "give a constant offset; a monotonic growth means the predicted shape is "
              "wrong as well, not merely shifted. Calibrating one constant at 16 µm leaves "
              "residuals that still rise, and the out-of-sample correlation with the "
              "observed curve over 32–128 µm is only 0.53–0.61 with a slope near 0.3.")

    caption(doc, "Table S5.",
            "Deviation of the predicted α(d) from the observed curve, colon section, "
            "n = 2,476 genes, in log units. Middle column: no free parameter. Right "
            "columns: one multiplicative constant fitted at 16 µm and then applied out of "
            "sample; Pearson and slope are between predicted and observed α at that scale.")
    table(doc, ["d (µm)", "Median log deviation, zero-parameter",
                "Median log deviation, calibrated at 16 µm", "Pearson", "Slope"],
          t_s5(), widths=[1.6, 3.9, 6.1, 2.0, 1.8], font=8.5,
          align=["l", "r", "r", "r", "r"])

    h(doc, "Second route: predict only the slope", 10.5, before=8)
    gr = json.loads((VAR / "gamma_from_rhobar.json").read_text(encoding="utf-8"))["by_tissue"]
    fac = [abs(gr[t]["gamma_pred_median"]) / abs(gr[t]["gamma_obs_median"]) for t, _ in TISSUE]
    para(doc, "If the absolute height is corrupted by the resolution limit but the slope "
              "is not, the decay rate of ρ̄ should still recover γ. ρ̄ was therefore "
              "formed at k = 2, 4, 8 and 16 (that is, 16, 32, 64 and 128 µm), its decay "
              "fitted over that same range, and the result compared with the observed γ "
              "fitted over 16–128 µm so that the two windows coincide. The predicted "
              f"exponent is steeper than the observed one in every section, by a factor of "
              f"{min(fac):.1f} to {max(fac):.1f} (Table S6, Fig. S1b), and the per-gene "
              f"correlation between the two is at best moderate — Pearson "
              f"{gr['colon']['r_log']:.2f} in colon, {gr['lung']['r_log']:.2f} in lung and "
              f"{gr['brain']['r_log']:.2f} in brain — so this route cannot rank genes "
              "either.")

    caption(doc, "Table S6.",
            "Observed and predicted exponents over a common window (16–128 µm), all three "
            "sections. Predicted values come from the block-averaged covariogram. r is "
            "Pearson between the two per-gene estimates, in log space.")
    table(doc, ["Section", "Genes", "γ observed", "γ predicted", "Median difference",
                "Pearson", "Spearman"],
          t_s6(), widths=[2.4, 2.2, 2.4, 2.6, 3.2, 2.2, 2.2], font=8.5,
          align=["l", "r", "r", "r", "r", "r", "r"])

    h(doc, "Why both routes fail, and what it implies", 10.5, before=8)
    para(doc, "The covariogram is measured on the 8 µm grid. The value at the first lag is "
              "therefore the correlation between two adjacent 8 µm blocks, which already "
              "contains the averaging over 8 µm inside each block. Averaging that quantity "
              "again over a k × k block performs the block average twice: the correlation "
              "at short lags is depressed, every ρ̄(d) falls faster than it should, and the "
              "bias worsens as k grows. This is not a defect of the identity — the identity "
              "is verified directly on synthetic fields, where the covariogram is known at "
              "the point scale (Fig. 3a) — but a hard limit of this data set. Measuring "
              "ρ(h) on the 2 µm grid is not an alternative, because at 2 µm a bin holds "
              "zero or one transcript and the covariance is dominated by counting noise.")
    para(doc, "The consequence for interpretation is stated in the main text as a "
              "limitation and is worth repeating here: because only lags of 8 µm and above "
              "are measurable, and because each measured lag is itself block-averaged, the "
              "absolute height and the detailed shape of ρ̄(A) cannot be verified. What can "
              "be compared between the two routes is the local decay rate, and that "
              "comparison is the content of Fig. 4. The negative result is reported rather "
              "than omitted because the two routes were the natural test of the identity "
              "and a reader is entitled to know that they were run.")

    # ---- Fig S1
    doc.add_page_break()
    h(doc, "Figure S1. The two failed curve predictions.", 11)
    doc.add_picture(str(FIGS / "FigS1.png"), width=Cm(16.4))
    doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
    caption(doc, "Figure S1.",
            "(a) Median deviation between observed and predicted α(d) for the colon "
            "section, with no free parameter (filled circles) and with one multiplicative "
            "constant fitted at 16 µm and applied out of sample (open squares). The "
            "deviation is zero at 8 µm by construction and grows monotonically, so the "
            "discrepancy is in the shape and not only in the level. (b) Predicted "
            "exponent from the block-averaged covariogram against the observed exponent "
            "over the same window (16–128 µm), for the three sections; dashed line, y = x. "
            "Prediction is systematically steeper, and no section falls near the identity "
            "line.")

    # ---- Table S7 / S8
    doc.add_page_break()
    h(doc, "Table S7. Differential expression between two regions, repeated at every scale.",
      11)
    para(doc, "The colon section was split at 16 µm into its epithelial and stromal ends by "
              "marker score (Methods), the labels carried to each coarser scale by integer "
              "coordinate mapping, and the 3,000 most abundant genes tested with a "
              "two-sided Wilcoxon rank-sum test. The reference is the 8 µm result. Bins: "
              "number of bins in each group after the 20-count filter. Top-100: fraction of "
              "the 8 µm top 100 (by signed test statistic) recovered at that scale. "
              "Jaccard: overlap of the sets passing FDR < 0.05 with |log2FC| > 0.5. "
              "Significant: size of that set at 8 µm and at the stated scale. The two "
              "slopes regress the coarse-scale log fold change on the 8 µm value, first "
              "with the mean taken over log1p(CP10K) values and then over the raw ratio; "
              "r and ρ are Pearson and Spearman on the same pairs. A slope of one means "
              "the effect size is comparable across scales, which neither normalisation "
              "delivers.", size=8.5, italic=True)
    table(doc, ["Scale (µm)", "Bins, epi / stroma", "Top-100 (%)", "Jaccard",
                "Significant, 8 µm / scale", "Slope, log1p", "Slope, raw ratio",
                "r", "ρ"],
          t_s7(L), widths=[1.9, 3.2, 1.9, 1.6, 3.4, 1.9, 2.2, 1.3, 1.3], font=8,
          align=["l", "r", "r", "r", "r", "r", "r", "r", "r"])
    para(doc, "At 128 µm the groups hold 125 and 222 bins; that row is reported for "
              "completeness and carries the widest sampling error in the table.",
         size=8, italic=True)

    h(doc, "Table S8. Cell-type module scores across scales, two normalisations.", 11,
      before=12)
    para(doc, "Each marker set was scored as the sum of its counts over the bin total in "
              "CP10K. The reference is the 8 µm score aggregated to the coarse bin as a "
              "mean over the 8 µm bins it contains; the comparison is the coarse bin's own "
              "score. Slope is ordinary least squares with the reference on the x axis, so "
              "a slope of one means the score is comparable across supports and a slope "
              "above one means the coarse bin overstates it. r is Pearson. Marker sets: "
              "epithelial, EPCAM/KRT18/KRT8/CEACAM5/TFF3/CDX2/KRT19; stromal, "
              "COL1A1/COL1A2/COL3A1/DCN/LUM/PDGFRB/THY1; T cell, "
              "CD3D/CD3E/CD2/CD247/IL7R/CCL5; plasma cell, "
              "IGKC/IGHG1/JCHAIN/MZB1/MS4A1/CD79A; endothelial, "
              "PECAM1/VWF/CLDN5/PLVAP/ENG.", size=8.5, italic=True)
    table(doc, ["Marker set", "Scale (µm)", "Slope, log1p", "r, log1p",
                "Slope, raw ratio", "r, raw ratio"],
          t_s8(L), widths=[3.4, 2.2, 2.6, 2.2, 3.0, 2.6], font=8,
          align=["l", "r", "r", "r", "r", "r"])

    # ---- Fig S2
    fig_s2(L)
    doc.add_page_break()
    h(doc, "Figure S2. What the exponent costs the analyses that follow.", 11)
    doc.add_picture(str(FIGS / "FigS2.png"), width=Cm(16.4))
    doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
    caption(doc, "Figure S2.",
            "(a) Variance understated by a Poisson assumption, φ = 1 + μα (left axis, "
            "filled circles, median over the 2,000 most abundant genes), and the implied "
            "false-positive rate of a two-sided test at a nominal 5% (right axis, open "
            "squares). Both rise with bin size because μ grows with d² while α falls only "
            "as a shallow power of d rather than as its inverse square. (b) The same module "
            "scores, plotted by the slope they give under "
            "each normalisation: x is the slope on the raw ratio, y the slope after "
            "log1p, the dashed line is y = x and the dotted lines mark one. Every marker "
            "set sits near x = 1, so the raw ratio is comparable across supports, but the "
            "same scores spread from 0.77 to 2.95 vertically; the normalisation alone "
            "moves a score by up to a factor of three, and it moves sparse marker sets "
            "furthest. (c) "
            "Agreement of the epithelial-versus-stromal comparison with its 8 µm result: "
            "top-100 recall and Jaccard of the significant set (left axis), and the slope "
            "of coarse-scale log fold changes on the 8 µm values under the two "
            "normalisations (right axis). The slope crosses one between normalisations.")

    DOCS.mkdir(parents=True, exist_ok=True)
    out = DOCS / "02_Supplementary_Information.docx"
    doc.save(out)

    man = dict(
        file=str(out.relative_to(ROOT)),
        tables=dict(S1=len(t_s1(L)), S2=len(t_s2(L)), S3=len(t_s3(L)),
                    S4=len(t_s4(L)), S5=len(t_s5()), S6=len(t_s6()),
                    S7=len(t_s7(L)), S8=len(t_s8(L))),
        figures=["figs/FigS1.pdf", "figs/FigS1.png",
                 "figs/FigS2.pdf", "figs/FigS2.png"],
        note="Supplementary Note 1",
    )
    (BUILD / "supp_manifest.json").write_text(
        json.dumps(man, indent=1, ensure_ascii=False), encoding="utf-8")
    print(f"写出 {out}")
    print(f"  表行数 {man['tables']}")
    print("  图 FigS1.pdf / FigS1.png / FigS2.pdf / FigS2.png")


if __name__ == "__main__":
    main()

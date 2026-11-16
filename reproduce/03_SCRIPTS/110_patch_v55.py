"""Build v5.5: the exponent range widens to 5- to 15-fold when four further
sections are added, so the title, the abstract, the Results and the Discussion
have to move together.

Why this is a separate script rather than an edit
-------------------------------------------------
The manuscript has no markdown source (the v5 text exists only as a .docx), so
the only disciplined way to change it is to re-run a patcher that (a) starts
from the delivered v5.4 bytes, (b) refuses to act if any anchor has moved, and
(c) recomputes every quoted number from `extra_tissues_20261116.json` — which
`120_aggregate_extras.py` writes from the raw 2 µm matrices — rather than from
anything typed here.

Budget
------
v5.4 body was 3,424 of 3,500 words and the abstract 148 of 150. Four new
sections have to fit in the remaining 76 body words, so the abstract edit is a
net *reduction* (three per-section exponents collapse into one range) and the
new Results sentences are paid for out of that slack. Both limits are asserted
after the edit and the remainder is printed.

Usage
-----
  python 110_patch_v55.py --dry-run    # measure, write nothing
  python 110_patch_v55.py              # write the three docx into the v5.5 package
"""
from __future__ import annotations

import argparse
import copy as _copy
import hashlib
import json
import re
import sys
from datetime import datetime
from pathlib import Path

import docx
from docx.shared import Cm
from docx.table import _Row
from docx.text.paragraph import Paragraph

FIXED_TS = datetime(2026, 11, 16, 0, 0, 0)
MINUS = "\u2212"
MICRO = "\u00b5"
ALPHA = "\u03b1"
PROP = "\u221d"
ENDASH = "\u2013"
EMDASH = "\u2014"

REPRO = Path(r"C:\Users\唐朝华\Desktop\scaling-report\reproduce")
SRC = Path(r"C:\Users\唐朝华\Desktop\23空间分辨率\A题_空间分辨率_投稿包_v5.4_20261116")
DST = Path(r"C:\Users\唐朝华\Desktop\23空间分辨率\A题_空间分辨率_投稿包_v5.5_20261116")
BUILD = REPRO / "07_SUBMISSION" / "_build_v55"
AUDIT = BUILD / "audit"
LOCK = REPRO / "02_OUT" / "NUMBERS_LOCK_v7.json"
SUMMARY = REPRO / "02_OUT" / "extra_tissues_20261116.json"

MS = "02_主稿/01_Manuscript.docx"
SI = "03_补充材料/02_Supplementary_Information.docx"
CL = "01_投稿信/00_Cover_Letter.docx"

BASE_SPECIES = {"colon": "Human", "lung": "Human", "brain": "Mouse"}


def freeze_zip(path: Path, when=(2026, 11, 16, 0, 0, 0)) -> None:
    """A .docx is a zip; unpinned entry stamps make two runs differ byte-wise."""
    import io
    import zipfile

    with zipfile.ZipFile(io.BytesIO(path.read_bytes())) as zin:
        items = [(i.filename, zin.read(i.filename)) for i in zin.infolist()]
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zout:
        for name, blob in items:
            zi = zipfile.ZipInfo(name, date_time=when)
            zi.compress_type = zipfile.ZIP_DEFLATED
            zi.external_attr = 0o644 << 16
            zout.writestr(zi, blob)
    path.write_bytes(buf.getvalue())


def wc(text: str) -> int:
    """The journal counts words by whitespace, so every check here does too."""
    return len(text.split())


class Num:
    """Every quoted value, formatted at the point of use, never typed in."""

    def __init__(self, summary: dict):
        self.s = summary
        self.sec = {x["key"]: x for x in summary["sections"]}
        self.rng = summary["range"]
        self.primary = [k for k in self.sec if k != "colon_p1"]

    def drop(self, key: str) -> float:
        return self.sec[key]["drop_fold"]

    def a2(self, key: str) -> float:
        return self.sec[key]["alpha_2um"]

    def umi(self, key: str) -> float:
        return self.sec[key]["umi_per_bin"]

    def ci(self, key: str) -> tuple[float, float]:
        b = self.sec[key].get("boot")
        if not b:
            raise SystemExit(f"{key} 没有 bootstrap 结果，不能报 CI")
        return b["ci95"][0], b["ci95"][1]

    def species(self, key: str) -> str:
        return self.sec[key].get("species") or BASE_SPECIES[key]

    def gamma_abs(self) -> tuple[float, float]:
        g = [abs(self.sec[k]["gamma"]) for k in self.primary]
        return min(g), max(g)

    def genes(self, species: str) -> int:
        vals = {self.sec[k]["n_genes_noMT"] for k in self.sec
                if self.species(k) == species}
        if len(vals) != 1:
            detail = {k: self.sec[k]["n_genes_noMT"] for k in self.sec
                      if self.species(k) == species}
            raise SystemExit(f"{species} 各片基因数不一致，Methods 那句话要重写: {detail}")
        return vals.pop()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--summary", default=str(SUMMARY))
    ap.add_argument("--dst", default=None)
    a = ap.parse_args()

    summary = json.loads(Path(a.summary).read_text(encoding="utf-8"))
    lock = json.loads(LOCK.read_text(encoding="utf-8"))
    n = Num(summary)
    rng = summary["range"]
    floor_max = rng["floor_max"]

    assert rng["floor_max"] != 16, (
        "floor 上界撞上 16 —— 与 16,384 在版面上混淆，标题方案需要单独拍板")

    si_k, ki_k, pa_k, em_k = ("mouse_small_intestine", "mouse_kidney",
                              "human_pancreas", "mouse_embryo")
    new_keys = (si_k, ki_k, pa_k, em_k)
    longest = max(new_keys, key=n.drop)
    lo_ci, hi_ci = n.ci(longest)
    colon_drop = n.drop("colon")
    ci_gap = colon_drop - lo_ci

    g_lo, g_hi = n.gamma_abs()
    floor_min = rng["floor_min"]
    WORDS = {2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven",
             8: "eight", 9: "nine", 10: "ten", 11: "eleven", 12: "twelve",
             13: "thirteen", 14: "fourteen", 15: "fifteen", 16: "sixteen",
             17: "seventeen", 18: "eighteen", 19: "nineteen", 20: "twenty"}
    w_min = WORDS.get(floor_min, str(floor_min))
    w_max = WORDS.get(floor_max, str(floor_max))

    # ---- texts -------------------------------------------------------------
    def tissue_phrase(k: str) -> str:
        """Lower-cased label minus the platform prefix, for mid-sentence use.

        Labels are Title Case ("Visium HD Mouse Kidney"), so lowering only the
        first character leaves "mouse Kidney" mid-sentence. No label here
        contains a proper noun, so the whole phrase is lowered; if one ever does,
        this is the line that has to change.
        """
        return n.sec[k]["label"].replace("Visium HD ", "").lower()

    title_v54 = ("Overdispersion in binned spatial transcriptomics falls 5- to "
                 "14-fold, not 16,384-fold, because it is set by tissue covariance")
    title_v55 = ("Overdispersion in binned spatial transcriptomics falls "
                 f"{floor_min}- to {floor_max}-fold, not 16,384-fold, because it is "
                 "set by tissue covariance")

    abs_old_mid = ("of human colon, human lung and mouse brain Visium HD sections "
                   "by exact summation and measure " + ALPHA + f" {PROP} d")
    abs_new_mid = ("of seven Visium HD sections from human and mouse by exact "
                   "summation and measure falls of "
                   f"{rng['min']:.1f}- to {rng['max']:.1f}-fold, with "
                   f"{ALPHA} {PROP} d")

    abs_old_drop = f" {EMDASH} falls of 14.3-, 13.8- and 5.8-fold."
    abs_old_tail = " (Pearson 0.78, 0.73 and 0.54)."
    abs_new_tail = " in three of them (Pearson 0.78, 0.73 and 0.54)."

    intro_old = "in three Visium HD sections " + EMDASH + " human colon cancer, human lung cancer and mouse brain " + EMDASH + " and identify what it measures."
    intro_new = ("in three Visium HD sections " + EMDASH + " human colon cancer, "
                 "human lung cancer and mouse brain " + EMDASH + " plus four further "
                 "sections (normal pancreas, small intestine, kidney, embryo), and "
                 "identify what it measures.")

    # "at most" is not decoration: once one section falls only 3.5-fold, a bare
    # "falls fifteen-fold" would be read as a claim about every section. The
    # heading states the bound, and the bound is what the data support.
    head_old = "Overdispersion falls fourteen-fold, not sixteen thousand-fold"
    head_new = (f"Overdispersion falls at most {w_max}-fold, not sixteen "
                "thousand-fold")

    # Results: the four new sections, their range, and a CI that forbids any
    # claim of separation from the colon section
    # If the interval reaches below the colon value the sentence must say so;
    # if it does not, no comparison is asserted at all. Which branch is taken is
    # decided by the numbers, never by a preference for a stronger claim.
    ci_clause = ""
    if lo_ci < colon_drop:
        ci_clause = ", which overlaps the colon value"
    res_add = (
        " Four further sections fall "
        f"{n.drop(si_k):.1f}-fold (small intestine), {n.drop(ki_k):.1f}-fold "
        f"(kidney), {n.drop(pa_k):.1f}-fold (pancreas) and "
        f"{n.drop(em_k):.1f}-fold (embryo). The largest has a bootstrap 95% "
        f"confidence interval of {lo_ci:.1f} to {hi_ci:.1f}"
        f"{ci_clause}."
    )

    # Results: depth is the obvious confounder; two new sections of near-equal
    # depth settle it without a model
    depth_ratio = n.a2("lung") / n.a2(si_k)
    drop_gap = 100 * (n.drop(si_k) / n.drop("lung") - 1)
    dep_add = (
        " At the same depth, lung and small intestine differ "
        f"{depth_ratio:.2f}-fold in {ALPHA}(2 {MICRO}m) but only {drop_gap:.0f}% in "
        "how far it falls."
    )

    disc_old = "overdispersion falls five- to fourteen-fold rather than 16,384-fold"
    disc_new = (f"overdispersion falls {w_min}- to {w_max}-fold rather than "
                "16,384-fold")

    meth_old = ("We analysed four Visium HD sections and one Stereo-seq section. "
                f"Visium HD human colon cancer (Visium_HD_Human_Colon_Cancer, 6.5 mm "
                "capture area), a second section of the same tumour block "
                "(Visium_HD_Human_Colon_Cancer_P1), Visium HD human lung cancer and "
                f"Visium HD mouse brain were obtained as Space Ranger 2 {MICRO}m "
                "output and used without further alignment.")
    meth_new = (
        "We analysed eight Visium HD sections and one Stereo-seq section. Visium HD "
        "human colon cancer (Visium_HD_Human_Colon_Cancer, 6.5 mm capture area), a "
        "second section of the same tumour block (Visium_HD_Human_Colon_Cancer_P1), "
        "Visium HD human lung cancer, human pancreas "
        "(Visium_HD_Human_Pancreas), and mouse brain, small intestine, "
        "kidney and embryo (Visium_HD_Mouse_Brain, "
        "Visium_HD_Mouse_Small_Intestine, Visium_HD_Mouse_Kidney, "
        f"Visium_HD_Mouse_Embryo) were obtained as Space Ranger 2 {MICRO}m output "
        "and used without further alignment. Each dataset contributes one capture "
        "area and one section. The second colon section is a replicate; the other "
        "seven enter the scaling analysis, one per tissue. The four sections added "
        "here were fixed before their exponents were computed and were analysed by "
        "the identical pipeline."
    )

    genes_h, genes_m = n.genes("Human"), n.genes("Mouse")
    genes_old = ("retaining 18,074 genes in the human sections and 19,046 in mouse "
                 "brain.")
    genes_new = (f"retaining {genes_h:,} genes in the human sections and "
                 f"{genes_m:,} in the mouse sections.")

    data_old = " The four Visium HD sections and the Stereo-seq section analysed here"
    data_new = " The eight Visium HD sections and the Stereo-seq section analysed here"

    # ---- manuscript --------------------------------------------------------
    src = SRC / MS
    sha_before = hashlib.sha256(src.read_bytes()).hexdigest()
    doc = docx.Document(str(src))
    ps = doc.paragraphs
    before = [p.text for p in ps]
    assert len(ps) == 165, f"v5.4 paragraph count changed: {len(ps)}"

    body_before = wc(" ".join(before[15:65]))
    abs_before = wc(before[13])
    assert body_before == 3424, f"v5.4 body wordcount drifted: {body_before}"
    assert abs_before == 148, f"v5.4 abstract wordcount drifted: {abs_before}"

    def sub_once(text: str, old: str, new: str, where: str) -> str:
        c = text.count(old)
        assert c == 1, f"{where}: anchor {old[:60]!r} count={c}"
        return text.replace(old, new)

    # [0] title
    assert ps[0].text.strip() == title_v54, f"title drift: {ps[0].text[:60]!r}"
    assert len(ps[0].runs) == 1
    ps[0].runs[0].text = title_v55

    # [13] abstract: three per-section exponents -> one range.
    # runs 3/5/7 hold the superscript exponents; 2 and 4/6 are plain text.
    p = ps[13]
    assert len(p.runs) == 9, f"abstract run count {len(p.runs)}"
    r = p.runs
    r[2].text = sub_once(r[2].text, abs_old_mid, abs_new_mid, "abstract r2")
    r[3].text = f"{MINUS}{g_lo:.2f}"
    r[4].text = " to d"
    r[5].text = f"{MINUS}{g_hi:.2f}"
    r[6].text = ""
    r[7].text = ""
    r[8].text = sub_once(r[8].text, abs_old_tail, abs_new_tail, "abstract r8")
    r[8].text = sub_once(r[8].text, abs_old_drop, ".", "abstract r8 drop clause")

    # [20] introduction: keep the three exponents, add the four sections
    p = ps[20]
    assert len(p.runs) == 7, f"intro run count {len(p.runs)}"
    p.runs[0].text = sub_once(p.runs[0].text, intro_old, intro_new, "intro r0")

    # [23] Results sub-heading
    p = ps[23]
    assert len(p.runs) == 1
    p.runs[0].text = sub_once(p.text, head_old, head_new, "results heading")

    # [25] Results: the four new sections
    p = ps[25]
    assert len(p.runs) == 1
    t25 = sub_once(p.text, "The counting argument predicts a fall of 16,384.",
                   "The counting argument predicts 16,384-fold.", "results counting sentence")
    p.runs[0].text = t25 + res_add

    # [26] Results: depth rebuttal, appended where depth is already discussed
    p = ps[26]
    assert len(p.runs) == 1
    assert p.text.endswith("is not an effect of sequencing depth.")
    p.runs[0].text = p.text + dep_add

    # [55] Discussion opening
    p = ps[55]
    assert len(p.runs) == 1
    p.runs[0].text = sub_once(p.text, disc_old, disc_new, "discussion")

    # [56] Discussion: three of the seven curves are visibly flat. The paragraph
    # already argues the finite-correlation-length case with brain as the example;
    # the sections that fall least have to be named in the same breath, or the
    # figure carries outliers that the text never mentions. Names and numbers come
    # from the lock, sorted by drop, never typed.
    p = ps[56]
    assert len(p.runs) == 1
    assert p.text.startswith("Brain behaves as a field with a correlation length")
    three = sorted(n.sec.values(), key=lambda s: s["drop_fold"])[:3]
    names = [tissue_phrase(s["key"]) for s in three]
    # Three of the seven labels start with "mouse"; repeating the genus three
    # times in one sentence reads like a list of datasets rather than of tissues.
    if names[-1].startswith("mouse ") and any(x.startswith("mouse ") for x in names[:-1]):
        names[-1] = names[-1][len("mouse "):]
    add56 = (f" The {', '.join(names[:-1])} and {names[-1]} sections fall only "
             f"{three[0]['drop_fold']:.1f}- to {three[-1]['drop_fold']:.1f}-fold.")
    anchor56 = "instead of staying flat."
    p.runs[0].text = sub_once(p.text, anchor56, anchor56 + add56, "discussion [56]")

    # ---- word budget -------------------------------------------------------
    # Naming three sections and reporting four values costs more than the range
    # they replace, and the ceiling is hard. Each cut below removes a filler word
    # and no claim; none of them touches a number. The last one also fixes a
    # sentence that the seven-tissue design has outgrown: with one section per
    # tissue the design still cannot carry a between-tissue inference, but it now
    # does carry a descriptive comparison, and the wording should not deny that.
    for idx, old, new, tag in (
        (60, "cannot carry a between-tissue comparison",
             "cannot carry a between-tissue inference", "between-tissue"),
        (61, "computed at one support is read as", "computed at one support reads as",
             "read as"),
        (62, "That is of the same order as", "That is the same order as", "order"),
        (63, "And any step that changes", "Any step that changes", "levers"),
    ):
        p = ps[idx]
        assert len(p.runs) == 1, f"[{idx}] runs {len(p.runs)}"
        p.runs[0].text = sub_once(p.text, old, new, tag)

    # [67] Methods: datasets
    p = ps[67]
    assert len(p.runs) == 1
    new67 = sub_once(p.text, meth_old, meth_new, "methods datasets")
    new67 = sub_once(new67, genes_old, genes_new, "methods genes")
    p.runs[0].text = new67

    # [154] Fig. 1 caption: Fig. 1b draws three of the seven, and must say so
    p = ps[154]
    assert p.text.startswith("Fig. 1.")
    assert len(p.runs) == 4
    p.runs[1].text = sub_once(
        p.runs[1].text,
        "for colon (dark blue), lung (mid blue) and mouse brain (ochre).",
        "for colon (dark blue), lung (mid blue) and mouse brain (ochre), three of "
        "the seven sections; all seven are in Fig. S3.",
        "fig 1 caption")

    # [102] Data availability
    p = ps[102]
    assert len(p.runs) == 2
    p.runs[1].text = sub_once(p.runs[1].text, data_old, data_new, "data availability")

    changed_expected = [0, 13, 20, 23, 25, 26, 55, 56, 60, 61, 62, 63, 67, 102, 154]
    after = [p.text for p in doc.paragraphs]
    assert len(after) == len(before)
    changed = [i for i, (x, y) in enumerate(zip(before, after)) if x != y]
    assert changed == changed_expected, f"changed {changed} != {changed_expected}"

    body_after = wc(" ".join(after[15:65]))
    abs_after = wc(after[13])
    print(f"body     {body_before} -> {body_after}   limit 3500, left {3500 - body_after}")
    print(f"abstract {abs_before} -> {abs_after}   limit 150, left {150 - abs_after}")
    print(f"title    {title_v55}")
    assert body_after <= 3500, f"body over limit: {body_after}"
    assert abs_after <= 150, f"abstract over limit: {abs_after}"

    m = re.search(r"falls of ([\d.]+)- to ([\d.]+)-fold", after[13])
    assert m, "abstract range not found"
    assert abs(float(m.group(1)) - rng["min"]) < 0.06, (m.group(1), rng["min"])
    assert abs(float(m.group(2)) - rng["max"]) < 0.06, (m.group(2), rng["max"])
    assert f"falls {floor_min}- to {floor_max}-fold, not 16,384-fold" in after[0]
    assert f"{w_min}- to {w_max}-fold rather than 16,384-fold" in after[55]
    assert "four Visium HD sections" not in " ".join(after)
    for k in new_keys:
        assert n.sec[k]["label"].split()[-1].lower() in after[67].lower(), k

    dst = Path(a.dst) if a.dst else DST

    if a.dry_run:
        print("\n--- DRY RUN, nothing written ---")
        for i in changed_expected:
            print(f"\n[{i}]\n{after[i]}")
        return 0

    (dst / "02_主稿").mkdir(parents=True, exist_ok=True)
    out = dst / MS
    doc.core_properties.created = FIXED_TS
    doc.core_properties.modified = FIXED_TS
    doc.save(str(out))
    freeze_zip(out)
    sha_after = hashlib.sha256(out.read_bytes()).hexdigest()
    back = [p.text for p in docx.Document(str(out)).paragraphs]
    assert back == after, "saved file does not round-trip"
    print(f"written {out}\n  sha256 {sha_after[:16]}")

    # ---- supplementary: cover title + four table rows ----------------------
    si_src = SRC / SI
    si_sha_before = hashlib.sha256(si_src.read_bytes()).hexdigest()
    sidoc = docx.Document(str(si_src))
    sps = sidoc.paragraphs
    si_before = [p.text for p in sps]
    assert sps[1].text.strip() == title_v54, f"SI cover title drift: {sps[1].text[:60]!r}"
    old_slice = "5- to 14-fold"
    new_slice = f"{floor_min}- to {floor_max}-fold"
    for run in sps[1].runs:
        if old_slice in run.text:
            run.text = run.text.replace(old_slice, new_slice)
    assert old_slice not in sps[1].text, "SI cover title still says 14-fold"
    assert sps[1].text.strip() == title_v55, "SI cover title differs from the paper title"

    tab = sidoc.tables[0]
    assert len(tab.rows) == 6, f"Table S1 rows {len(tab.rows)}"
    anchor_tr = tab.rows[3]._tr  # brain row; new sections go in front of Stereo-seq

    def fill(row, vals: list[str]) -> None:
        for cell, val in zip(row.cells, vals):
            runs = cell.paragraphs[0].runs
            assert len(runs) >= 1
            runs[0].text = val
            for extra in runs[1:]:
                extra.text = ""

    prev = anchor_tr
    for k in new_keys:
        s = n.sec[k]
        tr = _copy.deepcopy(anchor_tr)
        prev.addnext(tr)
        prev = tr
        fill(_Row(tr, tab), [
            s["label"], "Visium HD", n.species(k), s["tissue"],
            "6.5 \u00d7 6.5 mm", f"{s['n_bins_2um']:,}",
            f"{s['umi_per_bin']:.1f}", f"{s['n_genes_noMT']:,}",
            s.get("source", "10x Genomics (public)"),
        ])
    assert len(tab.rows) == 10, f"Table S1 rows after insert {len(tab.rows)}"

    # the opening sentence still described an early draft of the file
    p2 = sps[2]
    assert len(p2.runs) == 1
    p2.runs[0].text = (
        "This file contains the dataset table (Table S1), the window-sensitivity "
        "matrix (Table S2), the permutation, dose\u2013response, curve-deviation, "
        "exponent and downstream tables (Tables S3\u2013S8), one note recording two "
        "analyses that were attempted and failed (Supplementary Note 1), and three "
        "figures (Figures S1\u2013S3). Every number is taken from the numbers lock "
        "file (NUMBERS_LOCK_v7.json) or from the derived output of the script named "
        "in the corresponding Methods subsection; nothing is transcribed by hand."
    )

    # Figure S3: the seven-section curve, so the wider range has a carrier
    fig_png = BUILD / "figs" / "FigS3.png"
    fig_pdf = BUILD / "figs" / "FigS3.pdf"
    assert fig_png.exists() and fig_pdf.exists(), f"缺 FigS3: {fig_png}"

    title_copy = _copy.deepcopy(sps[33]._p)
    cap_copy = _copy.deepcopy(sps[35]._p)
    img_p = sidoc.add_paragraph()
    img_p._p.addprevious(title_copy)
    img_p._p.addnext(cap_copy)
    img_p.add_run().add_picture(str(fig_png), width=Cm(16.4))

    tp = Paragraph(title_copy, sidoc)
    tp.runs[0].text = "Figure S3. The same decay holds on all seven sections."
    cp = Paragraph(cap_copy, sidoc)
    cp.runs[0].text = "Figure S3."
    cp.runs[1].text = (
        " \u03b1(d) normalised to its value at 2 \u00b5m, on the same axes and with "
        "the same grey guides as Fig. 1b, for all seven sections. Exact integer "
        "re-aggregation to the same eight scales was applied to each, so the curves "
        "differ only in tissue and in depth, and the four added sections were fixed "
        "before their exponents were computed. Every curve stays above the counting "
        "limit at every scale."
    )

    si_after = [p.text for p in sidoc.paragraphs]
    si_changed = [i for i, (x, y) in enumerate(zip(si_before, si_after)) if x != y]
    assert si_changed == [1, 2], f"SI changed {si_changed}"
    assert len(si_after) == len(si_before) + 3, f"SI paragraphs {len(si_after)}"
    assert "Figure S3." in si_after[-1]

    (dst / "03_补充材料").mkdir(parents=True, exist_ok=True)
    fig_dir = dst / "04_图件"
    fig_dir.mkdir(parents=True, exist_ok=True)
    for f in (fig_pdf, fig_png):
        (fig_dir / f.name).write_bytes(f.read_bytes())

    si_out = dst / SI
    sidoc.core_properties.created = FIXED_TS
    sidoc.core_properties.modified = FIXED_TS
    sidoc.save(str(si_out))
    freeze_zip(si_out)
    si_sha_after = hashlib.sha256(si_out.read_bytes()).hexdigest()
    chk = docx.Document(str(si_out))
    assert len(chk.tables[0].rows) == 10, "table rows did not survive the save"
    assert [p.text for p in chk.paragraphs] == si_after, "SI does not round-trip"
    print(f"written {si_out}\n  sha256 {si_sha_after[:16]}")

    # ---- cover letter: title and section count -----------------------------
    cl_src = SRC / CL
    cl_sha_before = hashlib.sha256(cl_src.read_bytes()).hexdigest()
    cdoc = docx.Document(str(cl_src))
    cps = cdoc.paragraphs
    cl_before = [p.text for p in cps]
    assert title_v54 in cl_before[2] and len(cps[2].runs) == 1
    cps[2].runs[0].text = cps[2].text.replace(title_v54, title_v55)

    cl_old = ("Across four Visium HD sections and one Stereo-seq section, pooled "
              "overdispersion falls 14.3-fold in colon, 13.8-fold in lung and "
              f"5.8-fold in mouse brain between 2 and 256 {MICRO}m, against the "
              "16,384-fold the counting argument predicts; the exponents are 0.554, "
              "0.546 and 0.364.")
    k_lo = rng["min_key"]
    k_hi = rng["max_key"]
    cl_new = (f"Across seven Visium HD sections, one per tissue, and one "
              f"Stereo-seq section, pooled overdispersion between 2 and 256 {MICRO}m "
              f"falls {rng['min']:.1f}-fold in {tissue_phrase(k_lo)} and "
              f"{rng['max']:.1f}-fold in {tissue_phrase(k_hi)} {EMDASH} a "
              f"{floor_min}- to {floor_max}-fold range {EMDASH} against the "
              f"16,384-fold the counting argument predicts; the exponents run "
              f"{g_lo:.2f} to {g_hi:.2f}.")
    hit = [i for i, t in enumerate(cl_before) if cl_old in t]
    assert len(hit) == 1, f"cover section sentence found in {hit}"
    ci_ = hit[0]
    assert len(cps[ci_].runs) == 1
    cps[ci_].runs[0].text = cps[ci_].text.replace(cl_old, cl_new)

    # The opening paragraph claimed "an order of magnitude"; with a section that
    # falls only 2.1-fold that is a bound, not a description, so it has to be
    # written as one.
    cl_p4_old = "the noise falls by an order of magnitude where the argument predicts four"
    cl_p4_new = ("the noise falls by well under two orders of magnitude where the "
                 "argument predicts four")
    p4 = next(i for i, t in enumerate(cl_before) if cl_p4_old in t)
    assert len(cps[p4].runs) == 1
    cps[p4].runs[0].text = cps[p4].text.replace(cl_p4_old, cl_p4_new)

    cl_after = [p.text for p in cdoc.paragraphs]
    cl_changed = [i for i, (x, y) in enumerate(zip(cl_before, cl_after)) if x != y]
    assert cl_changed == sorted({2, p4, ci_}), f"cover changed {cl_changed}"

    (dst / "01_投稿信").mkdir(parents=True, exist_ok=True)
    cl_out = dst / CL
    cdoc.core_properties.created = FIXED_TS
    cdoc.core_properties.modified = FIXED_TS
    cdoc.save(str(cl_out))
    freeze_zip(cl_out)
    cl_sha_after = hashlib.sha256(cl_out.read_bytes()).hexdigest()
    assert [p.text for p in docx.Document(str(cl_out)).paragraphs] == cl_after
    print(f"written {cl_out}\n  sha256 {cl_sha_after[:16]}")

    AUDIT.mkdir(parents=True, exist_ok=True)
    rec = dict(
        change=f"v5.4 -> v5.5: four further sections; range widens to "
               f"5- to {floor_max}-fold",
        source=str(SRC),
        numbers_source="02_OUT/extra_tissues_20261116.json (120_aggregate_extras.py, "
                       "from the raw 2 um matrices)",
        range=rng,
        manuscript=dict(sha256_before=sha_before, sha256_after=sha_after,
                        paragraphs=len(after), changed=changed,
                        body_before=body_before, body_after=body_after,
                        abstract_before=abs_before, abstract_after=abs_after),
        supplementary=dict(sha256_before=si_sha_before, sha256_after=si_sha_after,
                           changed_paragraphs=si_changed,
                           table_s1_rows=len(chk.tables[0].rows),
                           rows_added=[n.sec[k]["label"] for k in new_keys]),
        cover=dict(sha256_before=cl_sha_before, sha256_after=cl_sha_after,
                   changed=cl_changed),
        text=dict(results_sentence=res_add, depth_sentence=dep_add,
                  ci_gap_to_colon=ci_gap, gamma_range=[g_lo, g_hi]),
    )
    (AUDIT / "patch_v5_5.json").write_text(
        json.dumps(rec, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"audit  {AUDIT / 'patch_v5_5.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""Guard for v5.5. Every item is computed here, on the delivered files, from the
numbers lock. Nothing is inherited from v5.4's 39 PASS.

The point of adding four sections is that a reader can check the new range
without leaving the package: the four new drops, the confidence interval, the
depth sentence, the abstract range and the supplementary table must all agree
with `extra_tissues_20261116.json`, and the old "14-fold" must be gone from
every carrier (title, abstract, Results heading, Discussion, SI cover, cover
letter, Table S1).

Exit code 0 only if no check fails.
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

import docx
import pymupdf

PKG = Path(r"C:\Users\唐朝华\Desktop\23空间分辨率\A题_空间分辨率_投稿包_v5.5_20261116")
REPRO = Path(r"C:\Users\唐朝华\Desktop\scaling-report\reproduce")
LOCK = REPRO / "02_OUT" / "NUMBERS_LOCK_v7.json"
SUMMARY = REPRO / "02_OUT" / "extra_tissues_20261116.json"
OUT = REPRO / "07_SUBMISSION" / "_build_v55" / "audit" / "guard_v5_5.json"
SUMMARY_P = SUMMARY
if "--test" in sys.argv:
    PKG = Path(r"E://A_spatial_20261116//WORK//_test_v55")
    SUMMARY_P = Path(r"E://A_spatial_20261116//WORK//fake_summary.json")
    OUT = Path(r"E://A_spatial_20261116//WORK//_test_v55//guard_test.json")

MAN = PKG / "02_主稿" / "01_Manuscript.docx"
MAN_PDF = PKG / "02_主稿" / "01_Manuscript.pdf"
SI = PKG / "03_补充材料" / "02_Supplementary_Information.docx"
SI_PDF = PKG / "03_补充材料" / "02_Supplementary_Information.pdf"
CL = PKG / "01_投稿信" / "00_Cover_Letter.docx"
CL_PDF = PKG / "01_投稿信" / "00_Cover_Letter.pdf"

results: list[dict] = []
MINUS = "\u2212"


def ck(cond: bool, name: str, detail: str = "") -> bool:
    results.append(dict(item=name, pass_=bool(cond), detail=detail))
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"  [{detail}]" if detail else ""))
    return bool(cond)


def flat_pdf(path: Path) -> str:
    """Text with Word's line numbers stripped; they land as standalone integers."""
    d = pymupdf.open(path)
    parts = []
    for page in d:
        t = page.get_text()
        t = re.sub(r"(?m)^[ \t]*\d{1,4}[ \t]*$", " ", t)
        parts.append(t)
    return " ".join("".join(parts).split())


def main() -> int:
    summary = json.loads(Path(SUMMARY_P).read_text(encoding="utf-8"))
    rng = summary["range"]
    sec = {x["key"]: x for x in summary["sections"]}
    floor_max = rng["floor_max"]
    floor_min = rng["floor_min"]
    WORDS = {2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven",
             8: "eight", 9: "nine", 10: "ten", 11: "eleven", 12: "twelve",
             13: "thirteen", 14: "fourteen", 15: "fifteen", 16: "sixteen",
             17: "seventeen", 18: "eighteen", 19: "nineteen", 20: "twenty"}
    w_min = WORDS.get(floor_min, str(floor_min))
    w_max = WORDS.get(floor_max, str(floor_max))

    man = docx.Document(str(MAN))
    paras = [p.text for p in man.paragraphs]
    title = paras[0].strip()

    si = docx.Document(str(SI))
    si_title = si.paragraphs[1].text.strip()
    cl = docx.Document(str(CL))
    cl_text = "\n".join(p.text for p in cl.paragraphs)

    body = " ".join(paras[15:65])
    abstract = paras[13]
    man_pdf = flat_pdf(MAN_PDF)
    si_flat = flat_pdf(SI_PDF)

    print("== 1. title carried by all three deliverables, with the new bound")
    ck(title == si_title, "SI 封面标题与主稿逐字一致", f"{len(title)} chars")
    ck(title in cl_text, "投稿信引用同一标题")
    nows = lambda x: re.sub(r"\s+", "", x)
    ck(nows(title) in nows(man_pdf), "主稿 PDF 首页标题与 docx 一致", "PDF 折行后拼接比较")
    ck(f"falls {floor_min}- to {floor_max}-fold" in title,
       "标题区间为锁的 floor 值", f"{floor_min}-{floor_max}")
    ck("not 16,384-fold" in title, "标题保留 16,384 对照")

    print("== 2. length limits")
    n_body = len(body.split())
    ck(n_body <= 3500, "正文词数 ≤ 3500", f"{n_body} words, {3500 - n_body} left")
    n_abs = len(abstract.split())
    ck(n_abs <= 150, "摘要词数 ≤ 150", f"{n_abs} words, {150 - n_abs} left")
    ck(len(paras) == 165, "主稿段落数未变", f"{len(paras)}")

    print("== 3. the four new sections, value by value")
    res = next((t for t in paras if t.startswith("Pooled overdispersion in the colon")), "")
    m = re.search(
        r"Four further sections fall ([\d.]+)-fold \(small intestine\), "
        r"([\d.]+)-fold \(kidney\), ([\d.]+)-fold \(pancreas\) and "
        r"([\d.]+)-fold \(embryo\)", res)
    ck(bool(m), "Results 出现四组织句式")
    if m:
        got = [float(x) for x in m.groups()]
        want = [round(sec[k]["drop_fold"], 1) for k in
                ("mouse_small_intestine", "mouse_kidney", "human_pancreas", "mouse_embryo")]
        ck(got == want, "Results 四个倍数与锁一致", f"{got} vs {want}")

    ci_m = re.search(r"confidence interval of (\d+\.\d) to (\d+\.\d)", res)
    longest = max(("mouse_small_intestine", "mouse_kidney", "human_pancreas",
                   "mouse_embryo"), key=lambda k: sec[k]["drop_fold"])
    want_ci = sec[longest]["boot"]["ci95"]
    ck(bool(ci_m), "Results 给出 bootstrap 区间")
    if ci_m:
        got_ci = [float(ci_m.group(1)), float(ci_m.group(2))]
        ck(abs(got_ci[0] - want_ci[0]) < 0.06 and abs(got_ci[1] - want_ci[1]) < 0.06,
           "区间两界与锁一致", f"{got_ci} vs {[round(v, 2) for v in want_ci]}")
    ck("overlaps the colon value" in res or want_ci[0] > sec["colon"]["drop_fold"],
       "区间与结肠的关系如实（无 '显著超出'）")
    ck("significantly" not in man_pdf.lower(), "全文无 'significantly' 类夸大")

    dep = next((t for t in paras if "the same depth, lung and small intestine" in t), "")
    dm = re.search(r"differ ([\d.]+)-fold in .\(2 .m\) but only (\d+)% in how far", dep)
    ck(bool(dm), "深度对照句存在")
    if dm:
        want_ratio = round(sec["lung"]["alpha_2um"] / sec["mouse_small_intestine"]["alpha_2um"], 2)
        want_gap = round(100 * (sec["mouse_small_intestine"]["drop_fold"] /
                                sec["lung"]["drop_fold"] - 1))
        ck(abs(float(dm.group(1)) - want_ratio) < 0.006 and int(dm.group(2)) == want_gap,
           "深度句两个数与锁一致", f"{dm.group(1)}/{dm.group(2)} vs {want_ratio}/{want_gap}")

    am = re.search(r"falls of ([\d.]+)- to ([\d.]+)-fold", abstract)
    ck(bool(am), "摘要给出区间")
    if am:
        ck(abs(float(am.group(1)) - rng["min"]) < 0.06 and
           abs(float(am.group(2)) - rng["max"]) < 0.06,
           "摘要区间与锁一致", f"{am.group(1)}-{am.group(2)} vs "
           f"{rng['min']:.2f}-{rng['max']:.2f}")

    # Three of the seven curves are visibly flat; the Discussion names them where
    # the finite-correlation-length argument is made. Names and numbers both come
    # from the lock, so a re-count cannot leave a stale label.
    print("== 3b. the flattest curves are named where the argument is made")
    flat_s = next((t for t in paras if "sections fall only" in t), "")
    three = sorted(sec.values(), key=lambda s: s["drop_fold"])[:3]
    tabs = [s["label"].replace("Visium HD ", "").lower() for s in three]
    if tabs[-1].startswith("mouse ") and any(x.startswith("mouse ") for x in tabs[:-1]):
        tabs[-1] = tabs[-1][len("mouse "):]
    want_s = (f"The {', '.join(tabs[:-1])} and {tabs[-1]} sections fall only "
              f"{three[0]['drop_fold']:.1f}- to {three[-1]['drop_fold']:.1f}-fold.")
    ck(bool(flat_s), "Discussion 点名三片下落最小的组织")
    ck(want_s in flat_s, "三片名称与区间均与锁一致",
       f"{three[0]['drop_fold']:.2f} / {three[-1]['drop_fold']:.2f}")

    print("== 4. the old bound is gone from every carrier")
    ck(f"Overdispersion falls at most {w_max}-fold" in " ".join(paras),
       "Results 小标题用新上界")
    ck(f"{w_min}- to {w_max}-fold rather than 16,384-fold" in man_pdf,
       "Discussion 用新上界")
    ck("14-fold" not in man_pdf, "主稿 PDF 不再出现 14-fold")
    ck("fourteen-fold" not in man_pdf, "主稿 PDF 不再出现 fourteen-fold")
    ck("5- to 14-fold" not in si_flat, "补充材料封面不再出现 5- to 14-fold")
    ck("5- to 14-fold" not in cl_text, "投稿信不再出现 5- to 14-fold")
    ck("four Visium HD sections" not in man_pdf, "不再写 four Visium HD sections")

    print("== 5. section count consistent across carriers")
    ck("eight Visium HD sections" in paras[67], "Methods 写 eight Visium HD sections")
    ck("eight Visium HD sections" in paras[102], "Data availability 写 eight")
    ck("seven Visium HD sections" in abstract, "摘要写 seven Visium HD sections")
    ck("four further sections" in paras[20], "引言承认新增四片")
    for k in ("mouse_small_intestine", "mouse_kidney", "human_pancreas", "mouse_embryo"):
        ck(sec[k]["dataset"] in man_pdf, f"Methods 点名 {sec[k]['dataset']}")

    print("== 6. supplementary table S1")
    tab = si.tables[0]
    ck(len(tab.rows) == 10, "Table S1 共 10 行", f"{len(tab.rows)}")
    rows = {r.cells[0].text.strip(): [c.text.strip() for c in r.cells] for r in tab.rows}
    for k in ("mouse_small_intestine", "mouse_kidney", "human_pancreas", "mouse_embryo"):
        lab = sec[k]["label"]
        row = rows.get(lab)
        ck(row is not None, f"Table S1 收录 {lab}")
        if row:
            ck(row[5] == f"{sec[k]['n_bins_2um']:,}", f"{lab} 2 µm 格数与锁一致", row[5])
            ck(abs(float(row[6]) - sec[k]["umi_per_bin"]) < 0.06,
               f"{lab} 每格转录本与锁一致", row[6])
            ck(row[7] == f"{sec[k]['n_genes_noMT']:,}", f"{lab} 保留基因数与锁一致", row[7])
            ck(row[2] == sec[k]["species"], f"{lab} 物种正确", row[2])

    print("== 7. figures and the fit window")
    cap_start = next(i for i, t in enumerate(paras) if t.strip().startswith("Fig. 1."))
    body_only = " ".join(paras[:cap_start])
    for n in range(1, 7):
        hits = len(re.findall(r"Fig\. ?%d(?![0-9])" % n, body_only))
        ck(hits >= 1, f"图 {n} 在正文被引用", f"{hits} 处")
    ck("β 2–24 grid units" in man_pdf, "图 4c 图注给出 β 2–24")
    ck("2–24" in next(p.text for p in man.paragraphs if "log ρ = c − β log h" in p.text),
       "Methods 的 β 窗口为 2–24")
    ck("β over 2–24 grid units" in si_flat, "Table S2 注的 β 窗口为 2–24")

    print("== 8. nothing left as a placeholder")
    for bad in ["TO BE INSERTED", "XXX", "TBD", "[URL", "PLACEHOLDER"]:
        ck(bad not in man_pdf, f"主稿无占位符 {bad!r}")

    print("== 9. the package is visible, with its repository")
    ck("https://github.com/tang-sir38/scaling-report" in nows(man_pdf),
       "Code availability 有仓库地址")
    ck("scaling_report" in man_pdf, "正文点名可安装包 scaling_report")
    ck("installs with pip" in man_pdf, "给出安装方式")
    ck("https://github.com/tang-sir38/scaling-report" in cl_text, "投稿信给出仓库地址")
    ck("will be deposited" not in cl_text, "投稿信不再写'将来存入'")
    ck("installable reporting package" in cl_text, "投稿信提到可安装包")

    print("== 10. the per-gene numbers still match the lock")
    lock = json.loads(LOCK.read_text(encoding="utf-8"))["spatial_proxy"]["by_tissue"]
    tis = ("colon", "lung", "brain")
    want_moran = [round(lock[t]["rho_gamma_moran"], 3) for t in tis]
    want_band = [round([b for b in lock[t]["band_rho"] if b["k"] == 10][0]["median"], 2)
                 for t in tis]
    want_semi = [round(lock[t]["semiparam"]["pearson"], 2) for t in tis]
    disc = next(p.text for p in man.paragraphs if "Moran's I at" in p.text)
    got_moran = [float(x.replace(MINUS, "-")) for x in re.findall(MINUS + r"\d\.\d{3}", disc)]
    ck(got_moran == want_moran, "Moran 三个值与锁一致", f"{got_moran} vs {want_moran}")
    got_band = [float(x.replace(MINUS, "-"))
                for x in re.findall(MINUS + r"\d\.\d{2}(?!\d)", disc)][:3]
    ck(got_band == want_band, "丰度分层三个值与锁一致", f"{got_band} vs {want_band}")
    meth = next(p.text for p in man.paragraphs
                if "semiparametric residual control gives" in p.text)
    got_semi = [float(x.replace(MINUS, "-"))
                for x in re.findall(MINUS + r"\d\.\d{2}(?!\d)", meth)][:3]
    ck(got_semi == want_semi, "半参数控制三个值与锁一致", f"{got_semi} vs {want_semi}")
    ck("0.194" not in man_pdf, "未引用被判定不适用的线性偏相关")

    print("== 11. references")
    ref_start = next(i for i, t in enumerate(paras) if t.strip() == "References")
    refs = [t for t in paras[ref_start:] if re.match(r"^\d+\.\s", t.strip())]
    ck(len(refs) <= 50, "参考文献 ≤ 50", f"{len(refs)} 条")
    cited = set()
    for mm in re.finditer(r"\[([\d,\s–-]+)\]", " ".join(paras[:ref_start])):
        for part in mm.group(1).split(","):
            part = part.strip()
            if re.fullmatch(r"\d+", part):
                cited.add(int(part))
            elif re.fullmatch(r"\d+\s*[–-]\s*\d+", part):
                a, b = re.split(r"[–-]", part)
                cited.update(range(int(a), int(b) + 1))
    missing = [i for i in range(1, len(refs) + 1) if i not in cited]
    ck(not missing, "参考文献全部被正文引用", f"缺 {missing}" if missing else "0 条未引")

    print("== 12. supplementary material intact")
    si_txt = "\n".join(p.text for p in si.paragraphs)
    found_tabs = sorted(re.findall(r"^Table (S\d+)\.", si_txt, re.M))
    ck(found_tabs == [f"S{i}" for i in range(1, 9)], "补充表 S1–S8 齐", str(found_tabs))
    ck(len(si.tables) == 8, "补充表对象 8 个", f"{len(si.tables)}")
    ck(len(si.inline_shapes) == 3, "补充图 3 张内嵌", f"{len(si.inline_shapes)}")
    ck("Figure S3." in si_txt, "Figure S3 图注存在")
    ck((PKG / "04_图件" / "FigS3.pdf").exists() and (PKG / "04_图件" / "FigS3.png").exists(),
       "FigS3 矢量与位图都在 04_图件")
    ck("Fig. S3" in man_pdf, "主稿图 1 图注指向 Fig. S3", "七组织曲线有载体")
    ck(len(si.paragraphs) == 39, "补充材料段落数 = 39", f"{len(si.paragraphs)}")

    print("== 13. page counts and no text loss in the PDF")
    ck(len(pymupdf.open(MAN_PDF)) >= 27, "主稿 PDF 页数 ≥ 27",
       f"{len(pymupdf.open(MAN_PDF))} 页")
    ck(len(pymupdf.open(SI_PDF)) >= 6, "补充材料 PDF 页数 ≥ 6",
       f"{len(pymupdf.open(SI_PDF))} 页")
    ck(len(pymupdf.open(CL_PDF)) == 1, "投稿信仍是 1 页",
       f"{len(pymupdf.open(CL_PDF))} 页")
    ck("suggest but not test" not in man_pdf, "旧的 'suggest but not test' 未回流")

    print("== 14. manifest covers every file that will be uploaded")
    man_path = PKG / "MANIFEST.json"
    if man_path.exists():
        listed = {e["file"] for e in json.loads(man_path.read_text(encoding="utf-8"))["files"]}
        present = {p.relative_to(PKG).as_posix() for p in PKG.rglob("*") if p.is_file()} - {"MANIFEST.json"}
        ck(listed == present, "MANIFEST 覆盖包内全部文件",
           f"{len(listed)} 条 / {len(present)} 文件")
    else:
        ck(False, "MANIFEST 存在")

    fails = [r for r in results if not r["pass_"]]
    status = "PASS" if not fails else "BLOCKED"
    rec = dict(
        status=status, pass_=len(results) - len(fails), fail=len(fails),
        total=len(results),
        fails=[f"{r['item']} :: {r['detail']}" for r in fails],
        checks=results,
        method="每项均在本次运行的 v5.5 文件上实算；数字一侧对 extra_tissues_20261116.json，"
               "不沿用 v5.4 的任何 carried-over 结论",
        range=rng,
        files=dict(
            manuscript=dict(sha256=hashlib.sha256(MAN.read_bytes()).hexdigest(),
                            pdf_sha256=hashlib.sha256(MAN_PDF.read_bytes()).hexdigest(),
                            pages=len(pymupdf.open(MAN_PDF)), body_words=n_body),
            supplementary=dict(sha256=hashlib.sha256(SI.read_bytes()).hexdigest(),
                               pages=len(pymupdf.open(SI_PDF)),
                               table_s1_rows=len(tab.rows)),
            cover=dict(sha256=hashlib.sha256(CL.read_bytes()).hexdigest(),
                       pages=len(pymupdf.open(CL_PDF))),
        ),
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(rec, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"\n{status}: {rec['pass_']}/{rec['total']} pass, {rec['fail']} fail")
    print(f"written {OUT}")
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())

"""Guard for v5.4. Runs against the files in the package directory, not against
a previous run's conclusions.

v5.3 carried 19 checks over from v5.2 because 103_guard_v5.py is not in the
tree. v5.4 changes the manuscript, so every check has to be re-established on
this byte image. This script does that: each item is computed here from the
docx, the PDF and the numbers lock, and nothing is inherited.

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

PKG = Path(r"C:\Users\唐朝华\Desktop\23空间分辨率\A题_空间分辨率_投稿包_v5.4_20261116")
LOCK = Path(r"C:\Users\唐朝华\Desktop\scaling-report\reproduce\02_OUT\NUMBERS_LOCK_v7.json")
OUT = Path(r"C://Users//唐朝华//Desktop//scaling-report//reproduce//07_SUBMISSION//_build_v54//audit//guard_v5_4.json")

MAN = PKG / "02_主稿" / "01_Manuscript.docx"
MAN_PDF = PKG / "02_主稿" / "01_Manuscript.pdf"
SI = PKG / "03_补充材料" / "02_Supplementary_Information.docx"
SI_PDF = PKG / "03_补充材料" / "02_Supplementary_Information.pdf"
CL = PKG / "01_投稿信" / "00_Cover_Letter.docx"
CL_PDF = PKG / "01_投稿信" / "00_Cover_Letter.pdf"

results: list[dict] = []


def ck(cond: bool, name: str, detail: str = "") -> bool:
    results.append(dict(item=name, pass_=bool(cond), detail=detail))
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"  [{detail}]" if detail else ""))
    return bool(cond)


def flat_pdf(path: Path) -> str:
    """Text with Word's line numbers removed.

    The manuscript carries line numbers (the journal asks for them), and they
    land in the text layer as standalone integers: "..., not 1 16,384-fold".
    Leaving them in breaks every substring test that spans a line, so they are
    dropped before the whitespace is collapsed.
    """
    d = pymupdf.open(path)
    parts = []
    for page in d:
        t = page.get_text()
        t = re.sub(r"(?m)^[ \t]*\d{1,4}[ \t]*$", " ", t)
        parts.append(t)
    return " ".join("".join(parts).split())


def main() -> int:
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

    # ---- 1. title consistency across the three deliverables
    ck(title == si_title, "SI 封面标题与主稿逐字一致", f"{len(title)} chars")
    ck(title in cl_text, "投稿信引用同一标题")
    nows = lambda x: re.sub(r"\s+", "", x)
    ck(nows(title) in nows(man_pdf), "主稿 PDF 首页标题与 docx 一致",
       "PDF 折行后拼接比较")

    # ---- 2. length limits
    n_body = len(body.split())
    ck(n_body <= 3500, "正文词数 ≤ 3500", f"{n_body} words, {3500 - n_body} left")
    n_abs = len(abstract.split())
    ck(n_abs <= 150, "摘要词数 ≤ 150", f"{n_abs} words, {150 - n_abs} left")
    ck(len([p for p in paras if p.strip()]) > 100, "主稿未截断", f"{len(paras)} paragraphs")

    # ---- 3. every main figure cited in the body (excluding captions)
    cap_start = next(i for i, t in enumerate(paras) if t.strip().startswith("Fig. 1."))
    body_only = " ".join(paras[:cap_start])
    for n in range(1, 7):
        hits = len(re.findall(r"Fig\. ?%d(?![0-9])" % n, body_only))
        ck(hits >= 1, f"图 {n} 在正文被引用", f"{hits} 处")

    # ---- 4. the fit window agrees in all three places
    ck("β 2–24 grid units" in man_pdf, "图 4c 图注给出 β 2–24")
    ck("2–24" in next(p.text for p in man.paragraphs if "log ρ = c − β log h" in p.text),
       "Methods 的 β 窗口为 2–24")
    ck("β over 2–24 grid units" in si_flat, "Table S2 注的 β 窗口为 2–24（已修）")

    # ---- 5. nothing left as a placeholder
    for bad in ["TO BE INSERTED", "XXX", "TBD", "[URL", "PLACEHOLDER"]:
        ck(bad not in man_pdf, f"主稿无占位符 {bad!r}")

    # ---- 6. the package is visible, with its repository
    ck("https://github.com/tang-sir38/scaling-report" in nows(man_pdf),
       "Code availability 有仓库地址", "PDF 在连字符处折行，去空白比较")
    ck("scaling_report" in man_pdf, "正文点名可安装包 scaling_report")
    ck("installs with pip" in man_pdf, "给出安装方式")
    ck("https://github.com/tang-sir38/scaling-report" in cl_text,
       "投稿信给出仓库地址")
    ck("will be deposited" not in cl_text,
       "投稿信不再写'将来存入'（时态与主稿一致）")
    ck("installable reporting package" in cl_text, "投稿信提到可安装包")

    # ---- 7. the new spatial numbers match the lock, value by value
    lock = json.loads(LOCK.read_text(encoding="utf-8"))["spatial_proxy"]["by_tissue"]
    tis = ("colon", "lung", "brain")
    want_moran = [round(lock[t]["rho_gamma_moran"], 3) for t in tis]
    want_band = [round([b for b in lock[t]["band_rho"] if b["k"] == 10][0]["median"], 2) for t in tis]
    want_semi = [round(lock[t]["semiparam"]["pearson"], 2) for t in tis]

    disc = next(p.text for p in man.paragraphs if "Moran's I at" in p.text)
    got_moran = [float(x.replace("\u2212", "-")) for x in re.findall("\u2212\\d\\.\\d{3}", disc)]
    ck(got_moran == want_moran, "Moran 三个值与锁一致", f"{got_moran} vs {want_moran}")

    def grab(block: str, k: int) -> list[float]:
        v = [float(x.replace("\u2212", "-")) for x in re.findall("\u2212\\d\\.\\d{2}(?!\\d)", block)]
        return v[:k]

    meth = next(p.text for p in man.paragraphs if "semiparametric residual control gives" in p.text)
    ck(grab(disc, 3) == want_band, "丰度分层三个值与锁一致", f"{grab(disc, 3)} vs {want_band}")
    semi_got = [float(x.replace("\u2212", "-")) for x in re.findall("\u2212\\d\\.\\d{2}(?!\\d)", meth)]
    ck(semi_got[:3] == want_semi, "半参数控制三个值与锁一致", f"{semi_got[:3]} vs {want_semi}")

    # the linear partial correlation must not appear anywhere as a claim
    ck("0.194" not in man_pdf and "+0.194" not in man_pdf,
       "未引用被判定不适用的线性偏相关")

    # ---- 8. references
    ref_start = next(i for i, t in enumerate(paras) if t.strip() == "References")
    refs = [t for t in paras[ref_start:] if re.match(r"^\d+\.\s", t.strip())]
    ck(len(refs) <= 50, "参考文献 ≤ 50", f"{len(refs)} 条")
    cited = set()
    for m in re.finditer(r"\[([\d,\s–-]+)\]", " ".join(paras[:ref_start])):
        for part in m.group(1).split(","):
            part = part.strip()
            if re.fullmatch(r"\d+", part):
                cited.add(int(part))
            elif re.fullmatch(r"\d+\s*[–-]\s*\d+", part):
                a, b = re.split(r"[–-]", part)
                cited.update(range(int(a), int(b) + 1))
    missing = [i for i in range(1, len(refs) + 1) if i not in cited]
    ck(not missing, "参考文献全部被正文引用", f"缺 {missing}" if missing else "0 条未引")

    # ---- 9. supplementary material intact
    si_txt = "\n".join(p.text for p in si.paragraphs)
    found_tabs = sorted(re.findall(r"^Table (S\d+)\.", si_txt, re.M))
    ck(found_tabs == [f"S{i}" for i in range(1, 9)], "补充表 S1–S8 齐", str(found_tabs))
    ck(len(si.tables) == 8, "补充表对象 8 个", f"{len(si.tables)}")
    ck(len(si.inline_shapes) == 2, "补充图 2 张内嵌", f"{len(si.inline_shapes)}")

    # ---- 10. page counts and no text loss in the PDF
    ck(len(pymupdf.open(MAN_PDF)) >= 27, "主稿 PDF 页数 ≥ 27", f"{len(pymupdf.open(MAN_PDF))} 页")
    ck(len(pymupdf.open(SI_PDF)) >= 6, "补充材料 PDF 页数 ≥ 6", f"{len(pymupdf.open(SI_PDF))} 页")
    ck("suggest but not test" not in man_pdf, "旧的 'suggest but not test' 已替换")

    # ---- 11. manifest covers every file that will be uploaded
    man_path = PKG / "MANIFEST.json"
    if man_path.exists():
        listed = {e["file"] for e in json.loads(man_path.read_text(encoding="utf-8"))["files"]}
        # the manifest does not list itself; it does list this guard's output,
        # which is deterministic (no timestamps), so two runs give the same bytes
        present = {p.relative_to(PKG).as_posix() for p in PKG.rglob("*") if p.is_file()} - {"MANIFEST.json"}
        ck(listed == present, "MANIFEST 覆盖包内全部文件", f"{len(listed)} 条 / {len(present)} 文件")

    fails = [r for r in results if not r["pass_"]]
    status = "PASS" if not fails else "BLOCKED"
    rec = dict(
        status=status,
        pass_=len(results) - len(fails),
        fail=len(fails),
        total=len(results),
        fails=[f"{r['item']} :: {r['detail']}" for r in fails],
        checks=results,
        method="每项均在本次运行的 v5.4 文件上实算，不沿用 v5.3 的 carried-over 结论",
        files=dict(
            manuscript=dict(sha256=hashlib.sha256(MAN.read_bytes()).hexdigest(),
                            pdf_sha256=hashlib.sha256(MAN_PDF.read_bytes()).hexdigest(),
                            pages=len(pymupdf.open(MAN_PDF)), body_words=n_body),
            supplementary=dict(sha256=hashlib.sha256(SI.read_bytes()).hexdigest(),
                               pages=len(pymupdf.open(SI_PDF))),
        ),
    )
    OUT.write_text(json.dumps(rec, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"\n{status}: {rec['pass_']}/{rec['total']} pass, {rec['fail']} fail")
    print(f"written {OUT}")
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())

# -*- coding: utf-8 -*-
"""111_guard_v56.py — v5.6 守卫。只对本次构建出的这两份文件有效。

不引用 110 的任何中间值：每个数字都从 E:/xen_p2/out 的 JSON 现场重算，
再与 docx 的 OOXML 文本和 PDF 的渲染文本比对。

用法：python 111_guard_v56.py
退出码 0 = 全过；1 = 有失败项。
"""
from __future__ import annotations

import hashlib
import json
import math
import re
import sys
from pathlib import Path

import docx
import pymupdf

OLD = Path(r"C:\Users\唐朝华\Desktop\23空间分辨率\A题_空间分辨率_投稿包_v5.5.1_20261117")
PKG = Path(r"C:\Users\唐朝华\Desktop\23空间分辨率\A题_空间分辨率_投稿包_v5.6_20261117")
XEN = Path(r"E:\xen_p2\out")
MAN = PKG / "02_主稿" / "01_Manuscript.docx"
SI = PKG / "03_补充材料" / "02_Supplementary_Information.docx"
MAN_PDF = PKG / "02_主稿" / "01_Manuscript.pdf"
SI_PDF = PKG / "03_补充材料" / "02_Supplementary_Information.pdf"

CHECKS: list[tuple[str, bool, str]] = []


def ck(name: str, ok: bool, detail: str = "") -> None:
    CHECKS.append((name, bool(ok), detail))


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def docx_text(p: Path) -> str:
    d = docx.Document(str(p))
    parts = [x.text for x in d.paragraphs]
    for t in d.tables:
        for r in t.rows:
            parts.extend(c.text for c in r.cells)
    return "\n".join(parts)


def pdf_text(p: Path) -> str:
    with pymupdf.open(str(p)) as d:
        return "".join(pg.get_text() for pg in d)


def main() -> int:
    # ---------------------------------------------------------- 0 原包未动
    old_man = json.loads((OLD / "MANIFEST.json").read_text(encoding="utf-8"))
    bad = []
    for e in old_man["files"]:
        f = OLD / e["file"]
        if not f.exists() or f.stat().st_size != e["bytes"] or sha(f) != e["sha256"]:
            bad.append(e["file"])
    ck("v5.5.1 原包 36 文件逐字节未动", not bad, f"changed={bad[:4]}")

    # ---------------------------------------------------------- 1 主稿结构
    dm = docx.Document(str(MAN))
    paras = [p.text for p in dm.paragraphs]
    body = len(" ".join(paras[15:65]).split())
    abstract = len(paras[13].split())
    ck("主稿段落数 165", len(paras) == 165, f"{len(paras)}")
    ck("正文词数 ≤3500", body <= 3500, f"{body}")
    ck("摘要词数 ≤150", abstract <= 150, f"{abstract}")

    mt = "\n".join(paras) + "\n" + "\n".join(
        c.text for t in dm.tables for r in t.rows for c in r.cells)
    for s, n in [("partly for that reason", 1),
                 ("Table S5, S6, S9", 1),
                 ("and one such section analysed here reproduces the failure (Note 1)", 1),
                 ("fail for that reason", 0),
                 ("in a 2 µm square", 0)]:
        ck(f"主稿 '{s[:44]}' 出现 {n} 次", mt.count(s) == n, f"{mt.count(s)}")

    # ---------------------------------------------------------- 2 SI 结构
    ds = docx.Document(str(SI))
    sp = [p.text for p in ds.paragraphs]
    st = docx_text(SI)
    ck("SI 段落 48", len(sp) == 48, f"{len(sp)}")
    ck("SI 表格 9", len(ds.tables) == 9, f"{len(ds.tables)}")
    ck("SI 新节标题存在", st.count("A third platform: single-molecule coordinates") == 1,
       f"{st.count('A third platform: single-molecule coordinates')}")
    ck("SI 旧串 'a hard limit' 已删", "a hard limit of this data set" not in st,
       f"{st.count('a hard limit of this data set')}")
    ck("SI 新串 'one limit ... not the only one' 存在",
       st.count("but one limit of this data set, and not the only one (below)") == 1, "")
    ck("SI 第 21 段新句存在",
       st.count("the single-molecule check below shows that the limit does not "
                "disappear when the grid does") == 1, "")
    ck("SI 无 [SCOTE] 占位符", "[SCOTE]" not in st, "")
    ck("SCOTE 引用写全", "doi:10.64898/2026.10.09.757915" in st, "")

    # ---------------------------------------------------------- 3 数字复算
    pred = json.loads((PRED := XEN / "xen_predict_closing10.json").read_text("utf-8"))
    rows = pred["rows"]
    ds_v = [r["d"] for r in rows]
    pv = [r["alpha_pred"] for r in rows]
    ov = [r["alpha_true"] for r in rows]
    # 独立重算斜率（log-log，1–256）
    x = [math.log(v) for v in ds_v]
    xm = sum(x) / len(x)
    def slope(yv):
        ym = sum(yv) / len(yv)
        num = sum((a - xm) * (b - ym) for a, b in zip(x, yv))
        den = sum((a - xm) ** 2 for a in x)
        return num / den
    s_pred = slope([math.log(v) for v in pv])
    s_obs = slope([math.log(v) for v in ov])
    ck("斜率 predicted −0.343", abs(s_pred - (-0.34346)) < 5e-4, f"{s_pred:.5f}")
    ck("斜率 observed −0.024", abs(s_obs - (-0.02447)) < 5e-4, f"{s_obs:.5f}")
    ck("斜率差 0.319", abs(abs(s_pred - s_obs) - 0.31899) < 5e-4,
       f"{abs(s_pred - s_obs):.5f}")
    # 2 µm 归一后的最大偏差
    i2 = ds_v.index(2.0)
    sc = ov[i2] / pv[i2]
    devs = [0.0 if d == 2.0 else math.log10(abs(p * sc / o))
            for d, p, o in zip(ds_v, pv, ov)]
    mx = max(abs(v) for v in devs)
    ck("最大 log10 偏差 0.831", abs(mx - 0.83122) < 5e-4, f"{mx:.5f}")
    ck("6.8 倍 = 10^0.831", abs(10 ** mx - 6.78) < 0.05, f"{10**mx:.2f}")

    # Table S9 逐格
    t9 = ds.tables[6]
    ck("Table S9 为 10 行 4 列", len(t9.rows) == 10 and len(t9.columns) == 4,
       f"{len(t9.rows)}x{len(t9.columns)}")
    head = [c.text.strip() for c in t9.rows[0].cells]
    ck("Table S9 表头", head == ["d (µm)", "α predicted", "α observed", "log10 deviation"],
       str(head))
    cell_bad = []
    for i, r in enumerate(rows, start=1):
        got = [c.text.strip() for c in t9.rows[i].cells]
        want_d = f"{r['d']:g}"
        want_p = f"{r['alpha_pred']:.3f}"
        want_o = f"{r['alpha_true']:.3f}"
        dv = devs[i - 1]
        want_v = "0.000" if abs(dv) < 5e-4 else f"{dv:+.3f}".replace("-", "\u2212")
        if [want_d, want_p, want_o, want_v] != got:
            cell_bad.append((want_d, got, [want_d, want_p, want_o, want_v]))
    ck("Table S9 九行与 JSON 复算一致", not cell_bad, str(cell_bad[:2]))

    a = json.loads((XEN / "xen_alpha_closing10.json").read_text("utf-8"))
    dsv = a["ds"]
    ap = a["alpha_pool"]
    ck("α(4 µm)=3.11", abs(ap[dsv.index(4.0)] - 3.1101) < 5e-3,
       f"{ap[dsv.index(4.0)]:.4f}")
    ck("α(0.5 µm)=−0.45", abs(ap[dsv.index(0.5)] - (-0.45353)) < 5e-3,
       f"{ap[dsv.index(0.5)]:.4f}")

    def seg(x0, x1):
        i, j = dsv.index(x0), dsv.index(x1)
        return (math.log(ap[j]) - math.log(ap[i])) / (math.log(x1) - math.log(x0))
    k18, k864 = seg(1.0, 8.0), seg(8.0, 64.0)
    ck("k(1–8)=+0.56", abs(k18 - 0.56434) < 5e-3, f"{k18:+.4f}")
    ck("k(8–64)=−0.19", abs(k864 - (-0.19042)) < 5e-3, f"{k864:+.4f}")

    g = json.loads((XEN / "gamma_compare.json").read_text("utf-8"))
    ck("per-gene ρ=0.871", abs(g["spearman_paired"] - 0.87136) < 5e-4,
       f"{g['spearman_paired']:.4f}")
    ck("n paired=403", g["n_paired"] == 403, str(g["n_paired"]))
    ck("median γ 0.54 / 0.41",
       abs(abs(g["visium_median_gamma"]) - 0.541) < 5e-3
       and abs(abs(g["xenium_median_gamma"]) - 0.412) < 5e-3,
       f"{abs(g['visium_median_gamma']):.3f}/{abs(g['xenium_median_gamma']):.3f}")
    # pooled 两个数（同 408 基因，各自权重）
    vv = abs(g["V_weight_on_visium_alpha"])
    xx = abs(g["X_weight_on_xenium_alpha"])
    ck("pooled 0.42 vs 0.10", abs(vv - 0.42) < 5e-3 and abs(xx - 0.10) < 5e-3,
       f"{vv:.4f}/{xx:.4f}")
    ck("gap 0.32", abs(abs(vv - xx) - 0.32) < 5e-3, f"{abs(vv-xx):.4f}")
    vv2 = abs(g["V_weight_on_xenium_alpha"])
    ck("gap 0.34（both Visium 权重）", abs(abs(vv - vv2) - 0.34) < 5e-3,
       f"{abs(vv-vv2):.4f}")

    t6 = json.loads((XEN / "xen_thin06b.json").read_text("utf-8"))
    ck("thinning 0.001", abs(t6["dgamma_S2"] - 0.0014) < 5e-4, f"{t6['dgamma_S2']:.5f}")
    ck("panel 0.089", abs(t6["dgamma_S4"] - 0.0893) < 5e-4, f"{t6['dgamma_S4']:.5f}")

    prep = json.loads((XEN / "xen_prep.json").read_text("utf-8"))
    ck("transcripts 66,248,117", prep["n_rows"] == 66248117, f"{prep['n_rows']:,}")
    ck("panel 422", prep["n_targets"] == 422, str(prep["n_targets"]))
    ck("pass 50,594,830", prep["n_pass_qv20_gene_inpanel"] == 50594830,
       f"{prep['n_pass_qv20_gene_inpanel']:,}")

    # ---------------------------------------------------------- 4 PDF 载体
    ck("主稿 PDF 页数 28", pymupdf.open(str(MAN_PDF)).page_count == 28,
       str(pymupdf.open(str(MAN_PDF)).page_count))
    ck("SI PDF 页数 8", pymupdf.open(str(SI_PDF)).page_count == 8,
       str(pymupdf.open(str(SI_PDF)).page_count))
    mpt = re.sub(r"\s+", " ", pdf_text(MAN_PDF))
    for s in ["partly for that reason", "Table S5, S6, S9",
              "one such section analysed here reproduces the failure"]:
        ck(f"主稿 PDF 含 '{s[:36]}'", mpt.count(s) >= 1, f"{mpt.count(s)}")
    spt = re.sub(r"\s+", " ", pdf_text(SI_PDF))
    for s in ["A third platform: single-molecule coordinates",
              "Table S9. The integral prediction repeated",
              "but one limit of this data set, and not the only one",
              "does not disappear when the grid does",
              "0.871", "66,248,117", "50,594,830"]:
        ck(f"SI PDF 含 '{s[:36]}'", spt.count(s) >= 1, f"{spt.count(s)}")
    # 七个关键数字在 SI PDF 里都在
    for num in ["0.42", "0.10", "0.32", "0.34", "0.343", "0.024", "0.83"]:
        ck(f"SI PDF 含数字 {num}", num in spt, "")

    # ---------------------------------------------------------- 5 输出
    fails = [c for c in CHECKS if not c[1]]
    w = max(len(c[0]) for c in CHECKS)
    for name, ok, detail in CHECKS:
        print(f"  {'PASS' if ok else 'FAIL'}  {name:<{w}}  {detail}")
    print(f"\n{len(CHECKS) - len(fails)}/{len(CHECKS)} PASS")

    rec = dict(generator="111_guard_v56.py", n_checks=len(CHECKS),
               n_pass=len(CHECKS) - len(fails),
               pass_all=not fails,
               files={str(p.relative_to(PKG)): sha(p) for p in
                      (MAN, SI, MAN_PDF, SI_PDF)},
               checks=[dict(name=n, ok=o, detail=d) for n, o, d in CHECKS])
    out = PKG / "05_Audit" / "guard_v5_6.json"
    out.write_text(json.dumps(rec, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"wrote {out}")
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())

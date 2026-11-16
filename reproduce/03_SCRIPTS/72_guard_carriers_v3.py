#!/usr/bin/env python
"""Guard 3 — carriers (v3.0).

Every Fig./Table/Supplementary identifier cited in the text must have a real
carrier: a figure file, a legend entry, or a table object. And every carrier
must be cited at least once (no orphan panels).
"""
import re
import sys
from pathlib import Path

from docx import Document

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "07_SUBMISSION" / "_build_v3.1"
DOCS = BUILD / "docs"
FIGS = BUILD / "figs"

ms = Document(DOCS / "01_Manuscript.docx")
supp = Document(DOCS / "04_Supplementary.docx")

ms_paras = [p.text for p in ms.paragraphs]
supp_paras = [p.text for p in supp.paragraphs]
ms_tables = len(ms.tables)
supp_tables = len(supp.tables)

# split the manuscript at the Figure legends heading
i_leg = next(i for i, t in enumerate(ms_paras) if t.strip() == "Figure legends")
i_dec = next(i for i, t in enumerate(ms_paras) if t.strip() == "Declarations")
body = "\n".join(ms_paras[:i_leg])
legends = "\n".join(ms_paras[i_leg:i_dec])

print("=" * 78)
print("GUARD 3  carriers")
print("=" * 78)

fail = 0

# ------------------------------------------------------------------ figures
fig_files = sorted(p.stem for p in FIGS.glob("Fig*.pdf"))
print(f"figure files: {len(fig_files)}  ->  {', '.join(fig_files)}")

main_files = {f for f in fig_files if re.fullmatch(r"Fig\d", f)}
supp_files = {f for f in fig_files if re.fullmatch(r"FigS\d", f)}

cited_main = set(re.findall(r"Fig\.?\s*(\d)(?![0-9])", body + "\n" + legends))
cited_supp = set(re.findall(r"Fig\.?\s*S(\d)", body + "\n" + legends))
cited_supp |= set(re.findall(r"Figure S(\d)", body + "\n" + legends))

print(f"\nmain figures cited: {sorted(cited_main, key=int)}")
print(f"supp figures cited: {sorted(cited_supp, key=int)}")

for n in sorted(cited_main, key=int):
    if f"Fig{n}" not in main_files:
        print(f"  MISSING FILE  Fig{n}.pdf")
        fail += 1
for n in sorted(cited_supp, key=int):
    if f"FigS{n}" not in supp_files:
        print(f"  MISSING FILE  FigS{n}.pdf")
        fail += 1
for f in sorted(main_files):
    if f[3:] not in cited_main:
        print(f"  ORPHAN FILE   {f}.pdf (never cited)")
        fail += 1
for f in sorted(supp_files):
    if f[4:] not in cited_supp:
        print(f"  ORPHAN FILE   {f}.pdf (never cited)")
        fail += 1

# legend entries vs files
leg_ids = re.findall(r"^Fig\.\s*(S?\d)\.", legends, flags=re.M)
print(f"\nlegend entries: {leg_ids}")
for lid in leg_ids:
    f = f"Fig{lid}"
    if f not in fig_files:
        print(f"  LEGEND WITHOUT FILE  {f}")
        fail += 1
missing_leg = [f for f in fig_files if f[3:] not in leg_ids]
for f in missing_leg:
    print(f"  FILE WITHOUT LEGEND  {f}")
    fail += 1

# ------------------------------------------------------------------ tables
print(f"\nmain table objects: {ms_tables}   supp table objects: {supp_tables}")
cited_tab = set(re.findall(r"Table\s*(\d)(?![0-9])", body + "\n" + legends))
cited_stab = set(re.findall(r"Table\s*S(\d+)", body + "\n" + legends))
print(f"main tables cited: {sorted(cited_tab, key=int)}")
print(f"supp tables cited: {sorted(cited_stab, key=int)}")
main_heads = set(re.findall(r"^Table\s*(\d)\.", "\n".join(ms_paras), flags=re.M))
supp_heads = set(re.findall(r"^Table\s*S(\d+)\.", "\n".join(supp_paras), flags=re.M))
print(f"main table headings: {sorted(main_heads, key=int)}   supp table headings: {sorted(supp_heads, key=int)}")
for n in sorted(cited_tab, key=int):
    if n not in main_heads:
        print(f"  CITED BUT NO HEADING  Table {n}")
        fail += 1
for n in sorted(main_heads, key=int):
    if n not in cited_tab:
        print(f"  HEADING NEVER CITED   Table {n}")
        fail += 1
for n in sorted(cited_stab, key=int):
    if n not in supp_heads:
        print(f"  CITED BUT NO HEADING  Table S{n}")
        fail += 1
for n in sorted(supp_heads, key=int):
    if n not in cited_stab:
        print(f"  HEADING NEVER CITED   Table S{n}")
        fail += 1

# ------------------------------------------------------- declared in metrics
metrics = [t for t in ms_paras if t.startswith("Main text:")]
if metrics:
    m = metrics[0]
    print(f"\nmetrics line: {m}")
    for pat, real in [(r"Figures: (\d) main", len(main_files)),
                      (r"Figures: \d main, (\d) supplementary", len(supp_files)),
                      (r"Tables: (\d) main", ms_tables),
                      (r"Tables: \d main, (\d) supplementary", supp_tables),
                      (r"References: (\d+)", None)]:
        mm = re.search(pat, m)
        if mm and real is not None and int(mm.group(1)) != real:
            print(f"  METRIC WRONG  '{pat}' says {mm.group(1)}, files/objects = {real}")
            fail += 1

# ------------------------------------------------------------------ panels
print("\npanel letters per figure legend")
panel_map = {}
for lid in leg_ids:
    m = re.search(r"^Fig\.\s*%s\.(.*?)(?=^Fig\.|\Z)" % lid, legends, flags=re.M | re.S)
    if not m:
        continue
    letters = sorted(set(re.findall(r"\(([a-f])\)(?=\s+[A-Z(\u0394\u0391-\u03a9\u00b5])", m.group(1))))
    panel_map[lid] = letters
    print(f"  Fig.{lid}: {letters if letters else 'single panel'}")

# 正文里每个面板都要被点到（Fig. 3a、Fig. 3a,b、Fig. 3a–c 三种写法）
print("\npanel citations in the body")
for lid, letters in panel_map.items():
    if not letters:
        continue
    cited = set()
    for mm in re.finditer(r"Fig\.\s*%s((?:[a-f](?:[,\u2013-]\s*[a-f])*)+(?![0-9]))" % re.escape(lid), body):
        cited |= set(re.findall(r"[a-f]", mm.group(1)))
    uncited = [x for x in letters if x not in cited]
    print(f"  Fig.{lid}: panels {letters}  cited {sorted(cited)}"
          + (f"   UNCITED {uncited}" if uncited else ""))
    if uncited:
        fail += 1

print("\n" + ("PASS" if fail == 0 else f"FAIL  ({fail} findings)"))
sys.exit(0 if fail == 0 else 1)

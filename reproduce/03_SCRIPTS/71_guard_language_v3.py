#!/usr/bin/env python
"""Guards 2, 4, 5, 6, 7 — language (v3.0).

2 assertions   : every sentence that makes a strong claim must carry a number
                 or a figure/table pointer in the same sentence.
4 idle         : a sentence that says what we did must end with a number or a
                 carrier, otherwise it is decoration.
5 anti-defensive: defensive hedges are counted; more than 6 means rewrite.
6 reader       : abstract length, sentence length, opening and closing punch.
7 english      : sentence-length distribution against the target journal.
"""
import re
import statistics as st
import sys
from pathlib import Path

from docx import Document

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "07_SUBMISSION" / "_build_v3.1" / "docs"

ms = Document(DOCS / "01_Manuscript.docx")
paras = [p.text.strip() for p in ms.paragraphs if p.text.strip()]

i_abs = paras.index("Abstract")
i_bg = paras.index("Background", i_abs)
i_res = paras.index("Results", i_bg)
i_dis = paras.index("Discussion", i_res)
i_met = paras.index("Methods", i_dis)
i_leg = paras.index("Figure legends", i_met)

abstract = " ".join(paras[i_abs + 1:i_bg])
body_paras = paras[i_bg:i_leg]
body = " ".join(body_paras)

SENT = re.compile(r"(?<=[.!?])\s+(?=[A-Z(\u201c])")
HAS_EVID = re.compile(r"\d|Fig\.|Table\s*S?\d")


def sentences(t):
    t = t.replace("e.g. ", "eg ").replace("i.e. ", "ie ").replace("et al. ", "et al ")
    return [s.strip() for s in SENT.split(t) if len(s.strip()) > 25]


def para_sentences(ps):
    out = []
    for p in ps:
        out.extend(sentences(p))
    return out


sentences_of_body = para_sentences(body_paras)


fail = 0

# ------------------------------------------------------------------ 2 assertions
print("=" * 78)
print("GUARD 2  assertions")
print("=" * 78)
STRONG = ["same", "identical", "prove", "proves", "confirm", "confirms", "demonstrate",
          "establishes", "always", "never", "all ", "none ", "cannot", "no evidence",
          "exactly", "precisely", "definitively", "conclusively"]
naked = []
for s in sentences_of_body:
    low = s.lower()
    hit = [w for w in STRONG if w in low]
    if hit and not HAS_EVID.search(s):
        naked.append((hit, s))
for hit, s in naked:
    print(f"  [{','.join(hit)}] {s[:190]}")
print(f"  candidates without in-sentence evidence: {len(naked)}")

# ------------------------------------------------------------------ 4 idle
print("\n" + "=" * 78)
print("GUARD 4  idle sentences")
print("=" * 78)
SELF = re.compile(r"\b(we (took|computed|used|applied|measured|fit|fitted|estimated|"
                  r"partitioned|simulated|split|ran)|to (do|test) this|our approach)\b", re.I)
idle = [s for s in sentences_of_body if SELF.search(s) and not HAS_EVID.search(s)]
for s in idle:
    print(f"  {s[:190]}")
print(f"  self-referential sentences without a number or carrier: {len(idle)}")
fail += len(idle)

# ------------------------------------------------------------------ 5 defensive
print("\n" + "=" * 78)
print("GUARD 5  anti-defensive")
print("=" * 78)
HEDGE = ["should be noted", "may suggest", "it is important to", "we cannot exclude",
         "it should be emphasized", "it is worth noting", "further studies are needed",
         "future work", "limitation", "caveat", "we acknowledge", "arguably",
         "to some extent", "it is possible that", "we believe", "in our opinion",
         "of note", "importantly"]
allt = (abstract + " " + body + " " + " ".join(paras[i_leg:])).lower()
total = 0
for h in HEDGE:
    c = allt.count(h)
    if c:
        print(f"  {h!r}: {c}")
        total += c
print(f"  total defensive markers: {total}  (threshold 6)")
if total > 6:
    fail += 1

# ------------------------------------------------------------------ 6 reader
print("\n" + "=" * 78)
print("GUARD 6  reader (abstract)")
print("=" * 78)
aw = len(re.findall(r"[A-Za-z][A-Za-z'\u2019\-]*", abstract))
asents = sentences(abstract)
print(f"  abstract words: {aw}  (limit 250)")
print(f"  abstract sentences: {len(asents)}  max length {max(len(s.split()) for s in asents)}")
print(f"  opening: {asents[0][:150]}")
print(f"  closing: {asents[-1][:150]}")
print(f"  closing words: {len(asents[-1].split())}")
if aw > 250:
    fail += 1
if len(asents[-1].split()) > 30:
    print("  closing sentence is long; the last line should land in one breath")
    fail += 1

# ------------------------------------------------------------------ 7 english
print("\n" + "=" * 78)
print("GUARD 7  english")
print("=" * 78)
bs = sentences_of_body
lens = [len(s.split()) for s in bs]
print(f"  body sentences: {len(bs)}")
print(f"  mean length {st.mean(lens):.1f}  median {st.median(lens):.0f}  max {max(lens)}")
print(f"  over 40 words: {sum(1 for x in lens if x > 40)}  "
      f"over 50 words: {sum(1 for x in lens if x > 50)}")
CONN = ["however", "because", "therefore", "but ", "so ", "yet ", "although", "while ",
        "instead", "rather", "which ", "that "]
cs = sum(1 for s in bs if any(c in s.lower() for c in CONN))
print(f"  sentences carrying an explicit logical connective: {cs}/{len(bs)} "
      f"({100 * cs / len(bs):.0f}%)")
longest = sorted(bs, key=lambda s: -len(s.split()))[:3]
for s in longest:
    print(f"    {len(s.split()):>3}w  {s[:150]}")
if max(lens) > 55:
    print("  a sentence exceeds 55 words; split it")
    fail += 1

print("\n" + ("PASS" if fail == 0 else f"FAIL  ({fail} findings)"))
sys.exit(0 if fail == 0 else 1)

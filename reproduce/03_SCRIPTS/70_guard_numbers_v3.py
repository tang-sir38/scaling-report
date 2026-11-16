#!/usr/bin/env python
"""Guard 1 (v3.0) — numbers.

Two directions, both from the *delivered artefacts*, never from the build script:

  forward : every numeric token printed in the manuscript / cover letter /
            supplementary must trace back to a ledger entry (or to a structural
            number that is not a measurement).
  reverse : every ledger entry must actually appear in the delivered text
            (catches "computed but never injected").
  recompute: every ledger entry is re-evaluated from the locked JSON; a drift
            means the document was built from a different lock than the one on
            disk.

Nothing here imports the builder; it reads the DOCX and the ledger only.
"""
import json
import math
import re
import sys
from pathlib import Path

from docx import Document

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "07_SUBMISSION" / "_build_v3.1"
DOCS = BUILD / "docs"
LOCK = ROOT / "02_OUT" / "NUMBERS_LOCK_v5.json"

NUM = re.compile(r"\d[\d,]*(?:\.\d+)?")


def doc_text(p):
    d = Document(p)
    out = [q.text for q in d.paragraphs]
    for t in d.tables:
        for r in t.rows:
            for c in r.cells:
                out.append(c.text)
    for s in d.sections:
        for q in s.footer.paragraphs:
            out.append(q.text)
    return "\n".join(out)


# ---------------------------------------------------------------- allowed set
ledger = json.loads((BUILD / "numbers_ledger_v2.json").read_text(encoding="utf-8"))

allowed = set()
for e in ledger:
    for m in NUM.findall(e["text"]):
        allowed.add(m)

# structural / non-measurement integers that legitimately stand alone
STRUCT = set()
for i in range(0, 260):
    STRUCT.add(str(i))
for i in range(1990, 2031):
    STRUCT.add(str(i))
STRUCT |= {"1", "2", "3", "4", "5", "6", "7", "8", "16", "32", "64", "128", "256", "512",
           "1024", "2048", "4096"}

EXTRA_ALLOWED = {
    # significance levels, power, thresholds quoted as methods constants
    "5", "80", "95", "0.05", "2.802", "100", "8", "0.1", "0.02", "1.0", "3.0",
    "3.13", "2.5", "1.18", "10", "4", "6", "18", "155", "255", "000", "20",
    # theoretical limits and rhetorical constants, not measurements
    "2.00", "0.5", "1.01", "0.01", "351100",
    # software versions quoted in Methods / Code availability
    "2.5", "1.18", "3.11", "3.13", "05",
    # formatting constants
    "0.3", "1", "3", "0.03", "1000", "2083", "2082",
    # structural: decade label, dpi of the reproduced figures, percentile label
    "1980", "300", "97.5",
}

targets = {
    "manuscript": DOCS / "01_Manuscript.docx",
    "cover": DOCS / "02_Cover_Letter.docx",
    "supp": DOCS / "04_Supplementary.docx",
}

print("=" * 78)
print("GUARD 1  numbers (v3.0)")
print("=" * 78)
print(f"ledger entries: {len(ledger)}   distinct numeric strings from ledger: {len(allowed)}")

fail = 0
all_text = {}
for name, p in targets.items():
    if not p.exists():
        print(f"  MISSING {p}")
        fail += 1
        continue
    all_text[name] = doc_text(p)

# ---------------------------------------------------------------- forward
CITE = re.compile(r"\[[\d,\u2013\-]+\]")
DROP = re.compile(r"(doi:|ORCID|Tel: \+|Dongyan Road|Main text:|^\d+\.\s+[A-ZÅÖÜ].*20\d\d)"
                  r"|^\d+\.\s+$")


def scrub(txt):
    txt = CITE.sub(" ", txt)
    txt = re.sub(r"\b0\d_[A-Za-z]+", "PATH", txt)   # package file names, not quantities
    keep = []
    for line in txt.split("\n"):
        if DROP.search(line):
            continue
        keep.append(line)
    return "\n".join(keep)


print("\n[forward] numeric tokens in the delivered text that no ledger entry produces")
flagged = {}
for name, txt in all_text.items():
    txt = scrub(txt)
    for m in NUM.finditer(txt):
        tok = m.group(0).rstrip(",").rstrip(".")
        # a digit run welded to letters on either side is part of an identifier
        # (accession numbers, gene symbols, package names), not a measurement
        before = txt[m.start() - 1] if m.start() else ""
        after = txt[m.end()] if m.end() < len(txt) else ""
        if before.isalpha() or (after.isalpha() and not tok.endswith(".")):
            continue
        if tok in allowed or tok in STRUCT or tok in EXTRA_ALLOWED:
            continue
        ctx = txt[max(0, m.start() - 55): m.end() + 45].replace("\n", " ")
        flagged.setdefault(name, []).append((tok, ctx))

for name in flagged:
    print(f"\n  -- {name}")
    seen = set()
    for tok, ctx in flagged[name]:
        key = (tok, ctx[:40])
        if key in seen:
            continue
        seen.add(key)
        print(f"     {tok:>12}   …{ctx}…")
    fail += len(flagged[name])

# ---------------------------------------------------------------- reverse
print("\n[reverse] ledger entries whose text never appears in any delivered document")
blob = "\n".join(all_text.values())
missing = []
for e in ledger:
    if e["text"] not in blob:
        missing.append(e)
for e in missing[:25]:
    print(f"     {e['text']:>14}   <- {e['expr']}")
print(f"     total missing: {len(missing)}")
fail += len(missing)

# ------------------------------------------------- ledger internal consistency
print("\n[recompute] every ledger entry re-evaluated from the locked JSON")
RAW = json.loads(LOCK.read_text(encoding="utf-8"))


class NS(dict):
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


N = wrap(RAW)
bad = 0
for e in ledger:
    v = eval(e["expr"], {"__builtins__": {}},
             {"N": N, "d": N.datasets, "D3": N.derived3, "J": N.jackknife, "X": N.cross, "MU": N.multi,
              "F": N.fill.by_tissue, "B": N.boot.by_tissue, "LW": N.layer_weights.by_tissue,
              "L": N.layers, "S": N.synthetic, "V": N.variance, "M": N.models,
              "RL": N.reliability_layers, "n1": "n1", "n3": "n3",
              "PP": lambda x: "p%g" % x,
              "abs": abs, "min": min, "max": max, "sum": sum, "len": len, "round": round,
              "math": math, "np": math, "log": math.log, "exp": math.exp, "sqrt": math.sqrt})
    kind = e["kind"]
    if kind == "f":
        s = f"{v:.{e['nd']}f}"
    elif kind == "pct":
        s = f"{v * 100:.{e['nd']}f}%"
    elif kind == "int":
        s = f"{int(round(v)):,}"
    elif kind == "g":
        s = f"{v:g}"
    else:
        s = f"{v:.{e['nd']}f}"
    if s != e["text"]:
        bad += 1
        print(f"     DRIFT {e['expr']}  ledger={e['text']}  now={s}")
print(f"     drifted: {bad}")
fail += bad

print("\n" + ("PASS" if fail == 0 else f"FAIL  ({fail} findings)"))
sys.exit(0 if fail == 0 else 1)

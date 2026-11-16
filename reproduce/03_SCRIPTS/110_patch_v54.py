"""Build v5.4: make the package visible in the manuscript and add the per-gene
spatial-structure result to the Discussion.

Every number that appears in the new text is computed here from
NUMBERS_LOCK_v7.json (`spatial_proxy.by_tissue`) at build time, not typed in.
A hand-rounded value had already slipped through once (-0.63 for -0.6248), so
the quotes are generated from the lock and then re-parsed out of the finished
paragraph and compared back to the lock.

The linear partial correlation (-0.084 / +0.194 / -0.494) is deliberately not
quoted: 52_spatial_proxy.py states that a log-linear partial correlation is not
adequate for these two statistics, and reports abundance bands and a
semiparametric residual control instead. Both of those are used.
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import docx
from datetime import datetime

# python-docx writes the wall clock into core.xml on save, which makes two runs
# of this script produce different bytes from identical text. Freezing it is
# what makes the build reproducible; the date is the day the text was fixed.
FIXED_TS = datetime(2026, 11, 16, 0, 0, 0)


def freeze_zip(path, when=(2026, 11, 16, 0, 0, 0)):
    """Rewrite the docx container with fixed entry timestamps.

    A .docx is a zip, and zipfile stamps every entry with the wall clock, so
    freezing core.xml alone still leaves the bytes different between two runs.
    Same entry order, same payloads, fixed timestamps.
    """
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

PKG = Path(r"C:\Users\唐朝华\Desktop\23空间分辨率\A题_空间分辨率_投稿包_v5.3_20261116")
BUILD = Path(r"C:\Users\唐朝华\Desktop\scaling-report\reproduce\07_SUBMISSION\_build_v54")
DOCS = BUILD / "docs"
AUDIT = BUILD / "audit"
LOCK = Path(r"C:\Users\唐朝华\Desktop\scaling-report\reproduce\02_OUT\NUMBERS_LOCK_v7.json")

MANUSCRIPT = PKG / "02_主稿" / "01_Manuscript.docx"
TISSUES = ("colon", "lung", "brain")


def wordcount(text: str) -> int:
    return len(re.findall(r"[A-Za-z][A-Za-z'-]*", text))


def load_numbers() -> dict:
    d = json.loads(LOCK.read_text(encoding="utf-8"))["spatial_proxy"]["by_tissue"]
    out = {}
    for key in ("rho_gamma_moran", "band_k10", "semipar_pearson", "semipar_spearman"):
        vals = []
        for t in TISSUES:
            if key == "band_k10":
                v = [b for b in d[t]["band_rho"] if b["k"] == 10][0]["median"]
            elif key == "semipar_pearson":
                v = d[t]["semiparam"]["pearson"]
            elif key == "semipar_spearman":
                v = d[t]["semiparam"]["spearman"]
            else:
                v = d[t][key]
            vals.append(v)
        out[key] = vals
    out["raw"] = d
    return out


def fmt_vals(vals, digits: int) -> str:
    # the manuscript uses U+2212 throughout, not the ASCII hyphen that f-strings emit
    return (
        ", ".join(f"{v:+.{digits}f}" for v in vals)
        .replace("+", "")
        .replace("-", "\u2212")
    )


def build_edits(n: dict) -> list[dict]:
    moran = fmt_vals(n["rho_gamma_moran"], 3)
    band = fmt_vals(n["band_k10"], 2)
    semi = fmt_vals(n["semipar_pearson"], 2)

    return [
        dict(
            idx=56,
            mode="replace_tail",
            anchor=(
                "If that reading is right the exponent should track the dominant "
                "length scale of the tissue \u2014 a prediction this data set can "
                "suggest but not test."
            ),
            new=(
                "If that reading is right the exponent should track the dominant "
                "length scale of the tissue. At the per-gene level it partly does, "
                "and the split is by tissue. The per-gene exponent correlates with "
                "Moran's I at " + moran + " in colon, lung and brain, and holding "
                "abundance fixed by stratifying genes into ten bins leaves " + band +
                ". The association survives in colon and in brain and nearly "
                "vanishes in lung, so the exponent carries a spatial component in "
                "the two sections whose covariance has a measurable structure, and "
                "mostly does not in the one where a few transcripts set the value. "
                "What this design still cannot do is test the prediction where it "
                "would bite, below 8 \u00b5m, which is a measurement we cannot make."
            ),
            label="Discussion: per-gene spatial result replaces 'suggest but not test'",
        ),
        dict(
            idx=81,
            mode="append",
            anchor=(
                "The local slope of \u03c1 is the negative finite difference of "
                "log \u03c1 between consecutive lags on the median profile."
            ),
            new=(
                " For the comparison in the Discussion, each gene's covariance "
                "range was fitted to its semivariogram on seven lags with an "
                "exponential model, Moran's I was computed on the 8 \u00b5m grid, "
                "and abundance was controlled by stratifying genes into ten "
                "equal-count bins and taking the median within-bin Spearman "
                "correlation. A log-linear partial correlation is not reported "
                "because both statistics are non-linear in log abundance; a "
                "semiparametric residual control gives " + semi + " for the same "
                "three sections."
            ),
            label="Methods: how the per-gene spatial quantities were computed",
        ),
        dict(
            idx=90,
            mode="append",
            anchor=(
                "code and derived tables are available at the URL in the Code "
                "availability statement."
            ),
            new=(
                " The same repository ships the reporting procedure as an "
                "installable package (scaling_report, Python \u22653.9, numpy only) "
                "that computes the pooled exponent and its per-gene decomposition "
                "from bin moments alone."
            ),
            label="Methods: name the installable package",
        ),
        dict(
            idx=103,
            mode="append",
            anchor=(
                "the derived tables are available at "
                "https://github.com/tang-sir38/scaling-report."
            ),
            new=(
                " The package installs with pip and ships tests and a worked "
                "example; running it on the bundled colon moments reproduces the "
                "pooled exponent, the per-gene distribution and the share of the "
                "exponent carried by the heaviest transcript."
            ),
            label="Code availability: package usage line",
        ),
    ]


def main() -> None:
    DOCS.mkdir(parents=True, exist_ok=True)
    AUDIT.mkdir(parents=True, exist_ok=True)

    numbers = load_numbers()
    edits = build_edits(numbers)

    sha_before = hashlib.sha256(MANUSCRIPT.read_bytes()).hexdigest()
    doc = docx.Document(str(MANUSCRIPT))
    paras = doc.paragraphs
    before = [p.text for p in paras]

    print(f"source: {MANUSCRIPT}")
    print(f"  sha256 {sha_before[:16]}  paragraphs {len(paras)}")
    body_before = wordcount(" ".join(before[15:65]))
    print(f"  body wordcount {body_before}")

    applied = []
    for e in edits:
        p = paras[e["idx"]]
        assert e["anchor"] in p.text, f"anchor missing at [{e['idx']}]: {e['anchor'][:60]!r}"
        assert p.text.count(e["anchor"]) == 1, f"anchor not unique at [{e['idx']}]"
        assert p.text.endswith(e["anchor"]), (
            f"anchor not at tail of [{e['idx']}]; refusing to splice mid-sentence"
        )
        if e["mode"] == "replace_tail":
            assert len(p.runs) == 1, f"expected 1 run at [{e['idx']}], got {len(p.runs)}"
            p.runs[0].text = p.text[: -len(e["anchor"])] + e["new"]
        else:
            p.runs[-1].text = p.runs[-1].text + e["new"]
        applied.append(dict(idx=e["idx"], label=e["label"], words_added=wordcount(e["new"])))

    after = [p.text for p in doc.paragraphs]
    assert len(after) == len(before), "paragraph count changed"

    changed = [i for i, (a, b) in enumerate(zip(before, after)) if a != b]
    expected = sorted(e["idx"] for e in edits)
    assert changed == expected, f"changed {changed} != expected {expected}"

    body_after = wordcount(" ".join(after[15:65]))
    assert body_after <= 3333, f"body over limit: {body_after}"
    print(f"  body wordcount after {body_after}  (limit 3333)")

    out = DOCS / "01_Manuscript.docx"
    doc.core_properties.created = FIXED_TS
    doc.core_properties.modified = FIXED_TS
    doc.save(str(out))
    freeze_zip(out)
    assert len(docx.Document(str(out)).paragraphs) == len(after), "container rewrite broke the docx"
    sha_after = hashlib.sha256(out.read_bytes()).hexdigest()
    print(f"written {out}\n  sha256 {sha_after[:16]}")

    check = docx.Document(str(out))
    text_after = [p.text for p in check.paragraphs]
    assert text_after == after, "saved file does not round-trip"

    MINUS = "\u2212"
    disc = text_after[56]
    found_moran = re.findall(MINUS + r"\d\.\d{3}", disc)
    found_band = re.findall(MINUS + r"\d\.\d{2}(?!\d)", disc)
    assert len(found_moran) == 3, f"expected 3 Moran values, found {found_moran}"
    assert len(found_band) == 3, f"expected 3 band values, found {found_band}"
    for got, want in zip(found_moran, numbers["rho_gamma_moran"]):
        assert abs(float(got.replace("\u2212", "-")) - round(want, 3)) < 5e-4, (
            f"Moran mismatch: page {got} vs lock {want:.6f}"
        )
    for got, want in zip(found_band, numbers["band_k10"]):
        assert abs(float(got.replace("\u2212", "-")) - round(want, 2)) < 5e-3, (
            f"band mismatch: page {got} vs lock {want:.6f}"
        )

    meth = text_after[81]
    found_semi = re.findall(MINUS + r"\d\.\d{2}(?!\d)", meth)
    assert len(found_semi) == 3, f"expected 3 semiparam values, found {found_semi}"
    for got, want in zip(found_semi, numbers["semipar_pearson"]):
        assert abs(float(got.replace("\u2212", "-")) - round(want, 2)) < 5e-3, (
            f"semiparam mismatch: page {got} vs lock {want:.6f}"
        )

    for frag in ["falls 5- to 14-fold", "16,384-fold", "\u03b2 2\u201324 grid units"]:
        assert any(frag in t for t in text_after), f"lost fragment: {frag}"

    # ---- cover letter. This is the one place that still promised the code
    # "will be deposited", in the future tense, while the manuscript already
    # gives a live URL. An editor reads the cover letter first; the two files
    # must not disagree about whether the code exists yet.
    cover_src = PKG / "01_投稿信" / "00_Cover_Letter.docx"
    cover_sha_before = hashlib.sha256(cover_src.read_bytes()).hexdigest()
    cover_before = [p.text for p in docx.Document(str(cover_src)).paragraphs]
    cover = docx.Document(str(cover_src))
    cps = cover.paragraphs
    cover_anchor = (
        "Analysis code and the numbers file from which every reported value is "
        "injected will be deposited publicly and linked in the Code availability "
        "statement."
    )
    cover_new = (
        "Analysis code, the installable reporting package and the numbers file "
        "from which every reported value is injected are available at "
        "https://github.com/tang-sir38/scaling-report."
    )
    hit = [i for i, p in enumerate(cps) if cover_anchor in p.text]
    assert len(hit) == 1, f"cover anchor found in {hit}"
    ci = hit[0]
    assert cps[ci].text.count(cover_anchor) == 1, "cover anchor not unique"
    assert len(cps[ci].runs) == 1, f"expected 1 run at cover [{ci}], got {len(cps[ci].runs)}"
    cps[ci].runs[0].text = cps[ci].text.replace(cover_anchor, cover_new)

    cover_after = [p.text for p in cover.paragraphs]
    assert len(cover_after) == len(cover_before), "cover paragraph count changed"
    cover_changed = [i for i, (a, b) in enumerate(zip(cover_before, cover_after)) if a != b]
    assert cover_changed == [ci], f"cover changed {cover_changed} != [{ci}]"

    cover_out = DOCS / "00_Cover_Letter.docx"
    cover.core_properties.created = FIXED_TS
    cover.core_properties.modified = FIXED_TS
    cover.save(str(cover_out))
    freeze_zip(cover_out)
    assert len(docx.Document(str(cover_out)).paragraphs) == len(cover_after), "container rewrite broke the cover letter"
    cover_sha_after = hashlib.sha256(cover_out.read_bytes()).hexdigest()
    cover_text = "\n".join(cover_after)
    assert "will be deposited" not in cover_text, "future tense survived"
    assert "https://github.com/tang-sir38/scaling-report" in cover_text, "URL missing"
    print(f"written {cover_out}\n  sha256 {cover_sha_after[:16]}")

    rec = dict(
        change="v5.3 -> v5.4: package made visible; per-gene spatial result added; "
               "cover letter no longer promises future deposition",
        source_docx=str(MANUSCRIPT),
        docx_sha256_before=sha_before,
        docx_sha256_after=sha_after,
        paragraphs=len(after),
        changed_paragraphs=changed,
        body_words_before=body_before,
        body_words_after=body_after,
        body_words_limit=3333,
        edits=applied,
        numbers_source="NUMBERS_LOCK_v7.json :: spatial_proxy.by_tissue",
        numbers_quoted=dict(
            on_page=dict(moran=found_moran, band=found_band, semiparam=found_semi),
            from_lock=dict(
                rho_gamma_moran=[round(v, 3) for v in numbers["rho_gamma_moran"]],
                band_k10=[round(v, 2) for v in numbers["band_k10"]],
                semiparam_pearson=[round(v, 2) for v in numbers["semipar_pearson"]],
            ),
            roundtrip="page values re-parsed and compared to lock, tolerance 5e-4 / 5e-3",
        ),
        cover_letter=dict(
            source=str(cover_src),
            docx_sha256_before=cover_sha_before,
            docx_sha256_after=cover_sha_after,
            paragraph=ci,
            changed_paragraphs=cover_changed,
            before="... will be deposited publicly and linked in the Code availability statement.",
            after="... are available at https://github.com/tang-sir38/scaling-report.",
            why="cover letter said the code would be deposited, in the future tense, "
                "while the manuscript already gives a live URL; the editor reads the "
                "cover letter first",
        ),
        numbers_not_quoted=dict(
            linear_partial_rho_moran=[
                round(numbers["raw"][t]["partial_rho_moran"], 3) for t in TISSUES
            ],
            why="52_spatial_proxy.py: log-linear partial correlation is not adequate; "
                "abundance bands and a semiparametric residual control are reported instead",
        ),
    )
    (AUDIT / "patch_v5_4.json").write_text(
        json.dumps(rec, ensure_ascii=False, indent=1), encoding="utf-8"
    )
    print(f"audit  {AUDIT / 'patch_v5_4.json'}")


if __name__ == "__main__":
    main()

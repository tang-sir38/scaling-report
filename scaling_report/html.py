"""A one-file HTML report: the answer first, the evidence below it.

The Markdown report states what was measured. This one states what it means for
the reader's next decision, and puts the two pictures that carry the argument
above the fold:

* left — the exponent against bin size, with the textbook expectation
  (alpha falling as d^-2) drawn as a dashed line from the same starting point.
  The gap between the two lines is the whole paper.
* right — what a Poisson test would actually do at each scale. This is the panel
  a reader acts on, so it is drawn in linear units with the nominal 5% marked.

No dependency, no JavaScript, no external font. The SVG is written by hand so the
file opens from a USB stick on a machine with nothing installed.
"""
from __future__ import annotations

import html
import math
from pathlib import Path

# --------------------------------------------------------------------------
# palette — light background, dark ink, prints in greyscale without ambiguity
# --------------------------------------------------------------------------
INK = "#1a1a1a"
MUTED = "#5c6670"
RULE = "#d8dde3"
ACCENT = "#1f4e79"
WARN = "#a4262c"
OK = "#2e6b34"
PANEL = "#f6f8fa"
SERIF = 'Georgia, "Times New Roman", serif'
SANS = '-apple-system, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif'


# --------------------------------------------------------------------------
# small SVG helpers
# --------------------------------------------------------------------------


def _scale(v: float, d0: float, d1: float, r0: float, r1: float) -> float:
    if d1 == d0:
        return 0.5 * (r0 + r1)
    return r0 + (v - d0) / (d1 - d0) * (r1 - r0)


def _fmt_tick(v: float) -> str:
    if v >= 1:
        return f"{v:g}"
    if v >= 0.01:
        return f"{v:g}"
    return f"{v:.0e}".replace("e-0", "e-")


def _polyline(pts: list[tuple[float, float]]) -> str:
    return " ".join(f"{x:.1f},{y:.1f}" for x, y in pts)


def _fig_alpha(rep: dict, width: int = 330, height: int = 300) -> str:
    """log-log alpha against bin size, with the d^-2 expectation."""
    ds = [float(d) for d in rep["ladder_um"]]
    alpha = [float(a) for a in rep["alpha_pool"]]
    if not ds or not alpha:
        return ""

    xs = [math.log10(d) for d in ds]
    ys = [math.log10(max(a, 1e-9)) for a in alpha]

    # the expectation: alpha falls as d^-2 from the finest measured value
    exp_y = [ys[0] - 2.0 * (x - xs[0]) for x in xs]

    y_lo = min(min(ys), min(exp_y)) - 0.25
    y_hi = max(max(ys), max(exp_y)) + 0.25
    x_lo, x_hi = min(xs) - 0.12, max(xs) + 0.12

    L, R, T, B = 46, 12, 14, 38
    pw, ph = width - L - R, height - T - B

    def px(x):
        return _scale(x, x_lo, x_hi, L, L + pw)

    def py(y):
        return _scale(y, y_lo, y_hi, T + ph, T)

    out = [f'<svg viewBox="0 0 {width} {height}" width="100%" role="img" '
           f'aria-label="Overdispersion against bin size on log-log axes">']
    out.append(f'<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>')

    # gridlines at whole decades of alpha
    lo_i, hi_i = math.floor(y_lo), math.ceil(y_hi)
    for e in range(lo_i, hi_i + 1):
        y = py(e)
        if T - 1 <= y <= T + ph + 1:
            out.append(f'<line x1="{L}" y1="{y:.1f}" x2="{L+pw}" y2="{y:.1f}" '
                       f'stroke="{RULE}" stroke-width="1"/>')
            out.append(f'<text x="{L-6}" y="{y+3.5:.1f}" font-size="9" fill="{MUTED}" '
                       f'text-anchor="end" font-family="{SANS}">'
                       f'{_fmt_tick(10.0**e)}</text>')

    # x ticks at the measured scale values
    for d, x in zip(ds, xs):
        xp = px(x)
        out.append(f'<line x1="{xp:.1f}" y1="{T+ph}" x2="{xp:.1f}" y2="{T+ph+3}" '
                   f'stroke="{MUTED}" stroke-width="1"/>')
        out.append(f'<text x="{xp:.1f}" y="{T+ph+14}" font-size="9" fill="{MUTED}" '
                   f'text-anchor="middle" font-family="{SANS}">{d:g}</text>')

    # expected line
    out.append(f'<polyline points="{_polyline([(px(x), py(y)) for x, y in zip(xs, exp_y)])}" '
               f'fill="none" stroke="{MUTED}" stroke-width="1.4" stroke-dasharray="4 3"/>')

    # measured line + points
    out.append(f'<polyline points="{_polyline([(px(x), py(y)) for x, y in zip(xs, ys)])}" '
               f'fill="none" stroke="{ACCENT}" stroke-width="2.2" '
               f'stroke-linejoin="round"/>')
    for x, y in zip(xs, ys):
        out.append(f'<circle cx="{px(x):.1f}" cy="{py(y):.1f}" r="3" fill="{ACCENT}"/>')

    # axis labels
    out.append(f'<text x="{L+pw/2:.1f}" y="{height-4}" font-size="10" fill="{INK}" '
               f'text-anchor="middle" font-family="{SANS}">bin edge (µm)</text>')
    out.append(f'<text x="12" y="{T+ph/2:.1f}" font-size="10" fill="{INK}" '
               f'text-anchor="middle" font-family="{SANS}" '
               f'transform="rotate(-90 12 {T+ph/2:.1f})">overdispersion α</text>')

    # legend, in the empty upper-right corner: the measured line starts high on
    # the left and ends low on the right, and the expectation dives away, so the
    # only clear space is top-right
    lx, ly = L + pw - 100, T + 10
    out.append(f'<line x1="{lx}" y1="{ly}" x2="{lx+22}" y2="{ly}" stroke="{ACCENT}" '
               f'stroke-width="2.2"/>')
    out.append(f'<text x="{lx+27}" y="{ly+3.5}" font-size="9" fill="{INK}" '
               f'font-family="{SANS}">measured</text>')
    ly2 = ly + 14
    out.append(f'<line x1="{lx}" y1="{ly2}" x2="{lx+22}" y2="{ly2}" stroke="{MUTED}" '
               f'stroke-width="1.4" stroke-dasharray="4 3"/>')
    out.append(f'<text x="{lx+27}" y="{ly2+3.5}" font-size="9" fill="{MUTED}" '
               f'font-family="{SANS}">if α ∝ d⁻²</text>')

    out.append("</svg>")
    return "".join(out)


def _fig_fpr(adv: dict, width: int = 330, height: int = 300) -> str:
    """False-positive rate of a Poisson test at each scale."""
    table = adv.get("table") or []
    if not table:
        return ""
    ds = [r["bin_um"] for r in table]
    fpr = [r["fpr_median"] * 100.0 for r in table]

    L, R, T, B = 46, 12, 14, 38
    pw, ph = width - L - R, height - T - B

    y_hi = max(100.0, max(fpr) * 1.05)

    def px(i):
        if len(ds) == 1:
            return L + pw / 2
        return _scale(math.log10(ds[i]), math.log10(ds[0]), math.log10(ds[-1]),
                      L, L + pw)

    def py(v):
        return _scale(v, 0, y_hi, T + ph, T)

    out = [f'<svg viewBox="0 0 {width} {height}" width="100%" role="img" '
           f'aria-label="False-positive rate of a Poisson test by bin size">']
    out.append(f'<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>')

    for v in (0, 20, 40, 60, 80, 100):
        if v > y_hi:
            continue
        y = py(v)
        out.append(f'<line x1="{L}" y1="{y:.1f}" x2="{L+pw}" y2="{y:.1f}" '
                   f'stroke="{RULE}" stroke-width="1"/>')
        out.append(f'<text x="{L-6}" y="{y+3.5:.1f}" font-size="9" fill="{MUTED}" '
                   f'text-anchor="end" font-family="{SANS}">{v}</text>')

    # nominal level and the 10% working target, labelled once in the empty
    # upper-left corner rather than at the line ends, where the two labels
    # collide with each other and with the x-axis ticks
    for v, colour in ((5, OK), (10, WARN)):
        y = py(v)
        out.append(f'<line x1="{L}" y1="{y:.1f}" x2="{L+pw}" y2="{y:.1f}" '
                   f'stroke="{colour}" stroke-width="1" stroke-dasharray="5 3" '
                   f'opacity="0.8"/>')

    lgx, lgy = L + 8, T + 8
    out.append(f'<line x1="{lgx}" y1="{lgy}" x2="{lgx+16}" y2="{lgy}" '
               f'stroke="{OK}" stroke-width="1" stroke-dasharray="5 3"/>')
    out.append(f'<text x="{lgx+21}" y="{lgy+3.5}" font-size="8" fill="{OK}" '
               f'font-family="{SANS}">nominal 5%</text>')
    lgy2 = lgy + 12
    out.append(f'<line x1="{lgx}" y1="{lgy2}" x2="{lgx+16}" y2="{lgy2}" '
               f'stroke="{WARN}" stroke-width="1" stroke-dasharray="5 3"/>')
    out.append(f'<text x="{lgx+21}" y="{lgy2+3.5}" font-size="8" fill="{WARN}" '
               f'font-family="{SANS}">target 10%</text>')

    pts = [(px(i), py(v)) for i, v in enumerate(fpr)]
    out.append(f'<polyline points="{_polyline(pts)}" fill="none" stroke="{WARN}" '
               f'stroke-width="2.2" stroke-linejoin="round"/>')
    for x, y in pts:
        out.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="3" fill="{WARN}"/>')

    for i, d in enumerate(ds):
        xp = px(i)
        out.append(f'<line x1="{xp:.1f}" y1="{T+ph}" x2="{xp:.1f}" y2="{T+ph+3}" '
                   f'stroke="{MUTED}" stroke-width="1"/>')
        out.append(f'<text x="{xp:.1f}" y="{T+ph+14}" font-size="9" fill="{MUTED}" '
                   f'text-anchor="middle" font-family="{SANS}">{d:g}</text>')

    out.append(f'<text x="{L+pw/2:.1f}" y="{height-4}" font-size="10" fill="{INK}" '
               f'text-anchor="middle" font-family="{SANS}">bin edge (µm)</text>')
    out.append(f'<text x="12" y="{T+ph/2:.1f}" font-size="10" fill="{INK}" '
               f'text-anchor="middle" font-family="{SANS}" '
               f'transform="rotate(-90 12 {T+ph/2:.1f})">false positives (%)</text>')
    out.append("</svg>")
    return "".join(out)


# --------------------------------------------------------------------------
# page
# --------------------------------------------------------------------------


def _bar(fraction: float, colour: str) -> str:
    frac = max(0.0, min(1.0, fraction))
    return (f'<div class="bar"><div class="barfill" style="width:{frac*100:.1f}%;'
            f'background:{colour}"></div></div>')


def to_html(rep: dict, adv: dict, *, meta: dict | None = None,
            version: str = "") -> str:
    meta = meta or {}
    c = rep["concentration"]
    pg = rep["per_gene"]
    table = adv.get("table") or []
    sentences = adv.get("sentences") or []
    checklist = adv.get("reporting_checklist") or []

    title = html.escape(str(rep.get("section") or "section"))
    headline = html.escape(sentences[0]) if sentences else "No usable scale."
    detail = " ".join(html.escape(s) for s in sentences[1:])

    coarsest = table[-1] if table else None
    finest = table[0] if table else None

    rows = []
    for r in table:
        # one decimal: a value printed as "10%" next to a "no" verdict reads as a
        # contradiction, because the verdict is decided at 10.1%
        rows.append(
            f'<tr><td class="num">{r["bin_um"]:.0f}</td>'
            f'<td class="num">{r["alpha_pool"]:.3g}</td>'
            f'<td class="num">{r["phi_median"]:.1f}</td>'
            f'<td class="num {"bad" if r["fpr_median"] > 0.10 else "good"}">'
            f'{r["fpr_median"]*100:.1f}%</td>'
            f'<td class="num muted">{r["fpr_q25"]*100:.1f}–'
            f'{r["fpr_q75"]*100:.1f}%</td>'
            f'<td>{"yes" if r["acceptable"] else "no"}</td></tr>')

    rule = rep.get("reporting_rule") or []
    rule_html = "".join(f"<li>{html.escape(s)}</li>" for s in rule)

    usable = adv.get("max_usable_bin_um")
    usable_txt = f"{usable:.0f} µm" if usable else "none"

    fam = ""
    if meta.get("native"):
        fam = (f'<p class="meta">Read from Space Ranger: native bins '
               f'{", ".join(str(x) for x in meta["native"])} µm; ladder rebuilt to '
               f'{max(r["bin_um"] for r in table):.0f} µm.</p>')

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>scaling report · {title}</title>
<style>
  :root {{ color-scheme: light; }}
  * {{ box-sizing: border-box; }}
  body {{ margin: 0; padding: 32px 20px 56px; background: #ffffff; color: {INK};
         font-family: {SANS}; font-size: 15px; line-height: 1.6;
         -webkit-font-smoothing: antialiased; }}
  .wrap {{ max-width: 760px; margin: 0 auto; }}
  h1 {{ font-family: {SERIF}; font-size: 25px; line-height: 1.3; margin: 0 0 6px;
        font-weight: 600; letter-spacing: -0.01em; }}
  .sub {{ color: {MUTED}; font-size: 13px; margin: 0 0 26px; }}
  .verdict {{ border-left: 3px solid {ACCENT}; padding: 2px 0 2px 18px;
              margin: 0 0 10px; }}
  .verdict p {{ font-family: {SERIF}; font-size: 19px; line-height: 1.5; margin: 0; }}
  .detail {{ color: #333; font-size: 14px; margin: 14px 0 0; }}
  .cards {{ display: grid; grid-template-columns: repeat(3, 1fr); gap: 12px;
            margin: 30px 0 34px; }}
  .card {{ background: {PANEL}; border: 1px solid {RULE}; border-radius: 4px;
           padding: 14px 14px 12px; }}
  .card .k {{ font-size: 11px; letter-spacing: 0.05em; color: {MUTED};
              margin: 0 0 6px; }}
  .card .v {{ font-family: {SERIF}; font-size: 24px; line-height: 1.1; margin: 0;
              font-weight: 600; }}
  .card .n {{ font-size: 11px; color: {MUTED}; margin: 6px 0 0; }}
  h2 {{ font-size: 13px; letter-spacing: 0.05em; color: {MUTED}; font-weight: 600;
        margin: 38px 0 14px; padding-bottom: 7px;
        border-bottom: 1px solid {RULE}; }}
  .figs {{ display: grid; grid-template-columns: 1fr 1fr; gap: 22px; }}
  .fig {{ margin: 0; }}
  .fig figcaption {{ font-size: 12px; color: {MUTED}; margin-top: 4px;
                     line-height: 1.45; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 13px;
           font-variant-numeric: tabular-nums; }}
  th {{ text-align: left; font-weight: 600; font-size: 12px; letter-spacing: 0.02em;
        color: {MUTED}; padding: 0 10px 7px 0;
        border-bottom: 1px solid {RULE}; white-space: nowrap; }}
  th.num {{ text-align: right; }}
  td {{ padding: 8px 10px 8px 0; border-bottom: 1px solid #eef1f4; }}
  td.num {{ text-align: right; }}
  td:last-child {{ padding-left: 18px; white-space: nowrap; }}
  .good {{ color: {OK}; font-weight: 600; }}
  .bad {{ color: {WARN}; font-weight: 600; }}
  .muted {{ color: {MUTED}; }}
  ol {{ margin: 0; padding-left: 22px; }}
  ol li {{ margin: 0 0 9px; }}
  ul {{ margin: 0; padding-left: 20px; }}
  ul li {{ margin: 0 0 5px; font-size: 14px; }}
  .meta {{ font-size: 12px; color: {MUTED}; margin: 10px 0 0; }}
  footer {{ margin-top: 46px; padding-top: 16px; border-top: 1px solid {RULE};
            font-size: 11.5px; color: {MUTED}; }}
  @media (max-width: 640px) {{
    .cards {{ grid-template-columns: 1fr; }}
    .figs {{ grid-template-columns: 1fr; }}
    body {{ padding: 22px 14px 40px; }}
  }}
  @media print {{
    body {{ padding: 0; font-size: 11pt; }}
    .cards {{ break-inside: avoid; }}
    .figs {{ break-inside: avoid; }}
    @page {{ margin: 16mm; }}
  }}
</style>
</head>
<body>
<div class="wrap">

  <h1>What a Poisson test costs at this bin size</h1>
  <p class="sub">{title} · γ<sub>pool</sub> = {rep["gamma_pool"]:.3f} ·
     {len(rep["ladder_um"])} scales ·
     {pg["n"]:,} genes fitted</p>

  <div class="verdict"><p>{headline}</p></div>
  <p class="detail">{detail}</p>

  <div class="cards">
    <div class="card">
      <p class="k">Largest safe bin</p>
      <p class="v">{usable_txt}</p>
      <p class="n">keeps false positives under 10%</p>
    </div>
    <div class="card">
      <p class="k">False positives at {coarsest["bin_um"]:.0f} µm</p>
      <p class="v" style="color:{WARN}">{coarsest["fpr_median"]*100:.0f}%</p>
      <p class="n">nominal 5% · {coarsest["phi_median"]:.0f}× variance ratio</p>
    </div>
    <div class="card">
      <p class="k">Exponent carries</p>
      <p class="v">{c["n_eff_exponent"]:.1f}</p>
      <p class="n">effective genes; {html.escape(str(c["top1_gene"]))} supplies
         {c["top1_share"]*100:.0f}%</p>
    </div>
  </div>

  <h2>The two numbers that disagree</h2>
  <div class="figs">
    <figure class="fig">
      {_fig_alpha(rep)}
      <figcaption>Measured overdispersion against bin size. The dashed line is
      what pooling predicts if the Poisson term drove the decay (α ∝ d⁻², a
      16,384-fold fall from 2 to 256 µm). The measured fall is
      {rep["alpha_pool"][0]/rep["alpha_pool"][-1]:.1f}-fold.</figcaption>
    </figure>
    <figure class="fig">
      {_fig_fpr(adv)}
      <figcaption>What a test that assumes Poisson actually rejects, for the
      median of the most abundant genes. This panel is the one to act on.</figcaption>
    </figure>
  </div>

  <h2>Scale by scale</h2>
  <table>
    <thead><tr>
      <th class="num" style="text-align:right">Bin (µm)</th>
      <th class="num">α pooled</th>
      <th class="num">Variance ratio</th>
      <th class="num">False pos.</th>
      <th class="num">IQR</th>
      <th>Safe</th>
    </tr></thead>
    <tbody>{"".join(rows)}</tbody>
  </table>
  <p class="meta">Variance ratio is 1 + α·μ for the median of the 2,000 most
  abundant genes at each scale, the definition used in the paper's S2 Fig; false
  positives follow from it at a two-sided 5% level. IQR spans the 25th to 75th
  percentile of that same set, not a confidence interval.</p>

  {"<h2>Report these four things</h2><ol>" + rule_html + "</ol>" if rule_html else ""}

  {"<h2>Checklist before you write the sentence</h2><ul>" + "".join(f"<li>{html.escape(x)}</li>" for x in checklist) + "</ul>" if checklist else ""}

  {fam}

  <footer>
    scaling_report {html.escape(version)} · generator: <code>scaling-report</code>.
    Every number above is recomputable from the moments table; the raw matrix is
    not needed after the table exists. The variance ratio is a dispersion
    statement, not a confidence interval: bins are spatially correlated, so a
    standard error needs a spatial block bootstrap or at least two sections per
    group.
  </footer>

</div>
</body>
</html>
"""


def write_html(text: str, path: str | Path) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding="utf-8")
    return p

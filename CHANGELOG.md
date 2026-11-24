# Changelog

## 1.1.1 — 2026-11-23

The launcher worked only on a machine whose `python` was the one with numpy.
It does not any more. No change to any number, any estimator, or any output.

### Fixed

- **`run_windows.bat` gave up with "numpy is missing".** It called whatever
  `python` came first on PATH. On this machine that is a bare interpreter
  without numpy, so the launcher tried to `pip install`, failed, and closed.
  It now walks a list of interpreters and uses the first one that can
  `import numpy`, falling back to PATH only after every candidate is exhausted.
- **The launcher had no menu.** It asked for a folder and nothing else, so a
  reader who typed `1` got `1` as a path. It now offers [1] self-test,
  [2] demo, [3] my own data, [0] exit, and accepts the same choices as an
  argument: `run_windows.bat selftest`, `run_windows.bat demo`,
  `run_windows.bat "D:/path/to/outs"`.
- **The launcher was written with LF line endings.** `cmd.exe` mis-parses
  multi-line `if (...)` blocks without CRLF and exits silently, which is what
  a reader sees as a window that flashes and closes. It is CRLF now, and every
  conditional is a `goto` rather than a parenthesised block.
- **`if ... set ... & goto` on one line ran the `goto` unconditionally.** `&`
  binds to the whole line, not to the `if`. Split into separate lines.

### Changed

- Version strings in `pyproject.toml`, `.zenodo.json`, `CITATION.cff` and
  `scaling_report/__init__.py` are all 1.1.1. `.zenodo.json` and `CITATION.cff`
  had been left at 1.0.1 when 1.1.0 shipped.

## 1.1.0 — 2026-11-23

The version that a reader can use without writing the aggregation themselves.
Interface additions only: every v1.0 call still works and still returns the same
numbers.

### Added

- **`--from-10x DIR`** — reads a Space Ranger run directly. Native
  `square_002um`, `square_008um` and `square_016um` matrices are taken as they
  are (no aggregation, no extra memory) and the ladder is extended upward from
  the coarsest native scale. This removes the step that used to stop a reader
  with real data: the tool previously wanted a table of per-bin means and
  variances, which is the *output* of the hardest part.
- **`--from-mtx DIR --base-um N`** — one Matrix Market directory, ladder built
  upward. Pure numpy and the standard library: no scipy, no h5py, no anndata.
- **`--demo DIR`** — writes a synthetic section with an exponent known by
  construction (power-law intensity field, Poisson counts on top), laid out
  exactly like a Space Ranger output.
- **`--selftest`** — seven end-to-end checks, including recovery of the
  constructed exponent and the false-positive arithmetic against the paper's
  value.
- **`scaling_report.html`** — one file, hand-written SVG, no JavaScript, no
  external font, prints to A4. Two figures (the measured exponent against the
  `d^-2` expectation; the false-positive rate by bin size) and a per-scale table.
- **`advice.json`** — the per-scale table and the sentences, for scripting.
- **`scaling_report/advice.py`** — turns the exponent into a decision:
  `phi = 1 + alpha*mu`, the implied false-positive rate at a nominal 5%, the
  largest bin that keeps it under a target, and the sample-size multiplier.
- **`快速开始.md`** — a Chinese quick-start for readers who are not going to read
  the English README.
- **`run_windows.bat`** — drag a folder onto it, get a report. English interface,
  because a batch file with Chinese text inside garbles on a stock console.

### Fixed

- **The gene set behind phi.** `phi` is now the median over the **2,000 most
  abundant genes at each scale**, which is the definition in the paper's S2 Fig.
  v1.1 drafts used the median over every fitted gene, which is about five times
  smaller at 256 µm and made the tool disagree with the paper it accompanies
  (56% where the manuscript prints 72%). `advice.json` carries both numbers.
  Verified: on colon P1 the tool returns φ = 1.083 at 8 µm and 29.5 at 256 µm,
  i.e. 6.0% and 71.8% false positives, against the manuscript's 1.08, 30, 6.0%
  and 72%.
- `_block_sum` carried bin-level coordinates and entry-level gene indices as if
  they were the same length; they are not, and the first aggregation raised a
  broadcast error. Caught by `--selftest` before it reached a real matrix.
- Table headers no longer apply `text-transform: uppercase`, which rendered
  `µm` as `ΜM` and `α` as `Α` (Greek capitals), and the numeric headers are now
  right-aligned to their columns.
- False-positive rates print to one decimal, so a value displayed as `10%` is no
  longer paired with a `no` verdict decided at 10.1%.

### Not changed

- Every number the v1.0 package reported for `examples/colon_moments.npz` is
  unchanged: γ_pool 0.5547, per-gene median 0.428, IGKC 87.7%, n_eff 1.30.
  `tests/test_core.py` still passes 12/12.
- No confidence interval, and no absolute exponent for a tissue. Both refusals
  stand; the reasons are in the README.

## 1.0.1 — 2026-11-20

Zenodo archive; removed internal planning notes from the repository and pointed
the package at the submitted manuscript.

## 1.0.0 — 2026-11-15

First release: pooled exponent, per-gene distribution, effective gene counts,
and the identifiability bound.

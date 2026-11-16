r"""Core estimators.

Every definition here is the one used in the paper, restated once, in one file:

    alpha_g(d)      = (var_g(d) - mu_g(d)) / mu_g(d)^2        per-gene overdispersion
    alpha_pool(d)   = sum_g [var_g(d) - mu_g(d)] / sum_g mu_g(d)^2
                    = sum_g w_g(d) alpha_g(d),  w_g = mu_g^2 / sum mu^2
    gamma           = -d log(alpha) / d log(d)                 over the bin ladder
    v_g             = w_g alpha_g / alpha_pool                 share of the exponent
    v_bar_g         = mean_d v_g(d),  renormalised to sum 1
    n_eff           = 1 / sum_g weight_g^2

The distinction between w and v is the point of the paper: w is the weight of the
*variance*, v is the weight of the *exponent*, and only v answers "how many genes
decide the number I am reporting".
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import numpy as np

DEFAULT_DS = (2, 4, 8, 16, 32, 64, 128, 256)


class MomentsError(ValueError):
    """Raised when a moments table cannot be interpreted."""


@dataclass
class Moments:
    """Mean and variance of per-bin counts, one column per bin size."""

    ds: np.ndarray
    genes: np.ndarray
    mu: np.ndarray
    var: np.ndarray

    def __post_init__(self) -> None:
        self.ds = np.asarray(self.ds, dtype=float).ravel()
        self.genes = np.asarray(self.genes).astype(str).ravel()
        self.mu = np.asarray(self.mu, dtype=np.float64)
        self.var = np.asarray(self.var, dtype=np.float64)
        if self.mu.ndim == 1:
            self.mu = self.mu[:, None]
        if self.var.ndim == 1:
            self.var = self.var[:, None]
        if self.mu.shape != self.var.shape:
            raise MomentsError(f"mu {self.mu.shape} and var {self.var.shape} differ")
        if self.mu.shape[0] != self.genes.size:
            raise MomentsError(
                f"{self.genes.size} gene names for {self.mu.shape[0]} rows")
        if self.mu.shape[1] != self.ds.size:
            raise MomentsError(
                f"{self.mu.shape[1]} scale columns for {self.ds.size} bin sizes")
        if np.any(self.ds <= 0) or np.any(np.diff(self.ds) <= 0):
            raise MomentsError("bin sizes must be positive and strictly increasing")
        if not np.all(np.isfinite(self.var)) or np.any(self.var < 0):
            raise MomentsError("variances must be finite and non-negative")
        if not np.all(np.isfinite(self.mu)) or np.any(self.mu < 0):
            raise MomentsError("means must be finite and non-negative")

    # ---------------------------------------------------------------- readers
    @classmethod
    def from_npz(cls, path: str | Path) -> "Moments":
        """Read the native layout: ds, genes, mu_<d>, var_<d>."""
        z = np.load(Path(path), allow_pickle=True)
        keys = set(z.files)
        if "ds" not in keys or "genes" not in keys:
            raise MomentsError("npz must contain 'ds' and 'genes'")
        ds = np.asarray(z["ds"], dtype=float)
        mu = np.column_stack([z[f"mu_{int(d)}"] for d in ds])
        var = np.column_stack([z[f"var_{int(d)}"] for d in ds])
        return cls(ds=ds, genes=z["genes"], mu=mu, var=var)

    @classmethod
    def from_tidy(cls, path: str | Path, *, sep: str | None = None) -> "Moments":
        """Read a long table with columns gene, scale, mu, var.

        The separator is inferred from the suffix: .csv -> comma, anything else
        -> tab. Header names are matched case-insensitively and 'bin'/'bin_size'
        are accepted for 'scale'.
        """
        p = Path(path)
        if sep is None:
            sep = "," if p.suffix.lower() == ".csv" else "\t"
        with p.open("r", encoding="utf-8") as fh:
            header = fh.readline().rstrip("\n").split(sep)
        names = [h.strip().lower() for h in header]
        alias = {"bin": "scale", "bin_size": "scale", "d_um": "scale",
                 "mean": "mu", "variance": "var"}
        names = [alias.get(n, n) for n in names]
        need = {"gene", "scale", "mu", "var"}
        if not need.issubset(names):
            raise MomentsError(
                f"tidy table needs {sorted(need)}, found {sorted(names)}")
        data = np.genfromtxt(p, delimiter=sep, names=names,
                             dtype=None, encoding="utf-8", skip_header=1)
        if data.size == 0:
            raise MomentsError("tidy table is empty")

        genes, index = [], {}
        rows = []
        for rec in np.atleast_1d(data):
            g = str(rec["gene"])
            if g not in index:
                index[g] = len(genes)
                genes.append(g)
            rows.append((index[g], float(rec["scale"]), float(rec["mu"]),
                         float(rec["var"])))
        ds = np.array(sorted({r[1] for r in rows}), dtype=float)
        gi = np.array([r[0] for r in rows])
        si = np.searchsorted(ds, [r[1] for r in rows])
        mu = np.full((len(genes), ds.size), np.nan)
        var = np.full((len(genes), ds.size), np.nan)
        mu[gi, si] = [r[2] for r in rows]
        var[gi, si] = [r[3] for r in rows]
        if np.isnan(mu).any():
            raise MomentsError("tidy table is not complete: some gene has no row "
                               "at some scale")
        return cls(ds=ds, genes=np.array(genes), mu=mu, var=var)

    @classmethod
    def load(cls, path: str | Path) -> "Moments":
        p = Path(path)
        if p.suffix.lower() == ".npz":
            return cls.from_npz(p)
        return cls.from_tidy(p)

    # ---------------------------------------------------------------- helpers
    @property
    def n_genes(self) -> int:
        return int(self.mu.shape[0])

    @property
    def n_scales(self) -> int:
        return int(self.ds.size)

    def subset(self, mask: np.ndarray) -> "Moments":
        return Moments(self.ds, self.genes[mask], self.mu[mask], self.var[mask])


def _slope(ds: np.ndarray, y: np.ndarray) -> float:
    """-d log y / d log d, the least-squares exponent."""
    return -float(np.polyfit(np.log(ds), np.log(y), 1)[0])


@dataclass
class Fit:
    """Result of one decomposition."""

    gamma_pool: float
    alpha_pool: np.ndarray
    w_mean: np.ndarray            # (G,) w averaged over scales
    v_bar: np.ndarray             # (G,) share of the exponent, sums to 1
    gamma_g: np.ndarray           # (G,) per-gene exponent, NaN where not fitted
    ok: np.ndarray                # (G,) gene entered the per-gene fit
    mu_min: float
    sum_w2: float = float("nan")  # mean over scales of sum_g w_g^2
    ds: np.ndarray = field(default_factory=lambda: np.array(DEFAULT_DS, float))
    alpha_pool_loo: np.ndarray | None = None

    # ------------------------------------------------------------ summaries
    @property
    def n_fitted(self) -> int:
        return int(self.ok.sum())

    @property
    def gamma_weighted(self) -> float:
        """sum_g v_bar_g gamma_g; equals gamma_pool up to weight drift."""
        return float(np.nansum(self.v_bar * self.gamma_g))

    @property
    def median(self) -> float:
        return float(np.nanmedian(self.gamma_g[self.ok]))

    @property
    def q25(self) -> float:
        return float(np.nanpercentile(self.gamma_g[self.ok], 25))

    @property
    def q75(self) -> float:
        return float(np.nanpercentile(self.gamma_g[self.ok], 75))

    @property
    def n_eff_w(self) -> float:
        """Effective number of genes behind the *variance*."""
        return 1.0 / float(self.sum_w2)

    @property
    def n_eff_v(self) -> float:
        """Effective number of genes behind the *exponent*."""
        return 1.0 / float((self.v_bar ** 2).sum())

    @property
    def top1(self) -> int:
        return int(np.argmax(np.where(self.ok, self.v_bar, 0.0)))

    @property
    def top1_share(self) -> float:
        return float(self.v_bar[self.top1])

    def top_table(self, n: int = 10) -> list[tuple[int, float, float]]:
        order = np.argsort(-np.where(self.ok, self.v_bar, 0.0))[:n]
        return [(int(i), float(self.v_bar[i]), float(self.gamma_g[i]))
                for i in order]

    def share_top(self, k: int) -> float:
        v = np.where(self.ok, self.v_bar, 0.0)
        return float(np.sort(v)[::-1][:k].sum())


def fit(m: Moments, *, mu_min: float = 0.5,
        min_scales: int = 3) -> Fit:
    """Decompose the pooled exponent.

    Parameters
    ----------
    mu_min
        Gene filter: a gene enters the per-gene summary only if its mean count at
        the *largest* bin size is at least this. 0.5 is the paper's value.
    min_scales
        Minimum number of usable scales for a per-gene fit.
    """
    if m.n_scales < min_scales:
        raise MomentsError(f"{m.n_scales} scales is fewer than min_scales={min_scales}")

    mu, var, ds = m.mu, m.var, m.ds
    with np.errstate(invalid="ignore", divide="ignore"):
        alpha_g = (var - mu) / mu ** 2
    # a gene with no counts at a scale contributes nothing to either sum; writing
    # 0 instead of 0/0 keeps the pooled value finite without dropping the row
    alpha_g = np.where(mu > 0, alpha_g, 0.0)
    mu2 = mu ** 2
    total = mu2.sum(axis=0)
    if np.any(total <= 0):
        raise MomentsError("a scale has no counts at all")
    w = mu2 / total
    alpha_pool = (w * alpha_g).sum(axis=0)
    if np.any(~np.isfinite(alpha_pool)) or np.any(alpha_pool <= 0):
        raise MomentsError("pooled overdispersion is not positive at every scale")
    gamma_pool = _slope(ds, alpha_pool)

    finite = np.all(np.isfinite(alpha_g) & (alpha_g > 0), axis=1)
    ok = finite & (mu[:, -1] >= mu_min)
    gamma_g = np.full(m.n_genes, np.nan)
    usable = int(ok.sum())
    if usable:
        gamma_g[ok] = -np.polyfit(
            np.log(ds), np.log(alpha_g[ok]).T, 1)[0]

    with np.errstate(invalid="ignore", divide="ignore"):
        v = w * alpha_g / alpha_pool
    v_bar = np.where(ok, v.mean(axis=1), 0.0)
    s = v_bar.sum()
    if s <= 0:
        raise MomentsError("no gene carries a positive share of the exponent")
    v_bar = v_bar / s

    return Fit(gamma_pool=gamma_pool, alpha_pool=alpha_pool,
               w_mean=w.mean(axis=1), v_bar=v_bar, gamma_g=gamma_g,
               ok=ok, mu_min=mu_min, ds=ds,
               sum_w2=float((w ** 2).sum(axis=0).mean()))


def leave_one_out(m: Moments, fit_: Fit, gene: int) -> float:
    """gamma_pool after dropping one gene; the honest sensitivity to one gene."""
    mu = np.delete(m.mu, gene, axis=0)
    var = np.delete(m.var, gene, axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        a = (var - mu) / mu ** 2
    a = np.where(mu > 0, a, 0.0)
    w = mu ** 2 / (mu ** 2).sum(axis=0)
    return _slope(m.ds, (w * a).sum(axis=0))


def identifiability(fit_: Fit, *, lo: float = 0.0,
                    hi: float = 1.0) -> dict:
    """How far the per-gene median can move while gamma_pool stays fixed.

    Existence proof, not an optimisation. Put the heaviest gene at a and every
    other gene at one constant c. Then

        gamma_pool = v1 a + (1 - v1) c   =>   c = (gamma_pool - v1 a) / (1 - v1)

    and c is the median because it is the value of the vast majority. Hence
    d(median)/da = -v1/(1-v1), and c spans [lo, hi] only while a stays in
    [lo, hi], which is what bounds the admissible range of the median.
    """
    v1 = float(fit_.top1_share)
    if not np.isfinite(v1) or v1 >= 1.0:
        return dict(defined=False, reason="the top gene carries the whole exponent",
                    top1_share=v1)
    g = float(fit_.gamma_pool)
    a_lo = max((g - (1 - v1) * hi) / v1, lo)
    a_hi = min((g - (1 - v1) * lo) / v1, hi)
    med_hi = (g - v1 * a_lo) / (1 - v1)
    med_lo = (g - v1 * a_hi) / (1 - v1)
    return dict(
        defined=True,
        gamma_pool=g, top1_share=v1,
        anchor_range=[float(a_lo), float(a_hi)],
        median_at_fixed_pool_min=float(med_lo),
        median_at_fixed_pool_max=float(med_hi),
        median_span=float(med_hi - med_lo),
        d_median_d_anchor=float(-v1 / (1 - v1)),
        observed_median=None if np.isnan(fit_.median) else float(fit_.median),
        bounds=[lo, hi],
    )


def warnings_for(m: Moments, f: Fit) -> list[str]:
    out: list[str] = []
    if m.n_scales < 4:
        out.append(f"only {m.n_scales} bin sizes: the exponent is barely constrained")
    step = np.unique(np.round(np.diff(np.log2(m.ds)), 6))
    if step.size > 1:
        out.append("bin sizes are not a geometric ladder; the OLS slope is still "
                   "defined but is no longer a scale-free exponent")
    frac = f.n_fitted / max(m.n_genes, 1)
    if frac < 0.5:
        out.append(f"only {frac*100:.0f}% of genes gave a usable per-gene fit; the "
                   f"median describes that subset, not the transcriptome")
    if f.n_fitted < 200:
        out.append(f"{f.n_fitted} genes fitted: too few for a distribution")
    if f.top1_share > 0.5:
        out.append(f"one transcript ({m.genes[f.top1]}) carries "
                   f"{f.top1_share*100:.0f}% of the exponent; report that share or "
                   f"the number will not be comparable with another tissue")
    if not np.array_equal(m.ds, np.array(DEFAULT_DS, float)):
        out.append("bin sizes differ from 2-256 um; absolute values are not "
                   "comparable with the paper's tables")
    return out

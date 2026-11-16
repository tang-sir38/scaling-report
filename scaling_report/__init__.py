r"""scaling_report -- report a binned spatial transcriptomics scaling exponent so
that two sections can be compared.

Minimal use:

    python -m scaling_report moments.npz --out report/

The package reports the pooled exponent, the per-gene distribution, how many
genes actually decide the pooled number, and what the pooled number does *not*
identify. It is deliberately small: numpy only, no dependency on a spatial
analysis framework, because the estimator needs nothing but mean and variance
per gene at a ladder of bin sizes.
"""
from .core import (DEFAULT_DS, Fit, Moments, MomentsError, fit,
                   identifiability, leave_one_out, warnings_for)
from .report import (Z_ALPHA, Z_POWER, compare, mde, rule_text, summarise,
                     to_markdown, write)

__version__ = "1.0.0"

__all__ = [
    "DEFAULT_DS", "Fit", "Moments", "MomentsError", "fit", "identifiability",
    "leave_one_out", "warnings_for", "Z_ALPHA", "Z_POWER", "compare", "mde",
    "rule_text", "summarise", "to_markdown", "write", "__version__",
]

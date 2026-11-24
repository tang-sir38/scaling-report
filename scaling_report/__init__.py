r"""scaling_report -- report a binned spatial transcriptomics scaling exponent so
that two sections can be compared, and say what a Poisson test would cost.

Minimal use, from a Space Ranger run:

    python -m scaling_report --from-10x run/binned_outputs --out report/

Or from a table you already have:

    python -m scaling_report moments.npz --out report/

The package reports the pooled exponent, the per-gene distribution, how many
genes actually decide the pooled number, and what the pooled number does *not*
identify. Since v1.1 it also reads Space Ranger binned output directly, so the
ladder no longer has to be built by hand, and turns the exponent into the one
number a reader acts on: the false-positive rate of a Poisson test at each bin
size. Still numpy only.
"""
from .advice import fpr_under_poisson, phi_for_fpr, scale_table, verdict
from .core import (DEFAULT_DS, Fit, Moments, MomentsError, fit,
                   identifiability, leave_one_out, warnings_for)
from .report import (Z_ALPHA, Z_POWER, compare, mde, rule_text, summarise,
                     to_markdown, write)
from .tenx import ladder_from_one_matrix, ladder_from_run

__version__ = "1.1.2"

__all__ = [
    "DEFAULT_DS", "Fit", "Moments", "MomentsError", "fit", "identifiability",
    "leave_one_out", "warnings_for", "Z_ALPHA", "Z_POWER", "compare", "mde",
    "rule_text", "summarise", "to_markdown", "write",
    "fpr_under_poisson", "phi_for_fpr", "scale_table", "verdict",
    "ladder_from_run", "ladder_from_one_matrix",
    "__version__",
]

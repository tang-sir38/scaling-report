"""Point the repository README at the v5.5 title.

The README quotes the submitted title, and the title now says 15-fold. A public
repository that still says 14 while the manuscript says 15 is exactly the kind of
two-carrier disagreement the guard looks for inside the package — the difference
is that this one is outside the package, so nothing would have caught it.

Only the quoted title changes. Push still has to be done by hand: the CLI has no
credentials for this remote.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SUMMARY = REPO / "reproduce" / "02_OUT" / "extra_tissues_20261116.json"

OLD = ("transcriptomics falls 5- to 14-fold, not 16,384-fold, because it is set by "
       "tissue\ncovariance* (Tang et al., submitted)")


def main() -> int:
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    rng = summary["range"]
    new = (f"transcriptomics falls {rng['floor_min']}- to {rng['floor_max']}-fold, "
           f"not 16,384-fold, because the exponent is the decay rate of\n"
           "tissue covariance* (Tang et al., submitted)")

    p = REPO / "README.md"
    txt = p.read_text(encoding="utf-8")
    assert txt.count(OLD) == 1, f"README 标题锚点命中 {txt.count(OLD)} 次"
    txt = txt.replace(OLD, new)
    p.write_text(txt, encoding="utf-8")
    print(f"README title -> {rng['floor_min']}- to {rng['floor_max']}-fold")

    # 顺带把这一版的做法写进 README 末尾，读者才知道区间是几个组织撑起来的
    return 0


if __name__ == "__main__":
    sys.exit(main())

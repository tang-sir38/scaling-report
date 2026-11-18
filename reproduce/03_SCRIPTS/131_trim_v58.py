r"""131_trim_v58.py — 把正文压回 3,500 词以内（v5.8 的九片改写多出 6 词）。

只动三处，信息不丢：
1 [25] 六片那句改回"19.9 收尾"，并删掉重复的 "largest has a bootstrap interval of"
2 [67] 方法那句去掉 "contribute one section each"（与 each dataset one section 重复）
3 [46] 删 "from different patients"（方法段已写明，重复）
"""
from __future__ import annotations

from pathlib import Path

import docx

MAN = Path(r"C:\Users\唐朝华\Desktop\23空间分辨率\A题_空间分辨率_投稿包_v5.8_20261118"
           r"\02_主稿\01_Manuscript.docx")


def replace_in_para(p, old: str, new: str) -> bool:
    runs = p.runs
    text = "".join(r.text for r in runs)
    idx = text.find(old)
    if idx < 0:
        return False
    end = idx + len(old)
    pos, spans = 0, []
    for r in runs:
        spans.append((pos, pos + len(r.text)))
        pos += len(r.text)
    first = next(i for i, (a, b) in enumerate(spans) if b > idx)
    last = next(i for i, (a, b) in enumerate(spans) if b >= end)
    a0, _ = spans[first]
    if first == last:
        r = runs[first]
        r.text = r.text[: idx - a0] + new + r.text[end - a0:]
    else:
        runs[first].text = runs[first].text[: idx - a0] + new
        for i in range(first + 1, last):
            runs[i].text = ""
        runs[last].text = runs[last].text[end - spans[last][0]:]
    return True


EDITS = [
    ("[25] 六片句压缩",
     "Six further sections fall 19.9- and 13.1- (two further colon tumours), "
     "15.6- (small intestine), 3.5- (kidney), 2.3- (embryo) and 2.1-fold "
     "(pancreas). The largest has a bootstrap interval of 16.2 to 23.4, "
     "which does not overlap the 14.3 above.",
     "Six further sections fall 15.6- (small intestine), 3.5- (kidney), "
     "2.3- (embryo), 2.1- (pancreas), 13.1- and 19.9-fold (two further colon "
     "tumours). The 19.9-fold section's interval, 16.2 to 23.4, does not "
     "overlap the 14.3 above."),
    ("[67] 方法句压缩",
     "Three are colon adenocarcinomas from different patients; the other six "
     "tissues contribute one section each, and all nine enter the scaling "
     "analysis.",
     "Three are colon adenocarcinomas from different patients, six tissues "
     "contribute one section each, and all nine enter the scaling analysis."),
    ("[46] 删重复的 from different patients",
     "Two further colon tumours from different patients fall 13.1- and 19.9-fold",
     "Two further colon tumours fall 13.1- and 19.9-fold"),
]


def main() -> int:
    d = docx.Document(str(MAN))
    for label, old, new in EDITS:
        hits = [p for p in d.paragraphs if old in "".join(r.text for r in p.runs)]
        if len(hits) != 1:
            raise SystemExit(f"{label}: 命中 {len(hits)} 段，应为 1")
        if not replace_in_para(hits[0], old, new):
            raise SystemExit(f"{label}: 替换失败")
        print("  ", label)
    d.save(str(MAN))

    ps = [p.text for p in docx.Document(str(MAN)).paragraphs]
    body = " ".join(ps[15:65]).split()
    print(f"\n  正文 15:65 = {len(body)} / 3500   余量 {3500 - len(body)}")
    print(f"  摘要 = {len(ps[13].split())} / 150")
    assert len(body) <= 3500, len(body)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

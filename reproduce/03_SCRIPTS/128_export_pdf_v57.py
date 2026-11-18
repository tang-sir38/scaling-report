"""Export the three v5.5.1 docx files to PDF through Word, read-only.

Word is only used as a renderer: each document is opened with ReadOnly=True and
never saved, so the OOXML layer of the delivered docx is exactly what the patch
script wrote. ExportAsFixedFormat(path, 17) is the same call used for v5.4/v5.5.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import win32com.client as win32

PKG = Path(r"C:\Users\唐朝华\Desktop\23空间分辨率\A题_空间分辨率_投稿包_v5.7_20261118")
DOCS = [
    PKG / "01_投稿信" / "00_Cover_Letter.docx",
    PKG / "02_主稿" / "01_Manuscript.docx",
    PKG / "03_补充材料" / "02_Supplementary_Information.docx",
]
WD_EXPORT_PDF = 17


def main() -> int:
    word = win32.Dispatch("Word.Application")   # late binding: 避开 gencache 缓存损坏
    word.Visible = False
    word.DisplayAlerts = 0
    rc = 0
    try:
        for doc_path in DOCS:
            if not doc_path.exists():
                print(f"missing {doc_path}")
                rc = 1
                continue
            pdf = doc_path.with_suffix(".pdf")
            if pdf.exists():
                pdf.unlink()
            t0 = time.time()
            doc = word.Documents.Open(str(doc_path), ReadOnly=True,
                                      AddToRecentFiles=False, Visible=False)
            try:
                doc.ExportAsFixedFormat(str(pdf), WD_EXPORT_PDF, False, 0, 0, 0, 0,
                                        0, True, True, 0, True, True, False)
            finally:
                doc.Close(SaveChanges=0)
            print(f"  {pdf.relative_to(PKG)}  {pdf.stat().st_size:,} B  "
                  f"{time.time() - t0:.1f}s")
    finally:
        word.Quit()
    return rc


if __name__ == "__main__":
    sys.exit(main())

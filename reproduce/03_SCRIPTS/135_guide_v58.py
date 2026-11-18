r"""135_guide_v58.py — 投稿指南升到 v5.8（状态表、变更说明、checklist、脚本表）。

原则：表里的"包内文件"由实际文件数现场算，不手写。
"""
from __future__ import annotations

from pathlib import Path

PKG = Path(r"C:\Users\唐朝华\Desktop\23空间分辨率\A题_空间分辨率_投稿包_v5.8_20261118")
G = PKG / "00_投稿指南.md"

V58 = """## 零之二、v5.8 相对 v5.7 改了什么（六处硬伤 + 九片全入图）

1. **两处 `fifteen-fold` 改 `twenty-fold`**（Results 小标题、Discussion）。
   标题与摘要写 20 倍，正文却写"最多十五倍"——审稿人第一页就能撞见。
2. **摘要指数上界 `0.55` → `0.63`**。P5 的拟合指数 0.629 已取代 colon 的 0.555
   成为上界；投稿信里本来就写 0.63，只有摘要没跟着改。
3. **Fig 1 与 Fig S3 改九条曲线**。v5.7 把 P5（19.9）画进图、把 P1（13.1）留在
   图外，两者唯一区别是"P5 最大"，等于挑最大的一片画。现在同组织三片同色不同
   marker（P2 空心圆、P1 空心方、P5 实心），图上直接看到组织内 13.1–19.9 的极差。
4. **Fig S3 图注**：`eight` → `nine`；"27-fold" 与 "3.8-fold" 现场重算
   （26.93 / 3.848），数值不变。
5. **SI Note 2 稀释句**：`0.1 transcripts per bin` → "a tenth of its counts"，
   `0.001` → `0.0004`，出处由 Fig. 3d（合成场）改为 Methods（结肠实跑）。
6. **Discussion `a halving of depth by 0.0004` → `0.00004`**。0.0004 是稀释到
   10% 那次，减半（46%）只动 0.00004。

连带同步：摘要/引言/Results/Methods/Fig 1 注/Fig S3 注/投稿信/本指南的片数一律
`nine`；正文字数 3,499 → 3,498（余量 2）；Table S1 的 13.0 与 27.6 复核后**未改**
——那一列整列都是含线粒体的原始口径（colon 30.5 = 266.6M/8.73M 同理）。

守卫升到 **154 项**，新增两类：①**docx 级词形黑名单**（`fifteen-fold`、
`the eight curves`、`0.1 transcripts per bin`……只要出现就 FAIL，不再依赖"我记得
它会怎么写"）；②Fig 1b 的 P1 端点（13.097）纳入检查。变异测试：故意把
`twenty-fold` 改回 `fifteen`、把三处 `nine` 改回 `eight`、把 SI 的 0.0004 改
`0.001`，守卫抓到 7 项 FAIL。

---

"""


def main() -> int:
    s = G.read_text(encoding="utf-8")

    reps = [
        ("# 投稿指南 · A 题（空间分辨率 / 分箱标度律）v5.7",
         "# 投稿指南 · A 题（空间分辨率 / 分箱标度律）v5.8"),
        ("版本 v5.7 ｜ 2026-11-18", "版本 v5.8 ｜ 2026-11-18"),
        ("| 正文字数（Introduction→Discussion，含小节标题） | 3,499 | 3,500 | 余 1 |",
         "| 正文字数（Introduction→Discussion，含小节标题） | 3,498 | 3,500 | 余 2 |"),
        ("| 分析切片 | **8 片入图**（7 个组织 + 第三片结肠；包内共 9 个 Visium HD 切片 "
         "+ 1 个 Stereo-seq，P1 只作重复不入图） | — | — |",
         "| 分析切片 | **9 片全部入图**（7 个组织；结肠三片 P2/P1/P5 同色不同 marker） | — | — |"),
        ("守卫状态：**见 `05_Audit/guard_v5_7.json`**（v5.6 为 59/59，见 `guard_v5_6.json`）",
         "守卫状态：**见 `05_Audit/guard_v5_8.json`（154/154）**"
         "（v5.7 为 142/142，见 `guard_v5_7.json`；v5.6 为 59/59）"),
        ("`111_guard_v56.py` 在 v5.6 的文件上逐项实算，每个数字都从 `E:/xen_p2/out`\n"
         "的 JSON 现场重算，不沿用上一版任何结论。",
         "`111_guard_v58.py` 在 v5.8 的文件上逐项实算：数字从 `02_OUT` 与 "
         "`_build_v58/audit` 的 JSON 现场重算，内嵌图按字节比对，PDF 与 docx 双载体各查一遍。"),
        ("- [ ] `01_Manuscript.pdf` 是**最后一次改 docx 之后**导出的（守卫只对跑过的那份文件有效）",
         "- [ ] `01_Manuscript.pdf` 是**最后一次改 docx 之后**导出的（守卫只对跑过的那份文件有效；"
         "v5.8 换过内嵌图，PDF 已重导）"),
        ("- [ ] 6 张主图逐张上传，`Fig1.pdf`…`Fig6.pdf`",
         "- [ ] 6 张主图逐张上传，`Fig1.pdf`…`Fig6.pdf`\n"
         "- [ ] Fig 1 与 Fig S3 都是**九条**曲线（P2 空心圆 / P1 空心方 / P5 实心圆，同色）"),
    ]
    for old, new in reps:
        n = s.count(old)
        if n != 1:
            raise SystemExit(f"命中 {n} 次（应 1）：{old[:70]!r}")
        s = s.replace(old, new)

    anchor = "## 零之一、v5.7 相对 v5.6 改了什么"
    if s.count(anchor) != 1:
        raise SystemExit("找不到 零之一 锚点")
    s = s.replace(anchor, V58 + anchor)

    # 脚本对应关系补四行
    tail = "| `120_aggregate_extras.py` | 四片汇成 `extra_tissues` 块，写进锁；倍数与 CI 现场算 |"
    if s.count(tail) == 1:
        s = s.replace(tail, tail + "\n"
                      "| `125_fig_v58.py` | Fig 1 / Fig S3 九片版；fold 与跨度现场算，写 `fig_v58.json` |\n"
                      "| `130_p1_alpha.py` | 从 P1 逐基因矩重建 α(d)，四个锁定锚点回验 |\n"
                      "| `133_figembed_v58.py` | 按字节换 docx 内嵌的 Fig1 / FigS3，并同步 extent |\n"
                      "| `134_build_guard_v58.py` | 由 v5.7 守卫生成 v5.8（含词形黑名单），改动逐条断言 |")

    G.write_text(s, encoding="utf-8")

    # 包内文件数：现场数，写回表格
    n_files = sum(1 for p in PKG.rglob("*") if p.is_file() and p.name != "MANIFEST.json")
    s = G.read_text(encoding="utf-8")
    import re
    s2, k = re.subn(r"\| 包内文件 \| [^|]+ \|", f"| 包内文件 | {n_files}（以 MANIFEST.json 为准） |", s)
    if k != 1:
        raise SystemExit(f"包内文件行命中 {k} 次")
    G.write_text(s2, encoding="utf-8")
    print(f"guide updated; 包内文件 = {n_files}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

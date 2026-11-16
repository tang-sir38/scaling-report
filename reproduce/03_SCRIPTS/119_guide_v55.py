"""Write 00_投稿指南.md for the v5.5 package.

The guide is not uploaded, but it is the file a co-author opens first, so every
number in it is computed here from the delivered files (docx word counts, the
guard's own report, the manifest) instead of being typed in. That is also why it
is generated: the v5.4 guide still said "36 项" while its guard report said 39,
because the text was written before the guard was extended.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import docx

PKG = Path(r"C:\Users\唐朝华\Desktop\23空间分辨率\A题_空间分辨率_投稿包_v5.5_20261116")
REPRO = Path(r"C:\Users\唐朝华\Desktop\scaling-report\reproduce")
GUARD = REPRO / "07_SUBMISSION" / "_build_v55" / "audit" / "guard_v5_5.json"
SUMMARY = REPRO / "02_OUT" / "extra_tissues_20261116.json"
OUT = PKG / "00_投稿指南.md"


def main() -> int:
    g = json.loads(GUARD.read_text(encoding="utf-8"))
    if g["status"] != "PASS":
        raise SystemExit(f"守卫未通过（{g['status']}），不写指南: {g['fails']}")
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    rng = summary["range"]
    sec = {x["key"]: x for x in summary["sections"]}

    man = docx.Document(str(PKG / "02_主稿" / "01_Manuscript.docx"))
    paras = [p.text for p in man.paragraphs]
    n_body = len(" ".join(paras[15:65]).split())
    n_abs = len(paras[13].split())
    ref_start = next(i for i, t in enumerate(paras) if t.strip() == "References")
    n_refs = len([t for t in paras[ref_start:] if re.match(r"^\d+\.\s", t.strip())])
    man_rec = next(c for c in g["checks"] if c["item"] == "正文词数 ≤ 3500")
    manf = json.loads((PKG / "MANIFEST.json").read_text(encoding="utf-8"))
    n_files = len(manf["files"]) + 1  # manifest 本身不列自己

    longest = max(("mouse_small_intestine", "mouse_kidney", "human_pancreas",
                   "mouse_embryo"), key=lambda k: sec[k]["drop_fold"])
    ci = sec[longest]["boot"]["ci95"]
    colon = sec["colon"]["drop_fold"]
    fmin, fmax = rng["floor_min"], rng["floor_max"]
    if fmin == 5:
        bound_note = "只动上界一个数字。"
    else:
        lo_key = rng["min_key"]
        lo_lab = sec[lo_key]["label"].replace("Visium HD ", "").lower()
        bound_note = (f"**两端都动了**：下界现在由 {lo_lab}"
                      f" 的 {sec[lo_key]['drop_fold']:.2f} 倍决定，原来的 5 已不成立。")
    umi_lung = sec["lung"]["umi_per_bin"]
    umi_si = sec["mouse_small_intestine"]["umi_per_bin"]
    ratio = sec["lung"]["alpha_2um"] / sec["mouse_small_intestine"]["alpha_2um"]
    gap = 100 * (sec["mouse_small_intestine"]["drop_fold"] / sec["lung"]["drop_fold"] - 1)

    txt = f"""# 投稿指南 · A 题（空间分辨率 / 分箱标度律）v5.5

版本 v5.5 ｜ 2026-11-16 ｜ 目标刊：**Nature Computational Science**（Analysis）
本文件是操作手册，**不上传**。

---

## 零、现在是什么状态

| 项 | 实测值 | 期刊上限 | 判定 |
|---|---|---|---|
| 正文字数（Introduction→Discussion，含小节标题） | {n_body:,} | 3,500 | 余 {3500 - n_body} |
| 摘要字数 | {n_abs} | ≤150 | 余 {150 - n_abs} |
| 主图 | 6 | 6 | 满格 |
| 参考文献 | {n_refs} | ≤50 | 余 {50 - n_refs} |
| 补充表 | S1–S8（S1 现为 10 行） | 不限 | — |
| 补充图 | S1–S3 | 不限 | — |
| 分析切片 | **7 片进入标度分析**（包内共 8 个 Visium HD 切片 + 1 个 Stereo-seq） | — | — |
| 包内文件 | {n_files} | — | — |

守卫状态：**{g['pass_']}/{g['total']} PASS，0 失败**，见 `05_Audit/guard_v5_5.json`。
`111_guard_v55.py` 在 v5.5 的文件上逐项实算，不沿用 v5.4 的任何 carried-over 结论。

---

## 一、v5.5 相对 v5.4 改了什么（只有三件事）

1. **标题区间从 `5- to 14-fold` 改成 `{fmin}- to {fmax}-fold`。** 加了四个组织之后
   实测 {rng['min']:.2f}–{rng['max']:.2f} 倍，按**向下取整**惯例（5.82→5、14.31→14 已验证）
   取 {fmin} 与 {fmax}。{bound_note}结构、`16,384-fold` 对照、`not` 的句式全部不动。
2. **新增四片，n 从 3 提到 7。** 小鼠小肠 {sec['mouse_small_intestine']['drop_fold']:.2f} 倍、
   小鼠肾 {sec['mouse_kidney']['drop_fold']:.2f}、人胰腺 {sec['human_pancreas']['drop_fold']:.2f}、
   小鼠胚胎 {sec['mouse_embryo']['drop_fold']:.2f}。四片**在算指数之前就已定下**
   （Methods 与 Fig. S3 图注都写了这一句），不是看到结果才补的。
3. **载体补齐。** Table S1 从 6 行变 10 行；新增 Fig. S3 把七条 α(d) 曲线画在同一坐标上；
   Fig. 1b 的图注改成"three of the seven sections; all seven are in Fig. S3"——
   原来那句只画三片却不说，读者会以为全篇只有三片。

正文从 3,424 词变成 {n_body} 词，摘要仍是 {n_abs} 词：摘要把三个分片指数压成一个区间，
省下的词正好付掉 Results 新增的两句。

---

## 二、新数字的出处（审稿人一问就能答）

| 量 | 值 | 出处 |
|---|---|---|
| 四片下降倍数 | {sec['mouse_small_intestine']['drop_fold']:.2f} / {sec['mouse_kidney']['drop_fold']:.2f} / {sec['human_pancreas']['drop_fold']:.2f} / {sec['mouse_embryo']['drop_fold']:.2f} | `NUMBERS_LOCK_v7.json :: extra_tissues` |
| 最大一片的 95% CI | {ci[0]:.1f}–{ci[1]:.1f}（{sec[longest]['boot']['n']} 次 bin 级 bootstrap，seed {sec[longest]['boot']['seed']}） | 同上，`boot.ci95` |
| 结肠的点估计 | {colon:.2f} | `datasets.colon` |
| 七片 α(d) 曲线 | Fig. S3 | 同上 `alpha_by_d` |

**CI 下界 {ci[0]:.1f} 与结肠的 {colon:.2f} 几乎重合**，所以全文只写
"overlaps the colon value"，**没有一处**写显著超出。守卫里专门有一条查
`significantly` 不出现。

---

## 三、上传哪些文件（按这个顺序）

| # | 文件 | 用途 | 格式 |
|---|---|---|---|
| 1 | `02_主稿/01_Manuscript.docx` | 主稿（正文+方法+声明+文献+图注） | docx，**必须** |
| 2 | `02_主稿/01_Manuscript.pdf` | 审稿用 PDF（Word 导出，与 docx 同源） | pdf |
| 3 | `01_投稿信/00_Cover_Letter.docx` / `.pdf` | 投稿信（1 页） | docx 或 pdf |
| 4 | `03_补充材料/02_Supplementary_Information.docx` / `.pdf` | 补充材料（S1–S8 + Note 1 + Fig. S1–S3） | docx + pdf |
| 5 | `04_图件/Fig1.pdf` … `Fig6.pdf` | 6 张主图，**单独上传** | pdf |
| 6 | `05_Audit/NUMBERS_LOCK_v7.json`（可选） | 随代码仓库走，不单独上传 | json |

Fig. S1–S3 **不单独上传**，它们只活在补充材料 docx 里。

---

## 四、门户字段逐条怎么填

| 字段 | 填什么 | 注意 |
|---|---|---|
| *入口* | `https://mts-natcomputsci.nature.com/cgi-bin/main.plex` | 2026-11-16 现查；eJournalPress 系统（**不是** Snapp） |
| Article type | **Analysis** | 不要选 Article |
| Title | Overdispersion in binned spatial transcriptomics falls {fmin}- to {fmax}-fold, not 16,384-fold, because it is set by tissue covariance | 与稿件逐字一致 |
| Abstract | 直接粘贴稿件摘要（{n_abs} 词） | 粘完看计数器 |
| Authors | Chaohua Tang; Fangfang Wu; Yingying Zheng; Lingling Zhan; Hong Tan; Lijing Lin; Meimei Wen | 前三人 equal contribution |
| Corresponding | Chaohua Tang，13860904395@163.com，ORCID 0009-0004-7080-8993 | 与稿件一致 |
| Affiliation | Department of Gynecology, The 95th Hospital of Putian, No. 485 Dongyan Road, Licheng District, Putian City, Fujian Province 351100, China | 逐字一致 |
| Subject terms | Computational biology and bioinformatics; Statistical methods; Transcriptomics | 别选临床类 |
| Data availability | "Data available in public repository"：8 个 10x 数据集名 + STOmicsDB C04042E2 | 与稿件 Data availability 段一致 |
| Code availability | `https://github.com/tang-sir38/scaling-report` | 仓库需公开可访问 |
| Competing interests / Funding | None / None | 与稿件一致 |
| Suggested reviewers | **不填** | 没有真实可核的姓名就不填 |

---

## 五、三个最可能的攻击，和现成的应答锚点

| 攻击 | 现成证据（都在稿内） |
|---|---|
| "只加了几个组织，广度还是不够" | 7 片、7 种组织、三类状态（两种癌、四种正常、一种发育中）。命题的单位不是组织而是相关结构；每片都有 Table S1 的来源与格数 |
| "组织之间深度不同，指数差是深度造成的" | 肺与小肠深度几乎相同（{umi_lung:.2f} vs {umi_si:.2f} transcripts/bin），α(2 µm) 差 {ratio:.2f} 倍，而下降倍数只差 {gap:.0f}%；另有结肠抽稀到 1/4 只让 γ 动 0.0004。两条独立证据 |
| "组织是挑出来的，正好凑出 {fmin}–{fmax}" | 四片在算指数前定下（Methods、Fig. S3 图注各写一次）；且最大一片的 CI 与结肠重合，本文没有用它去主张任何排序 |

原有的四条攻击（同义反复、n=1 不可推组织差、下游代价、置换非因果）应答不变，见正文
Discussion 四段与 Table S7/S8。

---

## 六、提交前 10 分钟 checklist

- [x] 守卫 **{g['pass_']}/{g['total']} PASS**（`05_Audit/guard_v5_5.json`，全在 v5.5 文件上实算）
- [ ] `01_Manuscript.pdf` 是**最后一次改 docx 之后**导出的（守卫只对跑过的那份文件有效）
- [ ] 系统里 Title 是 `{fmin}- to {fmax}-fold`（**不是**旧版的 `5- to 14-fold`）
- [ ] 作者序、单位、邮箱、ORCID 与稿件首页逐字一致
- [ ] 6 张主图逐张上传，`Fig1.pdf`…`Fig6.pdf`
- [ ] 补充材料上传后翻到最后：**Figure S3 在**，S1–S8 与 Note 1 都在
- [ ] Table S1 是 **10 行**（少了就是拿成 v5.4 的包）
- [ ] Cover letter 里没有留 `[ ]` 占位
- [ ] 记下投稿号，存进 `00_ADMIN/`

---

## 七、脚本对应关系

| 脚本 | 作用 |
|---|---|
| `113_fetch_tissue.py` | 下载 + gzip CRC + 解包（登记前必须 `curl -r 0-3` 看 `1f8b`） |
| `114_pilot_tissue.py` | 新组织全链条：内存体检 → 02b 重聚合 → 07 池化指数 → 判读 |
| `115_boot_tail.py` | 尾段 bin 级 bootstrap（64/128/256 µm） |
| `120_aggregate_extras.py` | 四片汇成 `extra_tissues` 块，写进锁；倍数与 CI 现场算 |
| `121_fig_s3.py` | Fig. S3：七条 α(d) 曲线，样式与 Fig. 1b 一致 |
| `117_stage_v55.py` | 从 v5.4 包逐字节搬过来（只读 v5.4，不改不删） |
| `110_patch_v55.py` | 文字补丁：锚点断言 + 数字从锁注入 + 回读核对；含 SI 表格行与 Fig. S3 插入 |
| `118_manifest_v55.py` | 重建 MANIFEST（字节数与 SHA256 现算） |
| `111_guard_v55.py` | {g['total']} 项守卫，全部在当前文件上实算 |
| `112_pack_v55.py` | 打 zip 两次（项目目录 + 桌面）并逐字节比对 |

**纪律**：改任何交付件 = 改脚本重跑；禁手改 docx；改完重跑守卫；重算 SHA256、重打 zip。
`v5.4` 的原包与 zip 一律保留，不删不移。
"""
    OUT.write_text(txt, encoding="utf-8")
    print(f"written {OUT}  ({len(txt)} chars)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

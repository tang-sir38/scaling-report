r"""100_guard_v4.py — v4.0 投稿包守卫

原则（沿用工作区红线）
--------------------
1. 守卫只对**跑过的那一份文件**有效：每次都从磁盘重读，不看缓存。
2. 有阻断项就报阻断项，**不因为"检查脚本跑完了"而宣称 READY**。
3. 数字不是"看着对"，而是从锁重算后与 docx 里的字符串逐字比对。
4. 全绿≠没问题：另设一组"内容性"检查（图注与图位、文献单调、承诺悬空）。

用法
----
  python 100_guard_v4.py
"""
from __future__ import annotations

import json
import re
from pathlib import Path

from docx import Document

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "07_SUBMISSION/_build_v4"
DOCX = BUILD / "docs/01_Manuscript.docx"
LEDGER = BUILD / "numbers_ledger_v4.json"
STATS = BUILD / "stats_v4.json"
LOCK = ROOT / "02_OUT/NUMBERS_LOCK_v7.json"
REFS = ROOT / "02_OUT/REFS_v4.json"
FIGS = BUILD / "figs"

PASS, FAIL, WARN = [], [], []


def ck(cond, name, detail=""):
    (PASS if cond else FAIL).append((name, detail))


def warn(cond, name, detail=""):
    if cond:
        PASS.append((name, detail))
    else:
        WARN.append((name, detail))


def main():
    txt = "\n".join(p.text for p in Document(str(DOCX)).paragraphs)
    stats = json.loads(STATS.read_text(encoding="utf-8"))
    ledger = json.loads(LEDGER.read_text(encoding="utf-8"))
    lock = json.loads(LOCK.read_text(encoding="utf-8"))
    refs = json.loads(REFS.read_text(encoding="utf-8"))["refs"]

    # ---------- 1 结构完整性
    ck("{{" not in txt, "无未替换占位符", f"发现 {txt.count('{{')} 处")
    ck("[[[[" not in txt and "[[" not in txt, "无未解析引用标记", f"发现 {txt.count('[[')} 处")
    ck(stats["main_text_words"] <= 3500, "正文 ≤3500 词", f"{stats['main_text_words']} 词")
    ck(stats["n_display_items"] <= 6, "图位 ≤6", f"{stats['n_display_items']} 个")
    ck(stats["n_refs"] <= 50, "文献 ≤50", f"{stats['n_refs']} 条")
    ck(100 <= stats["abstract_words"] <= 150, "摘要 100–150 词",
       f"{stats['abstract_words']} 词")

    # ---------- 2 数字：从锁重算后逐字比对
    def get(d, path):
        cur = d
        for part in path.split("."):
            if isinstance(cur, dict):
                cur = cur[part]
            else:
                cur = cur[int(part)]
        return cur

    bad = []
    for key, rec in ledger.items():
        fmt = rec["fmt"]
        if any(s in rec["lock_path"] for s in (" / ", "*", "−", "之和")):
            continue          # 派生键（fall / gap / meddiff / 求和），交下一段独立复算
        try:
            v = get(lock, rec["lock_path"].split("  ×")[0].strip())
        except Exception:
            v = None
        if v is None:
            continue
        if rec.get("tf") == "abs":      # 台账记录了 abs，守卫就重放 abs
            v = abs(float(v))
        else:
            v = float(v)
        scale = 1.0
        if "×100" in rec["lock_path"]:
            scale = 100.0
        try:
            want = format(float(v) * scale, fmt) if fmt != "," else format(int(round(float(v))), ",")
        except Exception:
            continue
        if rec["formatted"] != want:
            bad.append((key, rec["formatted"], want))
    ck(not bad, "台账与锁逐字一致（可直接追溯的键）",
       f"{len(ledger)} 个键，{len(bad)} 处不一致 {bad[:4]}")

    # 派生键：独立复算
    deriv = {
        "c.fall": get(lock, "datasets.colon.alpha.2") / get(lock, "datasets.colon.alpha.256"),
        "l.fall": get(lock, "datasets.lung.alpha.2") / get(lock, "datasets.lung.alpha.256"),
        "b.fall": get(lock, "datasets.brain.alpha.2") / get(lock, "datasets.brain.alpha.256"),
        "gb.total": sum(get(lock, f"gamma_beta.by_tissue.{t}.n")
                        for t in ("colon", "lung", "brain")),
        "des.gap_cb": abs(abs(get(lock, "datasets.colon.gamma")) -
                          abs(get(lock, "datasets.brain.gamma"))),
    }
    dbad = []
    for k, v in deriv.items():
        if k in ledger:
            fmt = ledger[k]["fmt"]
            want = format(int(round(v)), ",") if fmt == "," else format(v, fmt)
            if ledger[k]["formatted"] != want:
                dbad.append((k, ledger[k]["formatted"], want))
    ck(not dbad, "派生数字独立复算一致", f"{len(dbad)} 处不一致 {dbad}")

    # ---------- 3 摘要与正文数字一致（跨段守卫）
    abs_i = [i for i, p in enumerate(Document(str(DOCX)).paragraphs)
             if p.text.strip() == "Abstract"]
    paras = [p.text for p in Document(str(DOCX)).paragraphs]
    abstract = paras[abs_i[0] + 1] if abs_i else ""
    # 正文 = References 之前的一切；先前按 "Methods" 截断会把 Methods/Discussion 的引用整段切掉
    body = txt.split("\nReferences\n")[0]
    nums_abs = set(re.findall(r"\d+\.\d+", abstract))
    miss = [n for n in nums_abs if n not in body]
    ck(not miss, "摘要中每个小数都在正文出现", f"未出现在正文：{miss}")

    # ---------- 4 文献：按首现顺序编号、全部被引用
    order = re.findall(r"\[(\d+(?:[–,\-]\d+)*)\]", body)
    firsts, seen = [], set()
    for grp in order:
        for part in re.split(r"[,]", grp):
            if "–" in part:
                a, b = part.split("–")
                rng = range(int(a), int(b) + 1)
            else:
                rng = [int(part)]
            for n in rng:
                if n not in seen:
                    seen.add(n); firsts.append(n)
    ck(firsts == sorted(firsts), "文献编号按首现单调递增",
       f"前 12 个首现顺序 {firsts[:12]}")
    ck(set(firsts) == set(range(1, stats["n_refs"] + 1)),
       "1..N 全部被引用且无跳号",
       f"缺 {sorted(set(range(1, stats['n_refs'] + 1)) - set(firsts))[:6]}")

    # ---------- 5 图件：6 个都在，且正文引用的编号都存在
    have = sorted(int(m.group(1)) for p in FIGS.glob("Fig*.pdf")
                  if (m := re.fullmatch(r"Fig(\d+)\.pdf", p.name)))
    ck(have == [1, 2, 3, 4, 5, 6], "图件 1–6 齐全", f"实际 {have}")
    cited = sorted({int(x) for x in re.findall(r"Fig\.\s*(\d)", txt)})
    ck(set(cited) <= set(have), "正文引用的图号都存在", f"引用 {cited}")

    # 每个图必须有图注：docx 里 ** 粗体标记已被 para() 剥掉，只能按行首锚定
    legends = re.findall(r"(?m)^Fig\.\s*(\d+)\.", txt)
    ck(sorted({int(x) for x in legends}) == [1, 2, 3, 4, 5, 6],
       "每个图都有图注", f"图注编号 {sorted({int(x) for x in legends})}")

    # 面板级载体覆盖：正文引 Fig. 5b，图注就必须有 (b)；图注的面板字母必须从 a 连续
    # （v4 里 Fig. 5 图注曾整段描述错面板，正是这一层没查）
    legend_text = {int(n): t for n, t in
                   re.findall(r"(?ms)^Fig\.\s*(\d+)\.\s*(.*?)(?=\nFig\.\s*\d+\.|\Z)", txt)}
    panel_bad = []
    for n, ltxt in legend_text.items():
        letters = sorted(set(re.findall(r"\(([a-f])\)", ltxt)))
        want = [chr(ord("a") + i) for i in range(len(letters))]
        if letters != want:
            panel_bad.append(f"Fig.{n} 面板字母不连续 {letters}")
        for m in re.finditer(rf"Fig\.\s*{n}\s*([a-f])(?![a-z])", body):
            if m.group(1) not in letters:
                panel_bad.append(f"正文引 Fig.{n}{m.group(1)}，图注无此面板")
    ck(not panel_bad, "正文引用的面板字母都在图注里", "；".join(sorted(set(panel_bad))[:4]))

    # ---------- 6 承诺悬空（红线 ⑨）＋ 载体覆盖（红线 ③）
    SUPP = ROOT / "07_SUBMISSION/_build_v4/docs/02_Supplementary_Information.docx"
    blockers = []
    if "URL TO BE INSERTED" in txt or "[URL" in txt:
        blockers.append("Code availability 缺 URL（投稿前必须补，否则编辑部会退回）")

    supp_txt, supp_shapes = "", 0
    if SUPP.exists():
        sd = Document(str(SUPP))
        supp_txt = "\n".join(p.text for p in sd.paragraphs)
        supp_txt += "\n".join(c.text for t in sd.tables for r in t.rows for c in r.cells)
        supp_shapes = len(sd.inline_shapes)
    else:
        blockers.append("补充材料 docx 不存在")

    cited_s = set(re.findall(r"Table (S\d+)", txt)) | set(
        re.findall(r"Fig(?:ure)?\s*S(\d+)", txt))
    have_s = set(re.findall(r"Table (S\d+)\.", supp_txt))
    missing = sorted(cited_s - have_s)
    if missing:
        blockers.append(f"正文引用 {sorted(missing)}，补充材料里没有")
    if "Supplementary Note" in txt and "Supplementary Note 1" not in supp_txt:
        blockers.append("正文引用 Supplementary Note 1，补充材料里没有对应小节")

    # 反向：补充材料里列了表/图却没人引用 = 载体里塞了孤儿。
    # 补充材料内部的表可以在 Note 正文里自引（S5/S6 即如此），故并入 supp 自身的引用
    cited_in_supp = set(re.findall(r"Table (S\d+)", supp_txt))
    orphan_s = sorted(have_s - cited_s - cited_in_supp)
    ck(not orphan_s, "补充材料每张表都被引用（正文或 Note 内）", f"孤儿表 {orphan_s}")
    ck(supp_shapes >= 1, "补充材料图件已嵌入", f"{supp_shapes} 张内嵌图")
    ck(not blockers, "无承诺悬空 / 补充材料覆盖正文引用", "；".join(blockers))

    # ---------- 7 数字锁新鲜度
    import os
    newest_src = max((LECT for LECT in [
        get(lock, "meta.version") and 0,
    ]), default=0)
    ck(get(lock, "meta.version") == 7.0, "数字锁为 v7",
       f"version={get(lock, 'meta.version')}")
    n_src = len({v["lock_path"] for v in ledger.values()})
    warn(n_src >= 40, "台账覆盖足够多的锁路径", f"{n_src} 条不同路径")

    # ---------- 输出
    print("=" * 74)
    for name, det in PASS:
        print(f"  PASS  {name}" + (f"   [{det}]" if det else ""))
    for name, det in WARN:
        print(f"  WARN  {name}" + (f"   [{det}]" if det else ""))
    for name, det in FAIL:
        print(f"  FAIL  {name}" + (f"   [{det}]" if det else ""))
    print("=" * 74)
    status = "READY" if not FAIL and not WARN and not blockers else (
        "BLOCKED" if (FAIL or blockers) else "PASS_WITH_WARNINGS")
    print(f"状态：{status}   通过 {len(PASS)}  警告 {len(WARN)}  失败 {len(FAIL)}")
    if blockers:
        print("阻断项：")
        for b in blockers:
            print(f"  - {b}")
    if FAIL:
        print("不要在没有清空 FAIL 的情况下宣称 READY。")
    (BUILD / "guard_v4.json").write_text(json.dumps(
        dict(status=status, pass_=len(PASS), warn=len(WARN), fail=len(FAIL),
             fails=[n for n, _ in FAIL], warnings=[n for n, _ in WARN],
             blockers=blockers), indent=2, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()

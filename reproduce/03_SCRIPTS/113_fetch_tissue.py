#!/usr/bin/env python3
r"""
113_fetch_tissue.py — 取一片 10x Visium HD 公开数据，可断点续传。

为什么需要这个脚本
------------------
本机到 cf.10xgenomics.com 只能走代理，单连接 300-730 KB/s。实测 8 路并行
1.70 MB/s、16 路 2.17 MB/s，所以下载必须分块并发，否则一片 6 GB 要三小时。
分块文件留在磁盘上，中断后重跑只补没下完的块 —— 这一点决定了这件事能不
能跨会话做完。

用法
----
  python 113_fetch_tissue.py --list
  python 113_fetch_tissue.py --tissue mouse_small_intestine
  python 113_fetch_tissue.py --tissue mouse_small_intestine --jobs 16
  python 113_fetch_tissue.py --tissue mouse_small_intestine --verify-only

产物
----
  <root>/RAW/<tissue>/parts/part.00000 ...    分块（保留，用于续传）
  <root>/RAW/<tissue>/<name>_binned_outputs.tar.gz   拼接后的原始包
  <root>/RAW/<tissue>/extracted/             解开后的 binned_outputs/
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tarfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# 公开数据集登记表。只登记实测过的数据：每条都发过 range 请求，前四字节必须是
# gzip 魔数 1f8b。2026-11-16 实测，以下三个名字返回 HTTP 200 但内容是 XML
# AccessDenied —— Visium_HD_Human_Kidney / Visium_HD_Human_Breast_Cancer /
# Visium_HD_Human_Prostate，HEAD 会骗人，所以登记前必须看字节。
TISSUES = {
    "mouse_small_intestine": dict(
        name="Visium_HD_Mouse_Small_Intestine",
        species="Mus musculus",
        tissue="small intestine",
        condition="normal",
        bytes=6034152093,
        why="与结肠同源、正常组织，是判断 alpha 反映组织类型还是恶性程度的最短路径",
    ),
    "mouse_kidney": dict(
        name="Visium_HD_Mouse_Kidney",
        species="Mus musculus",
        tissue="kidney",
        condition="normal",
        bytes=9474823404,
        why="正常组织，上皮/血管/肾小球结构分明，检验 alpha 是否跟着结构复杂度走",
    ),
    "human_pancreas": dict(
        name="Visium_HD_Human_Pancreas",
        species="Homo sapiens",
        tissue="pancreas",
        condition="normal",
        bytes=8046375772,
        why="人正常实体器官，腺泡结构规则，作为结肠之外的第二个正常人组织",
    ),
    "mouse_embryo": dict(
        name="Visium_HD_Mouse_Embryo",
        species="Mus musculus",
        tissue="whole embryo",
        condition="developing",
        bytes=9808678337,
        why="发育中组织，结构梯度最强，是对 alpha 最不利的检验场景",
    ),
    "mouse_brain": dict(
        name="Visium_HD_Mouse_Brain",
        species="Mus musculus",
        tissue="brain",
        condition="normal",
        bytes=4616239264,
        why="已有一片，重下只为跨平台一致性对照，不是新组织",
    ),
}

BASE = "https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/{name}/{name}_binned_outputs.tar.gz"
DEFAULT_ROOT = Path(r"E:\A_spatial_20261116")
CHUNK = 32 * 1024 * 1024


def url_for(key: str) -> str:
    return BASE.format(name=TISSUES[key]["name"])


def head_size(url: str) -> int | None:
    """问服务端要文件大小。走 curl，与下载同一条路径，避免两套代理行为。"""
    p = subprocess.run(["curl", "-sIL", url], capture_output=True)
    if p.returncode != 0:
        print(f"  HEAD 失败: rc={p.returncode}")
        return None
    hits = re.findall(r"(?im)^content-length:\s*(\d+)\s*$",
                      p.stdout.decode("utf-8", "replace"))
    return int(hits[-1]) if hits else None


def fetch_range(url: str, start: int, end: int, dest: Path, tries: int = 5) -> int:
    """下载 [start, end] 闭区间到 dest，重试若干次。返回字节数。

    走 curl 子进程而不是 urllib：一块 32 MiB 若整块读进内存，16 路并发就是
    512 MiB，而这台机器常年只剩 1.7 GB 可用。curl 自己写盘，进程内存恒定。
    """
    want = end - start + 1
    if dest.exists() and dest.stat().st_size == want:
        return want
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    for k in range(1, tries + 1):
        if tmp.exists():
            tmp.unlink()
        # Timeouts matter more than they look. On 2026-11-16 three chunks of the
        # mouse kidney archive (32 MiB each) stopped receiving bytes and curl sat
        # on the dead socket indefinitely — the pool waited on them for 35
        # minutes with nothing on disk changing. --speed-limit/--speed-time abort
        # a transfer that has gone idle, --max-time caps one attempt, and the
        # retry loop above then re-asks for the same range.
        p = subprocess.run(
            ["curl", "-sS", "--fail",
             "--connect-timeout", "30",
             "--max-time", "1800",
             "--speed-limit", "2048", "--speed-time", "90",
             "-o", str(tmp), "-r", f"{start}-{end}", url],
            capture_output=True,
        )
        if p.returncode == 0 and tmp.exists() and tmp.stat().st_size == want:
            tmp.replace(dest)
            return want
        msg = (p.stderr or b"").decode("utf-8", "replace").strip()[:160]
        if k == tries:
            raise IOError(f"块 {dest.name} 下载失败（{p.returncode}）：{msg}")
        print(f"    重试 {k}/{tries - 1} 块 {dest.name}: rc={p.returncode} {msg}", flush=True)
    raise RuntimeError("unreachable")


def download(key: str, root: Path, jobs: int) -> Path:
    info = TISSUES[key]
    url = url_for(key)
    d = root / "RAW" / key
    parts = d / "parts"
    parts.mkdir(parents=True, exist_ok=True)

    total = head_size(url) or info.get("bytes")
    if not total:
        raise SystemExit("拿不到文件大小，无法分块")
    print(f"文件 {info['name']}")
    print(f"大小 {total:,} 字节 ({total / 2**30:.2f} GiB)")

    n = (total + CHUNK - 1) // CHUNK
    print(f"分块 {n} 块 × 32 MiB，并发 {jobs}")
    todo = []
    for i in range(n):
        s = i * CHUNK
        e = min(s + CHUNK - 1, total - 1)
        p = parts / f"part.{i:05d}"
        if not (p.exists() and p.stat().st_size == e - s + 1):
            todo.append((i, s, e, p))
    print(f"待下 {len(todo)} 块，已完成 {n - len(todo)} 块")

    if todo:
        done = 0
        with ThreadPoolExecutor(max_workers=jobs) as ex:
            futs = {ex.submit(fetch_range, url, s, e, p): i for i, s, e, p in todo}
            for f in as_completed(futs):
                f.result()
                done += 1
                if done % 10 == 0 or done == len(todo):
                    print(f"  进度 {n - len(todo) + done}/{n} 块", flush=True)

    # 拼接
    out = d / f"{info['name']}_binned_outputs.tar.gz"
    if out.exists() and out.stat().st_size == total:
        print(f"已拼接: {out}")
        return out
    print("拼接…")
    with out.open("wb") as w:
        for i in range(n):
            w.write((parts / f"part.{i:05d}").read_bytes())
    size = out.stat().st_size
    if size != total:
        raise SystemExit(f"拼接后大小不符: {size} != {total}")
    print(f"拼接完成: {out} ({size:,} 字节)")

    # gzip 完整性：读一遍校验尾部 CRC
    import gzip
    print("校验 gzip 完整性（读一遍全部解压流，不做落盘）…")
    with gzip.open(out, "rb") as g:
        while g.read(64 * 1024 * 1024):
            pass
    print("gzip CRC 通过")
    return out


def extract(key: str, root: Path, arc: Path) -> Path:
    d = root / "RAW" / key / "extracted"
    marker = d / ".done"
    if marker.exists():
        print(f"已解包: {d}")
        return d
    d.mkdir(parents=True, exist_ok=True)
    print("解包…（这一步读写量大，需要几分钟）")
    with tarfile.open(arc) as t:
        t.extractall(d)
    marker.write_text("ok", encoding="utf-8")
    print(f"解包完成: {d}")
    return d


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue")
    ap.add_argument("--root", default=str(DEFAULT_ROOT))
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--verify-only", action="store_true")
    ap.add_argument("--extract", action="store_true")
    a = ap.parse_args()

    if a.list:
        print(f"{'key':24s} {'name':34s} {'state'}")
        for k, v in TISSUES.items():
            print(f"{k:24s} {v['name']:34s} {v['species']}, {v['condition']}")
        print("\n用 --tissue <key> 取一片。")
        return 0

    if not a.tissue:
        print("需要 --tissue，或用 --list 看清单", file=sys.stderr)
        return 2
    if a.tissue not in TISSUES:
        print(f"未登记的组织: {a.tissue}", file=sys.stderr)
        return 2

    root = Path(a.root)
    if a.verify_only:
        d = root / "RAW" / a.tissue
        arc = d / f"{TISSUES[a.tissue]['name']}_binned_outputs.tar.gz"
        if not arc.exists():
            print(f"没有包: {arc}", file=sys.stderr)
            return 1
        import gzip
        with gzip.open(arc, "rb") as g:
            while g.read(64 * 1024 * 1024):
                pass
        print(f"OK {arc.stat().st_size:,} 字节")
        return 0

    arc = download(a.tissue, root, a.jobs)
    if a.extract:
        extract(a.tissue, root, arc)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

# -*- coding: utf-8 -*-
r"""
60_gef_pull.py -- pull one HDF5 dataset out of a remote Stereo-seq GEF.

The GEF is 16 GB but the bin1 expression records are a single contiguous block
of 2.25 GB. We resolve its byte offset with a couple of range requests and then
fetch only that window, in parallel chunks, resumably.

Usage:
  python 60_gef_pull.py --url <gef> --group bin1 \
      --out 00_RAW/stereo/C04042E2_bin1.bin --workers 16
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import threading
import time
import urllib.request
from pathlib import Path

UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"


def head_size(url: str) -> int:
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=60) as r:
        return int(r.headers["Content-Length"])


class HTTPFile:
    """Minimal seekable reader, blocks fetched on demand (HDF5 metadata only)."""

    def __init__(self, url, block=1 << 20, maxblocks=64):
        self.url, self.block, self.maxblocks = url, block, maxblocks
        self.pos = 0
        self.size = head_size(url)
        self.n = 0

    def read(self, n=-1):
        if n < 0:
            n = self.size - self.pos
        end = min(self.pos + n, self.size)
        out = bytearray()
        while self.pos < end:
            a = self.pos
            b = min((a // self.block + 1) * self.block, end) - 1
            req = urllib.request.Request(
                self.url, headers={"User-Agent": UA, "Range": f"bytes={a}-{b}"}
            )
            with urllib.request.urlopen(req, timeout=90) as r:
                out += r.read()
            self.pos = b + 1
            self.n += 1
        return bytes(out)

    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else (self.pos + off if whence == 1 else self.size + off)
        return self.pos

    def tell(self):
        return self.pos


def resolve(url: str, group: str, probe_cache: Path):
    if probe_cache.exists():
        return json.loads(probe_cache.read_text())
    import h5py

    h = h5py.File(HTTPFile(url), "r")
    d = h[f"geneExp/{group}/expression"]
    meta = dict(offset=int(d.id.get_offset()), nbytes=int(d.nbytes),
                n_records=int(d.shape[0]), itemsize=int(d.dtype.itemsize),
                dtype=d.dtype.descr if hasattr(d.dtype, "descr") else str(d.dtype))
    probe_cache.parent.mkdir(parents=True, exist_ok=True)
    probe_cache.write_text(json.dumps(meta, indent=2))
    return meta


def fetch(url, a, b, path, lock, prog):
    have = path.stat().st_size if path.exists() else 0
    total = b - a + 1
    if have >= total:
        with lock:
            prog["done"] += have
        return
    pos = a + have
    for attempt in range(200):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA,
                                                       "Range": f"bytes={pos}-{b}"})
            with urllib.request.urlopen(req, timeout=90) as r, open(path, "ab") as f:
                while True:
                    c = r.read(1 << 16)
                    if not c:
                        break
                    f.write(c)
                    pos += len(c)
                    with lock:
                        prog["done"] += len(c)
            if pos > b:
                return
        except Exception:                                    # noqa: BLE001
            time.sleep(min(2 ** min(attempt, 5), 30))
            pos = a + (path.stat().st_size if path.exists() else 0)
            if pos > b:
                return
    raise RuntimeError(f"chunk {a}-{b} failed")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", required=True)
    ap.add_argument("--group", default="bin1")
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--chunk-mb", type=int, default=32)
    a = ap.parse_args()

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    cache = out.with_suffix(".probe.json")
    meta = resolve(a.url, a.group, cache)
    print(f"dataset {a.group}: offset {meta['offset']:,} nbytes {meta['nbytes']:,} "
          f"records {meta['n_records']:,}", flush=True)

    parts = out.with_suffix(out.suffix + ".parts")
    parts.mkdir(exist_ok=True)
    step = a.chunk_mb * 1024 * 1024
    base = meta["offset"]
    total = meta["nbytes"]
    jobs, i = [], 0
    while i < total:
        n = min(step, total - i)
        jobs.append((i, n, parts / f"{i:012d}.part"))
        i += n

    prog = {"done": 0}
    lock = threading.Lock()
    t0 = time.time()
    stop = threading.Event()

    def reporter():
        while not stop.is_set():
            d = prog["done"]
            el = max(time.time() - t0, 1)
            sp = d / el / 1e6
            eta = (total - d) / max(d / el, 1) / 60
            print(f"  {d/1e9:.2f}/{total/1e9:.2f} GB  {100*d/total:5.1f}%  "
                  f"{sp:.2f} MB/s  eta {eta:.0f} min", flush=True)
            stop.wait(30)

    th = threading.Thread(target=reporter, daemon=True)
    th.start()
    idx = 0
    lock_job = threading.Lock()

    def worker():
        nonlocal idx
        while True:
            with lock_job:
                if idx >= len(jobs):
                    return
                off, n, path = jobs[idx]
                idx += 1
            fetch(a.url, base + off, base + off + n - 1, path, lock, prog)

    ts = [threading.Thread(target=worker) for _ in range(a.workers)]
    for t in ts:
        t.start()
    for t in ts:
        t.join()
    stop.set()

    with open(out, "wb") as w:
        for _, n, p in jobs:
            with open(p, "rb") as r:
                shutil.copyfileobj(r, w, 1 << 22)
    got = out.stat().st_size
    print(f"merged {got:,} bytes (expected {total:,})", flush=True)
    if got != total:
        return 1
    shutil.rmtree(parts, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

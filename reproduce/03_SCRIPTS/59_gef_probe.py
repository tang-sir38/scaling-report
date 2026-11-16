# -*- coding: utf-8 -*-
r"""
59_gef_probe.py -- read a remote Stereo-seq GEF over HTTP range requests.

A 16 GB GEF does not need to be downloaded to be inspected or, if it already
carries pre-binned groups, to be analysed. HDF5 needs random access; this wraps
the URL in a block-cached read-only file object so h5py can open it directly.

Usage:
  python 59_gef_probe.py --url <gef url> [--blocks 256]
"""
from __future__ import annotations

import argparse
import io
import sys
import time
import urllib.request
from collections import OrderedDict

import h5py

UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"
BLOCK = 1 << 20


class HTTPFile(io.RawIOBase):
    """Read-only, seekable view of a remote object with an LRU block cache."""

    def __init__(self, url: str, block: int = BLOCK, maxblocks: int = 256):
        self.url = url
        self.block = block
        self.maxblocks = maxblocks
        self._cache: OrderedDict[int, bytes] = OrderedDict()
        self.pos = 0
        self.nreq = 0
        self.nbytes = 0
        req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=60) as r:
            self.size = int(r.headers["Content-Length"])

    # ---------------------------------------------------------------- helpers
    def _fetch(self, idx: int) -> bytes:
        if idx in self._cache:
            self._cache.move_to_end(idx)
            return self._cache[idx]
        a = idx * self.block
        b = min(a + self.block, self.size) - 1
        req = urllib.request.Request(
            self.url, headers={"User-Agent": UA, "Range": f"bytes={a}-{b}"}
        )
        for attempt in range(4):
            try:
                with urllib.request.urlopen(req, timeout=120) as r:
                    data = r.read()
                break
            except Exception:                                  # noqa: BLE001
                if attempt == 3:
                    raise
                time.sleep(2 * (attempt + 1))
        self._cache[idx] = data
        self.nreq += 1
        self.nbytes += len(data)
        if len(self._cache) > self.maxblocks:
            self._cache.popitem(last=False)
        return data

    # ------------------------------------------------------------------- io
    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return True

    def tell(self) -> int:
        return self.pos

    def seek(self, off: int, whence: int = 0) -> int:
        if whence == 0:
            self.pos = off
        elif whence == 1:
            self.pos += off
        elif whence == 2:
            self.pos = self.size + off
        else:
            raise ValueError(whence)
        return self.pos

    def readinto(self, b) -> int:
        n = len(b)
        if self.pos >= self.size:
            return 0
        end = min(self.pos + n, self.size)
        got = 0
        while self.pos < end:
            idx = self.pos // self.block
            off = self.pos - idx * self.block
            chunk = self._fetch(idx)
            take = min(len(chunk) - off, end - self.pos)
            b[got:got + take] = chunk[off:off + take]
            got += take
            self.pos += take
        return got


def describe(name, obj):
    if isinstance(obj, h5py.Dataset):
        print(f"  {name:<44s} {str(obj.shape):>18s}  {obj.dtype}", flush=True)
    else:
        print(f"  {name:<44s} {'group':>18s}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", required=True)
    ap.add_argument("--blocks", type=int, default=256)
    ap.add_argument("--max-depth", type=int, default=3)
    a = ap.parse_args()

    f = HTTPFile(a.url, maxblocks=a.blocks)
    print(f"remote size {f.size:,} bytes", flush=True)
    t0 = time.time()
    h = h5py.File(f, "r")
    print(f"opened in {time.time() - t0:.1f}s, {f.nreq} range requests", flush=True)

    print("--- root attrs ---", flush=True)
    for k, v in h.attrs.items():
        s = str(v)
        print(f"  {k}: {s[:200]}", flush=True)

    print("--- tree ---", flush=True)
    h.visititems(describe)
    print(f"--- done, {f.nreq} requests, {f.nbytes/1e6:.1f} MB fetched, "
          f"{time.time()-t0:.1f}s ---", flush=True)


if __name__ == "__main__":
    sys.exit(main())

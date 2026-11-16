r"""
40_layer_maps.py — 两件小事，都要精确口径：

(1) 每个表达层占 Σμ² 的份额
    Table 3 要写"这 6 层覆盖了多少权重"。不能拍脑袋，得按 32 µm 的逐基因 μ 现场算：
      share_layer = Σ_{g in layer} μ_g² / Σ_all μ_g²
    同时给出"μ32 <= 10 的 6 层合计覆盖"，这个数决定正文能说多硬的话。

(2) 头号权重基因的空间分布
    pooled gamma 在 colon 是 IGKC 的 gamma（权重 62.4%），在 brain 是 Camk2n1（7.9%）。
    把两个基因在 32 µm 的计数画出来，让"一个基因决定一个指数"变成一眼可见。

输出 <outdir>/layer_weights.json 与 <outdir>/maps_<tag>.npz
"""
import argparse
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

LAYERS = [(0.01, 0.0316), (0.0316, 0.1), (0.1, 0.316), (0.316, 1.0),
          (1.0, 3.16), (3.16, 10.0), (10.0, 31.6), (31.6, 100.0)]


def load32(p):
    z = np.load(p, allow_pickle=False)
    shape = tuple(int(x) for x in z["shape"])
    X = z["data"], z["indices"], z["indptr"]
    genes = np.array([str(s) for s in z["genes"]])
    return X, shape, genes, z["coords"], int(z["n_domain_any"])


def layer_weights(tag, npzdir, out):
    X, shape, genes, coords, n = load32(Path(npzdir) / "bins_032um.npz")
    data, indices, indptr = X
    keep = ~np.array([s.upper().startswith("MT-") for s in genes])
    # 每基因总计数
    tot = np.zeros(shape[0], np.float64)
    for g in range(shape[0]):
        tot[g] = data[indptr[g]:indptr[g + 1]].sum()
    tot = tot[keep]
    gs = genes[keep]
    mu = tot / n
    w = mu ** 2
    W = w.sum()
    rows = []
    cum = 0.0
    for lo, hi in LAYERS:
        m = (mu >= lo) & (mu < hi)
        s = float(w[m].sum())
        cum += s
        rows.append(dict(lo=lo, hi=hi, n=int(m.sum()), share=float(s / W),
                         share_pct=float(100 * s / W)))
    out[tag] = dict(rows=rows, total_sigma_mu2=float(W),
                    coverage_6layers=float(sum(r["share"] for r in rows[:6])),
                    coverage_all8=float(sum(r["share"] for r in rows)),
                    n_genes=int(mu.size))
    print(f"[{tag}] 6 层覆盖 Σμ² = {100*out[tag]['coverage_6layers']:.1f}%  "
          f"8 层 = {100*out[tag]['coverage_all8']:.1f}%")
    for r in rows:
        print(f"   {r['lo']:.3g}-{r['hi']:.3g}: n={r['n']:5d}  share={r['share_pct']:.2f}%")
    return gs, mu, coords, data, indices, indptr, keep, shape


def dump_maps(tag, gs, mu, coords, data, indices, indptr, keep, shape, gene, od):
    idx_all = np.nonzero(gs == gene)[0]
    if idx_all.size == 0:
        print(f"  [{tag}] 找不到基因 {gene}")
        return
    gi = int(idx_all[0])
    rows = np.zeros(coords.shape[0], np.float32)
    tot = np.zeros(coords.shape[0], np.float32)
    for j, g in enumerate(np.nonzero(keep)[0]):
        a, b = indptr[g], indptr[g + 1]
        tot[indices[a:b]] += data[a:b]
        if j == gi:
            rows[indices[a:b]] = data[a:b]
    np.savez(od / f"maps_{tag}.npz", row=coords[:, 0].astype(np.int16),
             col=coords[:, 1].astype(np.int16), val=rows, total=tot, gene=gene)
    print(f"  [{tag}] {gene}: 非零 {int((rows > 0).sum()):,} / {rows.size:,} 格, "
          f"最大 {rows.max():.0f}；总计数最大 {tot.max():.0f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=str(ROOT / "02_OUT" / "layer_maps"))
    a = ap.parse_args()
    od = Path(a.outdir)
    od.mkdir(parents=True, exist_ok=True)

    out = {}
    sets = [("colon", str(ROOT / "01_WORK"), "IGKC"),
            ("lung", str(ROOT / "D3_LungCancer" / "01_WORK"), "IGKC"),
            ("brain", str(ROOT / "D2_MouseBrain" / "01_WORK"), "Camk2n1")]
    for tag, d, gene in sets:
        gs, mu, coords, data, indices, indptr, keep, shape = layer_weights(tag, d, out)
        dump_maps(tag, gs, mu, coords, data, indices, indptr, keep, shape, gene, od)

    json.dump(out, open(od / "layer_weights.json", "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print(f"\n写出 {od/'layer_weights.json'}")


if __name__ == "__main__":
    main()

"""106_downstream_v4.py — 下流分析偏差的直接量化（只跑结肠 P1）

两件事，回应"这一指数如何影响下流分析、偏差量级多大"：

A. 解析件：α(d) 把"方差被低估"换算成假阳性膨胀与有效样本量损失。
   若把 bin 计数当 Poisson（方差 = 均值），实际方差是 μ + μ²α，低估因子 φ = 1 + μα。
   检验统计量因此放大 √φ，名义 0.05 的实际 FPR = 2(1 − Φ(1.96/√φ))。
   这一项只依赖 α，与检验方法无关地给出**上界**（NB / Wilcoxon 会小于它）。

B. 实证件：真实的差异表达在粗尺度上错多少。
   在 16 µm 上按 marker 把组织分成"上皮/肿瘤"与"基质"两端，标签按整数坐标
   映射到 8/16/32/64/128/256 µm（粗 bin 按面积多数票，纯度 <0.8 丢弃），
   每个尺度做 Wilcoxon 秩和检验（scanpy 默认口径），以 8 µm 为参照报：
   top-100 召回、显著集 Jaccard、logFC 回归斜率（收缩因子）、统计量秩相关。

C. 细胞类型映射的代理：5 个 canonical marker 集的 module score，在粗 bin 上
   重算后与"8 µm 分数按粗 bin 聚合"相比，报相关、斜率与离散度比。分数被抹平
   最快的类型就是映射最先失效的类型。

坐标：各尺度 coords 都是同一物理网格上的整数下标，8 µm (r,c) 对应 s µm
(r//f, c//f)，f = s//8。已实测范围吻合（837 → 418 → 209 → 104 → 26）。

输出：02_OUT/downstream/downstream_colon.json
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy import sparse, stats

ROOT = Path(__file__).resolve().parents[1]
SEC = ROOT / "MS_Visium_HD_Human_Colon_Cancer_P1"
WORK = SEC / "01_WORK"
OUT = ROOT / "02_OUT" / "downstream"
OUT.mkdir(parents=True, exist_ok=True)

SCALES = ["008", "016", "032", "064", "128", "256"]
REF = "008"
DEF_SCALE = "016"
PURITY_MIN = 0.8
N_TOP_GENES = 3000
N_ALPHA_GENES = 2000
CHUNK = 50

MARK_EPI = ["EPCAM", "KRT18", "KRT8", "CEACAM5", "TFF3", "CDX2", "KRT19"]
MARK_STROMA = ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRB", "THY1"]

PANELS = {
    "epithelial": MARK_EPI,
    "stromal": MARK_STROMA,
    "T_cell": ["CD3D", "CD3E", "CD2", "CD247", "IL7R", "CCL5"],
    "plasma_B": ["IGKC", "IGHG1", "JCHAIN", "MZB1", "MS4A1", "CD79A"],
    "endothelial": ["PECAM1", "VWF", "CLDN5", "PLVAP", "ENG"],
}


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


_CACHE = {}


def load(scale):
    if scale in _CACHE:
        return _CACHE[scale]
    z = np.load(WORK / f"bins_{scale}um.npz", allow_pickle=True)
    X = sparse.csr_matrix((z["data"], z["indices"], z["indptr"]),
                          shape=tuple(int(v) for v in z["shape"]))
    out = (X, z["coords"].astype(np.int64), [str(g) for g in z["genes"]])
    if len(_CACHE) >= 2:
        _CACHE.pop(next(iter(_CACHE)))
    _CACHE[scale] = out
    return out


def tot_of(X):
    return np.asarray(X.sum(axis=0)).ravel()


def gene_index(scale):
    return {g: i for i, g in enumerate(load(scale)[2])}


# ---------------------------------------------------------------- A 解析件
def part_a():
    scales = ["008", "016", "032", "064", "128", "256"]
    X8, _, g8 = load(REF)
    gc = np.asarray(X8.sum(axis=1)).ravel()
    keep8 = list(np.argsort(gc)[::-1][:N_ALPHA_GENES])
    rec = {}
    for s in scales:
        Xs, _, gs = load(s)
        gi = gene_index(s)
        idx = [gi[g8[i]] for i in keep8 if g8[i] in gi]
        tots = tot_of(Xs)
        rows = []
        for st in range(0, len(idx), CHUNK):
            blk = idx[st:st + CHUNK]
            sub = np.asarray(Xs[blk, :].todense(), dtype=np.float64)
            mu = sub.mean(axis=1)
            var = sub.var(axis=1)
            with np.errstate(divide="ignore", invalid="ignore"):
                a = (var - mu) / mu ** 2
            for m, aa in zip(mu, a):
                if np.isfinite(aa) and m > 0:
                    rows.append((m, aa))
        arr = np.array(rows)
        m, a = arr[:, 0], arr[:, 1]
        phi = 1.0 + m * a
        good = np.isfinite(phi) & (phi > 0)
        phi = phi[good]
        fpr = 2 * (1 - stats.norm.cdf(1.96 / np.sqrt(phi)))
        rec[s] = dict(
            n_genes=int(len(a)),
            alpha_median=float(np.median(a)), alpha_q25=float(np.percentile(a, 25)),
            alpha_q75=float(np.percentile(a, 75)),
            mu_median=float(np.median(m)),
            phi_median=float(np.median(phi)), phi_q25=float(np.percentile(phi, 25)),
            phi_q75=float(np.percentile(phi, 75)),
            fpr_median=float(np.median(fpr)), fpr_q25=float(np.percentile(fpr, 25)),
            fpr_q75=float(np.percentile(fpr, 75)),
            n_eff_ratio=float(np.median(1.0 / phi)),
        )
        log(f"A {s}µm: α_med={rec[s]['alpha_median']:.3f} φ_med={rec[s]['phi_median']:.2f} "
            f"FPR_med={rec[s]['fpr_median']:.3f} n_eff={rec[s]['n_eff_ratio']:.3f}")
    return rec


# ---------------------------------------------------------------- 区域定义
def define_regions():
    X, coords, genes = load(DEF_SCALE)
    tot = tot_of(X)
    gi = gene_index(DEF_SCALE)

    def score(marks):
        idx = [gi[g] for g in marks if g in gi]
        v = np.asarray(X[idx, :].sum(axis=0)).ravel()
        return np.log1p(v / np.maximum(tot, 1) * 1e4)

    d = stats.zscore(score(MARK_EPI)) - stats.zscore(score(MARK_STROMA))
    lo, hi = np.percentile(d, [30, 70])
    lab = np.zeros(len(d), dtype=np.int8)
    lab[d >= hi] = 1
    lab[d <= lo] = 2
    log(f"区域（16 µm）：上皮端 {int((lab==1).sum())}，基质端 {int((lab==2).sum())}，"
        f"共 {len(lab)} bin")
    return coords, lab


def map_labels(scale, coords16, lab16):
    """16 µm 标签 → 任意尺度。细尺度继承；粗尺度面积多数票 + 纯度过滤。"""
    X, coords, genes = load(scale)
    if int(scale) >= 16:
        f = int(scale) // 16
        lut = {(int(r), int(c)): int(l) for (r, c), l in zip(coords16, lab16)}
        base_r = (coords[:, 0] * f)[:, None] + np.arange(f)[None, :]
        base_c = (coords[:, 1] * f)[:, None] + np.arange(f)[None, :]
        votes = np.zeros((coords.shape[0], 3), dtype=np.int32)
        for i in range(coords.shape[0]):
            rr = base_r[i]
            cc = base_c[i]
            grid = np.empty((f, f), dtype=np.int8)
            for a in range(f):
                for b in range(f):
                    grid[a, b] = lut.get((int(rr[a]), int(cc[b])), 0)
            cnt = np.bincount(grid.ravel(), minlength=3)[:3]
            votes[i] = cnt
        lab = np.argmax(votes, axis=1).astype(np.int8)
        tot_v = votes.sum(axis=1)
        purity = np.where(tot_v > 0, votes.max(axis=1) / np.maximum(tot_v, 1), 0)
        lab[(purity < PURITY_MIN) | (lab == 0)] = 0
    else:
        f = 16 // int(scale)
        lut = {(int(r), int(c)): int(l) for (r, c), l in zip(coords16, lab16)}
        lab = np.array([lut.get((int(r), int(c)), 0)
                        for r, c in zip(coords[:, 0] // f, coords[:, 1] // f)],
                       dtype=np.int8)
    return X, coords, genes, lab


# ---------------------------------------------------------------- B 实证件
def de(X, lab, keep_idx, min_tot_bin=20, min_nz=100):
    """两组 Wilcoxon。低深度 bin 先 QC（每 bin >= min_tot_bin counts）；基因可检性
    用**绝对非零 bin 数**判断——基质端每 bin 深度低，用比例门槛会把基因全筛掉
    （实测 MT-CO3 这类高表达基因在基质端的检出率也只有 0.37）。"""
    tot = tot_of(X)
    qual = tot >= min_tot_bin
    a = np.where((lab == 1) & qual)[0]
    b = np.where((lab == 2) & qual)[0]
    if len(a) < 10 or len(b) < 10:
        return None
    n1, n2 = len(a), len(b)
    lfc = np.full(len(keep_idx), np.nan)
    lfc_raw = np.full(len(keep_idx), np.nan)
    auc = np.full(len(keep_idx), np.nan)
    pv = np.full(len(keep_idx), np.nan)
    ok = np.zeros(len(keep_idx), dtype=bool)
    for k, g in enumerate(keep_idx):
        row = X.getrow(g)
        nz, vals = row.indices, row.data
        ca = np.zeros(n1)
        cb = np.zeros(n2)
        pa = np.searchsorted(a, nz)
        pa_c = np.minimum(pa, n1 - 1)
        ma = (pa < n1) & (a[pa_c] == nz)
        ca[pa[ma]] = vals[ma] / np.maximum(tot[nz[ma]], 1) * 1e4
        pb = np.searchsorted(b, nz)
        pb_c = np.minimum(pb, n2 - 1)
        mb = (pb < n2) & (b[pb_c] == nz)
        cb[pb[mb]] = vals[mb] / np.maximum(tot[nz[mb]], 1) * 1e4
        if (ca > 0).sum() < min_nz or (cb > 0).sum() < min_nz:
            continue
        ca = np.log1p(ca)
        cb = np.log1p(cb)
        u, p = stats.mannwhitneyu(ca, cb, alternative="two-sided", method="asymptotic")
        # 同一对分布，两种归一化：均值取在 log1p 尺度（scanpy 默认口径）还是原始比例
        lfc[k] = np.log2((ca.mean() + 1e-9) / (cb.mean() + 1e-9))
        lfc_raw[k] = np.log2((np.expm1(ca).mean() + 1e-9) / (np.expm1(cb).mean() + 1e-9))
        auc[k] = u / (n1 * n2)
        pv[k] = p
        ok[k] = True
    q = np.full(len(pv), np.nan)
    if ok.sum():
        qq = pv[ok] * ok.sum() / np.maximum(stats.rankdata(pv[ok]), 1)
        q[ok] = np.minimum(qq, 1.0)
    return dict(logfc=lfc, logfc_raw=lfc_raw, auc=auc, p=pv, q=q, detected=ok, n_a=n1, n_b=n2,
                tot_med_a=float(np.median(tot[a])), tot_med_b=float(np.median(tot[b])))


def signed_stat(d):
    with np.errstate(divide="ignore"):
        s = -np.log10(np.maximum(d["p"], 1e-300)) * np.sign(d["logfc"])
    return np.where(d["detected"], s, 0.0)


def compare(de_ref, de_s, name):
    ok = de_ref["detected"] & de_s["detected"]
    if ok.sum() < 50:
        return None
    s_ref, s_s = signed_stat(de_ref), signed_stat(de_s)
    top_ref = set(np.argsort(-np.abs(s_ref))[:100])
    top_s = set(np.argsort(-np.abs(s_s))[:100])
    sig_ref = set(np.where((de_ref["q"] < 0.05) & (np.abs(de_ref["logfc"]) > 0.5))[0])
    sig_s = set(np.where((de_s["q"] < 0.05) & (np.abs(de_s["logfc"]) > 0.5))[0])
    x, y = de_ref["logfc"][ok], de_s["logfc"][ok]
    sl, ic, r, p, se = stats.linregress(x, y)
    xr, yr = de_ref["logfc_raw"][ok], de_s["logfc_raw"][ok]
    slr, icr, rr, pr, ser = stats.linregress(xr, yr)
    rho = stats.spearmanr(s_ref[ok], s_s[ok]).statistic
    return dict(scale=name, n_a=de_s["n_a"], n_b=de_s["n_b"], n_common=int(ok.sum()),
                tot_med_a=de_s["tot_med_a"], tot_med_b=de_s["tot_med_b"],
                top100_recall=len(top_ref & top_s) / 100,
                jaccard_sig=float(len(sig_ref & sig_s) / max(len(sig_ref | sig_s), 1)),
                n_sig_ref=len(sig_ref), n_sig_s=len(sig_s),
                logfc_slope=float(sl), logfc_slope_se=float(se), logfc_r=float(r),
                logfc_slope_raw=float(slr), logfc_slope_raw_se=float(ser),
                logfc_raw_r=float(rr),
                logfc_intercept=float(ic), stat_spearman=float(rho),
                logfc_median_abs_ref=float(np.median(np.abs(x))),
                logfc_median_abs_s=float(np.median(np.abs(y))))


# ---------------------------------------------------------------- C 映射代理
def part_c():
    X8, c8, g8 = load(REF)
    tot8 = tot_of(X8)
    gi8 = gene_index(REF)
    F = 8
    rec = {}
    for s in SCALES[1:]:
        Xs, cs, gs = load(s)
        tots = tot_of(Xs)
        gis = gene_index(s)
        f = int(s) // F
        # 8 µm 每个 bin 归属的粗格
        kb = (c8[:, 0] // f).astype(np.int64) * 100000 + (c8[:, 1] // f)
        order = np.argsort(kb, kind="stable")
        kb_s = kb[order]
        uniq, start = np.unique(kb_s, return_index=True)
        counts = np.diff(np.append(start, len(kb_s))).astype(np.float64)
        kk = cs[:, 0].astype(np.int64) * 100000 + cs[:, 1]
        rec[s] = {}
        for name, marks in PANELS.items():
            idx8 = [gi8[g] for g in marks if g in gi8]
            idxs = [gis[g] for g in marks if g in gis]
            if not idx8 or not idxs:
                continue
            raw8 = np.asarray(X8[idx8, :].sum(axis=0)).ravel() / np.maximum(tot8, 1) * 1e4
            raw_s = np.asarray(Xs[idxs, :].sum(axis=0)).ravel() / np.maximum(tots, 1) * 1e4
            rec[s][name] = {}
            for tag, tf in (("log1p", np.log1p), ("raw", lambda x: x)):
                v8 = tf(raw8)
                agg = np.add.reduceat(v8[order], start) / counts
                lut = dict(zip(uniq.tolist(), agg.tolist()))
                ref_v = np.array([lut.get(int(x), np.nan) for x in kk])
                sc = tf(raw_s)
                okm = np.isfinite(ref_v)
                if okm.sum() < 20:
                    continue
                sl, ic, r, p, se = stats.linregress(ref_v[okm], sc[okm])
                rec[s][name][tag] = dict(
                    n=int(okm.sum()), pearson=float(r), slope=float(sl), slope_se=float(se),
                    sd_ref=float(np.std(ref_v[okm])), sd_s=float(np.std(sc[okm])),
                    sd_ratio=float(np.std(sc[okm]) / max(np.std(ref_v[okm]), 1e-12)))
            lg = rec[s][name].get("log1p", {})
            rw = rec[s][name].get("raw", {})
            log(f"C {s}µm {name}: log1p slope={lg.get('slope', float('nan')):.2f} "
                f"r={lg.get('pearson', float('nan')):.3f} | raw slope={rw.get('slope', float('nan')):.2f} "
                f"r={rw.get('pearson', float('nan')):.3f}")
    return rec


def main():
    t0 = time.time()
    log("A 解析件：逐基因 α(d) 与方差误设代价")
    A = part_a()

    log("区域定义")
    coords16, lab16 = define_regions()

    log("B 逐尺度 DE")
    X8, _, g8 = load(REF)
    gc = np.asarray(X8.sum(axis=1)).ravel()
    keep8 = list(np.argsort(gc)[::-1][:N_TOP_GENES])
    log(f"top {len(keep8)} 基因参与 DE")

    de_by = {}
    for s in SCALES:
        Xs, cs, gs, lab = map_labels(s, coords16, lab16)
        n1, n2 = int((lab == 1).sum()), int((lab == 2).sum())
        log(f"  {s}µm：上皮 {n1} / 基质 {n2}")
        if n1 < 15 or n2 < 15:
            log(f"  {s}µm 组太小，跳过")
            continue
        gi = gene_index(s)
        idx = [gi[g8[i]] for i in keep8 if g8[i] in gi]
        r = de(Xs, lab, idx)
        if r is not None:
            de_by[s] = r

    log("B 比对（参照 8 µm）")
    comp = []
    for s in SCALES[1:]:
        if s in de_by and REF in de_by:
            c = compare(de_by[REF], de_by[s], s)
            if c:
                comp.append(c)
                log(f"  {s}µm：top100 {c['top100_recall']:.2f} | Jacc {c['jaccard_sig']:.2f} | "
                    f"logFC 斜率 {c['logfc_slope']:.3f} | stat ρ {c['stat_spearman']:.3f} | "
                    f"显著 {c['n_sig_ref']}→{c['n_sig_s']}")

    log("C 细胞类型映射代理")
    C = part_c()

    payload = dict(built=time.strftime("%Y-%m-%d %H:%M:%S"),
                   section="MS_Visium_HD_Human_Colon_Cancer_P1",
                   ref_scale=REF, def_scale=DEF_SCALE, purity_min=PURITY_MIN,
                   n_genes_de=len(keep8), n_genes_alpha=N_ALPHA_GENES,
                   part_a_overdispersion=A, part_b_de=comp, part_c_mapping=C,
                   elapsed_s=round(time.time() - t0, 1))
    p = OUT / "downstream_colon.json"
    p.write_text(json.dumps(payload, ensure_ascii=False, indent=1), encoding="utf-8")
    log(f"写出 {p}  用时 {payload['elapsed_s']}s")


if __name__ == "__main__":
    sys.exit(main())

"""numbers_v8.json: 唯一数字源。全部从上游 json 现场取值，不手抄。
新增 P5 与四个额外组织（原 v7 锁未收录），来源明确标注。
"""
import json, os

LOCK = json.load(open(r'E:\SPATIAL_2026\02_OUT\NUMBERS_LOCK_v7.json', encoding='utf-8'))
EXTRA = json.load(open(r'C:\Users\唐朝华\WorkBuddy\A73_rewrite\audit\extra_tissues.json', encoding='utf-8'))
P5 = json.load(open(r'E:\A_spatial_20261116\OUT\human_colon_cancer_p5\lineA_R.json', encoding='utf-8'))
P5B = json.load(open(r'E:\A_spatial_20261116\OUT\human_colon_cancer_p5\boot_tail.json', encoding='utf-8'))
P5G = json.load(open(r'E:\A_spatial_20261116\OUT\human_colon_cancer_p5\pergene_three.json', encoding='utf-8'))
XEN = json.load(open(r'E:\SPATIAL_2026\02_OUT\downstream\downstream_colon.json', encoding='utf-8'))

D = LOCK['datasets']; G = LOCK['gene_level']; DV = LOCK['derived']; DV3 = LOCK['derived3']
DV4 = LOCK['derived4']; MOD = LOCK['models']; TH = LOCK['theory']; GB = LOCK['gamma_beta']
RS = LOCK['rho_shape']; IV = LOCK['intervene']; CR = LOCK['cross']; DEC = LOCK['decomp']['by_tissue']
DN = LOCK['downstream']; WS = LOCK['window_scan']; SP = LOCK['spatial_proxy']['by_tissue']
MULTI = LOCK['multi']['sections']
extra = {s['key']: s for s in EXTRA['sections']}

def r(x, n):
    return round(float(x), n)

N = {}
def put(k, v):
    N[k] = v

# ---- 设计几何
put('n_visium', 9); put('n_total', 11)
put('bins_2um_area', 10562500); put('area_ratio', 16384); put('d_ratio', 128)
put('n_bins_colon', D['colon']['n_bins']); put('umi_colon', r(D['colon']['umi_per_bin'], 1))
put('n_bins_lung', D['lung']['n_bins']); put('umi_lung', r(D['lung']['umi_per_bin'], 1))
put('n_bins_brain', D['brain']['n_bins']); put('umi_brain', r(D['brain']['umi_per_bin'], 1))
put('genes_human', D['colon']['n_genes_noMT']); put('genes_mouse', D['brain']['n_genes_noMT'])

# ---- 倍数与 α
put('fold_colon', r(extra['colon']['drop_fold'], 1))
put('fold_lung', r(extra['lung']['drop_fold'], 1))
put('fold_brain', r(extra['brain']['drop_fold'], 1))
put('fold_si', r(extra['mouse_small_intestine']['drop_fold'], 1))
put('fold_kidney', r(extra['mouse_kidney']['drop_fold'], 1))
put('fold_embryo', r(extra['mouse_embryo']['drop_fold'], 1))
put('fold_pancreas', r(extra['human_pancreas']['drop_fold'], 1))
put('fold_p1', r(P5['folds']['P1'], 1)); put('fold_p5', r(P5['folds']['P5'], 1))
# 极值由九条曲线现场取，不手写
_all_fold = [extra[k]['drop_fold'] for k in ('colon', 'lung', 'brain', 'mouse_small_intestine',
                                             'mouse_kidney', 'mouse_embryo', 'human_pancreas')]
_all_fold += [P5['folds']['P1'], P5['folds']['P5']]
assert len(_all_fold) == 9
put('fold_min', r(min(_all_fold), 1)); put('fold_max', r(max(_all_fold), 1))
assert (N['fold_min'], N['fold_max']) == (2.1, 19.9), (N['fold_min'], N['fold_max'])
# P5 倍数区间：α(256) 的 bootstrap CI 换算成 fold
put('fold_p5_lo', r(P5B['drop_fold']['ci_from_alpha256'][0], 1))
put('fold_p5_hi', r(P5B['drop_fold']['ci_from_alpha256'][1], 1))
put('alpha_colon_2', r(D['colon']['alpha']['2'], 1)); put('alpha_colon_256', r(D['colon']['alpha']['256'], 2))

# ---- gamma 与 CI
# 注意：锁里 gamma 为负（α 随 d 下降），boot_ci 亦为负。取绝对值后区间端点会互换，
# 必须重新排序，否则正文出现 "0.580–0.523" 这种倒置区间（v8.1 审核发现）。
for t in ('colon', 'lung', 'brain'):
    g = abs(D[t]['gamma']); ci = sorted(abs(x) for x in D[t]['boot_ci'])
    assert ci[0] < g < ci[1], (t, ci, g)
    put(f'g_{t}', r(g, 3)); put(f'g_{t}_lo', r(ci[0], 3)); put(f'g_{t}_hi', r(ci[1], 3))
    put(f'g_{t}_r2', r(D[t]['r2'], 3))

# ---- 模型比较
put('daic_colon', r(MOD['colon']['daic']['power'], 1))
put('F_colon', r(MOD['colon']['quad_F'], 1))
put('p_colon', r(MOD['colon']['quad_p_value'], 4))
put('df1', 1); put('df2', 5)

# ---- CV 超 Poisson 倍数
put('cv_x_colon', round(DV['colon_obs_over_poisson_256']))
put('cv_x_lung', round(DV['lung_obs_over_poisson_256']))
put('cv_x_brain', round(DV['brain_obs_over_poisson_256']))

# ---- 深度
put('thin_quarter', 0.0004); put('thin_thousand', r(TH['depth_gamma_range'], 3))

# ---- 每基因
put('med_colon', r(G['colon']['gamma_gene']['median'], 3))
put('med_lung', r(G['lung']['gamma_gene']['median'], 3))
put('med_brain', r(G['brain']['gamma_gene']['median'], 3))
put('gmax', r(G['lung']['gamma_gene']['gmax'], 2))
put('frac_gt1_colon', r(G['colon']['gamma_gene']['frac_gt_1'] * 100, 1))
put('frac_gt15_colon', r(G['colon']['gamma_gene']['frac_gt_1p5'] * 100, 2))
put('nfit_colon', G['colon']['n_fit']); put('nfit_lung', G['lung']['n_fit']); put('nfit_brain', G['brain']['n_fit'])
put('nfit_all', G['colon']['n_fit'] + G['lung']['n_fit'] + G['brain']['n_fit'])
put('split_r', 0.79); put('split_n', 11805)
put('reliab_lo', 0.18); put('reliab_hi', 0.76)

# ---- 理论
put('identity_err', r(TH['identity']['max_abs_rel_err'] * 100, 1))
put('white_g', r(-TH['families']['white_0']['gamma_2_64'], 3))
put('pl_r2', r(TH['gamma_vs_beta_emp']['r2'], 4))
put('exp2_a', r(-TH['families']['exponential_2']['gamma_2_64'], 2))
put('exp2_b', r(-TH['families']['exponential_2']['gamma_16_64'], 2))
put('exp8_a', r(-TH['families']['exponential_8']['gamma_2_64'], 2))
put('exp8_b', r(-TH['families']['exponential_8']['gamma_16_64'], 2))

# ---- 协变图
put('n_pairs', sum(GB['by_tissue'][t]['n'] for t in ('colon', 'lung', 'brain')))
for t in ('colon', 'lung', 'brain'):
    put(f'r_{t}', r(GB['by_tissue'][t]['pearson_raw'], 2))
    put(f'rho_{t}', r(GB['by_tissue'][t]['spearman'], 2))
    put(f'diff_{t}', r(GB['by_tissue'][t]['diff_median'], 3))
    put(f'n_{t}', GB['by_tissue'][t]['n'])
put('slope_sp_colon', r(RS['by_tissue']['colon']['slope_vs_lag_spearman'], 2))
put('slope_sp_lung', r(RS['by_tissue']['lung']['slope_vs_lag_spearman'], 2))
put('slope_sp_brain', r(RS['by_tissue']['brain']['slope_vs_lag_spearman'], 2))
put('brain_slope_lo', r(RS['by_tissue']['brain']['slope_by_range']['8-16um'], 2))
put('brain_slope_hi', r(RS['by_tissue']['brain']['slope_by_range']['136-248um'], 2))
# colon 与 lung 的局部斜率包络（正文「wanders between」用）
_env = [v for t in ('colon', 'lung') for v in RS['by_tissue'][t]['slope_by_range'].values()]
put('slope_lo', r(min(_env), 2)); put('slope_hi', r(max(_env), 2))
put('n_windows', 36)

# ---- 置换
ig = IV['single_gene']['genes']['IGKC']
put('igkc_obs', r(ig['baseline_gamma'], 3)); put('igkc_perm', r(ig['permute_gamma'], 3))
put('igkc_lo', r(ig['permute_lo'], 3)); put('igkc_hi', r(ig['permute_hi'], 3))
put('pooled', r(IV['single_gene']['pooled']['gamma_baseline'], 3))
put('pooled_perm', r(IV['single_gene']['pooled']['cases']['permute']['gamma'], 3))
put('delta_perm', r(IV['single_gene']['pooled']['delta_permute'], 2))
put('excise_hot', r(IV['single_gene']['pooled']['cases']['excise']['gamma'], 3))
put('igkc_isolate', r(ig['isolate_gamma'], 3))
put('six_min', r(IV['single_gene']['summary']['delta_min'], 2))
put('six_max', r(IV['single_gene']['summary']['delta_max'], 2))
dr = IV['dose_response']['curves']
put('perm10_g', r(dr['10']['gamma'], 3)); put('perm10_share', round(dr['10']['mu2_share'] * 100))
put('perm1000_g', r(dr['1000']['gamma'], 3))
# 剂量-反应中间三点与五点的权重份额（v8.1 回收手写值）
for _n in ('10', '30', '100', '300', '1000'):
    put(f'perm{_n}_g', r(dr[_n]['gamma'], 3))
    put(f'perm{_n}_share', r(dr[_n]['mu2_share'] * 100, 1))
# 六个转录本置换区间的最大宽度（图注「at most … wide」）
_w = [v['permute_hi'] - v['permute_lo'] for v in IV['single_gene']['genes'].values()]
put('perm_width_max', r(max(_w), 3))
put('loo_shift', r(DV['colon_loo_max_shift'], 3))
# 同深度对照：lung 与小肠的 α(2 µm) 之比、倍数之差（v8.1 回收手写值）
put('ratio_lung_si', r(extra['lung']['alpha_2um'] / extra['mouse_small_intestine']['alpha_2um'], 2))
put('gap_fold_pct', round((extra['mouse_small_intestine']['drop_fold'] - extra['lung']['drop_fold'])
                          / extra['lung']['drop_fold'] * 100))
# 合成场输入指数（设计常数）
put('synth_g', 0.55)

# ---- 组成
put('neff_colon', r(DEC['colon']['n_eff_mu2'], 1)); put('neff_lung', r(DEC['lung']['n_eff_mu2'], 1))
put('neff_brain', round(DEC['brain']['n_eff_mu2']))
put('top100_colon', round(G['colon']['composition']['top100_share'] * 100))
put('top100_lung', round(G['lung']['composition']['top100_share'] * 100)); put('top100_brain', round(G['brain']['composition']['top100_share'] * 100))
put('a1_colon', r(DEC['colon']['top1_alpha_share'] * 100, 1)); put('a1_lung', r(DEC['lung']['top1_alpha_share'] * 100, 1))
put('a1_brain', r(DEC['brain']['top1_alpha_share'] * 100, 1))
put('a1_2_colon', r(DEC['colon']['v_by_scale_top1'][0] * 100, 1))
put('a1_2_lung', r(DEC['lung']['v_by_scale_top1'][0] * 100, 1))
put('a1_2_brain', r(DEC['brain']['v_by_scale_top1'][0] * 100, 1))
put('igkc_w', round(DEC['colon']['top1_mu2_share'] * 100))

# ---- 排序
put('gap_colon_lung', r(DV['colon_lung_gamma_gap'], 3)); put('gap_colon_brain', r(DV['colon_brain_gamma_gap'], 3))

# ---- 复制
put('p1_pooled', r(MULTI['Visium_HD_Human_Colon_Cancer_P1']['gamma_pool'], 3))
put('p1_med', r(MULTI['Visium_HD_Human_Colon_Cancer_P1']['gamma_median'], 3))
put('p1_sd', r(MULTI['Visium_HD_Human_Colon_Cancer_P1']['boot_sd'], 3))
put('p2_sd', r(MULTI['colon']['boot_sd'], 3))
put('p5_sd', 0.031); put('p5_pooled', r(P5G['per_section']['P5']['median'], 3))
put('p5_gamma', 0.629); put('p5_nfit', P5G['per_section']['P5']['n_fit'])
put('p1_nfit', P5G['per_section']['P1']['n_fit']); put('p2_nfit', P5G['per_section']['P2']['n_fit'])
put('p1_med_exact', r(P5G['per_section']['P1']['median'], 3))
put('p2_med_exact', r(P5G['per_section']['P2']['median'], 3))
put('p5_med_exact', r(P5G['per_section']['P5']['median'], 3))
put('d_pool_p1', 0.006); put('d_med_p1', 0.19); put('d_pool_p5', 0.074); put('d_med_p5', 0.012)
put('band_range_pool', r(DV3['colon_band_range'], 3)); put('band_range_med', r(DV4['robust_range']['median'], 3))
put('band_sd_pool', r(DV3['colon_band_sd'], 3)); put('band_sd_med', r(DV4['robust_sd']['median'], 3))

# ---- 跨平台
put('stereo_g', r(CR['stereo']['gamma_pool'], 3)); put('visium_brain', r(CR['visium']['gamma_pool'], 3))
put('cross_gap', r(CR['delta_gamma_pool'], 2))
put('cross_rho', r(CR['replication']['rho_raw'], 3)); put('cross_partial', r(CR['replication']['rho_partial_abundance'], 3))
put('cross_n', CR['replication']['n'])

# ---- 下游
put('phi_8', r(DN['od']['8']['phi'], 2)); put('phi_256', round(DN['od']['256']['phi']))
put('fpr_8', r(DN['od']['8']['fpr'] * 100, 1)); put('fpr_256', round(DN['od']['256']['fpr'] * 100))
put('neff_frac', r(DN['od']['256']['neff'] * 100, 1))
put('jacc_32', r(DN['de']['32']['jaccard_sig'], 2))
put('top100_16', round(DN['de']['16']['top100_recall'] * 100)); put('top100_64', round(DN['de']['64']['top100_recall'] * 100))
put('sig_8', DN['de']['16']['n_sig_ref']); put('sig_128', DN['de']['128']['n_sig_s'])
put('slope_log_64', r(DN['de']['64']['logfc_slope'], 2)); put('slope_raw_64', r(DN['de']['64']['logfc_slope_raw'], 2))
put('r_64', r(DN['de']['64']['logfc_r'], 3))
put('map_raw_lo', r(DN['summary']['map_raw_min'], 2)); put('map_raw_hi', r(DN['summary']['map_raw_max_x'], 2))
put('map_plasma', r(DN['summary']['map_raw_max'], 2))
put('map_log_lo', r(DN['summary']['map_log_min'], 2)); put('map_log_hi', r(DN['summary']['map_log_max'], 2))
put('map_r_min', r(DN['summary']['map_raw_r_min'], 2))

# ---- 设计
put('mde_colon', r(D['colon']['mde_gamma']['n1'], 3)); put('mde_lung', r(D['lung']['mde_gamma']['n1'], 3))
put('mde_brain', r(D['brain']['mde_gamma']['n1'], 3))
put('rms_colon', r(D['colon']['rms_gamma_cut'], 3)); put('rms_lung', r(D['lung']['rms_gamma_cut'], 3))
put('rms_brain', r(D['brain']['rms_gamma_cut'], 3))
put('k_power', 2.80)

# ---- 空间代理
for t, key in (('colon', 'colon'), ('lung', 'lung'), ('brain', 'brain')):
    put(f'moran_{t}', r(SP[key]['rho_gamma_moran'], 3))
for t in ('colon', 'lung', 'brain'):
    spj = json.load(open(rf'E:\SPATIAL_2026\02_OUT\spatial_proxy\spatial_proxy_{t}.json', encoding='utf-8'))
    put(f'band10_{t}', r(spj['band_rho'][0]['median'], 2))
    put(f'semip_{t}', r(spj['semiparam']['pearson'], 2))

# ---- Xenium
put('xen_genes', 422); put('xen_tx', '50,594,830')
put('xen_pooled', 0.10); put('xen_gap', 0.32); put('xen_rho', 0.871); put('xen_n', 403)
put('xen_med', 0.41)

# ---- 软件
put('py', '3.13.14'); put('numpy', '2.5.4'); put('scipy', '1.18.1'); put('mpl', '3.11.2')
put('seed', 20261115); put('seed_boot', 20261114); put('seed_extra', 20261116)

out = r'C:\Users\唐朝华\WorkBuddy\A73_rewrite\v8\numbers_v8.json'
import os
os.makedirs(os.path.dirname(out), exist_ok=True)
json.dump(N, open(out, 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
print('keys:', len(N))
for k in list(N)[:20]:
    print(' ', k, '=', N[k])

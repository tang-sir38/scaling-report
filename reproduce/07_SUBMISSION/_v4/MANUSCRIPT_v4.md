# TITLE
Binning is set by spatial covariance and by which transcripts dominate

# AUTHORS
Chaohua Tang¹,*,**, Fangfang Wu¹,*, Yingying Zheng¹,*, Lingling Zhan¹, Hong Tan¹, Lijing Lin¹, Meimei Wen¹

* These authors contributed equally to this work.
** Corresponding author.

# AFFILIATION
¹ Department of Gynecology, The 95th Hospital of Putian, No. 485 Dongyan Road, Licheng District, Putian City, Fujian Province 351100, China

# CORRESPONDENCE
Chaohua Tang, Department of Gynecology, The 95th Hospital of Putian, No. 485 Dongyan Road, Licheng District, Putian City, Fujian Province 351100, China. Telephone: +86 13860904395. Email: 13860904395@163.com. ORCID: 0009-0004-7080-8993.

# RUNNING_TITLE
Binning noise is a covariance and composition statistic

# KEYWORDS
spatial transcriptomics; Visium HD; overdispersion; change of support; scaling; experimental design

# ABSTRACT
Binning buys counts with resolution, on the argument that overdispersion falls as the inverse bin area, α ∝ d^−2. Re-aggregating native 2 µm Visium HD bins to 256 µm in human colon, human lung and mouse brain, we measure α ∝ d^−{{c.gamma:.3f}}, d^−{{l.gamma:.3f}} and d^−{{b.gamma:.3f}}: {{c.fall:.0f}}-fold in colon and lung and {{b.fall:.0f}}-fold in brain, instead of 16,384-fold. We show that α(d) equals the block-averaged spatial correlation function to within {{th.idmaxpct:.1f}}% in synthetic fields, so the exponent is a covariance scaling exponent. Per gene it agrees with an independently measured covariogram across {{gb.total:,}} genes (Pearson {{gb.c.r:.2f}}, {{gb.l.r:.2f}}, {{gb.b.r:.2f}}), and permuting transcript positions — which leaves every abundance unchanged — raises the pooled exponent from {{iv.pool.base:.2f}} to {{iv.pool.perm:.2f}}. Because the pooled exponent is a weighted average in which a handful of transcripts carry almost all the weight, it is a composition statistic of the sample, not a constant of the tissue.

# BODY

## Introduction

Spatial transcriptomics measures expression without dissociating the tissue, and pays for that with sparsity [[stahl2016,stickels2020,chen2022,eng2019,chen2015,he2022]]. A Visium HD section tiles its capture area with a contiguous 2 µm grid, so one 6.5 mm section is {{c.nbins:,}} bins holding {{c.umi:.0f}} transcripts each [[janesick2023]]. Most bins are empty for most genes, and every downstream step — clustering, differential expression, cell-type mapping — inherits that noise [[moses2022,rao2021,luecken2019,heumos2023,li2022]].

The field's standard response is to merge neighbouring bins [[sun2025,polanski2024]], or to keep the native grid and deconvolve mixed pixels instead [[zhao2021,kleshchevnikov2022,cable2021,biancalani2021]]. The merged unit is then treated as a more precise measurement of the same quantity, and the justification is arithmetic. For a Poisson count with mean μ the coefficient of variation of the per-bin count is CV² = 1/μ; if μ grows with bin area A then CV² falls as 1/A, so the overdispersion α = (Var − E)/E² falls as d⁻², in the same mean–variance framework that governs count modelling elsewhere [[grun2014,hafemeister2019,ahlmann2023,sarkar2021,brennecke2013,eling2019,svensson2020]]. Coarsening from 2 to 256 µm multiplies the area by 16,384, so α should fall 16,384-fold and the loss of resolution should cost almost nothing. That expectation is why bin sizes are usually chosen for convenience rather than measured.

Whether α falls as d⁻² has not been measured over the range in which it is used. The quantity is easy to compute and the data to compute it are already on disk: any section can be re-aggregated exactly at every dyadic scale, without interpolation and without dropping reads. What has been missing is a reason to expect a particular value, and therefore a way to recognise the answer.

Two adjacent literatures supply the reason. Ecology has measured how the variance of a field depends on the grain at which it is sampled since the 1980s; the variance–grain relation and the autocorrelation of the same field are two views of one quantity [[wiens1989,legendre1993]]. Geostatistics formalised the same relation as change of support: the variance of a block is the integral of the covariance over the block, so regularising to a coarser support removes within-block variance at a rate fixed by the variogram [[matheron1963,journel1978]]. In geography the resulting dependence of a summary on the reporting unit is the modifiable areal unit problem [[jelinski1996]]. In none of these is the exponent a free parameter, and in none is it 2 unless the field has no spatial structure. Binned transcriptomics sits inside this framework without having noticed.

Here we measure the exponent in three Visium HD sections — human colon cancer, human lung cancer and mouse brain — and give it a mechanism. We show that α(d) equals the block-averaged spatial correlation function up to a constant, so that the exponent is a covariance scaling exponent: 2 for an uncorrelated field, equal to the covariogram exponent for a power-law field, and undefined for a field with a finite correlation length. We verify the identity in synthetic fields, test the per-gene prediction against covariograms measured along an independent path, and confirm it causally by permuting transcript positions while leaving every abundance unchanged. We then show what the exponent's dependence on a handful of transcripts implies, namely that the pooled value is a composition statistic of the sample: reproducible between sections, but not a property of the tissue.

## Results

### Noise falls fourteen-fold, not sixteen-thousand-fold

We re-aggregated the native 2 µm bins of each section by exact integer summation to 4, 8, 16, 32, 64, 128 and 256 µm on the same grid (Fig. 1a). No read is dropped and no value is interpolated: a 256 µm bin is the sum of the 16,384 native bins inside it.

Pooled overdispersion falls {{c.fall:.1f}}-fold in colon, from {{c.a2:.1f}} at 2 µm to {{c.a256:.2f}} at 256 µm (Fig. 1b). The counting argument predicts 16,384. A single power law fits the eight points with r² = {{c.r2:.3f}} and exponent γ = {{c.gamma:.3f}}. The other two sections behave the same way: lung falls {{l.fall:.1f}}-fold (γ = {{l.gamma:.3f}}, r² = {{l.r2:.3f}}) and mouse brain {{b.fall:.1f}}-fold (γ = {{b.gamma:.3f}}, r² = {{b.r2:.3f}}).

The gap is not a small correction to an otherwise correct argument. At 256 µm the colon value is still {{c.pois:.0f}} times the Poisson value for the same counts, and the lung and brain sections are {{l.pois:.0f}} and {{b.pois:.0f}} times. Binning removes the Poisson component and nothing else, so what remains at every scale is the variance of the underlying intensity field, which decays slowly. Thinning the colon section to a quarter of its counts moves γ by {{des.depth:.4f}}, so this is not a depth effect.

### The exponent is the decay rate of the correlation function

The arithmetic that motivates binning can be carried one step further. Let the intensity field have mean λ̄ and covariance C(h) = σ²ρ(h). For a block A, the count N_A satisfies E[N_A] = λ̄|A| and Var(N_A) = λ̄|A| + Var(Λ_A) with Λ_A = ∫_A λ. Hence

  α(A) ≡ (Var N_A − E N_A)/(E N_A)² = Var(Λ_A)/(λ̄|A|)² = (σ²/λ̄²) · ρ̄(A),

where ρ̄(A) = |A|⁻²∫∫_A ρ(x − y) dx dy is the correlation function averaged over pairs inside the block (Methods). The Poisson term cancels identically, so α is the block-averaged correlation function multiplied by a constant that does not depend on support.

Three consequences follow, and they separate the cases the field has conflated. For an uncorrelated field ρ̄ = 1/k² at block size k, so α ∝ d^−2: the counting limit is one special case, not the rule. For a power-law covariance ρ(h) ∝ h^−β in two dimensions ρ̄(A) ∝ A^−β/2 = d^−β, so γ = β — the noise–resolution exponent equals the covariance scaling exponent, with no free parameter. For a covariance with a finite correlation length ξ, ρ̄ has no power-law range, so no single γ exists and any fitted value moves with the fit window.

We verified all three in synthetic log-normal fields on a periodic 1,024² grid (Methods). The identity holds to a maximum relative error of {{th.idmaxpct:.1f}}% (median {{th.idmedpct:.1f}}%) across coefficients of variation spanning two orders of magnitude. An uncorrelated field returns γ = {{th.wn:.3f}}. A power-law field returns γ equal to the covariogram exponent measured on the same realisations (r² = {{th.r2:.4f}}), although both sit above the nominal continuum β because a discrete grid truncates the spectrum. Exponential fields with ξ = 2 and 8 grid units return γ = {{th.ex2a:.2f}} and {{th.ex8a:.2f}} over the full range but {{th.ex2b:.2f}} and {{th.ex8b:.2f}} over the upper half: the fitted exponent moves with the window, exactly as the identity requires. Thinning a synthetic field a thousand-fold changes γ by {{th.depth:.3f}}.

The identity also makes a prediction that does not use the α curve at all. If γ is a covariance exponent, gene g's exponent should equal the exponent of its covariogram measured directly on the tissue — a different estimator applied to different objects, lagged products of the 8 µm field rather than block variances. We measured both for every gene passing the filters, {{gb.total:,}} genes across the three sections. The two agree in rank (Pearson {{gb.c.r:.2f}}, {{gb.l.r:.2f}} and {{gb.b.r:.2f}}; Spearman {{gb.c.rho:.2f}}, {{gb.l.rho:.2f}} and {{gb.b.rho:.2f}}) and in level: the median difference γ − β is {{gb.c.diff:+.3f}} in colon, {{gb.l.diff:+.3f}} in lung and {{gb.b.diff:+.3f}} in brain (Fig. 4a).

The residual disagreement is informative rather than embarrassing, because γ and β are different functionals. β is the local slope of ρ(h); γ is the slope of its block average, which weights every separation inside the block roughly equally and is therefore smoother. The two coincide if and only if ρ is a power law. In colon and lung they are: the local slope of the median ρ(h) wanders between {{rho.c.lo:.1f}} and {{rho.c.hi:.1f}} with no trend across the whole 8–256 µm range (Spearman against lag {{rho.c.trend:+.2f}}). In mouse brain it climbs monotonically, from {{rho.b.lo:.2f}} at 8–16 µm to {{rho.b.hi:.2f}} at 136–248 µm (Spearman {{rho.b.trend:+.2f}}), so there is no single β (Fig. 4b). Accordingly there is no single number for γ to equal, and the value it takes, {{gb.b.g:.3f}}, sits at the low end of the range spanned by the local slopes — as a smoothed functional must.

### What the identity does and does not imply

A reader's first thought is that the pooled exponent is a weighted average and therefore contains no information. It is a weighted average, and the arithmetic has to be separated from the empirical content. Because α_pool(d) = Σ_g μ_g²α_g/Σ_g μ_g², three statements are true by construction: the pooled exponent lies between the per-gene exponents; it can be rebuilt exactly from them, which we verify to machine precision; and genes with vanishing μ² cannot influence it. None of that is evidence about tissue.

What the arithmetic does not supply is the following. That the exponent is between 0.2 and 1.9 rather than 2. That the three sections order differently under the pooled and the per-gene summary. That one transcript supplies {{g.c.share1:.0%}} of the pooled α in colon and {{g.b.share1:.0%}} in brain. That shuffling positions moves the pooled exponent by {{iv.pool.delta:.2f}}. Each of these is a fact about a particular tissue that could have come out otherwise, and one of them — the brain section's covariance — did. The decomposition tells us where to look; it does not tell us what we will find.

### Moving transcripts moves the exponent

If γ is a covariogram exponent then it is a statement about where transcripts are, and that can be tested by moving them. For a chosen gene we kept its multiset of counts exactly — same total, same values, same occupancy — and reassigned them to uniformly random positions inside the tissue domain, recomputing all eight scales (Methods). Nothing about abundance, depth or chemistry changes; only the arrangement does.

For IGKC, the transcript that carries {{g.c.top1share:.0%}} of the colon section's μ² weight, the per-gene exponent rises from {{iv.igkc.base:.3f}} to {{iv.igkc.perm:.3f}} (95% range over 20 independent permutations {{iv.igkc.lo:.3f}}–{{iv.igkc.hi:.3f}}), close to the uncorrelated limit (Fig. 5a). The pooled exponent rises from {{iv.pool.base:.3f}} to {{iv.pool.perm:.3f}}. Removing IGKC's high-density regions instead of permuting them lowers the pooled exponent to {{iv.pool.excise:.3f}}, and taking IGKC as the only gene gives {{iv.pool.isolate:.3f}}: the arrangement, not the abundance, sets the exponent. Five further transcripts spanning three orders of magnitude in abundance behave identically — IGHG1, JCHAIN, MZB1, EPCAM and COL1A1 rise by {{iv.dmin:.3f}} to {{iv.dmax:.3f}}, every one of the six in the same direction (Fig. 5a; Table S3).

The effect is graded. Permuting the top N genes by μ² weight, with N = 10, 30, 100, 300 and 1,000, raises the pooled exponent monotonically from {{iv.dn0:.3f}} to {{iv.dn1:.3f}}, {{iv.dn2:.3f}}, {{iv.dn3:.3f}}, {{iv.dn4:.3f}} and {{iv.dn5:.3f}} (Fig. 5b; Table S4). Ten transcripts, {{iv.sh1:.0%}} of the weight, take the exponent from {{iv.dn0:.2f}} to {{iv.dn1:.2f}}. The pooled exponent is not a property of the assay, the platform or the depth; it is a property of how a small number of transcripts are arranged.

### A composition statistic

The effective number of genes behind the pooled exponent, 1/Σw² with w the μ² weights, is {{g.c.neff:.1f}} in colon, {{g.l.neff:.1f}} in lung and {{g.b.neff:.0f}} in brain; the top 100 genes carry {{g.c.top100:.0%}}, {{g.l.top100:.0%}} and {{g.b.top100:.0%}} of the weight. Weighting by the share of the exponent itself concentrates it further: the heaviest transcript supplies {{g.c.share1:.0%}} of the pooled α in colon and {{g.l.share1:.0%}} in lung, but {{g.b.share1:.0%}} in brain (Fig. 6).

Because the weights differ between samples, one statistic ranks the three sections in opposite orders. Pooled, colon ({{g.c.pooled:.3f}}) exceeds lung ({{g.l.pooled:.3f}}), which exceeds brain ({{g.b.pooled:.3f}}). Per gene, the medians run colon {{g.c.med:.3f}}, lung {{g.l.med:.3f}}, brain {{g.b.med:.3f}}: brain moves from last to first. Both numbers are correct and they answer different questions — how much of this sample's noise coarsening removes, and how a typical transcript's noise is removed — and a literature that quotes one while arguing about the other will not converge.

The pooled value is also reproducible where the per-gene summary is not. A second colon section, captured separately and carried through the same pipeline, gives {{c2.gamma:.3f}} ± {{c2.sd:.3f}} against {{c.gamma:.3f}}; its per-gene median is {{c2.med:.3f}} against {{g.c.med:.3f}}, a difference of {{c2.meddiff:.2f}} that no abundance correction removes. Across a change of chemistry the same tissue gives {{st.gamma:.3f}} on Stereo-seq against {{b.gamma:.3f}} on Visium HD, {{st.delta:.2f}} apart, while the per-gene values still replicate in rank across platforms (ρ = {{st.rho:.3f}} over {{st.n:,}} shared genes, {{st.rhopart:.3f}} after controlling for abundance). Spatial-statistics tools that operate on a binned grid inherit the same weights [[svensson2018,edsgard2018,sun2020,weber2023,dries2021,palla2022]]. Within one section, cutting the colon into four column bands moves the pooled exponent over {{rob.pooled.range:.3f}} and the per-gene median over {{rob.median.range:.3f}}; the between-band standard deviations are {{rob.pooled.sd:.3f}} and {{rob.median.sd:.3f}}.

### What the exponent costs the analyses that follow

The exponent can be used before any analysis is run, because it prices the downstream steps directly. A test that treats the counts in a bin as Poisson assumes Var = μ where the variance is μ + μ²α, so it understates it by a factor φ = 1 + μα. In the colon section φ is {{ds.phi_lo:.2f}} at 8 µm and {{ds.phi_hi:.0f}} at 256 µm. The growth is not caused by noisier data: μ rises with the square of the bin width while α falls only as d^−γ with γ below 2, so the product rises. For a nominal 5% test the implied false-positive rate goes from {{ds.fpr_lo:.1%}} to {{ds.fpr_hi:.0%}}, at which point the effective sample size is {{ds.neff:.1%}} of the bin count (Fig. S2a; Table S7). A negative-binomial or rank-based test does not pay this price in full, so this is an upper bound — but it is the bound that the counting argument implicitly sets to one.

We then asked the same question of a real comparison. Separating the section into an epithelial end and a stromal end by marker score at 16 µm and repeating a Wilcoxon test at every coarser scale, the direction of every effect survives: log fold changes correlate with the 8 µm values at r = {{ds.r64:.3f}} even at 64 µm. Rank does not. The significant set is almost unchanged at 32 µm (Jaccard {{ds.jacc32:.2f}}), and the top 100 genes are {{ds.top16:.0%}} recoverable at 16 µm but only {{ds.top64:.0%}} at 64 µm, while the number passing FDR < 0.05 with |log2FC| > 0.5 falls from {{ds.sigref:,.0f}} to {{ds.sig128:,.0f}} at 128 µm.

Magnitude is worse than unstable: it reverses with the normalisation. Regressing coarse-scale log fold changes on the 8 µm values gives a slope of {{ds.sl64:.2f}} when the mean is taken after log1p, the scanpy default, and {{ds.sr64:.2f}} when it is taken on the raw ratio, from the same counts. The two differ in sign, not in size. A pipeline that reports attenuation and one that reports inflation can both be right about their own numbers and both wrong about the tissue.

Cell-type mapping shows the same asymmetry. Module scores computed at 8 µm and aggregated to a coarse bin are reproduced by that bin's own counts with slopes of {{ds.rawmin:.2f}}–{{ds.rawmaxx:.2f}} on the raw ratio, reaching {{ds.rawmax:.2f}} for the plasma-cell set that contains the transcript carrying most of the weight, and the spatial pattern is preserved throughout (r ≥ {{ds.rawrmin:.2f}}). After log1p the same scores span {{ds.logmin:.2f}}–{{ds.logmax:.2f}} (Fig. S2b; Table S8). Sparsely expressed marker sets move furthest. The pattern survives; the number attached to it does not, and a module score compared across scales at absolute value is reporting the normalisation rather than the tissue.

### Consequences for design

The design consequences follow from an inequality that does not depend on the mechanism. Every transcript measured here has γ ≤ {{des.gmax:.2f}}; {{g.c.gt1:.1%}} of colon genes exceed 1 and {{g.c.gt15:.2%}} exceed 1.5. Whatever the arrangement, α cannot fall faster than d^−1.91 over the range studied, so the variance reduction bought by coarsening is bounded well below the counting limit and the resolution spent to obtain it is not recoverable from the counts.

Four things follow. First, a sample's noise is not summarised by its exponent: report the per-gene distribution, the effective number of genes and the between-region spread, because a single number can be moved by one transcript. Second, comparisons between conditions need sections rather than bins, for the same reason that pseudoreplication inflates single-cell studies [[squair2021,zimmerman2021,murphy2022]]: with one section per group the minimum detectable difference in γ at 80% power is {{des.mde.c:.3f}}, {{des.mde.l:.3f}} and {{des.mde.b:.3f}}, larger than the colon–lung gap of {{des.gap_cl:.3f}} and comparable to the colon–brain gap of {{des.gap_cb:.3f}}, and of the same order as the sample sizes power analyses for other spatial assays already demand [[svensson2017,ryaboshapkina2023,vallejos2014]]. Third, the available levers are unequal in a way that can be read off the decomposition: halving depth moves γ by {{des.depth:.4f}}, dropping the heaviest transcript by {{des.drop:.3f}}, and changing region by up to {{rob.pooled.range:.3f}}. Fourth, when a pipeline changes which transcripts dominate — panel selection, ambient removal, deconvolution [[eraslan2019,wu2021]] — it will change the pooled exponent without changing the tissue, and that change should be reported as part of the result rather than interpreted.

## Discussion

Binned spatial transcriptomics does not behave the way its defaults assume [[lahnemann2020]], and the reason is not that the counts are noisier than expected. It is that the quantity being measured is not a counting statistic at all. Once the Poisson term is removed, what is left is the variance of the intensity field, and the rate at which that variance falls with support is the rate at which the field's correlation function decays. In geostatistics this is the regularisation relation [[matheron1963,journel1978]]; in ecology the variance–grain relation [[wiens1989,legendre1993]]; in geography the modifiable areal unit problem [[jelinski1996]]. The exponent is therefore not a property of the assay — and not a constant of a tissue either, because a weighted average over transcripts is not a property of anything except the weights.

That reframing changes what a binned analysis should report. An exponent quoted for a sample is a composition statistic: reproducible for that sample, at {{c.gamma:.3f}} and {{c2.gamma:.3f}} for two sections of one tumour, and not transportable to another, at {{b.gamma:.3f}} and {{st.gamma:.3f}} for one tissue on two chemistries. What does transport is the per-gene exponent, which replicates in rank across platforms (ρ = {{st.rho:.3f}}) and whose level is predicted by the covariogram. The four sections were not chosen to span tissues but to span covariance structures — two power-law and one not — because that is what the exponent depends on.

Three limitations bound these claims. The covariogram is measured on the 8 µm grid, so correlations below 8 µm are invisible and only the decay rate of ρ̄ can be checked, not its absolute height; the two analyses that fail for exactly this reason are reported in the Supplementary Note. The causal experiment permutes positions within the observed domain, which destroys spatial structure without reproducing a biologically realistic alternative, so it establishes that arrangement matters rather than what a particular arrangement would do. And the sections are whole-transcriptome on two platforms but a single capture-chemistry family, so a targeted panel may concentrate the weights differently, although it will obey the same identity.

One prediction is worth stating because it is falsifiable. Any step that changes which transcripts dominate the μ² weight — a panel, a filter, a deconvolution that redistributes counts — should move the pooled exponent by an amount computable from the weight it removes. We expect that prediction to hold, and we expect the resulting exponent to be reported alongside the step that produced it, because otherwise the number will be read as a property of the tissue.

The downstream prices are computable in advance too, and that changes what a methods section owes its reader. A study that reports an exponent can also report what its Poisson assumption costs at the scale it used, what fraction of a top-100 list survives a change of scale, and which normalisation its effect sizes live on. None of these requires new data: all three follow from the identity and from counts already in hand. A reader who is not told them will read a number computed at one support as a property of the tissue, which is the error this paper set out to name.

# METHODS

## Datasets

Four Visium HD sections and one Stereo-seq section were analysed. Visium HD human colon cancer (Visium_HD_Human_Colon_Cancer, 6.5 mm capture area), a second section of the same tumour block (Visium_HD_Human_Colon_Cancer_P1), Visium HD human lung cancer and Visium HD mouse brain were obtained as Space Ranger 2 µm output and used without further alignment. Stereo-seq mouse whole brain (C04042E2) was used as a cross-chemistry comparison. Mitochondrial genes were removed before analysis (18,074 genes retained in the human sections, 19,046 in mouse brain). No sample size calculation was performed; the number of sections was fixed by data availability.

## Exact re-aggregation and the analysis domain

Each section is represented as a sparse count matrix on its native 2 µm grid. A block of side d contains k = d/2 native cells per axis, and its count vector is the exact integer sum of the k² native count vectors inside it, obtained as a sparse matrix product with a 0/1 aggregation operator. No interpolation, no smoothing and no read filtering are applied, so the total number of transcripts is identical at every scale. The analysis domain at scale d is the set of blocks containing at least one in-tissue 2 µm bin; blocks lying entirely outside the tissue are excluded. Because the aggregation is exact, α(d) at every scale refers to the same transcripts, which is what makes the eight scales comparable.

## Pooled overdispersion

For each gene g and scale d, let μ_gd and v_gd be the mean and the unbiased variance of the per-block count over the analysis domain. The overdispersion is α_gd = (v_gd − μ_gd)/μ_gd². The pooled value is

  α_pool(d) = Σ_g μ_gd² α_gd / Σ_g μ_gd²,

which is the value obtained by computing (Var − E)/E² on the summed count field, because the Poisson term is additive in the same weights. The exponent γ is the negative slope of log α against log d over all eight scales (2–256 µm), fitted by ordinary least squares. This window was fixed before any comparison and is not varied between sections.

## Per-gene exponents, filtering and reliability

A gene enters the per-gene analysis if α_gd is finite and positive at all eight scales and if its mean count per 256 µm block is at least 0.5. This retains 13,115 of 18,074 genes in colon, 11,716 in lung and 11,571 in mouse brain. Reliability was assessed by splitting each section at the midpoint of the capture area into two halves, aggregating and fitting each half independently, and correlating the two per-gene estimates (Pearson r = 0.79 over 11,805 genes in colon). Reliability rises with abundance, from r = 0.18 in the lowest expression layer to r = 0.76 in the highest.

## Theory: the identity and synthetic fields

The identity α(A) = (σ²/λ̄²)ρ̄(A) follows from conditional independence of the counting step given the intensity field, which gives Var(N_A) = E[N_A] + Var(Λ_A); the first term is exactly E[N_A] and cancels in (Var − E)/E². It was verified on synthetic log-normal fields on a periodic 1,024² grid with eight realisations per configuration, a log-scale standard deviation of 0.8 and a mean count of 1.0 per cell. Three covariance families were used, defined by their spectra: a power law, S(k) ∝ k^(β−2); an exponential covariance with correlation length ξ, S(k) ∝ (1 + (2πξk)²)^(−3/2); and white noise. For each realisation the block-averaged correlation ρ̄ was computed from the exact autocovariance by summing over within-block lags with the correct pair multiplicities, and compared with the measured α at block sizes k = 1, 2, 4, 8, 16, 32 and 64. The reference exponent for the power-law family is the covariogram exponent measured on the same realisations, not the nominal spectral parameter, because discretisation and Nyquist truncation raise the measured value above the nominal one. Depth invariance was tested by binomial thinning of a realised field to mean counts of 30, 10, 3, 1, 0.3 and 0.1 per cell.

## Covariogram and the independent exponent

For each gene the covariogram was estimated on the 8 µm grid, the coarsest scale at which individual transcripts are not quantised. The masked autocovariance was computed by FFT with the mask autocorrelation in the denominator, so that only pairs of in-domain cells contribute, and the mean was subtracted before the transform to avoid the lag-zero bias. The two-dimensional covariance was reduced to an isotropic profile by averaging over rings of integer lag, normalised by the lag-zero value, and fitted by log ρ = c − β log h over lags 2–24 (16–192 µm). A gene was retained if this fit had r² > 0.5, which left 2,476 genes in colon, 1,726 in lung and 809 in mouse brain. The local slope of ρ was computed as the negative finite difference of log ρ between consecutive lags on the median profile.

## Permutation and dose–response experiments

For a chosen gene, the multiset of non-zero per-bin counts was retained and reassigned to uniformly random in-domain 2 µm positions without replacement; the aggregation and the α curve were then recomputed exactly as for the observed field. This leaves the total count, the count distribution and the occupancy unchanged, so any change in α is attributable to spatial arrangement alone. Twenty independent permutations were performed per gene, and the reported range is the 2.5th–97.5th percentile of the resulting per-gene exponent. For the dose–response experiment, genes were ranked by their μ² weight at 2 µm and the top N were permuted simultaneously and independently for N = 10, 30, 100, 300 and 1,000, covering 78.4%, 83.1%, 87.2%, 90.7% and 94.8% of the total weight; the pooled exponent was then rebuilt by replacing only the permuted genes' contributions in the weighted sum.

## Window sensitivity

To establish that no conclusion depends on the fit window, γ was refitted over six windows — lower bound fixed at 2 µm with upper bounds of 32, 64, 128 and 256 µm, then upper bound fixed at 256 µm with lower bounds of 8 and 16 µm — and β over six lag windows (2–8, 2–12, 2–16, 2–24, 2–32 and 3–16 grid units). The full 6 × 6 matrix of correlation between the two per-gene estimates is reported (Fig. 4c and Table S2); the window used in the main text is the one fixed a priori, and the matrix is reported as a sensitivity analysis rather than used to select a window.

## Decomposition, regional bands, replication and cross-platform comparison

The pooled exponent was decomposed by gene: the share of the pooled numerator contributed by gene g is v_g = w_g α_g/α_pool with w_g = μ_g²/Σμ², and the effective number of genes is 1/Σ_g w_g². Rebuilding α_pool from the per-gene exponents reproduces the measured value to machine precision at every scale, which was verified as a consistency check rather than treated as a result. Regional instability was assessed by cutting the colon section into four column bands on the 256 µm grid, so that no coarse block is divided, and recomputing five candidate summaries inside each band. Section-level uncertainty was obtained from a block bootstrap over 512 µm superblocks with 200 resamples. The two colon sections were processed by the identical pipeline. For the cross-platform comparison, per-gene exponents were matched by gene symbol between the Stereo-seq and Visium HD mouse brain sections and compared by Spearman correlation, raw and after partialling out mean abundance.

## Design calculations

The minimum detectable difference in γ was computed as 2.80 σ_region, where σ_region is the standard deviation of γ across the four column bands and 2.80 = z_0.975 + z_0.80, for a two-sided test at α = 0.05 and 80% power with one section per group. Depth sensitivity was assessed by binomial thinning of the colon section to 46%, 25% and 10% of its counts and refitting. The effect of dropping the heaviest transcript was taken from the leave-one-gene-out analysis. All calculations were performed in Python 3.13.14 with numpy 2.5.4, scipy 1.18.1 and matplotlib 3.11.2; code and derived tables are available at the URL in the Code availability statement.

## Downstream cost calculations

The cost of the Poisson assumption was computed from the per-gene overdispersion at each scale. For a gene with mean μ per bin, treating counts as Poisson sets the variance to μ where the tissue gives μ + μ²α, so the variance is understated by φ = 1 + μα; the implied false-positive rate of a two-sided test at a nominal 5% is 2[1 − Φ(1.96/√φ)], and 1/φ is the fraction of the bin count that behaves as independent replicates. Values are medians over the 2,000 most abundant genes, computed on all bins of the section without further filtering.

For the differential-expression comparison the two ends of the section were defined at 16 µm by the difference between the mean epithelial marker score (EPCAM, KRT18, KRT8, CEACAM5, TFF3, CDX2, KRT19) and the stromal score (COL1A1, COL1A2, COL3A1, DCN, LUM, PDGFRB, THY1) on log1p(CP10K) values, taking the top and bottom 30% of bins. Labels were carried to other scales through the integer coordinate mapping of the grid (a bin at s µm covers (s/16)² bins at 16 µm); a coarse bin was retained when at least 80% of the 16 µm bins it covers carried the same label, and bins with fewer than 20 counts were excluded at every scale. The 3,000 most abundant genes were tested with a two-sided Wilcoxon rank-sum test on log1p(CP10K) values, with Benjamini–Hochberg q values. Log fold changes are reported both as the ratio of log1p means and as the ratio of raw means, because across scales the two have opposite signs.

For the mapping proxy, five marker sets (epithelial, stromal, T cell, plasma cell and endothelial, listed in Table S8) were scored as the sum of their counts divided by the bin total, in CP10K, and the 8 µm scores were aggregated to each coarse bin as the mean over the 8 µm bins it contains. The slope is an ordinary least-squares regression of the coarse bin's own score on that aggregated reference; a slope of one means the score is comparable across supports. Both quantities were computed with and without log1p.

# DECLARATIONS

**Funding.** This work received no external funding.

**Competing interests.** The authors declare no competing interests.

**Ethics approval.** Not applicable. All data analysed here are previously published, de-identified, publicly available datasets; no new human or animal material was collected.

**Consent for publication.** Not applicable.

**Data availability.** The four Visium HD sections and the Stereo-seq section analysed here are public data, distributed by 10x Genomics and by STOmicsDB respectively and used unmodified; the dataset names and the STOmicsDB identifier are listed in Table S1. No new data were generated.

**Code availability.** All analysis code, the numbers lock file from which every number in this manuscript is injected, and the derived tables are available at [URL TO BE INSERTED BEFORE SUBMISSION].

**Author contributions.** C.T. conceived the study, performed the analysis and wrote the manuscript. F.W., Y.Z. and L.Z. contributed to data curation and interpretation. H.T., L.L. and M.W. contributed to interpretation and revision. All authors approved the final version.

**Acknowledgements.** We thank the authors of the original datasets for making their data public.

# REFERENCES

（此节由构建脚本按正文首现顺序自动生成，模板内无需手写）

# FIGURE_LEGENDS

**Fig. 1. Noise falls one order of magnitude, not four.** (a) Exact integer re-aggregation: native 2 µm bins are summed into 4, 8, 16, 32, 64, 128 and 256 µm blocks on the same grid, with no interpolation and no dropped reads. (b) Pooled overdispersion normalised to its value at 2 µm, against bin size, for colon (dark blue), lung (mid blue) and mouse brain (ochre). Grey dashed, α ∝ d^−2, the counting-statistics limit that merging is meant to approach; grey dotted, α ∝ d⁰. If the field's working assumption held, the data would lie on the steep grey line.

**Fig. 2. The exponent belongs to transcripts.** (a) Per-gene exponent against mean counts per 32 µm bin for the three sections; lines are medians within abundance bins. (b) Distribution of per-gene exponents; white dot, median; bar, interquartile range; dashed, pooled value. (c) Per-gene exponent in the two halves of the colon section; r, Pearson correlation. (d) Median exponent by abundance decade for the three sections (solid) and for a synthetic control with a known input exponent of 0.55 (dashed); the control is flat, the data are not.

**Fig. 3. The exponent is the decay rate of the correlation function.** (a) The identity tested on synthetic fields: measured α against the prediction CV²₀ρ̄, for seven power-law covariances (circles) and six coefficients of variation (squares); dashed, the identity line. (b) Recovered exponent against the covariogram exponent measured on the same realisations; dashed, y = x. (c) Recovered exponent for exponential covariances of four correlation lengths, fitted over 2–64 µm (open) and 16–64 µm (filled); the exponent moves with the window, and only the uncorrelated field returns 2 (diamond). (d) Recovered exponent under binomial thinning over three orders of magnitude of depth; the vertical axis is magnified to show the size of the variation.

**Fig. 4. Two independent paths give the same per-gene exponent.** (a) Per-gene exponent from aggregation against the covariogram exponent measured on the 8 µm grid, for the three sections; lines are Deming fits, and the numbers give Pearson and Spearman correlations and the median difference. (b) Local slope of the median covariogram against lag: flat in colon and lung, monotonically rising in mouse brain. (c) Pearson correlation between the two per-gene exponents for every combination of γ fit window and β lag window; the box marks the window fixed a priori.

**Fig. 5. Moving transcripts moves the exponent.** (a) Per-gene exponent of the six transcripts, observed (open circles) and after 20 independent random permutations of their positions (filled circles); the shaded bar is the 95% range over permutations, and the six genes are ordered by their observed value. Permutation leaves each gene's total, count distribution and occupancy unchanged, and moves all six towards the uncorrelated limit. (b) Pooled exponent after permuting the top N genes by μ² weight, N = 0, 10, 30, 100, 300 and 1,000; the percentages beside the points give the share of total μ² weight permuted.

**Fig. 6. The pooled exponent is carried by a few transcripts.** (a) Cumulative share of Σμ² against genes ranked by μ², for colon (dark blue), lung (mid blue) and mouse brain (ochre); dots mark the number of genes needed for half of the weight. (b) Exponent recomputed from the top-k genes only (solid circles) against the random-k 95% range (shaded), in the same colours, with the full-gene value as a dashed line; the first point of each curve is the heaviest transcript alone. (c) Largest change in exponent produced by leaving one gene out; the gene responsible is named below each bar.

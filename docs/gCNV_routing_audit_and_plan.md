# Batch-Aware gCNV Routing: Design Audit & Implementation Plan

Scope: WGS, WES, Illumina TruSight Hereditary Cancer (TSHC) · GATK 4.6 `GermlineCNVCaller` · Snakemake + Slurm
Status: design review, written before production code (2026-11-22)

---

## 0. Verdict at a glance

| # | Assumption | Verdict | Main change |
|---|-----------|---------|-------------|
| 1a | Autosome-only, segdup/PAR-excluded features | **Keep, extend** | Also exclude common polymorphic CNVs (GSTM1, UGT2B17, LCE3B/C, RHD, APOBEC3B…) and paralog/pseudogene regions. Otherwise they drive PCs and create "clusters" that are really genotypes. |
| 1b | Median norm → log2 → clip ±1.2/±1.5 → StandardScaler | **Change** | The clip doesn't do its job: het del = −1.0 and dup = +0.58 in log2 both sit inside the bounds. Replace with robust per-feature centering + z-winsorizing + a floored robust scale, and use robust (IRLS) projection at inference. |
| 2a | PCA (~8 PCs) + GMM, linear projection | **Keep, constrain** | Pick k_PC by parallel analysis. Use `tied`/`diag` covariance with `reg_covar`, or a Bayesian GMM. Enforce a minimum cluster size. |
| 2b | Mahalanobis 99th pct as OOD gate | **Insufficient** | Add the Q-statistic (SPE, the reconstruction error). A novel batch often varies *orthogonally* to the training PCs, where T² can't see it. Use Hotelling/F thresholds or leave-one-out (LOO) empirical thresholds, not in-sample χ² percentiles. |
| 3a | One genome-wide model per sample | **Keep** | chrX stays in the cluster model, since both sexes inform it. chrY stays in the cluster model when the cluster has enough males. Otherwise: pooled-male chrY model for WGS (not kit-adjusted), chrY ploidy only for WES/panel. |
| 3b | Mixed-sex cohorts + ploidy priors | **Keep, with guards** | Set minimum male/female counts per cluster. Split `FilterIntervals` into autosomes+X (all samples) and chrY (males only). |
| 4 | Fixed 24-chromosome scatter | **Keep for WGS/WES, reject for panel** | Panel: one shard. Per-chromosome panel shards hold a handful of targets, can't learn bias factors, and chrY yields an empty shard that crashes. The shard plan must come from data, not a constant. |
| 5 | Adaptive QS by CN × nInt | **Change direction** | QS grows with NP, so stricter thresholds for small NP wipe out single-exon calls. Calibrate thresholds empirically per (assay, DEL/DUP, NP bin) and add a review tier for reportable genes instead of hard filtering. |
| 5b | PAR1/PAR2 handled in postprocessing | **Move upstream** | Neither DGCP nor gCNV nor `PostprocessGermlineCNVCalls` has a PAR option, so handle PAR on the interval/count level. WGS/WES (SHOX required): relabel PAR1/PAR2 as diploid pseudo-contigs `chrPAR1`/`chrPAR2` and map them back to chrX after postprocessing (§1.3a). TSHC: exclude PAR. chrY PAR and XTR are always excluded. |
| — | Implicit: hard clusters beat a monolithic model | **Unproven** | gCNV's own bias factors are a linear batch model. The baseline to beat is monolithic **with more bias factors**, not the default. |

**Confirmed decisions (2026-11-23):**
- CASE calling is **single-sample**: each sample runs DGCP CASE + gCNV CASE alone against its routed model.
- **Bias factors are enabled** for WES and TSHC, in both COHORT and CASE (`--enable-bias-factors true`). WGS stays configurable.
- The **TSHC model is genome-wide, one shard**, with no chromosome scatter. WGS/WES keep per-chromosome shards.
- **SHOX calls are required** (WGS/WES). PAR1/PAR2 are modelled as in-silico pseudo-autosomes and reconstructed onto chrX after postprocessing (§1.3a).

---

## Part 1 — Critical Audit

### 1.1 The most important blind spots (read these first)

**(A) gCNV already contains a batch model, so the clustering partly duplicates it.**
In COHORT mode gCNV learns per-interval bias loadings (`--max-bias-factors`, default 5). Each sample gets latent factor scores, inferred again in CASE mode. A new batch whose bias lies in the span of the learned loadings is handled already. Routing helps only when:
- (i) batch effects are non-linear, or have more dimensions than the factor budget, or
- (ii) per-interval dispersion (ψ) and mean bias really differ between batches.

Consequences:
- The comparison arm in validation must be **monolithic with `--max-bias-factors` 10–20** (and ARD), not monolithic at the defaults.
- The current WGS configuration runs with `--enable-bias-factors false`. That WGS model has **no** latent batch correction at all, so in WGS the routing matters more than you may assume. Alternatively, turn bias factors on (memory cost).
- In [gatk_gCNV.cohort.smk](../snakemake/gatk_gCNV.cohort.smk), the comment `#--enable-bias-factors false #For WES data` is backwards. For WES and panels, bias factors are the main denoiser and must stay **on**.

**(B) Cluster size vs. homogeneity is the central trade-off, and the design doesn't model it.**
Every COHORT parameter per interval (mean bias, ψ, bias loadings) is estimated from the cluster's N samples. Going from N = 200 (heterogeneous) to N = 30 (homogeneous) cuts batch bias but raises estimation variance, which shows up as FP calls at noisy intervals. Clustering wins only where batch bias is larger than that added variance. This has to be **measured per assay**, not assumed (§2.4). An alternative that avoids the trade-off is *soft clusters*: train cluster k on its members plus its nearest neighbours from other clusters, padded up to N_min.

**(C) Cluster-specific `FilterIntervals` changes the callable territory per sample, and that is a clinical problem.**
If BRCA1 exon 2 survives filtering in cluster A but not in cluster B, B's patients get a silent no-call. The design needs:
- a **reportable-target BED** of must-call exons
- a per-cluster "callable fraction of reportable targets" metric, plus a release gate
- per-sample report text listing reportable exons that were not assessed

Never let a reportable exon disappear without trace.

**(D) The OOD gate is blind to new variation directions.**
Mahalanobis distance in the top-k PC space measures distance *within* the reference subspace. A new extraction kit typically adds a new signature (e.g., a different GC-dropout shape) that lies largely in the PCA residual. Standard multivariate process control uses **both** Hotelling T² (inside the subspace) and **Q/SPE** (‖x − x̂‖², outside the subspace). Without Q, the gate will pass exactly the samples it exists to catch.

**(E) Common CNV polymorphisms and paralogs will take over the PCA.**
In exomes, GSTM1/UGT2B17/LCE3B-C/RHD/APOBEC3B and similar loci are bi- or tri-modal at high allele frequency. Per-feature scaling (StandardScaler) gives them large weight, and GMM then splits on deletion genotype. Exclude intervals overlapping common CNVs (gnomAD-SV / DGV gold standard, AF > 1%), segdups, and known paralog pairs (PMS2/PMS2CL, CHEK2 pseudogenes, SMN1/2…) from the **clustering** features. As a QC step, check that cluster labels are independent of sex and of known polymorphic genotypes (χ² test).

### 1.2 Statistical & mathematical edge cases

**GMM covariance at N = 30–60, p = 8.**
- `full`: 36 covariance + 8 mean parameters per component. With N = 30 the estimate is unstable, and a single outlier can make it near-singular. Use `tied` (one shared covariance, sensible when batches differ mainly in location) or `diag`, with `reg_covar` ≈ 1e-3–1e-2 of the median PC variance.
- Alternatively use `BayesianGaussianMixture` (Dirichlet-process prior) so unsupported components collapse instead of overfitting.
- Choose K by BIC **subject to** min cluster size ≥ N_min(assay). Reject any solution with a cluster below N_min.
- HDBSCAN is good for *discovery* (it finds noise points and drift continua). Its `approximate_predict` is not a stable deployment-time projection, so deploy with GMM or kNN.
- Batches are often a **continuum** (reagent-lot drift), not discrete blobs. Require posterior ≥ 0.9 **and** a margin between the top two posteriors, else mark `AMBIGUOUS` (§2.2).

**Thresholds — in-sample percentiles are biased.**
- In-sample Mahalanobis distances are optimistically small because the sample being scored helped estimate the mean and covariance.
- For a new observation, T² follows `p(n+1)(n−1) / (n(n−p)) · F(p, n−p)`. With p = 8, n = 30, α = 0.01 the threshold is ≈ **37.6**. The naive χ²₈ 99% threshold is 20.1. A χ² or in-sample cutoff will therefore flag far too many routine samples as OOD.
- Use the F-based threshold, or LOO / cross-fitted empirical distances, with robust location/scatter (MCD, `sklearn.covariance.MinCovDet`) per cluster.
- Derive Q thresholds from the Jackson–Mudholkar approximation, or from LOO residuals.

**Feature pipeline.**
- *Clipping:* in log2, het del ≈ −1.0, hom del ≈ −6.6 (with the +0.01 pseudocount), dup ≈ +0.58, CN4 ≈ +1.0. A ±1.2/±1.5 clip only touches homozygous deletions and extreme outliers. It does **not** stop multi-exon CNVs from shifting PC scores.
  - Fix part 1: residualize first (subtract the per-feature reference median), then winsorize the robust z-score at ±3.
  - Fix part 2, the important one: **robust projection** at inference. Fit PC scores by IRLS (Huber/Tukey) regression of the sample's residual vector on the loadings. Contiguous blocks of deviating features (a CNV) are then down-weighted automatically. Report the down-weighted fraction as a QC metric.
- *StandardScaler:* scaling each feature to unit variance inflates near-constant, low-noise features and makes the clip boundary a discrete artifact. Use center + robust scale (MAD) **with a floor** (e.g., 25th percentile of feature MADs), or no scaling at all.
- *Median normalization* is robust to a 30 Mb event in WGS (≈1% of the genome) and to a whole-gene BRCA1 deletion on a panel (≈1% of targets). The median isn't the issue. PC scores are, which is what robust projection fixes.
- *Feature selection by assay:*
  - **Panel** (a few thousand targets): use all eligible targets. Expect 3–5 useful PCs at this N, not 8.
  - **WES** (~200k targets): stratified random subset of ~20k targets, balanced across chromosomes and GC/length strata. Fix the seed and store the subset in the registry.
  - **WGS** (~3M 1 kb bins): don't aggregate into large windows, because that blurs local GC-bias signatures. Instead combine (i) a stratified random subset of ~20k bins with (ii) a **technical fingerprint**: median residual log-depth per GC bin (e.g., 20 bins), per mappability bin and, for WES/panel, per target-length bin. Batch effects are mostly GC/length/fragment-size curves, and the fingerprint captures them almost free of CNV influence.
  - **Metadata features** (Picard insert size, duplication rate, fold-80, AT/GC dropout, % off-target) are cheap, CNV-free, and help explain clusters. Use them to interpret clusters and as a secondary OOD check, not as the primary embedding.
- *Choosing k_PC:* use parallel analysis (Horn) or a broken-stick test on the training matrix. A fixed "8" is arbitrary.
- *Large events / mosaic aneuploidy:* in **training**, a sample with a 30 Mb CNV adds a block of correlated features that can capture a PC. Exclude training samples with known large CNVs or aneuploidy (your karyotype filter already does this for WGS). Also run a *loading-locality check*: if more than X% of a PC's loading mass sits on one chromosome arm, it is biology, so reject that PC. In **inference**, robust projection plus the down-weighted-fraction QC handles it. Mosaic events are partial shifts, which the IRLS weights handle gracefully.

**Determinism.**
- Randomized SVD and GMM initialization need a fixed `random_state`. PC signs flip between refits (canonicalize the sign, e.g., max-abs loading positive).
- `joblib`/pickle is tied to library versions. Serialize the router as **plain arrays (NPZ) + JSON**: center, scale, loadings, component means/covariances, thresholds. The projection is a few lines of NumPy that any future version can reproduce and an auditor can read.
- Pin thread counts (`OMP_NUM_THREADS=1` for the router) so BLAS reductions are bitwise reproducible.

**Train/test leakage.** A training sample re-run in CASE mode against a model trained on it will look artificially good. Every benchmark must use samples held out from the cohort they are called against.

### 1.3 Sex chromosomes & cohort composition

**`FilterIntervals` and chrY.** The low-count filter (defaults: count < 5 in > 90% of samples) removes **all of chrY** once females exceed ~90% of a cluster. Below that point chrY survives, but its interval parameters are effectively learned from the males alone (females carry no chrY information). Five males in a cluster means a five-sample chrY model.

Fix (keeps chrX/chrY kit-adjusted wherever sample size allows):
1. Run `FilterIntervals` on **autosomes + chrX with all cluster samples**. chrX needs no separate model. gCNV models each sample's counts relative to its contig ploidy, so males (ploidy 1) and females (ploidy 2) both inform the chrX interval parameters. A separate non-kit-adjusted chrX model would only lose the kit correction.
2. Run `FilterIntervals` on **chrY with the cluster's males only**, then concatenate the two lists into the cluster's single filtered interval list. DGCP and gCNV then see chrY in every cluster whatever its sex ratio, so chrY ploidy stays callable (including XXY vs. XX).
3. Require a minimum number of males per cluster for chrY *segment* calling (e.g., ≥ 20). Below it, decide per assay:
   - **WGS:** a pooled-male chrY model across clusters. PCR-free WGS has small kit effects, so losing the kit adjustment costs little.
   - **WES/panel:** report chrY ploidy only, no chrY segments. Diagnostic yield of chrY CNV calls is low anyway, because the AZF regions are ampliconic and mostly removed by mappability filtering.

   If a pooled chrY model is used, run CASE + postprocessing for chrY separately and merge the VCFs.

**PAR1/PAR2: where the correction goes.**
None of the gCNV tools has a PAR option. `PostprocessGermlineCNVCalls --allosomal-contig` only sets each contig's baseline from the ploidy calls, for the whole contig. So the correction has to happen in the **interval list**, before any model sees PAR counts.

- **Where:** the interval list. For TSHC (and any assay with `par_mode: exclude`), pass a PAR (+ XTR) BED as `-XL` to `FilterIntervals`, or to `PreprocessIntervals` when the base list is built. For WGS/WES, chrX PAR intervals are *relabeled* to pseudo-contigs instead of excluded (§1.3a), while chrY PAR and XTR are still excluded. DGCP, gCNV COHORT and gCNV CASE all take the filtered list via `-L`, so PAR never enters the ploidy or CNV models. CASE mode must use the model's intervals anyway, so routed samples inherit the exclusion automatically. Use a plain `.bed` or `.interval_list`. In the current Snakefile, `par_bed_gz` is declared as an input of `determine_germline_contig_ploidy` but never used. The minimal fix is `-XL par.bed` in `filter_intervals`.
- **Why:**
  - GRCh38 no-alt analysis sets typically hard-mask chrY PARs, so all PAR reads align to chrX PAR. Males then show 2 copies of chrX PAR against a chrX ploidy of 1, giving a false DUP in every male. PAR reads also inflate the male chrX depth that DGCP uses to call ploidy.
  - If chrY PAR is *not* masked, PAR reads get MAPQ 0 and the default MAPQ ≥ 30 count filter drops them, giving false DELs or low-count filtering instead.
- **Postprocessing** keeps only a safety net: `filter-and-qc` suppresses any call overlapping PAR/XTR, in case an interval slipped through.
- **SHOX (PAR1):** excluding PAR1 would make SHOX deletions uncallable. WGS/WES therefore use the pseudo-autosome design in §1.3a.
- GRCh38 PAR coordinates (BED, 0-based): `chrX 10000 2781479`, `chrX 155701382 156030895`, `chrY 10000 2781479`, `chrY 56887902 57217415`.
- Also mask the **X-transposed region** (XTR, Xq21.3/Yp11.2, ~3.4 Mb, ~99% X–Y identity), a frequent source of male artifacts.

### 1.3a PAR1/PAR2 as in-silico pseudo-autosomes (SHOX)

**Decision:** for WGS and WES, PAR1 and PAR2 are modelled as two extra diploid contigs, `chrPAR1` and `chrPAR2`, and mapped back to chrX after postprocessing. TSHC has no PAR targets and keeps plain exclusion. Controlled by `par_mode: pseudo_contig | exclude` in `assay.json`.

**Precondition: PAR reads must land in one place.** The design works directly only if **chrY PAR is hard-masked** in the reference. Every PAR read then aligns to chrX PAR with normal MAPQ, and chrX PAR depth reflects X + Y PAR copies: 2 in both 46,XX and 46,XY.
- Check the FASTA:
  ```bash
  samtools faidx REF.fasta chrY:10001-2781479 | grep -v '>' | tr -d 'Nn\n' | wc -c   # 0 = masked
  ```
- Check empirically: median PAR1 bin depth ratio male/female ≈ 1.0 (X-specific bins ≈ 0.5).
- If chrY PAR is *not* masked, PAR reads get MAPQ 0 and the MAPQ ≥ 30 count filter drops them. You would then have to count both PARs at MAPQ ≥ 0 and sum X + Y per bin. Avoid this unless forced.

**Step 1: custom sequence dictionary.**
- Take the original dict and insert two contigs **after chr22 and before chrX**: `chr1 … chr22, chrPAR1, chrPAR2, chrX, chrY, …`.
- Keep the original coordinates, so this is a pure contig rename with no offset arithmetic. Set `LN` = chrX length for both.
- The position after chr22 matters: shards must be concatenated in dictionary order in postprocessing, and the PAR contigs share a shard with chr22 (step 4).

**Step 2: relabel intervals and counts.** Do this after `PreprocessIntervals`/`AnnotateIntervals`/`CollectReadCounts`, which run on the original reference, and before `FilterIntervals`.
- chrX intervals fully inside PAR1 become `chrPAR1`; fully inside PAR2 become `chrPAR2`. Intervals straddling a PAR boundary are dropped (with 1 kb bins one bin straddles chrX:2,781,479).
- chrY PAR intervals are dropped (masked, zero counts). XTR stays excluded.
- Apply the same relabel to the preprocessed interval list, the annotated intervals and **every sample's count file**. GATK CNV tools require identical sequence dictionaries across count files and annotated intervals. The relabel works on TSV; convert HDF5 counts to TSV first or write GATK-compatible HDF5.
- Everything downstream (`FilterIntervals`, DGCP, gCNV COHORT/CASE) works only on counts and intervals and never reads the FASTA.

**Step 3: ploidy priors.** PAR copy number = X count + Y count: 1 in 45,X; 3 in 47,XXY, XXX and XYY. Use a slightly wider prior than for autosomes:
```
chrPAR1	0.01	0.02	0.95	0.02
chrPAR2	0.01	0.02	0.95	0.02
```

**Step 4: sharding.** PAR1 is ~2.7 Mb; PAR2 is ~330 kb (a few hundred 1 kb bins, only a handful of WES exons). That is too few intervals to learn bias factors from on their own in WES. Put both pseudo-contigs **into the chr22 shard** for WGS and WES. A gCNV shard may span several contigs.

**Step 5: postprocessing.**
- Do **not** list `chrPAR1`/`chrPAR2` as `--allosomal-contig`. They then get the reference copy number from `--autosomal-ref-copy-number 2` in every sample, so a 45,X sample's single PAR copy shows as CN1, which is biologically correct.
- Pass the custom dictionary (`--sequence-dictionary pseudo.dict`) instead of `--reference`, because the FASTA has no `chrPAR1/2`. Verify on a test sample that the VCF header comes out as expected.

**Step 6: reconstruction** (`gcnvr par-reconstruct`, pysam).
- For the segments VCF, intervals VCF and denoised copy ratios: rename `chrPAR1/2 → chrX` (coordinates are already original), replace the header contigs with the original reference dictionary, re-sort and index.
- Add INFO `PAR=PAR1|PAR2` and `REFCN` (2 for PAR records, the chrX ploidy call for X-specific records), so downstream tools and analysts know the baseline changes inside chrX.
- **Do not merge segments across a PAR boundary automatically.** A deletion crossing chrX:2,781,479 comes out as two adjacent records, e.g. PAR CN1 vs. ref 2 plus X-specific CN0 vs. ref 1 in a male. Flag both `PAR_BOUNDARY_ADJACENT` for review; the analyst reports it as one event.
- For plots, rescale denoised copy ratios to the record's `REFCN`, so PAR and X-specific regions share one axis.

**Step 7: whole-PAR events, the ploidy trap.** DGCP calls a ploidy per contig per sample. A deletion covering most of PAR1 (terminal Xp/Yp deletions, X;Y rearrangements) can make DGCP call `chrPAR1` ploidy 1. gCNV then works from that baseline, and the loss may not appear as a normal segment. Two guards:
- **QC rule:** expected ploidy(chrPAR1) = ploidy(chrX) + ploidy(chrY), and the same for chrPAR2. A mismatch produces an explicit "whole-PAR copy-number change" finding in the report, not a pipeline failure. It is also a sample-swap and sex-check cross-check.
- Typical SHOX and enhancer deletions span tens to a few hundred kb of the 2.7 Mb PAR1. With a 0.95 prior on ploidy 2 they should not flip the contig ploidy. Verify with the positive controls below, and run a 45,X sample to confirm how postprocessing labels a contig whose ploidy is 1.

**Limits.**
- WES sees only SHOX exons. Enhancer-only deletions downstream of SHOX need WGS (or MLPA).
- Validation set: MLPA-confirmed SHOX gene deletions, enhancer-only deletions, SHOX duplications, 45,X, 47,XXY, 47,XYY.
- Negative check: false-call rate in PAR must not differ between sexes.

**Karyotype edge cases in routine diagnostics.**
- 47,XXY / 45,X / 47,XYY: the priors allow them (`chrX` ploidy 3 = 0.01; consider a small non-zero prior for `chrY` ploidy 2).
- **Mosaic loss of Y** is common in blood DNA from older men, and hereditary-cancer patients skew older. It produces ambiguous chrY ploidy and noisy chrY training data. Exclude males with low chrY depth ratio from chrY training, and don't gate sample QC on chrY.
- Sex-check discordance (ploidy call vs. LIMS sex) should block the report. It is also your best sample-swap detector.

**Panel specifics.**
- `DetermineGermlineContigPloidy` is unstable for contigs with few targets.
- If TSHC has few or no chrX/chrY targets, derive the priors table **from the contigs actually present** (don't feed priors for absent contigs). Use high-confidence autosomal priors, since aneuploidy detection isn't the goal of a panel.
- Get sex from LIMS / an orthogonal check (off-target reads, SNP heterozygosity) rather than relying on DGCP for X/Y.

### 1.4 Adaptive QS vs. single-exon sensitivity

- **QS adds up evidence across intervals.** An NP = 1 event has a structurally capped QS. Thresholds that *rise* as nInt falls (the usual "be strict with small events" reflex) systematically remove true single-exon deletions, which are exactly the clinically decisive BRCA1/2, MLH1, MSH2 calls.
- **DUP and DEL are asymmetric:** CR 1.5 vs. 0.5 means about half the log-fold change at equal depth. DUP needs its own threshold curve, and single-exon dups on a panel will usually stay in a review tier.
- **Pseudogene/paralog interference:** `CollectReadCounts` drops reads below MAPQ 30 by default.
  - PMS2 exons 12–15 (and 9, 11) vs. PMS2CL, CHEK2 3′ exons vs. its pseudogenes, and BRCA1 exons 1–2 vs. BRCA1P1 lose counts in a mapping-dependent way. The result is either filtering by `FilterIntervals` or spurious deletions.
  - These regions should be classified **"not assessable by gCNV → orthogonal method (MLPA / long-range PCR)"** in the reportable-target BED.
  - Optional advanced step: paralog-aware summed counting (MAPQ ≥ 0, PMS2+PMS2CL combined), reported as a combined copy number.
- **Recommended policy, three tiers, nothing deleted:**
  1. *Report candidate:* passes a threshold calibrated per (assay, DEL/DUP, NP bin ∈ {1, 2, 3–5, 6–20, > 20}), chosen to hit a pre-set sensitivity (e.g., ≥ 95% for NP = 1 DEL in reportable genes) on the truth set.
  2. *Review:* overlaps a reportable target, NP ≤ 2, QS between a floor and the tier-1 threshold, **and** denoised CR is consistent (e.g., < 0.65 / > 1.35), **and** the target is not paralog-flagged. Confirm by MLPA.
  3. *Suppressed:* logged with its reason.
- Add a **reference-cohort site-frequency filter**: an event seen at the same intervals in > X% of the *training* cohort is an artifact or polymorphism. Whitelist known founder CNVs (e.g., BRCA1 exon 13 dup, CHEK2 del exons 9–10) so a recurrent true event isn't filtered.
- `--p-alt` (default 1e-6) is a second sensitivity knob for single-target events. Tune it together with QS on the truth set, not on its own.

### 1.5 24-chromosome scattering

- **WGS/WES:** per-chromosome is reasonable. `PostprocessGermlineCNVCalls` concatenates shards and segments across them, so the "mid-arm HMM reset" mainly affects posteriors near boundaries, not final segment stitching. The real benefit is simpler bookkeeping and no arbitrary cut points.
- gCNV is not a 30 Mb event caller. Sample-level factor scores can partially absorb very large events. Back it up with an arm-level / ploidy-level depth check (or a dedicated large-CNV tool) and compare raw vs. denoised CR for events > 5 Mb.
- **Panel: use one shard.** Per-chromosome panel shards contain anywhere from 0 to a few hundred targets. Bias factors and ψ cannot be learned from a few targets, and an empty chrY shard crashes the run. The current `scatter_filtered_intervals` with a hard-coded `CHROMS` list will write header-only interval lists for absent contigs.
- **Memory:** scale it with the shard's interval count from a manifest, using a fitted model (`mem = a + b·n_intervals [+ c·n_samples·n_intervals]` for COHORT). Calibrate it from `sacct` MaxRSS of past runs instead of hard-coding chromosome names. Use `attempt`-based escalation on OOM.

### 1.6 Operational & lifecycle gaps

**Cold start (new kit/sequencer, < N_min samples).** Never silently route to the "least-bad" cluster. Policy ladder:
1. **Change control first** (ISO 15189 spirit): a wet-lab change is announced, and a validation batch (≥ 24–30 samples, including known positives) is run before go-live.
2. Samples that pass T² **and** Q for an existing cluster → route normally. The change was harmless.
3. OOD samples → `PROVISIONAL`: call with the nearest cluster **and** with the pooled monolithic high-factor model, report only concordant tier-1 calls, and send everything to review / orthogonal confirmation.
4. Once ≥ N_min new-batch samples accumulate → train a new cluster model and run a release (§2.1).

**Drift.** Track per-run T², Q, posterior entropy, AMBIGUOUS/OOD rate, calls/sample and denoised-CR MAD on **control charts** (EWMA/CUSUM). Alarms trigger review, not automatic retraining.

**Retraining & reproducibility.**
- Releases are **immutable**. Cluster IDs are never reused. Map new to old clusters by Hungarian matching of centroids, for continuity reports only.
- Every result carries `release_id`, `router_id`, `cohort_id`, interval-list SHA-256, GATK image **digest**, and the parameter-profile hash.
- Retain old models for the clinical retention period so reanalysis can reproduce the original call.

### 1.7 Findings in the current repo code

- [gatk_gCNV.cohort.smk](../snakemake/gatk_gCNV.cohort.smk): `--enable-bias-factors false` is fine for WGS memory but must not carry over to WES/panel. The WES comment says the opposite.
- `CHROMS` is hard-coded to 24 and used for scattering, which breaks for panels (empty shards).
- `par_bed_gz` is an unused input. PAR handling currently depends on the preprocessed interval list.
- Memory lambda keys on chromosome names, not interval counts.
- Parameters (`class/cnv-coherence-length 100000`) are WGS-1 kb choices. WES/panel should start from GATK defaults (10000). Make parameters an **assay profile** in the registry, not literals in the rule.
- `srun` + `docker run -u root` inside `shell:` hides jobs from Snakemake's scheduler and needs root on compute nodes. Consider the Snakemake ≥ 8 Slurm executor plugin with Apptainer images pinned by digest.

---

## Part 2 — Implementation Plan

### 2.1 Artifact & model registry

```
registry/
└── assays/
    └── TSHC_v1/                                  # assay = enrichment design + reference build
        ├── assay.json                            # build, kit, profile name, N_min, shard strategy, par_mode
        ├── intervals/
        │   ├── base.preprocessed.interval_list   # PAR/XTR excluded at PreprocessIntervals
        │   ├── base.annotated.tsv                # GC, mappability, segdup
        │   ├── reportable_targets.bed            # must-call exons; col: assessable={gcnv,orthogonal}
        │   ├── clustering_features.bed           # autosomal, −segdup −commonCNV −paralog −PAR
        │   ├── masks/{par,xtr,common_cnv,paralog,founder_whitelist}.bed
        │   └── par/{pseudo.dict,par_map.tsv}     # WGS/WES only (par_mode=pseudo_contig); par_map: pseudo-contig ↔ chrX start/end
        ├── profiles/
        │   └── gcnv_params.v3.yaml               # COHORT/CASE flags, QS calibration table ref
        ├── routers/
        │   └── R-2026.11.0/
        │       ├── router.json                   # metadata + thresholds (schema below)
        │       ├── router.npz                    # center, scale, loadings, means, covs
        │       ├── feature_index.tsv             # exact ordered feature list (+ sha256)
        │       ├── training_samples.tsv          # sample, cluster, posterior, T2, Q, sex, batch metadata
        │       └── report.html                   # scree, PC plots, cluster×metadata tables
        ├── cohorts/
        │   ├── C-2026.11.0-k01/
        │   │   ├── cohort.json
        │   │   ├── samples.tsv
        │   │   ├── filtered.autosomes_chrX.interval_list
        │   │   ├── shards/shard_manifest.tsv     # shard_id, contigs, n_intervals, mem_mb_case, mem_mb_cohort
        │   │   ├── shards/shard_000.interval_list ...
        │   │   ├── ploidy-model/
        │   │   ├── gcnv-model/shard_000/ ...
        │   │   ├── calibration/qs_thresholds.tsv # assay × DEL/DUP × NP-bin → QS threshold
        │   │   ├── reference_site_freq.tsv.gz    # per-interval event freq in training calls
        │   │   └── qc/callable_reportable.tsv    # per reportable exon: kept/filtered
        │   └── C-2026.11.0-chrY-pooledMale/      # optional: WGS clusters with too few males
        └── releases/
            ├── 2026.11.0.json                    # pins router + cohorts + images + profile
            └── CURRENT -> 2026.11.0.json
```

**`router.json` contract (abridged):**
```json
{
  "schema": "gcnv-router/1",
  "router_id": "R-2026.11.0", "assay": "TSHC_v1", "created": "2026-11-22",
  "code_version": "gcnvrouter 0.4.1 (git 4719c9d)",
  "features": {"bed": "clustering_features.bed", "sha256": "…", "n": 2841,
               "normalization": "median_autosomal", "pseudocount": 0.01,
               "center": "reference_median", "scale": "mad_floored", "winsor_z": 3.0},
  "pca": {"k": 4, "selection": "parallel_analysis", "sign_convention": "maxabs_positive"},
  "mixture": {"type": "gmm", "covariance": "tied", "reg_covar": 1e-3, "K": 3, "seed": 1729},
  "projection": {"method": "irls_huber", "c": 1.345, "max_downweighted_frac": 0.15},
  "gates": {"posterior_min": 0.90, "margin_min": 0.20,
            "t2": {"method": "hotelling_F", "alpha": 0.01, "per_cluster": {"k01": 31.2}},
            "q":  {"method": "loo_empirical", "quantile": 0.99, "value": 412.7}},
  "clusters": {"k01": {"cohort_id": "C-2026.11.0-k01", "n": 58, "n_male": 21, "n_female": 37}},
  "files": {"router.npz": "sha256:…", "feature_index.tsv": "sha256:…"}
}
```

**`cohort.json`** pins: sample list hash, `FilterIntervals` arguments, filtered-list SHA-256, shard manifest hash, ploidy priors hash, profile hash, GATK image digest, `callable_reportable_fraction`, and `release_gate` results. **`releases/*.json`** contains only pointers plus hashes. `registry validate` re-hashes everything and refuses to run on a mismatch.

### 2.2 Python package `gcnvrouter`

```
src/gcnvrouter/
├── io.py          # read CollectReadCounts TSV/HDF5 → (sample, intervals, counts); interval-key alignment
├── features.py    # build_feature_mask(), normalize(), residualize(), winsorize(), fingerprint_gc()
├── router.py      # fit_router(), project_robust(), score(T2,Q,posterior), save_npz/json, load
├── cohort.py      # select_cohort(), sex_balance_check(), plan_filter_intervals(), plan_shards()
├── calibrate.py   # QS threshold fitting per (DEL/DUP, NP bin) from truth calls
├── postprocess.py # tiering, site-freq filter, PAR/XTR suppression, reportable coverage, VCF annotate
├── qc.py          # sample gates; run-level control-chart stats
├── registry.py    # manifests, hashing, schema validation (pydantic/jsonschema), release/CURRENT
└── cli.py         # Typer entry point `gcnvr`
```

**CLI commands (offline):**

| Command | Input → Output | Key behaviour |
|---|---|---|
| `gcnvr build-features --assay A --counts DIR` | counts → `clustering_features.bed`, `feature_index.tsv` | Autosomes only; removes segdup/common-CNV/paralog/PAR. Stratified subset for WES/WGS with fixed seed. |
| `gcnvr build-router --assay A --samples S.tsv` | → `routers/R-*/` | Excludes known large-CNV/aneuploid samples. Parallel analysis → k. Loading-locality check. GMM/BGMM with `tied`/`diag`, BIC with N_min constraint. MCD per cluster. Hotelling-F T² and LOO Q thresholds. Sex/genotype independence test. HTML report. |
| `gcnvr train-cohort --router R --cluster k01 [--knn-pad N_min]` | → `cohorts/C-*/` + Snakemake config | Writes the `FilterIntervals` plan: autosomes+X on all samples, chrY on the cluster's males, PAR/XTR excluded via `-XL`. Checks the male count for chrY segment calling (else pooled-male chrY for WGS, ploidy only for WES/panel). Sets the assay profile (bias factors on for WES/panel). Computes `callable_reportable`. `plan_shards` (WGS/WES per contig, panel single shard, empty contigs dropped). Fills in memory estimates. |
| `gcnvr calibrate-qs --cohort C --truth T.vcf` | → `qs_thresholds.tsv` | Per (DEL/DUP, NP bin): smallest QS that keeps target sensitivity. Bootstrap CI. |
| `gcnvr release --assay A --router R --cohorts …` | → `releases/X.json` | Runs release gates (callable fraction = 100% of gCNV-assessable reportable targets, benchmark pass) and flips CURRENT. |

**CLI commands (online):**

| Command | Output |
|---|---|
| `gcnvr par-relabel --assay A --counts s.tsv` | relabeled counts on `pseudo.dict` (no-op for `par_mode: exclude`). The same function relabels the interval lists offline. |
| `gcnvr route-sample --release X --counts s.tsv` | `route.json` (below) |
| `gcnvr par-reconstruct --assay A --segments … --intervals … --dcr …` | chrX-reconstructed VCFs/TSV with original header, `PAR`/`REFCN`/`PAR_BOUNDARY_ADJACENT` tags, PAR ploidy consistency check |
| `gcnvr filter-and-qc --route route.json --segments seg.vcf.gz --intervals int.vcf.gz --dcr dcr.tsv` | annotated VCF (`TIER`, `FILTER` reason, `REPORTABLE`, `ASSESSABLE`), `sample_qc.json`, `not_assessed_targets.tsv` |
| `gcnvr qc-run --run RUN` | run-level control-chart record; alarms |

```json
{"sample": "S1", "release_id": "2026.11.0", "router_id": "R-2026.11.0",
 "status": "ROUTED",                     // ROUTED | AMBIGUOUS | OOD | QC_FAIL
 "cluster": "k01", "cohort_id": "C-2026.11.0-k01",
 "posterior": 0.97, "margin": 0.91,
 "t2": 12.4, "t2_threshold": 31.2, "q": 188.0, "q_threshold": 412.7,
 "downweighted_frac": 0.012,
 "fallback": null}                        // for OOD/AMBIGUOUS: {"cohort_id": "C-…-monolithic", "mode": "PROVISIONAL"}
```

`filter-and-qc` gates per sample:
- route status
- ploidy-vs-LIMS sex
- PLOIDY_GQ
- calls/sample ≤ cohort 99th percentile
- denoised-CR MAD ≤ cohort 99th percentile
- IRLS down-weighted fraction
- callable reportable fraction

Any failure → `QC_FAIL` with reasons, and no automatic report.

### 2.3 Snakemake + Slurm architecture

**Recommendation: two phases with a frozen routing table, instead of checkpoints.**

1. **Phase A, `route.smk`:** `CollectReadCounts` → `route-sample` per sample → aggregate `run/routing_table.tsv` (sample, status, cohort_id, release_id). This is a small, fast DAG. The table is a reviewable, signable artifact: you can see before calling which model each patient will get.
2. **Phase B, `call.smk`:** reads `routing_table.tsv` at parse time, so the DAG is fully static and `--dry-run` shows every job. Shard lists come from each cohort's `shard_manifest.tsv`.

This is easier to audit and debug than checkpoint re-evaluation. If you need a single invocation, the same input functions work behind a `checkpoint aggregate_routes` (sketch at the end of this section).

```python
# call.smk (sketch)
import pandas as pd, json
from gcnvrouter.registry import Registry

REG    = Registry(config["registry"])
ROUTES = pd.read_csv(config["routing_table"], sep="\t").set_index("sample")
CALLABLE = ROUTES[ROUTES.status.isin(["ROUTED", "AMBIGUOUS", "OOD"])]   # QC_FAIL stops here

def cohort(wc):  return REG.cohort(ROUTES.loc[wc.sample, "cohort_id"])
def shards(wc):  return cohort(wc).shard_ids                 # ['000'] for panel; per-contig for WGS/WES

rule all:
    input: expand("results/{sample}/{sample}.gcnv.final.vcf.gz", sample=CALLABLE.index)

rule ploidy_case:
    input:  counts="counts/{sample}.counts.hdf5",
            model=lambda wc: cohort(wc).ploidy_model,
            intervals=lambda wc: cohort(wc).filtered_intervals
    output: directory("work/{sample}/ploidy-calls")
    resources: mem_mb=8000, runtime=30
    container: config["gatk_image"]          # pinned by digest
    shell: "gatk DetermineGermlineContigPloidy --run-mode CASE --model {input.model} "
           "-I {input.counts} -L {input.intervals} --interval-merging-rule OVERLAPPING_ONLY "
           "--output work/{wildcards.sample} --output-prefix ploidy"

rule gcnv_case:
    input:  counts="counts/{sample}.counts.hdf5",
            ploidy="work/{sample}/ploidy-calls",
            model=lambda wc: cohort(wc).shard_model(wc.shard)
    output: directory("work/{sample}/gcnv/shard_{shard}-calls")
    threads: 4
    resources:
        mem_mb=lambda wc, attempt: cohort(wc).shard_mem_case(wc.shard) * (1 + 0.5 * (attempt - 1)),
        runtime=lambda wc: cohort(wc).shard_runtime_case(wc.shard)
    retries: 2
    container: config["gatk_image"]
    shell: "gatk GermlineCNVCaller --run-mode CASE --model {input.model} -I {input.counts} "
           "--contig-ploidy-calls {input.ploidy} --output work/{wildcards.sample}/gcnv "
           "--output-prefix shard_{wildcards.shard} {params.profile_flags}"

rule postprocess:
    input:  calls=lambda wc: expand("work/{s}/gcnv/shard_{sh}-calls", s=wc.sample, sh=shards(wc)),
            models=lambda wc: [cohort(wc).shard_model(sh) for sh in shards(wc)],
            ploidy="work/{sample}/ploidy-calls"
    output: seg="work/{sample}/segments.vcf.gz", ints="work/{sample}/intervals.vcf.gz",
            dcr="work/{sample}/dcr.tsv"
    # only if the cohort uses a pooled-male chrY model (WGS, too few males):
    # separate CASE + postprocess for chrY, merged in `finalize`

rule filter_and_qc:
    input:  route="routes/{sample}.json", seg=rules.postprocess.output.seg,
            ints=rules.postprocess.output.ints, dcr=rules.postprocess.output.dcr
    output: vcf="results/{sample}/{sample}.gcnv.final.vcf.gz",
            qc="results/{sample}/{sample}.qc.json"
    shell:  "gcnvr filter-and-qc --route {input.route} --segments {input.seg} "
            "--intervals {input.ints} --dcr {input.dcr} -o {output.vcf} --qc {output.qc}"
```

Notes:
- **Slurm:** `snakemake --executor slurm` (plugin) with `default-resources: slurm_partition=all`. Remove `srun`/`docker` from `shell:` so Snakemake sees real resource requests, retries, and accounting. Use Apptainer (`--sdm apptainer`) to avoid root containers.
- **Shard memory** comes from the manifest. `plan_shards` writes it from a fitted model, and `retries` + `attempt` escalate on OOM. Record MaxRSS after each run to refit.
- **PAR rules (WGS/WES):** `CollectReadCounts` → `par_relabel` → `route_sample` / `ploidy_case` / `gcnv_case` → `postprocess` (with `--sequence-dictionary pseudo.dict`, PAR contigs not allosomal) → `par_reconstruct` → `filter_and_qc`. In the shard manifest, the PAR contigs belong to the chr22 shard.
- **Single-sample CASE (decided):** every sample runs DGCP CASE and gCNV CASE alone against its routed model. Routing is per sample, a failing sample never blocks others, and `--sample-index` in postprocessing is always 0. Cost: the model is loaded once per sample × shard. That is noticeable for WGS (24 shards) and negligible for TSHC (1 shard).
- **COHORT training** gets its own workflow, `train.smk`, driven by `train-cohort` output: FilterIntervals → DGCP COHORT → gCNV COHORT per shard → calibration → release. This replaces today's `gatk_gCNV.cohort.smk` and reuses its rule bodies with assay profiles.
- **Checkpoint variant** (single invocation): `checkpoint aggregate_routes` writes the routing table, and `cohort(wc)` reads `checkpoints.aggregate_routes.get().output[0]` instead of the parse-time table. Everything else stays the same.

### 2.4 Validation & benchmarking

**Arms (same samples, same truth):**
- **A0** monolithic COHORT, default bias factors
- **A1** monolithic, `--max-bias-factors` 10–20 (the real competitor)
- **B** cluster-routed hard clusters
- **B′** cluster-routed with kNN padding to N_min

All arms call single-sample CASE against their model.

**Truth sets:**

| Assay | Real truth | Synthetic |
|---|---|---|
| WGS | GIAB HG002 SV/CNV benchmarks (deletions); in-house samples with chromosomal microarray (CMA) results (large CNVs, aneuploidy, mosaics); Coriell lines with known events | Binomial thinning (DEL) / count augmentation (DUP) of 30 kb–30 Mb blocks |
| WES | In-house CMA/MLPA-confirmed exon CNVs | Same, exon-level, NP 1–20 |
| WGS/WES PAR | MLPA-confirmed SHOX gene and enhancer-only deletions, SHOX dups; 45,X / 47,XXY / 47,XYY | Thinning/augmentation of PAR1 blocks in both sexes |
| TSHC | In-house MLPA positives (BRCA1/2, MLH1, MSH2, MSH6, PMS2, EPCAM, CHEK2), incl. founders (BRCA1 ex13 dup, CHEK2 del ex9–10); **ICR96** (MLPA-validated exon CNVs on TruSight Cancer data, controlled access) for method-level comparison | Read-level spike-ins (downsample/duplicate reads in target exons of real BAMs, then rerun `CollectReadCounts`), which preserve real batch noise |

**Designs:**
- **Leave-batch-out:** train on all batches but one, call the held-out batch. This tests both OOD detection (is it flagged?) and degradation (how much worse when wrongly routed?).
- **Cross-fitting:** cases are never members of the cohort they are called against.
- **Prospective shadow mode:** run B alongside the current pipeline for ≥ 3 months, and review discordant calls blinded to arm.

**Metrics, stratified by assay × DEL/DUP × NP bin {1, 2, 3–5, 6–20, > 20} × reportable/non-reportable:**
- Event-level sensitivity and PPV. Matching: ≥ 50% reciprocal overlap, or ≥ 1 shared target for NP ≤ 2.
- **FP calls per sample.** This drives review workload and is often the decisive clinical metric.
- **Large CNV retention (WGS):** detection, size ratio, and breakpoint distance for ≥ 1 Mb events. Also check for fragmentation (one true event called as many segments).
- **Denoised-CR MAD** in diploid regions per sample. It is continuous, so it is far more powerful than event counts for detecting a batch-noise effect with modest N.
- Calibration: empirical precision vs. QS curves per arm.
- Callable reportable fraction per arm (B can lose territory through cluster-specific filtering).

**Statistics & acceptance (pre-register before unblinding):**
- **Sensitivity:** paired McNemar per event (A1 vs. B), with sample-level cluster bootstrap CIs. Non-inferiority margin e.g. −2 pp overall, and **zero loss** of reportable-gene true positives.
- **FP/sample:** Wilcoxon signed-rank on paired counts. B must reduce FP/sample with CI excluding 0.
- **Sample size:** with zero misses, a one-sided 95% lower bound ≥ 90% needs ≥ 30 true events per stratum (two-sided: ≥ 36). For NP = 1 strata this will mostly need read-level spike-ins on top of real positives.
- **Reproducibility:** rerun the same input twice and on two nodes; the VCF body hash must be identical.

**Decision rule per assay:**
- Adopt B (or B′) only if it beats **A1** on FP/sample and denoised-CR MAD without losing sensitivity.
- Expected outcome to test: likely a clear win on panels/WES with heterogeneous wet-lab history; possibly no gain on PCR-free WGS once bias factors are on.

### 2.5 Suggested build order

1. Registry + `io` + `registry validate`. Turn `gatk_gCNV.cohort.smk` into `train.smk` with assay profiles and a data-driven shard plan (fixes the panel/chrY breakage now).
2. `build-features` + `build-router` + HTML report on historical data. Inspect clusters against metadata **before** training any cluster model.
3. `train-cohort` for A1 and B on TSHC (the cheapest compute, the most clinical truth).
4. `calibrate-qs`, `filter-and-qc`, and the benchmark harness. Run §2.4 on TSHC, then WES, then WGS.
5. `route.smk` + `call.smk` in shadow mode, with control charts.
6. Release gate, SOP for cold start/drift, go-live per assay.

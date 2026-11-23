# Plan: Batch Clustering, Router & Visualization

Status: **Phase 1 (exploration) implemented** in `src/gcnvrouter/` (2026-11-23). Phases 2–4 are still a plan.
Parent design: [gCNV_routing_audit_and_plan.md](gCNV_routing_audit_and_plan.md), §1.1, §1.2, §2.1 and §2.2.
Scope: everything up to and including `route-sample`. Cohort training per cluster, CASE calling and QS filtering are out of scope here.

---

## Usage (phase 1)

```bash
mamba env create -f envs/gcnvrouter.yaml && mamba run -n gcnvrouter pip install -e .
gcnvr explore --counts-dir gCNV/03_read_counts/tshc --pattern '*.hg38.tsv' --profile panel -o explore_tshc \
    [--metadata sample_metadata.tsv] [--annotated hg38.annotated.tsv] [--mask-bed common_cnv.bed] \
    [--exclude SAMPLE] [--n-min 3] [--n-boot 100]
# without installing: PYTHONPATH=src python -m gcnvrouter.cli explore ...
```

Outputs:
- `report.html`: self-contained, works offline.
- `samples.tsv`: per sample, the preflight result, PC scores, cluster, posterior/margin, T²/Q against their limits, IRLS down-weighting and status.
- `sweep.tsv` (BIC per K and covariance type), `pca.tsv` (eigenvalues + parallel-analysis threshold), `locality.tsv`, `features.tsv` (features + loadings) and `summary.json`.

Tests: `python -m pytest tests/router` (synthetic data only).

**Implementation differences from the plan text below:**
- **Stability uses 80% subsampling without replacement** (consensus clustering), not the classic bootstrap. With small cohorts a bootstrap resample holds only ~63% unique samples, so small real clusters often cannot be refitted.
- **Exploration minimum cluster size is 3 for panels** (WES/WGS 10). Trainability is reported separately against the production minimum (30/50/30).
- **Preflight adds an assay-compatibility check.** A sample with more than 5% of targets at < 5% of median depth is marked `ASSAY_MISMATCH` and excluded from the fit. It typically means a different design or probe pool, not a batch.
- **Down-weighting** is the fraction of features with Huber weight < 0.5 (|u| > 2c), about 0.7% under Gaussian noise. Q is computed over the remaining features, so a large CNV does not inflate it.
- The **loading-locality check** works per chromosome, not per arm, because no centromere table is bundled yet.
- **Catch-all warning:** if a cluster is > 3× more spread out than the tightest one, the report warns that it is probably several batches or outliers merged.

---

## 0. Guiding decisions

1. **Explore before building the router.** Phase 1 answers one question with plots: *is there real wet-lab batch structure in our historical counts, and does it line up with known metadata?* If the answer is "no structure beyond run-to-run noise" (K = 1), stop. A monolithic model with more bias factors is then the right answer, and the routing machinery isn't needed.
2. **K = 1 is a valid result.** The model selection must be allowed to choose a single cluster.
3. **One code path for WGS, WES and TSHC.** Only the feature-selection and threshold parameters differ, and they live in an assay profile.
4. **Plain, auditable artifacts.** The router is stored as NPZ + JSON. No pickle.
5. **Reports are self-contained HTML** (plotly with plotly.js embedded, so they work offline on the clinical network), plus static PNGs for SOPs and validation reports.
6. **UMAP is for exploration only.** It is an optional plot and never used for routing.

---

## 1. Inputs

| Input | Source in this repo / infrastructure | Needed for |
|---|---|---|
| Read counts | `CollectReadCounts` TSV (`*.hg38.tsv`, `@` header + `CONTIG START END COUNT`) | features |
| Annotated intervals | `assets/{protocol}/hg38.annotated.interval_list` (columns `GC_CONTENT`, `MAPPABILITY`, `SEGMENTAL_DUPLICATION_CONTENT`) | feature eligibility, GC fingerprint, loading interpretation |
| Masks | `common_cnv.bed` (gnomAD-SV / DGV gold standard, AF > 1%), `paralog.bed`, PAR coordinates (from `scripts/par/par_common.py`) | feature eligibility |
| Run metadata | parsed from run folder names, e.g. `240913_A01917_0057_AHN3M5DSXC` → date, instrument, run number, flowcell | colouring, contingency tests, drift |
| QC metrics | GRZ QC files (`quality_control/GRZ/G*.hg38.final.txt`): mean depth plus any other columns present | colouring, secondary interpretation |
| Sex | karyotype / selection step (`sex` column), or DGCP ploidy calls | independence check |
| **Wet-lab metadata (LIMS)** | user-provided TSV: extraction kit, library-prep kit/lot, DNA source (blood/saliva), input amount, operator, sequencer | **the key validation input**: clusters can only be shown to be wet-lab batches if these exist |

The metadata contract is a single `sample_metadata.tsv`, keyed by `sample_id`. Columns are free-form, but these names are recognised: `run`, `run_date`, `instrument`, `flowcell`, `sex`, `coverage`, `extraction_kit`, `library_kit`, `library_lot`, `dna_source`. A helper builds the automatic part (run, QC, sex), and LIMS columns are merged in from a user TSV.

---

## 2. Package layout

Located in `src/gcnvrouter/` as in the parent plan. It runs in its own conda env, `envs/gcnvrouter.yaml`: python 3.11, numpy, scipy, pandas, scikit-learn, plotly, matplotlib, pyyaml, jinja2, pytest, and optionally umap-learn. It is run with the existing `mamba run -n …` pattern, **not** in the GATK container (that container's Python lacks sklearn and plotly). A Dockerfile can come later.

```
src/gcnvrouter/
├── profiles.py      # assay profiles (wgs / wes / panel): feature caps, thresholds, k range
├── io_counts.py     # fast TSV reader (pandas, comment='@', usecols, float32); interval alignment check
├── metadata.py      # run-name parser, GRZ QC reader, LIMS merge → sample_metadata.tsv
├── features.py      # eligibility mask, stratified subsample, normalization, GC fingerprint
├── cache.py         # per-sample feature vectors (.npy) keyed by feature-index hash → incremental
├── pca.py           # SVD, parallel analysis, sign convention, loading-locality check, IRLS projection
├── cluster.py       # GMM/BGMM sweep, constrained BIC selection, bootstrap stability, MCD
├── gates.py         # Hotelling-F T², cross-fitted Q, posterior/margin → status
├── router.py        # fit pipeline, save/load NPZ+JSON, route_sample()
├── report/
│   ├── plots.py     # one function per figure (plotly fig + matplotlib PNG)
│   ├── build.py     # assemble HTML from jinja2 template
│   └── template.html.j2
├── simulate.py      # synthetic count generator for tests & demos
└── cli.py           # `gcnvr` entry point
tests/router/        # pytest, synthetic data only
```

---

## 3. Algorithms (concrete defaults)

### 3.1 Feature construction

1. **Alignment check.** Every count file must have identical `CONTIG/START/END` to the first one, checked via a hash of the coordinate columns. Files that don't match are rejected with a clear error.
2. **Eligibility.** Keep a feature only if it passes all of these:
   - on an autosome
   - `MAPPABILITY ≥ 0.9`
   - `SEGMENTAL_DUPLICATION_CONTENT == 0`
   - `0.25 ≤ GC ≤ 0.75`
   - not in common-CNV, paralog or PAR masks
   - reference median count ≥ profile minimum
3. **Subsample.** Stratified random selection by chromosome × GC decile, with a fixed seed. Caps: WGS 20,000 bins, WES 20,000 targets, panel all eligible targets. The ordered feature list is saved as `feature_index.tsv` with its SHA-256.
4. **Per-sample normalization.**
   - `x = log2(count / median_eligible_autosomal + 0.01)`.
   - **Residualize** against the reference: `r = x − ref_median_feature`.
   - **Robust scale:** `z = r / max(MAD_feature, floor)`, where the floor is the 25th percentile of all feature MADs.
   - **Winsorize** `z` at ±3. This replaces the ±1.2/±1.5 clip, which didn't work.
5. **GC fingerprint** (interpretation and secondary feature). For each sample, the median residual `r` within each of 20 GC bins over all eligible bins, not just the subsample. It is cheap and almost free of CNV influence.
6. **Cache.** One `.npy` per sample per feature-index hash. Routing a new sample reads one count file. Rebuilding the router reuses the cache.

### 3.2 PCA

- Full SVD of the N × p matrix (N ≤ ~2,000, p ≤ 20k fits in memory as float32).
- **k selection: parallel analysis.** Permute each column independently 20 times and keep the PCs whose eigenvalue exceeds the 95th percentile of the permuted eigenvalues. The result is capped by the profile range (panel 2–5, WES/WGS 2–10).
- **Sign convention:** the largest-|loading| entry of each PC is positive.
- **Loading-locality check:** for each PC, compare its squared-loading mass per chromosome arm with that arm's share of features. If any arm holds more than 3× its expected share, the PC is flagged as likely biological and the report shows it in red. Flagged PCs are excluded from clustering unless overridden.

### 3.3 Clustering

- Sweep K = 1…K_max (profile, default 8), with covariance types `tied` and `diag`, `reg_covar = 1e-3 × median PC variance`, and 10 initializations. The seed is fixed.
- Select by **BIC, subject to** min cluster size ≥ N_min (profile: panel 30, WES 50, WGS 30). Solutions that violate it are discarded.
- **Stability:** 100 bootstrap resamples, refit at the chosen K, and compute the ARI against the full fit. Also build a co-clustering matrix. Require median ARI ≥ 0.8. If it fails, step down to K−1, which may reach K = 1.
- Per cluster: **MinCovDet** location and scatter in PC space. Use these for T²; the GMM covariance is used only for posteriors.
- **Exploratory only:** HDBSCAN noise labels shown in the report, as a hint of drift or a continuum.

### 3.4 Gates

- **Posterior:** status `AMBIGUOUS` if the top posterior < 0.90 or the margin between the top two < 0.20.
- **Hotelling T²** against the assigned cluster's robust centre: threshold `p(n+1)(n−1)/(n(n−p)) · F₀.₉₉(p, n−p)`.
- **Q / SPE** = ‖z − ẑ‖² outside the PC subspace. Threshold = 99th percentile of **10-fold cross-fitted** Q on the training set, because an in-sample Q is optimistically small.
- **IRLS robust projection** (Huber, c = 1.345, ≤ 20 iterations) for every scored sample. The down-weighted feature fraction is reported, and a fraction > 0.15 → `QC_FAIL` (probably a very large CNV or a bad sample).
- **Statuses:** `ROUTED | AMBIGUOUS | OOD (T² or Q) | QC_FAIL`.

### 3.5 Independence checks (report section "Is it biology?")

χ² / Fisher tests with Cramér's V of cluster against sex. Optionally also against genotype at the strongest polymorphic CNV loci, measured as raw depth at UGT2B17 and GSTM1 even though those loci are excluded from the features. A significant association is a red flag.

---

## 4. Visualization

### 4.1 Reference report (`router_report.html` + `png/`)

Each figure has a one-line caption saying what to look for.

| # | Figure | What it answers |
|---|--------|-----------------|
| 1 | **Overview panel**: N samples, k PCs, K clusters, cluster sizes, sex per cluster, stability ARI, warnings | Is this router usable? |
| 2 | **Scree + parallel-analysis curve** | How many real technical axes are there? |
| 3 | **PC scatter matrix** (PC1–PC4). A dropdown recolours by cluster / run / instrument / flowcell / run date / extraction kit / library lot / DNA source / sex / coverage / status. Hover shows sample ID, run and metrics. Cluster ellipses come from MCD | Do clusters coincide with known wet-lab factors? |
| 4 | **Cluster × metadata heatmaps**, one per categorical column, with Cramér's V in the title | Which factor explains each cluster? |
| 5 | **Timeline**: PC1/PC2 and cluster label against run date, one point per sample, runs as vertical bands | Drift vs. discrete switches (kit change date) |
| 6 | **GC-bias fingerprint curves** per sample, coloured by cluster, with cluster medians in bold | Physical interpretation of clusters |
| 7 | **Loadings genome track** per PC (Manhattan-style by chromosome), locality flag highlighted. Includes loading vs. GC scatter | Is a PC technical (GC-driven, genome-wide) or biological (one arm)? |
| 8 | **BIC vs. K** per covariance type, with the N_min-violating region shaded | Why this K |
| 9 | **Bootstrap co-clustering heatmap** + ARI histogram | Are the clusters real or arbitrary splits? |
| 10 | **Sample–sample correlation heatmap** (features), ordered by cluster and then hierarchically | Batch structure at a glance, independent of PCA |
| 11 | **T² vs. Q control chart** (log-log), thresholds as lines, coloured by cluster, training outliers labelled | The OOD gate and how tight it is |
| 12 | **Posterior / margin histograms**, AMBIGUOUS samples listed | How sharp the boundaries are |
| 13 | **Independence table** (cluster vs. sex / polymorphic genotypes) | Is it biology? |
| 14 | *(optional, labelled "exploration only")* UMAP coloured like #3 | Non-linear structure PCA might miss |

### 4.2 Per-sample routing view (`route_report.html`, small)

- The sample as a star on the reference PC1/PC2 (and PC3/PC4) cloud, coloured by cluster.
- The sample's position on the T² vs. Q chart.
- Its posterior bar chart.
- The IRLS weight track along the genome, which shows *where* a big CNV was down-weighted.

These figures are also embeddable in the later per-sample QC report.

### 4.3 Run-level drift (later, after go-live)

T², Q, AMBIGUOUS/OOD rate, and EWMA per sequencing run against time, with control limits. This reuses figure 5 and figure 11 code.

---

## 5. CLI

```
gcnvr build-metadata  --runs-root DIR --grz-glob PATTERN [--lims lims.tsv] [--sex sex.tsv] -o sample_metadata.tsv
gcnvr build-features  --profile {wgs,wes,panel} --annotated A.tsv --masks DIR --counts-list counts.txt -o features/
gcnvr build-router    --features features/ --metadata sample_metadata.tsv --profile P -o router/R-YYYY.MM.N/
gcnvr report          --router router/R-…/ --metadata sample_metadata.tsv -o router/R-…/report/
gcnvr route-sample    --router router/R-…/ --counts sample.tsv -o route.json [--html route_report.html]
gcnvr explore         (= build-features + build-router + report with --no-gates; phase 1 shortcut)
```

`router/` contents (matches parent plan §2.1): `router.json`, `router.npz`, `feature_index.tsv`, `training_samples.tsv` (sample, cluster, posterior, margin, T², Q, CV-Q, down-weighted fraction, status, plus metadata), `report/`.

---

## 6. Phases & acceptance criteria

### Phase 1: exploration (decision gate)
Build: `io_counts`, `metadata`, `features`, `cache`, `pca`, `cluster` (sweep + stability), report figures 1–10 and 13–14, `gcnvr explore`, `simulate.py` with tests.

Run on: historical WGS (`wgs.1k`) and TSHC cohorts, plus WES if available.

**Decision per assay**, written up from the report:
- (a) K = 1 or clusters unstable → no routing for this assay.
- (b) Stable clusters explained by wet-lab metadata → proceed to phase 2.
- (c) Stable clusters *not* explained by metadata → investigate first (sex, polymorphic CNVs, sample quality) before trusting them.

### Phase 2: router + gates
Build: `gates`, IRLS projection, NPZ/JSON serialization, figures 11–12, `build-router`, `report`.

Accept when:
- **Held-out calibration:** in 10-fold CV, the false-OOD rate on in-distribution samples is ≤ 2%.
- **Leave-batch-out:** a held-out known batch (e.g., a known kit change) is flagged OOD at ≥ 90%. Report whether T² or Q caught it.
- **Robustness:** a sample with a simulated 30 Mb deletion (count-level binomial thinning) keeps its cluster assignment, with an elevated but sub-threshold down-weighted fraction.
- **Determinism:** two builds on the same input give byte-identical `router.npz`.

### Phase 3: `route-sample` + per-sample view
Build: `route-sample`, the route report, a `router.smk` (per-sample cached features → build → report) and a `route.smk` stub writing `routes/{sample}.json` + `routing_table.tsv`.

Accept when:
- routing a single new sample takes < 30 s for WGS, from one count file
- `route.json` validates against its schema

### Phase 4 (later): run-level drift charts, registry/release integration.

---

## 7. Tests (synthetic, `tests/router/`)

`simulate.py` generates negative-binomial counts for N samples × p intervals:
- with per-batch GC-bias curves (quadratic in GC, different per batch)
- with library-size variation
- with sex (no effect on autosomes)
- with polymorphic CNV loci (genotypes 0/1/2 at AF 0.3)
- with optional large CNVs and novel-batch samples

| Test | Expectation |
|---|---|
| 3 batches, N = 3 × 50 | K = 3, ARI ≥ 0.9 against truth |
| No batch effect | K = 1 selected |
| Polymorphic loci **not** masked | the test shows the failure (genotype clusters) |
| Polymorphic loci masked | clusters = batches |
| Novel batch with a signature orthogonal to the training PCs | flagged by Q, while T² alone misses ≥ 50% (demonstrates the need for Q) |
| In-distribution held-out samples | false-OOD ≈ 1–2% |
| 30 Mb het deletion sample | same cluster as without the CNV; down-weighted fraction > 0 |
| PC driven by one chromosome (injected) | loading-locality flag raised |
| Router round trip | save → load → identical `route_sample` output |
| Determinism | identical NPZ bytes across two runs |
| Interval mismatch between count files | clear error |

---

## 8. Compute notes

- WGS count TSVs have ~2.9M rows. Read them with pandas `usecols=[3]`, `dtype=float32`, `comment='@'` after one alignment check. That is about 1–2 s per file, and a few minutes for 500 files using a multiprocessing pool. Keep only the eligible rows and fingerprint aggregates in memory.
- The PCA/GMM/bootstrap step for N = 500, p = 20k takes seconds to minutes on one node. A Slurm job with 8 cores and 32 GB is ample.
- Per-sample feature caching makes routine routing and router rebuilds incremental.

---

## 9. Open inputs needed from the lab

1. **LIMS export** with extraction kit, library-prep kit/lot, DNA source, and sequencer per sample. Without it, phase 1 can only relate clusters to run, instrument, flowcell and date.
2. **Historical cohorts per assay:** sample counts and time span for WGS, WES and TSHC, including the dates of known wet-lab changes. Those dates are the ground truth for leave-batch-out tests.
3. **Common-CNV / paralog BED** source and version (proposal: gnomAD-SV v4 sites with AF > 1%, lifted to the interval list's build).
4. Whether a dedicated conda env (`gcnvrouter`) next to the existing `snakemake` env is acceptable on the cluster.

# build_model.sh

## Overview

Executes GATK gCNV model building using the final selected cohort. Determines germline contig ploidy and builds the cohort CNV model using GermlineCNVCaller.

## Usage

```bash
build_model.sh <sample_list> <read_counts_dir> <intervals_file> <reference_fasta> <output_dir>
```

### Parameters

- `sample_list`: TSV file with final selected samples
- `read_counts_dir`: Directory containing read count TSV files
- `intervals_file`: Preprocessed interval list
- `reference_fasta`: Reference genome FASTA
- `output_dir`: Output directory for model files

### Example

```bash
./scripts/build_model.sh \
    20241228/final_samples/final_samples.txt \
    20241228/filtered_by_read_depth/gCNV/03_read_counts \
    20241228/preprocessed.interval_list \
    /mnt/storage/db/references/GRCh38.fasta \
    20241228/model/
```

## Input Requirements

### Sample List Format

Tab-delimited with columns:
```
sample_id    protocol    run    coverage    sex    karyotype_status
```

### Read Count Files

Expected location:
```
${read_counts_dir}/${protocol}/${sample_id}.hg38.tsv
```

Example:
```
read_counts_dir/wgs.1k/G12345.hg38.tsv
```

### Read Count TSV Format

GATK read count format (TSV):
```
CONTIG    START    END      COUNT
chr1      1        1000     42
chr1      1001     2000     38
...
```

## Pipeline Steps

### Step 1: Prepare Input Arguments

Collects read count files:
```bash
for sample in final_samples:
    READ_COUNT_ARGS="${READ_COUNT_ARGS} -I ${read_count_file}"
```

### Step 2: Create Contig Ploidy Priors

Generates prior expectations:
```tsv
CONTIG    PLOIDY_PRIOR_0    PLOIDY_PRIOR_1    PLOIDY_PRIOR_2    PLOIDY_PRIOR_3    PLOIDY_PRIOR_4
chr1      0.01              0.01              0.97              0.01              0.00
chr2      0.01              0.01              0.97              0.01              0.00
...
chrX      0.01              0.01              0.97              0.01              0.00
chrY      0.50              0.49              0.01              0.00              0.00
```

**Interpretation**:
- **Autosomes (chr1-22)**: 97% prior for ploidy=2 (diploid)
- **chrX**: 97% prior for ploidy=2 (accounts for both XX and XY)
- **chrY**: 50% for ploidy=0 (females), 49% for ploidy=1 (males)

### Step 3: Determine Germline Contig Ploidy

Runs GATK DetermineGermlineContigPloidy:
```bash
gatk DetermineGermlineContigPloidy \
    -L intervals.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    -I sample1.tsv -I sample2.tsv ... \
    --contig-ploidy-priors priors.tsv \
    --output ploidy/ \
    --output-prefix ploidy
```

**Output**: Per-sample contig ploidy calls in `ploidy/ploidy-calls/`

### Step 4: Build gCNV Model

Runs GATK GermlineCNVCaller in COHORT mode:
```bash
gatk GermlineCNVCaller \
    --run-mode COHORT \
    -L intervals.interval_list \
    -I sample1.tsv -I sample2.tsv ... \
    --contig-ploidy-calls ploidy/ploidy-calls \
    --annotated-intervals intervals.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    --output cohort-model/ \
    --output-prefix cohort
```

**Output**: Cohort model in `cohort-model/cohort-model/`

## Output Structure

```
output_dir/
├── temp/
│   └── contig_ploidy_priors.tsv
├── ploidy/
│   ├── ploidy-calls/
│   │   ├── sample_1/           # Per-sample ploidy
│   │   ├── sample_2/
│   │   └── ...
│   ├── ploidy-model/           # Learned ploidy model
│   └── ploidy_determination.log
├── cohort-model/
│   ├── cohort-model/           # **FINAL MODEL**
│   │   ├── model_*_0/         # Model shards
│   │   ├── model_*_1/
│   │   └── ...
│   └── cohort-calls/          # Optional: calls on training cohort
├── determine_ploidy.log
└── build_model.log
```

## GATK Parameters

### DetermineGermlineContigPloidy

Key parameters:
- `--interval-merging-rule OVERLAPPING_ONLY`: Preserves interval boundaries
- `--contig-ploidy-priors`: Uses custom priors for human genome
- `--output-prefix`: Prefix for output files

### GermlineCNVCaller

Key parameters:
- `--run-mode COHORT`: Builds model (vs CASE for calling)
- `--interval-merging-rule OVERLAPPING_ONLY`: Consistent with ploidy determination
- `--contig-ploidy-calls`: Uses determined ploidy
- `--annotated-intervals`: Uses preprocessed intervals

## Model Components

### Ploidy Model

Per-contig ploidy distributions:
- Learned from cohort data
- Refined from priors
- Sex-stratified for sex chromosomes

### CNV Model

Statistical model including:
- **Baseline coverage**: Mean coverage per interval
- **Variance model**: Coverage variance structure
- **CNV priors**: Prior probability of CNV events
- **Length distribution**: CNV size distribution

## Performance

### Computational Requirements

**DetermineGermlineContigPloidy**:
- Memory: ~4-8 GB
- Time: ~10-30 minutes for 200 samples
- Parallelization: Per-sample (automatic)

**GermlineCNVCaller COHORT mode**:
- Memory: ~16-32 GB
- Time: ~1-3 hours for 200 samples
- Parallelization: Interval shards (automatic)

### Optimization

For faster execution:
```bash
# Increase memory
gatk --java-options "-Xmx32g" GermlineCNVCaller ...

# Adjust sharding
--interval-merging-rule OVERLAPPING_ONLY \
--number-of-eigensamples 20
```

## Log Files

### determine_ploidy.log

Contains:
- Ploidy determination progress
- Per-sample ploidy calls
- Model convergence information

### build_model.log

Contains:
- Model building progress
- Interval processing status
- Model parameters and convergence

## Error Handling

### Common Errors

**Missing read count files**:
```
WARNING: Read count file not found: {file}
```
- Skips sample
- Continues with available samples

**GATK errors**:
- Logged to respective log files
- Pipeline exits with error code

## Integration

Called by `build_gCNVModel.sh` in Step 4:

```bash
bash "${SCRIPT_DIR}/scripts/build_model.sh" \
    "${FINAL_SAMPLES_DIR}/final_samples.txt" \
    "${GCNV_DIR}" \
    "${BASE_DIR}/preprocessed.interval_list" \
    "${REFERENCE_FASTA}" \
    "${MODEL_DIR}"
```

## Using the Model

### For CNV Calling

Use the built model for CASE mode analysis:
```bash
gatk GermlineCNVCaller \
    --run-mode CASE \
    -I new_sample.tsv \
    --contig-ploidy-calls model/ploidy/ploidy-calls \
    --model model/cohort-model/cohort-model \
    --output new_sample_calls/
```

### Model Validation

Check model quality:
```bash
# Review logs for convergence
grep -i "converged" build_model.log

# Check model files exist
ls model/cohort-model/cohort-model/model_*

# Validate ploidy calls
ls model/ploidy/ploidy-calls/
```

## Customization

### Advanced GATK Parameters

**Memory optimization**:
```bash
--interval-merging-rule OVERLAPPING_ONLY \
--learning-rate 0.05 \
--num-training-epochs 50
```

**Coverage normalization**:
```bash
--enable-bias-factors \
--minimum-contig-length 1000000
```

### Custom Ploidy Priors

Modify for different organisms or models:
```bash
# Example: Higher autosomal aneuploidy rate
chr1    0.01    0.05    0.88    0.05    0.01
```

## Troubleshooting

### Model fails to converge

**Problem**: Model building doesn't converge

**Solutions**:
- Increase training epochs: `--num-training-epochs 100`
- Check sample quality and coverage uniformity
- Review interval list for problematic regions

### Poor model performance

**Problem**: Model performs badly on new samples

**Solutions**:
- Increase cohort size (>200 samples recommended)
- Ensure cohort diversity (multiple batches)
- Check for batch effects in read counts
- Validate ploidy determination

### Memory errors

**Problem**: Out of memory during model building

**Solutions**:
- Increase Java heap: `--java-options "-Xmx64g"`
- Reduce interval count
- Process subsets and combine models

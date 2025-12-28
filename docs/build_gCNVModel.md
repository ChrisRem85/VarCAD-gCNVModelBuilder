# build_gCNVModel.sh - Main Pipeline Script

## Overview

The main orchestrator script that executes the complete gCNV model building pipeline. It coordinates all filtering steps, manages directory structure, and ensures proper execution flow.

## Configuration Variables

### Required Paths
- `BASE_PATH`: Base directory for WGS data (default: `/mnt/storage/genetic_data/WGS`)
- `INTERVALS_FILE`: Preprocessed interval list file
- `REFERENCE_FASTA`: Reference genome FASTA file
- `GATK_PATH`: Path to GATK executable (default: `gatk`)

### Pipeline Parameters
- `MIN_COVERAGE`: Minimum read depth filter (default: 25X)
- `MAX_COVERAGE`: Maximum read depth filter (default: 60X)
- `N_MALES`: Number of male samples to select (default: 100)
- `N_FEMALES`: Number of female samples to select (default: 100)
- `PROTOCOL`: Sequencing protocol name (default: `wgs.1k`)

### Auto-Generated
- `MODEL_VERSION`: Auto-generated date stamp (YYYYMMDD format)
- `BASE_DIR`: Output directory for this model version

## Directory Structure

The pipeline creates the following directory structure:

```
YYYYMMDD/
├── filtered_by_coverage.txt          # Samples passing coverage filter
├── samples_filtered_by_coverage.txt  # Samples with hard-linked read counts
├── preprocessed.interval_list        # Preprocessed intervals for GATK
├── filtered_by_karyotype/
│   └── samples_normal_karyotype.txt  # Samples with normal karyotype
├── filtered_unrelated/
│   ├── somalier_relate.*             # Somalier output files
│   ├── related_pairs.txt             # Identified related sample pairs
│   ├── samples_to_remove.txt         # Samples to remove for relatedness
│   └── samples_unrelated.txt         # Final unrelated samples
├── final_samples/
│   └── final_samples.txt             # Selected 100M + 100F samples
├── filtered_by_read_depth/
│   └── gCNV/
│       └── 03_read_counts/
│           └── wgs.1k/               # Hard-linked read count files
└── model/
    ├── temp/                         # Temporary files
    ├── ploidy/                       # Contig ploidy calls
    └── cohort-model/                 # Final gCNV model
```

## Pipeline Flow

### 0. Initialization
```bash
check_dependencies()  # Validates GATK and Somalier are available
create_directories()  # Creates output directory structure
```

### 1. Coverage Filtering
```bash
filter_by_read_depth()
```
- Calls `filter_by_GRZ_MeanDepthOfCoverage.sh` to parse GRZ QC files
- Filters samples by coverage range (25-60X)
- Hard links read count TSV files from run directories
- Output: `samples_filtered_by_coverage.txt`

### 2. Interval Preprocessing
```bash
prepare_intervals_and_ploidy()
```
- Runs GATK PreprocessIntervals
- Prepares intervals for gCNV analysis
- Output: `preprocessed.interval_list`

### 4. Karyotype Filtering
```bash
filter_by_karyotype()
```
- Calls `filter_karyotype.sh`
- Filters for normal diploid karyotypes
- Output: `filtered_by_karyotype/samples_normal_karyotype.txt`

### 5. Relatedness Filtering
```bash
filter_related_samples()
```
- Calls `filter_related.sh`
- Uses Somalier to identify related samples
- Removes minimal samples to eliminate relatedness
- Output: `filtered_unrelated/samples_unrelated.txt`

### 6. Sample Selection
```bash
select_final_samples()
```
- Calls `select_samples.sh`
- Randomly selects 100 males and 100 females
- Output: `final_samples/final_samples.txt`

### 7. Model Building
```bash
build_gcnv_model()
```
- Calls `build_model.sh`
- Determines contig ploidy
- Builds final gCNV cohort model
- Output: `model/cohort-model/`

## Logging

All operations are logged to `YYYYMMDD/build_log_YYYYMMDD.log` with timestamps.

Log format:
```
[YYYY-MM-DD HH:MM:SS] Message
```

## Error Handling

- Uses `set -euo pipefail` for strict error handling
- Exits on any command failure
- Missing dependencies cause immediate exit
- All errors are logged with timestamps

## Execution

```bash
# Standard execution
./build_gCNVModel.sh

# The pipeline will:
# 1. Create a new dated directory
# 2. Execute all steps sequentially
# 3. Log all operations
# 4. Exit with status 0 on success
```

## Output Files

Key output files:
- **Sample Lists**: Tab-delimited files with columns: `sample_id`, `protocol`, `run`, `coverage`, `sex`, `karyotype_status`
- **Log File**: Complete execution log with timestamps
- **Model Files**: GATK gCNV model ready for use in variant calling

## Customization

To modify pipeline behavior:

1. **Coverage Range**: Adjust `MIN_COVERAGE` and `MAX_COVERAGE`
2. **Sample Counts**: Modify `N_MALES` and `N_FEMALES`
3. **Protocol**: Change `PROTOCOL` variable for different sequencing types
4. **Skip Steps**: Comment out function calls in `main()`

## Dependencies

The pipeline checks for:
- `gatk`: GATK toolkit
- `somalier`: For relatedness analysis

Additional implicit dependencies:
- bash, awk, sed, grep
- Python 3 (for relatedness filtering)
- Standard Unix utilities

# select_samples.sh

## Overview

Randomly selects a balanced cohort of male and female samples from the unrelated sample pool. Ensures representation from both sexes for robust model training.

## Usage

```bash
select_samples.sh <input_samples> <n_males> <n_females> <output_file>
```

### Parameters

- `input_samples`: TSV file with unrelated samples
- `n_males`: Number of male samples to select
- `n_females`: Number of female samples to select
- `output_file`: Output file with selected samples

### Example

```bash
./scripts/select_samples.sh \
    20241228/filtered_unrelated/samples_unrelated.txt \
    100 \
    100 \
    20241228/final_samples/final_samples.txt
```

## Input Format

Tab-delimited file with columns:
```
sample_id    protocol    run    coverage    sex    karyotype_status
```

**Required**: `sex` column must contain "male" or "female"

## Output Format

Same format as input:
```
sample_id    protocol    run    coverage    sex    karyotype_status
```

Contains randomly selected samples (N_MALES males + N_FEMALES females)

## Selection Algorithm

### Step 1: Separate by Sex

Samples are split into two groups based on the `sex` column (column 5):
```bash
males=$(awk '$5=="male"' input)
females=$(awk '$5=="female"' input)
```

### Step 2: Check Availability

Verifies sufficient samples are available:
```bash
n_males_available=$(wc -l males.txt)
n_females_available=$(wc -l females.txt)
```

If insufficient samples:
- Logs warning
- Selects all available samples for that sex

### Step 3: Random Selection

Uses `shuf` for random sampling:
```bash
shuf males.txt | head -n ${N_MALES}
shuf females.txt | head -n ${N_FEMALES}
```

**Properties**:
- Unbiased random selection
- No replacement (each sample selected at most once)
- Different result on each run (use `shuf --random-source` for reproducibility)

### Step 4: Combine Results

Merges male and female selections into single output file.

## Randomization

### Default Behavior

Random seed from system entropy (non-reproducible).

### Reproducible Selection

For reproducible results, modify script:
```bash
# Add random seed
shuf --random-source=<(echo ${RANDOM_SEED}) males.txt | head -n ${N_MALES}
```

Or use a fixed seed file:
```bash
shuf --random-source=/dev/urandom males.txt | head -n ${N_MALES}
```

## Sex Balance Rationale

### Why 100:100?

- **Balanced sex representation** for model training
- **Prevents sex bias** in CNV calling
- **Standard practice** for gCNV models

### Alternative Ratios

Modify `N_MALES` and `N_FEMALES` in main script for different ratios:
```bash
N_MALES=150
N_FEMALES=50   # 3:1 ratio
```

## Error Handling

### Insufficient Samples

**Male shortage**:
```
WARNING: Only N males available, requested M
```
- Selects all N available males
- Continues with available samples

**Female shortage**:
```
WARNING: Only N females available, requested M
```
- Selects all N available females
- Continues with available samples

### Output

```
Available: N males, M females
Selected X samples (A males, B females)
```

## Temporary Files

Uses `mktemp -d` for temporary directory:
```
/tmp/tmp.XXXXXXXXXX/
├── males.txt
├── females.txt
├── selected_males.txt
└── selected_females.txt
```

Automatically cleaned up via `trap`:
```bash
trap "rm -rf ${TEMP_DIR}" EXIT
```

## Integration

Called by `build_gCNVModel.sh` in Step 3:

```bash
bash "${SCRIPT_DIR}/scripts/select_samples.sh" \
    "${UNRELATED_DIR}/samples_unrelated.txt" \
    "${N_MALES}" \
    "${N_FEMALES}" \
    "${FINAL_SAMPLES_DIR}/final_samples.txt"
```

## Customization

### Different Selection Criteria

To select by additional criteria:

**Coverage-weighted selection**:
```bash
# Prefer higher coverage samples
sort -t$'\t' -k4 -n -r males.txt | head -n ${N_MALES}
```

**Run-balanced selection**:
```bash
# Ensure samples from multiple runs
# Group by run, select proportionally
```

### Stratified Selection

Select from different strata:
```bash
# Example: Select 50 from each of two age groups
young_males=$(awk '$7=="young" && $5=="male"' | shuf | head -50)
old_males=$(awk '$7=="old" && $5=="male"' | shuf | head -50)
```

## Reproducibility

### For Reproducible Builds

1. **Set seed in environment**:
```bash
export RANDOM_SEED=12345
```

2. **Modify script**:
```bash
shuf --random-source=<(echo $RANDOM_SEED) input.txt
```

3. **Document seed**:
```bash
echo "Random seed: ${RANDOM_SEED}" >> build_log.txt
```

### For Auditing

Log selected sample IDs:
```bash
cut -f1 final_samples.txt > final_sample_ids.txt
```

## Troubleshooting

### No males/females

**Problem**: Empty sex category

**Solutions**:
- Check sex assignment in previous steps
- Verify karyotype filtering worked correctly
- Review input file: `cut -f5 input.txt | sort | uniq -c`

### Unexpected selection

**Problem**: Same samples selected repeatedly

**Solutions**:
- Check `shuf` is working: `echo -e "1\n2\n3" | shuf`
- Verify input file is not corrupted
- Check for duplicate entries: `sort input.txt | uniq -d`

### Wrong column for sex

**Problem**: Sex column not in position 5

**Solutions**:
- Update awk command: `awk '$CORRECT_COLUMN=="male"'`
- Verify column positions: `head -1 input.txt | tr '\t' '\n' | nl`

# filter_related.sh

## Overview

Uses Somalier to identify and remove related samples from the cohort. Implements an iterative algorithm to remove the minimum number of samples needed to eliminate all relatedness relationships above a threshold.

## Usage

```bash
filter_related.sh <input_samples> <base_path> <output_dir>
```

### Parameters

- `input_samples`: TSV file with karyotype-filtered samples
- `base_path`: Base WGS data directory (contains somalier files)
- `output_dir`: Output directory for results

### Example

```bash
./scripts/filter_related.sh \
    20241228/filtered_by_karyotype/samples_normal_karyotype.txt \
    /mnt/storage/genetic_data/WGS \
    20241228/filtered_unrelated/
```

## Somalier Files

### Expected Location

```
${BASE_PATH}/somalier/${sample_id}.somalier
```

### File Format

Binary Somalier format containing:
- Genotype information at informative sites
- Used for relatedness inference

### Generation

Somalier files should be pre-generated for all samples:
```bash
somalier extract \
    -d output_dir/ \
    --sites sites.vcf.gz \
    -f reference.fasta \
    sample.bam
```

## Relatedness Threshold

**Default**: 0.25

This corresponds to the relatedness coefficient threshold:
- **< 0.25**: Unrelated or distant relatives
- **≥ 0.25**: Second-degree relatives or closer
  - 0.25: Second-degree (grandparent, aunt/uncle, half-sibling)
  - 0.50: First-degree (parent, child, full sibling)

## Algorithm

### Step 1: Collect Somalier Files

```bash
for sample_id in input_samples:
    find ${BASE_PATH}/somalier/${sample_id}.somalier
```

### Step 2: Run Somalier Relate

```bash
somalier relate \
    --infer \
    -o output_prefix \
    file1.somalier file2.somalier ...
```

**Output**: `output_prefix.pairs.tsv` with pairwise relationships

### Step 3: Identify Related Pairs

Parse Somalier output for relationships > threshold:
```bash
awk -v threshold=0.25 'NR>1 && $3>threshold {print $1"\t"$2"\t"$3}'
```

### Step 4: Iterative Removal

**Goal**: Remove minimum samples to break all relationships

**Algorithm**:
1. Count relationships per sample
2. Remove sample with most relationships
3. Update relationship graph
4. Repeat until no relationships remain

**Implementation**: Python script

```python
while related_pairs_exist:
    count_relationships_per_sample()
    remove_sample_with_most_relationships()
    update_remaining_pairs()
```

### Step 5: Create Unrelated Sample List

Filter input samples by removing identified samples.

## Input Format

Tab-delimited with columns:
```
sample_id    protocol    run    coverage    sex    karyotype_status
```

## Output Files

### related_pairs.txt

Related sample pairs above threshold:
```
sample1_id    sample2_id    relatedness_coefficient
G12345       G67890        0.48
G11111       G22222        0.26
```

### samples_to_remove.txt

Samples identified for removal:
```
G12345
G22222
```

### samples_unrelated.txt

Final unrelated sample list (same format as input):
```
sample_id    protocol    run    coverage    sex    karyotype_status
```

### somalier_relate.* files

Somalier output files:
- `.pairs.tsv`: Pairwise relatedness
- `.samples.tsv`: Per-sample QC
- `.html`: Interactive visualization

## Optimization

The iterative removal algorithm minimizes sample loss by:
1. **Greedy approach**: Always removes most connected sample
2. **Graph-based**: Treats relationships as edges, samples as nodes
3. **Iterative**: Recalculates after each removal

This generally produces near-optimal results but is not guaranteed to find the absolute minimum.

## Performance

- **Time**: O(n²) for Somalier relate, O(n·e) for removal where:
  - n = number of samples
  - e = number of relationships
- **Memory**: Loads all relationships into memory
- **Disk**: Stores intermediate Somalier files

## Error Handling

### Warnings

- Missing Somalier file: `WARNING: Somalier file not found for {sample_id}`
- Sample skipped from analysis

### Errors

- No Somalier files found: `ERROR: No somalier files found`, exit 1

## Output Messages

```
Running somalier relate on N samples...
Identifying related samples...
Removed M samples to eliminate relatedness
N unrelated samples remaining
```

## Integration

Called by `build_gCNVModel.sh` in Step 2b:

```bash
bash "${SCRIPT_DIR}/scripts/filter_related.sh" \
    "${KARYOTYPE_DIR}/samples_normal_karyotype.txt" \
    "${BASE_PATH}" \
    "${UNRELATED_DIR}"
```

## Customization

### Adjust Relatedness Threshold

Modify the threshold variable:
```bash
RELATEDNESS_THRESHOLD=0.125  # For stricter filtering
RELATEDNESS_THRESHOLD=0.5    # For looser filtering (only first-degree)
```

### Alternative Removal Strategies

Replace the Python removal algorithm:

**Random removal**:
```bash
# Remove one from each pair randomly
awk '{print $1}' related_pairs.txt | sort -u > samples_to_remove.txt
```

**Manual curation**:
```bash
# Review related_pairs.txt and manually select samples to remove
```

## Troubleshooting

### Somalier fails

**Problem**: Somalier command fails

**Solutions**:
- Check Somalier is installed: `which somalier`
- Verify Somalier files are valid: `somalier --version`
- Check file permissions

### Too many samples removed

**Problem**: Excessive sample removal

**Solutions**:
- Increase relatedness threshold
- Check for systematic issues (e.g., sample duplication)
- Review Somalier QC metrics in `.samples.tsv`

### Python script fails

**Problem**: Python removal script error

**Solutions**:
- Check Python 3 is available
- Verify file paths in heredoc
- Debug with: `python3 -c "import sys; print(sys.version)"`

### No relationships found

**Problem**: All samples appear unrelated

**Solutions**:
- Verify Somalier was run on correct reference
- Check sites VCF used in extraction
- Review Somalier QC in `.samples.tsv`

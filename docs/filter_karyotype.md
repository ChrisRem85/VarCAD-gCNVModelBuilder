# filter_karyotype.sh

## Overview

Filters samples to include only those with normal karyotypes: diploid autosomes (chr1-22 with 2 copies) and standard sex chromosomes (XX or XY).

## Usage

```bash
filter_karyotype.sh <input_samples> <output_samples>
```

### Parameters

- `input_samples`: Input TSV file with coverage-filtered samples
- `output_samples`: Output TSV file with karyotype-filtered samples

### Example

```bash
./scripts/filter_karyotype.sh \
    20241228/samples_filtered_by_coverage.txt \
    20241228/filtered_by_karyotype/samples_normal_karyotype.txt
```

## Input Format

Expects tab-delimited file with columns:
```
sample_id    protocol    run    coverage
```

## Output Format

Tab-delimited file with additional columns:
```
sample_id    protocol    run    coverage    sex         karyotype_status
G12345       WGS         2409...  42.5       female      normal
G67890       WGS         2409...  38.2       male        normal
```

### New Columns

- **sex**: "male" or "female" based on karyotype
- **karyotype_status**: "normal" for samples passing filter

## Karyotype Detection

### Expected Data Source

The script looks for karyotype files at:
```
${KARYOTYPE_DATA}/${sample_id}.karyotype.txt
```

**Note**: `KARYOTYPE_DATA` path must be configured in the script.

### Karyotype Format

Expected content in karyotype files:
```
46,XX    # Normal female
46,XY    # Normal male
```

### Accepted Patterns

- **46,XX**: Normal female karyotype
- **46,XY**: Normal male karyotype

### Rejected Patterns

Any deviation from 46,XX or 46,XY is filtered out:
- Aneuploidies: 45,X or 47,XXY
- Structural abnormalities
- Mosaicism

## Implementation Notes

### Current Status

⚠️ **Configuration Required**: The script requires implementation-specific karyotype data access.

Set the `KARYOTYPE_DATA` variable:
```bash
KARYOTYPE_DATA="/path/to/karyotype/data"
```

### Alternative Implementations

If karyotype files are not available, consider:

1. **X/Y Coverage Ratio**: Infer sex from read depth on sex chromosomes
2. **gCNV Ploidy Analysis**: Use prior gCNV ploidy calls
3. **Clinical Data**: Import from sample metadata
4. **Skip Step**: If population screening isn't needed

## Sex Determination

```bash
if karyotype matches "XX":
    sex = "female"
else if karyotype matches "XY":
    sex = "male"
```

## Error Handling

### Warnings

- Missing karyotype file: `WARNING: Karyotype data not found for {sample_id}, skipping`
- Sample is skipped (not included in output)

### Output

```
Filtered to N samples with normal karyotype
```

## Integration

Called by `build_gCNVModel.sh` in Step 2:

```bash
bash "${SCRIPT_DIR}/scripts/filter_karyotype.sh" \
    "${BASE_DIR}/samples_filtered_by_coverage.txt" \
    "${KARYOTYPE_DIR}/samples_normal_karyotype.txt"
```

## Customization

### Adding Data Sources

To integrate with different karyotype data:

1. **Modify file lookup**:
```bash
karyotype_file="${YOUR_DATA_PATH}/${sample_id}.your_format"
```

2. **Update parsing logic**:
```bash
karyotype=$(awk 'YOUR_PARSING_LOGIC' "${karyotype_file}")
```

3. **Adjust validation**:
```bash
if [[ "${karyotype}" =~ YOUR_PATTERN ]]; then
    # validate
fi
```

### Inferring from Coverage

Example X/Y ratio method:
```bash
# Get X and Y chromosome coverage
x_coverage=$(get_x_coverage "${sample_id}")
y_coverage=$(get_y_coverage "${sample_id}")

# Calculate ratio
ratio=$(echo "scale=2; $y_coverage / $x_coverage" | bc)

# Determine sex (example thresholds)
if (( $(echo "$ratio > 0.1" | bc -l) )); then
    sex="male"     # Y chromosome present
else
    sex="female"   # Y chromosome absent
fi
```

## Troubleshooting

### No samples output

**Problem**: All samples filtered out

**Solutions**:
- Check `KARYOTYPE_DATA` path is correct
- Verify karyotype files exist and are readable
- Check file format matches expected pattern
- Review karyotype validation logic

### Wrong sex assignments

**Problem**: Incorrect sex determination

**Solutions**:
- Validate karyotype file content
- Check regex patterns for XX/XY detection
- Review parsing logic for file format

### Missing karyotype data

**Problem**: Karyotype files not available

**Solutions**:
1. Implement alternative method (coverage-based)
2. Import from clinical/metadata database
3. Run preliminary ploidy analysis with GATK
4. Skip filtering step if not critical

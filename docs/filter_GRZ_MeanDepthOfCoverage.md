# filter_GRZ_MeanDepthOfCoverage.sh

## Overview

Parses GRZ QC files to extract sample coverage information and filters samples by read depth. This is the first filtering step in the pipeline.

## Usage

```bash
filter_GRZ_MeanDepthOfCoverage.sh <base_path> <min_coverage> <max_coverage> <output_file>
```

### Parameters

- `base_path`: Base directory containing sequencer run folders
- `min_coverage`: Minimum mean depth of coverage (e.g., 25)
- `max_coverage`: Maximum mean depth of coverage (e.g., 60)
- `output_file`: Output TSV file path

### Example

```bash
./scripts/filter_GRZ_MeanDepthOfCoverage.sh \
    /mnt/storage/genetic_data/WGS \
    25 \
    60 \
    filtered_samples.txt
```

## Input Data Structure

The script expects GRZ QC files in the following structure:

```
${BASE_PATH}/
├── 230303_A01917_0010_AHJ3KKDMXY/
│   └── quality_control/
│       └── GRZ/
│           ├── G12345.hg38.final.txt
│           └── G67890.hg38.final.txt
├── 240913_A01917_0057_AHN3M5DSXC/
│   └── quality_control/
│       └── GRZ/
│           └── G*.hg38.final.txt
...
```

### Run Directory Pattern

The script searches for runs matching the pattern:
- Format: `YYMMDD_SEQUENCER_RUNNUMBER_FLOWCELL`
- Example: `240913_A01917_0057_AHN3M5DSXC`

### QC File Pattern

- Pattern: `G*.hg38.final.txt`
- Location: `{run}/quality_control/GRZ/`
- Sample ID: Extracted from filename (e.g., `G12345` from `G12345.hg38.final.txt`)

## GRZ QC File Format

The script extracts the `MeanDepthOfCoverage` field from GRZ QC files.

### Extraction Methods

1. **Primary**: `grep "MeanDepthOfCoverage" | awk '{print $NF}'`
2. **Fallback**: Tab-delimited parsing for column-based format

### Example QC File Content

```
SampleID: G12345
MeanDepthOfCoverage: 42.5
...
```

## Output Format

Tab-delimited file with header:

```
sample_id    protocol    run                               coverage
G12345       WGS         240913_A01917_0057_AHN3M5DSXC    42.5
G67890       WGS         240913_A01917_0057_AHN3M5DSXC    38.2
```

### Columns

1. **sample_id**: Sample identifier extracted from filename
2. **protocol**: Always "WGS" (from directory structure)
3. **run**: Full run directory name
4. **coverage**: Mean depth of coverage value

## Filtering Logic

```bash
if coverage >= min_coverage AND coverage <= max_coverage:
    include sample
```

The script uses `awk` for float-safe comparison.

## Error Handling

### Warnings

- Missing coverage value: `WARNING: Could not extract coverage from {file}`
- Continues processing other files

### Errors

- No QC files found: Exits with code 1
- Invalid paths: Bash error from file operations

## Output Messages

```
Searching for GRZ QC files ${BASE_PATH}/{24,25}*/quality_control/GRZ/G*.hg38.final.txt
Found N GRZ QC files
Filtered M samples with coverage between XXX and YYY
```

## Performance Considerations

- Uses file globbing for efficient directory traversal
- Processes files sequentially
- Memory usage: O(1) per file (streams output)
- Time complexity: O(n) where n = number of QC files

## Integration

Called by `build_gCNVModel.sh` in Step 1:

```bash
bash "${SCRIPT_DIR}/scripts/filter_GRZ_MeanDepthOfCoverage.sh" \
    "${BASE_PATH}" \
    "${MIN_COVERAGE}" \
    "${MAX_COVERAGE}" \
    "${BASE_DIR}/samples_filtered_by_coverage.txt"
```

## Troubleshooting

### No files found

**Problem**: `ERROR: No GRZ QC files found`

**Solutions**:
- Verify `BASE_PATH` is correct
- Check directory structure matches expected pattern
- Ensure QC files exist: `ls ${BASE_PATH}/*/quality_control/GRZ/*.final.txt`

### Coverage not extracted

**Problem**: `WARNING: Could not extract coverage`

**Solutions**:
- Check GRZ file format
- Verify `MeanDepthOfCoverage` field exists
- Examine sample file manually: `cat {qc_file} | grep -i coverage`

### Wrong protocol

Currently hardcoded to "WGS". To support other protocols:
- Modify protocol detection logic
- Add protocol parameter or detection from directory structure

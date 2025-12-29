#!/bin/bash

# Script to copy/link read count files for filtered samples
# Usage: copy_read_count_files.sh <filtered_samples_file> <base_path> <protocol> <output_dir>

set -euo pipefail

INPUT_FILE="$1"
BASE_PATH="$2"
PROTOCOL="$3"
OUTPUT_DIR="$4"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo "Copying read count files for input samples..."

# Read header and identify column positions
header=$(head -n1 "${INPUT_FILE}")
sample_col=$(echo "${header}" | awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="sample_id") print i}')
run_col=$(echo "${header}" | awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="run") print i}')

if [ -z "${sample_col}" ] || [ -z "${run_col}" ]; then
    echo "ERROR: Could not find required columns (sample_id, run) in header"
    exit 1
fi

echo "Using columns: sample_id=${sample_col}, run=${run_col}"

input_samples=0
copied_samples=0

tail -n +2 "${INPUT_FILE}" | while IFS=$'\t' read -r -a fields; do
    sample_id="${fields[$((sample_col-1))]}"
    run="${fields[$((run_col-1))]}"
    
    filtered_samples=$((filtered_samples + 1))
    
    # Use run directory to locate read count file directly
    tsv_file="${BASE_PATH}/${run}/gCNV/03_read_counts/${PROTOCOL}/${sample_id}.hg38.tsv"
    
    if [ -f "${tsv_file}" ]; then
        dst_file="${OUTPUT_DIR}/${sample_id}.hg38.tsv"
        if cp -l "${tsv_file}" "${dst_file}" 2>/dev/null || cp "${tsv_file}" "${dst_file}"; then
            copied_samples=$((copied_samples + 1))
        fi
    else
        echo "WARNING: Read count file not found: ${tsv_file}" >&2
    fi
done

echo "Read count files copied: ${copied_samples}/${input_samples} samples"
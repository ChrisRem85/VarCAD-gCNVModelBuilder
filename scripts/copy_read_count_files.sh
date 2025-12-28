#!/bin/bash

# Script to copy/link read count files for filtered samples
# Usage: copy_read_count_files.sh <filtered_samples_file> <base_path> <protocol> <output_dir>

set -euo pipefail

FILTERED_FILE="$1"
BASE_PATH="$2"
PROTOCOL="$3"
OUTPUT_DIR="$4"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo "Copying read count files for filtered samples..."

filtered_samples=0
copied_files=0

while IFS=$'\t' read -r sample_id run coverage; do
    # Skip header
    if [ "${sample_id}" = "sample_id" ]; then
        continue
    fi
    
    filtered_samples=$((filtered_samples + 1))
    
    # Use run directory to locate read count file directly
    tsv_file="${BASE_PATH}/${run}/gCNV/03_read_counts/${PROTOCOL}/${sample_id}.hg38.tsv"
    
    if [ -f "${tsv_file}" ]; then
        dst_file="${OUTPUT_DIR}/${sample_id}.hg38.tsv"
        if cp -l "${tsv_file}" "${dst_file}" 2>/dev/null || cp "${tsv_file}" "${dst_file}"; then
            copied_files=$((copied_files + 1))
        fi
    else
        echo "WARNING: Read count file not found: ${tsv_file}" >&2
    fi
done < "${FILTERED_FILE}"

echo "Read count files copied: ${copied_files}/${filtered_samples} samples"
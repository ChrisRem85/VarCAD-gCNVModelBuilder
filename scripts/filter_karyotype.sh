#!/bin/bash

# Script to filter samples by normal karyotype
# Normal karyotype: chr1-22 diploid (2n), sex chromosomes XX or XY
# Usage: filter_karyotype.sh <input_samples> <output_samples>

set -euo pipefail

INPUT_SAMPLES="$1"
OUTPUT_SAMPLES="$2"

# Path to karyotype data - update this
KARYOTYPE_DATA="/path/to/karyotype/data"

# Create output directory if needed
mkdir -p "$(dirname "${OUTPUT_SAMPLES}")"

# Create header
echo -e "sample_id\tprotocol\trun\tcoverage\tsex\tkaryotype_status" > "${OUTPUT_SAMPLES}"

# Check each sample's karyotype
while IFS=$'\t' read -r sample_id protocol run coverage; do
    # Skip header
    if [ "${sample_id}" = "sample_id" ]; then
        continue
    fi
    
    # Check karyotype - adjust based on your data format
    # This could come from:
    # 1. Ploidy calls from previous GATK runs
    # 2. Clinical data
    # 3. Sex chromosome ratio analysis
    
    karyotype_file="${KARYOTYPE_DATA}/${sample_id}.karyotype.txt"
    
    if [ -f "${karyotype_file}" ]; then
        # Parse karyotype file
        karyotype=$(cat "${karyotype_file}")
        
        # Check for normal karyotype patterns
        if [[ "${karyotype}" =~ ^46,XX$ ]] || [[ "${karyotype}" =~ ^46,XY$ ]]; then
            # Determine sex
            if [[ "${karyotype}" =~ XX ]]; then
                sex="female"
            else
                sex="male"
            fi
            
            echo -e "${sample_id}\t${protocol}\t${run}\t${coverage}\t${sex}\tnormal" >> "${OUTPUT_SAMPLES}"
        fi
    else
        # Alternative: infer from X/Y coverage ratio if karyotype file not available
        # This is a simplified approach - adjust as needed
        
        # For now, we'll use a placeholder that checks sex from filename or metadata
        # You should implement actual karyotype checking based on your data
        
        # Example: Check if sample has normal autosomal ploidy (2) and either XX or XY
        # This would typically involve analyzing the gCNV read counts or other ploidy data
        
        echo "WARNING: Karyotype data not found for ${sample_id}, skipping" >&2
    fi
done < "${INPUT_SAMPLES}"

n_filtered=$(tail -n +2 "${OUTPUT_SAMPLES}" | wc -l)
echo "Filtered to ${n_filtered} samples with normal karyotype"

#!/bin/bash

# Script to get GRZ QC data and filter samples by read depth
# Searches for GRZ QC files matching pattern: {base_path}/25*/quality_control/GRZ/G*.hg38.final.txt
# Usage: get_GRZ_QC.sh <base_path> <min_coverage> <max_coverage> <output_file>

set -euo pipefail

BASE_PATH="$1"
MIN_COVERAGE="$2"
MAX_COVERAGE="$3"
OUTPUT_FILE="$4"

# Create header
echo -e "sample_id\trun\tcoverage" > "${OUTPUT_FILE}"

echo "Searching for GRZ QC files ${BASE_PATH}/{24,25}*/quality_control/GRZ/G*.hg38.final.txt"

# Counter for found files
files_found=0
samples_passed=0

# Search for GRZ QC files matching the pattern (all subdirectories)
for qc_file in "${BASE_PATH}/{24,25}*/quality_control/GRZ/G*.hg38.final.txt"; do
    if [ -f "${qc_file}" ]; then
        files_found=$((files_found + 1))
        
        # Extract full run directory name (e.g., 240913_A01917_0057_AHN3M5DSXC)
        run=$(echo "${qc_file}" | sed -E 's|.*/([^/]+)/quality_control/.*|\1|')
        
        # Extract sample ID from filename (G*.hg38.final.txt)
        filename=$(basename "${qc_file}")
        sample_id=$(echo "${filename}" | sed 's/\.hg38\.final\.txt$//')
        
        # Parse GRZ QC file to extract MeanDepthOfCoverage
        coverage=$(grep "MeanDepthOfCoverage" "${qc_file}" | awk '{print $NF}')
        
        if [ -z "${coverage}" ]; then
            # Fallback: try to extract from tab-delimited format
            coverage=$(awk -F'\t' '/MeanDepthOfCoverage/ {for(i=1;i<=NF;i++) if($i=="MeanDepthOfCoverage") print $(i+1)}' "${qc_file}")
        fi
        
        if [ -n "${coverage}" ]; then
            # Check if coverage is within range using awk for float comparison
            if awk -v cov="${coverage}" -v min="${MIN_COVERAGE}" -v max="${MAX_COVERAGE}" \
                'BEGIN {exit !(cov >= min && cov <= max)}'; then
                echo -e "${sample_id}\t${run}\t${coverage}" >> "${OUTPUT_FILE}"
                samples_passed=$((samples_passed + 1))
            fi
        else
            echo "WARNING: Could not extract coverage from ${qc_file}" >&2
        fi
    fi
done

if [ ${files_found} -eq 0 ]; then
    echo "ERROR: No GRZ QC files found matching pattern ${BASE_PATH}/{24,25}*/quality_control/GRZ/G*.hg38.final.txt"
    exit 1
fi

# Report results
echo "Found ${files_found} GRZ QC files"
echo "Filtered ${samples_passed} samples with coverage between ${MIN_COVERAGE}X and ${MAX_COVERAGE}X"


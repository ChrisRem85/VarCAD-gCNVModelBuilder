#!/bin/bash

# Script to filter related samples using somalier
# Usage: filter_related.sh <input_samples> <base_path> <output_dir>

set -euo pipefail

INPUT_SAMPLES="$1"
BASE_PATH="$2"
OUTPUT_DIR="$3"

RELATEDNESS_THRESHOLD=0.1  # Threshold for considering samples related

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}/somalier/extracted"

# Extract sample IDs
tail -n +2 "${INPUT_SAMPLES}" | cut -f1 > "${OUTPUT_DIR}/sample_ids.txt"

# Collect and copy somalier files for these samples
SOMALIER_FILES=()
while read -r sample_id run coverage sex; do
    # Read the full line to get run information
    line=$(grep "^${sample_id}" "${INPUT_SAMPLES}" | head -n1)
    run=$(echo "${line}" | cut -f2)
    
    somalier_file="${BASE_PATH}/${run}/somalier/extracted/${sample_id}.hg38.somalier"
    if [ -f "${somalier_file}" ]; then
        # Copy to output directory
        cp -l "${somalier_file}" "${OUTPUT_DIR}/somalier/extracted/" 2>/dev/null || cp "${somalier_file}" "${OUTPUT_DIR}/somalier/extracted/"
        SOMALIER_FILES+=("${OUTPUT_DIR}/somalier/extracted/${sample_id}.hg38.somalier")
    else
        echo "WARNING: Somalier file not found for ${sample_id}: ${somalier_file}" >&2
    fi
done < "${OUTPUT_DIR}/sample_ids.txt"

if [ ${#SOMALIER_FILES[@]} -eq 0 ]; then
    echo "ERROR: No somalier files found"
    exit 1
fi

# Run somalier relate
echo "Running somalier relate on ${#SOMALIER_FILES[@]} samples..."
/mnt/storage/groups/heinz/tools/somalier-latest relate \
    --infer \
    -o "${OUTPUT_DIR}/somalier_relate" \
    "${SOMALIER_FILES[@]}"

# Parse results and identify related pairs
echo "Identifying related samples..."
awk -v threshold="${RELATEDNESS_THRESHOLD}" \
    'NR>1 && $3>threshold {print $1"\t"$2"\t"$3}' \
    "${OUTPUT_DIR}/somalier_relate.pairs.tsv" > "${OUTPUT_DIR}/related_pairs.txt"

# Remove minimum number of samples to eliminate all relationships
# Strategy: iteratively remove the sample with most relationships
python3 - <<'EOF' > "${OUTPUT_DIR}/samples_to_remove.txt"
import sys

# Read related pairs
related_pairs = []
with open("${OUTPUT_DIR}/related_pairs.txt", 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            related_pairs.append((parts[0], parts[1]))

# Count relationships per sample
from collections import defaultdict
relationship_counts = defaultdict(int)
for s1, s2 in related_pairs:
    relationship_counts[s1] += 1
    relationship_counts[s2] += 1

# Iteratively remove samples with most relationships
removed = set()
remaining_pairs = set(related_pairs)

while remaining_pairs:
    # Count relationships for non-removed samples
    current_counts = defaultdict(int)
    for s1, s2 in remaining_pairs:
        if s1 not in removed and s2 not in removed:
            current_counts[s1] += 1
            current_counts[s2] += 1
    
    if not current_counts:
        break
    
    # Find sample with most relationships
    max_sample = max(current_counts.items(), key=lambda x: x[1])[0]
    removed.add(max_sample)
    
    # Remove pairs involving this sample
    remaining_pairs = {(s1, s2) for s1, s2 in remaining_pairs 
                      if s1 != max_sample and s2 != max_sample}

# Output samples to remove
for sample in removed:
    print(sample)
EOF

# Create list of unrelated samples
echo -e "sample_id\trun\tcoverage\tsex" > "${OUTPUT_DIR}/samples_unrelated.txt"
while IFS=$'\t' read -r sample_id run coverage sex; do
    if [ "${sample_id}" = "sample_id" ]; then
        continue
    fi
    
    if ! grep -q "^${sample_id}$" "${OUTPUT_DIR}/samples_to_remove.txt" 2>/dev/null; then
        echo -e "${sample_id}\t${run}\t${coverage}\t${sex}" >> "${OUTPUT_DIR}/samples_unrelated.txt"
    fi
done < "${INPUT_SAMPLES}"

n_removed=$(wc -l < "${OUTPUT_DIR}/samples_to_remove.txt" 2>/dev/null || echo 0)
n_unrelated=$(tail -n +2 "${OUTPUT_DIR}/samples_unrelated.txt" | wc -l)

echo "Removed ${n_removed} samples to eliminate relatedness"
echo "${n_unrelated} unrelated samples remaining"

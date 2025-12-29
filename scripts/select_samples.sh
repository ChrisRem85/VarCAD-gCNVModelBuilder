#!/bin/bash

# Script to randomly select N males and N females
# Usage: select_samples.sh <input_samples> <n_males> <n_females> <output_file>

set -euo pipefail

INPUT_SAMPLES="$1"
N_MALES="$2"
N_FEMALES="$3"
OUTPUT_DIR="$4"


# Read header and identify column positions
header=$(head -n1 "${INPUT_SAMPLES}")
sex_col=$(echo "${header}" | awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="sex") print i}')

if [ -z "${sex_col}" ]; then
    echo "ERROR: Could not find 'sex' column in header"
    exit 1
fi

echo "Using column ${sex_col} for sex"

# Separate by sex using the identified column
tail -n +2 "${INPUT_SAMPLES}" | awk -F'\t' -v col="${sex_col}" '$col=="male"' > "${OUTPUT_DIR}/males.txt"
tail -n +2 "${INPUT_SAMPLES}" | awk -F'\t' -v col="${sex_col}" '$col=="female"' > "${OUTPUT_DIR}/females.txt"

n_males_available=$(wc -l < "${OUTPUT_DIR}/males.txt")
n_females_available=$(wc -l < "${OUTPUT_DIR}/females.txt")
echo "Available: ${n_males_available} males, ${n_females_available} females"

# Check if we have enough samples
if [ "${n_males_available}" -lt "${N_MALES}" ]; then
    echo "WARNING: Only ${n_males_available} males available, requested ${N_MALES}"
    N_MALES="${n_males_available}"
fi

if [ "${n_females_available}" -lt "${N_FEMALES}" ]; then
    echo "WARNING: Only ${n_females_available} females available, requested ${N_FEMALES}"
    N_FEMALES="${n_females_available}"
fi

# Randomly select samples
shuf "${OUTPUT_DIR}/males.txt" | head -n "${N_MALES}" > "${OUTPUT_DIR}/selected_males.txt"
shuf "${OUTPUT_DIR}/females.txt" | head -n "${N_FEMALES}" > "${OUTPUT_DIR}/selected_females.txt"

# Combine and create output with same header as input
head -n1 "${INPUT_SAMPLES}" > "${OUTPUT_DIR}/final_samples.txt"
cat "${OUTPUT_DIR}/selected_males.txt" "${OUTPUT_DIR}/selected_females.txt" >> "${OUTPUT_DIR}/final_samples.txt"

n_selected=$(tail -n +2 "${OUTPUT_DIR}/final_samples.txt" | wc -l)
echo "Selected ${n_selected} samples (${N_MALES} males, ${N_FEMALES} females)"


# Copy read count files for selected samples
log "Copying read count files for selected samples..."

mkdir -p "${BASE_DIR}/selected_samples/gCNV/03_read_counts/${PROTOCOL}"
    
bash "${SCRIPT_DIR}/scripts/copy_read_count_files.sh" \
    "${BASE_DIR}/selected_samples/gCNV/03_read_counts/final_samples.txt" \
    "${BASE_PATH}" \
    "${PROTOCOL}" \
    "${BASE_DIR}/selected_samples/gCNV/03_read_counts/${PROTOCOL}"
    
log "Read count files copied to ${BASE_DIR}/selected_samples/gCNV/03_read_counts/${PROTOCOL}/"
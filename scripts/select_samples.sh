#!/bin/bash

# Script to randomly select N males and N females
# Usage: select_samples.sh <input_samples> <n_males> <n_females> <output_file>

set -euo pipefail

INPUT_SAMPLES="$1"
N_MALES="$2"
N_FEMALES="$3"
OUTPUT_FILE="$4"

TEMP_DIR=$(mktemp -d)
trap "rm -rf ${TEMP_DIR}" EXIT

# Separate by sex (sex is in column 5 now with run added)
tail -n +2 "${INPUT_SAMPLES}" | awk '$5=="male"' > "${TEMP_DIR}/males.txt"
tail -n +2 "${INPUT_SAMPLES}" | awk '$5=="female"' > "${TEMP_DIR}/females.txt"

n_males_available=$(wc -l < "${TEMP_DIR}/males.txt")
n_females_available=$(wc -l < "${TEMP_DIR}/females.txt")

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
shuf "${TEMP_DIR}/males.txt" | head -n "${N_MALES}" > "${TEMP_DIR}/selected_males.txt"
shuf "${TEMP_DIR}/females.txt" | head -n "${N_FEMALES}" > "${TEMP_DIR}/selected_females.txt"

# Combine and create output
echo -e "sample_id\tprotocol\trun\tcoverage\tsex\tkaryotype_status" > "${OUTPUT_FILE}"
cat "${TEMP_DIR}/selected_males.txt" "${TEMP_DIR}/selected_females.txt" >> "${OUTPUT_FILE}"

n_selected=$(tail -n +2 "${OUTPUT_FILE}" | wc -l)
echo "Selected ${n_selected} samples (${N_MALES} males, ${N_FEMALES} females)"

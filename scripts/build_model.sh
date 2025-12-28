#!/bin/bash

# Script to build GATK gCNV model
# Usage: build_model.sh <sample_list> <read_counts_dir> <intervals_file> <reference_fasta> <output_dir>

set -euo pipefail

SAMPLE_LIST="$1"
READ_COUNTS_DIR="$2"
INTERVALS_FILE="$3"
REFERENCE_FASTA="$4"
OUTPUT_DIR="$5"

GATK_PATH="gatk"
TEMP_DIR="${OUTPUT_DIR}/temp"
PLOIDY_DIR="${OUTPUT_DIR}/ploidy"
MODEL_OUTPUT="${OUTPUT_DIR}/cohort-model"

mkdir -p "${TEMP_DIR}" "${PLOIDY_DIR}" "${MODEL_OUTPUT}"

# Prepare read count input arguments
READ_COUNT_ARGS=""
while IFS=$'\t' read -r sample_id protocol run coverage sex karyotype; do
    if [ "${sample_id}" = "sample_id" ]; then
        continue
    fi
    
    tsv_file="${READ_COUNTS_DIR}/${protocol}/${sample_id}.hg38.tsv"
    if [ -f "${tsv_file}" ]; then
        READ_COUNT_ARGS="${READ_COUNT_ARGS} -I ${tsv_file}"
    else
        echo "WARNING: Read count file not found: ${tsv_file}" >&2
    fi
done < "${SAMPLE_LIST}"

# Step 1: Determine germline contig ploidy
echo "Step 1: Determining germline contig ploidy..."

# Create contig ploidy priors (standard human: autosomes=2, X/Y as appropriate)
cat > "${TEMP_DIR}/contig_ploidy_priors.tsv" <<EOF
CONTIG	PLOIDY_PRIOR_0	PLOIDY_PRIOR_1	PLOIDY_PRIOR_2	PLOIDY_PRIOR_3	PLOIDY_PRIOR_4
chrX	0.01	0.01	0.97	0.01	0.00
chrY	0.50	0.49	0.01	0.00	0.00
EOF

for chr in {1..22}; do
    echo -e "chr${chr}\t0.01\t0.01\t0.97\t0.01\t0.00" >> "${TEMP_DIR}/contig_ploidy_priors.tsv"
done

${GATK_PATH} DetermineGermlineContigPloidy \
    -L "${INTERVALS_FILE}" \
    --interval-merging-rule OVERLAPPING_ONLY \
    ${READ_COUNT_ARGS} \
    --contig-ploidy-priors "${TEMP_DIR}/contig_ploidy_priors.tsv" \
    --output "${PLOIDY_DIR}" \
    --output-prefix "ploidy" \
    --verbosity DEBUG \
    2>&1 | tee "${OUTPUT_DIR}/determine_ploidy.log"

# Step 2: Build gCNV model
echo "Step 2: Building gCNV cohort model..."

${GATK_PATH} GermlineCNVCaller \
    --run-mode COHORT \
    -L "${INTERVALS_FILE}" \
    ${READ_COUNT_ARGS} \
    --contig-ploidy-calls "${PLOIDY_DIR}/ploidy-calls" \
    --annotated-intervals "${INTERVALS_FILE}" \
    --interval-merging-rule OVERLAPPING_ONLY \
    --output "${MODEL_OUTPUT}" \
    --output-prefix "cohort" \
    --verbosity DEBUG \
    2>&1 | tee "${OUTPUT_DIR}/build_model.log"

echo "gCNV model building completed successfully!"
echo "Ploidy calls: ${PLOIDY_DIR}"
echo "Cohort model: ${MODEL_OUTPUT}"

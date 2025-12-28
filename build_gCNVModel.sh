#!/bin/bash

# GATK gCNV Model Builder
# Main pipeline script for building gCNV models with QC filtering

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODEL_VERSION=$(date +%Y%m%d)
BASE_DIR="${SCRIPT_DIR}/${MODEL_VERSION}"
PROTOCOL="wgs.1k"  # Protocol name (e.g., wgs.1k, wes)
MIN_COVERAGE=25
MAX_COVERAGE=60
N_MALES=100
N_FEMALES=100

# Paths
BASE_PATH="/mnt/storage/genetic_data/WGS"  # Base path for WGS data
INTERVALS_FILE="${SCRIPT_DIR}/assets/hg38.preprocessed.interval_list"
REFERENCE_FASTA="/mnt/storage/db/references/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta"
GATK_VERSION="4.6.0.0"  # GATK version

# Log file
LOG_FILE="${BASE_DIR}/build_log_${MODEL_VERSION}.log"

# Function to log messages
log() {
    # Ensure log directory exists before writing
    if [ ! -d "$(dirname "${LOG_FILE}")" ]; then
        mkdir -p "$(dirname "${LOG_FILE}")"
    fi
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG_FILE}"
}

# Function to check if required tools are available
check_dependencies() {
    log "Checking dependencies..."
    local missing_deps=()
    
    for cmd in gatk somalier; do
        if ! command -v $cmd &> /dev/null; then
            missing_deps+=($cmd)
        fi
    done
    
    if [ ${#missing_deps[@]} -ne 0 ]; then
        log "ERROR: Missing dependencies: ${missing_deps[*]}"
        exit 1
    fi
    
    log "All dependencies found"
}


# Step 1: Get GRZ Mean Depth of Coverage
filter_by_GRZ_MeanDepthOfCoverage() {
    log "STEP 1: Filtering by GRZ Mean Depth of Coverage..."

    mkdir -p "${BASE_DIR}/filtered_by_GRZ_MeanDepthOfCoverage"
    
    bash "${SCRIPT_DIR}/scripts/filter_by_GRZ_MeanDepthOfCoverage.sh" \
        "${BASE_PATH}" \
        "${MIN_COVERAGE}" \
        "${MAX_COVERAGE}" \
        "${BASE_DIR}/filtered_by_GRZ_MeanDepthOfCoverage/filtered_by_GRZ_MeanDepthOfCoverage.txt"
    
    log "Step 1 completed. Coverage data saved to ${BASE_DIR}/filtered_by_GRZ_MeanDepthOfCoverage/filtered_by_GRZ_MeanDepthOfCoverage.txt"
}

# Copy read count files for filtered samples
copy_read_count_files() {
    log "Copying read count files for filtered samples..."

    mkdir -p "${BASE_DIR}/filtered_by_GRZ_MeanDepthOfCoverage/${PROTOCOL}/03_read_counts/"
    
    bash "${SCRIPT_DIR}/scripts/copy_read_count_files.sh" \
        "${BASE_DIR}/filtered_by_GRZ_MeanDepthOfCoverage/filtered_by_GRZ_MeanDepthOfCoverage.txt" \
        "${BASE_PATH}" \
        "${PROTOCOL}" \
        "${BASE_DIR}/filtered_by_GRZ_MeanDepthOfCoverage/${PROTOCOL}/03_read_counts/"
    
    log "Read count files copied to ${BASE_DIR}/filtered_by_GRZ_MeanDepthOfCoverage/${PROTOCOL}/03_read_counts/"
}


gatk_filter_intervals() {
    log "Filtering intervals for gCNV model..."
    
    local threads=128
    local memory=256
    local read_counts_dir="${BASE_DIR}/filtered_by_GRZ_MeanDepthOfCoverage/${PROTOCOL}/03_read_counts"
    local output_dir="${BASE_DIR}/filtered_by_GRZ_MeanDepthOfCoverage/${PROTOCOL}/04_filtered_intervals"
    
    mkdir -p "${output_dir}"
    
    # Build input arguments for all TSV files
    local input_args=""
    for tsv_file in "${read_counts_dir}"/*.hg38.tsv; do
        if [ -f "${tsv_file}" ]; then
            input_args="${input_args} -I ${tsv_file}"
        fi
    done
    
    if [ -z "${input_args}" ]; then
        log "ERROR: No TSV files found in ${read_counts_dir}"
        exit 1
    fi
    
    log "Found $(echo ${input_args} | grep -o ' -I ' | wc -l) read count files"
    
    srun -p all -c ${threads} --mem=${memory}GB \
    docker run --cpus ${threads} -m ${memory}g -u $UID:1002 --rm \
        -v ${BASE_DIR}:${BASE_DIR} \
        -v ${SCRIPT_DIR}:${SCRIPT_DIR}:ro \
        broadinstitute/gatk:${GATK_VERSION} /bin/bash -c " \
            umask 0027; \
            gatk FilterIntervals \
                --java-options '-Xmx${memory}G' \
                --intervals ${SCRIPT_DIR}/assets/${PROTOCOL}/hg38.preprocessed.interval_list \
                --annotated-intervals ${SCRIPT_DIR}/assets/${PROTOCOL}/hg38.annotated.interval_list \
                ${input_args} \
                --interval-merging-rule OVERLAPPING_ONLY \
                -O ${output_dir}/cohort.hg38.filtered.interval_list;" \
    
    log "Intervals filtered and saved to ${output_dir}/cohort.hg38.filtered.interval_list"
}

# # Step 1: Filter intervals and determine contig ploidy
# prepare_intervals_and_ploidy() {
#     log "Filtering intervals and determining contig ploidy..."
    
#     # Filter intervals (example - customize as needed)
#     ${GATK_PATH} PreprocessIntervals \
#         -R "${REFERENCE_FASTA}" \
#         -L "${INTERVALS_FILE}" \
#         --bin-length 0 \
#         --interval-merging-rule OVERLAPPING_ONLY \
#         -O "${BASE_DIR}/preprocessed.interval_list" \
#         2>&1 | tee -a "${LOG_FILE}"
    
#     # Determine contig ploidy - will be done with cohort
#     log "Intervals prepared. Contig ploidy will be determined in cohort analysis"
# }

# Step 2: Filter by karyotype
filter_by_karyotype() {
    log "STEP 2: Filtering samples by normal karyotype..."
    
    bash "${SCRIPT_DIR}/scripts/filter_karyotype.sh" \
        "${BASE_DIR}/samples_filtered_by_coverage.txt" \
        "${KARYOTYPE_DIR}/samples_normal_karyotype.txt"
    
    log "Step 2a completed. $(wc -l < ${KARYOTYPE_DIR}/samples_normal_karyotype.txt) samples with normal karyotype"
}

# Step 2: Filter by relatedness using somalier
filter_related_samples() {
    log "STEP 2b: Filtering related samples using somalier..."
    
    bash "${SCRIPT_DIR}/scripts/filter_related.sh" \
        "${KARYOTYPE_DIR}/samples_normal_karyotype.txt" \
        "${BASE_PATH}" \
        "${UNRELATED_DIR}"
    
    log "Step 2b completed. $(wc -l < ${UNRELATED_DIR}/samples_unrelated.txt) unrelated samples remain"
}

# Step 3: Randomly select final samples
select_final_samples() {
    log "STEP 3: Randomly selecting ${N_MALES} males and ${N_FEMALES} females..."
    
    bash "${SCRIPT_DIR}/scripts/select_samples.sh" \
        "${UNRELATED_DIR}/samples_unrelated.txt" \
        "${N_MALES}" \
        "${N_FEMALES}" \
        "${FINAL_SAMPLES_DIR}/final_samples.txt"
    
    log "Step 3 completed. Selected $(wc -l < ${FINAL_SAMPLES_DIR}/final_samples.txt) samples for model building"
}

# Step 4: Build gCNV model
build_gcnv_model() {
    log "STEP 4: Building gCNV model..."
    
    bash "${SCRIPT_DIR}/scripts/build_model.sh" \
        "${FINAL_SAMPLES_DIR}/final_samples.txt" \
        "${GCNV_DIR}" \
        "${BASE_DIR}/preprocessed.interval_list" \
        "${REFERENCE_FASTA}" \
        "${MODEL_DIR}"
    
    log "Step 4 completed. Model saved to ${MODEL_DIR}"
}

# Main pipeline
main() {
    log "=========================================="
    log "Starting gCNV Model Builder Pipeline"
    log "Version: ${MODEL_VERSION}"
    log "=========================================="
    
    #check_dependencies
    
    # Step 1
    filter_by_GRZ_MeanDepthOfCoverage
    
    copy_read_count_files
    
    gatk_filter_intervals
    
    #filter_by_read_depth
    #prepare_intervals_and_ploidy
    
    # Step 2
    #filter_by_karyotype
    #filter_related_samples
    
    # Step 3
    #select_final_samples
    
    # Step 4
    #build_gcnv_model
    
    log "=========================================="
    log "Pipeline completed successfully!"
    log "Model version: ${MODEL_VERSION}"
    log "Model location: ${MODEL_DIR}"
    log "=========================================="
}

# Run main pipeline
main "$@"

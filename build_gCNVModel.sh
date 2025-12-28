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
GATK_PATH="gatk"  # Update if needed

# Output directories
FILTERED_DIR="${BASE_DIR}/filtered_by_read_depth"
GCNV_DIR="${FILTERED_DIR}/gCNV/03_read_counts"
KARYOTYPE_DIR="${BASE_DIR}/filtered_by_karyotype"
UNRELATED_DIR="${BASE_DIR}/filtered_unrelated"
FINAL_SAMPLES_DIR="${BASE_DIR}/final_samples"
MODEL_DIR="${BASE_DIR}/model"

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

# Create directory structure
create_directories() {
    log "Creating directory structure for version ${MODEL_VERSION}..."
    mkdir -p "${BASE_DIR}"
    #mkdir -p "${GCNV_DIR}"/{wgs.1k}
    #mkdir -p "${KARYOTYPE_DIR}"
    #mkdir -p "${UNRELATED_DIR}"
    #mkdir -p "${FINAL_SAMPLES_DIR}"
    #mkdir -p "${MODEL_DIR}"
    log "Directory structure created"
}

# Step 1: Get GRZ Mean Depth of Coverage
filter_by_GRZ_MeanDepthOfCoverage() {
    log "STEP 1: Filtering by GRZ Mean Depth of Coverage..."
    bash "${SCRIPT_DIR}/scripts/filter_by_GRZ_MeanDepthOfCoverage.sh" \
        "${BASE_PATH}" \
        "${MIN_COVERAGE}" \
        "${MAX_COVERAGE}" \
        "${BASE_DIR}/filtered_by_GRZ_MeanDepthOfCoverage.txt"
    log "Step 1 completed. Coverage data saved to ${BASE_DIR}/filtered_by_GRZ_MeanDepthOfCoverage.txt"
}

# Step 1: Get QC data and filter by read depth
filter_by_read_depth() {
    log "STEP 1: Filtering samples by read depth (${MIN_COVERAGE}X - ${MAX_COVERAGE}X)..."
    
    # Run the QC script
    bash "${SCRIPT_DIR}/scripts/filter_by_GRZ_MeanDepthOfCoverage.sh" \
        "${BASE_PATH}" \
        "${MIN_COVERAGE}" \
        "${MAX_COVERAGE}" \
        "${BASE_DIR}/samples_filtered_by_coverage.txt"
    
    # Hard link the read count files
    log "Hard linking read count files..."
    while IFS=$'\t' read -r sample_id protocol run coverage; do
        # Skip header
        if [ "${sample_id}" = "sample_id" ]; then
            continue
        fi
        
        # Use run directory to locate read count file directly
        tsv_file="${BASE_PATH}/${run}/gCNV/03_read_counts/wgs.1k/${sample_id}.hg38.tsv"
        
        if [ -f "${tsv_file}" ]; then
            dst_file="${GCNV_DIR}/${protocol}/${sample_id}.hg38.tsv"
            ln "${tsv_file}" "${dst_file}" 2>/dev/null || cp "${tsv_file}" "${dst_file}"
        else
            log "WARNING: Read count file not found: ${tsv_file}"
        fi
    done < "${BASE_DIR}/samples_filtered_by_coverage.txt"
    
    log "Step 1 completed. $(wc -l < ${BASE_DIR}/samples_filtered_by_coverage.txt) samples passed coverage filter"
}

# Step 1: Filter intervals and determine contig ploidy
prepare_intervals_and_ploidy() {
    log "Filtering intervals and determining contig ploidy..."
    
    # Filter intervals (example - customize as needed)
    ${GATK_PATH} PreprocessIntervals \
        -R "${REFERENCE_FASTA}" \
        -L "${INTERVALS_FILE}" \
        --bin-length 0 \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O "${BASE_DIR}/preprocessed.interval_list" \
        2>&1 | tee -a "${LOG_FILE}"
    
    # Determine contig ploidy - will be done with cohort
    log "Intervals prepared. Contig ploidy will be determined in cohort analysis"
}

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
    create_directories
    
    # Step 1
    filter_by_GRZ_MeanDepthOfCoverage
    
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

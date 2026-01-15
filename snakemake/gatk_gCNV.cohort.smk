# ************************************************************************************************
# * Snakefile for SNP calling (test version)
# ************************************************************************************************

import os
import math


#Define variables
GDA_VERSION = "3.0"
GATK_VERSION = "4.6.0.0"
GITC_VERSION = "2.3.1-1512499786"
PROTOCOL = config.get("protocol", "wgs.1k")  # Get from config, default to wgs.1k
DB_DIR = "/mnt/storage/db"
FASTA = "references/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta"
INTERVAL_LIST = "regions/wgs_covered.hg38.bed"


#Set FASTA_PATH
FASTA_PATH = DB_DIR + "/" + FASTA


#Set INTERVAL_LIST_PATH
INTERVAL_LIST_PATH = DB_DIR + "/" + INTERVAL_LIST


#Get current working directory
CWD = os.getcwd()

#Get script directory (parent of snakemake directory)
SCRIPT_DIR = os.path.dirname(workflow.basedir)

#Chromosomes
CHROMS = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]


# Create wildcards
DATASETS, = glob_wildcards(CWD + "/gCNV/03_read_counts/" + PROTOCOL + "/{dataset, [A-Za-z0-9\-\_]+}.hg38.tsv")


#Set SAMPLE_INDICES
SAMPLE_INDICES = range(0, len(DATASETS))
#SAMPLE_INDICES = [ "0" ]


SAMPLE_DATASETS_INDICES = {}
SAMPLE_INDICES_DATASETS = {}
for i in range(0, len(DATASETS)):

	SAMPLE_DATASETS_INDICES[sorted(DATASETS)[i]] = str(i)
	SAMPLE_INDICES_DATASETS[str(i)] = sorted(DATASETS)[i]

print(SAMPLE_DATASETS_INDICES)
print(SAMPLE_INDICES_DATASETS)



#Logging
print("Version of Docker image 'genetic_data_analysis': " + GDA_VERSION)
print("Version of Docker image 'gatk': " + GATK_VERSION)
print("Version of Docker image 'gitc': " + GITC_VERSION)
print("Database directory: " + DB_DIR)
print("Path to reference genome: " + FASTA_PATH)
print("Path to interval list: " + INTERVAL_LIST_PATH)
print("Protocol: " + PROTOCOL)
print("Chromosomes: " + ','.join(CHROMS))

print("Current workding directory: " + CWD)
print("Datasets: " + ','.join(DATASETS))
print("Sample indices: " + ','.join(map(str, SAMPLE_INDICES)))


# *** Define Output
OUTPUT = []

OUTPUT.append(CWD + "/gCNV/04_filtered_intervals/" + PROTOCOL + "/cohort.hg38.filtered.interval_list")

OUTPUT.append(CWD + "/gCNV/05_contig_ploidy/" + PROTOCOL + "/cohort.hg38.ploidy-model/ploidy_config.json")
OUTPUT = OUTPUT + expand(CWD + "/gCNV/05_contig_ploidy/" + PROTOCOL + "/cohort.hg38.ploidy-calls/SAMPLE_{sample_index}/contig_ploidy.tsv", sample_index=SAMPLE_INDICES)

OUTPUT = OUTPUT + expand(CWD + "/gCNV/06_scattered_intervals/" + PROTOCOL + "/cohort.hg38.filtered.{chrom}.interval_list", chrom=CHROMS)

OUTPUT = OUTPUT + expand(CWD + "/gCNV/07_raw_cnv/" + PROTOCOL + "/cohort.hg38.{chrom}.cnv-tracking/main_elbo_history.tsv", chrom=CHROMS)
OUTPUT = OUTPUT + expand(CWD + "/gCNV/07_raw_cnv/" + PROTOCOL + "/cohort.hg38.{chrom}.cnv-model/calling_config.json", chrom=CHROMS)
OUTPUT = OUTPUT + expand(CWD + "/gCNV/07_raw_cnv/" + PROTOCOL + "/cohort.hg38.{chrom}.cnv-calls/calling_config.json", chrom=CHROMS)

OUTPUT = OUTPUT + expand(CWD + "/gCNV/08_postprocessed_cnv/" + PROTOCOL + "/cohort.SAMPLE_{sample_index}/denoised_copy_ratios.tsv", sample_index=SAMPLE_INDICES)
OUTPUT = OUTPUT + expand(CWD + "/gCNV/08_postprocessed_cnv/" + PROTOCOL + "/cohort.SAMPLE_{sample_index}/genotyped_intervals.vcf.gz", sample_index=SAMPLE_INDICES)
OUTPUT = OUTPUT + expand(CWD + "/gCNV/08_postprocessed_cnv/" + PROTOCOL + "/cohort.SAMPLE_{sample_index}/genotyped_segments.vcf.gz", sample_index=SAMPLE_INDICES)

# OUTPUT = OUTPUT + expand(CWD + "/gCNV/final/" + PROTOCOL + "/{dataset}.hg38.gCNV." + PROTOCOL + ".cohort.denoised_copy_ratios.bedGraph", dataset=DATASETS)
# OUTPUT = OUTPUT + expand(CWD + "/gCNV/final/" + PROTOCOL + "/{dataset}.hg38.gCNV." + PROTOCOL + ".cohort.genotyped_intervals.vcf.gz", dataset=DATASETS)
OUTPUT = OUTPUT + expand(CWD + "/gCNV/final/" + PROTOCOL + "/{dataset}.hg38.gCNV." + PROTOCOL + ".cohort.genotyped_segments.vcf.gz", dataset=DATASETS)


# ************************************************************************************************

rule all:
    input: OUTPUT

rule test:
    shell:  print(OUTPUT) # print(DATASETS),


# ************************************************************************************************	
# Rules
# ************************************************************************************************		

rule filter_intervals:
    input:   tsv=sorted(expand("{{cwd}}/gCNV/03_read_counts/{{protocol}}/{dataset}.hg38.tsv", dataset=DATASETS)), preprocessed_interval_list=SCRIPT_DIR + "/assets/{protocol}/hg38.preprocessed.interval_list", annotated_interval_list=SCRIPT_DIR + "/assets/{protocol}/hg38.annotated.interval_list"
    output:  interval_list="{cwd}/gCNV/04_filtered_intervals/{protocol}/cohort.hg38.filtered.interval_list"
    params:  input=sorted(expand("-I {{cwd}}/gCNV/03_read_counts/{{protocol}}/{dataset}.hg38.tsv", dataset=DATASETS))
    message: "executing {rule} with output {output} and input {input}"
    threads: 128
    resources:
	    mem_gb=512
    shell:   "umask 0027; \
				mkdir -p $(dirname {output.interval_list}); \
                srun -p all -c {threads} --mem={resources.mem_gb}GB \
				docker run --cpus {threads} -m {resources.mem_gb}g -u $UID:1002 --rm -v {CWD}:{CWD} -v {SCRIPT_DIR}:{SCRIPT_DIR}:ro -v {DB_DIR}:{DB_DIR}:ro broadinstitute/gatk:{GATK_VERSION} /bin/bash -c \" \
					printf 'Container ID:\\t'; hostname; \
					printf 'Start time:\\t'; date; \
					umask 0027; \
					gatk FilterIntervals \
						--java-options '-Xmx{resources.mem_gb}G' \
						--intervals {input.preprocessed_interval_list} \
						--annotated-intervals {input.annotated_interval_list} \
						{params.input} \
						--interval-merging-rule OVERLAPPING_ONLY \
						-O {output.interval_list}; \
					printf 'End time:\\t'; date; \" \
				&> {output.interval_list}.log;"


rule determine_germline_contig_ploidy:
    input:   tsv=sorted(expand("{{cwd}}/gCNV/03_read_counts/{{protocol}}/{dataset}.hg38.tsv", dataset=DATASETS)), interval_list="{cwd}/gCNV/04_filtered_intervals/{protocol}/cohort.hg38.filtered.interval_list", par_bed_gz="/mnt/storage/db/GATK_resources/hg38/par.bed.gz", priors=SCRIPT_DIR + "/assets/contig_ploidy_prior.tsv"
    output:  ploidy_model="{cwd}/gCNV/05_contig_ploidy/{protocol}/cohort.hg38.ploidy-model/ploidy_config.json", ploidy_calls=sorted(expand("{{cwd}}/gCNV/05_contig_ploidy/{{protocol}}/cohort.hg38.ploidy-calls/SAMPLE_{sample_index}/contig_ploidy.tsv", sample_index=SAMPLE_INDICES))
    params:  input=sorted(expand("-I {{cwd}}/gCNV/03_read_counts/{{protocol}}/{dataset}.hg38.tsv", dataset=DATASETS)), ploidy_model_dir="{cwd}/gCNV/05_contig_ploidy/{protocol}/cohort.hg38.ploidy-model", ploidy_calls_dir="{cwd}/gCNV/05_contig_ploidy/{protocol}/cohort.hg38.ploidy-calls"
    message: "executing {rule} with output {output} and input {input}"
    threads: 128
    resources:
	    mem_gb=512
    shell:  "umask 0027; \
			mkdir -p {params.ploidy_model_dir}; \
			mkdir -p {params.ploidy_calls_dir}; \
			rm -r {params.ploidy_model_dir}; \
			rm -r {params.ploidy_calls_dir}; \
			mkdir -p {params.ploidy_model_dir}; \
			mkdir -p {params.ploidy_calls_dir}; \
			srun -p all -c {threads} --mem={resources.mem_gb}GB \
			docker run --cpus {threads} -m {resources.mem_gb}g -u root:1002 --rm -v {CWD}:{CWD} -v {SCRIPT_DIR}:{SCRIPT_DIR}:ro -v {DB_DIR}:{DB_DIR}:ro broadinstitute/gatk:{GATK_VERSION} /bin/bash -c \" \
				printf 'Container ID:\\t'; hostname; \
				printf 'Start time:\\t'; date; \
				umask 0027; \
				gatk DetermineGermlineContigPloidy \
					--java-options '-Xmx{resources.mem_gb}G' \
					--intervals {input.interval_list} \
					--interval-merging-rule OVERLAPPING_ONLY \
					{params.input} \
					--contig-ploidy-priors {input.priors} \
					--output $(dirname {params.ploidy_model_dir}) \
					--output-prefix cohort.hg38.ploidy \
					--verbosity DEBUG; \
				chown -R $UID:1002 {params.ploidy_model_dir}; \
				chown -R $UID:1002 {params.ploidy_calls_dir}; \
				printf 'End time:\\t'; date; \" \
			&> $(dirname {params.ploidy_model_dir})/cohort.hg38.ploidy.log;"


rule scatter_filtered_intervals:
    input:   interval_list="{cwd}/gCNV/04_filtered_intervals/{protocol}/cohort.hg38.filtered.interval_list"
    output:  interval_list="{cwd}/gCNV/06_scattered_intervals/{protocol}/cohort.hg38.filtered.{chrom}.interval_list"
    params:  
    message: "executing {rule} with output {output} and input {input}"
    threads: 2
    resources:
            mem_gb=4
    shell:  "umask 0027; \
			mkdir -p $(dirname {output.interval_list}); \
			srun -p all -c {threads} --mem={resources.mem_gb}GB \
			docker run --cpus {threads} -m {resources.mem_gb}g -u $UID:1002 --rm -v {CWD}:{CWD} -v {SCRIPT_DIR}:{SCRIPT_DIR}:ro -v {DB_DIR}:{DB_DIR}:ro broadinstitute/gatk:{GATK_VERSION} /bin/bash -c \" \
				printf 'Container ID:\\t'; hostname; \
				printf 'Start time:\\t'; date; \
				umask 0027; \
				cat {input.interval_list} | \
				grep '^@SQ' | \
				cat > {output.interval_list}; \
				cat {input.interval_list} | \
				grep -w '^{wildcards.chrom}' | \
				cat >> {output.interval_list}; \
				printf 'End time:\\t'; date; \" \
			&> {output.interval_list}.log;"


rule run_GermlineCNVCaller:
    input:   tsv=sorted(expand("{{cwd}}/gCNV/03_read_counts/{{protocol}}/{dataset}.hg38.tsv", dataset=DATASETS)), ploidy_calls=sorted(expand("{{cwd}}/gCNV/05_contig_ploidy/{{protocol}}/cohort.hg38.ploidy-calls/SAMPLE_{sample_index}/contig_ploidy.tsv", sample_index=SAMPLE_INDICES)), scattered_interval_list="{cwd}/gCNV/06_scattered_intervals/{protocol}/cohort.hg38.filtered.{chrom}.interval_list", annotated_interval_list=SCRIPT_DIR + "/assets/{protocol}/hg38.annotated.interval_list"
    output:  cnv_tracking="{cwd}/gCNV/07_raw_cnv/{protocol}/cohort.hg38.{chrom}.cnv-tracking/main_elbo_history.tsv", cnv_model="{cwd}/gCNV/07_raw_cnv/{protocol}/cohort.hg38.{chrom}.cnv-model/calling_config.json", cnv_calls="{cwd}/gCNV/07_raw_cnv/{protocol}/cohort.hg38.{chrom}.cnv-calls/calling_config.json"
    params:  input=sorted(expand("-I {{cwd}}/gCNV/03_read_counts/{{protocol}}/{dataset}.hg38.tsv", dataset=DATASETS)), ploidy_calls_dir="{cwd}/gCNV/05_contig_ploidy/{protocol}/cohort.hg38.ploidy-calls", cnv_tracking_dir="{cwd}/gCNV/07_raw_cnv/{protocol}/cohort.hg38.{chrom}.cnv-tracking", cnv_model_dir="{cwd}/gCNV/07_raw_cnv/{protocol}/cohort.hg38.{chrom}.cnv-model", cnv_calls_dir="{cwd}/gCNV/07_raw_cnv/{protocol}/cohort.hg38.{chrom}.cnv-calls"
    message: "executing {rule} with output {output} and input {input}"
    threads: 64
    resources:
	    mem_gb=lambda wildcards: 320 if wildcards.chrom in ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chrX"] else 160
    shell:  "umask 0027; \
			mkdir -p {params.cnv_tracking_dir}; \
			mkdir -p {params.cnv_model_dir}; \
			mkdir -p {params.cnv_calls_dir}; \
			rm -rf {params.cnv_tracking_dir}; \
			rm -rf {params.cnv_model_dir}; \
			rm -rf {params.cnv_calls_dir}; \
			mkdir -p {params.cnv_tracking_dir}; \
			mkdir -p {params.cnv_model_dir}; \
			mkdir -p {params.cnv_calls_dir}; \
			srun -p all -c {threads} --mem={resources.mem_gb}GB \
			docker run --cpus {threads} -m {resources.mem_gb}g -u root:1002 --rm -v {CWD}:{CWD} -v {SCRIPT_DIR}:{SCRIPT_DIR}:ro -v {DB_DIR}:{DB_DIR}:ro broadinstitute/gatk:{GATK_VERSION} /bin/bash -c \" \
				printf 'Container ID:\\t'; hostname; \
				printf 'Start time:\\t'; date; \
				umask 0027; \
				export MKL_NUM_THREADS={threads}; \
				export OMP_NUM_THREADS={threads}; \
				gatk GermlineCNVCaller \
					--java-options '-Xmx{resources.mem_gb}G' \
					--run-mode COHORT \
					--intervals {input.scattered_interval_list} \
					{params.input} \
					--contig-ploidy-calls {params.ploidy_calls_dir} \
					--annotated-intervals {input.annotated_interval_list} \
					--interval-merging-rule OVERLAPPING_ONLY \
					--class-coherence-length 100000.0 \
					--cnv-coherence-length 100000.0 \
						--enable-bias-factors false \
						--interval-psi-scale 0.001 \
						--log-mean-bias-standard-deviation 0.1 \
						--sample-psi-scale 0.0001 \
					--depth-correction-tau	10000.0 \
					--p-active 0.01 \
					--p-alt 0.000001 \
					--output $(dirname {params.cnv_tracking_dir}) \
					--output-prefix cohort.hg38.{wildcards.chrom}.cnv \
					--verbosity DEBUG \
					--max-copy-number 5; \
				[[ \$(cat {output.cnv_tracking} | wc -l) -lt 1 ]] && exit 101 || echo 'File size: OK'; \
				[[ \$(cat {output.cnv_model} | wc -l) -lt 1 ]] && exit 101 || echo 'File size: OK'; \
				[[ \$(cat {output.cnv_calls} | wc -l) -lt 1 ]] && exit 101 || echo 'File size: OK'; \
				chown -R $UID:1002 {params.cnv_tracking_dir}; \
				chown -R $UID:1002 {params.cnv_model_dir}; \
				chown -R $UID:1002 {params.cnv_calls_dir}; \
				printf 'End time:\\t'; date; \" \
			&> $(dirname {params.cnv_tracking_dir})/cohort.hg38.{wildcards.chrom}.cnv.log;"


#--disable-caller true \
#--enable-bias-factors false #For WES data


rule postprocess_GermlineCNVCalls:
    input:   ploidy_calls=sorted(expand("{{cwd}}/gCNV/05_contig_ploidy/{{protocol}}/cohort.hg38.ploidy-calls/SAMPLE_{sample_index}/contig_ploidy.tsv", sample_index=SAMPLE_INDICES)), cnv_model=sorted(expand("{{cwd}}/gCNV/07_raw_cnv/{{protocol}}/cohort.hg38.{chrom}.cnv-model/calling_config.json", chrom=CHROMS)), cnv_calls=sorted(expand("{{cwd}}/gCNV/07_raw_cnv/{{protocol}}/cohort.hg38.{chrom}.cnv-calls/calling_config.json", chrom=CHROMS)), fasta=FASTA_PATH
    output:  denoised_copy_ratios_tsv="{cwd}/gCNV/08_postprocessed_cnv/{protocol}/cohort.SAMPLE_{sample_index}/denoised_copy_ratios.tsv", genotyped_intervals_vcf_gz="{cwd}/gCNV/08_postprocessed_cnv/{protocol}/cohort.SAMPLE_{sample_index}/genotyped_intervals.vcf.gz", genotyped_segments_vcf_gz="{cwd}/gCNV/08_postprocessed_cnv/{protocol}/cohort.SAMPLE_{sample_index}/genotyped_segments.vcf.gz"
    params:  ploidy_calls_dir="{cwd}/gCNV/05_contig_ploidy/{protocol}/cohort.hg38.ploidy-calls", model_shard_path=sorted(expand("--model-shard-path {{cwd}}/gCNV/07_raw_cnv/{{protocol}}/cohort.hg38.{chrom}.cnv-model", chrom=CHROMS)), calls_shard_path=sorted(expand("--calls-shard-path {{cwd}}/gCNV/07_raw_cnv/{{protocol}}/cohort.hg38.{chrom}.cnv-calls", chrom=CHROMS)), postprocessed_cnv_dir="{cwd}/gCNV/08_postprocessed_cnv/{protocol}/cohort.SAMPLE_{sample_index}"
    wildcard_constraints:
             sample_index="\d+"
    message: "executing {rule} with output {output} and input {input}"
    threads: 16
    resources:
	    mem_gb=32
    shell:   "umask 0027; \
			mkdir -p {params.postprocessed_cnv_dir}; \
			srun -p all -c {threads} --mem={resources.mem_gb}GB \
			docker run --cpus {threads} -m {resources.mem_gb}g -u root:1002 --rm -v {CWD}:{CWD} -v {SCRIPT_DIR}:{SCRIPT_DIR}:ro -v {DB_DIR}:{DB_DIR}:ro broadinstitute/gatk:{GATK_VERSION} /bin/bash -c \" \
				printf 'Container ID:\\t'; hostname; \
				printf 'Start time:\\t'; date; \
				umask 0027; \
				export MKL_NUM_THREADS={threads}; \
				export OMP_NUM_THREADS={threads}; \
				gatk PostprocessGermlineCNVCalls \
					{params.model_shard_path} \
					{params.calls_shard_path} \
					--contig-ploidy-calls {params.ploidy_calls_dir} \
					--allosomal-contig chrX \
					--allosomal-contig chrY \
					--sample-index {wildcards.sample_index} \
					--output-denoised-copy-ratios {output.denoised_copy_ratios_tsv} \
					--output-genotyped-intervals {output.genotyped_intervals_vcf_gz} \
					--output-genotyped-segments {output.genotyped_segments_vcf_gz} \
					--reference {input.fasta} \
					--verbosity DEBUG; \
				[[ \$(cat {output.denoised_copy_ratios_tsv} | wc -l) -lt 1 ]] && exit 101 || echo 'File size: OK'; \
				[[ \$(bcftools view -H {output.genotyped_intervals_vcf_gz} | wc -l) -lt 1 ]] && exit 101 || echo 'File size: OK'; \
				[[ \$(bcftools view -H {output.genotyped_segments_vcf_gz} | wc -l) -lt 1 ]] && exit 101 || echo 'File size: OK'; \
				chown -R $UID:1002 {params.postprocessed_cnv_dir}; \
				printf 'End time:\\t'; date; \" \
			&> {params.postprocessed_cnv_dir}.log;"


# def get_denoised_copy_ratios(wildcards):
# 	return "{cwd}/gCNV/08_postprocessed_cnv/{protocol}/cohort.SAMPLE_" + SAMPLE_DATASETS_INDICES[wildcards.dataset] + "/denoised_copy_ratios.tsv"


# rule finalize_denoised_copy_ratios:
#     input:   denoised_copy_ratios_tsv=get_denoised_copy_ratios
#     output:  denoised_copy_ratios_bedgraph="{cwd}/gCNV/final/{protocol}/{dataset}.hg38.gCNV.{protocol}.cohort.denoised_copy_ratios.bedGraph"
#     params:  
#     wildcard_constraints:
#              sample_index="\d+"
#     message: "executing {rule} with output {output} and input {input}"
#     threads: 2
#     resources:
#             mem_gb=4
#     shell:   "umask 0027; \
#                 mkdir -p $(dirname {output.denoised_copy_ratios_bedgraph}); \
#                 srun -p all -c {threads} --mem={resources.mem_gb}GB \
#                 docker run --cpus {threads} -m {resources.mem_gb}g -u $UID:1002 --rm -v {CWD}:{CWD} -v {SCRIPT_DIR}:{SCRIPT_DIR}:ro -v {DB_DIR}:{DB_DIR}:ro storage-node:5000/own/genetic_data_analysis:{GDA_VERSION} /bin/bash -c " \
#                         printf 'Container ID:\\t'; hostname; \
#                         printf 'Start time:\\t'; date; \
#                         umask 0027; \
# 			echo 'track type=bedGraph name=\\\"{wildcards.dataset} CN\\\" description=\\\"BedGraph format\\\" visibility=full graphType=points viewLimits=0:5 yLineMark=2 yLineOnOff=on windowingFunction=median color=200,100,0 altColor=0,100,200 priority=20' > {output.denoised_copy_ratios_bedgraph}; \
# 			grep '^chr' {input.denoised_copy_ratios_tsv} | \
# 			cat >> {output.denoised_copy_ratios_bedgraph}; \
# 			[[ \$(cat {output.denoised_copy_ratios_bedgraph} | wc -l) -lt 2 ]] && exit 101 || echo 'File size: OK'; \
#                         printf 'End time:\\t'; date; \" \
#                 &> {output.denoised_copy_ratios_bedgraph}.log;"


# def get_genotyped_intervals(wildcards):
# 	return "{cwd}/gCNV/08_postprocessed_cnv/{protocol}/cohort.SAMPLE_" + SAMPLE_DATASETS_INDICES[wildcards.dataset] + "/genotyped_intervals.vcf.gz"


# rule finalize_genotyped_intervals:
#     input:   genotyped_intervals_vcf_gz=get_genotyped_intervals
#     output:  genotyped_intervals_vcf_gz="{cwd}/gCNV/final/{protocol}/{dataset}.hg38.gCNV.{protocol}.cohort.genotyped_intervals.vcf.gz"
#     params:  
#     wildcard_constraints:
# 	    sample_index="\d+"
#     message: "executing {rule} with output {output} and input {input}"
#     threads: 2
#     resources:
# 	    mem_gb=4
#     shell:   "umask 0027; \
# 		mkdir -p $(dirname {output.genotyped_intervals_vcf_gz}); \
# 		srun -p all -c {threads} --mem={resources.mem_gb}GB \
# 		docker run --cpus {threads} -m {resources.mem_gb}g -u $UID:1002 --rm -v {CWD}:{CWD} -v {SCRIPT_DIR}:{SCRIPT_DIR}:ro -v {DB_DIR}:{DB_DIR}:ro storage-node:5000/own/genetic_data_analysis:{GDA_VERSION} /bin/bash -c " \
# 			printf 'Container ID:\\t'; hostname; \
# 			printf 'Start time:\\t'; date; \
# 			umask 0027; \
# 			cp {input.genotyped_intervals_vcf_gz} {output.genotyped_intervals_vcf_gz}; \
# 			[[ \$(bcftools view -H {output.genotyped_intervals_vcf_gz} | wc -l) -lt 1 ]] && exit 101 || echo 'File size: OK'; \
# 			tabix {output.genotyped_intervals_vcf_gz}; \
# 			printf 'End time:\\t'; date; \" \
# 		&> {output.genotyped_intervals_vcf_gz}.log;"


def get_genotyped_segments(wildcards):
	return "{cwd}/gCNV/08_postprocessed_cnv/{protocol}/cohort.SAMPLE_" + SAMPLE_DATASETS_INDICES[wildcards.dataset] + "/genotyped_segments.vcf.gz"


rule finalize_genotyped_segments:
    input:   genotyped_segments_vcf_gz=get_genotyped_segments
    output:  genotyped_segments_vcf_gz="{cwd}/gCNV/final/{protocol}/{dataset}.hg38.gCNV.{protocol}.cohort.genotyped_segments.vcf.gz"
    params:
    wildcard_constraints:
            sample_index="\d+"
    message: "executing {rule} with output {output} and input {input}"
    threads: 2
    resources:
            mem_gb=4
    shell:   "umask 0027; \
			mkdir -p $(dirname {output.genotyped_segments_vcf_gz}); \
			srun -p all -c {threads} --mem={resources.mem_gb}GB \
			docker run --cpus {threads} -m {resources.mem_gb}g -u $UID:1002 --rm -v {CWD}:{CWD} -v {SCRIPT_DIR}:{SCRIPT_DIR}:ro -v {DB_DIR}:{DB_DIR}:ro storage-node:5000/own/genetic_data_analysis:{GDA_VERSION} /bin/bash -c \" \
				printf 'Container ID:\\t'; hostname; \
				printf 'Start time:\\t'; date; \
				umask 0027; \
				bcftools view \
					--threads {threads} \
					-i 'GT!=\\\"ref\\\" && ALT!=\\\".\\\"' \
					{input.genotyped_segments_vcf_gz} | \
				sed 's&\\./\\.:&0/1:&g' | \
				sed 's&\\.:&1:&g' | \
				bcftools view \
					--threads {threads} \
					-O z \
					-o {output.genotyped_segments_vcf_gz}; \
				[[ \$(bcftools view -H {output.genotyped_segments_vcf_gz} | wc -l) -lt 1 ]] && exit 101 || echo 'File size: OK'; \
				tabix {output.genotyped_segments_vcf_gz}; \
				printf 'End time:\\t'; date; \" \
			&> {output.genotyped_segments_vcf_gz}.log;"


# ************************************************************************************************

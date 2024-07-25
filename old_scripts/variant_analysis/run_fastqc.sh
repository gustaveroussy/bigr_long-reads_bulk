#!/bin/bash
################################################################
## Various Tools Testing
## Using : sbatch run_fastqc.sh
################################################################

#SBATCH --job-name=3582_65NG_dorado_filtered_fastqc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=4
#SBATCH --mem=20G
#SBATCH --partition=mediumq

# Input information
WDIR="/mnt/beegfs/scratch/bioinfo_core/B23043_NADR_02"

#SAMPLE_NAME="3683_CD14plus_ADNg_01082022"
#SAMPLE_NAME="3700_CD14plus_ADNg_25072022"
#SAMPLE_NAME="3700_CD14plus_ADNg_28072022"
#SAMPLE_NAME="3700_R10"

SAMPLE_NAME="3582_R10_65NG"
#SAMPLE_NAME="3582_R10_250NG"
#SAMPLE_NAME="3582_R10_1000NG"
#SAMPLE_NAME="3666_ultra_long_reads"
#SAMPLE_NAME="BC35_control"
#SAMPLE_NAME="2572_CD14"

#FASTQ="${WDIR}/data_input/${SAMPLE_NAME}/${SAMPLE_NAME}/20220801_1123_6C_PAM59295_dd30a530/test/${SAMPLE_NAME}.fastq.gz"
#FASTQ="${WDIR}/data_input/${SAMPLE_NAME}/${SAMPLE_NAME}/20220725_1340_6A_PAM62320_0eacc38a/concat_fastq/${SAMPLE_NAME}.fastq.gz"
#FASTQ="${WDIR}/data_input/${SAMPLE_NAME}/${SAMPLE_NAME}/20220728_1123_6B_PAM62320_b6f6a5e4/concat_fastq/${SAMPLE_NAME}.fastq.gz"
#FASTQ="${WDIR}/data_input/${SAMPLE_NAME}/dorado_rebasecalling/fastq/${SAMPLE_NAME}_dorado.fastq.gz"
#FASTQ="${WDIR}/data_output/${SAMPLE_NAME}/seqkit/3582_R10_65NG_dorado_extracted_from_guppy_pass_readsID.fastq.gz"
FASTQ="${WDIR}/data_output/${SAMPLE_NAME}_dorado/chopper/${SAMPLE_NAME}_merged_noadapter_minlen200_q10.fastq.gz"
#FASTQ="${WDIR}/data_output/${SAMPLE_NAME}/chopper/${SAMPLE_NAME}_filtered_reads.fastq.gz"

THREADS_FASTQC="4"
#4 threads for 71G FastQC and 2 for 50G FastQC

# fastQC - FASTQ QC

mkdir -p ${WDIR}/data_output/${SAMPLE_NAME}/fastqc/

# module load fastqc 
# with the old fastqc version I obtain this error : Exception in thread "Thread-1" java.lang.OutOfMemoryError: Java heap space
# They increased the memory allocation by default in the next versions from FastQC v0.12.1

module load java/17.0.4.1
/mnt/beegfs/userdata/n_rabearivelo/FastQC/fastqc --threads ${THREADS_FASTQC} --memory 4096 --outdir ${WDIR}/data_output/${SAMPLE_NAME}/fastqc/ ${FASTQ}

# default memory option: 512MB - not enough for 50G fastq files --> --memory 1024
# for 71G --> --memory 4096
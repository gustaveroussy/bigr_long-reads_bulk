#!/bin/bash
################################################################
## Various Tools Testing
## Using : sbatch run_pycoqc.sh
################################################################

#SBATCH --job-name=3700_R9_pycoqc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=2
#SBATCH --mem=16G
#SBATCH --partition=shortq

# Input information
WDIR="/mnt/beegfs/scratch/bioinfo_core/B23043_NADR_02"

#SAMPLE_NAME="3683_CD14plus_ADNg_01082022"
#SAMPLE_NAME="3700_R10"
SAMPLE_NAME="3700_CD14plus_ADNg_25072022"

BAM="${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_hg38.sorted.bam"

#sequencing_summary="/mnt/beegfs/scratch/bioinfo_core/B23043_NADR_02/data_input/3700_R10/sequencing_summary_PAQ48746_1980738a_73ad1478.txt"
sequencing_summary="/mnt/beegfs/scratch/bioinfo_core/B23043_NADR_02/data_input/3700_CD14plus_ADNg_25072022/3700_CD14plus_ADNg_25072022/20220725_1340_6A_PAM62320_0eacc38a/sequencing_summary_PAM62320_5af8d2bb.txt"

#FASTQ="${WDIR}/data_input/${SAMPLE_NAME}/${SAMPLE_NAME}/20220801_1123_6C_PAM59295_dd30a530/test/${SAMPLE_NAME}.fastq.gz"
#FASTQ="${WDIR}/data_input/${SAMPLE_NAME}/${SAMPLE_NAME}/20220725_1340_6A_PAM62320_0eacc38a/concat_fastq/${SAMPLE_NAME}.fastq.gz"
#FASTQ="${WDIR}/data_input/${SAMPLE_NAME}/${SAMPLE_NAME}/20220728_1123_6B_PAM62320_b6f6a5e4/concat_fastq/${SAMPLE_NAME}.fastq.gz"
#FASTQ="${WDIR}/data_input/${SAMPLE_NAME}/${SAMPLE_NAME}.fastq.gz"
#FASTQ=${WDIR}/data_input/${SAMPLE_NAME}_cat/${SAMPLE_NAME}.fastq.gz

source /mnt/beegfs/userdata/n_rabearivelo/mambaforge/etc/profile.d/conda.sh

#mkdir -p ${WDIR}/data_output/${SAMPLE_NAME}/nanoqc/
#conda activate nanoQC
#nanoQC -o ${WDIR}/data_output/${SAMPLE_NAME}/nanoqc/ ${FASTQ}
#conda deactivate

mkdir -p ${WDIR}/data_output/${SAMPLE_NAME}/pycoqc/
conda activate pycoqc
pycoQC -f ${sequencing_summary} -a ${BAM} -o ${WDIR}/data_output/${SAMPLE_NAME}/nanoqc/${SAMPLE_NAME}_pycoqc_report.html
conda deactivate 
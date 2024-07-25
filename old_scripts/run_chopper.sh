#!/bin/bash
################################################################
## Various Tools Testing
## Using : sbatch run_chopper.sh
################################################################

#SBATCH --job-name=3700_R10_chopper
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=8
#SBATCH --mem=20G
#SBATCH --partition=mediumq

# Input information
#WDIR="/mnt/beegfs/scratch/bioinfo_core/B23043_NADR_02"

LAMBDA_PHAGE="/mnt/beegfs/userdata/y_mesloub/ressources/references/lambda_phage/enterobacteria_phage_lambda.fasta"

#SAMPLE_NAME="3683_CD14plus_ADNg_01082022"
#SAMPLE_NAME="3700_CD14plus_ADNg_25072022"
#SAMPLE_NAME="3700_R10"
#SAMPLE_NAME="3582_R10_65NG"

#FASTQ="${WDIR}/data_input/${SAMPLE_NAME}/${SAMPLE_NAME}/20220801_1123_6C_PAM59295_dd30a530/test/${SAMPLE_NAME}.fastq.gz"
#FASTQ="${WDIR}/data_input/${SAMPLE_NAME}/${SAMPLE_NAME}/20220725_1340_6A_PAM62320_0eacc38a/concat_fastq/${SAMPLE_NAME}.fastq.gz"
#FASTQ="${WDIR}/data_input/${SAMPLE_NAME}/${SAMPLE_NAME}/20220728_1123_6B_PAM62320_b6f6a5e4/concat_fastq/${SAMPLE_NAME}.fastq.gz"
#FASTQ="${WDIR}/data_input/${SAMPLE_NAME}/${SAMPLE_NAME}.fastq.gz"

THREADS="8"

# Create output directory
mkdir -p ${WDIR}/data_output/${SAMPLE_NAME}_dorado/chopper/

# Chopper - Rust implementation of NanoLyse+NanoFilt : https://github.com/wdecoster/chopper
source /mnt/beegfs/userdata/n_rabearivelo/mambaforge/etc/profile.d/conda.sh
conda activate chopper

echo "Remove Lambda Phage reads from FASTQ file..."
gunzip -c ${FASTQ} | chopper --threads ${THREADS} --minlength 200 --quality 10 --contam ${LAMBDA_PHAGE} | gzip > ${WDIR}/data_output/${SAMPLE_NAME}_dorado/chopper/${SAMPLE_NAME}_noadapter_minlen200_q10.fastq.gz

conda deactivate
#!/bin/bash
################################################################
## Various Tools Testing
## Using : sbatch porechop_abi.sh
################################################################

#SBATCH --job-name=250NG_chopper_filtering
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=6
#SBATCH --mem=10G
#SBATCH --partition=longq
#SBATCH --output=TRIMMING-%x.%j.out
#SBATCH --error=TRIMMING-%x.%j.err

# Input information
WDIR="/mnt/beegfs/scratch/bioinfo_core/B23043_NADR_02"
SAMPLE_NAME="3582_R10_250NG"
#FASTQ="${WDIR}/data_input/${SAMPLE_NAME}/${SAMPLE_NAME}.fastq.gz"
FASTQ_DIR="${WDIR}/data_input/${SAMPLE_NAME}/dorado_rebasecalling/fastq"

LAMBDA_PHAGE="/mnt/beegfs/userdata/y_mesloub/ressources/references/lambda_phage/enterobacteria_phage_lambda.fasta"

# Create output directory
mkdir -p ${WDIR}/data_output/${SAMPLE_NAME}_dorado/porechop_abi/
mkdir -p ${WDIR}/data_output/${SAMPLE_NAME}_dorado/chopper/

# Use Porechop_ABI & Chopper
source /mnt/beegfs/userdata/n_rabearivelo/mambaforge/etc/profile.d/conda.sh

for FASTQ in ${FASTQ_DIR}/${SAMPLE_NAME}_batch*.fastq.gz; do
    BATCH=$(echo ${FASTQ} | cut -d"_" -f11 | cut -d"." -f1)
    #conda activate porechop_abi
    #porechop_abi -t 8 --discard_database --custom_adapters /mnt/beegfs/userdata/n_rabearivelo/references/porechop_adapters.txt --input ${FASTQ} --output ${WDIR}/data_output/${SAMPLE_NAME}_dorado/porechop_abi/${SAMPLE_NAME}_noadapter_${BATCH}.fastq.gz
    #conda deactivate
    
    FASTQ2="${WDIR}/data_output/${SAMPLE_NAME}_dorado/porechop_abi/${SAMPLE_NAME}_noadapter_${BATCH}.fastq.gz"
    
    conda activate chopper
    gunzip -c ${FASTQ2} | chopper --threads 6 --minlength 200 --quality 10 --contam ${LAMBDA_PHAGE} | gzip > ${WDIR}/data_output/${SAMPLE_NAME}_dorado/chopper/${SAMPLE_NAME}_noadapter_minlen200_q10_${BATCH}.fastq.gz
    conda deactivate
done
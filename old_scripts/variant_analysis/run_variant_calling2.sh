#!/bin/bash
################################################################
## Using : sbatch run_variant_calling2.sh
################################################################
 
#SBATCH --job-name=chr20_variant_calling2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=10
#SBATCH --mem=20G
#SBATCH --partition=shortq

# Input information
WDIR="/mnt/beegfs/scratch/bioinfo_core/B23043_NADR_02"
SAMPLE_NAME="3700_R10"
#REF="/mnt/beegfs/userdata/n_rabearivelo/references/Ensembl/T2T-CHM13v2.0/Homo_sapiens-GCA_009914755.4-unmasked.fa"
REF="${WDIR}/script/Homo_sapiens-GCA_009914755.4-unmasked.fa"
BAM="${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_T2T-CHM13.chr20.bam"

# Set the number of CPUs to use
THREADS="10"

# Set up output directory
OUTPUT_DIR="${WDIR}/data_output/${SAMPLE_NAME}/PEPPER_Margin_DeepVariant/chr20/"
OUTPUT_PREFIX="${SAMPLE_NAME}_T2T-CHM13.chr20"
OUTPUT_VCF="${SAMPLE_NAME}_T2T-CHM13.chr20.vcf.gz"

## Create local directory structure
mkdir -p "${OUTPUT_DIR}"

# Pull the docker images
#singularity pull docker://jmcdani20/hap.py:v0.3.12
#singularity pull docker://kishwars/pepper_deepvariant:r0.8

# Run PEPPER-Margin-DeepVariant
module load singularity
singularity exec --bind ${WDIR} /mnt/beegfs/userdata/n_rabearivelo/containers/pepper_deepvariant_r0.8.sif run_pepper_margin_deepvariant call_variant \
--bam "${BAM}" \
--fasta "${REF}" \
--output_dir "${OUTPUT_DIR}" \
--output_prefix "${OUTPUT_PREFIX}" \
--threads ${THREADS} \
--sample_name ${SAMPLE_NAME} \
--ont_r10_q20

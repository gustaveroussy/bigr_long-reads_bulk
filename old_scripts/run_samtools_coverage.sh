#!/bin/bash
################################################################
## Using : sbatch run_samtools_coverage.sh
################################################################

#SBATCH --job-name=3700_R10_hg38_coverage
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=2
#SBATCH --mem=10G
#SBATCH --partition=shortq

# Input information
#WDIR="/mnt/beegfs/scratch/bioinfo_core/B23043_NADR_02"

#SAMPLE_NAME="3683_CD14plus_ADNg_01082022"
#SAMPLE_NAME="3700_CD14plus_ADNg_25072022"
#SAMPLE_NAME="3700_CD14plus_ADNg_28072022"
#SAMPLE_NAME="3700_R10"
#SAMPLE_NAME="3700_R9"
#SAMPLE_NAME="3582_R10_250NG"
#SAMPLE_NAME="3666_ultra_long_reads"
#SAMPLE_NAME="2572_CD14"
#SAMPLE_NAME="BC35_control"
#SAMPLE_NAME="3582_R10_1000NG"

#BAM="${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_hg38.sorted.bam"
#BAM="${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_T2T-CHM13.sorted.bam"

THREADS="2"

# Create output directory
mkdir -p ${WDIR}/data_output/${SAMPLE_NAME}_dorado/samtools/
#mkdir -p ${WDIR}/data_output/${SAMPLE_NAME}/bedtools/

# Samtools Coverage
module load samtools
samtools coverage --output ${WDIR}/data_output/${SAMPLE_NAME}_dorado/samtools/${SAMPLE_NAME}_${REF}_samtools_coverage.txt ${BAM}
#samtools coverage --output ${WDIR}/data_output/${SAMPLE_NAME}/samtools/${SAMPLE_NAME}_T2T_samtools_coverage.txt ${BAM}

# Bedtools genomecov
#module load bedtools
#bedtools genomecov -ibam ${BAM} 
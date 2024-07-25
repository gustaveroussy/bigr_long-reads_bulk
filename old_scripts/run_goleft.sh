#!/bin/bash
################################################################
## Various Tools Testing
## Using : sbatch run_goleft.sh
################################################################

#SBATCH --job-name=3700_R10_goleft
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --partition=shortq

# Input information
WDIR="/mnt/beegfs/scratch/bioinfo_core/B23043_NADR_02"
REF="hg38"

#SAMPLE_NAME="3683_CD14plus_ADNg_01082022"
#SAMPLE_NAME="3700_CD14plus_ADNg_25072022"
#SAMPLE_NAME="3700_CD14plus_ADNg_28072022"
SAMPLE_NAME="3700_R10"
#SAMPLE_NAME="3582_R10_250NG"
#SAMPLE_NAME="3582_R10_1000NG"
#SAMPLE_NAME="3666_ultra_long_reads"
#SAMPLE_NAME="BC35_control"
#SAMPLE_NAME="2572_CD14"

BAM="${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_${REF}.sorted.bam"

# Create output directory
mkdir -p ${WDIR}/data_output/${SAMPLE_NAME}/goleft/

# Run NanoPlot
source /mnt/beegfs/userdata/n_rabearivelo/mambaforge/etc/profile.d/conda.sh
conda activate goleft

goleft indexcov --directory ${WDIR}/data_output/${SAMPLE_NAME}/goleft/ ${BAM}
goleft covstats --directory ${WDIR}/data_output/${SAMPLE_NAME}/goleft/ ${BAM}

conda deactivate
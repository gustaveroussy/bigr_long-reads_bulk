#!/bin/bash
################################################################
## Various Tools Testing
## Using : sbatch run_mosdepth.sh
################################################################

#SBATCH --job-name=3666_ultra_long_reads_mosdepth
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=4
#SBATCH --mem=20G
#SBATCH --partition=shortq 

# Input information
WDIR="/mnt/beegfs/scratch/bioinfo_core/B23043_NADR_02"

#SAMPLE_NAME="3683_CD14plus_ADNg_01082022"
#SAMPLE_NAME="3700_CD14plus_ADNg_25072022"
#SAMPLE_NAME="3700_CD14plus_ADNg_28072022"
#SAMPLE_NAME="3700_R10"
#SAMPLE_NAME="3700_R9"
#SAMPLE_NAME="3582_R10_250NG"
SAMPLE_NAME="3666_ultra_long_reads"

BAM="${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_hg38.sorted.bam"

THREADS_MOSDEPTH="3"
#4 threads for 71G FastQC and 2 for 50G FastQC

# Create output directory
mkdir -p ${WDIR}/data_output/${SAMPLE_NAME}/mosdepth/

# Mosdepth - BAM depth calculator
source /mnt/beegfs/userdata/n_rabearivelo/mambaforge/etc/profile.d/conda.sh
conda activate mosdepth

## wgs mode recommanded on mosdepth github page
cd ${WDIR}/data_output/${SAMPLE_NAME}/mosdepth/
mosdepth -n --fast-mode --by 500 -t ${THREADS_MOSDEPTH} ${SAMPLE_NAME}_wgs_mode ${BAM} && 
mosdepth -t ${THREADS_MOSDEPTH} ${SAMPLE_NAME} ${BAM}

conda deactivate 
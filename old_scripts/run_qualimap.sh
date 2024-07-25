#!/bin/bash
################################################################
## QUALIMAP : Quality control of alignment sequencing data and its derivatives like feature counts. 
## http://qualimap.conesalab.org/doc_html/command_line.html
## Using : sbatch run_qualimap.sh
################################################################

#SBATCH --job-name=2572_T2T_qualimap
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=8
#SBATCH --mem=50G
#SBATCH --partition=mediumq

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

# Create output directory
mkdir -p ${WDIR}/data_output/${SAMPLE_NAME}_dorado/qualimap/${REF}/

# Run Qualimap
source /mnt/beegfs/userdata/n_rabearivelo/mambaforge/etc/profile.d/conda.sh
conda activate qualimap

qualimap bamqc -bam ${BAM} -outdir ${WDIR}/data_output/${SAMPLE_NAME}_dorado/qualimap/${REF}/ --paint-chromosome-limits -nt 10 --java-mem-size=40G

conda deactivate
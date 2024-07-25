#!/bin/bash
################################################################
## Various Tools Testing
## Using : sbatch run_nanoplot.sh
################################################################

#SBATCH --job-name=3582_1000NG_nanoplot_FASTQ
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=2
#SBATCH --mem=16G
#SBATCH --partition=shortq

# Input information
#WDIR="/mnt/beegfs/scratch/bioinfo_core/B23043_NADR_02"

#SAMPLE_NAME="3683_CD14plus_ADNg_01082022"
#SAMPLE_NAME="3700_CD14plus_ADNg_25072022"
#SAMPLE_NAME="3700_CD14plus_ADNg_28072022"
#SAMPLE_NAME="3700_R10"
#SAMPLE_NAME="3582_R10_250NG"
#SAMPLE_NAME="3582_R10_1000NG"
#SAMPLE_NAME="3666_ultra_long_reads"
#SAMPLE_NAME="BC35_control"
#SAMPLE_NAME="2572_CD14"

#FASTQ="${WDIR}/data_input/${SAMPLE_NAME}/${SAMPLE_NAME}/20220801_1123_6C_PAM59295_dd30a530/test/${SAMPLE_NAME}.fastq.gz"
#FASTQ="${WDIR}/data_input/${SAMPLE_NAME}/${SAMPLE_NAME}/20220725_1340_6A_PAM62320_0eacc38a/concat_fastq/${SAMPLE_NAME}.fastq.gz"
#FASTQ="${WDIR}/data_input/${SAMPLE_NAME}/${SAMPLE_NAME}/20220728_1123_6B_PAM62320_b6f6a5e4/concat_fastq/${SAMPLE_NAME}.fastq.gz"
#FASTQ="${WDIR}/data_input/${SAMPLE_NAME}/${SAMPLE_NAME}.fastq.gz"
#FASTQ="${WDIR}/data_output/${SAMPLE_NAME}/chopper/${SAMPLE_NAME}_filtered_reads.fastq.gz"

#SEQUENCING_SUMMARY="${WDIR}/data_input/${SAMPLE_NAME}/sequencing_summary_PAQ48746_1980738a_73ad1478.txt"
#BAM="${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_hg38.sorted.bam"
#BAM="${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_T2T-CHM13.sorted.bam"

# Create output directory
#mkdir -p ${WDIR}/data_output/${SAMPLE_NAME}/nanoplot/

# Run NanoPlot
source /mnt/beegfs/userdata/n_rabearivelo/mambaforge/etc/profile.d/conda.sh
conda activate nanoplot

NanoPlot --threads 2 --outdir ${WDIR}/data_output/${SAMPLE_NAME}/nanoplot/fastq/ --prefix ${SAMPLE_NAME}_unfiltered_fastq --tsv_stats --info_in_report --N50 --fastq ${FASTQ}
#NanoPlot --threads 2 --outdir ${WDIR}/data_output/${SAMPLE_NAME}/nanoplot/sequencing_summary/ --prefix ${SAMPLE_NAME}_seq_summary --N50 --tsv_stats --info_in_report --summary ${SEQUENCING_SUMMARY}
#NanoPlot --threads 10 --outdir ${WDIR}/data_output/${SAMPLE_NAME}/nanoplot/bam/ --prefix ${SAMPLE_NAME}_sorted_bam --N50 --tsv_stats --info_in_report --bam ${BAM}

conda deactivate
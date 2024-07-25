#!/bin/bash
################################################################
## SAMBAMBA : Convert SAM to BAM and SORT + INDEX BAM file
## Using : sbatch run_sambamba.sh
################################################################

#SBATCH --job-name=sort_index_3700_R10_chr20
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=2
#SBATCH --mem=24G
#SBATCH --partition=longq

# Input information
WDIR="/mnt/beegfs/scratch/bioinfo_core/B23043_NADR_02"
#SAMPLE_NAME="3683_CD14plus_ADNg_01082022"
#SAMPLE_NAME="3683_partial"
#INPUT="${WDIR}/data_input/3683_CD14plus_ADNg_01082022/3683_CD14plus_ADNg_01082022/20220801_1123_6C_PAM59295_dd30a530/test/${SAMPLE_NAME}.fastq.gz"

#SAMPLE_NAME="3700_CD14plus_ADNg_25072022"
#INPUT="${WDIR}/data_input/${SAMPLE_NAME}/${SAMPLE_NAME}/20220725_1340_6A_PAM62320_0eacc38a/concat_fastq/${SAMPLE_NAME}.fastq.gz"

SAMPLE_NAME="3700_R10"

#SAMPLE_NAME="3582_R10_250NG"
#SAMPLE_NAME="3666_ultra_long_reads"
#SAMPLE_NAME="2572_CD14"
#SAMPLE_NAME="BC35_control"
#SAMPLE_NAME="3582_R10_1000NG"

# Load needed modules
module load sambamba

# Create output dir
mkdir -p ${WDIR}/data_output/${SAMPLE_NAME}/sambamba/

# sambamba
# VS HG38
# view (2G, 1 CPU/task)
echo "Sorting..."
#sambamba view -S -f bam ${WDIR}/data_output/${SAMPLE_NAME}/minimap2/${SAMPLE_NAME}_${REF}.sam > ${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_${REF}.bam 

# sort (16G, 1 CPU/task)
#sambamba sort ${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_${REF}.bam
sambamba sort /mnt/beegfs/scratch/bioinfo_core/B23043_NADR_02/data_output/3700_R10/whatshap/chr20/3700_R10_hg38_variants.chr20.haplotagged.bam
echo "Indexing..."
#sambamba index ${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_${REF}.sorted.bam
sambamba index /mnt/beegfs/scratch/bioinfo_core/B23043_NADR_02/data_output/3700_R10/whatshap/chr20/3700_R10_hg38_variants.chr20.haplotagged.sorted.bam

#Run on multiple BAM files / batchs
#FASTQ_DIR="/mnt/beegfs/scratch/n_rabearivelo/long_read_NADR/P34DEV3582_1000NG/fastq_pass"
#FASTQ_LIST=$(ls ${FASTQ_DIR} | sed 's,.fastq.gz,,g')
#for FASTQ in ${FASTQ_LIST}; do
    # view (2G, 1 CPU/task)
    #sambamba view -S -f bam ${WDIR}/data_output/${SAMPLE_NAME}/minimap2/${SAMPLE_NAME}_${FASTQ}_hg38.sam > ${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_${FASTQ}_hg38.bam
    # sort (16G, 1 CPU/task)
    #sambamba sort ${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_${FASTQ}_hg38.bam
    #sambamba index ${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_${FASTQ}_hg38.sorted.bam
#done;

# ---- No Need to re-run ---- #

# Create index for both references
#samtools faidx ${HG38}
#samtools faidx ${T2T}

# --------Was not Used ----------

#samtools sort - TO COMPARE ? Sambamba should be better.
#samtools sort ${WDIR}/data_output/${SAMPLE_NAME}/minimap2/PAM59295_7c170434_hg38.sam -o ${WDIR}/data_output/${SAMPLE_NAME}/samtools/${SAMPLE_NAME}_hg38.bam
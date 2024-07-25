#!/bin/bash
################################################################
##CONCAT ALL BAM, FILTER, SORT, INDEX AND SPLIT BY CHROM
################################################################

#SBATCH --job-name=concat_sort_splitByChr
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=12
#SBATCH --partition=mediumq


#LOAD MODULE
module load samtools

#PROJECT FOLDER
WKDIR="/path/to/project/directory"
DATASET=""
INPUT="${WKDIR}/data_output/${DATASET}"

#Moving in dataset folder
cd ${INPUT}

#Concat all bam
mkdir 3_concat_sort_index
echo "Concatenating all aligned bam files"
samtools cat 2_custom_align/b*.bam > 3_concat_sort_index/all_bam_concat_filtered.bam

#Sort and index
echo "Sorting and indexing concatenated bam"
samtools sort -@ 12 -o 3_concat_sort_index/all_bam_concat_filtered.bam_sorted.bam 3_concat_sort_index/all_bam_concat_filtered.bam
samtools index 3_concat_sort_index/all_bam_concat_filtered.bam_sorted.bam

#Separate by chromosome
mkdir 4_chromosome
echo "Separating concatenated by chromosome and indexing"
for i in {1..22} MT X Y;do echo "Chromosome $i";samtools view -b -@ 12 -o 4_chromosome/chr_${i}_filtered.bam 3_concat_sort_index/all_bam_concat_filtered.bam_sorted.bam ${i};samtools index 4_chromosome/chr_${i}_filtered.bam;done

echo "Finished"

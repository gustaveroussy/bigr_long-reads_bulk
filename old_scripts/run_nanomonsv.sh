#!/bin/bash
################################################################
## Usage : sbatch run_nanomonsv.sh
## https://github.com/friend1ws/nanomonsv
################################################################

#SBATCH --job-name=nanomonSV
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=4
#SBATCH --mem=20G
#SBATCH --partition=mediumq
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

# Input parameters
WDIR="/mnt/beegfs/scratch/bioinfo_core/B24018_OLBE_01"

TUMOR="${WDIR}/data_output/somatic_variant_calling/concat_chr/WM99_CD19_all_chr_concat_filtered_sorted.bam"
NORMAL="${WDIR}/data_output/somatic_variant_calling/concat_chr/WM99_CD3_all_chr_concat_filtered_sorted.bam"
REF="/mnt/beegfs/userdata/n_rabearivelo/references/Ensembl/GRCh38.109/homo_sapiens.GRCh38.109.fasta"

# Create output directory
mkdir -p ${WDIR}/data_output/somatic_variant_calling/SV_calling/

# NanomonSV - Somatic Structural Variant Calling
source /mnt/beegfs/userdata/n_rabearivelo/mambaforge/etc/profile.d/conda.sh
conda activate nanomonsv

# 1- Parses all the supporting reads of putative somatic SVs
#nanomonsv parse --reference_fasta ${REF} ${TUMOR} ${WDIR}/data_output/somatic_variant_calling/SV_calling/WM99_CD19_Tumor
#nanomonsv parse --reference_fasta ${REF} ${NORMAL} ${WDIR}/data_output/somatic_variant_calling/SV_calling/WM99_CD3_Normal

# 2- Gets the SV result from the parsed supporting reads data
#nanomonsv get --control_prefix ${WDIR}/data_output/somatic_variant_calling/SV_calling/WM99_CD3_Normal --control_bam ${NORMAL} --processes 10 --max_memory_minimap2 20 --qv15 --use_racon --single_bnd ${WDIR}/data_output/somatic_variant_calling/SV_calling/WM99_CD19_Tumor ${TUMOR} ${REF}

#misc/post_fileter.py to filter the results

# 3- Classifies the long insertions into several mobile element insertions
nanomonsv insert_classify --genome_id hg38 ${WDIR}/data_output/somatic_variant_calling/SV_calling/WM99_CD19_Tumor.nanomonsv.result.txt ${WDIR}/data_output/somatic_variant_calling/SV_calling/WM99_CD19_Tumor.nanomonsv.insert_classify.txt ${REF}

conda deactivate
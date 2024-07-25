#!/bin/bash
################################################################
## Using : sbatch run_structural_variant_calling.sh
################################################################

#SBATCH --job-name=SV_calling
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=4
#SBATCH --mem=20G
#SBATCH --partition=shortq
#SBATCH --output=SV_calling-%x.%j.out
#SBATCH --error=SV_calling-%x.%j.err
 
# Input information
WDIR="/mnt/beegfs/scratch/bioinfo_core/B23043_NADR_02"
SAMPLE_NAME="3582_R10_250NG"
REF="hg38"
BAM="${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_${REF}.sorted.bam"

#Create output directory
mkdir -p ${WDIR}/data_output/dorado/${SAMPLE_NAME}/sniffles/

#Run Sniffles2
source /mnt/beegfs/userdata/n_rabearivelo/mambaforge/etc/profile.d/conda.sh
conda activate sniffles

sniffles -i ${BAM} -v ${WDIR}/data_output/dorado/${SAMPLE_NAME}/sniffles/${SAMPLE_NAME}_${REF}_SV.vcf

conda deactivate

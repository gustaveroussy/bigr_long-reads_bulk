#!/bin/bash
################################################################
## Using : sbatch run_annotsv.sh
################################################################

#SBATCH --job-name=AnnotSV_generate_vcf
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=6G
#SBATCH --partition=shortq

# Input information
INPUT_PATH="${WDIR}/data_output/${SAMPLE_NAME}/sniffles"
OUTDIR="${WDIR}/data_output/${SAMPLE_NAME}/annotsv_vcf"

mkdir -p ${OUTDIR}

source /mnt/beegfs/userdata/n_rabearivelo/mambaforge/etc/profile.d/conda.sh
conda activate annotsv

AnnotSV -genomeBuild GRCh38 -annotationsDir /mnt/beegfs/userdata/n_rabearivelo/AnnotSV_annotations -SVinputFile ${INPUT_PATH}/${SAMPLE_NAME}_SV.trf.vcf -outputDir ${OUTDIR} -vcf 1
#-variantconvertDir /mnt/beegfs/userdata/n_rabearivelo/mambaforge/envs/annotsv/share/python3/variantconvert/src/variantconvert/

conda deactivate
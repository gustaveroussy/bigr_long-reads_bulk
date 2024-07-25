#!/bin/bash
################################################################
## Using : sbatch run_bcftools_stats.sh
################################################################

#SBATCH --job-name=vcfstats
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --partition=shortq

# Input parameters
WDIR="/mnt/beegfs/scratch/bioinfo_core/B24018_OLBE_01"
SAMPLE_NAME="WM99"
REF="hg38"
#VCF="${WDIR}/data_output/somatic_variant_calling/SNV_annotation/snpSift/${SAMPLE_NAME}_${REF}_SNV_PASS.ann.dbnsfp.clinvar.vcf"
#VCF="${WDIR}/data_output/somatic_variant_calling/SNV_annotation/snpSift/${SAMPLE_NAME}_${REF}_SNV_indel_PASS.ann.dbnsfp.clinvar.vcf"
#VCF="${WDIR}/data_output/somatic_variant_calling/SNV_annotation/snpEff/${SAMPLE_NAME}_${REF}_SNV_indel.ann.vcf"
#VCF="${WDIR}/data_output/somatic_variant_calling/SV_calling/WM99_CD19_Tumor.nanomonsv.result.vcf"
VCF="${WDIR}/data_output/somatic_variant_calling/SNV_annotation/snpEff/${SAMPLE_NAME}_${REF}_SNV.ann.vcf"

# Create output directory
mkdir -p ${WDIR}/data_output/somatic_variant_calling/QC/bcftools/

# BCF tools vcf stats
module load bcftools
#bcftools stats ${VCF} > ${WDIR}/data_output/somatic_variant_calling/QC/bcftools/VCF-STATS_${SAMPLE_NAME}_${REF}_SNV_PASS.ann.dbnsfp.clinvar.txt
#bcftools stats ${VCF} > ${WDIR}/data_output/somatic_variant_calling/QC/bcftools/VCF-STATS_${SAMPLE_NAME}_${REF}_SNV_indel_PASS.ann.dbnsfp.clinvar.txt
#bcftools stats ${VCF} > ${WDIR}/data_output/somatic_variant_calling/QC/bcftools/VCF-STATS_${SAMPLE_NAME}_${REF}_SNV_indel.ann.txt
#bcftools stats ${VCF} > ${WDIR}/data_output/somatic_variant_calling/QC/bcftools/VCF-STATS_${SAMPLE_NAME}_${REF}_SV.txt
bcftools stats ${VCF} > ${WDIR}/data_output/somatic_variant_calling/QC/bcftools/VCF-STATS_${SAMPLE_NAME}_${REF}_SNV.ann.txt


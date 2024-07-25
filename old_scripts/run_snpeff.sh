#!/bin/bash
################################################################
## Using : sbatch run_snpeff.sh
################################################################

#SBATCH --job-name=annotation_BC35
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --partition=mediumq

# Input parameters
WDIR="/mnt/beegfs/scratch/bioinfo_core/B23043_NADR_02"
SAMPLE_NAME="BC35_control"
REF="hg38"
#VCF="${WDIR}/data_output/somatic_variant_calling/SNV_calling/WM99_SNV_indel.vcf.gz"
VCF="${WDIR}/data_output/dorado/${SAMPLE_NAME}/clair3_HAC_${REF}/${SAMPLE_NAME}_trimmed_dorado_${REF}_SNV.vcf.gz"

# Create output directory
mkdir -p ${WDIR}/data_output/dorado/${SAMPLE_NAME}/snpEff/
mkdir -p ${WDIR}/data_output/dorado/${SAMPLE_NAME}/snpSift/

module load java

# ------ Execute snpEff
SNPEFF="/mnt/beegfs/userdata/n_rabearivelo/snpEff/snpEff.jar"
java -jar ${SNPEFF} GRCh38.105 ${VCF} -csvStats ${WDIR}/data_output/dorado/${SAMPLE_NAME}/snpEff/${SAMPLE_NAME}_${REF}_SNV_snpEff_stats.csv -stats ${WDIR}/data_output/dorado/${SAMPLE_NAME}/snpEff/${SAMPLE_NAME}_${REF}_SNV_snpEff_summary.html > ${WDIR}/data_output/dorado/${SAMPLE_NAME}/snpEff/${SAMPLE_NAME}_${REF}_SNV.ann.vcf
#java -jar ${SNPEFF} GRCh38.105 ${VCF} -csvStats ${WDIR}/data_output/somatic_variant_calling/SNV_annotation/snpEff/${SAMPLE_NAME}_${REF}_SNV_indel_snpEff_stats.csv -stats ${WDIR}/data_output/somatic_variant_calling/SNV_annotation/snpEff/${SAMPLE_NAME}_${REF}_SNV_indel_snpEff_summary.html > ${WDIR}/data_output/somatic_variant_calling/SNV_annotation/snpEff/${SAMPLE_NAME}_${REF}_SNV_indel.ann.vcf

#Get only the summary and genes files
java -jar ${SNPEFF} GRCh38.105 ${VCF} -stats ${WDIR}/data_output/dorado/${SAMPLE_NAME}/snpEff/${SAMPLE_NAME}_${REF}_SNV_snpEff_summary.html
#java -jar ${SNPEFF} GRCh38.105 ${VCF} -stats ${WDIR}/data_output/${SAMPLE_NAME}/SNV_annotation/snpEff/${SAMPLE_NAME}_${REF}_SNV_indel_snpEff_summary.html

# ------ Filter annotated VCF
SNPSIFT="/mnt/beegfs/userdata/n_rabearivelo/snpEff/SnpSift.jar"
#VCF_ANNOT="${WDIR}/data_output/somatic_variant_calling/SNV_annotation/snpEff/${SAMPLE_NAME}_${REF}_SNV.ann.vcf"
VCF_ANNOT="${WDIR}/data_output/dorado/${SAMPLE_NAME}/snpEff/${SAMPLE_NAME}_${REF}_SNV.ann.vcf"

#cat ${VCF_ANNOT} | java -jar ${SNPSIFT} filter "( FILTER = 'PASS' )" > ${WDIR}/data_output/somatic_variant_calling/SNV_annotation/snpSift/${SAMPLE_NAME}_${REF}_SNV_PASS.ann.vcf
#cat ${VCF_ANNOT} | java -jar ${SNPSIFT} filter "( FILTER = 'PASS' )" > ${WDIR}/data_output/somatic_variant_calling/SNV_annotation/snpSift/${SAMPLE_NAME}_${REF}_SNV.ann.vcf

#-------- SNPSIFT (same as run_snpsift.sh)

# 1 - Annotation (of the VCF previously filtered to keep only the SNV with "PASS" value in the Filter column) using DBNSFP
#VCF_ANNOT="${WDIR}/data_output/somatic_variant_calling/SNV_annotation/snpSift/${SAMPLE_NAME}_${REF}_SNV_PASS.ann.vcf"
default_fields="1000Gp3_AC,1000Gp3_AF,1000Gp3_AFR_AC,1000Gp3_AFR_AF,1000Gp3_AMR_AC,1000Gp3_AMR_AF,1000Gp3_EAS_AC,1000Gp3_EAS_AF,1000Gp3_EUR_AC,1000Gp3_EUR_AF,1000Gp3_SAS_AC,1000Gp3_SAS_AF,ESP6500_AA_AC,ESP6500_AA_AF,ESP6500_EA_AC,ESP6500_EA_AF,ExAC_AC,ExAC_AF,ExAC_AFR_AC,ExAC_AFR_AF,ExAC_AMR_AC,ExAC_AMR_AF,ExAC_Adj_AC,ExAC_Adj_AF,ExAC_EAS_AC,ExAC_EAS_AF,ExAC_FIN_AC,ExAC_FIN_AF,ExAC_NFE_AC,ExAC_NFE_AF,ExAC_SAS_AC,ExAC_SAS_AF,FATHMM_pred,GERP++_NR,GERP++_RS,Interpro_domain,LRT_pred,MetaSVM_pred,MutationAssessor_pred,MutationTaster_pred,PROVEAN_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,SIFT_pred,Uniprot_acc,phastCons100way_vertebrate,CADD_phred"
supp_fields="CADD_raw,CADD_raw_rankscore,DANN_score,DANN_rankscore,REVEL_score,REVEL_rankscore,ClinPred_score,ClinPred_rankscore,ClinPred_pred,MutPred_score,MutPred_protID,MutPred_AAchange,VEP_canonical,gnomAD_genomes_flag,gnomAD_genomes_AC,gnomAD_genomes_AN,gnomAD_genomes_AF,GTEx_V8_gene,GTEx_V8_tissue"
#java -jar ${SNPSIFT} dbnsfp -v -db /mnt/beegfs/database/bioinfo/Index_DB/dbNSFP/4.1/GRCh38/dbNSFP4.1a.txt.gz -f ${default_fields},${supp_fields} ${VCF_ANNOT} > ${WDIR}/data_output/somatic_variant_calling/SNV_annotation/snpSift/${SAMPLE_NAME}_${REF}_SNV_PASS.ann.dbnsfp.vcf
java -jar ${SNPSIFT} dbnsfp -v -db /mnt/beegfs/database/bioinfo/Index_DB/dbNSFP/4.1/GRCh38/dbNSFP4.1a.txt.gz -f ${default_fields},${supp_fields} ${VCF_ANNOT} > ${WDIR}/data_output/dorado/${SAMPLE_NAME}/snpSift/${SAMPLE_NAME}_${REF}_SNV.ann.dbnsfp.vcf

# 2 - Annotation (of the VCF previously annotated with DBNSFP) using CLINVAR
#VCF_ANNOT="${WDIR}/data_output/somatic_variant_calling/SNV_annotation/snpSift/${SAMPLE_NAME}_${REF}_SNV_PASS.ann.dbnsfp.vcf"
VCF_ANNOT="${WDIR}/data_output/dorado/${SAMPLE_NAME}/snpSift/${SAMPLE_NAME}_${REF}_SNV.ann.dbnsfp.vcf"
#java -jar ${SNPSIFT} annotate -v /mnt/beegfs/database/bioinfo/Index_DB/ClinVar/GRCh38/clinvar_20230710.vcf.gz ${VCF_ANNOT} > ${WDIR}/data_output/somatic_variant_calling/SNV_annotation/snpSift/${SAMPLE_NAME}_${REF}_SNV_PASS.ann.dbnsfp.clinvar.vcf
java -jar ${SNPSIFT} annotate -v /mnt/beegfs/database/bioinfo/Index_DB/ClinVar/GRCh38/clinvar_20230710.vcf.gz ${VCF_ANNOT} > ${WDIR}/data_output/dorado/${SAMPLE_NAME}/snpSift/${SAMPLE_NAME}_${REF}_SNV.ann.dbnsfp.clinvar.vcf
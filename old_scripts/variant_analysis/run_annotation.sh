#!/bin/bash
################################################################
## Using : sbatch run_annotation.sh
################################################################

#SBATCH --job-name=3582_65NG_filter_VCF
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=40G
#SBATCH --partition=mediumq

# Input information
#WDIR="/mnt/beegfs/scratch/bioinfo_core/B23043_NADR_02"
#SAMPLE_NAME="3582_R10_65NG"
#VCF="/mnt/beegfs/scratch/bioinfo_core/B23043_NADR_02/data_output/3700_R10/clair3_HAC_hg38/3700_R10_hg38_merge_output_Q15.vcf"
REF="hg38"

# Reference genomes (FASTA)
#HG38="/mnt/beegfs/userdata/n_rabearivelo/references/Ensembl/GRCh38.109/homo_sapiens.GRCh38.109.fasta"
#T2T="/mnt/beegfs/userdata/n_rabearivelo/references/Ensembl/T2T-CHM13v2.0/Homo_sapiens-GCA_009914755.4-unmasked.fa"

# ------ Extract chr20
#module load samtools
#module load sambamba

#samtools view -bo ${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_T2T-CHM13.chr20.bam ${BAM} 20
#sambamba sort ${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_T2T-CHM13.chr20.bam
#sambamba index ${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_T2T-CHM13.chr20.bam

#samtools view -bo ${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_hg38.chr20.bam ${BAM} 20
#sambamba sort ${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_hg38.chr20.bamis
#sambamba index ${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_hg38.chr20.bam

# ----- Variant calling + Haplotyping
# R10 model - Guppy 6.5.7 - HAC / Rerio model (considered deprecated since it's using Guppy), use the latest when analyzing Dorado basecalled data
#CLAIR3_MODEL="/mnt/beegfs/userdata/n_rabearivelo/rerio/clair3_models/r1041_e82_400bps_hac_g632"
#CLAIR3_THREADS="10"
#CONTIGS_LIST="20"
#source /mnt/beegfs/userdata/n_rabearivelo/mambaforge/etc/profile.d/conda.sh
#conda activate clair3
#run_clair3.sh --bam_fn=${BAM} --ref_fn=${T2T} --threads=${CLAIR3_THREADS} --platform="ont" --model_path="${CLAIR3_MODEL}" --output="${WDIR}/data_output/${SAMPLE_NAME}/clair3_HAC_T2T/chr20/" --ctg_name=${CONTIGS_LIST} --use_whatshap_for_final_output_haplotagging --enable_long_indel --enable_phasing
#BAM="${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_hg38.chr20.sorted.bam"
#run_clair3.sh --bam_fn=${BAM} --ref_fn=${T2T} --threads=${CLAIR3_THREADS} --platform="ont" --model_path="${CLAIR3_MODEL}" --output="${WDIR}/data_output/${SAMPLE_NAME}/clair3_HAC_hg38/chr20/" --use_whatshap_for_final_output_haplotagging --enable_long_indel --enable_phasing
#conda deactivate

# --------

# Create output directory
mkdir -p ${WDIR}/data_output/dorado/${SAMPLE_NAME}/snpEff/
mkdir -p ${WDIR}/data_output/dorado/${SAMPLE_NAME}/snpSift/

module load java

# ------ Execute snpEff
SNPEFF="/mnt/beegfs/userdata/n_rabearivelo/snpEff/snpEff.jar"
#gunzip ${VCF}
mv ${WDIR}/data_output/dorado/${SAMPLE_NAME}/clair3_HAC_${REF}/merge_output.vcf.gz ${WDIR}/data_output/dorado/${SAMPLE_NAME}/clair3_HAC_${REF}/${SAMPLE_NAME}_trimmed_dorado_${REF}_SNV.vcf.gz
VCF="${WDIR}/data_output/dorado/${SAMPLE_NAME}/clair3_HAC_${REF}/${SAMPLE_NAME}_trimmed_dorado_${REF}_SNV.vcf.gz"
java -jar ${SNPEFF} GRCh38.105 ${VCF} -csvStats ${WDIR}/data_output/dorado/${SAMPLE_NAME}/snpEff/${SAMPLE_NAME}_${REF}_merged_snpEff_stats.csv -stats ${WDIR}/data_output/dorado/${SAMPLE_NAME}/snpEff/${SAMPLE_NAME}_${REF}_merged_snpEff_summary.html > ${WDIR}/data_output/dorado/${SAMPLE_NAME}/snpEff/${SAMPLE_NAME}_${REF}_merged.ann.vcf
#java -jar ${SNPEFF} GRCh38.105 ${VCF} -csvStats ${WDIR}/data_output/${SAMPLE_NAME}/snpEff/${SAMPLE_NAME}_${REF}_merged_snpEff_summary.csv
#cat ${WDIR}/data_output/${SAMPLE_NAME}_dorado/snpEff/${SAMPLE_NAME}_${REF}_merged.ann.vcf | java -jar /mnt/beegfs/userdata/n_rabearivelo/snpEff/SnpSift.jar filter "(FILTER = 'PASS')" > ${WDIR}/data_output/${SAMPLE_NAME}_dorado/snpEff/${SAMPLE_NAME}_merge_output_PASS.vcf
#cat ${WDIR}/data_output/${SAMPLE_NAME}_dorado/snpEff/${SAMPLE_NAME}_${REF}_merged.ann.vcf| java -jar /mnt/beegfs/userdata/n_rabearivelo/snpEff/SnpSift.jar filter "( QUAL >= 10 )" > ${WDIR}/data_output/${SAMPLE_NAME}_dorado/snpEff/${SAMPLE_NAME}_merge_output_Q10.vcf
#cat ${WDIR}/data_output/${SAMPLE_NAME}_dorado/snpEff/${SAMPLE_NAME}_${REF}_merged.ann.vcf | java -jar /mnt/beegfs/userdata/n_rabearivelo/snpEff/SnpSift.jar filter "( QUAL >= 15 )" > ${WDIR}/data_output/${SAMPLE_NAME}_dorado/snpEff/${SAMPLE_NAME}_merge_output_Q15.vcf

#Get only the summary and genes files !
#java -jar ${SNPEFF} GRCh38.105 ${VCF} -stats ${WDIR}/data_output/${SAMPLE_NAME}/snpEff/${SAMPLE_NAME}_${REF}_merged_snpEff_summary.html

# ------ Execute snpSift
#SNPSIFT="/mnt/beegfs/userdata/n_rabearivelo/snpEff/SnpSift.jar"
#VCF_ANNOT="${WDIR}/data_output/${SAMPLE_NAME}/snpEff/${SAMPLE_NAME}_${REF}_merged.ann.vcf"

#dbNSFP v4.1
#java -jar ${SNPSIFT} dbnsfp -v -db /mnt/beegfs/database/bioinfo/Index_DB/dbNSFP/4.1/GRCh38/dbNSFP4.1a.txt.gz ${VCF_ANNOT}

#dbNSFP v4.5
#java -jar ${SNPSIFT} dbnsfp -v -db /mnt/beegfs/userdata/n_rabearivelo/snpEff/dbNSFP4.5a.txt.gz ${VCF_ANNOT} > ${WDIR}/data_output/${SAMPLE_NAME}/snpSift/${SAMPLE_NAME}_hg38_merged.ann.dbnsfp.vcf
# http://database.liulab.science/dbNSFP#version

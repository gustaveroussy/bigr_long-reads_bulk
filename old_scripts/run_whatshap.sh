#!/bin/bash
################################################################
## Using : sbatch run_whatshap.sh
## https://whatshap.readthedocs.io/en/latest/guide.html#subcommands
################################################################

#SBATCH --job-name=3700_R10_hg38_whatshap_phasing
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=4
#SBATCH --mem=20G
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

# Input
#BAM="${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_hg38.sorted.bam"
#VCF="${WDIR}/data_output/${SAMPLE_NAME}/clair3_HAC/merge_output.vcf.gz"

# Ensembl Reference genomes
#HG38="/mnt/beegfs/userdata/n_rabearivelo/references/Ensembl/GRCh38.109/homo_sapiens.GRCh38.109.fasta"
#HG38_minimap2_Index="/mnt/beegfs/userdata/n_rabearivelo/references/Ensembl/GRCh38.109/homo_sapiens.GRCh38.109.mmi"
#T2T="/mnt/beegfs/userdata/n_rabearivelo/references/Ensembl/T2T-CHM13v2.0/Homo_sapiens-GCA_009914755.4-unmasked.fa"
#T2T_minimap2_Index="/mnt/beegfs/userdata/n_rabearivelo/references/Ensembl/T2T-CHM13v2.0/Homo_sapiens-GCA_009914755.4-unmasked.mmi"

# Create output directory
mkdir -p ${WDIR}/data_output/dorado/${SAMPLE_NAME}/whatshap/

# Use WhatsHap
source /mnt/beegfs/userdata/n_rabearivelo/mambaforge/etc/profile.d/conda.sh
conda activate whatshap

# Phasing
#whatshap phase --output ${WDIR}/data_output/${SAMPLE_NAME}/whatshap/${SAMPLE_NAME}_hg38_variants.phased.vcf.gz --reference ${HG38} --mapping-quality 20 --ignore-read-groups ${VCF} ${BAM}
#BAM="/mnt/beegfs/scratch/bioinfo_core/B23043_NADR_02/data_output/3700_R10/sambamba/chr20/3700_R10_hg38.chr20.sorted.bam"
#VCF="/mnt/beegfs/scratch/bioinfo_core/B23043_NADR_02/data_output/3700_R10/clair3_HAC_hg38/chr20/merge_output.vcf.gz"
#bgzip ${VCF}

#whatshap phase --output ${WDIR}/data_output/dorado/${SAMPLE_NAME}/whatshap/${SAMPLE_NAME}_${REF}_SNV.phased.vcf.gz --reference ${ref_fasta} --mapping-quality 20 --ignore-read-groups ${VCF} ${BAM}
#vcf file should be in vcf.gz format (bgzip file.vcf)

#tabix ${WDIR}/data_output/dorado/${SAMPLE_NAME}/whatshap/${SAMPLE_NAME}_${REF}_SNV.phased.vcf.gz

# Stats
# TSV file can be used by MultiQC
#whatshap stats --tsv ${WDIR}/data_output/${SAMPLE_NAME}/whatshap/chr20/${SAMPLE_NAME}_hg38_phasing_stats.tsv --block-list=${WDIR}/data_output/${SAMPLE_NAME}/whatshap/chr20/${SAMPLE_NAME}_hg38_haplotype_blocks.tsv --gtf=${WDIR}/data_output/${SAMPLE_NAME}/whatshap/chr20/${SAMPLE_NAME}_hg38_haplotype_blocks.gtf ${WDIR}/data_output/${SAMPLE_NAME}/whatshap/chr20/${SAMPLE_NAME}_hg38_variants.phased.chr20.vcf.gz > ${WDIR}/data_output/${SAMPLE_NAME}/whatshap/chr20/${SAMPLE_NAME}_hg38_general_phasing_stats.txt

# Haplotagging : Tagging reads by haplotype for visualization
# Show the reads along with the variants
# It tags each read in a BAM file with HP:i:1 or HP:i:2 depending on which haplotype it belongs to + adds a PS tag describing in which haplotype block the read is.
#whatshap haplotag --output-threads=4 --ignore-read-groups --output-haplotag-list ${WDIR}/data_output/dorado/${SAMPLE_NAME}/whatshap/${SAMPLE_NAME}_hg38_SNV_meth_haplotag_list.tsv  --output ${WDIR}/data_output/dorado/${SAMPLE_NAME}/whatshap/${SAMPLE_NAME}_hg38_SNV_meth.haplotagged.bam --reference ${ref_fasta} ${WDIR}/data_output/dorado/${SAMPLE_NAME}/whatshap/${SAMPLE_NAME}_${REF}_SNV.phased.vcf.gz ${BAM}
# output haplotag list : to use if we want to use whatshap split after

conda deactivate
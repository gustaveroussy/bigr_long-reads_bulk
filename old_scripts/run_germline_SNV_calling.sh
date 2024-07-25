#!/bin/bash
################################################################
## Variant Calling Tools Testing
## Using : sbatch run_variant_calling.sh
################################################################
 
#SBATCH --job-name=SNV_calling_clair3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=10
#SBATCH --mem=16G
#SBATCH --partition=longq

# Input information
#WDIR="/mnt/beegfs/scratch/bioinfo_core/B23043_NADR_02"

#SAMPLE_NAME="3683_CD14plus_ADNg_01082022"
#SAMPLE_NAME="3700_R10"
#SAMPLE_NAME="3700_CD14plus_ADNg_25072022"
#SAMPLE_NAME="3582_R10_250NG"
SAMPLE_NAME="3582_R10_1000NG"

BAM="${WDIR}/data_output/${SAMPLE_NAME}/sambamba/${SAMPLE_NAME}_hg38.sorted.bam"

# R9 model - Guppy 6.5.7 - SUP & HAC : https://github.com/HKU-BAL/Clair3/blob/main/docs/guppy5_20220113.md "use the sup model for hac data. But use sup data for the best performance." + You can use the Guppy5 model r941_prom_sup_g5014 for Guppy6 data.
# Clair3 included model
#CLAIR3_MODEL="r941_prom_sup_g5014"

# R10 model - Guppy 6.5.7 - HAC
# Rerio model (considered deprecated since it's using Guppy), use the latest when analyzing Dorado basecalled data
#CLAIR3_MODEL="/mnt/beegfs/userdata/n_rabearivelo/rerio/clair3_models/r1041_e82_400bps_hac_g632"

# Latest R10 model from Rerio for SUP / downloaded on 2024-02-25
CLAIR3_MODEL="/mnt/beegfs/userdata/n_rabearivelo/rerio/clair3_models/r1041_e82_400bps_sup_v430"
CLAIR3_THREADS="10"

# Reference genomes (FASTA)
#HG38="/mnt/beegfs/userdata/n_rabearivelo/references/Ensembl/GRCh38.109/homo_sapiens.GRCh38.109.fasta"
#T2T="/mnt/beegfs/userdata/n_rabearivelo/references/Ensembl/T2T-CHM13v2.0/Homo_sapiens-GCA_009914755.4-unmasked.fa"

# Clair3 - SNP & Long indels Calling
source /mnt/beegfs/userdata/n_rabearivelo/mambaforge/etc/profile.d/conda.sh
conda activate clair3

run_clair3.sh --bam_fn=${BAM} --ref_fn=${ref_fasta} --threads=${CLAIR3_THREADS} --platform="ont" --model_path="${CLAIR3_MODEL}" --output="${WDIR}/data_output/dorado/${SAMPLE_NAME}/clair3_HAC_${REF}/" --include_all_ctgs
# By default, Clair3 searches for chr{1..22,X,Y} and {1..22,X,Y}. If you are using a non human genome reference or a human ref with different chromosomes/contigs names, add the option "--include_all_ctgs".
# Add long indels detection (EXPERIMENTAL) : --enable_long_indel

conda deactivate

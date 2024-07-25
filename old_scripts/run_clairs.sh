#!/bin/bash
################################################################
## Usage : sbatch run_clairs.sh
## https://github.com/HKU-BAL/ClairS#installation
################################################################

#SBATCH --job-name=clairS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=10
#SBATCH --mem=20G
#SBATCH --partition=mediumq
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

# Input parameters
WDIR="/mnt/beegfs/scratch/bioinfo_core/B24018_OLBE_01"
CLAIRS_THREADS="10"
PLATFORM="ont_r10_dorado_sup_5khz" # options: {ont_r10_dorado_sup_4khz, ont_r10_dorado_sup_5khz, ont_r10_guppy, ont_r9_guppy, ilmn, hifi_sequel2, hifi_revio, etc}

BAM_PATH="${WDIR}/data_output/somatic_variant_calling/concat_chr"
REF="/mnt/beegfs/userdata/n_rabearivelo/references/Ensembl/GRCh38.109"

OUTPUT_DIR="${WDIR}/data_output/somatic_variant_calling/SNV_calling/"

# Create output directory
mkdir -p ${OUTPUT_DIR}

# ClairS - Somatic Variant Calling
module load singularity

singularity exec \
  -B ${WDIR},${BAM_PATH},${REF},${OUTPUT_DIR} \
  /mnt/beegfs/userdata/n_rabearivelo/containers/clairs_latest.sif \
  /opt/bin/run_clairs \
  --tumor_bam_fn ${BAM_PATH}/WM99_CD19_all_chr_concat_filtered_sorted.bam \
  --normal_bam_fn ${BAM_PATH}/WM99_CD3_all_chr_concat_filtered_sorted.bam \
  --ref_fn ${REF}/homo_sapiens.GRCh38.109.fasta \
  --threads ${CLAIRS_THREADS} \
  --platform ${PLATFORM} \
  --output_dir ${OUTPUT_DIR} \
  --output_prefix WM99_SNV_indel \
  --ctg_name 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT \
  --conda_prefix /opt/conda/envs/clairs \
  --enable_indel_calling

#!/bin/bash
################################################################
## MultiQC
## Using : sbatch run_multiqc.sh
################################################################

#SBATCH --job-name=multiqc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --partition=shortq

# Input parameters
WDIR="/mnt/beegfs/scratch/bioinfo_core/B24018_OLBE_01"
SAMPLE_NAME="WM99"
REF="hg38"

source /mnt/beegfs/userdata/n_rabearivelo/mambaforge/etc/profile.d/conda.sh
conda activate multiqc

multiqc ${WDIR}/data_output/somatic_variant_calling/QC/1_For_MultiQC/ --outdir ${WDIR}/data_output/somatic_variant_calling/ --force --title "Quality Control - After alignment" --comment "Reads were aligned against hg38 human reference sequence (release 109)." --config /mnt/beegfs/pipelines/long-reads-bulk/v0_test/config/multiqc_config.yaml

conda deactivate

#!/bin/bash
################################################################
## Fastq Screen - Search for contamination
## Using : sbatch run_fastq_screen.sh
################################################################

#SBATCH --job-name=fastq_screen_3582_65NG
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=4
#SBATCH --mem=20G
#SBATCH --partition=mediumq

# Input information
WDIR="/mnt/beegfs/scratch/bioinfo_core/B23043_NADR_02"
SAMPLE_NAME="3582_R10_65NG"
FASTQ="${WDIR}/data_input/${SAMPLE_NAME}/${SAMPLE_NAME}.fastq.gz"

#BARCODES_LIST=$(ls ${WDIR}/data_input/concat/*.fastq.gz | sed 's,/mnt/beegfs/scratch/n_rabearivelo/NAPE_QC/data_input/concat/,,g' | sed 's,.fastq.gz,,g')
#BARCODES_LIST=$(ls ${WDIR}/data_input/GUPPY/fastq_pass/)
FASTQ_SCREEN="/mnt/beegfs/userdata/n_rabearivelo/FastQ-Screen/fastq_screen"
CONF="/mnt/beegfs/userdata/n_rabearivelo/references/FastQ_Screen_Genomes/fastq_screen.conf"

mkdir -p ${WDIR}/data_output/${SAMPLE_NAME}/fastq_screen/hg38/
${FASTQ_SCREEN} --aligner minimap2 --threads 8 --outdir ${WDIR}/data_output/${SAMPLE_NAME}/fastq_screen/hg38/ --force --conf ${CONF} ${FASTQ}

#for BARCODE in ${BARCODES_LIST}; do
    # Create output dir
    #mkdir -p ${WDIR}/data_output/DORADO/fastq_screen/${BARCODE}/
    # Use fastq_screen/v0.15.3 branch "minimap"
    #${FASTQ_SCREEN} --aligner minimap2 --threads 4 --outdir ${WDIR}/data_output/DORADO/fastq_screen/${BARCODE}/ --force --conf ${CONF} ${WDIR}/data_input/DORADO/fastq_from_ubam/SQK-NBD114-24_${BARCODE}.fastq.gz
#done;

#################################

# ---- No Need to re-run ---- #

## Create all indexes for minimap2
#minimap2 -d /mnt/beegfs/userdata/n_rabearivelo/references/FastQ_Screen_Genomes/minimap2/Adapters/Contaminants.mmi /mnt/beegfs/userdata/n_rabearivelo/references/FastQ_Screen_Genomes/minimap2/Adapters/Contaminants.fasta
#minimap2 -d /mnt/beegfs/userdata/n_rabearivelo/references/FastQ_Screen_Genomes/minimap2/Vectors/UniVec.mmi /mnt/beegfs/userdata/n_rabearivelo/references/FastQ_Screen_Genomes/minimap2/Vectors/UniVec.fasta
#minimap2 -d /mnt/beegfs/userdata/n_rabearivelo/references/FastQ_Screen_Genomes/minimap2/Lambda_Phage/lambda_phage.mmi /mnt/beegfs/userdata/n_rabearivelo/references/FastQ_Screen_Genomes/minimap2/Lambda_Phage/lambda_phage.fasta
#minimap2 -d /mnt/beegfs/userdata/n_rabearivelo/references/FastQ_Screen_Genomes/minimap2/PhiX/PhiX.mmi /mnt/beegfs/userdata/n_rabearivelo/references/FastQ_Screen_Genomes/minimap2/PhiX/PhiX.fasta
#minimap2 -d /mnt/beegfs/userdata/n_rabearivelo/references/FastQ_Screen_Genomes/minimap2/E_coli/Ecoli.mmi /mnt/beegfs/userdata/n_rabearivelo/references/FastQ_Screen_Genomes/minimap2/E_coli/Ecoli.fasta
#minimap2 -d /mnt/beegfs/userdata/n_rabearivelo/references/FastQ_Screen_Genomes/minimap2/Arabidopsis/Arabidopsis_thaliana.TAIR10.mmi /mnt/beegfs/userdata/n_rabearivelo/references/FastQ_Screen_Genomes/minimap2/Arabidopsis/Arabidopsis_thaliana.TAIR10.fasta
#minimap2 -d /mnt/beegfs/userdata/n_rabearivelo/references/FastQ_Screen_Genomes/minimap2/Yeast/S288C_R64-4-1_20230830.mmi /mnt/beegfs/userdata/n_rabearivelo/references/FastQ_Screen_Genomes/minimap2/Yeast/S288C_R64-4-1_20230830.fasta
#minimap2 -d /mnt/beegfs/userdata/n_rabearivelo/references/FastQ_Screen_Genomes/minimap2/Worm/Caenorhabditis_elegans.WBcel235.mmi /mnt/beegfs/userdata/n_rabearivelo/references/FastQ_Screen_Genomes/minimap2/Worm/Caenorhabditis_elegans.WBcel235.fasta
#minimap2 -d /mnt/beegfs/userdata/n_rabearivelo/references/FastQ_Screen_Genomes/minimap2/Drosophila/Drosophila_melanogaster.BDGP6.mmi /mnt/beegfs/userdata/n_rabearivelo/references/FastQ_Screen_Genomes/minimap2/Drosophila/Drosophila_melanogaster.BDGP6.fasta
#minimap2 -d /mnt/beegfs/userdata/n_rabearivelo/references/FastQ_Screen_Genomes/minimap2/Rat/Rattus_norvegicus.mRatBN7.mmi /mnt/beegfs/userdata/n_rabearivelo/references/FastQ_Screen_Genomes/minimap2/Rat/Rattus_norvegicus.mRatBN7.fasta
#minimap2 -d /mnt/beegfs/userdata/n_rabearivelo/references/FastQ_Screen_Genomes/minimap2/Mouse/Mus_musculus.GRCm39.mmi /mnt/beegfs/userdata/n_rabearivelo/references/FastQ_Screen_Genomes/minimap2/Mouse/Mus_musculus.GRCm39.fasta
#minimap2 -d /mnt/beegfs/scratch/n_rabearivelo/NAPE_QC/data_input/NanoForkSpeed_ref_genome/S288CwExtrarDNA_ROMAN.mmi /mnt/beegfs/scratch/n_rabearivelo/NAPE_QC/data_input/NanoForkSpeed_ref_genome/S288CwExtrarDNA_ROMAN.fa.gz

#################################

#Perl module GD::Graph::bars not installed, skipping charts
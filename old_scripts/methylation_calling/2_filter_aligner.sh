#!/bin/bash
################################################################
##DORADO ALIGN AFTER DORADO BASECALLING
################################################################

#SBATCH --job-name=filter_alignement_LR
#SBATCH --nodes=1
#SBATCH --mem=25G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=8
#SBATCH --partition=mediumq

#should change to longq

#MODULE LOAD
module load samtools/1.11
module load python

#WORKKING DIRECTORY AND RESOURCES PATH
WKDIR="/mnt/beegfs/scratch/bioinfo_core/B24018_OLBE_01" #no slash at the end of path
DATASET="WM100_CD19" #empty if only one sample


DORADO="/mnt/beegfs/pipelines/dorado"
REF="/mnt/beegfs/userdata/n_rabearivelo/references/Ensembl/T2T-CHM13v2.0/Homo_sapiens-GCA_009914755.4-unmasked.fa"
INPUT=${WKDIR}/data_output/${DATASET}

echo "Relocating in data_output and creating subfolder"
cd ${INPUT}
mkdir 2_custom_align

echo "Alignement to custom ref for each bam"
echo $REF

for bam in 1_meth_call/batch{85..89}*.bam
do
	echo $bam
	time python3 /mnt/beegfs/userdata/y_mesloub/script/filter_bam.py $bam
	${DORADO}/0.5.3/dorado-0.5.3-linux-x64/bin/dorado aligner -t 8 --bandwidth 500,20000 --secondary=no ${REF} ${bam}_filtered.bam > 2_custom_align/${bam##*/}_custom_align.bam
done

echo "Finished!!"



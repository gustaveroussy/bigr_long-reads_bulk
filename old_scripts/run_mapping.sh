#!/bin/bash
################################################################
## MINIMAP2 : Mapping step
## Using : sbatch run_mapping.sh
################################################################

#SBATCH --job-name=250NG_trimmed_split_Merge_Index
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=2
#SBATCH --mem=20G
#SBATCH --partition=mediumq
#SBATCH --output=MAPPING-%x.%j.out
#SBATCH --error=MAPPING-%x.%j.err

#Mem usually = 16G

# Input information
#WDIR="/mnt/beegfs/scratch/bioinfo_core/B23043_NADR_02"
#SAMPLE_NAME="3683_CD14plus_ADNg_01082022"
#SAMPLE_NAME="3683_partial"
#FASTQ="${WDIR}/data_input/3683_CD14plus_ADNg_01082022/3683_CD14plus_ADNg_01082022/20220801_1123_6C_PAM59295_dd30a530/test/${SAMPLE_NAME}.fastq.gz"

#SAMPLE_NAME="3700_CD14plus_ADNg_25072022"
#FASTQ="${WDIR}/data_input/${SAMPLE_NAME}/${SAMPLE_NAME}/20220725_1340_6A_PAM62320_0eacc38a/concat_fastq/${SAMPLE_NAME}.fastq.gz"

#SAMPLE_NAME="3700_R10"
#FASTQ="${WDIR}/data_input/${SAMPLE_NAME}/${SAMPLE_NAME}.fastq.gz"

#SAMPLE_NAME="3582_R10_250NG"
#SAMPLE_NAME="3666_ultra_long_reads"
#SAMPLE_NAME="BC35_control"
#SAMPLE_NAME="2572_CD14"
#SAMPLE_NAME="3582_R10_1000NG"

#FASTQ="${WDIR}/data_input/${SAMPLE_NAME}/${SAMPLE_NAME}.fastq.gz"

# Ensembl Reference genomes
HG38="/mnt/beegfs/userdata/n_rabearivelo/references/Ensembl/GRCh38.109/homo_sapiens.GRCh38.109.fasta"
HG38_minimap2_Index="/mnt/beegfs/userdata/n_rabearivelo/references/Ensembl/GRCh38.109/homo_sapiens.GRCh38.109.mmi"
T2T="/mnt/beegfs/userdata/n_rabearivelo/references/Ensembl/T2T-CHM13v2.0/Homo_sapiens-GCA_009914755.4-unmasked.fa"
T2T_minimap2_Index="/mnt/beegfs/userdata/n_rabearivelo/references/Ensembl/T2T-CHM13v2.0/Homo_sapiens-GCA_009914755.4-unmasked.mmi"

# Create output dir
#mkdir -p ${WDIR}/data_output/${SAMPLE_NAME}/minimap2/

# Launch Minimap2 to perform alignment (16-->14G, 8-->2 CPU per task)
#minimap2 --secondary=no -x map-ont -a ${mmi_index} ${FILTERED_FASTQ} > ${WDIR}/data_output/${SAMPLE_NAME}/minimap2/${SAMPLE_NAME}_${REF}.sam

# Launch Minimap2 to perform alignment (16-->14G, 8-->2 CPU per task)
#minimap2 --secondary=no -x map-ont -a ${mmi_index} ${FASTQ_DIR}/*.part_001.fastq.gz > ${WDIR}/data_output/${SAMPLE_NAME}_dorado/minimap2/${SAMPLE_NAME}_${REF}.part_001.sam

# Input information
WDIR="/mnt/beegfs/scratch/bioinfo_core/B23043_NADR_02"
SAMPLE_NAME="3582_R10_250NG"
REF="hg38"
#FASTQ_DIR="${WDIR}/data_input/${SAMPLE_NAME}/dorado_rebasecalling/fastq"
FASTQ_DIR=" ${WDIR}/data_output/${SAMPLE_NAME}_dorado/chopper"

# Create output directory
mkdir -p ${WDIR}/data_output/${SAMPLE_NAME}_dorado/minimap2/
mkdir -p ${WDIR}/data_output/${SAMPLE_NAME}_dorado/sambamba/

# Use Minimap2 & Sambamba
#for i in {10..78}; do
    #for FASTQ in ${FASTQ_DIR}/${SAMPLE_NAME}_noadapter_minlen200_q10_batch${i}.fastq.gz; do
        #BATCH=$(echo ${FASTQ} | cut -d"_" -f14 | cut -d"." -f1)
        #BATCH="batch${i}"
        
        #jobid1=$(sbatch --job-name=250NG-trimmed-vs-${REF}-${BATCH}-minimap2 --ntasks=4 --mem=20G --partition=shortq --output=MAPPING-%x.%j.out --error=MAPPING-%x.%j.err --wrap "minimap2 --secondary=no -x map-ont -a ${HG38_minimap2_Index} ${FASTQ} > ${WDIR}/data_output/${SAMPLE_NAME}_dorado/minimap2/${SAMPLE_NAME}_trimmed_${REF}_${BATCH}.sam")
        #jobid1=$(echo ${jobid1} | cut -d" " -f4)
        
        #SAM="${WDIR}/data_output/${SAMPLE_NAME}_dorado/minimap2/${SAMPLE_NAME}_trimmed_${REF}_${BATCH}.sam"
        #jobid2=$(sbatch --job-name=250NG-trimmed-vs-HG38-${BATCH}-sam2bam --ntasks=4 --mem=5G --partition=shortq --output=MAPPING-%x.%j.out --error=MAPPING-%x.%j.err --dependency=afterok:${jobid1} --wrap "module load sambamba; sambamba view -S -f bam ${SAM} > ${WDIR}/data_output/${SAMPLE_NAME}_dorado/sambamba/${SAMPLE_NAME}_trimmed_${REF}_${BATCH}.bam")
        #jobid2=$(echo ${jobid2} | cut -d" " -f4)
        
        #BAM="${WDIR}/data_output/${SAMPLE_NAME}_dorado/sambamba/${SAMPLE_NAME}_trimmed_${REF}_${BATCH}.bam"
        #module load sambamba; sambamba sort ${BAM}
        #sbatch --job-name=250NG-trimmed-vs-${REF}-${BATCH}-sort --ntasks=4 --mem=10G --partition=shortq --output=MAPPING-%x.%j.out --error=MAPPING-%x.%j.err --dependency=afterok:${jobid2} --wrap "module load sambamba; sambamba sort ${BAM}"
    #done
#done

jobid_merge=$(sbatch --job-name=1000NG-trimmed-vs-${REF}-${BATCH}-merge --ntasks=4 --mem=24G --partition=mediumq --output=MAPPING-%x.%j.out --error=MAPPING-%x.%j.err --wrap "module load samtools; samtools merge -f ${WDIR}/data_output/${SAMPLE_NAME}_dorado/sambamba/${SAMPLE_NAME}_trimmed_${REF}_merged.sorted.bam ${WDIR}/data_output/${SAMPLE_NAME}_dorado/sambamba/${SAMPLE_NAME}_trimmed_${REF}_*.sorted.bam")
jobid_merge=$(echo ${jobid_merge} | cut -d" " -f4)
sbatch --job-name=1000NG-trimmed-MERGED-vs-${REF}-${BATCH}-index --ntasks=4 --mem=24G --partition=mediumq --output=MAPPING-%x.%j.out --error=MAPPING-%x.%j.err --dependency=afterok:${jobid_merge} --wrap "module load sambamba; sambamba index ${WDIR}/data_output/${SAMPLE_NAME}_dorado/sambamba/${SAMPLE_NAME}_trimmed_${REF}_merged.sorted.bam"


#module load sambamba; sambamba index ${WDIR}/data_output/${SAMPLE_NAME}_dorado/sambamba/${SAMPLE_NAME}_trimmed_${REF}_merged.sorted.bam

# ---- No Need to re-run ---- #

# Create index for both references
#samtools faidx ${HG38}
#samtools faidx ${T2T}

# Indexing reference files for Minimap2 : 12G, 1CPU/task
# Save Index file for HG38 (before using it for the alignment step, replacing the fasta reference file)
#minimap2 -d ${HG38_minimap2_Index} ${HG38}

# Save Index file for T2T (before using it for the alignment step, replacing the fasta reference file)
#minimap2 -d ${T2T_minimap2_Index} ${T2T}
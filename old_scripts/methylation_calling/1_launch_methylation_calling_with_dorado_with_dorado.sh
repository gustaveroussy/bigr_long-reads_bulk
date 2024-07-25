#!/bin/bash
################################################################
##BASECALLING AND MODIFIED BASE DETECTION ON DORADO SUP R10
################################################################

#SBATCH --job-name=BC_R10_ONT_CG
#SBATCH --mem=30gb
#SBATCH --time=50:00:00
#SBATCH --cpus-per-task=2
#SBATCH --partition=gpgpuq
#SBATCH --gres=gpu:a100:1


#export cuda location
export LD_LIBRARY_PATH=/usr/local/cuda-11.1/targets/x86_64-linux/lib:$LD_LIBRARY_PATH
export PATH=/usr/local/cuda/bin:$PATH

echo $1

#WORKING DIRECTORY AND RESOURCES
WKDIR=""    #no slash at the end of path
DORADO="/mnt/beegfs/pipelines/dorado"  #no slash at the end of path
DATASET_FOLDER="" #name of dataset in data_input #no slash at the end

if [[ -z $1 ]];then echo $1; echo "$1 is empty";fi

#MODELS AND KIT 
MODEL="${DORADO}/model/model_r10/dna_r10.4.1_e8.2_400bps_sup@v4.2.0"
METH_MOD="${DORADO}/model/model_r10/dna_r10.4.1_e8.2_400bps_sup@v4.2.0_5mCG_5hmCG@v2"


#CHANGING LOCATION TO PROJECT PATH AND CREATING SUBDIRECTORY
cd ${WKDIR}
mkdir -m 771 -p ${WKDIR}/data_output/${DATASET_FOLDER}/1_meth_call


#creation des batchs d'input (de taille d'environ 40G)
#echo "Creating batches of pod5 file in tmp folder"
mkdir -m 771 ${WKDIR}/tmp/${DATASET_FOLDER}
find ${WKDIR}/data_input/${DATASET_FOLDER}/ -type f -name "*.pod5" |split -d -l 1 - ${WKDIR}/tmp/${DATASET_FOLDER}/batch
for batch in `ls ${WKDIR}/tmp/${DATASET_FOLDER}/batch*`
do
	echo ${batch}; mkdir ${batch}_folder
	for pod5 in `cat ${batch}`;do ln -s -t ${batch}_folder $pod5;done
done

#INPUT FOLDERS
INPUT="${WKDIR}/tmp/${DATASET_FOLDER}/b*_folder"

echo -e "Launching of methylation calling for 5mCG and 5hmCG...\n"
#LAUNCHING METHYLATION CALLING (SEQUENTIAL)
for folder in ${INPUT}; do echo -e "\n\n${folder}";date;${DORADO}/0.5.3/dorado-0.5.3-linux-x64/bin/dorado basecaller ${MODEL} ${folder} --modified-bases-models ${METH_MOD} -x cuda:0 --emit-moves --no-trim -b 1728 > ${WKDIR}/data_output/${DATASET_FOLDER}/1_meth_call/${folder##*/}_meth_call.bam;done

echo "Finished"

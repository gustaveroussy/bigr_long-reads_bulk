# bigr_long-reads_bulk
This pipeline allows to analyse bulk long-read DNA-seq from PromethION (Oxford Nanopore Technologies) to identify SNV and SV for germline and somatic data.
For now, the input files are bam with their index (.bai) files.

## Installation
### Flamingo
The pipeline and environments are already installed on the Flamingo cluster of Gustave Roussy.  
The pipeline is localized here: /mnt/beegfs/pipelines/bigr_long-reads_bulk/<version>

### Installation outside GR or for new version
#### Download pipeline
```
VERSION="1.0.0"
cd  /mnt/beegfs/pipelines/bigr_long-reads_bulk/
git clone https://github.com/gustaveroussy/bigr_long-reads_bulk.git ${VERSION}
```
#### Download environments
```
source /mnt/beegfs/software/miniconda/24.3.0/etc/profile.d/conda.sh
conda create --prefix="/mnt/beegfs/userdata/<user_id>/.environnement_conda/git_lfs" git-lfs
conda activate /mnt/beegfs/userdata/<user_id>/.environnement_conda/git_lfs
cd ${VERSION}
git lfs install
git lfs pull
```
#### Install snakemake environment
```
source /mnt/beegfs/software/miniconda/24.3.0/etc/profile.d/conda.sh
conda env create -f /mnt/beegfs/pipelines/bigr_long-reads_bulk/${VERSION}/envs/conda/snakemake.yaml --prefix=/mnt/beegfs/pipelines/bigr_long-reads_bulk/${VERSION}/envs/conda/snakemake -y
```

## Using
You need to make 2 files: a design file and a configuration file. 

### Configuration file
You can copy and modify the example from config/config.yaml.  

### Design file
It must be a comma separated file (.csv where comma is ",").
If you wante a germline analysis:
- **sample_id**: the sample name of you sample (it could be different that your fastq files).
- **bam_file**: absolute path to the bam file.
Example:
```
sample_id,bam_file
3700_R10,/mnt/beegfs/userdata/m_aglave/long-reads-bulk/test/data_input/3700_R10_chr_22_filtered.bam
2572_CD14,/mnt/beegfs/userdata/m_aglave/long-reads-bulk/test/data_input/2572_CD14_chr_22_filtered.bam
```
If you wante a somatic analysis:
- **sample_id**: the sample name of you sample (it could be different that your fastq files).
- **bam_file_tumor**: absolute path to the tumor bam file.
- **bam_file_normal**: absolute path to the normal bam file.

Example:
```
sample_id,bam_file_tumor,bam_file_normal
somatic_test_data_bam,/mnt/beegfs/userdata/m_aglave/long-reads-bulk/test/data_input/3700_R10_chr_22_filtered.bam,/mnt/beegfs/userdata/m_aglave/long-reads-bulk/test/data_input/2572_CD14_chr_22_filtered.bam
```

> Notes:
> - sample names mustn't contain special characters or spaces.
> - each bam file must be sorted, and have its bai index in the same directory than them.

### Run
You need snakemake (via conda) and singularity (via module load). They are already installed for you on Flamingo, just follow the example below.  
Don't forget to change the version of the pipeline and the path to your configuration file.  
Example of script:
```
#!/bin/bash
#using: sbatch run.sh
#SBATCH --job-name=LR_analysis
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=250M
#SBATCH --partition=longq

source /mnt/beegfs/software/miniconda/24.3.0/etc/profile.d/conda.sh
conda activate /mnt/beegfs/pipelines/bigr_long-reads_bulk/<version>/envs/conda/snakemake
module load singularity

LR_pipeline="/mnt/beegfs/pipelines/bigr_long-reads_bulk/<version>/"

snakemake --profile ${LR_pipeline}/profiles/slurm \
          -s ${LR_pipeline}/Snakefile \
          --configfile <path_to/my_configuration_file.yaml>
```

import pandas as pd
import os

sys.stderr.write("\n############################################################ \n\n")

sys.stderr.write("                           ++++++               \n")
sys.stderr.write("                        ++++++++++++            \n")
sys.stderr.write("                     ++++++++++++++++++         \n")
sys.stderr.write("                  +++++++++++++++++++++++++     \n")
sys.stderr.write("               ++++++++++++++++++++++++++++++   \n")
sys.stderr.write("               ++++|  _ \(_)/ ____|  __ \++++   \n")
sys.stderr.write("               ++++| |+) |_| |++++| |++) |+++   \n")
sys.stderr.write("               ++++|  _ <| | |+|_ |  _  /++++   \n")
sys.stderr.write("               ++++| |+) | | |++| | |+\ \++++   \n")
sys.stderr.write("               ++++|____/|_|\_____|_|++\_\+++   \n")
sys.stderr.write("               ++++++++++++++++++++++++++++++   \n")
sys.stderr.write("                  +++++++++++++++++++++++++     \n")
sys.stderr.write("                     +++++++++++++++++++        \n")
sys.stderr.write("                        +++++++++++++           \n")
sys.stderr.write("                           +++++++              \n")

sys.stderr.write("\n                Bulk long reads pipeline      \n")

sys.stderr.write("\nFor any question, sent an email to bigr@gustaveroussy.fr")
sys.stderr.write("\n############################################################ \n\n")

sys.stderr.write("\n#################### Setting Parameters ####################\n\n")

#pipeline directory
PIPELINE_DIR = workflow.snakefile
PIPELINE_DIR = PIPELINE_DIR.replace("/Snakefile", "")
c
#output directory (working directory)
if "output_dir" in config:
    OUTPUT_DIR = config["output_dir"]
    sys.stderr.write("Output directory set to " + OUTPUT_DIR + ".\n")
else:
    OUTPUT_DIR = os.getcwd()
    sys.stderr.write("No 'output_dir' found in config file, output directory set to " + OUTPUT_DIR + ".\n")


#set enviroements and parameters of steps
if config["steps"]["bam_qc"]:
  CONDA_ENV_QUALITYMAP = PIPELINE_DIR + "/envs/conda/qualimap.yaml"
  CONDA_ENV_MOSDEPTH = PIPELINE_DIR + "/envs/conda/mosdepth.yaml"
  SING_ENV_NANOPLOT = PIPELINE_DIR + "/envs/singularity/nanoplot_1.42.0.simg"
  CONDA_ENV_MULTIQC = PIPELINE_DIR + "/envs/conda/multiqc.yaml"
  
#if config["steps"]["fastq_qc"]:

if config["steps"]["snv_calling"]:
  CONDA_ENV_CLAIR3 = PIPELINE_DIR + "/envs/conda/clair3.yaml"
  SING_ENV_PEPPER_DEEPVARIANT = PIPELINE_DIR + "/envs/singularity/pepper_deepvariant_r0.8.sif"
  NAME_CLAIR3_MODEL = [os.path.basename(model) for model in config["clair3"]["model"]]
  SING_ENV_SNPEFF = PIPELINE_DIR + "/envs/singularity/snpeff_5.2c.simg"
  SNPSIFT_FILTERS =  [filter.replace("'", "\'").replace('"', '\"') for filter in config["snpsift"]["filters"]]
  SNPSIFT_FILTERS_NAMES = [filter.replace("'", "").replace("(", "").replace(")", "").replace(" ", "").replace(">", "sup").replace("<", "inf").replace("=", "eq") for filter in config["snpsift"]["filters"]] #to do: find a way to keep "()" if therer are more than 2 parentheses. Ex: "((QUAL >= 10) && (QUAL <= 30)) || (FILTER = 'PASS')"
  SING_ENV_CLAIRS = PIPELINE_DIR + "/envs/singularity/clairs_0.1.7.sif"

if config["steps"]["sv_calling"]:
  CONDA_ENV_SNIFFLES = PIPELINE_DIR + "/envs/conda/sniffles.yaml"
  CONDA_ENV_NANOMONSV = PIPELINE_DIR + "/envs/conda/nanomonsv.yaml"

if config["steps"]["phasing"]:
  CONDA_ENV_WHATSHAP = PIPELINE_DIR + "/envs/conda/whatshap.yaml"


sys.stderr.write("Parameters validated.\n")



sys.stderr.write("\n################### Checking Design File ###################\n\n")

#read design file
design=pd.read_table(config["design"],sep=",")
#design=pd.read_table("/mnt/beegfs/userdata/m_aglave/long-reads-bulk/test/script/design_bam_germline.tsv",sep=",")

#check if input_format and variant_calling_mode are agree with design format
if config["input_format"] == "bam" and config["variant_calling_mode"] == "germline" : format_design = ['sample_id', 'bam_file']
if config["input_format"] == "bam" and config["variant_calling_mode"] == "somatic" : format_design = ['sample_id', 'bam_file_tumor', 'bam_file_normal']
if config["input_format"] == "fastq" and config["variant_calling_mode"] == "germline" : format_design = ['sample_id', 'upstream_fastq_file', 'downstream_fastq_file']
if config["input_format"] == "fastq" and config["variant_calling_mode"] == "somatic" : format_design = ['sample_id', 'upstream_fastq_file_tumor', 'downstream_fastq_file_tumor', 'upstream_fastq_file_normal', 'downstream_fastq_file_normal']
if set(format_design).issubset(design.columns):
  sys.stderr.write("Design file well formated.\n")
  design=design[format_design]
else: sys.exit("Error in the format of the design file: missing column(s).")

#check if all sample_id are differents
if not len(set(design["sample_id"])) == len(design["sample_id"]): sys.exit("Error: All sample_id have to be different! Check your design file.")

#check if all files are differents on each line
for i in range(0, len(design["sample_id"]), 1):
  if not len(set(design.iloc[i])) == len(design.iloc[i]): sys.exit("Error: All files have to be different! Check your design file.")

#make the data structure
SAMPLE_NAME =  []
BAM_NAME =  []
ORIG_FILE = []
SYMLINK_FILES = []
for col in design.drop(columns="sample_id").columns.tolist():
  for line in range(0,len(design[col].tolist()),1):
    SAMPLE_NAME.append(design["sample_id"].iloc[line])
    ORIG_FILE.append(design[col].iloc[line])
    if col == "bam_file" : SUPPL_NAME=""
    if col == "bam_file_tumor": SUPPL_NAME="_tumor"
    if col == "bam_file_normal": SUPPL_NAME="_normal"
    BAM_NAME.append(design["sample_id"].iloc[line] + SUPPL_NAME)
    SYMLINK_FILES.append(os.path.normpath(OUTPUT_DIR + "/symlink_input/" + design["sample_id"].iloc[line] + SUPPL_NAME + ".bam"))
    
#Note: for germline SAMPLE_NAME = BAM_NAME


sys.stderr.write("\n########################### Run ############################\n\n")


include: "rules/rule_all.smk"
rule all:
    input:
        **get_targets()
    message:
        "Pipeline finished!"



# Include rules files
if config["steps"]["bam_qc"]:
    include: "rules/bam_qc.smk"
    include: "rules/multiqc.smk"

if config["steps"]["snv_calling"] and config["variant_calling_mode"] == "germline":
    include: "rules/germline_snv_calling.smk"
    include: "rules/snv_annotation.smk"

if config["steps"]["sv_calling"] and config["variant_calling_mode"] == "germline":
    include: "rules/germline_sv_calling.smk"

if config["steps"]["snv_calling"] and config["variant_calling_mode"] == "somatic":
    include: "rules/somatic_snv_calling.smk"
    include: "rules/snv_annotation.smk"

if config["steps"]["sv_calling"] and config["variant_calling_mode"] == "somatic":
    include: "rules/somatic_sv_calling.smk"

if config["steps"]["phasing"]:
    include: "rules/phasing.smk"



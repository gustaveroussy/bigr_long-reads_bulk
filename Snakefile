import pandas as pd
import os
import glob
import itertools
import re
import sys

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

# Pipeline directory
PIPELINE_DIR = workflow.snakefile
PIPELINE_DIR = PIPELINE_DIR.replace("/Snakefile", "")

# Output directory (working directory)
if "output_dir" in config:
    OUTPUT_DIR = config["output_dir"]
    sys.stderr.write("Output directory set to " + OUTPUT_DIR + ".\n")
else:
    OUTPUT_DIR = os.getcwd()
    sys.stderr.write("No 'output_dir' found in config file, output directory set to " + OUTPUT_DIR + ".\n")


###################  Get fasta reference and chromosomes list ###################

if "references" not in config or "genome" not in config["references"]:
    sys.exit("No reference fasta found for 'genome'. Check your configuration file.")
#else:
    #sys.stderr.write("The reference used is:")
    #sys.stderr.write(config["references"]["genome"] + "\n\n")
if "references" in config and "chromosomes" in config["references"]:
    CHR_NUMBER = config["references"]["chromosomes"]
else: sys.exit("No chromosomes list found for 'chromosomes'. Check your configuration file.")

# Function to get chromosomes list from fasta file (but problem if reference contains alternatives chromosomes)
def chrom_extract(fasta_file):
    fasta_sequences=open(fasta_file,"r")
    chrom=[i.strip(">").split(" ")[0] for i in fasta_sequences if "chromosome:" in i]
    return chrom
#CHR_NUMBER=['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '2', '20', '21', '22', '3', '4', '5', '6', '7', '8', '9', 'MT', 'X', 'Y']
#CHR_NUMBER=chrom_extract(config["references"]["genome"])

###################  Set environments and parameters of the various pipeline steps ###################

if "input_format" not in config: sys.exit("Error: 'input_format' have to be set. Check your configuration file.")
if config["input_format"] not in ["pod5","ubam","bam"]: sys.exit("Error: 'input_format' unknown. Check your configuration file.")

if "basecalling_mode" not in config: sys.exit("Error: 'basecalling_mode' have to be set. Check your configuration file.")
if config["basecalling_mode"] not in ["methylation","basic"]: sys.exit("Error: 'basecalling_mode' unknown. Check your configuration file.")

if config["basecalling_mode"] == "basic":
    if config["steps"]["differential_methylation_sample"] or config["steps"]["differential_methylation_condition"]: sys.exit("Error: 'differential_methylation' can't be done if 'basecalling_mode' is 'basic'. Check your configuration file.")

# to do: check steps -for example, if 'phasing:true', 'snv_calling' has to be true too.

# BASECALLING
if config["steps"]["basecalling"]:

    # environments
    TOOL_DORADO = "/mnt/beegfs/pipelines/dorado/0.5.3/dorado-0.5.3-linux-x64/bin/dorado"
    
    # parameters
    if 'basecalling_mode' in config:
        if config['basecalling_mode'] == "basic":
            if 'dorado_basecaller' in config and 'model' in config['dorado_basecaller']:
                DORADO_MODEL = [config['dorado_basecaller']['model'], ""]
            else:
                DORADO_MODEL = ["/mnt/beegfs/pipelines/dorado/model/model_r10/dna_r10.4.1_e8.2_400bps_sup@v4.2.0", ""]
        elif config['basecalling_mode'] == "methylation":
            if 'dorado_basecaller' in config and 'model' in config['dorado_basecaller']:
                MODE_BASIC = config['dorado_basecaller']['model']
            else:
                MODE_BASIC = "/mnt/beegfs/pipelines/dorado/model/model_r10/dna_r10.4.1_e8.2_400bps_sup@v4.2.0"
            if 'dorado_basecaller' in config and 'model_meth' in config['dorado_basecaller']:
                MODEL_METH = config['dorado_basecaller']['model_meth']
            else:
                MODEL_METH = "/mnt/beegfs/pipelines/dorado/model/model_r10/dna_r10.4.1_e8.2_400bps_sup@v4.2.0_5mCG_5hmCG@v2"
            DORADO_MODEL = [MODE_BASIC,MODEL_METH]
        else: sys.exit("Error: 'basecalling_mode' unknown. Check your configuration file.")
    else: sys.exit("Error: 'basecalling_mode' have to be set. Check your configuration file.")


# ALIGNMENT
if config["steps"]["alignment"]:

    # environments
    TOOL_DORADO = "/mnt/beegfs/pipelines/dorado/0.5.3/dorado-0.5.3-linux-x64/bin/dorado"

    # parameters
    # to do: check reference

# QUALITY-CONTROL
if config["steps"]["alignment"] or config["input_format"] == "bam" :

    # environments
    #PATH_PYTHON3 = "/usr/bin/python3"
    CONDA_ENV_PYTHON = PIPELINE_DIR + "/envs/conda/python_3.9.1.yaml"
    CONDA_ENV_SAMTOOLS = PIPELINE_DIR + "/envs/conda/samtools_1.11.yaml"
    CONDA_ENV_QUALITYMAP = PIPELINE_DIR + "/envs/conda/qualimap.yaml"
    CONDA_ENV_MOSDEPTH = PIPELINE_DIR + "/envs/conda/mosdepth.yaml"
    SING_ENV_NANOPLOT = PIPELINE_DIR + "/envs/singularity/nanoplot_1.42.0.simg"
    SING_ENV_FASTQ_SCREEN = PIPELINE_DIR + "/envs/singularity/fastq_screen_v0.15.3_minimap2.simg"
    SING_ENV_FASTQC = PIPELINE_DIR + "/envs/singularity/fastqc_1.20.simg"
    CONDA_ENV_MULTIQC = PIPELINE_DIR + "/envs/conda/multiqc.yaml"

    # methylation
    if config['basecalling_mode'] == "methylation":
        # environments
        SING_ENV_GENEDMR = "/mnt/beegfs/pipelines/bigr_long-reads_bulk/dev_test/envs/singularity/GeneDMR.simg"
        # parameters
        def list_meth_type(list_of_model):
        	type_meth=[ re.split("@v[0-9\.]+",x)[-2].split("_")[1:] for x in list_of_model[1:]]
        	return list(chain(*type_meth))
        if 'dorado_basecaller' in config and 'meth_type' in config['dorado_basecaller'] :
            METH_TYPE = config['dorado_basecaller']['meth_type']
        elif MODEL_METH in locals() :
            METH_TYPE = list_meth_type(DORADO_MODEL)
        else: sys.exit("Error: 'meth_type' in 'dorado_basecaller' have to be set. Check your configuration file.")
        print(METH_TYPE)
    

# PREPROCESSING AND METHYLATION ANALYSIS
if config["steps"]["differential_methylation_sample"] or config["steps"]["differential_methylation_condition"]:

    # environments
    TOOL_MODKIT = "/mnt/beegfs/pipelines/dorado/tools/modkit_v0.2.6/dist/modkit"
    SING_ENV_GENEDMR = "/mnt/beegfs/pipelines/bigr_long-reads_bulk/dev_test/envs/singularity/GeneDMR.simg"
    CONDA_ENV_SAMTOOLS = PIPELINE_DIR + "/envs/conda/samtools_1.11.yaml"
    
    # parameters
    # Get methylation analysis type(s) from the list of models
    def list_meth_type(list_of_model):
    	type_meth=[ re.split("@v[0-9\.]+",x)[-2].split("_")[1:] for x in list_of_model[1:]]
    	return list(chain(*type_meth))
    if 'dorado_basecaller' in config and 'meth_type' in config['dorado_basecaller'] :
        METH_TYPE = config['dorado_basecaller']['meth_type']
    elif MODEL_METH in locals() :
        METH_TYPE = list_meth_type(DORADO_MODEL)
    else: sys.exit("Error: 'meth_type' in 'dorado_basecaller' have to be set. Check your configuration file.")
    print(METH_TYPE)
    REF_TYPE=["Alu","Transcript","CpG"]
    STRAND=["fwd","rev"]
    
    
# VARIANT ANALYSIS
if 'variant_calling_mode' in config:
    if config["steps"]["snv_calling"]:
    
        # germline
        if config['variant_calling_mode'] == "germline":
            # environments
            CONDA_ENV_CLAIR3 = PIPELINE_DIR + "/envs/conda/clair3.yaml"
            SING_ENV_PEPPER_DEEPVARIANT = PIPELINE_DIR + "/envs/singularity/pepper_deepvariant_r0.8.sif"
            # parameters
            if 'clair3' in config and 'model' in config['clair3']:
                CLAIR3_MODEL = config["clair3"]["model"]
            else:
                CLAIR3_MODEL = ["/mnt/beegfs/database/bioinfo/bigr_long-reads_bulk/MODELS/rerio/clair3_models/r1041_e82_400bps_sup_v430"]
            NAME_CLAIR3_MODEL = [os.path.basename(model) for model in CLAIR3_MODEL]
        # somatic
        elif config['variant_calling_mode'] == "somatic":
            # environments
            SING_ENV_CLAIRS = PIPELINE_DIR + "/envs/singularity/clairs_latest.sif"
            
        else: sys.exit("Error: 'variant_calling_mode' unknown. Check your configuration file.")
        
        # germline & somatic
        # environments
        SING_ENV_SNPEFF = PIPELINE_DIR + "/envs/singularity/snpeff_5.2c.simg"
        SING_ENV_VCF2MAF = PIPELINE_DIR + "/envs/singularity/vcf2maf_1.6.22.simg"
        CONDA_ENV_MAFTOOLS = PIPELINE_DIR + "/envs/conda/maftools.yaml"
        # parameters
        if 'snpsift' in config and 'filters' in config['snpsift']:
            FILTERS = config["snpsift"]["filters"]
        else:
            FILTERS = ["(FILTER = 'PASS')"]
        SNPSIFT_FILTERS = [
            filter.replace("'", "'").replace('"', '"')
            for filter in FILTERS
        ]
        SNPSIFT_FILTERS_NAMES = [
            filter.replace("'", "")
            .replace("(", "")
            .replace(")", "")
            .replace(" ", "")
            .replace(">", "sup")
            .replace("<", "inf")
            .replace("=", "eq")
            for filter in FILTERS
        ]  # to do: find a way to keep "()" if therer are more than 2 parentheses. Ex: "((QUAL >= 10) && (QUAL <= 30)) || (FILTER = 'PASS')"
        if 'maftools' in config and 'genes_of_interest' in config['maftools']:
            GENES_FILE = config["maftools"]["genes_of_interest"]
        else:
            print("No genes file found.")
            GENES_FILE = ""

        SNPEFF_SUFFIX = ""
        DBNSFP_SUFFIX = ""
        CLINVAR_SUFFIX = ""
        if "references" in config:
            if "genome_snpEff" in config["references"] and config["references"]["genome_snpEff"] != "": SNPEFF_SUFFIX = "_annotated"
            if "dbnsfp" in config["references"] and config["references"]["dbnsfp"] != "": DBNSFP_SUFFIX = "_dbnsfp"
            if "clinvar" in config["references"] and config["references"]["clinvar"] != "": CLINVAR_SUFFIX = "_clinvar"

    if config["steps"]["sv_calling"]:
        
        # germline
        if config['variant_calling_mode'] == "germline":
            # environments
            CONDA_ENV_SNIFFLES = PIPELINE_DIR + "/envs/conda/sniffles.yaml"
            CONDA_ENV_CUTESV = PIPELINE_DIR + "/envs/conda/cutesv.yaml"
            CONDA_ENV_ANNOTSV = PIPELINE_DIR + "/envs/conda/annotsv.yaml"
        
        #somatic
        elif config['variant_calling_mode'] == "somatic":
            # environments
            CONDA_ENV_NANOMONSV = PIPELINE_DIR + "/envs/conda/nanomonsv.yaml"
            CONDA_ENV_ANNOTSV = PIPELINE_DIR + "/envs/conda/annotsv.yaml"
        else: sys.exit("Error: 'variant_calling_mode' unknown. Check your configuration file.")
        
else: sys.exit("Error: 'variant_calling_mode' have to be set. Check your configuration file.")

if config["steps"]["phasing"]:
    # environments
    CONDA_ENV_WHATSHAP = PIPELINE_DIR + "/envs/conda/whatshap.yaml"


# CNV ANALYSIS
if config["steps"]["cnv_calling"]:
    # environments
    CONDA_ENV_MOSDEPTH = PIPELINE_DIR + "/envs/conda/mosdepth.yaml"
    CONDA_ENV_SPECTRE = PIPELINE_DIR + "/envs/conda/spectre.yaml"
    # parameters
    EXTRA_PARAMS_SPECTRE = ""
    if "spectre" in config:
        if "blacklist" in config["spectre"] and config["spectre"]["blacklist"] != "":
                if os.path.exists(config["spectre"]["blacklist"]):
                    EXTRA_PARAMS_SPECTRE = EXTRA_PARAMS_SPECTRE + "--blacklist " + config["spectre"]["blacklist"]
                else: sys.exit("Error: 'blacklist' file not found. Check your configuration file.")
        if config["spectre"]["cancer_sample"]:
            EXTRA_PARAMS_SPECTRE = EXTRA_PARAMS_SPECTRE + " --cancer "


sys.stderr.write("Parameters validated.\n")

sys.stderr.write("\n################### Checking Design File ###################\n\n")
# Check if input_format, differential_methylation_condition and variant_calling_mode are consistent with design format

# Read design file
if "design" in config:
    design = pd.read_table(config["design"], sep=",")
else: sys.exit("Error: 'design' have to be set. Check your configuration file.")

# FASTQ as input file
if config["input_format"] == "fastq":
    if config["steps"]["differential_methylation_condition"]:
        if config["variant_calling_mode"] == "somatic":
            format_design = ["sample_id","upstream_fastq_file","downstream_fastq_file","methyl_group","somatic_ctrl"]
        else:
            format_design = ["sample_id","upstream_fastq_file","downstream_fastq_file","methyl_group"]
    else:
        if config["variant_calling_mode"] == "somatic":
            format_design = ["sample_id","upstream_fastq_file","downstream_fastq_file","somatic_ctrl"]
        else:
            format_design = ["sample_id", "upstream_fastq_file", "downstream_fastq_file"]

#POD5, BAM, UBAM as input file
else:
    if config["steps"]["differential_methylation_condition"]:
        if config["variant_calling_mode"] == "somatic":
            format_design = ["sample_id", "path_file", "methyl_group", "somatic_ctrl"]
        else:
            format_design = ["sample_id","path_file","methyl_group"]
    else:
        if config["variant_calling_mode"] == "somatic":
            format_design = ["sample_id", "path_file", "somatic_ctrl"]
        else:
            format_design = ["sample_id","path_file"]

# Check design file:
if set(format_design).issubset(design.columns):
    sys.stderr.write("Design file well formated.\n")
    design = design[format_design]
else:
    sys.exit("Error in the format of the design file: missing column(s) or wrong column name(s). Check your design file.")

# Check if all sample_id are differents:
if not len(set(design["sample_id"])) == len(design["sample_id"]):
    sys.exit("Error: All sample_id have to be different! Check your design file.")

# Check if all files are differents and if they exist:
if config["input_format"] == "fastq":
    for i in range(0, len(design["sample_id"]), 1):
        if design["upstream_fastq_file"].iloc[i] == design["downstream_fastq_file"].iloc[i]:
            print("Error: upstream_fastq_file and downstream_fastq_file have to be different:")
            print("Line " + str(i+1) + " of design file: " + design["upstream_fastq_file"].iloc[i] + " same as " + design["downstream_fastq_file"].iloc[i])
            sys.exit("Check your design file!")
        if not os.path.isfile(design["upstream_fastq_file"].iloc[i]): sys.exit(design["upstream_fastq_file"].iloc[i] + " file doesn't exist. Check your design file.")
        if not os.path.isfile(design["downstream_fastq_file"].iloc[i]): sys.exit(design["upstream_fastq_file"].iloc[i] + " file doesn't exist. Check your design file.")
    if not len(set(design["upstream_fastq_file"])) == len(design["upstream_fastq_file"]):
        sys.exit("Error: All upstream_fastq_file have to be different! Check your design file.")
    if not len(set(design["downstream_fastq_file"])) == len(design["downstream_fastq_file"]):
        sys.exit("Error: All downstream_fastq_file have to be different! Check your design file.")
else:
    if not len(set(design["path_file"])) == len(design["path_file"]):
        sys.exit("Error: All path_file have to be different! Check your design file.")
    for i in range(0, len(design["sample_id"]), 1):
        if config["input_format"] == "pod5" and not os.path.isdir(design["path_file"].iloc[i]): sys.exit(design["path_file"].iloc[i] + " pod5 folder doesn't exist. Check your design file.")
        if config["input_format"] == "bam" and not os.path.isfile(design["path_file"].iloc[i]): sys.exit(design["path_file"].iloc[i] + " bam file doesn't exist. Check your design file.")
        if config["input_format"] == "ubam" and not os.path.isfile(design["path_file"].iloc[i]): sys.exit(design["path_file"].iloc[i] + " ubam file doesn't exist. Check your design file.")

#Check if somatic_ctrl is not the same as its sample_id
if config["variant_calling_mode"] == "somatic":
    for i in range(0, len(design["sample_id"]), 1):
        if design["sample_id"].iloc[i] == design["somatic_ctrl"].iloc[i]:
            print("Error: somatic_ctrl and sample_id have to be different on the same line.")
            print("Line " + str(i+1) + " of design file: " + design["sample_id"].iloc[i] + " same as " + design["somatic_ctrl"].iloc[i])
            sys.exit("Check your design file.")
        if not pd.isna(design["somatic_ctrl"].iloc[i]) and not design["somatic_ctrl"].iloc[i] in design["sample_id"].tolist(): 
            sys.exit("somatic_ctrl " + str(design["somatic_ctrl"].iloc[i]) + " is not in sample_id column. Check your design file")

sys.stderr.write("Design information appears correct.\n")
#to do: check if there is no "_vs_" in sample_id and methyl_group.

################### Creating needed data structures depending of data input ###################

# Function to create POD5 links
def split_pod5_size(indir, outdir, sample, max_size):
    n = 0
    size = 0
    file_name = []
    batches = []
    pod5_file = [i for i in glob.glob(indir + "/*.pod5", recursive=True)]
    for pod5 in pod5_file:
        file_name.append(pod5)
        size = size + os.stat(pod5).st_size
        if (size >= max_size) or (pod5_file.index(pod5) == (len(pod5_file) - 1)):
            batch_name = "batch_" + str(n)
            if (os.path.exists(str(outdir + "/tmp/" + sample + "/" + sample + "_" + batch_name))== False):
                os.makedirs(str(outdir + "/tmp/" + sample + "/" + sample + "_" + batch_name))
            for j in file_name:
                if (os.path.islink(str(outdir + "/tmp/"+ sample+ "/"+ sample+ "_"+ batch_name+ "/"+ os.path.basename(j)))== False):
                    os.symlink(str(j),str(outdir + "/tmp/"+ sample+ "/"+ sample+ "_"+ batch_name+ "/"+ os.path.basename(j)))
            n = n + 1
            size = 0
            file_name = []
            batches.append(sample + "_" + batch_name)
    return batches

# Make the data structure for pod5 files as input
if config["steps"]["basecalling"]:
    if config["input_format"] == "pod5":
        print("input format is POD5 files.")
        SAMPLE_NAME = []
        BATCH_NAME = []
        for line in range(0, len(design["path_file"]), 1):
            BATCH_FOLDER = split_pod5_size(design["path_file"].iloc[line],OUTPUT_DIR,design["sample_id"].iloc[line],20000000000)
            SAMPLE_NAME = SAMPLE_NAME + ([design["sample_id"].iloc[line]] * len(BATCH_FOLDER))
            BATCH_NAME = BATCH_NAME + BATCH_FOLDER
        print("BATCH_NAME:")
        print(BATCH_NAME)
        print("SAMPLE_NAME:")
        print(SAMPLE_NAME)
    else: sys.exit("To make 'basecalling', 'input_format' have to be 'pod5'. Check your design file")

# Make the data structure for UBAM files as input
if config["steps"]["alignment"]:
    if config["input_format"] == "ubam":
        print("input format is UBAM file.")
    # to to
    elif config["input_format"] == "bam": sys.exit("To make 'alignment', 'input_format' have to be 'pod5' or 'ubam'. Check your design file")

# Make the data structure for BAM files as input
if config["steps"]["differential_methylation_sample"] or config["steps"]["differential_methylation_condition"] or config["steps"]["snv_calling"] or config["steps"]["sv_calling"]:
    if config["input_format"] == "bam":
        print("input format is BAM file.")
        SAMPLE_NAME = []
        ORIG_FILE = []
        SYMLINK_FILES = []
        for line in range(0, len(design["path_file"]), 1):
            SAMPLE_NAME.append(design["sample_id"].iloc[line])
            ORIG_FILE.append(design["path_file"].iloc[line])
            SYMLINK_FILES.append(os.path.normpath(OUTPUT_DIR + "/concat_sort/" + design["sample_id"].iloc[line] + "/" + design["sample_id"].iloc[line] + "_sorted.bam"))


################### Creating needed data structures for paired analyses ###################

if config["steps"]["differential_methylation_sample"]:
    PAIR_TMP=list(itertools.permutations(list(design['sample_id']),2))
    PAIR_METHYL=["_vs_".join(elms) for elms in PAIR_TMP]
    print("PAIR_METHYL: ")
    print(PAIR_METHYL)

if config["steps"]["differential_methylation_condition"]:
    CASE=list(design[design.methyl_group == 'case']['sample_id'])
    CONTROL=list(design[design.methyl_group == 'control']['sample_id'])
    print("CASE:")
    print(CASE)
    print("CONTROL:")
    print(CONTROL)

if config["variant_calling_mode"] == "somatic":
    PAIR_SOMATIC = []
    PAIR_SOMATIC_duplicate = []
    for i in range(0, len(design["sample_id"]), 1):
        if not pd.isna(design["somatic_ctrl"].iloc[i]):
            PAIR_SOMATIC = PAIR_SOMATIC + [ design["sample_id"].iloc[i] + "_vs_" + design["somatic_ctrl"].iloc[i] ]
        else:
            PAIR_SOMATIC = PAIR_SOMATIC + [ design[design.somatic_ctrl == design["sample_id"].iloc[i]]["sample_id"].iloc[0] + "_vs_" + design["sample_id"].iloc[i] ]
    print("PAIR_SOMATIC:")
    print(PAIR_SOMATIC)


################### Wilcards Constraints ###################

wildcard_constraints:
    sample_name = '|'.join([x for x in SAMPLE_NAME]),
    path_calling_tool_params = "clair3.+|pepper_margin_deepvariant|clairs",
    compl = "_merge_output||_snv|_indel",
    path_calling_tool_params_SV = "sniffles|cuteSV|nanomonsv",
    compl_SV = "|_SV|.nanomonsv.result",

if config["steps"]["snv_calling"] and config['variant_calling_mode'] == "germline":
    wildcard_constraints:
        clair3_model = '|'.join([x for x in NAME_CLAIR3_MODEL]),

if config["steps"]["snv_calling"]:
    wildcard_constraints:
        filter = '|'.join([x for x in SNPSIFT_FILTERS_NAMES])

if config["variant_calling_mode"] == "somatic":
    wildcard_constraints:
        pair_somatic = '|'.join([x for x in PAIR_SOMATIC]),
        sample_name_and_all_samples = ("|".join(SAMPLE_NAME) + "|all_samples|" + "|".join(PAIR_SOMATIC)),
        variant_calling_mode = "Somatic"
        

if config["variant_calling_mode"] == "germline":
    wildcard_constraints:
        sample_name_and_all_samples = ("|".join(SAMPLE_NAME) + "|all_samples"),
        variant_calling_mode = "Germline"


if config["steps"]["differential_methylation_condition"]:
    wildcard_constraints:
        pair_methyl = '|'.join([x for x in PAIR_METHYL]),

#    chromos = "|".join(CHR_NUMBER)
#    batch_name = "|".join(BATCH_NAME),

################### Functions needed into rules ###################

# to get input normal bam files
def get_input_normal_bam(wildcards):
        return os.path.normpath(OUTPUT_DIR + "/reconcat/" + str(wildcards.pair_somatic).split("_vs_")[0] + "/" + str(wildcards.pair_somatic).split("_vs_")[0] + "_sorted.bam")

# to get input tumor bam files
def get_input_tumor_bam(wildcards):
        return os.path.normpath(OUTPUT_DIR + "/reconcat/" + str(wildcards.pair_somatic).split("_vs_")[1] + "/" + str(wildcards.pair_somatic).split("_vs_")[1] + "_sorted.bam")

sys.stderr.write("\n########################### Run ############################\n\n")

include: "rules/rule_all.smk"
rule all:
    input:
        **get_targets(),
    message:
        "Pipeline finished!"

# BASECALLING
if config["steps"]["basecalling"]:
    include: "rules/basecalling.smk"

# ALIGNMENT
if config["steps"]["alignment"]:
    include: "rules/filter_align.smk"

# BAM PREPARATION
if config["steps"]["alignment"] or config["input_format"] == "bam":
    include: "rules/split_concat_bam.smk"
    include: "rules/bam_qc.smk"
    include: "rules/bam2fq_qc.smk"
    include: "rules/multiqc.smk"
    if 'methylation' in config['basecalling_mode']:
        include: "rules/methylation_extract.smk"
        include: "rules/methylation_qc.smk"

# METHYLATION ANALYSIS
if config["steps"]["differential_methylation_sample"] or config["steps"]["differential_methylation_condition"]:
    include: "rules/dmr.smk"

# VARIANT ANALYSIS
if config["steps"]["snv_calling"] and config["variant_calling_mode"] == "germline":
    include: "rules/germline_snv_calling.smk"
    include: "rules/snv_annotation.smk"
    include: "rules/snv_graphs.smk"

if config["steps"]["sv_calling"] and config["variant_calling_mode"] == "germline":
    include: "rules/germline_sv_calling.smk"
    include: "rules/sv_annotation.smk"
    include: "rules/sv_graphs.smk"

if config["steps"]["snv_calling"] and config["variant_calling_mode"] == "somatic":
    include: "rules/somatic_snv_calling.smk"
    include: "rules/snv_annotation.smk"
    include: "rules/snv_graphs.smk"

if config["steps"]["sv_calling"] and config["variant_calling_mode"] == "somatic":
    include: "rules/somatic_sv_calling.smk"
    include: "rules/sv_annotation.smk"

if config["steps"]["phasing"]:
    include: "rules/phasing.smk"

if config["steps"]["cnv_calling"]:
    include: "rules/cnv_calling.smk"
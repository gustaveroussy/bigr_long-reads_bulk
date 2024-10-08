"""
##########################################################################
These rules launches ubam preprocessing, alignment and following steps
##########################################################################
"""

if config["basecalling_mode"] == "methylation" and config["input_format"] == "ubam":
    """
    This rules checks if the UBAM files contain methylation information before proceeding to the following steps
    """
    rule check_ubam_methylation:
        input: 
            ubam = os.path.normpath(OUTPUT_DIR + "/tmp/ubam/{sample_name}/{batch_name}.bam"),
            index = os.path.normpath(OUTPUT_DIR + "/tmp/ubam/{sample_name}/{batch_name}.bam.bai")
        output: 
            flag = temp(os.path.normpath(OUTPUT_DIR + "/tmp/ubam/{sample_name}/{batch_name}_check_bam_methylation_OK.txt"))
        params:
            script = os.path.normpath(PIPELINE_DIR + "/script")
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(1024 + 5120 * (attempt - 1 ),2048),
            time_min = (lambda wildcards, attempt: attempt * 60)
        conda:
            CONDA_ENV_PYTHON
        shell:
            """
            flag=$(python3 {params.script}/check_methylation_in_bam.py {input.ubam})

            if [ "$flag" = "True" ]; then
                echo "Ubam file seems to contain methylation information." > {output.flag}
            elif [ "$flag" = "False" ]; then
                echo "Ubam file does not seem to contain methylation information as MM read tag was not found for the first read from the input BAM file."
            else
                echo "It looks like the script encountered an issue while running."
            fi
            """

"""
This rule launches ubam filtering
"""
def ubam_filtering_input(wildcards):
    input = []
    if config["input_format"] == "pod5":
        input.append(os.path.normpath(OUTPUT_DIR + "/tmp/calling/{sample_name}/{batch_name}.bam"))
    elif config["input_format"] == "ubam" and config['basecalling_mode'] == "methylation":
        input.append(os.path.normpath(OUTPUT_DIR + "/tmp/ubam/{sample_name}/{batch_name}.bam"))
        input.append(os.path.normpath(OUTPUT_DIR + "/tmp/ubam/{sample_name}/{batch_name}.bam.bai"))
        input.append(os.path.normpath(OUTPUT_DIR + "/tmp/ubam/{sample_name}/{batch_name}_check_bam_methylation_OK.txt"))
    elif config["input_format"] == "ubam" and config['basecalling_mode'] == "basic":
        input.append(os.path.normpath(OUTPUT_DIR + "/tmp/ubam/{sample_name}/{batch_name}.bam"))
        input.append(os.path.normpath(OUTPUT_DIR + "/tmp/ubam/{sample_name}/{batch_name}.bam.bai"))
    return input

rule ubam_filtering:
    input: 
        ubam_filtering_input
    output: 
        ubam = temp(os.path.normpath(OUTPUT_DIR + "/tmp/filtered/{sample_name}/{batch_name}_filtered.bam"))
    params:
        script = os.path.normpath(PIPELINE_DIR + "/script"),
        min_length = MIN_READ_LENGTH,
        min_qual = MIN_READ_QUALITY_SCORE
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: min(1024 + 5120 * (attempt - 1 ),20480),
        time_min = (lambda wildcards, attempt: attempt * 300)
    conda:
        CONDA_ENV_PYTHON
    shell:
        """
        python3 {params.script}/filter_bam.py --input_bam {input[0]} --output_bam {output.ubam} --min_length {params.min_length} --min_qual {params.min_qual}
        """

"""
This rule launches the alignment step
"""
rule alignment:
    input:
        ubam = os.path.normpath(OUTPUT_DIR + "/tmp/filtered/{sample_name}/{batch_name}_filtered.bam")
    output:
        bam = temp(os.path.normpath(OUTPUT_DIR + "/tmp/alignment/{sample_name}/{batch_name}_aligned.bam"))
    params:
        reference = config["references"]["genome"],
        ref_path = os.path.normpath(config["references"]["genome"])
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: min(25600 + 5120 * (attempt - 1),35840),
        time_min = (lambda wildcards, attempt: attempt * 300)
    shell:
       	"""
       	TMP_DIR=$(mktemp -d -t lr_pipeline-XXXXXXXXXX) && \
        singularity exec --contain -B {OUTPUT_DIR},{params.ref_path} -B ${{TMP_DIR}}:/tmp \
        {SING_ENV_DORADO} dorado aligner --threads {threads} --mm2-opts "-r 500,20000 --secondary=no" {params.reference} {input.ubam} > {output.bam} 
       	"""

"""
This rule concatenates all BAM per sample
"""
def input_concat_bam(wildcards):
    indices = [i for i in range(len(SAMPLE_NAME)) if SAMPLE_NAME[i] == wildcards.sample_name]
    input = expand(os.path.normpath(OUTPUT_DIR + "/tmp/alignment/"+ wildcards.sample_name +"/{batch_name}_aligned.bam"), batch_name=[BATCH_NAME[i] for i in indices])
    return input

rule concat_bam:
    input:
        bams = input_concat_bam
    output:
        temp(os.path.normpath(OUTPUT_DIR + "/Bam/{sample_name}/{sample_name}_concat.bam"))
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: min(2048 + 5120 * (attempt - 1),10240),
        time_min = (lambda wildcards, attempt: attempt * 300)
    conda:
        CONDA_ENV_SAMTOOLS
    shell:
        """
       	samtools cat {input.bams} > {output}
        """

"""
This rule sorts the concatenated BAM
"""
rule sort:
    input:
        bam = os.path.normpath(OUTPUT_DIR + "/Bam/{sample_name}/{sample_name}_concat.bam")
    output:
        os.path.normpath(OUTPUT_DIR + "/Bam/{sample_name}/{sample_name}_sorted.bam")
    threads: 12
    resources:
        mem_mb=lambda wildcards, attempt: min(20480 + 5120 * (attempt - 1 ),35840),
        time_min = (lambda wildcards, attempt: attempt * 300)
    conda:
        CONDA_ENV_SAMTOOLS
    shell:
        """
       	samtools sort -@ 12 -o {output} {input.bam}
        """

"""
This rule indexes the sorted BAM
"""
rule index:
    input:
        bam = os.path.normpath(OUTPUT_DIR + "/Bam/{sample_name}/{sample_name}_sorted.bam")
    output:
        os.path.normpath(OUTPUT_DIR + "/Bam/{sample_name}/{sample_name}_sorted.bam.bai")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: min(1024 + 5120 * (attempt - 1 ),10240),
        time_min = (lambda wildcards, attempt: attempt * 300)
    conda:
        CONDA_ENV_SAMTOOLS
    shell:
        """
       	samtools index {input.bam}
       	"""

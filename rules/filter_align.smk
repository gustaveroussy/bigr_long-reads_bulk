"""
##########################################################################
These rules launches the alignment and following steps
##########################################################################

"""

"""
This rule launches ubam filtering
"""

rule ubam_filtering:
    input: 
        bam = os.path.normpath(OUTPUT_DIR + "/calling/{sample_name}/{batch_name}.bam")
    output: 
        temp(os.path.normpath(OUTPUT_DIR + "/filtered/{sample_name}/{batch_name}_filtered.bam"))
    params:
        script = os.path.normpath(PIPELINE_DIR + "/script")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: min(1024 + 5120 * (attempt - 1 ),20480),
        time_min = (lambda wildcards, attempt: attempt * 300)
    conda:
        CONDA_ENV_PYTHON
    shell:
        """
        python3 {params.script}/filter_bam2.py {input.bam} {output}
        """

"""
This rule launches the alignment step
"""
rule alignment:
    input:
        bam = os.path.normpath(OUTPUT_DIR + "/filtered/{sample_name}/{batch_name}_filtered.bam")
    output:
        temp(os.path.normpath(OUTPUT_DIR + "/alignment/{sample_name}/{batch_name}_aligned.bam"))
    params:
        reference = config["references"]["genome"]
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: min(25600 + 5120 * (attempt - 1),35840),
        time_min = (lambda wildcards, attempt: attempt * 300)
    shell:
       	"""
       	{TOOL_DORADO} aligner -t 8 --bandwidth 500,20000 --secondary=no {params.reference} {input.bam} > {output} 
       	"""


"""
This rule concatenates all BAM per sample
"""
def input_concat_bam(wildcards):
    indices = [i for i in range(len(SAMPLE_NAME)) if SAMPLE_NAME[i] == wildcards.sample_name]
    input = expand(os.path.normpath(OUTPUT_DIR + "/alignment/"+ wildcards.sample_name +"/{batch_name}_aligned.bam"), batch_name=[BATCH_NAME[i] for i in indices])
    return input

rule concat_bam:
    input:
        bams = input_concat_bam
    output:
        temp(os.path.normpath(OUTPUT_DIR + "/concat_sort/{sample_name}/{sample_name}_concat.bam"))
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: min(2048 + 5120 * (attempt - 1),10240),
        time_min = (lambda wildcards, attempt: attempt * 300)
    conda:
        CONDA_ENV_SAMTOOLS
    shell:
        """
        module load samtools/1.11
       	samtools cat {input.bams} > {output}
        """

"""
This rule sorts the concatenated BAM
"""
rule sort:
    input:
        bam = os.path.normpath(OUTPUT_DIR + "/concat_sort/{sample_name}/{sample_name}_concat.bam")
    output:
        os.path.normpath(OUTPUT_DIR + "/concat_sort/{sample_name}/{sample_name}_sorted.bam")
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
        bam = os.path.normpath(OUTPUT_DIR + "/concat_sort/{sample_name}/{sample_name}_sorted.bam")
    output:
        os.path.normpath(OUTPUT_DIR + "/concat_sort/{sample_name}/{sample_name}_sorted.bam.bai")
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

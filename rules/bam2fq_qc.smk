"""
################################################################################
These rules make the BAM to FASTQ file conversion and generate Quality Control
################################################################################
"""

"""
This rule concatenates all UBAM per sample
"""
def input_concat_ubam(wildcards):
    indices = [i for i in range(len(SAMPLE_NAME)) if SAMPLE_NAME[i] == wildcards.sample_name]
    if config["input_format"] == "pod5":
        input = expand(os.path.normpath(OUTPUT_DIR + "/tmp/calling/"+ wildcards.sample_name +"/{batch_name}.bam"), batch_name=[BATCH_NAME[i] for i in indices])
    elif config["input_format"] == "ubam":
        input = expand(os.path.normpath(OUTPUT_DIR + "/tmp/ubam/"+ wildcards.sample_name +"/{batch_name}.bam"), batch_name=[BATCH_NAME[i] for i in indices])
    return input

rule concat_ubam:
    input:
        ubams = input_concat_ubam
    output:
        temp(os.path.normpath(OUTPUT_DIR + "/Ubam/{sample_name}/{sample_name}_concat.bam"))
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: min(2048 + 5120 * (attempt - 1),10240),
        time_min = (lambda wildcards, attempt: attempt * 300)
    conda:
        CONDA_ENV_SAMTOOLS
    shell:
        """
       	samtools cat {input.ubams} > {output}
        """

"""
This rule sorts the concatenated UBAM
"""
rule sort_ubam:
    input:
        ubam = os.path.normpath(OUTPUT_DIR + "/Ubam/{sample_name}/{sample_name}_concat.bam")
    output:
        os.path.normpath(OUTPUT_DIR + "/Ubam/{sample_name}/{sample_name}_sorted.bam")
    threads: 12
    resources:
        mem_mb=lambda wildcards, attempt: min(20480 + 5120 * (attempt - 1 ),35840),
        time_min = (lambda wildcards, attempt: attempt * 300)
    conda:
        CONDA_ENV_SAMTOOLS
    shell:
        """
       	samtools sort -@ 12 -o {output} {input.ubam}
        """

"""
This rule indexes the sorted UBAM
"""
rule index_ubam:
    input:
        ubam = os.path.normpath(OUTPUT_DIR + "/Ubam/{sample_name}/{sample_name}_sorted.bam")
    output:
        os.path.normpath(OUTPUT_DIR + "/Ubam/{sample_name}/{sample_name}_sorted.bam.bai")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: min(1024 + 5120 * (attempt - 1 ),10240),
        time_min = (lambda wildcards, attempt: attempt * 300)
    conda:
        CONDA_ENV_SAMTOOLS
    shell:
        """
       	samtools index {input.ubam}
       	"""


"""
This rule converts UBAM file to FASTQ file
"""
rule samtools_ubam_to_fastq:
    input:
        ubam_file = os.path.normpath(OUTPUT_DIR + "/Ubam/{sample_name}/{sample_name}_sorted.bam")
    output:
        fastq_file = os.path.normpath(OUTPUT_DIR + "/Fastq/{sample_name}.fastq.gz")
    threads:
        8
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 24576),
        time_min = (lambda wildcards, attempt: attempt * 360)
    conda:
        CONDA_ENV_SAMTOOLS
    shell:
        """
        samtools fastq --threads {threads} {input.ubam_file} > {OUTPUT_DIR}/Fastq/{wildcards.sample_name}.fastq && \
        gzip -v {OUTPUT_DIR}/Fastq/{wildcards.sample_name}.fastq
        """

"""
This rule generates QC for fastq file
/!\ default java memory option = 512MB, which is not enough for 50G fastq files, hence the use of "--memory 1024" (or more)
"""

rule fastqc:
    input:
        fastq_file = os.path.normpath(OUTPUT_DIR + "/Fastq/{sample_name}.fastq.gz")
    output:
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/fastq_QC/fastqc/{sample_name}/{sample_name}_fastqc.html"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/fastq_QC/fastqc/{sample_name}/{sample_name}_fastqc.zip")
    threads:
        4
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 20480),
        time_min = (lambda wildcards, attempt: attempt * 240)
    shell:
        """
        TMP_DIR=$(mktemp -d -t lr_pipeline-XXXXXXXXXX) && \
        _JAVA_OPTIONS="-Djava.io.tmpdir=${{TMP_DIR}}" && \
        export _JAVA_OPTION && \
        singularity exec --contain -B {OUTPUT_DIR},${{TMP_DIR}} {SING_ENV_FASTQC} \
        fastqc --threads {threads} --memory 4096 --dir ${{TMP_DIR}} --outdir {OUTPUT_DIR}/Quality_Control/fastq_QC/fastqc/{wildcards.sample_name}/ {input.fastq_file}
        """

"""
This rule generates QC for fastq file to check contaminations
"""

rule fastq_screen:
    input:
        fastq_file = os.path.normpath(OUTPUT_DIR + "/Fastq/{sample_name}.fastq.gz"),
        config = config["references"]["fastq_screen_conf"]
    output:
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/fastq_QC/fastq_screen/{sample_name}/{sample_name}_screen.txt"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/fastq_QC/fastq_screen/{sample_name}/{sample_name}_screen.html")
    threads:
        8
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 20480),
        time_min = (lambda wildcards, attempt: attempt * 240)
    params:
        path_config = os.path.dirname(config["references"]["fastq_screen_conf"])
    shell:
        """
        TMP_DIR=$(mktemp -d -t lr_pipeline-XXXXXXXXXX) && \
        singularity exec --contain -B {OUTPUT_DIR},{params.path_config} -B ${{TMP_DIR}}:/tmp {SING_ENV_FASTQ_SCREEN} \
        fastq_screen --aligner minimap2 --threads {threads} --outdir {OUTPUT_DIR}/Quality_Control/fastq_QC/fastq_screen/{wildcards.sample_name}/ --force --conf {input.config} {input.fastq_file}
        """
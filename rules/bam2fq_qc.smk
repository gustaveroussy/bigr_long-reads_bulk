"""
################################################################################
These rules make the BAM to FASTQ file conversion and generate Quality Control
################################################################################
"""

"""
This rule converts BAM file to FASTQ file
"""

rule samtools_bam_to_fastq:
    input:
        bam_file = os.path.normpath(OUTPUT_DIR + "/reconcat/{sample_name}/{sample_name}_sorted.bam")
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
        samtools fastq --threads {threads} {input.bam_file} > {OUTPUT_DIR}/Fastq/{wildcards.sample_name}.fastq && \
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
        singularity exec --contain -B {OUTPUT_DIR} -B ${{TMP_DIR}}:/tmp {SING_ENV_FASTQC} \
        fastqc --threads {threads} --memory 4096 --outdir {OUTPUT_DIR}/Quality_Control/fastq_QC/fastqc/{wildcards.sample_name}/ {input.fastq_file}
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
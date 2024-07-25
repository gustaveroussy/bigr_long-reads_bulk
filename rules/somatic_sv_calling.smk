"""
##########################################################################
These rules make the SV Calling for somatic variants
##########################################################################
"""

wildcard_constraints:
    sample_name = '|'.join([x for x in SAMPLE_NAME]),
    bam_name = '|'.join([x for x in BAM_NAME])

"""
This rule makes the parsing of all the supporting reads of putative somatic SVs by nanomonsv
"""

def nanomonsv_parsing_input_bam(wildcards):
    index = BAM_NAME.index(wildcards.bam_name)
    return SYMLINK_FILES[index]

rule nanomonsv_parsing:
    input:
        bam_file = nanomonsv_parsing_input_bam,
        fa_ref = config["references"]["genome"],
    output:
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{bam_name}.deletion.sorted.bed.gz"),
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{bam_name}.deletion.sorted.bed.gz.tbi"),
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{bam_name}.insertion.sorted.bed.gz"),
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{bam_name}.insertion.sorted.bed.gz.tbi"),
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{bam_name}.rearrangement.sorted.bedpe.gz"),
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{bam_name}.rearrangement.sorted.bedpe.gz.tbi"),
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{bam_name}.bp_info.sorted.bed.gz"),
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{bam_name}.bp_info.sorted.bed.gz.tbi")
    threads:
        10
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 20480),
        time_min = (lambda wildcards, attempt: attempt * 720)
    params:
        output_path = os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{bam_name}")
    conda:
        CONDA_ENV_NANOMONSV
    shell:
        """
        nanomonsv parse --reference_fasta {input.fa_ref} {input.bam_file} {params.output_path}
  
        """

"""
This rule gets the SV result from the parsed supporting reads data by nanomonsv
"""
def nanomonsv_SV_normal_bam(wildcards):
    index = BAM_NAME.index(wildcards.sample_name + "_normal")
    return SYMLINK_FILES[index]

def nanomonsv_SV_tumor_bam(wildcards):
    index = BAM_NAME.index(wildcards.sample_name + "_tumor")
    return SYMLINK_FILES[index]

rule nanomonsv_SV:
    input:
        normal_bam_file = nanomonsv_SV_normal_bam,
        tumor_bam_file = nanomonsv_SV_tumor_bam,
        res_parsing_normal = os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{sample_name}_normal.deletion.sorted.bed.gz"),
        res_parsing_tumor = os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{sample_name}_tumor.deletion.sorted.bed.gz"),
        fa_ref = config["references"]["genome"],
    output:
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{sample_name}.nanomonsv.result.txt"),
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{sample_name}.nanomonsv.result.vcf"),
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{sample_name}.nanomonsv.sbnd.result.txt"),
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{sample_name}.nanomonsv.supporting_read.txt")
    threads:
        10
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 20480),
        time_min = (lambda wildcards, attempt: attempt * 720)
    params:
        normal_path = os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{sample_name}_normal"),
        tumor_path = os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{sample_name}_tumor")
    conda:
        CONDA_ENV_NANOMONSV
    shell:
        """
        res=$(({resources.mem_mb}/1000))
        nanomonsv get \
        {params.tumor_path} {input.tumor_bam_file} {input.fa_ref} \
        --control_prefix {params.normal_path} \
        --control_bam {input.normal_bam_file} \
        --processes {threads} --max_memory_minimap2 $res --qv15 --use_racon \
        --single_bnd && \
        cd {OUTPUT_DIR}/SV_Calling/nanomonsv/{wildcards.sample_name}/ && \
        mv {wildcards.sample_name}_tumor.nanomonsv.result.txt {wildcards.sample_name}.nanomonsv.result.txt && \
        mv {wildcards.sample_name}_tumor.nanomonsv.result.vcf {wildcards.sample_name}.nanomonsv.result.vcf && \
        mv {wildcards.sample_name}_tumor.nanomonsv.sbnd.result.txt {wildcards.sample_name}.nanomonsv.sbnd.result.txt && \
        mv {wildcards.sample_name}_tumor.nanomonsv.supporting_read.txt {wildcards.sample_name}.nanomonsv.supporting_read.txt
        
        """


"""
This rule classifies the long insertions into several mobile element insertions by nanomonsv
"""
"""
Error:
[E::bwa_idx_load_from_disk] fail to locate the index files
bwa mem -h 200 /mnt/beegfs/database/bioinfo/Index_DB/Fasta/Ensembl/GRCh38.109/homo_sapiens.GRCh38.109.fasta /mnt/beegfs/userdata/m_aglave/long-reads-bulk/test/data_output_somatic/SV_Calling/nanomonsv/somatic_test_data_bam/somatic_test_data_bam.nanomonsv.insert_classify.txt.tmp.fasta
Traceback (most recent call last):
  File "/mnt/beegfs/userdata/m_aglave/long-reads-bulk/envs/conda/f46ba2c9/bin/nanomonsv", line 10, in <module>
    sys.exit(main())
  File "/mnt/beegfs/userdata/m_aglave/long-reads-bulk/envs/conda/f46ba2c9/lib/python3.10/site-packages/nanomonsv/__init__.py", line 13, in main
    args.func(args)
  File "/mnt/beegfs/userdata/m_aglave/long-reads-bulk/envs/conda/f46ba2c9/lib/python3.10/site-packages/nanomonsv/run.py", line 571, in insert_classify_main
    subprocess.check_call(["bwa", "mem", "-h", "200", args.reference_fasta, args.output_file + ".tmp.fasta"], stdout = hout)
  File "/mnt/beegfs/userdata/m_aglave/long-reads-bulk/envs/conda/f46ba2c9/lib/python3.10/subprocess.py", line 369, in check_call
    raise CalledProcessError(retcode, cmd)
subprocess.CalledProcessError: Command '['bwa', 'mem', '-h', '200', '/mnt/beegfs/database/bioinfo/Index_DB/Fasta/Ensembl/GRCh38.109/homo_sapiens.GRCh38.109.fasta', '/mnt/beegfs/userdata/m_aglave/long-reads-bulk/test/data_output_somatic/SV_Calling/nanomonsv/somatic_test_data_bam/somatic_test_data_bam.nanomonsv.insert_classify.txt.tmp.fasta']' returned non-zero exit status 1.

"""
rule nanomonsv_classifier:
    input:
        sv_txt_file = os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{sample_name}.nanomonsv.result.txt"),
        fa_ref = config["references"]["genome"]
    output:
        sv_txt_file = os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{sample_name}.nanomonsv.insert_classify.txt")
    threads:
        10
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 20480),
        time_min = (lambda wildcards, attempt: attempt * 720)
    params:
        genome_name = "hg38" if re.search("GRCh38", config["references"]["genome_name"]) else "hg19"
    conda:
        CONDA_ENV_NANOMONSV
    shell:
        """
        nanomonsv insert_classify \
        --genome_id {params.genome_name} \
        {input.sv_txt_file} \
        {output.sv_txt_file} \
        {input.fa_ref}

        """

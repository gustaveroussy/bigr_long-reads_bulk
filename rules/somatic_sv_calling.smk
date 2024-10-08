"""
##########################################################################
These rules make the SV Calling for somatic variants
##########################################################################
"""

"""
This rule makes the parsing of all the supporting reads of putative somatic SVs by nanomonsv
"""

rule nanomonsv_parsing:
    input:
        bam_file = os.path.normpath(OUTPUT_DIR + "/tmp/reconcat/{sample_name}/{sample_name}_sorted.bam"),
        bai_file = os.path.normpath(OUTPUT_DIR + "/tmp/reconcat/{sample_name}/{sample_name}_sorted.bam.bai"),
        fa_ref = config["references"]["genome"],
    output:
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/{sample_name}.deletion.sorted.bed.gz"),
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/{sample_name}.deletion.sorted.bed.gz.tbi"),
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/{sample_name}.insertion.sorted.bed.gz"),
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/{sample_name}.insertion.sorted.bed.gz.tbi"),
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/{sample_name}.rearrangement.sorted.bedpe.gz"),
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/{sample_name}.rearrangement.sorted.bedpe.gz.tbi"),
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/{sample_name}.bp_info.sorted.bed.gz"),
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/{sample_name}.bp_info.sorted.bed.gz.tbi")
    threads:
        10
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 20480),
        time_min = (lambda wildcards, attempt: attempt * 720)
    params:
        output_path = os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/{sample_name}")
    conda:
        CONDA_ENV_NANOMONSV
    shell:
        """
        nanomonsv parse --reference_fasta {input.fa_ref} {input.bam_file} {params.output_path}
  
        """

"""
This rule gets the SV result from the parsed supporting reads data by nanomonsv
"""
def get_res_parsing_normal(wildcards):
    return os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/" + wildcards.pair_somatic + "/" + str(wildcards.pair_somatic).split("_vs_")[0] + ".deletion.sorted.bed.gz")

def get_res_parsing_tumor(wildcards):
    return os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/" + wildcards.pair_somatic + "/" + str(wildcards.pair_somatic).split("_vs_")[1] + ".deletion.sorted.bed.gz")

def get_res_parsing_normal_path(wildcards):
    return os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/" + wildcards.pair_somatic + "/" + str(wildcards.pair_somatic).split("_vs_")[0])

def get_res_parsing_tumor_path(wildcards):
    return os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/" + wildcards.pair_somatic + "/" + str(wildcards.pair_somatic).split("_vs_")[1])

def get_tumor_sample_name(wildcards):
    return str(wildcards.pair_somatic).split("_vs_")[1]
    
rule nanomonsv_SV:
    input:
        normal_bam_file = get_input_normal_bam,
        tumor_bam_file = get_input_tumor_bam,
        normal_bai_file = get_input_normal_bai,
        tumor_bai_file = get_input_tumor_bai,
        res_parsing_normal = get_res_parsing_normal,
        res_parsing_tumor = get_res_parsing_tumor,
        fa_ref = config["references"]["genome"]
    output:
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/{pair_somatic}.nanomonsv.result.txt"),
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/{pair_somatic}.nanomonsv.result.vcf"),
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/{pair_somatic}.nanomonsv.sbnd.result.txt"),
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/{pair_somatic}.nanomonsv.supporting_read.txt")
    threads:
        10
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 20480),
        time_min = (lambda wildcards, attempt: attempt * 720)
    params:
        normal_path = get_res_parsing_normal_path,
        tumor_path = get_res_parsing_tumor_path,
        tumor_sample = get_tumor_sample_name
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
        cd {OUTPUT_DIR}/SV_Calling/Somatic/nanomonsv/{wildcards.pair_somatic}/ && \
        mv {params.tumor_sample}.nanomonsv.result.txt {wildcards.pair_somatic}.nanomonsv.result.txt && \
        mv {params.tumor_sample}.nanomonsv.result.vcf {wildcards.pair_somatic}.nanomonsv.result.vcf && \
        mv {params.tumor_sample}.nanomonsv.sbnd.result.txt {wildcards.pair_somatic}.nanomonsv.sbnd.result.txt && \
        mv {params.tumor_sample}.nanomonsv.supporting_read.txt {wildcards.pair_somatic}.nanomonsv.supporting_read.txt
        
        """


"""
This rule classifies the long insertions into several mobile element insertions by nanomonsv
"""
"""
Error:
[E::bwa_idx_load_from_disk] fail to locate the index files
bwa mem -h 200 /mnt/beegfs/database/bioinfo/Index_DB/Fasta/Ensembl/GRCh38.109/homo_sapiens.GRCh38.109.fasta /mnt/beegfs/userdata/m_aglave/long-reads-bulk/test/data_output_somatic/SV_Calling/Somatic/nanomonsv/somatic_test_data_bam/somatic_test_data_bam.nanomonsv.insert_classify.txt.tmp.fasta
Traceback (most recent call last):
  File "/mnt/beegfs/userdata/m_aglave/long-reads-bulk/envs/conda/f46ba2c9/bin/nanomonsv", line 10, in <module>
    sys.exit(main())
  File "/mnt/beegfs/userdata/m_aglave/long-reads-bulk/envs/conda/f46ba2c9/lib/python3.10/site-packages/nanomonsv/__init__.py", line 13, in main
    args.func(args)
  File "/mnt/beegfs/userdata/m_aglave/long-reads-bulk/envs/conda/f46ba2c9/lib/python3.10/site-packages/nanomonsv/run.py", line 571, in insert_classify_main
    subprocess.check_call(["bwa", "mem", "-h", "200", args.reference_fasta, args.output_file + ".tmp.fasta"], stdout = hout)
  File "/mnt/beegfs/userdata/m_aglave/long-reads-bulk/envs/conda/f46ba2c9/lib/python3.10/subprocess.py", line 369, in check_call
    raise CalledProcessError(retcode, cmd)
subprocess.CalledProcessError: Command '['bwa', 'mem', '-h', '200', '/mnt/beegfs/database/bioinfo/Index_DB/Fasta/Ensembl/GRCh38.109/homo_sapiens.GRCh38.109.fasta', '/mnt/beegfs/userdata/m_aglave/long-reads-bulk/test/data_output_somatic/SV_Calling/Somatic/nanomonsv/somatic_test_data_bam/somatic_test_data_bam.nanomonsv.insert_classify.txt.tmp.fasta']' returned non-zero exit status 1.

#SOLUTION (proposée par Marine mais à tester): faire un index bwa mem "cp /mnt/beegfs/database/bioinfo/Index_DB/Fasta/Ensembl/GRCh38.109/homo_sapiens.GRCh38.109.fasta /mnt/beegfs/database/bioinfo/<PathtoTheFuturBWAIndex>/; bwa index -p /mnt/beegfs/database/bioinfo/<PathtoTheFuturBWAIndex>/homo_sapiens.GRCh38.109 /mnt/beegfs/database/bioinfo/Index_DB/Fasta/Ensembl/GRCh38.109/homo_sapiens.GRCh38.109.fasta"


"""
"""
rule nanomonsv_classifier:
    input:
        sv_txt_file = os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/{pair_somatic}.nanomonsv.result.txt"),
        fa_ref = config["references"]["genome"]
    output:
        sv_txt_file = os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/{pair_somatic}.nanomonsv.insert_classify.txt")
    threads:
        10
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 20480),
        time_min = (lambda wildcards, attempt: attempt * 720)
    params:
        genome_name = config["references"]["genome_nanomonsv_class"]
    conda:
        CONDA_ENV_NANOMONSV
    shell:
        \"""
        nanomonsv insert_classify \
        --genome_id {params.genome_name} \
        {input.sv_txt_file} \
        {output.sv_txt_file} \
        {input.fa_ref}

        \"""
"""
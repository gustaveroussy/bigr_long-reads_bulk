"""
##########################################################################
These rules make the SNV Phasing
##########################################################################
"""
wildcard_constraints:
    sample_name = '|'.join([x for x in SAMPLE_NAME]),
    path_calling_tool_params = "clair3.+|pepper_margin_deepvariant|clairs",
    compl = "_merge_output|_XXX|_snv|_indel"

"""
This rule makes the phasing of SNV by whatshap
"""
def phasing_input_bam(wildcards):
    if config["variant_calling_mode"] == "germline":
        index = BAM_NAME.index(wildcards.sample_name)
        return SYMLINK_FILES[index]
    if config["variant_calling_mode"] == "somatic":
        index_n = BAM_NAME.index(wildcards.sample_name + "_normal")
        index_t = BAM_NAME.index(wildcards.sample_name + "_tumor")
        return [SYMLINK_FILES[index_n], SYMLINK_FILES[index_t]]

def phasing_input_vcf(wildcards):
    return os.path.normpath(OUTPUT_DIR + "/SNV_Calling/" + str(wildcards.path_calling_tool_params) + "/" + str(wildcards.sample_name) + "/" + str(wildcards.sample_name) + str(wildcards.compl) + ".vcf.gz")

rule phasing:
    input:
        bam_file = phasing_input_bam,
        vcf_file = phasing_input_vcf,
        fa_ref = config["references"]["genome"]
    output:
        phased_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{sample_name}/whatshap/{sample_name}{compl}_phased.vcf.gz")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 720)
    conda:
        CONDA_ENV_WHATSHAP
    shell:
        """
        whatshap phase --output {output.phased_vcf_file} --reference {input.fa_ref} --mapping-quality 20 --ignore-read-groups {input.vcf_file} {input.bam_file}
   
        """


"""
This rule makes the vcf index file
"""
def input_vcf_gz(wildcards):
    return os.path.normpath(OUTPUT_DIR + "/SNV_Calling/" + str(wildcards.path_calling_tool_params) + "/" + str(wildcards.sample_name) + "/whatshap/" + str(wildcards.sample_name) + str(wildcards.compl) + "_phased.vcf.gz")

rule tabix_vcf:
    input:
        vcf_gz_file = input_vcf_gz
    output:
        vcf_index = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{sample_name}/whatshap/{sample_name}{compl}_phased.vcf.gz.tbi")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 720)
    conda:
        CONDA_ENV_WHATSHAP
    shell:
        """
        tabix {input.vcf_gz_file}
        """


"""
This rule makes the stat of phasing of SNV by whatshap
"""
def input_vcf_tbi(wildcards):
    return os.path.normpath(OUTPUT_DIR + "/SNV_Calling/" + str(wildcards.path_calling_tool_params) + "/" + str(wildcards.sample_name) + "/whatshap/" + str(wildcards.sample_name) + str(wildcards.compl) + "_phased.vcf.gz.tbi")

rule phasing_stat:
    input:
        vcf_file = input_vcf_gz,
        vcf_tbo_file = input_vcf_tbi
    output:
        phased_stat_txt = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{sample_name}/whatshap/{sample_name}{compl}_phasing_stats.txt"),
        phased_stat_tsv = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{sample_name}/whatshap/{sample_name}{compl}_phasing_stats.tsv"),
        phased_block_tsv = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{sample_name}/whatshap/{sample_name}{compl}_phasing_haplotype_blocks.tsv"),
        phased_block_gtf = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{sample_name}/whatshap/{sample_name}{compl}_phasing_haplotype_blocks.gtf")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 720)
    conda:
        CONDA_ENV_WHATSHAP
    shell:
        """
        whatshap stats --tsv={output.phased_stat_tsv} --block-list={output.phased_block_tsv} --gtf={output.phased_block_gtf} {input.vcf_file} > {output.phased_stat_txt}
   
        """

"""
This rule makes the haplotagging of phasing of SNV by whatshap
"""
rule phasing_haplotagging:
    input:
        bam_file = phasing_input_bam,
        vcf_gz_file = input_vcf_gz,
        vcf_index_file = input_vcf_tbi
    output:
        haplotag_list_tsv = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{sample_name}/whatshap/{sample_name}{compl}_phasing_haplotag_list.tsv"),
        haplotag_bam = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{sample_name}/whatshap/{sample_name}{compl}_haplotagged.bam")
    threads:
        4
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 720)
    conda:
        CONDA_ENV_WHATSHAP
    params:
        fa_ref = config["references"]["genome"]
    shell:
        """
        whatshap haplotag --output-threads={threads} --ignore-read-groups --output-haplotag-list {output.haplotag_list_tsv} --output {output.haplotag_bam} --reference {params.fa_ref} {input.vcf_gz_file} {input.bam_file}

        """


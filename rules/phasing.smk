"""
##########################################################################
These rules make the SNV Phasing
##########################################################################
"""

"""
This rule makes the phasing of SNV by whatshap
"""
def phasing_input_bam(wildcards):
    if config["variant_calling_mode"] == "germline":
        input = os.path.normpath(OUTPUT_DIR + "/tmp/reconcat/" + wildcards.sample_name_or_pair_somatic + "/" + wildcards.sample_name_or_pair_somatic + "_sorted.bam")
        return input
    if config["variant_calling_mode"] == "somatic":
        input_n = os.path.normpath(OUTPUT_DIR + "/tmp/reconcat/" + str(wildcards.sample_name_or_pair_somatic).split("_vs_")[0] + "/" + str(wildcards.sample_name_or_pair_somatic).split("_vs_")[0] + "_sorted.bam")
        input_t = os.path.normpath(OUTPUT_DIR + "/tmp/reconcat/" + str(wildcards.sample_name_or_pair_somatic).split("_vs_")[1] + "/" + str(wildcards.sample_name_or_pair_somatic).split("_vs_")[1] + "_sorted.bam")
        return [input_n, input_t]

rule phasing:
    input:
        bam_file = phasing_input_bam,
        vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/{sample_name_or_pair_somatic}{compl}.vcf.gz"),
        fa_ref = config["references"]["genome"]
    output:
        phased_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/whatshap/{sample_name_or_pair_somatic}{compl}_phased.vcf.gz")
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

rule tabix_vcf:
    input:
        vcf_gz_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/whatshap/{sample_name_or_pair_somatic}{compl}_phased.vcf.gz")
    output:
        vcf_tbi_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/whatshap/{sample_name_or_pair_somatic}{compl}_phased.vcf.gz.tbi")
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

rule phasing_stat:
    input:
        vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/whatshap/{sample_name_or_pair_somatic}{compl}_phased.vcf.gz"),
        vcf_tbi_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/whatshap/{sample_name_or_pair_somatic}{compl}_phased.vcf.gz.tbi")
    output:
        phased_stat_txt = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/whatshap/{sample_name_or_pair_somatic}{compl}_phasing_stats.txt"),
        phased_stat_tsv = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/whatshap/{sample_name_or_pair_somatic}{compl}_phasing_stats.tsv"),
        phased_block_tsv = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/whatshap/{sample_name_or_pair_somatic}{compl}_phasing_haplotype_blocks.tsv"),
        phased_block_gtf = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/whatshap/{sample_name_or_pair_somatic}{compl}_phasing_haplotype_blocks.gtf")
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
This rule haplotags reads (separates reads between Haplotype 1 and Haplotype 2), based on the SNV phasing, by WhatsHap (only for germline data)
"""

rule phasing_haplotagging:
    input:
        bam_file = os.path.normpath(OUTPUT_DIR + "/tmp/reconcat/{sample_name}/{sample_name}_sorted.bam"),
        vcf_gz_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name}/whatshap/{sample_name}{compl}_phased.vcf.gz"),
        vcf_index_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name}/whatshap/{sample_name}{compl}_phased.vcf.gz.tbi")
    output:
        haplotag_list_tsv = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name}/whatshap/{sample_name}{compl}_phasing_haplotag_list.tsv"),
        haplotag_bam = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name}/whatshap/{sample_name}{compl}_haplotagged.bam")
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


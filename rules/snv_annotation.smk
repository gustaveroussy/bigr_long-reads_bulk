"""
##########################################################################
These rules make the SNV Annotation
##########################################################################
"""
wildcard_constraints:
    sample_name = '|'.join([x for x in SAMPLE_NAME]),
    claire3_model = '|'.join([x for x in NAME_CLAIR3_MODEL]),
    path_calling_tool_params = "clair3.+|pepper_margin_deepvariant|clairs",
    filter = '|'.join([x for x in SNPSIFT_FILTERS_NAMES]),
    compl = "_merge_output|_XXX|_snv|_indel"

"""
This rule makes the annotation of SNV by snpEff
"""

def snpeff_annotation_input_vcf(wildcards):
    return os.path.normpath(OUTPUT_DIR + "/SNV_Calling/" + str(wildcards.path_calling_tool_params) + "/" + str(wildcards.sample_name) + "/" + str(wildcards.sample_name) + str(wildcards.compl) + ".vcf.gz")

rule snpeff_annotation:
    input:
        vcf_file = snpeff_annotation_input_vcf
    output:
        annotated_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{sample_name}/snpEff/{sample_name}{compl}_annotated.vcf"),
        snpEff_stat_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{sample_name}/snpEff/{sample_name}{compl}_annotated_stats.csv"),
        snpEff_summary_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{sample_name}/snpEff/{sample_name}{compl}_annotated_summary.html")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 350)
    params:
        genome_name = config["references"]["genome_snpEff"]
    shell:
        """
        TMP_DIR=$(mktemp -d -t lr_pipeline-XXXXXXXXXX) && \
        TMP_DIR2=$(mktemp -d -t lr_pipeline-XXXXXXXXXX) && \
        singularity exec --contain -B {OUTPUT_DIR} -B ${{TMP_DIR2}}:/snpEff/data -B ${{TMP_DIR}}:/tmp {SING_ENV_SNPEFF} \
        java -jar /snpEff/snpEff.jar \
        {params.genome_name} \
        {input.vcf_file} \
        -csvStats {output.snpEff_stat_file} \
        -stats {output.snpEff_summary_file} > {output.annotated_vcf_file}
   
        """


"""
This rule makes the annotation of SNV by snpSift with dbsnp
"""

rule snpsift_annotation_dbsnp:
    input:
        annotated_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{sample_name}/snpEff/{sample_name}{compl}_annotated.vcf"),
        database = config["references"]["dbsnp"]
    output:
        annotated_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{sample_name}/snpEff/{sample_name}{compl}_annotated_dbnsfp.vcf")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 350)
    params:
        database = os.path.dirname(config["references"]["dbsnp"])
    shell:
        """
        TMP_DIR=$(mktemp -d -t lr_pipeline-XXXXXXXXXX) && \
        singularity exec --contain -B {OUTPUT_DIR},{params.database} -B ${{TMP_DIR}}:/tmp {SING_ENV_SNPEFF} \
        java -jar /snpEff/SnpSift.jar \
        dbnsfp -v -db {input.database} \
        {input.annotated_vcf_file} > {output.annotated_vcf_file}

        """

"""
This rule makes the annotation of SNV by snpSift with clinvar
"""

rule snpsift_annotation_clinvar:
    input:
        annotated_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{sample_name}/snpEff/{sample_name}{compl}_annotated_dbnsfp.vcf"),
        database = config["references"]["clinvar"]
    output:
        annotated_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{sample_name}/snpEff/{sample_name}{compl}_annotated_dbnsfp_clinvar.vcf")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 350)
    params:
        database = os.path.dirname(config["references"]["clinvar"])
    shell:
        """
        TMP_DIR=$(mktemp -d -t lr_pipeline-XXXXXXXXXX) && \
        singularity exec --contain -B {OUTPUT_DIR},{params.database} -B ${{TMP_DIR}}:/tmp {SING_ENV_SNPEFF} \
        java -jar /snpEff/SnpSift.jar \
        annotate -v {input.database} \
        {input.annotated_vcf_file} > {output.annotated_vcf_file}

        """


"""
This rule makes the filtering of SNV by snpSift
"""
def snpsift_filter_params(wildcards):
    index = SNPSIFT_FILTERS_NAMES.index(wildcards.filter)
    return SNPSIFT_FILTERS[index]

rule snpsift_filtering:
    input:
        annotated_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{sample_name}/snpEff/{sample_name}{compl}_annotated.vcf")
    output:
        filtered_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{path_calling_tool_params}/{sample_name}/snpEff/{sample_name}{compl}_annotated_{filter}.vcf")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 350)
    params:
        filter = snpsift_filter_params
    shell:
        """
        TMP_DIR=$(mktemp -d -t lr_pipeline-XXXXXXXXXX)
        echo "cat {input.annotated_vcf_file} | \
        java -jar /snpEff/SnpSift.jar \
        filter \\\"{params.filter}\\\" > {output.filtered_vcf_file}" | singularity exec --contain -B {OUTPUT_DIR} -B ${{TMP_DIR}}:/tmp {SING_ENV_SNPEFF} bash

        """

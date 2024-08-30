"""
##########################################################################
These rules make the SNV Annotation
##########################################################################
"""

"""
This rule makes the annotation of SNV by snpEff
"""

rule snpeff_annotation:
    input:
        vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/{sample_name_or_pair_somatic}{compl}.vcf.gz")
    output:
        annotated_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/snpEff/{sample_name_or_pair_somatic}{compl}_annotated.vcf"),
        snpEff_stat_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/snpEff/{sample_name_or_pair_somatic}{compl}_annotated_stats.csv"),
        snpEff_summary_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/snpEff/{sample_name_or_pair_somatic}{compl}_annotated_summary.html")
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
This rule makes the annotation of SNV by snpSift with dbnsfp
"""
def input_snpsift_annotation_dbnsfp(wildcards):
    if SNPEFF_SUFFIX != "":
        return os.path.normpath(OUTPUT_DIR + "/SNV_Calling/" + str(wildcards.variant_calling_mode) + "/" + str(wildcards.path_calling_tool_params) + "/" + str(wildcards.sample_name_or_pair_somatic) + "/snpEff/" + str(wildcards.sample_name_or_pair_somatic) + str(wildcards.compl) + SNPEFF_SUFFIX + ".vcf"),
    else:
        return os.path.normpath(OUTPUT_DIR + "/SNV_Calling/" + str(wildcards.variant_calling_mode) + "/" + str(wildcards.path_calling_tool_params) + "/" + str(wildcards.sample_name_or_pair_somatic) + "/" + str(wildcards.sample_name_or_pair_somatic) + str(wildcards.compl) + ".vcf.gz"),

rule snpsift_annotation_dbnsfp:
    input:
        annotated_vcf_file = input_snpsift_annotation_dbnsfp,
        database = config["references"]["dbnsfp"]
    output:
        annotated_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/snpEff/{sample_name_or_pair_somatic}{compl}" + SNPEFF_SUFFIX + "_dbnsfp.vcf")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 350)
    params:
        database = os.path.dirname(config["references"]["dbnsfp"])
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
def input_snpsift_annotation_clinvar(wildcards):
    if SNPEFF_SUFFIX != "" or DBNSFP_SUFFIX != "":
        return os.path.normpath(OUTPUT_DIR + "/SNV_Calling/" + str(wildcards.variant_calling_mode) + "/" + str(wildcards.path_calling_tool_params) + "/" + str(wildcards.sample_name_or_pair_somatic) + "/snpEff/" + str(wildcards.sample_name_or_pair_somatic) + str(wildcards.compl) + SNPEFF_SUFFIX + DBNSFP_SUFFIX + ".vcf"),
    else:
        return os.path.normpath(OUTPUT_DIR + "/SNV_Calling/" + str(wildcards.variant_calling_mode) + "/" + str(wildcards.path_calling_tool_params) + "/" + str(wildcards.sample_name_or_pair_somatic) + "/" + str(wildcards.sample_name_or_pair_somatic) + str(wildcards.compl) + ".vcf.gz"),

rule snpsift_annotation_clinvar:
    input:
        annotated_vcf_file = input_snpsift_annotation_clinvar,
        database = config["references"]["clinvar"]
    output:
        annotated_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/snpEff/{sample_name_or_pair_somatic}{compl}" + SNPEFF_SUFFIX + DBNSFP_SUFFIX + "_clinvar.vcf")
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

def input_snpsift_filtering(wildcards):
    if SNPEFF_SUFFIX != "" or DBNSFP_SUFFIX != "" or CLINVAR_SUFFIX != "":
        return os.path.normpath(OUTPUT_DIR + "/SNV_Calling/" + str(wildcards.variant_calling_mode) + "/" + str(wildcards.path_calling_tool_params) + "/" + str(wildcards.sample_name_or_pair_somatic) + "/snpEff/" + str(wildcards.sample_name_or_pair_somatic) + str(wildcards.compl) + SNPEFF_SUFFIX + DBNSFP_SUFFIX + CLINVAR_SUFFIX + ".vcf"),
    else:
        return os.path.normpath(OUTPUT_DIR + "/SNV_Calling/" + str(wildcards.variant_calling_mode) + "/" + str(wildcards.path_calling_tool_params) + "/" + str(wildcards.sample_name_or_pair_somatic) + "/" + str(wildcards.sample_name_or_pair_somatic) + str(wildcards.compl) + ".vcf.gz"),

def snpsift_filter_params(wildcards):
    index = SNPSIFT_FILTERS_NAMES.index(wildcards.filter)
    return SNPSIFT_FILTERS[index]

rule snpsift_filtering:
    input:
        annotated_vcf_file = input_snpsift_filtering
    output:
        filtered_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/snpEff/{sample_name_or_pair_somatic}{compl}" + SNPEFF_SUFFIX + DBNSFP_SUFFIX + CLINVAR_SUFFIX + "_{filter}.vcf")
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

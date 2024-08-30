"""
##########################################################################
These rules make the SV Annotation
##########################################################################
"""

"""
This rule makes the annotation of SV by AnnotSV
"""

rule annotsv_annotation:
    input:
        vcf_file = os.path.normpath(OUTPUT_DIR + "/SV_Calling/{variant_calling_mode}/{path_calling_tool_params_SV}/{sample_name_and_all_samples}/{sample_name_and_all_samples}{compl_SV}.vcf")
    output:
        annotated_vcf_file = os.path.normpath(OUTPUT_DIR + "/SV_Calling/{variant_calling_mode}/{path_calling_tool_params_SV}/{sample_name_and_all_samples}/AnnotSV/{sample_name_and_all_samples}{compl_SV}.annotated.vcf"),
        annotated_tsv_file = os.path.normpath(OUTPUT_DIR + "/SV_Calling/{variant_calling_mode}/{path_calling_tool_params_SV}/{sample_name_and_all_samples}/AnnotSV/{sample_name_and_all_samples}{compl_SV}.annotated.tsv")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 4096),
        time_min = (lambda wildcards, attempt: attempt * 60)
    params:
        genome_annotsv = config["references"]["genome_annotsv"],
        annot_path = config["references"]["annotsv"]
    conda:
        CONDA_ENV_ANNOTSV
    shell:
        """
        AnnotSV \
        -genomeBuild {params.genome_annotsv} \
        -annotationsDir {params.annot_path} \
        -SVinputFile {input.vcf_file} \
        -outputDir {OUTPUT_DIR}/SV_Calling/{wildcards.variant_calling_mode}/{wildcards.path_calling_tool_params_SV}/{wildcards.sample_name_and_all_samples}/AnnotSV/ \
        -vcf 1
        """
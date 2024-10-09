"""
##########################################################################
These rules make the SNV Annotation
##########################################################################
"""

if 'snpsift' not in config or 'keep_only_pass' not in config["snpsift"] or config["snpsift"]["keep_only_pass"]: 

    """
    This rule makes the filtering of "PASS" SNV by snpSift
    """
    
    rule snpsift_filtering_PASS:
        input:
            vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/{sample_name_or_pair_somatic}{compl}.vcf.gz")
        output:
            filtered_vcf_file = temp(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/snpEff/{sample_name_or_pair_somatic}{compl}.vcf.gz")),
            index = temp(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/snpEff/{sample_name_or_pair_somatic}{compl}.vcf.gz.tbi"))
        threads:
            3
        resources:
            mem_mb = (lambda wildcards, attempt: attempt * 40960),
            time_min = (lambda wildcards, attempt: attempt * 350)
        params:
            filter = "(FILTER = 'PASS')"
        conda:
            CONDA_ENV_SAMTOOLS
        shell:
            """
            TMP_DIR=$(mktemp -d -t lr_pipeline-XXXXXXXXXX) && \
            singularity exec --contain -B {OUTPUT_DIR} -B ${{TMP_DIR}}:/tmp {SING_ENV_SNPEFF} \
            java -jar /snpEff/SnpSift.jar \
            filter "{params.filter}" {input.vcf_file} > {OUTPUT_DIR}/SNV_Calling/{wildcards.variant_calling_mode}/{wildcards.path_calling_tool_params}/{wildcards.sample_name_or_pair_somatic}/snpEff/{wildcards.sample_name_or_pair_somatic}{wildcards.compl}.vcf
            
            bgzip {OUTPUT_DIR}/SNV_Calling/{wildcards.variant_calling_mode}/{wildcards.path_calling_tool_params}/{wildcards.sample_name_or_pair_somatic}/snpEff/{wildcards.sample_name_or_pair_somatic}{wildcards.compl}.vcf
            tabix {output.filtered_vcf_file}
            """


"""
This rule makes the annotation of SNV ("PASS" or not filtered) by snpEff
"""
def get_input_snpeff(wildcards):
    if 'snpsift' not in config or 'keep_only_pass' not in config["snpsift"] or config["snpsift"]["keep_only_pass"]:
        input = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/snpEff/{sample_name_or_pair_somatic}{compl}.vcf.gz") #Keep only PASS variants
    else:
        input = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/{sample_name_or_pair_somatic}{compl}.vcf.gz") #Keep all variants
    return input

rule snpeff_annotation:
    input:
        vcf_file = get_input_snpeff
    output:
        annotated_vcf_file = temp(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/snpEff/{sample_name_or_pair_somatic}{compl}_tmp_annot1.vcf.gz"))
    threads:
        3
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 350)
    params:
        genome_name = config["references"]["genome_snpEff"]
    conda:
        CONDA_ENV_SAMTOOLS
    shell:
        """
        if [ "{params.genome_name}" == "" ]
        then
            cp {input.vcf_file} {OUTPUT_DIR}/SNV_Calling/{wildcards.variant_calling_mode}/{wildcards.path_calling_tool_params}/{wildcards.sample_name_or_pair_somatic}/snpEff/{wildcards.sample_name_or_pair_somatic}{wildcards.compl}_tmp_annot1.vcf
        else
            TMP_DIR=$(mktemp -d -t lr_pipeline-XXXXXXXXXX) && \
            TMP_DIR2=$(mktemp -d -t lr_pipeline-XXXXXXXXXX) && \
            TMP_DIR3=$(mktemp -d -t lr_pipeline-XXXXXXXXXX) && \
            _JAVA_OPTIONS="-Djava.io.tmpdir=${{TMP_DIR3}}" && \
            export _JAVA_OPTION && \
            singularity exec --contain -B {OUTPUT_DIR},${{TMP_DIR3}} -B ${{TMP_DIR}}:/tmp -B ${{TMP_DIR2}}:/snpEff/data {SING_ENV_SNPEFF} \
            java -jar /snpEff/snpEff.jar \
            {params.genome_name} \
            {input.vcf_file} \
            -csvStats {OUTPUT_DIR}/SNV_Calling/{wildcards.variant_calling_mode}/{wildcards.path_calling_tool_params}/{wildcards.sample_name_or_pair_somatic}/snpEff/{wildcards.sample_name_or_pair_somatic}{wildcards.compl}_annotated_stats.csv \
            -stats {OUTPUT_DIR}/SNV_Calling/{wildcards.variant_calling_mode}/{wildcards.path_calling_tool_params}/{wildcards.sample_name_or_pair_somatic}/snpEff/{wildcards.sample_name_or_pair_somatic}{wildcards.compl}_annotated_summary.html > {OUTPUT_DIR}/SNV_Calling/{wildcards.variant_calling_mode}/{wildcards.path_calling_tool_params}/{wildcards.sample_name_or_pair_somatic}/snpEff/{wildcards.sample_name_or_pair_somatic}{wildcards.compl}_tmp_annot1.vcf
        fi
        
        bgzip {OUTPUT_DIR}/SNV_Calling/{wildcards.variant_calling_mode}/{wildcards.path_calling_tool_params}/{wildcards.sample_name_or_pair_somatic}/snpEff/{wildcards.sample_name_or_pair_somatic}{wildcards.compl}_tmp_annot1.vcf
        tabix {output.annotated_vcf_file}
        
        """

"""
This rule makes the annotation of SNV by snpSift with dbnsfp
"""

rule snpsift_annotation_dbnsfp:
    input:
        annotated_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/snpEff/{sample_name_or_pair_somatic}{compl}_tmp_annot1.vcf.gz")
    output:
        annotated_vcf_file = temp(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/snpEff/{sample_name_or_pair_somatic}{compl}_tmp_annot2.vcf.gz"))
    threads:
        3
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 350)
    params:
        database = config["references"]["dbnsfp"],
        database_path = os.path.dirname(config["references"]["dbnsfp"])
    conda:
        CONDA_ENV_SAMTOOLS
    shell:
        """
        if [ "{params.database_path}" == "" ]
        then
            cp {input.annotated_vcf_file} {OUTPUT_DIR}/SNV_Calling/{wildcards.variant_calling_mode}/{wildcards.path_calling_tool_params}/{wildcards.sample_name_or_pair_somatic}/snpEff/{wildcards.sample_name_or_pair_somatic}{wildcards.compl}_tmp_annot2.vcf
        else
        TMP_DIR=$(mktemp -d -t lr_pipeline-XXXXXXXXXX) && \
            singularity exec --contain -B {OUTPUT_DIR},{params.database_path} -B ${{TMP_DIR}}:/tmp {SING_ENV_SNPEFF} \
            java -jar /snpEff/SnpSift.jar \
            dbnsfp -v -db {params.database} \
            {input.annotated_vcf_file} > {OUTPUT_DIR}/SNV_Calling/{wildcards.variant_calling_mode}/{wildcards.path_calling_tool_params}/{wildcards.sample_name_or_pair_somatic}/snpEff/{wildcards.sample_name_or_pair_somatic}{wildcards.compl}_tmp_annot2.vcf
        fi
        
        bgzip {OUTPUT_DIR}/SNV_Calling/{wildcards.variant_calling_mode}/{wildcards.path_calling_tool_params}/{wildcards.sample_name_or_pair_somatic}/snpEff/{wildcards.sample_name_or_pair_somatic}{wildcards.compl}_tmp_annot2.vcf
        tabix {output.annotated_vcf_file}
        """

"""
This rule makes the annotation of SNV by snpSift with clinvar
"""

rule snpsift_annotation_clinvar:
    input:
        annotated_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/snpEff/{sample_name_or_pair_somatic}{compl}_tmp_annot2.vcf.gz")
    output:
        annotated_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/snpEff/{sample_name_or_pair_somatic}{compl}_annotated.vcf.gz")
    threads:
        3
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 350)
    params:
        database = config["references"]["clinvar"],
        database_path = os.path.dirname(config["references"]["clinvar"])
    conda:
        CONDA_ENV_SAMTOOLS
    shell:
        """
        if [ "{params.database_path}" == ""]
        then
            cp {input.annotated_vcf_file} {OUTPUT_DIR}/SNV_Calling/{wildcards.variant_calling_mode}/{wildcards.path_calling_tool_params}/{wildcards.sample_name_or_pair_somatic}/snpEff/{wildcards.sample_name_or_pair_somatic}{wildcards.compl}_annotated.vcf
        else
            TMP_DIR=$(mktemp -d -t lr_pipeline-XXXXXXXXXX) && \
            singularity exec --contain -B {OUTPUT_DIR},{params.database_path} -B ${{TMP_DIR}}:/tmp {SING_ENV_SNPEFF} \
            java -jar /snpEff/SnpSift.jar \
            annotate -v {params.database} \
            {input.annotated_vcf_file} > {OUTPUT_DIR}/SNV_Calling/{wildcards.variant_calling_mode}/{wildcards.path_calling_tool_params}/{wildcards.sample_name_or_pair_somatic}/snpEff/{wildcards.sample_name_or_pair_somatic}{wildcards.compl}_annotated.vcf
        fi
        
        bgzip {OUTPUT_DIR}/SNV_Calling/{wildcards.variant_calling_mode}/{wildcards.path_calling_tool_params}/{wildcards.sample_name_or_pair_somatic}/snpEff/{wildcards.sample_name_or_pair_somatic}{wildcards.compl}_annotated.vcf
        tabix {output.annotated_vcf_file}
        """

"""
This rule makes the filtering of SNV by snpSift
"""

def snpsift_filter_params(wildcards):
    index = SNPSIFT_FILTERS_NAMES.index(wildcards.filter)
    filter = SNPSIFT_FILTERS[index]
    filter_format = filter.replace('SAMPLE', str(wildcards.sample_name_or_pair_somatic))
    return filter_format


rule snpsift_filtering:
    input:
        annotated_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/snpEff/{sample_name_or_pair_somatic}{compl}_annotated.vcf.gz")
    output:
        filtered_vcf_file = temp(os.path.normpath(OUTPUT_DIR + "/tmp/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/snpEff/{sample_name_or_pair_somatic}{compl}_annotated_{filter}.vcf")),
        filtered_vcf_file_gz = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name_or_pair_somatic}/snpEff/{sample_name_or_pair_somatic}{compl}_annotated_{filter}.vcf.gz"),
    threads:
        3
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 350)
    params:
        filter = snpsift_filter_params
    conda:
        CONDA_ENV_SAMTOOLS
    shell:
        """
        TMP_DIR=$(mktemp -d -t lr_pipeline-XXXXXXXXXX) && \
        singularity exec --contain -B {OUTPUT_DIR} -B ${{TMP_DIR}}:/tmp {SING_ENV_SNPEFF} \
        java -jar /snpEff/SnpSift.jar \
        filter "{params.filter}" {input.annotated_vcf_file} > {OUTPUT_DIR}/SNV_Calling/{wildcards.variant_calling_mode}/{wildcards.path_calling_tool_params}/{wildcards.sample_name_or_pair_somatic}/snpEff/{wildcards.sample_name_or_pair_somatic}{wildcards.compl}_annotated_{wildcards.filter}.vcf

        cp {OUTPUT_DIR}/SNV_Calling/{wildcards.variant_calling_mode}/{wildcards.path_calling_tool_params}/{wildcards.sample_name_or_pair_somatic}/snpEff/{wildcards.sample_name_or_pair_somatic}{wildcards.compl}_annotated_{wildcards.filter}.vcf {output.filtered_vcf_file}
        bgzip {OUTPUT_DIR}/SNV_Calling/{wildcards.variant_calling_mode}/{wildcards.path_calling_tool_params}/{wildcards.sample_name_or_pair_somatic}/snpEff/{wildcards.sample_name_or_pair_somatic}{wildcards.compl}_annotated_{wildcards.filter}.vcf
        tabix {output.filtered_vcf_file_gz}
        """

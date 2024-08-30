"""
##########################################################################
These rules make the SV Calling for germline variants
##########################################################################
"""

"""
This rule makes the SV Calling by sniffles
"""

rule sniffles:
    input:
        bam_file = os.path.normpath(OUTPUT_DIR + "/reconcat/{sample_name}/{sample_name}_sorted.bam"),
        fa_ref = config["references"]["genome"]
    output:
        vcf_file = os.path.normpath(OUTPUT_DIR + "/SV_Calling/Germline/sniffles/{sample_name}/{sample_name}_SV.vcf"),
        snf_file = os.path.normpath(OUTPUT_DIR + "/SV_Calling/Germline/sniffles/{sample_name}/{sample_name}_SV.snf")
    threads:
        4
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 20480),
        time_min = (lambda wildcards, attempt: attempt * 720)
    conda:
        CONDA_ENV_SNIFFLES
    shell:
        """
        sniffles -i {input.bam_file} --reference {input.fa_ref} -v {output.vcf_file} --snf {output.snf_file}

        """

"""
This rule merges the SV Calling by sniffles for all samples
"""

rule merge_sniffles:
    input:
        snf_files = expand(os.path.normpath(OUTPUT_DIR + + "/SV_Calling/Germline/sniffles/" + {sample_name} + "/" + {sample_name} + "_SV.snf"), sample_name = SAMPLE_NAME)
    output:
        vcf_file = os.path.normpath(OUTPUT_DIR + "/SV_Calling/Germline/sniffles/all_samples/all_samples_SV.vcf")
    threads:
        4
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 20480),
        time_min = (lambda wildcards, attempt: attempt * 720)
    conda:
        CONDA_ENV_SNIFFLES
    shell:
        """
        sniffles --input {input.snf_files} --vcf {output.vcf_file} && \
        line_number=$(grep -n "CHROM" {output.vcf_file} | cut -f1 -d":") && \
        sed -i "${{line_number}}s/_SV//g" {output.vcf_file} #Rename samples by removing the _SV suffix
        """
        
"""
This rule makes the SV Calling by cuteSV
"""

rule cuteSV:
    input:
        bam_file = os.path.normpath(OUTPUT_DIR + "/reconcat/{sample_name}/{sample_name}_sorted.bam"),
        fa_ref = config["references"]["genome"]
    output:
        vcf_file = os.path.normpath(OUTPUT_DIR + "/SV_Calling/Germline/cuteSV/{sample_name}/{sample_name}_SV.vcf")
    threads:
        2
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 20480),
        time_min = (lambda wildcards, attempt: attempt * 720)
    conda:
        CONDA_ENV_CUTESV
    shell:
        """
        cuteSV {input.bam_file} {input.fa_ref} {output.vcf_file} {OUTPUT_DIR}/SV_Calling/Germline/cuteSV/{wildcards.sample_name}/ \
        --max_cluster_bias_INS 100 \
        --diff_ratio_merging_INS 0.3 \
        --max_cluster_bias_DEL 100 \
        --diff_ratio_merging_DEL 0.3
        """
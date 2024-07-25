"""
##########################################################################
These rules make the SV Calling for germline variants
##########################################################################
"""
wildcard_constraints:
    bam_name = '|'.join([x for x in BAM_NAME]),

"""
This rule makes the SV Calling by sniffles
"""

def sniffles_input_bam(wildcards):
    index = BAM_NAME.index(wildcards.bam_name)
    return SYMLINK_FILES[index]

rule sniffles:
    input:
        bam_file = sniffles_input_bam,
        fa_ref = config["references"]["genome"]
    output:
        vcf_file = os.path.normpath(OUTPUT_DIR + "/SV_Calling/sniffles/{bam_name}/{bam_name}_SV.vcf"),
        snf_file = os.path.normpath(OUTPUT_DIR + "/SV_Calling/sniffles/{bam_name}/{bam_name}_SV.snf")
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

def merge_sniffles_input_snf(wildcards):
    return [os.path.normpath(OUTPUT_DIR + "/SV_Calling/sniffles/" + x + "/" + x + "_SV.snf") for x in BAM_NAME]
    
rule merge_sniffles:
    input:
        snf_files = merge_sniffles_input_snf
    output:
        vcf_file = os.path.normpath(OUTPUT_DIR + "/SV_Calling/sniffles/all_samples_SV.vcf")
    threads:
        4
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 20480),
        time_min = (lambda wildcards, attempt: attempt * 720)
    conda:
        CONDA_ENV_SNIFFLES
    shell:
        """
        sniffles --input {input.snf_files} --vcf {output.vcf_file}

        """
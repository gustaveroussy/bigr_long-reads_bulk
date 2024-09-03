"""
##########################################################################
These rules make the SNV Calling for germline variants
##########################################################################
"""

"""
This rule makes the mosdepth coverage calculation as recommanded by Spectre for CNV calling
"""

rule cnv_mosdepth:
    input:
        bam_file = os.path.normpath(OUTPUT_DIR + "/reconcat/{sample_name}/{sample_name}_sorted.bam")
    output:
        os.path.normpath(OUTPUT_DIR + "/CNV_Calling/mosdepth/{sample_name}/{sample_name}_Q20.regions.bed.gz"),
        os.path.normpath(OUTPUT_DIR + "/CNV_Calling/mosdepth/{sample_name}/{sample_name}_Q20.regions.bed.gz.csi"),
        os.path.normpath(OUTPUT_DIR + "/CNV_Calling/mosdepth/{sample_name}/{sample_name}_Q20.mosdepth.global.dist.txt"),
        os.path.normpath(OUTPUT_DIR + "/CNV_Calling/mosdepth/{sample_name}/{sample_name}_Q20.mosdepth.region.dist.txt"),
        os.path.normpath(OUTPUT_DIR + "/CNV_Calling/mosdepth/{sample_name}/{sample_name}_Q20.mosdepth.summary.txt")
    threads:
        2
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 20480),
        time_min = (lambda wildcards, attempt: attempt * 120)
    conda:
        CONDA_ENV_MOSDEPTH
    shell:
        """
        cd {OUTPUT_DIR}/CNV_Calling/mosdepth/{wildcards.sample_name}/
        # Parameters recommanded for Spectre usage
        mosdepth -t {threads} --fast-mode --by 1000 --mapq 20 {wildcards.sample_name}_Q20 {input.bam_file}
        """

"""
This rule makes the CNV Calling by Spectre
"""
def input_cancer(wildcards):
    index = CNV_SAMPLE_NAME.index(wildcards.sample_name)
    if CNV_CANCER_BOOL[index]:
        return "--cancer"
    else:
        return ""

rule spectre_cnv:
    input:
        regions = os.path.normpath(OUTPUT_DIR + "/CNV_Calling/mosdepth/{sample_name}/{sample_name}_Q20.regions.bed.gz"),
        fa_ref = config["references"]["genome"]
    output:
        vcf_file = os.path.normpath(OUTPUT_DIR + "/CNV_Calling/spectre/{sample_name}/{sample_name}.vcf"),
        bed_file = os.path.normpath(OUTPUT_DIR + "/CNV_Calling/spectre/{sample_name}/{sample_name}_cnv.bed.gz"),
        bed_tbi = os.path.normpath(OUTPUT_DIR + "/CNV_Calling/spectre/{sample_name}/{sample_name}_cnv.bed.gz.tbi"),
        img = expand(os.path.normpath(OUTPUT_DIR + "/CNV_Calling/spectre/{{sample_name}}/img/{{sample_name}}_plot_cnv_chr_{chromos}.png"), chromos = [x for x in CHR_NUMBER if x not in ["MT","M"]])
    threads:
        2
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 4096),
        time_min = (lambda wildcards, attempt: attempt * 60)
    params:
        chromos_list = ",".join([x for x in CHR_NUMBER if x not in ["MT","M"]]),
        cancer = input_cancer
    conda:
        CONDA_ENV_SPECTRE
    shell:
        """
        spectre CNVCaller \
        --coverage {input.regions} \
        --sample-id {wildcards.sample_name} \
        --output-dir {OUTPUT_DIR}/CNV_Calling/spectre/{wildcards.sample_name}/ \
        --reference {input.fa_ref} {EXTRA_PARAMS_SPECTRE} {params.cancer} \
        --only-chr {params.chromos_list}
        
        cd {OUTPUT_DIR}/CNV_Calling/spectre/{wildcards.sample_name}/img/
        for file in *.png
        do
        new_name=$(echo ${{file}} | sed 's/plot_cnv/plot_cnv_chr/')
        mv ${{file}} ${{new_name}}
        done
        
        """

    
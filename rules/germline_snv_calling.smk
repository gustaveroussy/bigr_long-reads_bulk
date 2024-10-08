"""
##########################################################################
These rules make the SNV Calling for germline variants
##########################################################################
"""


"""
This rule makes the SNV Calling by clair3 with various models
"""

def clair3_input_model_path(wildcards):
    index = NAME_CLAIR3_MODEL.index(wildcards.clair3_model)
    return config["clair3"]["model"][index]

rule clair3:
    input:
        bam_file = os.path.normpath(OUTPUT_DIR + "/tmp/reconcat/{sample_name}/{sample_name}_sorted.bam"),
        index = os.path.normpath(OUTPUT_DIR + "/tmp/reconcat/{sample_name}/{sample_name}_sorted.bam.bai"),
        fa_ref = config["references"]["genome"],
        clair3_path = clair3_input_model_path
    output:
        vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/clair3/{clair3_model}/{sample_name}/{sample_name}_merge_output.vcf.gz"),
        vcf_tbi = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/clair3/{clair3_model}/{sample_name}/{sample_name}_merge_output.vcf.gz.tbi"),
        full_al_vcf_file = temp(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/clair3/{clair3_model}/{sample_name}/full_alignment.vcf.gz")),
        full_al_vcf_tbi = temp(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/clair3/{clair3_model}/{sample_name}/full_alignment.vcf.gz.tbi")),
        pil_vcf_file = temp(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/clair3/{clair3_model}/{sample_name}/pileup.vcf.gz")),
        pil_vcf_tbi = temp(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/clair3/{clair3_model}/{sample_name}/pileup.vcf.gz.tbi"))
    threads:
        10
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 51200),
        time_min = (lambda wildcards, attempt: attempt * 720)
    conda:
        CONDA_ENV_CLAIR3
    shell:
        """
        run_clair3.sh --bam_fn={input.bam_file} --ref_fn={input.fa_ref} --threads={threads} --platform="ont" --model_path={input.clair3_path} --output={OUTPUT_DIR}/SNV_Calling/Germline/clair3/{wildcards.clair3_model}/{wildcards.sample_name}/ --include_all_ctgs && \
        mv {OUTPUT_DIR}/SNV_Calling/Germline/clair3/{wildcards.clair3_model}/{wildcards.sample_name}/merge_output.vcf.gz {output.vcf_file} && \
        mv {OUTPUT_DIR}/SNV_Calling/Germline/clair3/{wildcards.clair3_model}/{wildcards.sample_name}/merge_output.vcf.gz.tbi {output.vcf_tbi}

        """


"""
This rule removes duplicate alignments with same read IDs for pepper_margin_deepvariant
The other solution could be adding a suffix to keep every primary alignments: https://github.com/kishwarshafin/pepper/issues/67
"""

rule sambamba_markdup:
    input:
        bam_file = os.path.normpath(OUTPUT_DIR + "/tmp/reconcat/{sample_name}/{sample_name}_sorted.bam"),
        index = os.path.normpath(OUTPUT_DIR + "/tmp/reconcat/{sample_name}/{sample_name}_sorted.bam.bai")
    output:
        markdup_bam = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/pepper_margin_deepvariant/{sample_name}/markdup/{sample_name}.bam")
    threads:
        4
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 4096),
        time_min = (lambda wildcards, attempt: attempt * 720)
    conda:
        CONDA_ENV_SAMBAMBA
    shell:
        """
        sambamba markdup -r -t {threads} {input.bam_file} {output.markdup_bam}
        """


"""
This rule makes the SNV Calling by pepper_margin_deepvariant
"""

rule pepper_margin_deepvariant:
    input:
        bam_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/pepper_margin_deepvariant/{sample_name}/markdup/{sample_name}.bam"),
        index = os.path.normpath(OUTPUT_DIR + "/tmp/reconcat/{sample_name}/{sample_name}_sorted.bam.bai"),
        fa_ref = config["references"]["genome"]
    output:
        os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/pepper_margin_deepvariant/{sample_name}/{sample_name}.vcf.gz"),
        os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/pepper_margin_deepvariant/{sample_name}/{sample_name}.visual_report.html")
    threads:
        10
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 1440)
    params:
        path_fa_ref = os.path.dirname(config["references"]["genome"])
    shell:
        """
        temp=$(mktemp -d -t lr_pipeline-XXXXXXXXXX) && \
        singularity exec --contain -B {OUTPUT_DIR},{params.path_fa_ref} -B ${{temp}}:${{TEMPDIR}} {SING_ENV_PEPPER_DEEPVARIANT} run_pepper_margin_deepvariant call_variant \
        --bam {input.bam_file} \
        --fasta {input.fa_ref} \
        --output_dir {OUTPUT_DIR}/SNV_Calling/Germline/pepper_margin_deepvariant/{wildcards.sample_name}/ \
        --output_prefix {wildcards.sample_name} \
        --threads {threads} \
        --sample_name {wildcards.sample_name} \
        --ont_r10_q20
 
        """

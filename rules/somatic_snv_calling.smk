"""
##########################################################################
These rules make the SNV Calling for somatic variants 
##########################################################################
"""

"""
This rule makes the SNV Calling by clairS with various models
"""

rule clairs:
    input:
        normal_bam_file = get_input_normal_bam,
        tumor_bam_file = get_input_tumor_bam,
        fa_ref = config["references"]["genome"],
    output:
        snv_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Somatic/clairs/{pair_somatic}/{pair_somatic}_snv.vcf.gz"),
        indel_vcf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Somatic/clairs/{pair_somatic}/{pair_somatic}_indel.vcf.gz")
    threads:
        10
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 20480),
        time_min = (lambda wildcards, attempt: attempt * 720)
    params:
        path_fa_ref = os.path.dirname(config["references"]["genome"]),
        model = config["clairs"]["model"]
    shell:
        """
        TMP_DIR=$(mktemp -d -t lr_pipeline-XXXXXXXXXX) && \
        singularity exec --contain \
          -B {OUTPUT_DIR},{params.path_fa_ref} -B ${{TMP_DIR}}:${{TMPDIR}} \
          {SING_ENV_CLAIRS} \
          /opt/bin/run_clairs \
          --tumor_bam_fn {input.tumor_bam_file} \
          --normal_bam_fn {input.normal_bam_file} \
          --ref_fn {input.fa_ref} \
          --threads {threads} \
          --platform {params.model} \
          --output_dir {OUTPUT_DIR}/SNV_Calling/Somatic/clairs/{wildcards.pair_somatic} \
          --output_prefix {wildcards.pair_somatic}_snv \
          --indel_output_prefix {wildcards.pair_somatic}_indel \
          --sample_name {wildcards.pair_somatic} \
          --include_all_ctgs \
          --remove_intermediate_dir \
          --conda_prefix /opt/conda/envs/clairs \
          --enable_indel_calling
  
        """

 
"""
##########################################################################
These rules make the SNV graphs
##########################################################################
"""

"""
This rule converts VCF files to MAF files
"""

rule vcf_to_maf:
    input:
        vcf_file = os.path.normpath(OUTPUT_DIR + "/tmp/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name}/snpEff/{sample_name}{compl}_annotated_{filter}.vcf"),
        fa_ref = config["references"]["genome"]
    output:
        maf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name}/maftools/{sample_name}{compl}_annotated_{filter}.maf")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 4096),
        time_min = (lambda wildcards, attempt: attempt * 60)
    params:
        path_fa_ref = os.path.dirname(config["references"]["genome"])
    shell:
        """
        TMP_DIR=$(mktemp -d -t lr_pipeline-XXXXXXXXXX) && \
        singularity run --no-home -B {OUTPUT_DIR},{params.path_fa_ref} -B ${{TMP_DIR}}:/tmp {SING_ENV_VCF2MAF} \
        perl /vcf2maf-1.6.22/vcf2maf.pl \
        --input-vcf {input.vcf_file} \
        --output-maf {output.maf_file} \
        --inhibit-vep \
        --ref-fasta {input.fa_ref}
        """

"""
These rule make graphs from MAF files using MAFTOOLS
"""

rule maftools_graphs:
    input:
        maf_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name}/maftools/{sample_name}{compl}_annotated_{filter}.maf")
    output:
        flag_file = os.path.normpath(OUTPUT_DIR + "/SNV_Calling/{variant_calling_mode}/{path_calling_tool_params}/{sample_name}/maftools/{sample_name}{compl}_annotated_{filter}_maftools_graphs_DONE.txt")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 4096),
        time_min = (lambda wildcards, attempt: attempt * 60)
    conda:
        CONDA_ENV_MAFTOOLS
    params: 
        genes_file = GENES_FILE,
        variantType = config["variant_calling_mode"]
    shell:
        """
        Rscript {PIPELINE_DIR}/script/maftools_graphs.R --sampleName {wildcards.sample_name}{wildcards.compl}_annotated_{wildcards.filter} --inputMAF {input.maf_file} --outputDir {OUTPUT_DIR}/SNV_Calling/{wildcards.variant_calling_mode}/{wildcards.path_calling_tool_params}/{wildcards.sample_name}/maftools/ --genesFile {params.genes_file} --variant_type {params.variantType}
        """

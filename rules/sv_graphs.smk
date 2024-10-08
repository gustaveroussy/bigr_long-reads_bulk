"""
##########################################################################
These rules make the SV graphs
##########################################################################
"""

"""
This rule makes some graphs of SV with sniffles2_plot for a single vcf file
"""

rule sniffles2_plot:
    input:
        vcf_file = os.path.normpath(OUTPUT_DIR + "/SV_Calling/{variant_calling_mode}/{path_calling_tool_params_SV}/{sample_name}/{sample_name}_SV.vcf")
    output:
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/{variant_calling_mode}/{path_calling_tool_params_SV}/{sample_name}/sniffles2_plot/{sample_name}.png"),
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/{variant_calling_mode}/{path_calling_tool_params_SV}/{sample_name}/sniffles2_plot/del_ins_genotype.jpg"),
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/{variant_calling_mode}/{path_calling_tool_params_SV}/{sample_name}/sniffles2_plot/del_ins_type_size.jpg"),
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/{variant_calling_mode}/{path_calling_tool_params_SV}/{sample_name}/sniffles2_plot/dup_inv_type_size.jpg"),
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/{variant_calling_mode}/{path_calling_tool_params_SV}/{sample_name}/sniffles2_plot/inv_dup_genotype.jpg"),
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/{variant_calling_mode}/{path_calling_tool_params_SV}/{sample_name}/sniffles2_plot/length_variant.jpg"),
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/{variant_calling_mode}/{path_calling_tool_params_SV}/{sample_name}/sniffles2_plot/variant_count.jpg")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 4096),
        time_min = (lambda wildcards, attempt: attempt * 60)
    conda:
        CONDA_ENV_SNIFFLES
    shell:
        """
        python3 -m sniffles2_plot -i {input.vcf_file} -o {OUTPUT_DIR}/SV_Calling/{wildcards.variant_calling_mode}/{wildcards.path_calling_tool_params_SV}/{wildcards.sample_name}/sniffles2_plot/ && \
        mv {OUTPUT_DIR}/SV_Calling/{wildcards.variant_calling_mode}/{wildcards.path_calling_tool_params_SV}/{wildcards.sample_name}/sniffles2_plot/*.png {OUTPUT_DIR}/SV_Calling/{wildcards.variant_calling_mode}/{wildcards.path_calling_tool_params_SV}/{wildcards.sample_name}/sniffles2_plot/{wildcards.sample_name}.png #Remove new line character introduced by Sniffles2_plot at the end of the sample name
        """

"""
This rule makes some graphs of SV with sniffles2_plot for a multiple vcf files inside a given directory
! Not in the pipeline output for now
! does not seem to work for sv_site_per_genome.jpg and sample_upset.jpg for some samples
"""

rule sniffles2_plot_multisamples:
    input:
        vcf_multisamples_file = os.path.normpath(OUTPUT_DIR + "/SV_Calling/{variant_calling_mode}/{path_calling_tool_params_SV}/all_samples/all_samples_SV.vcf")
    output:
        os.path.normpath(OUTPUT_DIR + "/SV_Calling/{variant_calling_mode}/sniffles/all_samples/sniffles2_plot/heatmap.jpg")
        #os.path.normpath(OUTPUT_DIR + "/SV_Calling/{variant_calling_mode}/sniffles/all_samples/sniffles2_plot/sample_upset.jpg"),
        #os.path.normpath(OUTPUT_DIR + "/SV_Calling/{variant_calling_mode}/sniffles/all_samples/sniffles2_plot/sv_site_per_genome.jpg"),
        #os.path.normpath(OUTPUT_DIR + "/SV_Calling/{variant_calling_mode}/sniffles/all_samples/sniffles2_plot/variant_count.jpg")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 4096),
        time_min = (lambda wildcards, attempt: attempt * 60)
    conda:
        CONDA_ENV_SNIFFLES
    shell:
        """
        python3 -m sniffles2_plot -i {input.vcf_multisamples_file} -o {OUTPUT_DIR}/SV_Calling/{wildcards.variant_calling_mode}/sniffles/all_samples/sniffles2_plot/ && \
        for filename in {OUTPUT_DIR}/SV_Calling/{wildcards.variant_calling_mode}/sniffles/all_samples/sniffles2_plot/*$'\\n'.png; do mv "${{filename}}" "${{filename//[$'\\n']/}}"; done #Remove new line character introduced by Sniffles2_plot at the end of the last sample name of the VCF multisample
        """
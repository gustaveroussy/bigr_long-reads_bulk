"""
This rule launches methylation QC
"""


rule meth_QC_ratio_fwd_rev:
    input:
        uncomb_bed = os.path.normpath(OUTPUT_DIR + "/tmp/bed_uncombined_strands/{sample_name}/{meth_type}/{sample_name}_chr{chr_number}_{meth_type}_uncomb.bed")
    output:
        per_pos_per_strand_gam_ratio = os.path.normpath(OUTPUT_DIR + "/Quality_Control/methylation_QC/per_pos_per_strand_gam_ratio_{sample_name}_chr{chr_number}_{meth_type}_mqc.png")
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: min(5120 + 5120 * (attempt - 1 ),20480),
        time_min = (lambda wildcards, attempt: attempt * 300)
    params:
        script = os.path.normpath(PIPELINE_DIR + "/script")
    shell:
        """
        #module load r/3.6.1
        singularity exec --contain -B {params.script},{OUTPUT_DIR} {SING_ENV_GENEDMR} \
        Rscript {params.script}/QC_methylation_ratio_fwd_rev.R --input_bed {input.uncomb_bed} --output_path {OUTPUT_DIR}/Quality_Control/methylation_QC/
        """


"""
This rule generate motif CG bed in reference
"""


rule motif_cg_gref:
    input:
        ref_fa = config["references"]["genome"],
        ref_path = os.path.normpath(config["references"]["genome"])
    output:
        cg_motif = os.path.normpath(OUTPUT_DIR + "/tmp/resources/motif_cg.bed")
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt,
        time_min = lambda wildcards, attempt: attempt * 300
    shell:
        """
        TMP_DIR=$(mktemp -d -t lr_pipeline-XXXXXXXXXX) && \
        singularity exec --contain -B {OUTPUT_DIR},{input.ref_path} -B ${{TMP_DIR}}:/tmp \
        {SING_ENV_MODKIT} modkit motif-bed {input.ref_fa} CG 0 1> {output.cg_motif}
        """


"""
This rule concatenate all chromosomes of a sample, from strand combined bed
"""

rule concat_all_chromosome_per_sample:
    input:
        uncomb_bed = expand(os.path.normpath(OUTPUT_DIR + "/tmp/bed_uncombined_strands/{{sample_name}}/{{meth_type}}/{{sample_name}}_chr{chr_number}_{{meth_type}}_uncomb.bed"), chr_number = CHR_NUMBER),
    output:
        bed_uncomb_concat_formated = os.path.normpath(OUTPUT_DIR + "/tmp/bed_uncombined_strands/{sample_name}/{meth_type}/{sample_name}_{meth_type}_uncomb.bed")
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt,
        time_min = lambda wildcards, attempt: attempt * 60
    shell:
        """
        sed "s/ /\t/g" {input.uncomb_bed} | awk  '{{print $1,$2,$3,$11,$12,$10-$12}}'> {output.bed_uncomb_concat_formated}
        """

"""
This rule launch methylation QC barplot of methylated GC
"""

def params_bed_uncomb_concat_formated_all_samples(wildcards):
        return  ",".join(list(dict.fromkeys(expand(os.path.normpath(OUTPUT_DIR + "/tmp/bed_uncombined_strands/{sample_name}/"+ wildcards.meth_type  + "/{sample_name}_" + wildcards.meth_type + "_uncomb.bed"), sample_name = SAMPLE_NAME))))

rule meth_QC_barplot_CG:
    input:
        bed_comb_concat_formated_all_samples = expand(os.path.normpath(OUTPUT_DIR + "/tmp/bed_uncombined_strands/{sample_name}/{{meth_type}}/{sample_name}_{{meth_type}}_uncomb.bed"), sample_name = SAMPLE_NAME),
        cg_motif = os.path.normpath(OUTPUT_DIR + "/tmp/resources/motif_cg.bed"),
    output:
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/methylation_QC/barplot_methylated_CG_{meth_type}_mqc.png")
    params:
        bed_uncomb_concat_formated_all_samples = params_bed_uncomb_concat_formated_all_samples,
        script = os.path.normpath(PIPELINE_DIR + "/script")
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: min(5120 + 5120 * (attempt - 1 ),20480),
        time_min = (lambda wildcards, attempt: attempt * 300)
    shell:
        """
        #module load r/3.6.1
        singularity exec --contain -B {params.script},{OUTPUT_DIR} {SING_ENV_GENEDMR} \
        Rscript {params.script}/QC_methylation_barplot_methylated_CG.R --list_input_bed_uncomb_concat {params.bed_uncomb_concat_formated_all_samples} --input_bed_motif_cg {input.cg_motif} --output_path {OUTPUT_DIR}/Quality_Control/methylation_QC/ 
        """





"""\ #je dois encore ecrire cette partie là
rule meth_QC_percentage_methylation_per_sample:
    input:
        bed=
        nanoplot=
    output:
        table
    threads: 2
    resources:
        mem_mb=1024,
        time_min = (lambda wildcards, attempt: attempt * 300)
    shell:
        
        module load r/3.6.1
        time Rscript /mnt/beegfs/userdata/y_mesloub/script/graph_QC_meth.R ${input.bed}


"""




"""
##########################################################################
These rules launch methylation analysis
##########################################################################

"""

def ignore_meth(wildcards):
    if wildcards.meth_type=="5mCG":
        return "m"
    elif wildcards.meth_type=="5hmCG":
        return "h"
    else:
        return None

"""
This # launches methylations separation#

"""
# modkit_separate_mod:#
#    input:
#       bam = os.path.normpath(OUTPUT_DIR + "/tmp/split_chr/{sample_name}/{sample_name}_chr_{chr_number}.bam")
#    output:
#        bam = temp(os.path.normpath(OUTPUT_DIR + "/tmp/separate_mod/{sample_name}/{meth_type}/{sample_name}_chr_{chr_number}_{meth_type}.bam"))
#    params:
#        meth_to_ignore = ignore_meth
#    threads: 12
#    resources:
#        mem_mb=lambda wildcards, attempt: min(5120 + 5120 * (attempt - 1 ), 20480),
#        time_min = (lambda wildcards, attempt: attempt * 300)
#    shell:
#        """
#        echo {wildcards.meth_type}
#        TMP_DIR=$(mktemp -d -t lr_pipeline-XXXXXXXXXX) && \
#        singularity exec --contain -B {OUTPUT_DIR} -B ${{TMP_DIR}}:/tmp \
#        {SING_ENV_MODKIT} modkit adjust-mods -t 12 --ignore {params.meth_to_ignore} {input.bam} {output.bam}
#        """

"""
This rule indexes the BAM file from modkit
"""
#rule bam_index:
#    input:
#        bam = os.path.normpath(OUTPUT_DIR + "/tmp/separate_mod/{sample_name}/{meth_type}/{sample_name}_chr_{chr_number}_{meth_type}.bam")
#    output:
#        index = temp(os.path.normpath(OUTPUT_DIR + "/tmp/separate_mod/{sample_name}/{meth_type}/{sample_name}_chr_{chr_number}_{meth_type}.bam.bai"))
#    threads: 1
#    resources:
#        mem_mb=lambda wildcards, attempt: min(1024 + 5120 * (attempt - 1 ), 20480),
#        time_min = (lambda wildcards, attempt: attempt * 300)
#    conda:
#        CONDA_ENV_SAMTOOLS
#    shell:
#        """
#       	samtools index {input.bam}
#       	"""

"""
These 2 rules launch BED files generation (strand uncombined and combined)
"""

rule modkit_pileup_uncomb:
    input:
        bam = os.path.normpath(OUTPUT_DIR + "/tmp/split_chr/{sample_name}/{sample_name}_chr_{chr_number}.bam"),
        index = os.path.normpath(OUTPUT_DIR + "/tmp/split_chr/{sample_name}/{sample_name}_chr_{chr_number}.bam.bai")
    output:
        uncomb_bed = temp(os.path.normpath(OUTPUT_DIR + "/tmp/bed_uncombined_strands/{sample_name}/{sample_name}_chr{chr_number}_uncomb.bed"))
    params:
        reference = config["references"]["genome"],
        ref_path = os.path.normpath(config["references"]["genome"])
    log:
        os.path.normpath(OUTPUT_DIR + "/tmp/logs/modkit_pileup/{sample_name}_chr{chr_number}_uncombined.log")
    threads: 12
    resources:
        mem_mb=lambda wildcards, attempt: min(5120 + 10240 * (attempt - 1 ), 102400),
        time_min = (lambda wildcards, attempt: attempt * 300)
    shell:
        """
        TMP_DIR=$(mktemp -d -t lr_pipeline-XXXXXXXXXX) && \
        singularity exec --contain -B {OUTPUT_DIR},{params.ref_path} -B ${{TMP_DIR}}:/tmp \
        {SING_ENV_MODKIT} modkit pileup --ref {params.reference} --cpg --threads {threads} --log-filepath {log} {input.bam} {output.uncomb_bed}
        """


rule modkit_pileup_comb:
    input:
        bam = os.path.normpath(OUTPUT_DIR + "/tmp/split_chr/{sample_name}/{sample_name}_chr_{ nchr_number}.bam"),
        index = os.path.normpath(OUTPUT_DIR + "/tmp/split_chr/{sample_name}/{sample_name}_chr_{chr_number}.bam.bai")
    output:
        comb_bed = temp(os.path.normpath(OUTPUT_DIR + "/bed_combined_strands/{sample_name}/{sample_name}_chr{chr_number}_comb.bed"))
    params:
        reference = config["references"]["genome"],
        ref_path = os.path.normpath(config["references"]["genome"])
    log:
        os.path.normpath(OUTPUT_DIR + "/tmp/logs/modkit_pileup/{sample_name}_chr_{chr_number}_combined.log")
    threads: 12
    resources:
        mem_mb=lambda wildcards, attempt: min(5120 + 10240 * (attempt - 1 ), 102400),
        time_min = (lambda wildcards, attempt: attempt * 300)
    shell:
        """
        TMP_DIR=$(mktemp -d -t lr_pipeline-XXXXXXXXXX) && \
        singularity exec --contain -B {OUTPUT_DIR},{params.ref_path} -B ${{TMP_DIR}}:/tmp \
        {SING_ENV_MODKIT} modkit pileup --ref {params.reference} --cpg --threads {threads}  --combine-strands --log-filepath {log} {input.bam} {output.comb_bed}
        
        """


rule split_by_meth_uncomb:
    input:
        uncomb_bed = os.path.normpath(OUTPUT_DIR + "/tmp/bed_uncombined_strands/{sample_name}/{sample_name}_chr{chr_number}_uncomb.bed") 
    output:
        uncomb_meth_bed = temp(os.path.normpath(OUTPUT_DIR + "/tmp/bed_uncombined_strands/{sample_name}/{meth_type}/{sample_name}_chr{chr_number}_{meth_type}_uncomb.bed"))
    params:
        meth_to_keep = ignore_meth
    threads:1
    resources:
        mem_mb=2048,
        time_min = (lambda wildcards, attempt: attempt * 300)
    shell:
        """
        cat {input.uncomb_bed}|awk -v meth={params.meth_to_keep} '$4==meth' >{output.uncomb_meth_bed}
        
        """







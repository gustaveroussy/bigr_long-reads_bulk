"""
##########################################################################
These rules launches the splitting, filtering and concatenation of bam files
##########################################################################

"""


"""
This rule makes the symbolic links of bam files with the good sample name.
"""

def symlink_rename_input_bam(wildcards):
    index = SAMPLE_NAME.index(wildcards.sample_name)
    return ORIG_FILE[index]

rule symlink_rename_bam:
    input:
        bam = symlink_rename_input_bam
    output:
        bam_link = temp(os.path.normpath(OUTPUT_DIR + "/concat_sort/{sample_name}/{sample_name}_sorted.bam")),
        index_link = temp(os.path.normpath(OUTPUT_DIR + "/concat_sort/{sample_name}/{sample_name}_sorted.bam.bai"))
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 2048, 20480)),
        time_min = (lambda wildcards, attempt: attempt * 300)
    run:
        #sys.stderr.write("\t Create symbolic link: \n")
        #sys.stderr.write("\t From :" + "\t" + str(input.bam) + "\n")
        #sys.stderr.write("\t To :" + "\t" + str(output.bam_link) + "\n")
        #os.symlink(str(input.bam), str(output.bam_link))
        #os.symlink(str(input.index_link), str(output.index_link))
        import shutil
        sys.stderr.write("\t Copy file: \n")
        sys.stderr.write("\t From :" + "\t" + str(input.bam) + "\n")
        sys.stderr.write("\t To :" + "\t" + str(output.bam_link) + "\n")
        shutil.copy(str(input.bam), str(output.bam_link))
        shutil.copy(str(input.bam) + ".bai", str(output.index_link ))

"""
This rule splits sorted BAM by chromosome
"""

rule split_sorted_bam:
    input:
        bam = os.path.normpath(OUTPUT_DIR + "/concat_sort/{sample_name}/{sample_name}_sorted.bam"),
        index = os.path.normpath(OUTPUT_DIR + "/concat_sort/{sample_name}/{sample_name}_sorted.bam.bai")
    output:
        bam = temp(os.path.normpath(OUTPUT_DIR + "/split_chr/{sample_name}/{sample_name}_chr_{chr_number}.bam")),
        index = temp(os.path.normpath(OUTPUT_DIR + "/split_chr/{sample_name}/{sample_name}_chr_{chr_number}.bam.bai"))
    threads: 12
    resources:
        mem_mb = lambda wildcards, attempt: min(10240 + 5120 * (attempt - 1 ), 102400),
        time_min = (lambda wildcards, attempt: attempt * 300)
    conda:
        CONDA_ENV_SAMTOOLS
    shell:
        """
       	samtools view -b -@ 12 -o {output.bam} {input.bam} {wildcards.chr_number}
       	samtools index {output.bam}
        """

"""
This rule re-concatenates all split BAM files per sample
"""

rule reconcat_split_bam:
    input:
        bams = expand(os.path.normpath(OUTPUT_DIR + "/split_chr/{{sample_name}}/{{sample_name}}_chr_{chr_number}.bam"),chr_number=CHR_NUMBER)
    output:
        bam = temp(os.path.normpath(OUTPUT_DIR + "/reconcat/{sample_name}/{sample_name}_reconcat.bam")),
        bam_sorted = temp(os.path.normpath(OUTPUT_DIR + "/reconcat/{sample_name}/{sample_name}_sorted.bam")),
        index = temp(os.path.normpath(OUTPUT_DIR + "/reconcat/{sample_name}/{sample_name}_sorted.bam.bai"))
    params:
        tmp = os.path.normpath(OUTPUT_DIR + "/tmp/reconcat/")
    threads: 12
    resources:
        mem_mb = lambda wildcards, attempt: min(20480 + 51200 * (attempt - 1 ), 204800),
        time_min = (lambda wildcards, attempt: attempt * 300)
    conda:
        CONDA_ENV_SAMTOOLS
    shell:
        """
       	samtools cat {input.bams} > {output.bam}
       	mkdir -m 771 -p {params.tmp}
       	samtools sort -@ 12 -T {params.tmp} -o {output.bam_sorted} {output.bam}
       	samtools index {output.bam_sorted}
       	"""
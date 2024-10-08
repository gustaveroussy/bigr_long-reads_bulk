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
        bam_link = temp(os.path.normpath(OUTPUT_DIR + "/Bam/{sample_name}/{sample_name}_sorted.bam")),
        index_link = temp(os.path.normpath(OUTPUT_DIR + "/Bam/{sample_name}/{sample_name}_sorted.bam.bai"))
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 2048, 20480)),
        time_min = (lambda wildcards, attempt: attempt * 300)
    run:
        import shutil
        sys.stderr.write("\t Copy file: \n")
        sys.stderr.write("\t From :" + "\t" + str(input.bam) + "\n")
        sys.stderr.write("\t To :" + "\t" + str(output.bam_link) + "\n")
        shutil.copy(str(input.bam), str(output.bam_link))
        shutil.copy(str(input.bam) + ".bai", str(output.index_link ))

"""
If basecalling mode is set to "methylation":
This rules checks if the bam file contains methylation information before proceeding to the following steps
"""
if config["basecalling_mode"] == "methylation" and config["input_format"] == "bam":
    rule check_bam_methylation:
        input: 
            bam = os.path.normpath(OUTPUT_DIR + "/Bam/{sample_name}/{sample_name}_sorted.bam")
        output: 
            flag = temp(os.path.normpath(OUTPUT_DIR + "/tmp/bam/{sample_name}/{sample_name}_check_bam_methylation_OK.txt"))
        params:
            script = os.path.normpath(PIPELINE_DIR + "/script")
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(1024 + 5120 * (attempt - 1 ),2048),
            time_min = (lambda wildcards, attempt: attempt * 60)
        conda:
            CONDA_ENV_PYTHON
        shell:
            """
            flag = $( python3 {params.script}/check_methylation_in_bam.py {input.bam} )

            if [ "$flag" = "True" ]; then
                echo "Bam file seems to contain methylation information." > {output.flag}
            elif [ "$flag" = "False" ]; then
                echo "Bam file does not seem to contain methylation information as MM read tag was not found for the first read from the input BAM file."
            else
                echo "It looks like the script encountered an issue while running."
            fi
            """

"""
This rule splits sorted BAM by chromosome
"""
def split_sorted_bam_input(wildcards):
    input = []
    if config["basecalling_mode"] == "methylation" and config["input_format"] == "bam":
        input.append(os.path.normpath(OUTPUT_DIR + "/tmp/bam/{sample_name}/{sample_name}_check_bam_methylation_OK.txt"))
    return input
    
rule split_sorted_bam:
    input:
        bam = os.path.normpath(OUTPUT_DIR + "/Bam/{sample_name}/{sample_name}_sorted.bam"),
        index = os.path.normpath(OUTPUT_DIR + "/Bam/{sample_name}/{sample_name}_sorted.bam.bai"),
        flag = split_sorted_bam_input
    output:
        bam = temp(os.path.normpath(OUTPUT_DIR + "/tmp/split_chr/{sample_name}/{sample_name}_chr_{chr_number}.bam")),
        index = temp(os.path.normpath(OUTPUT_DIR + "/tmp/split_chr/{sample_name}/{sample_name}_chr_{chr_number}.bam.bai"))
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
        bams = expand(os.path.normpath(OUTPUT_DIR + "/tmp/split_chr/{{sample_name}}/{{sample_name}}_chr_{chr_number}.bam"),chr_number=CHR_NUMBER)
    output:
        bam = temp(os.path.normpath(OUTPUT_DIR + "/tmp/reconcat/{sample_name}/{sample_name}_reconcat.bam")),
        bam_sorted = temp(os.path.normpath(OUTPUT_DIR + "/tmp/reconcat/{sample_name}/{sample_name}_sorted.bam")),
        index = temp(os.path.normpath(OUTPUT_DIR + "/tmp/reconcat/{sample_name}/{sample_name}_sorted.bam.bai"))
    params:
        tmp = os.path.normpath(OUTPUT_DIR + "/tmp/tmp/reconcat/")
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

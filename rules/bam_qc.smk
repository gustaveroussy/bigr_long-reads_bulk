"""
##########################################################################
These rules make the control-quality of the alignment
##########################################################################
"""
wildcard_constraints:
    bam_name = '|'.join([x for x in BAM_NAME])

"""
This rule makes the symbolic links of fastq files with the good sample name.
"""

def symlink_rename_input_bam(wildcards):
    index = BAM_NAME.index(wildcards.bam_name)
    return ORIG_FILE[index]

rule symlink_rename_bam:
    input:
        bam = symlink_rename_input_bam
    output:
        bam_link = temp(os.path.normpath(OUTPUT_DIR + "/symlink_input/" + "/{bam_name}.bam"))
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
        #os.symlink(str(input.bam + ".bai"), str(output.bam_link + ".bai"))
        import shutil
        sys.stderr.write("\t Copy file: \n")
        sys.stderr.write("\t From :" + "\t" + str(input.bam) + "\n")
        sys.stderr.write("\t To :" + "\t" + str(output.bam_link) + "\n")
        shutil.copy(str(input.bam), str(output.bam_link))
        shutil.copy(str(input.bam + ".bai"), str(output.bam_link + ".bai"))


"""
This rule makes the qualimap
"""

def qualimap_input_bam(wildcards):
    index = BAM_NAME.index(wildcards.bam_name)
    return SYMLINK_FILES[index]

rule qualimap:
    input:
        bam_file = qualimap_input_bam
    output:
        os.path.normpath(OUTPUT_DIR + "/bam_QC/qualimap/{bam_name}/qualimapReport.html"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/qualimap/{bam_name}/raw_data_qualimapReport/genome_fraction_coverage.txt"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/qualimap/{bam_name}/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/qualimap/{bam_name}/genome_results.txt"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/qualimap/{bam_name}/raw_data_qualimapReport/coverage_histogram.txt")
        #add css & images_qualimapReport ??
    threads:
        10
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 720)
    conda:
        CONDA_ENV_QUALITYMAP
    shell:
        """
        res=$(({resources.mem_mb}/1000)) && \
        qualimap bamqc -bam {input.bam_file} -outdir {OUTPUT_DIR}/bam_QC/qualimap/{wildcards.bam_name}/ --paint-chromosome-limits -nt {threads} --java-mem-size=${{res}}G
        
        """


"""
This rule makes the mosdepth
"""

def mosdepth_input_bam(wildcards):
    index = BAM_NAME.index(wildcards.bam_name)
    return SYMLINK_FILES[index]

rule mosdepth:
    input:
        bam_file = mosdepth_input_bam
    output:
        os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}_wgs_mode.mosdepth.global.dist.txt"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}_wgs_mode.mosdepth.region.dist.txt"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}_wgs_mode.mosdepth.summary.txt"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}_wgs_mode.regions.bed.gz"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}_wgs_mode.regions.bed.gz.csi")
    threads:
        3
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 720)
    conda:
        CONDA_ENV_MOSDEPTH
    shell:
        """
        cd {OUTPUT_DIR}/bam_QC/mosdepth/{wildcards.bam_name}/
        mosdepth -n --fast-mode --by 500 -t {threads} {wildcards.bam_name}_wgs_mode {input.bam_file} && \
        mosdepth -t {threads} {wildcards.bam_name} {input.bam_file}
        
        """


"""
This rule makes the mosdepth
"""

def nanoplot_input_bam(wildcards):
    index = BAM_NAME.index(wildcards.bam_name)
    return SYMLINK_FILES[index]

rule nanoplot_bam:
    input:
        bam_file = nanoplot_input_bam
    output:
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_AlignedReadlengthvsSequencedReadLength_dot.html"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_AlignedReadlengthvsSequencedReadLength_dot.png"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_AlignedReadlengthvsSequencedReadLength_kde.html"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_AlignedReadlengthvsSequencedReadLength_kde.png"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_LengthvsQualityScatterPlot_dot.html"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_LengthvsQualityScatterPlot_dot.png"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_LengthvsQualityScatterPlot_kde.html"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_LengthvsQualityScatterPlot_kde.png"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_MappingQualityvsAverageBaseQuality_dot.html"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_MappingQualityvsAverageBaseQuality_dot.png"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_MappingQualityvsAverageBaseQuality_kde.html"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_MappingQualityvsAverageBaseQuality_kde.png"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_MappingQualityvsReadLength_dot.html"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_MappingQualityvsReadLength_dot.png"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_MappingQualityvsReadLength_kde.html"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_MappingQualityvsReadLength_kde.png"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_NanoPlot-report.html"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_NanoStats.txt"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_Non_weightedHistogramReadlength.html"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_Non_weightedHistogramReadlength.png"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_Non_weightedLogTransformed_HistogramReadlength.html"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_Non_weightedLogTransformed_HistogramReadlength.png"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityHistogramDynamic_Histogram_percent_identity.html"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityHistogramDynamic_Histogram_percent_identity.png"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityvsAlignedReadLength_dot.html"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityvsAlignedReadLength_dot.png"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityvsAlignedReadLength_kde.html"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityvsAlignedReadLength_kde.png"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityvsAverageBaseQuality_dot.html"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityvsAverageBaseQuality_dot.png"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityvsAverageBaseQuality_kde.html"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityvsAverageBaseQuality_kde.png"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_WeightedHistogramReadlength.html"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_WeightedHistogramReadlength.png"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_WeightedLogTransformed_HistogramReadlength.html"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_WeightedLogTransformed_HistogramReadlength.png"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_Yield_By_Length.html"),
        os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_Yield_By_Length.png")
    threads:
        10
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 720)
    shell:
        """
        TMP_DIR=$(mktemp -d -t lr_pipeline-XXXXXXXXXX) && \
        singularity exec --contain -B {OUTPUT_DIR} -B ${{TMP_DIR}}:/tmp -B ${{TMP_DIR}}:${{TMPDIR}} {SING_ENV_NANOPLOT} \
        NanoPlot --threads {threads} --outdir {OUTPUT_DIR}/bam_QC/nanoplot/{wildcards.bam_name}/ --prefix {wildcards.bam_name}_ --N50 --tsv_stats --info_in_report --bam {input.bam_file}

        """



"""
##########################################################################
This ruls summarise the control-quality of the alignment
##########################################################################
"""
#wildcard_constraints:
#    sample_name = '|'.join([x for x in SAMPLE_NAME])

"""
This rule agglomerates qc into one html file thanks to multiqc
"""

rule multiqc:
    input:
        #qualimap
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/qualimap/{sample_name}/qualimapReport.html"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/qualimap/{sample_name}/raw_data_qualimapReport/genome_fraction_coverage.txt"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/qualimap/{sample_name}/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/qualimap/{sample_name}/genome_results.txt"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/qualimap/{sample_name}/raw_data_qualimapReport/coverage_histogram.txt"), sample_name=SAMPLE_NAME),
        #mosdepth
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/mosdepth/{sample_name}/{sample_name}_wgs_mode.mosdepth.global.dist.txt"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/mosdepth/{sample_name}/{sample_name}_wgs_mode.mosdepth.region.dist.txt"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/mosdepth/{sample_name}/{sample_name}_wgs_mode.mosdepth.summary.txt"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/mosdepth/{sample_name}/{sample_name}_wgs_mode.regions.bed.gz"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/mosdepth/{sample_name}/{sample_name}_wgs_mode.regions.bed.gz.csi"), sample_name=SAMPLE_NAME),
        #nanoplot
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_AlignedReadlengthvsSequencedReadLength_dot.html"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_AlignedReadlengthvsSequencedReadLength_dot.png"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_AlignedReadlengthvsSequencedReadLength_kde.html"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_AlignedReadlengthvsSequencedReadLength_kde.png"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_LengthvsQualityScatterPlot_dot.html"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_LengthvsQualityScatterPlot_dot.png"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_LengthvsQualityScatterPlot_kde.html"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_LengthvsQualityScatterPlot_kde.png"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_MappingQualityvsAverageBaseQuality_dot.html"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_MappingQualityvsAverageBaseQuality_dot.png"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_MappingQualityvsAverageBaseQuality_kde.html"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_MappingQualityvsAverageBaseQuality_kde.png"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_MappingQualityvsReadLength_dot.html"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_MappingQualityvsReadLength_dot.png"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_MappingQualityvsReadLength_kde.html"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_MappingQualityvsReadLength_kde.png"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_NanoPlot-report.html"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_NanoStats.txt"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_Non_weightedHistogramReadlength.html"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_Non_weightedHistogramReadlength.png"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_Non_weightedLogTransformed_HistogramReadlength.html"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_Non_weightedLogTransformed_HistogramReadlength.png"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_PercentIdentityHistogramDynamic_Histogram_percent_identity.html"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_PercentIdentityHistogramDynamic_Histogram_percent_identity.png"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_PercentIdentityvsAlignedReadLength_dot.html"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_PercentIdentityvsAlignedReadLength_dot.png"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_PercentIdentityvsAlignedReadLength_kde.html"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_PercentIdentityvsAlignedReadLength_kde.png"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_PercentIdentityvsAverageBaseQuality_dot.html"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_PercentIdentityvsAverageBaseQuality_dot.png"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_PercentIdentityvsAverageBaseQuality_kde.html"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_PercentIdentityvsAverageBaseQuality_kde.png"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_WeightedHistogramReadlength.html"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_WeightedHistogramReadlength.png"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_WeightedLogTransformed_HistogramReadlength.html"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_WeightedLogTransformed_HistogramReadlength.png"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_Yield_By_Length.html"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_Yield_By_Length.png"), sample_name=SAMPLE_NAME),
        #fastqc
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/fastq_QC/fastqc/{sample_name}/{sample_name}_fastqc.zip"), sample_name=SAMPLE_NAME),
        #fastq_screen
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/fastq_QC/fastq_screen/{sample_name}/{sample_name}_screen.txt"), sample_name=SAMPLE_NAME),
        #methylation qc
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/methylation_QC/per_pos_per_strand_gam_ratio_{sample_name}_chr{chr_number}_{meth_type}_mqc.png"), chr_number = CHR_NUMBER, sample_name = SAMPLE_NAME, meth_type = METH_TYPE),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/methylation_QC/barplot_methylated_CG_{meth_type}_mqc.png"), meth_type = METH_TYPE)
    output:
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/multiqc_report.html"),
        temp(directory(os.path.normpath(OUTPUT_DIR + "/Quality_Control/multiqc_data/")))
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 720)
    conda:
        CONDA_ENV_MULTIQC
    shell:
        """
        cd {OUTPUT_DIR}/Quality_Control/
        multiqc {OUTPUT_DIR}/Quality_Control/bam_QC/qualimap/ {OUTPUT_DIR}/Quality_Control/bam_QC/mosdepth/ {OUTPUT_DIR}/Quality_Control/bam_QC/nanoplot/ {OUTPUT_DIR}/Quality_Control/fastq_QC/fastqc/ {OUTPUT_DIR}/Quality_Control/fastq_QC/fastq_screen/ {OUTPUT_DIR}/Quality_Control/methylation_QC/ --config {PIPELINE_DIR}/config/multiqc_config.yaml
        
        """

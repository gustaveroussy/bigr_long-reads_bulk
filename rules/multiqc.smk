"""
##########################################################################
This ruls summarise the control-quality of the alignment
##########################################################################
"""
wildcard_constraints:
    bam_name = '|'.join([x for x in BAM_NAME])

"""
This rule agglomerates bam qc into one html file thanks to multiqc
"""

rule multiqc:
    input:
        #qualimap
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/qualimap/{bam_name}/qualimapReport.html"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/qualimap/{bam_name}/raw_data_qualimapReport/genome_fraction_coverage.txt"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/qualimap/{bam_name}/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/qualimap/{bam_name}/genome_results.txt"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/qualimap/{bam_name}/raw_data_qualimapReport/coverage_histogram.txt"), bam_name=BAM_NAME),
        #mosdepth
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}_wgs_mode.mosdepth.global.dist.txt"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}_wgs_mode.mosdepth.region.dist.txt"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}_wgs_mode.mosdepth.summary.txt"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}_wgs_mode.regions.bed.gz"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/mosdepth/{bam_name}/{bam_name}_wgs_mode.regions.bed.gz.csi"), bam_name=BAM_NAME),
        #nanoplot
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_AlignedReadlengthvsSequencedReadLength_dot.html"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_AlignedReadlengthvsSequencedReadLength_dot.png"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_AlignedReadlengthvsSequencedReadLength_kde.html"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_AlignedReadlengthvsSequencedReadLength_kde.png"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_LengthvsQualityScatterPlot_dot.html"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_LengthvsQualityScatterPlot_dot.png"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_LengthvsQualityScatterPlot_kde.html"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_LengthvsQualityScatterPlot_kde.png"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_MappingQualityvsAverageBaseQuality_dot.html"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_MappingQualityvsAverageBaseQuality_dot.png"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_MappingQualityvsAverageBaseQuality_kde.html"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_MappingQualityvsAverageBaseQuality_kde.png"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_MappingQualityvsReadLength_dot.html"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_MappingQualityvsReadLength_dot.png"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_MappingQualityvsReadLength_kde.html"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_MappingQualityvsReadLength_kde.png"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_NanoPlot-report.html"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_NanoStats.txt"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_Non_weightedHistogramReadlength.html"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_Non_weightedHistogramReadlength.png"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_Non_weightedLogTransformed_HistogramReadlength.html"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_Non_weightedLogTransformed_HistogramReadlength.png"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityHistogramDynamic_Histogram_percent_identity.html"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityHistogramDynamic_Histogram_percent_identity.png"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityvsAlignedReadLength_dot.html"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityvsAlignedReadLength_dot.png"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityvsAlignedReadLength_kde.html"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityvsAlignedReadLength_kde.png"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityvsAverageBaseQuality_dot.html"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityvsAverageBaseQuality_dot.png"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityvsAverageBaseQuality_kde.html"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_PercentIdentityvsAverageBaseQuality_kde.png"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_WeightedHistogramReadlength.html"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_WeightedHistogramReadlength.png"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_WeightedLogTransformed_HistogramReadlength.html"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_WeightedLogTransformed_HistogramReadlength.png"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_Yield_By_Length.html"), bam_name=BAM_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_Yield_By_Length.png"), bam_name=BAM_NAME)
    output:
        os.path.normpath(OUTPUT_DIR + "/bam_QC/multiqc_report.html"),
        temp(directory(os.path.normpath(OUTPUT_DIR + "/bam_QC/multiqc_data/")))
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 720)
    conda:
        CONDA_ENV_MULTIQC
    shell:
        """
        cd {OUTPUT_DIR}/bam_QC/
        multiqc {OUTPUT_DIR}/bam_QC/qualimap/ {OUTPUT_DIR}/bam_QC/mosdepth/ {OUTPUT_DIR}/bam_QC/nanoplot/ --config {PIPELINE_DIR}/config/multiqc_config.yaml
        
        """

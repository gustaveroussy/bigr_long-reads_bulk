"""
##########################################################################
This function make the input of rule all
##########################################################################
"""

def get_targets():
  targets = {}
  if config["steps"]["bam_qc"]:
      targets["bam_qc"]=[
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
        expand(os.path.normpath(OUTPUT_DIR + "/bam_QC/nanoplot/{bam_name}/{bam_name}_Yield_By_Length.png"), bam_name=BAM_NAME),
        #multiqc
        os.path.normpath(OUTPUT_DIR + "/bam_QC/multiqc_report.html")
        ]
  if config["steps"]["snv_calling"]:
      if config["variant_calling_mode"] == "germline":
        targets["snv_calling"]=[
            #clair3
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{sample_name}/{sample_name}_merge_output.vcf.gz"), sample_name=SAMPLE_NAME, claire3_model=NAME_CLAIR3_MODEL),
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{sample_name}/{sample_name}_merge_output.vcf.gz.tbi"), sample_name=SAMPLE_NAME, claire3_model=NAME_CLAIR3_MODEL)
            #pepper_margin_deepvariant
        ]
        targets["snv_annotation"]=[
            #clair3 & snpEff & snpSift
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{sample_name}/snpEff/{sample_name}{compl}_annotated.vcf"), sample_name=SAMPLE_NAME, claire3_model=NAME_CLAIR3_MODEL, compl="_merge_output"),
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{sample_name}/snpEff/{sample_name}{compl}_annotated_stats.csv"), sample_name=SAMPLE_NAME, claire3_model=NAME_CLAIR3_MODEL, compl="_merge_output"),
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{sample_name}/snpEff/{sample_name}{compl}_annotated_summary.html"), sample_name=SAMPLE_NAME, claire3_model=NAME_CLAIR3_MODEL, compl="_merge_output"),
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{sample_name}/snpEff/{sample_name}{compl}_annotated_dbnsfp.vcf"), sample_name=SAMPLE_NAME, claire3_model=NAME_CLAIR3_MODEL, compl="_merge_output"),
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{sample_name}/snpEff/{sample_name}{compl}_annotated_dbnsfp_clinvar.vcf"), sample_name=SAMPLE_NAME, claire3_model=NAME_CLAIR3_MODEL, compl="_merge_output"),
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{sample_name}/snpEff/{sample_name}{compl}_annotated_{filter}.vcf"), sample_name=SAMPLE_NAME, claire3_model=NAME_CLAIR3_MODEL, filter=SNPSIFT_FILTERS_NAMES, compl="_merge_output")
            #pepper_margin_deepvariant & snpEff & snpSift
        ]
      if config["variant_calling_mode"] == "somatic":
        targets["snv_calling"]=[
            #clairs
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clairs/{sample_name}/{sample_name}{compl}.vcf.gz"), sample_name=SAMPLE_NAME, compl=["_snv","_indel"])
        ]
        targets["snv_annotation"]=[
            #clairs & snpEff & snpSift
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clairs/{sample_name}/snpEff/{sample_name}{compl}_annotated.vcf"), sample_name=SAMPLE_NAME, compl=["_snv","_indel"]),
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clairs/{sample_name}/snpEff/{sample_name}{compl}_annotated_stats.csv"), sample_name=SAMPLE_NAME, compl=["_snv","_indel"]),
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clairs/{sample_name}/snpEff/{sample_name}{compl}_annotated_summary.html"), sample_name=SAMPLE_NAME, compl=["_snv","_indel"]),
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clairs/{sample_name}/snpEff/{sample_name}{compl}_annotated_dbnsfp.vcf"), sample_name=SAMPLE_NAME, compl=["_snv","_indel"]),
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clairs/{sample_name}/snpEff/{sample_name}{compl}_annotated_dbnsfp_clinvar.vcf"), sample_name=SAMPLE_NAME, compl=["_snv","_indel"]),
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clairs/{sample_name}/snpEff/{sample_name}{compl}_annotated_{filter}.vcf"), sample_name=SAMPLE_NAME, filter=SNPSIFT_FILTERS_NAMES, compl=["_snv","_indel"])
        ]
  if config["steps"]["phasing"]:
      if config["variant_calling_mode"] == "germline":
        targets["phasing"]=[
          #clair3
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{sample_name}/whatshap/{sample_name}{compl}_phased.vcf.gz"), sample_name=SAMPLE_NAME, claire3_model=NAME_CLAIR3_MODEL, compl="_merge_output"),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{sample_name}/whatshap/{sample_name}{compl}_phasing_stats.txt"), sample_name=SAMPLE_NAME, claire3_model=NAME_CLAIR3_MODEL, compl="_merge_output"),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{sample_name}/whatshap/{sample_name}{compl}_phasing_stats.tsv"), sample_name=SAMPLE_NAME, claire3_model=NAME_CLAIR3_MODEL, compl="_merge_output"),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{sample_name}/whatshap/{sample_name}{compl}_phasing_haplotype_blocks.tsv"), sample_name=SAMPLE_NAME, claire3_model=NAME_CLAIR3_MODEL, compl="_merge_output"),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{sample_name}/whatshap/{sample_name}{compl}_phasing_haplotype_blocks.gtf"), sample_name=SAMPLE_NAME, claire3_model=NAME_CLAIR3_MODEL, compl="_merge_output"),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{sample_name}/whatshap/{sample_name}{compl}_phased.vcf.gz.tbi"), sample_name=SAMPLE_NAME, claire3_model=NAME_CLAIR3_MODEL, compl="_merge_output"),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{sample_name}/whatshap/{sample_name}{compl}_phasing_haplotag_list.tsv"), sample_name=SAMPLE_NAME, claire3_model=NAME_CLAIR3_MODEL, compl="_merge_output"),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clair3/{claire3_model}/{sample_name}/whatshap/{sample_name}{compl}_haplotagged.bam"), sample_name=SAMPLE_NAME, claire3_model=NAME_CLAIR3_MODEL, compl="_merge_output")
          #pepper_margin_deepvariant
        ]
      if config["variant_calling_mode"] == "somatic":
        targets["phasing"]=[
          #clairs
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clairs/{sample_name}/whatshap/{sample_name}{compl}_phased.vcf.gz"), sample_name=SAMPLE_NAME, compl=["_snv","_indel"]),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clairs/{sample_name}/whatshap/{sample_name}{compl}_phasing_stats.txt"), sample_name=SAMPLE_NAME, compl=["_snv","_indel"]),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clairs/{sample_name}/whatshap/{sample_name}{compl}_phasing_stats.tsv"), sample_name=SAMPLE_NAME, compl=["_snv","_indel"]),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clairs/{sample_name}/whatshap/{sample_name}{compl}_phasing_haplotype_blocks.tsv"), sample_name=SAMPLE_NAME, compl=["_snv","_indel"]),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clairs/{sample_name}/whatshap/{sample_name}{compl}_phasing_haplotype_blocks.gtf"), sample_name=SAMPLE_NAME, compl=["_snv","_indel"]),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clairs/{sample_name}/whatshap/{sample_name}{compl}_phased.vcf.gz.tbi"), sample_name=SAMPLE_NAME, compl=["_snv","_indel"]),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clairs/{sample_name}/whatshap/{sample_name}{compl}_phasing_haplotag_list.tsv"), sample_name=SAMPLE_NAME, compl=["_snv","_indel"]),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/clairs/{sample_name}/whatshap/{sample_name}{compl}_haplotagged.bam"), sample_name=SAMPLE_NAME, compl=["_snv","_indel"])
        ]
  if config["steps"]["sv_calling"]:
      if config["variant_calling_mode"] == "germline":
          targets["germline_sv_calling"]=[
            #sniffles
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/sniffles/{bam_name}/{bam_name}_SV.vcf"), bam_name=BAM_NAME),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/sniffles/{bam_name}/{bam_name}_SV.snf"), bam_name=BAM_NAME)
            #sniffles2-plot
            ]
          if len(BAM_NAME) > 1 :
            targets["germline_sv_calling"].append(os.path.normpath(OUTPUT_DIR + "/SV_Calling/sniffles/all_samples_SV.vcf"))
      if config["variant_calling_mode"] == "somatic":
          targets["germline_sv_calling"]=[
            #nanomonsv parsing
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{sample_name}{sample_type}.deletion.sorted.bed.gz"), sample_name=SAMPLE_NAME, sample_type=["_tumor", "_normal"]),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{sample_name}{sample_type}.deletion.sorted.bed.gz.tbi"), sample_name=SAMPLE_NAME, sample_type=["_tumor", "_normal"]),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{sample_name}{sample_type}.insertion.sorted.bed.gz"), sample_name=SAMPLE_NAME, sample_type=["_tumor", "_normal"]),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{sample_name}{sample_type}.insertion.sorted.bed.gz.tbi"), sample_name=SAMPLE_NAME, sample_type=["_tumor", "_normal"]),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{sample_name}{sample_type}.rearrangement.sorted.bedpe.gz"), sample_name=SAMPLE_NAME, sample_type=["_tumor", "_normal"]),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{sample_name}{sample_type}.rearrangement.sorted.bedpe.gz.tbi"), sample_name=SAMPLE_NAME, sample_type=["_tumor", "_normal"]),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{sample_name}{sample_type}.bp_info.sorted.bed.gz"), sample_name=SAMPLE_NAME, sample_type=["_tumor", "_normal"]),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{sample_name}{sample_type}.bp_info.sorted.bed.gz.tbi"), sample_name=SAMPLE_NAME, sample_type=["_tumor", "_normal"]),
            #nanomonsv SV
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{sample_name}.nanomonsv.result.txt"), sample_name=SAMPLE_NAME),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{sample_name}.nanomonsv.result.vcf"), sample_name=SAMPLE_NAME),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{sample_name}.nanomonsv.sbnd.result.txt"), sample_name=SAMPLE_NAME),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{sample_name}.nanomonsv.supporting_read.txt"), sample_name=SAMPLE_NAME)
          ]
          #if config["references"]["species"] == "human":
            #nanomonsv classifier
            #targets["germline_sv_calling"].append(expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/nanomonsv/{sample_name}/{sample_name}.nanomonsv.insert_classify.txt"), sample_name=SAMPLE_NAME))
  return targets


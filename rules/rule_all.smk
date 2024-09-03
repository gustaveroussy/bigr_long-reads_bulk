"""
##########################################################################
This function make the input of rule all
##########################################################################
"""

def get_targets():
  targets = {}
  
  #BASECALLING
  if config["steps"]["basecalling"]:
      targets["basecalling"]=[
        #basecalled ubam
        expand(os.path.normpath(OUTPUT_DIR + "/calling/{sample_name}/{batch_name}.bam"),zip,batch_name=BATCH_NAME,sample_name=SAMPLE_NAME),
      ]
  #ALIGNMENT
  if config["steps"]["alignment"]:
      targets["filter_align"]=[
        #filter
#        expand(os.path.normpath(OUTPUT_DIR + "/filtered/{sample_name}/{batches}_filtered.bam"),zip,batches=BATCH_NAME,sample_name=SAMPLE_NAME),
        #alignment
#        expand(os.path.normpath(OUTPUT_DIR + "/alignment/{sample_name}/{batches}_aligned.bam"),zip,batches=BATCH_NAME,sample_name=SAMPLE_NAME),
        #concat
#        expand(os.path.normpath(OUTPUT_DIR + "/concat_sort/{sample_name}/{sample_name}_concat.bam"),sample_name=SAMPLE_NAME),
        #sort
        expand(os.path.normpath(OUTPUT_DIR + "/concat_sort/{sample_name}/{sample_name}_sorted.bam"),sample_name=SAMPLE_NAME),
        #index
        expand(os.path.normpath(OUTPUT_DIR + "/concat_sort/{sample_name}/{sample_name}_sorted.bam.bai"),sample_name=SAMPLE_NAME),
      ]
  if config["steps"]["alignment"] or config["input_format"] == "bam":
      targets["split_concat_bam"]=[
        #split
#        expand(os.path.normpath(OUTPUT_DIR + "/split_chr/{sample_name}/{sample_name}_chr_{chr_number}.bam"),chr_number=CHR_NUMBER,sample_name=SAMPLE_NAME),
#        expand(os.path.normpath(OUTPUT_DIR + "/split_chr/{sample_name}/{sample_name}_chr_{chr_number}.bam.bai"),chr_number=CHR_NUMBER,sample_name=SAMPLE_NAME),
        #reconcat
#        expand(os.path.normpath(OUTPUT_DIR + "/reconcat/{sample_name}/{sample_name}_reconcat.bam"),sample_name=SAMPLE_NAME),
#        expand(os.path.normpath(OUTPUT_DIR + "/reconcat/{sample_name}/{sample_name}_sorted.bam"),sample_name=SAMPLE_NAME),
#        expand(os.path.normpath(OUTPUT_DIR + "/reconcat/{sample_name}/{sample_name}_sorted.bam.bai"),sample_name=SAMPLE_NAME)
      ]
      #QUALITY-CONTROL
      targets["bam_qc"]=[
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
      ]
      targets["fastq_qc"]=[
        #bam_to_fastq
        expand(os.path.normpath(OUTPUT_DIR + "/Fastq/{sample_name}.fastq.gz"), sample_name=SAMPLE_NAME),
        #fastqc
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/fastq_QC/fastqc/{sample_name}/{sample_name}_fastqc.html"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/fastq_QC/fastqc/{sample_name}/{sample_name}_fastqc.zip"), sample_name=SAMPLE_NAME),
        #fastq_screen
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/fastq_QC/fastq_screen/{sample_name}/{sample_name}_screen.txt"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/fastq_QC/fastq_screen/{sample_name}/{sample_name}_screen.html"), sample_name=SAMPLE_NAME),
      ]
      if config["basecalling_mode"] == "methylation":
        targets["methylation_qc"]=[
            #meth_QC_ratio_fwd_rev
            expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/methylation_QC/per_pos_per_strand_gam_ratio_{sample_name}_chr{chr_number}_{meth_type}_mqc.png"),chr_number=CHR_NUMBER,sample_name=SAMPLE_NAME,meth_type=METH_TYPE),
            #QC_barplot_CG
            expand(os.path.normpath(OUTPUT_DIR + "/Quality_Control/methylation_QC/barplot_methylated_CG_{meth_type}_mqc.png"),meth_type=METH_TYPE)
        ]
      targets["multiqc"]=[
        #multiqc
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/multiqc_report.html")
      ]
      
  #METHYLATION
  if config["steps"]["differential_methylation_sample"] or config["steps"]["differential_methylation_condition"]:
      targets["methylation"]=[
        #modkit separate mod
#        expand(os.path.normpath(OUTPUT_DIR + "/separate_mod/{sample_name}/{meth_type}/{sample_name}_chr_{chr_number}_{meth_type}.bam"),chr_number=CHR_NUMBER,sample_name=SAMPLE_NAME,meth_type=METH_TYPE),
#        expand(os.path.normpath(OUTPUT_DIR + "/separate_mod/{sample_name}/{meth_type}/{sample_name}_chr_{chr_number}_{meth_type}.bam.bai"),chr_number=CHR_NUMBER,sample_name=SAMPLE_NAME,meth_type=METH_TYPE),
        #modkit pileup uncombined strand
#        expand(os.path.normpath(OUTPUT_DIR + "/bed_uncombined_strands/{sample_name}/{meth_type}/{sample_name}_chr{chr_number}_{meth_type}_uncomb.bed"),chr_number=CHR_NUMBER,sample_name=SAMPLE_NAME,meth_type=METH_TYPE),
        #modkit pileup combined strand
        #expand(os.path.normpath(OUTPUT_DIR + "/bed_combined_strands/{sample_name}/{meth_type}/{sample_name}_chr{chr_number}_{meth_type}_comb.bed"),chr_number=CHR_NUMBER,sample_name=SAMPLE_NAME,meth_type=METH_TYPE),
        #expand(os.path.normpath(OUTPUT_DIR + "/bed_combined_strands/{sample_name}/{meth_type}/{sample_name}_chr{chr_number}_{meth_type}_comb_format.bed"),chr_number=CHR_NUMBER,sample_name=SAMPLE_NAME,meth_type=METH_TYPE),
        #QC methylation: prend une des sorties du QC general
        #expand(output du QC)
        #split bed
        #expand(os.path.normpath(OUTPUT_DIR + "/Methylation_Analysis/bed_files/chr{chr_number}/{meth_type}/{sample_name}-chr{chr_number}-{meth_type}-{strands}.bed"),chr_number=CHR_NUMBER,sample_name=SAMPLE_NAME,meth_type=METH_TYPE,strands=STRAND),
        #split Alu
#        expand(os.path.normpath(OUTPUT_DIR + "/resources/Alu/refseq.bed_Alu_chr{chr_number}_{strands}.txt"),chr_number=CHR_NUMBER,strands=STRAND),
        #split Transcript
#        expand(os.path.normpath(OUTPUT_DIR + "/resources/Transcript/refseq.bed_Transcript_chr{chr_number}_{strands}.txt"),chr_number=CHR_NUMBER,strands=STRAND),
        #split CpG
#        expand(os.path.normpath(OUTPUT_DIR + "/resources/CpG/cpgi.bed_CpG_chr{chr_number}.txt"),chr_number=CHR_NUMBER)
      ]
  if config["steps"]["differential_methylation_sample"]:
      targets["differential_methylation_sample"]=[
        #mean_table
#        expand(os.path.normpath(OUTPUT_DIR + "/Methylation_Analysis/{ref_types}/{meth_type}/{pair_methyl}/mean_table_{pair_methyl}_chr{chr_number}_{strands}.tsv"), pair_methyl=PAIR_METHYL, chr_number=CHR_NUMBER, meth_type=METH_TYPE, strands=STRAND, ref_types=REF_TYPE),
        #dmr_stat
        expand(os.path.normpath(OUTPUT_DIR + "/Methylation_Analysis/{ref_types}/{meth_type}/{pair_methyl}/mean_table_{pair_methyl}_all_chr.tsv"), pair_methyl = PAIR_METHYL, meth_type = METH_TYPE, ref_types = REF_TYPE),
        expand(os.path.normpath(OUTPUT_DIR + "/Methylation_Analysis/{ref_types}/{meth_type}/{pair_methyl}/meth_res_table_{pair_methyl}_all_chr.csv"), pair_methyl = PAIR_METHYL, meth_type = METH_TYPE, ref_types = REF_TYPE),
      ]
  if config["steps"]["differential_methylation_condition"]:
      targets["differential_methylation_condition"]=[
        #mean_table
#        expand(os.path.normpath(OUTPUT_DIR + "/Methylation_Analysis/{ref_types}/{meth_type}/cases_vs_controls/mean_table_cases_vs_controls_chr{chr_number}_{strands}.tsv"), chr_number = CHR_NUMBER, meth_type = METH_TYPE, strands = STRAND, ref_types = REF_TYPE),
        #dmr_stat
        expand(os.path.normpath(OUTPUT_DIR + "/Methylation_Analysis/{ref_types}/{meth_type}/cases_vs_controls/mean_table_cases_vs_controls_all_chr.tsv"), meth_type = METH_TYPE, ref_types = REF_TYPE),
        expand(os.path.normpath(OUTPUT_DIR + "/Methylation_Analysis/{ref_types}/{meth_type}/cases_vs_controls/meth_res_table_cases_vs_controls_all_chr.csv"), meth_type = METH_TYPE, ref_types = REF_TYPE),
      ]

  #SINGLE NUCLEOTIDE VARIANT CALLING
  if config["steps"]["snv_calling"]:
      #GERMLINE
      if config["variant_calling_mode"] == "germline":
        targets["snv_calling"]=[
            #clair3
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/clair3/{clair3_model}/{sample_name}/{sample_name}_merge_output.vcf.gz"), sample_name=SAMPLE_NAME, clair3_model=NAME_CLAIR3_MODEL),
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/clair3/{clair3_model}/{sample_name}/{sample_name}_merge_output.vcf.gz.tbi"), sample_name=SAMPLE_NAME, clair3_model=NAME_CLAIR3_MODEL),
            #pepper_margin_deepvariant
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/pepper_margin_deepvariant/{sample_name}/{sample_name}.vcf.gz"), sample_name=SAMPLE_NAME),
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/pepper_margin_deepvariant/{sample_name}/{sample_name}.visual_report.html"), sample_name=SAMPLE_NAME),
        ]
        targets["snv_annotation"]=[]
        #clair3 & snpEff & snpSift
        if "genome_snpEff" in config["references"] and config["references"]["genome_snpEff"] != "": targets["snv_annotation"].append(expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/clair3/{clair3_model}/{sample_name}/snpEff/{sample_name}{compl}_annotated.vcf"), sample_name=SAMPLE_NAME, clair3_model=NAME_CLAIR3_MODEL, compl="_merge_output")),
        if "genome_snpEff" in config["references"] and config["references"]["genome_snpEff"] != "": targets["snv_annotation"].append(expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/clair3/{clair3_model}/{sample_name}/snpEff/{sample_name}{compl}_annotated_stats.csv"), sample_name=SAMPLE_NAME, clair3_model=NAME_CLAIR3_MODEL, compl="_merge_output")),
        if "genome_snpEff" in config["references"] and config["references"]["genome_snpEff"] != "": targets["snv_annotation"].append(expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/clair3/{clair3_model}/{sample_name}/snpEff/{sample_name}{compl}_annotated_summary.html"), sample_name=SAMPLE_NAME, clair3_model=NAME_CLAIR3_MODEL, compl="_merge_output")),
        if "dbnsfp" in config["references"] and config["references"]["dbnsfp"] != "": targets["snv_annotation"].append(expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/clair3/{clair3_model}/{sample_name}/snpEff/{sample_name}{compl}" + SNPEFF_SUFFIX + "_dbnsfp.vcf"), sample_name=SAMPLE_NAME, clair3_model=NAME_CLAIR3_MODEL, compl="_merge_output")),
        if "clinvar" in config["references"] and config["references"]["clinvar"] != "": targets["snv_annotation"].append(expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/clair3/{clair3_model}/{sample_name}/snpEff/{sample_name}{compl}" + SNPEFF_SUFFIX + DBNSFP_SUFFIX + "_clinvar.vcf"), sample_name=SAMPLE_NAME, clair3_model=NAME_CLAIR3_MODEL, compl="_merge_output")),
        targets["snv_annotation"].append(expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/clair3/{clair3_model}/{sample_name}/snpEff/{sample_name}{compl}" + SNPEFF_SUFFIX + DBNSFP_SUFFIX + CLINVAR_SUFFIX + "_{filter}.vcf"), sample_name=SAMPLE_NAME, clair3_model=NAME_CLAIR3_MODEL, filter=SNPSIFT_FILTERS_NAMES, compl="_merge_output")),
        #pepper_margin_deepvariant & snpEff & snpSift
        if "genome_snpEff" in config["references"] and config["references"]["genome_snpEff"] != "": targets["snv_annotation"].append(expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/pepper_margin_deepvariant/{sample_name}/snpEff/{sample_name}_annotated.vcf"), sample_name=SAMPLE_NAME)),
        if "genome_snpEff" in config["references"] and config["references"]["genome_snpEff"] != "": targets["snv_annotation"].append(expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/pepper_margin_deepvariant/{sample_name}/snpEff/{sample_name}_annotated_stats.csv"), sample_name=SAMPLE_NAME)),
        if "genome_snpEff" in config["references"] and config["references"]["genome_snpEff"] != "": targets["snv_annotation"].append(expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/pepper_margin_deepvariant/{sample_name}/snpEff/{sample_name}_annotated_summary.html"), sample_name=SAMPLE_NAME)),
        if "dbnsfp" in config["references"] and config["references"]["dbnsfp"] != "": targets["snv_annotation"].append(expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/pepper_margin_deepvariant/{sample_name}/snpEff/{sample_name}" + SNPEFF_SUFFIX + "_dbnsfp.vcf"), sample_name=SAMPLE_NAME)),
        if "clinvar" in config["references"] and config["references"]["clinvar"] != "": targets["snv_annotation"].append(expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/pepper_margin_deepvariant/{sample_name}/snpEff/{sample_name}" + SNPEFF_SUFFIX + DBNSFP_SUFFIX + "_clinvar.vcf"), sample_name=SAMPLE_NAME)),
        targets["snv_annotation"].append(expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/pepper_margin_deepvariant/{sample_name}/snpEff/{sample_name}" + SNPEFF_SUFFIX + DBNSFP_SUFFIX + CLINVAR_SUFFIX + "_{filter}.vcf"), sample_name=SAMPLE_NAME, filter=SNPSIFT_FILTERS_NAMES))
        
        targets["snv_graphs"]=[
            #clair3 & maftools
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/clair3/{clair3_model}/{sample_name}/maftools/{sample_name}{compl}" + SNPEFF_SUFFIX + DBNSFP_SUFFIX + CLINVAR_SUFFIX + "_{filter}.maf"), sample_name=SAMPLE_NAME, compl="_merge_output", clair3_model=NAME_CLAIR3_MODEL, filter=SNPSIFT_FILTERS_NAMES),
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/clair3/{clair3_model}/{sample_name}/maftools/{sample_name}{compl}" + SNPEFF_SUFFIX + DBNSFP_SUFFIX + CLINVAR_SUFFIX + "_{filter}_mafSummary_plot.pdf"), sample_name=SAMPLE_NAME, clair3_model=NAME_CLAIR3_MODEL, compl="_merge_output", filter=SNPSIFT_FILTERS_NAMES),
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/clair3/{clair3_model}/{sample_name}/maftools/{sample_name}{compl}" + SNPEFF_SUFFIX + DBNSFP_SUFFIX + CLINVAR_SUFFIX + "_{filter}_mafbarplot.pdf"), sample_name=SAMPLE_NAME, clair3_model=NAME_CLAIR3_MODEL, compl="_merge_output", filter=SNPSIFT_FILTERS_NAMES),
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/clair3/{clair3_model}/{sample_name}/maftools/{sample_name}{compl}" + SNPEFF_SUFFIX + DBNSFP_SUFFIX + CLINVAR_SUFFIX + "_{filter}_oncoplot.pdf"), sample_name=SAMPLE_NAME, clair3_model=NAME_CLAIR3_MODEL, compl="_merge_output", filter=SNPSIFT_FILTERS_NAMES),
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/clair3/{clair3_model}/{sample_name}/maftools/{sample_name}{compl}" + SNPEFF_SUFFIX + DBNSFP_SUFFIX + CLINVAR_SUFFIX + "_{filter}_titv_Transitions_Transversions.pdf"), sample_name=SAMPLE_NAME, clair3_model=NAME_CLAIR3_MODEL, compl="_merge_output", filter=SNPSIFT_FILTERS_NAMES),
            #pepper_margin_deepvariant & maftools
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/pepper_margin_deepvariant/{sample_name}/maftools/{sample_name}" + SNPEFF_SUFFIX + DBNSFP_SUFFIX + CLINVAR_SUFFIX + "_{filter}.maf"), sample_name=SAMPLE_NAME, filter=SNPSIFT_FILTERS_NAMES),
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/pepper_margin_deepvariant/{sample_name}/maftools/{sample_name}" + SNPEFF_SUFFIX + DBNSFP_SUFFIX + CLINVAR_SUFFIX + "_{filter}_mafSummary_plot.pdf"), sample_name=SAMPLE_NAME, filter=SNPSIFT_FILTERS_NAMES),
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/pepper_margin_deepvariant/{sample_name}/maftools/{sample_name}" + SNPEFF_SUFFIX + DBNSFP_SUFFIX + CLINVAR_SUFFIX + "_{filter}_mafbarplot.pdf"), sample_name=SAMPLE_NAME, filter=SNPSIFT_FILTERS_NAMES),
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/pepper_margin_deepvariant/{sample_name}/maftools/{sample_name}" + SNPEFF_SUFFIX + DBNSFP_SUFFIX + CLINVAR_SUFFIX + "_{filter}_oncoplot.pdf"), sample_name=SAMPLE_NAME, filter=SNPSIFT_FILTERS_NAMES),
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/pepper_margin_deepvariant/{sample_name}/maftools/{sample_name}" + SNPEFF_SUFFIX + DBNSFP_SUFFIX + CLINVAR_SUFFIX + "_{filter}_titv_Transitions_Transversions.pdf"), sample_name=SAMPLE_NAME, filter=SNPSIFT_FILTERS_NAMES)
        ]
      #SOMATIC
      if config["variant_calling_mode"] == "somatic":
        targets["snv_calling"]=[
            #clairs
            expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Somatic/clairs/{pair_somatic}/{pair_somatic}{compl}.vcf.gz"), pair_somatic=PAIR_SOMATIC, compl=["_snv","_indel"])
        ]
        targets["snv_annotation"]=[]
        #clairs & snpEff & snpSift
        if "genome_snpEff" in config["references"] and config["references"]["genome_snpEff"] != "": targets["snv_annotation"].append(expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Somatic/clairs/{pair_somatic}/snpEff/{pair_somatic}{compl}_annotated.vcf"), pair_somatic=PAIR_SOMATIC, compl=["_snv","_indel"])),
        if "genome_snpEff" in config["references"] and config["references"]["genome_snpEff"] != "": targets["snv_annotation"].append(expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Somatic/clairs/{pair_somatic}/snpEff/{pair_somatic}{compl}_annotated_stats.csv"), pair_somatic=PAIR_SOMATIC, compl=["_snv","_indel"])),
        if "genome_snpEff" in config["references"] and config["references"]["genome_snpEff"] != "": targets["snv_annotation"].append(expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Somatic/clairs/{pair_somatic}/snpEff/{pair_somatic}{compl}_annotated_summary.html"), pair_somatic=PAIR_SOMATIC, compl=["_snv","_indel"])),
        if "dbnsfp" in config["references"] and config["references"]["dbnsfp"] != "": targets["snv_annotation"].append(expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Somatic/clairs/{pair_somatic}/snpEff/{pair_somatic}{compl}" + SNPEFF_SUFFIX + "_dbnsfp.vcf"), pair_somatic=PAIR_SOMATIC, compl=["_snv","_indel"])),
        if "clinvar" in config["references"] and config["references"]["clinvar"] != "": targets["snv_annotation"].append(expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Somatic/clairs/{pair_somatic}/snpEff/{pair_somatic}{compl}" + SNPEFF_SUFFIX + DBNSFP_SUFFIX + "_clinvar.vcf"), pair_somatic=PAIR_SOMATIC, compl=["_snv","_indel"])),
        targets["snv_annotation"].append(expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Somatic/clairs/{pair_somatic}/snpEff/{pair_somatic}{compl}" + SNPEFF_SUFFIX + DBNSFP_SUFFIX + CLINVAR_SUFFIX + "_{filter}.vcf"), pair_somatic=PAIR_SOMATIC, filter=SNPSIFT_FILTERS_NAMES, compl=["_snv","_indel"]))

  #PHASING
  if config["steps"]["phasing"]:
      #GERMLINE
      if config["variant_calling_mode"] == "germline":
        targets["phasing"]=[
          #clair3
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/clair3/{clair3_model}/{sample_name}/whatshap/{sample_name}{compl}_phased.vcf.gz"), sample_name=SAMPLE_NAME, clair3_model=NAME_CLAIR3_MODEL, compl="_merge_output"),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/clair3/{clair3_model}/{sample_name}/whatshap/{sample_name}{compl}_phasing_stats.txt"), sample_name=SAMPLE_NAME, clair3_model=NAME_CLAIR3_MODEL, compl="_merge_output"),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/clair3/{clair3_model}/{sample_name}/whatshap/{sample_name}{compl}_phasing_stats.tsv"), sample_name=SAMPLE_NAME, clair3_model=NAME_CLAIR3_MODEL, compl="_merge_output"),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/clair3/{clair3_model}/{sample_name}/whatshap/{sample_name}{compl}_phasing_haplotype_blocks.tsv"), sample_name=SAMPLE_NAME, clair3_model=NAME_CLAIR3_MODEL, compl="_merge_output"),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/clair3/{clair3_model}/{sample_name}/whatshap/{sample_name}{compl}_phasing_haplotype_blocks.gtf"), sample_name=SAMPLE_NAME, clair3_model=NAME_CLAIR3_MODEL, compl="_merge_output"),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/clair3/{clair3_model}/{sample_name}/whatshap/{sample_name}{compl}_phased.vcf.gz.tbi"), sample_name=SAMPLE_NAME, clair3_model=NAME_CLAIR3_MODEL, compl="_merge_output"),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/clair3/{clair3_model}/{sample_name}/whatshap/{sample_name}{compl}_phasing_haplotag_list.tsv"), sample_name=SAMPLE_NAME, clair3_model=NAME_CLAIR3_MODEL, compl="_merge_output"),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/clair3/{clair3_model}/{sample_name}/whatshap/{sample_name}{compl}_haplotagged.bam"), sample_name=SAMPLE_NAME, clair3_model=NAME_CLAIR3_MODEL, compl="_merge_output"),
          #pepper_margin_deepvariant
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/pepper_margin_deepvariant/{sample_name}/whatshap/{sample_name}{compl}_phased.vcf.gz"), sample_name=SAMPLE_NAME, compl=""),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/pepper_margin_deepvariant/{sample_name}/whatshap/{sample_name}{compl}_phasing_stats.txt"), sample_name=SAMPLE_NAME, compl=""),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/pepper_margin_deepvariant/{sample_name}/whatshap/{sample_name}{compl}_phasing_stats.tsv"), sample_name=SAMPLE_NAME, compl=""),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/pepper_margin_deepvariant/{sample_name}/whatshap/{sample_name}{compl}_phasing_haplotype_blocks.tsv"), sample_name=SAMPLE_NAME, compl=""),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/pepper_margin_deepvariant/{sample_name}/whatshap/{sample_name}{compl}_phasing_haplotype_blocks.gtf"), sample_name=SAMPLE_NAME, compl=""),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/pepper_margin_deepvariant/{sample_name}/whatshap/{sample_name}{compl}_phased.vcf.gz.tbi"), sample_name=SAMPLE_NAME, compl=""),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/pepper_margin_deepvariant/{sample_name}/whatshap/{sample_name}{compl}_phasing_haplotag_list.tsv"), sample_name=SAMPLE_NAME, compl=""),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Germline/pepper_margin_deepvariant/{sample_name}/whatshap/{sample_name}{compl}_haplotagged.bam"), sample_name=SAMPLE_NAME, compl=""),
        ]
      #SOMATIC
      if config["variant_calling_mode"] == "somatic":
        targets["phasing"]=[
          #clairs & whatshap
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Somatic/clairs/{pair_somatic}/whatshap/{pair_somatic}{compl}_phased.vcf.gz"), pair_somatic=PAIR_SOMATIC, compl=["_snv","_indel"]),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Somatic/clairs/{pair_somatic}/whatshap/{pair_somatic}{compl}_phasing_stats.txt"), pair_somatic=PAIR_SOMATIC, compl=["_snv","_indel"]),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Somatic/clairs/{pair_somatic}/whatshap/{pair_somatic}{compl}_phasing_stats.tsv"), pair_somatic=PAIR_SOMATIC, compl=["_snv","_indel"]),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Somatic/clairs/{pair_somatic}/whatshap/{pair_somatic}{compl}_phasing_haplotype_blocks.tsv"), pair_somatic=PAIR_SOMATIC, compl=["_snv","_indel"]),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Somatic/clairs/{pair_somatic}/whatshap/{pair_somatic}{compl}_phasing_haplotype_blocks.gtf"), pair_somatic=PAIR_SOMATIC, compl=["_snv","_indel"]),
          expand(os.path.normpath(OUTPUT_DIR + "/SNV_Calling/Somatic/clairs/{pair_somatic}/whatshap/{pair_somatic}{compl}_phased.vcf.gz.tbi"), pair_somatic=PAIR_SOMATIC, compl=["_snv","_indel"])
        ]
  #STRUCTURAL VARIANT CALLING
  if config["steps"]["sv_calling"]:
      #GERMLINE
      if config["variant_calling_mode"] == "germline":
        targets["germline_sv_calling"]=[
            #sniffles
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Germline/sniffles/{sample_name}/{sample_name}_SV.vcf"), sample_name=SAMPLE_NAME),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Germline/sniffles/{sample_name}/{sample_name}_SV.snf"), sample_name=SAMPLE_NAME),
            #cuteSV
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Germline/cuteSV/{sample_name}/{sample_name}_SV.vcf"), sample_name=SAMPLE_NAME)
        ]
        if config["references"]["genome_annotsv"] != "" :
            targets["sv_annotation"]=[
                #sniffles & annotSV
                expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Germline/sniffles/{sample_name}/AnnotSV/{sample_name}{compl_SV}.annotated.vcf"), sample_name=SAMPLE_NAME, compl_SV=["_SV"]),
                expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Germline/sniffles/{sample_name}/AnnotSV/{sample_name}{compl_SV}.annotated.tsv"), sample_name=SAMPLE_NAME, compl_SV=["_SV"]),
                #cuteSV & annotSV
                expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Germline/cuteSV/{sample_name}/AnnotSV/{sample_name}{compl_SV}.annotated.tsv"), sample_name=SAMPLE_NAME, compl_SV=["_SV"]),
                expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Germline/cuteSV/{sample_name}/AnnotSV/{sample_name}{compl_SV}.annotated.vcf"), sample_name=SAMPLE_NAME, compl_SV=["_SV"])
            ]
        targets["germline_sv_graphs"]=[
            #sniffles & sniffles2_plot
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Germline/sniffles/{sample_name}/sniffles2_plot/{sample_name}.png"), sample_name=SAMPLE_NAME),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Germline/sniffles/{sample_name}/sniffles2_plot/del_ins_genotype.jpg"), sample_name=SAMPLE_NAME),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Germline/sniffles/{sample_name}/sniffles2_plot/del_ins_type_size.jpg"), sample_name=SAMPLE_NAME),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Germline/sniffles/{sample_name}/sniffles2_plot/dup_inv_type_size.jpg"), sample_name=SAMPLE_NAME),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Germline/sniffles/{sample_name}/sniffles2_plot/inv_dup_genotype.jpg"), sample_name=SAMPLE_NAME),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Germline/sniffles/{sample_name}/sniffles2_plot/length_variant.jpg"), sample_name=SAMPLE_NAME),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Germline/sniffles/{sample_name}/sniffles2_plot/variant_count.jpg"), sample_name=SAMPLE_NAME),
            #cuteSV & sniffles2_plot
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Germline/cuteSV/{sample_name}/sniffles2_plot/{sample_name}.png"), sample_name=SAMPLE_NAME),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Germline/cuteSV/{sample_name}/sniffles2_plot/del_ins_genotype.jpg"), sample_name=SAMPLE_NAME),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Germline/cuteSV/{sample_name}/sniffles2_plot/del_ins_type_size.jpg"), sample_name=SAMPLE_NAME),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Germline/cuteSV/{sample_name}/sniffles2_plot/dup_inv_type_size.jpg"), sample_name=SAMPLE_NAME),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Germline/cuteSV/{sample_name}/sniffles2_plot/inv_dup_genotype.jpg"), sample_name=SAMPLE_NAME),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Germline/cuteSV/{sample_name}/sniffles2_plot/length_variant.jpg"), sample_name=SAMPLE_NAME),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Germline/cuteSV/{sample_name}/sniffles2_plot/variant_count.jpg"), sample_name=SAMPLE_NAME)
        ]
        if len(SAMPLE_NAME) > 1 :
            targets["germline_sv_calling"].append(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Germline/sniffles/all_samples/all_samples_SV.vcf")),
            if config["references"]["genome_annotsv"] != "" : targets["sv_annotation"].append(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Germline/sniffles/all_samples/AnnotSV/all_samples_SV.annotated.vcf")),
            if config["references"]["genome_annotsv"] != "" : targets["sv_annotation"].append(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Germline/sniffles/all_samples/AnnotSV/all_samples_SV.annotated.tsv"))
      
      #SOMATIC
      if config["variant_calling_mode"] == "somatic":
        targets["somatic_sv_calling"]=[
            #nanomonsv
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/{sample_name}.deletion.sorted.bed.gz"), zip, pair_somatic=PAIR_SOMATIC, sample_name=list(dict.fromkeys(SAMPLE_NAME))), #"list(dict.fromkeys(" removes duplicate SAMPLE_NAME du to BATCH_NAME duplication  in  methylation part; needed beacause of "zip" using for expand()
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/{sample_name}.deletion.sorted.bed.gz.tbi"), zip, pair_somatic=PAIR_SOMATIC, sample_name=list(dict.fromkeys(SAMPLE_NAME))),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/{sample_name}.insertion.sorted.bed.gz"), zip, pair_somatic=PAIR_SOMATIC, sample_name=list(dict.fromkeys(SAMPLE_NAME))),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/{sample_name}.insertion.sorted.bed.gz.tbi"), zip, pair_somatic=PAIR_SOMATIC, sample_name=list(dict.fromkeys(SAMPLE_NAME))),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/{sample_name}.rearrangement.sorted.bedpe.gz"), zip, pair_somatic=PAIR_SOMATIC, sample_name=list(dict.fromkeys(SAMPLE_NAME))),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/{sample_name}.rearrangement.sorted.bedpe.gz.tbi"), zip, pair_somatic=PAIR_SOMATIC, sample_name=list(dict.fromkeys(SAMPLE_NAME))),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/{sample_name}.bp_info.sorted.bed.gz"), zip, pair_somatic=PAIR_SOMATIC, sample_name=list(dict.fromkeys(SAMPLE_NAME))),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/{sample_name}.bp_info.sorted.bed.gz.tbi"), zip, pair_somatic=PAIR_SOMATIC, sample_name=list(dict.fromkeys(SAMPLE_NAME))),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/{pair_somatic}.nanomonsv.result.txt"), pair_somatic=PAIR_SOMATIC),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/{pair_somatic}.nanomonsv.result.vcf"), pair_somatic=PAIR_SOMATIC),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/{pair_somatic}.nanomonsv.sbnd.result.txt"), pair_somatic=PAIR_SOMATIC),
            expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/{pair_somatic}.nanomonsv.supporting_read.txt"), pair_somatic=PAIR_SOMATIC)
        ]
        #nanomonsv & AnnotSV
        if config["references"]["genome_annotsv"] != "" : 
            targets["somatic_sv_annotation"]=[
                expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/AnnotSV/{pair_somatic}{compl_SV}.annotated.vcf"), pair_somatic=PAIR_SOMATIC, compl_SV=[".nanomonsv.result"]),
                expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/AnnotSV/{pair_somatic}{compl_SV}.annotated.tsv"), pair_somatic=PAIR_SOMATIC, compl_SV=[".nanomonsv.result"])
            ]
        #if config["references"]["genome_annotsv"] == "hg38":
            #nanomonsv classifier
            #targets["somatic_sv_calling"].append(expand(os.path.normpath(OUTPUT_DIR + "/SV_Calling/Somatic/nanomonsv/{pair_somatic}/{pair_somatic}.nanomonsv.insert_classify.txt"), pair_somatic=PAIR_SOMATIC))
  
  #COPY NUMBER VARIANT CALLING
  if config["steps"]["cnv_calling"]:
    targets["cnv_calling"]=[
        #spectre
        expand(os.path.normpath(OUTPUT_DIR + "/CNV_Calling/spectre/{sample_name}/{sample_name}.vcf.gz"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/CNV_Calling/spectre/{sample_name}/{sample_name}_cnv.bed.gz"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/CNV_Calling/spectre/{sample_name}/{sample_name}_cnv.bed.gz.tbi"), sample_name=SAMPLE_NAME),
        expand(os.path.normpath(OUTPUT_DIR + "/CNV_Calling/spectre/{sample_name}/img/{sample_name}_plot_cnv_chr_{chromos}.png"), sample_name=SAMPLE_NAME, chromos=[x for x in CHR_NUMBER if x not in ["MT","M"]])
    ]

  return targets
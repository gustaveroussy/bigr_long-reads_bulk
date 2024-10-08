"""
##########################################################################
These rules make the control-quality of the alignment
##########################################################################
"""

"""
This rule makes the qualimap QC
"""

rule qualimap:
    input:
        bam_file = os.path.normpath(OUTPUT_DIR + "/tmp/reconcat/{sample_name}/{sample_name}_sorted.bam"),
        index = os.path.normpath(OUTPUT_DIR + "/tmp/reconcat/{sample_name}/{sample_name}_sorted.bam.bai")
    output:
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/qualimap/{sample_name}/qualimapReport.html"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/qualimap/{sample_name}/raw_data_qualimapReport/genome_fraction_coverage.txt"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/qualimap/{sample_name}/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/qualimap/{sample_name}/genome_results.txt"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/qualimap/{sample_name}/raw_data_qualimapReport/coverage_histogram.txt")
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
        qualimap bamqc -bam {input.bam_file} -outdir {OUTPUT_DIR}/Quality_Control/bam_QC/qualimap/{wildcards.sample_name}/ --paint-chromosome-limits -nt {threads} --java-mem-size=${{res}}G
        
        """


"""
This rule makes the mosdepth QC
"""

rule mosdepth:
    input:
        bam_file = os.path.normpath(OUTPUT_DIR + "/tmp/reconcat/{sample_name}/{sample_name}_sorted.bam"),
        index = os.path.normpath(OUTPUT_DIR + "/tmp/reconcat/{sample_name}/{sample_name}_sorted.bam.bai")
    output:
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/mosdepth/{sample_name}/{sample_name}_wgs_mode.mosdepth.global.dist.txt"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/mosdepth/{sample_name}/{sample_name}_wgs_mode.mosdepth.region.dist.txt"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/mosdepth/{sample_name}/{sample_name}_wgs_mode.mosdepth.summary.txt"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/mosdepth/{sample_name}/{sample_name}_wgs_mode.regions.bed.gz"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/mosdepth/{sample_name}/{sample_name}_wgs_mode.regions.bed.gz.csi")
    threads:
        3
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 720)
    conda:
        CONDA_ENV_MOSDEPTH
    shell:
        """
        cd {OUTPUT_DIR}/Quality_Control/bam_QC/mosdepth/{wildcards.sample_name}/
        mosdepth -n --fast-mode --by 500 -t {threads} {wildcards.sample_name}_wgs_mode {input.bam_file} && \
        mosdepth -t {threads} {wildcards.sample_name} {input.bam_file}
        
        """


"""
This rule makes the nanoplot QC
"""

rule nanoplot_bam:
    input:
        bam_file = os.path.normpath(OUTPUT_DIR + "/tmp/reconcat/{sample_name}/{sample_name}_sorted.bam"),
        index = os.path.normpath(OUTPUT_DIR + "/tmp/reconcat/{sample_name}/{sample_name}_sorted.bam.bai")
    output:
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_AlignedReadlengthvsSequencedReadLength_dot.html"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_AlignedReadlengthvsSequencedReadLength_dot.png"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_AlignedReadlengthvsSequencedReadLength_kde.html"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_AlignedReadlengthvsSequencedReadLength_kde.png"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_LengthvsQualityScatterPlot_dot.html"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_LengthvsQualityScatterPlot_dot.png"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_LengthvsQualityScatterPlot_kde.html"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_LengthvsQualityScatterPlot_kde.png"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_MappingQualityvsAverageBaseQuality_dot.html"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_MappingQualityvsAverageBaseQuality_dot.png"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_MappingQualityvsAverageBaseQuality_kde.html"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_MappingQualityvsAverageBaseQuality_kde.png"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_MappingQualityvsReadLength_dot.html"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_MappingQualityvsReadLength_dot.png"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_MappingQualityvsReadLength_kde.html"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_MappingQualityvsReadLength_kde.png"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_NanoPlot-report.html"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_NanoStats.txt"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_Non_weightedHistogramReadlength.html"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_Non_weightedHistogramReadlength.png"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_Non_weightedLogTransformed_HistogramReadlength.html"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_Non_weightedLogTransformed_HistogramReadlength.png"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_PercentIdentityHistogramDynamic_Histogram_percent_identity.html"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_PercentIdentityHistogramDynamic_Histogram_percent_identity.png"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_PercentIdentityvsAlignedReadLength_dot.html"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_PercentIdentityvsAlignedReadLength_dot.png"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_PercentIdentityvsAlignedReadLength_kde.html"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_PercentIdentityvsAlignedReadLength_kde.png"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_PercentIdentityvsAverageBaseQuality_dot.html"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_PercentIdentityvsAverageBaseQuality_dot.png"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_PercentIdentityvsAverageBaseQuality_kde.html"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_PercentIdentityvsAverageBaseQuality_kde.png"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_WeightedHistogramReadlength.html"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_WeightedHistogramReadlength.png"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_WeightedLogTransformed_HistogramReadlength.html"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_WeightedLogTransformed_HistogramReadlength.png"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_Yield_By_Length.html"),
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/nanoplot/{sample_name}/{sample_name}_Yield_By_Length.png")
    threads:
        10
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 40960),
        time_min = (lambda wildcards, attempt: attempt * 720)
    shell:
        """
        TMP_DIR=$(mktemp -d -t lr_pipeline-XXXXXXXXXX) && \
        singularity exec --contain -B {OUTPUT_DIR} -B ${{TMP_DIR}}:/tmp -B ${{TMP_DIR}}:${{TMPDIR}} {SING_ENV_NANOPLOT} \
        NanoPlot --threads {threads} --outdir {OUTPUT_DIR}/Quality_Control/bam_QC/nanoplot/{wildcards.sample_name}/ --prefix {wildcards.sample_name}_ --N50 --tsv_stats --info_in_report --bam {input.bam_file}

        """


"""
These rules make various samtools stats files: stats, flagstat, idxstats and coverage
"""

rule samtools_stats:
    input:
        bam_file = os.path.normpath(OUTPUT_DIR + "/tmp/reconcat/{sample_name}/{sample_name}_sorted.bam")
    output:
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/samtools/{sample_name}/{sample_name}_stats.txt")
    threads:
        2
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 4096),
        time_min = (lambda wildcards, attempt: attempt * 60)
    conda:
        CONDA_ENV_SAMTOOLS
    shell:
        """
        samtools stats --threads {threads} {input.bam_file} > {output}
        """

rule samtools_flagstat:
    input:
        bam_file = os.path.normpath(OUTPUT_DIR + "/tmp/reconcat/{sample_name}/{sample_name}_sorted.bam")
    output:
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/samtools/{sample_name}/{sample_name}_flagstat.txt")
    threads:
        2
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 4096),
        time_min = (lambda wildcards, attempt: attempt * 60)
    conda:
        CONDA_ENV_SAMTOOLS
    shell:
        """
        samtools flagstat --threads {threads} {input.bam_file} > {output}
        """

rule samtools_idxstats:
    input:
        bam_file = os.path.normpath(OUTPUT_DIR + "/tmp/reconcat/{sample_name}/{sample_name}_sorted.bam"),
        index_file = os.path.normpath(OUTPUT_DIR + "/tmp/reconcat/{sample_name}/{sample_name}_sorted.bam.bai")
    output:
        os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/samtools/{sample_name}/{sample_name}_idxstats.txt")
    threads:
        2
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 1024),
        time_min = (lambda wildcards, attempt: attempt * 60)
    conda:
        CONDA_ENV_SAMTOOLS
    shell:
        """
        samtools idxstats {input.bam_file} > {output}
        """
        
rule samtools_coverage:
    input:
        bam_file = os.path.normpath(OUTPUT_DIR + "/tmp/reconcat/{sample_name}/{sample_name}_sorted.bam")
    output:
        all_chr_coverage = temp(os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/samtools/{sample_name}/{sample_name}_all_chr_coverage.txt")),
        coverage = os.path.normpath(OUTPUT_DIR + "/Quality_Control/bam_QC/samtools/{sample_name}/{sample_name}_coverage.txt")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 4096),
        time_min = (lambda wildcards, attempt: attempt * 60)
    conda:
        CONDA_ENV_SAMTOOLS
    shell:
        """
        samtools coverage --output {output.all_chr_coverage} {input.bam_file}
        
        #Get only chromosomes of interest
        while IFS= read -r line; do chr_input=$(echo $line | awk '{{print $1}}'); for chr in {CHR_NUMBER}; do if [ $chr = $chr_input ]; then echo $line >> {output.coverage}; fi; done; done < {output.all_chr_coverage}

        #Sort chr numerically, with alphabetical chr names (X, Y...) at the end
        sort -k1,1V {output.coverage}

        #Insert header
        sed -i '1i #rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq' {output.coverage}

        #replace " " separator introduced by previous operations by \t        
        sed -i 's, ,\t,g' {output.coverage}
        """

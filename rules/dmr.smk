"""
##########################################################################
These rules launch Differential Methylation Regions (DMR) analysis
##########################################################################
"""

"""
This rule computes chromosome size file from fasta reference file
"""
rule get_chromosome_size:
    input:
        fai_file = config["references"]["genome_fai"]
    output:
        chromSize = temp(os.path.normpath(OUTPUT_DIR + "/tmp/resources/chromSize.txt")),
    threads:1
    resources:
        mem_mb=lambda wildcards, attempt: min(1024 + 2040 * (attempt - 1 ),8192),
        time_min = (lambda wildcards, attempt: attempt * 60)
    conda:
        CONDA_ENV_SAMTOOLS
    shell:
        """
        cat {input.fai_file} | cut -f1,2 > {output.chromSize}
        """


"""
This rule separate strands
"""
rule split_bed_by_strand:
    input:
        os.path.normpath(OUTPUT_DIR + "/tmp/bed_uncombined_strands/{samples_name}/{meth_type}/{samples_name}_chr{chr_number}_{meth_type}_uncomb.bed")
    output:
        fwd = temp(os.path.normpath(OUTPUT_DIR + "/Methylation_Analysis/bed_files/chr{chr_number}/{meth_type}/{samples_name}-chr{chr_number}-{meth_type}-fwd.bed")),
        rev = temp(os.path.normpath(OUTPUT_DIR + "/Methylation_Analysis/bed_files/chr{chr_number}/{meth_type}/{samples_name}-chr{chr_number}-{meth_type}-rev.bed"))
    threads:1
    resources:
        mem_mb=lambda wildcards, attempt: min(1024 + 2040 * (attempt - 1 ),8192),
        time_min = (lambda wildcards, attempt: attempt * 60)
    shell:
        """
        sed "s/ /\t/g" {input} |awk '$6 == "+"'|awk  '{{print $1,$2,$3,$11,$12,$10-$12}}' > {output.fwd}
        sed "s/ /\t/g" {input} |awk '$6 == "-"'|awk  '{{print $1,$2,$3,$11,$12,$10-$12}}' > {output.rev}
        """


"""
These rules separate references by chromosome and by strand
"""

rule split_alu_by_chr_strand:
    input:
        os.path.normpath(config["Alu"])
    output:
        fwd = temp(os.path.normpath(OUTPUT_DIR + "/tmp/resources/Alu/refseq.bed_Alu_chr{chr_number}_fwd.txt")),
        rev = temp(os.path.normpath(OUTPUT_DIR + "/tmp/resources/Alu/refseq.bed_Alu_chr{chr_number}_rev.txt"))
    threads:1
    resources:
        mem_mb=lambda wildcards, attempt: min(1024 + 2040 * (attempt - 1 ),8192),
        time_min = (lambda wildcards, attempt: attempt * 60)
    shell:
        """
        echo {wildcards.chr_number}
        chromos_to_search="{wildcards.chr_number}"
        cat {input} | awk 'BEGIN{{FS=OFS="\\t"}} {{ gsub("chrM", "chrMT", $1); print }}' | tr -d "^chr" | awk -v chromos_to_search="$chromos_to_search" '$1 == chromos_to_search' | awk '$4 == "+"'|cut -f1,2,3,5 > {output.fwd}
        cat {input} | awk 'BEGIN{{FS=OFS="\\t"}} {{ gsub("chrM", "chrMT", $1); print }}' | tr -d "^chr" | awk -v chromos_to_search="$chromos_to_search" '$1 == chromos_to_search' | awk '$4 == "-"'|cut -f1,2,3,5 > {output.rev}
        """

rule split_transcript_by_chr_strand:
    input:
        os.path.normpath(config["Transcript"])
    output:
        fwd = temp(os.path.normpath(OUTPUT_DIR + "/tmp/resources/Transcript/refseq.bed_Transcript_chr{chr_number}_fwd.txt")),
        rev = temp(os.path.normpath(OUTPUT_DIR + "/tmp/resources/Transcript/refseq.bed_Transcript_chr{chr_number}_rev.txt"))
    threads:1
    resources:
        mem_mb=lambda wildcards, attempt: min(1024 + 2040 * (attempt - 1 ),8192),
        time_min = (lambda wildcards, attempt: attempt * 60)
    shell:
        """
        echo {wildcards.chr_number}
        chromos_to_search="{wildcards.chr_number}"
        cat {input} | awk 'BEGIN{{FS=OFS="\\t"}} {{ gsub("chrM", "chrMT", $1); print }}' | tr -d "^chr" | awk -v chromos_to_search="$chromos_to_search" '$1 == chromos_to_search' | awk '$6 == "+"' > {output.fwd}
        cat {input} | awk 'BEGIN{{FS=OFS="\\t"}} {{ gsub("chrM", "chrMT", $1); print }}' | tr -d "^chr" | awk -v chromos_to_search="$chromos_to_search" '$1 == chromos_to_search' | awk '$6 == "-"' > {output.rev}
        """


rule split_cpg_by_chr:
    input:
        os.path.normpath(config["CpG"])
    output:
        temp(os.path.normpath(OUTPUT_DIR + "/tmp/resources/CpG/cpgi.bed_CpG_chr{chr_number}.txt"))
    threads:1
    resources:
        mem_mb=lambda wildcards, attempt: min(1024 + 2040 * (attempt - 1 ),8192),
        time_min = (lambda wildcards, attempt: attempt * 60)
    shell:
        """
        echo {wildcards.chr_number}
        chromos_to_search="{wildcards.chr_number}"
        cat {input} | awk 'BEGIN{{FS=OFS="\\t"}} {{ gsub("chrM", "chrMT", $1); print }}' | tr -d "^chr" | awk -v chromos_to_search="$chromos_to_search" '$1 == chromos_to_search' > {output}
   
        """


"""
These rules produce the mean table
"""

###MODIFIER LES CODE GENEDMR POUR RECEVOIR LES 2 FICHIERS

if config["steps"]["differential_methylation_sample"]:
    rule DMR_mean_table_combination:
        input:
            control = lambda wildcards : os.path.normpath(OUTPUT_DIR + "/Methylation_Analysis/bed_files/chr{chr_number}/{meth_type}/" + wildcards.pair_methyl.split("_vs_")[1] + "-chr{chr_number}-{meth_type}-{strand}.bed"),
            case = lambda wildcards : os.path.normpath(OUTPUT_DIR + "/Methylation_Analysis/bed_files/chr{chr_number}/{meth_type}/" + wildcards.pair_methyl.split("_vs_")[0] + "-chr{chr_number}-{meth_type}-{strand}.bed"),
            alu = os.path.normpath(OUTPUT_DIR + "/tmp/resources/Alu/refseq.bed_Alu_chr{chr_number}_{strand}.txt"),
            trans = os.path.normpath(OUTPUT_DIR + "/tmp/resources/Transcript/refseq.bed_Transcript_chr{chr_number}_{strand}.txt"),
            cpg = os.path.normpath(OUTPUT_DIR + "/tmp/resources/CpG/cpgi.bed_CpG_chr{chr_number}.txt")
        output:
            os.path.normpath(OUTPUT_DIR + "/Methylation_Analysis/{ref_type}/{meth_type}/{pair_methyl}/mean_table_{pair_methyl}_chr{chr_number}_{strand}.tsv")
        params:
            ref_folder = os.path.normpath(OUTPUT_DIR + "/tmp/resources/{ref_type}/"),
            script = os.path.normpath(PIPELINE_DIR + "/script")
        threads:1
        resources:
            mem_mb=lambda wildcards, attempt: min(10240 + 2040 * (attempt - 1 ),81920),
            time_min = (lambda wildcards, attempt: attempt * 300)
        shell:
            """
            echo {params.ref_folder}
            singularity exec --contain -B {params.script},{params.ref_folder},{OUTPUT_DIR} {SING_ENV_GENEDMR} \
            Rscript {params.script}/GeneDMRs_code{wildcards.ref_type}_final.R {input.control} {input.case} {wildcards.chr_number} {wildcards.meth_type} {wildcards.strand} {params.ref_folder} {output}
            """

if config["steps"]["differential_methylation_condition"]:
    rule DMR_mean_table_condition:
        input:
            control = expand(os.path.normpath(OUTPUT_DIR + "/Methylation_Analysis/bed_files/chr{{chr_number}}/{{meth_type}}/{controls}-chr{{chr_number}}-{{meth_type}}-{{strand}}.bed"),controls=CONTROL),
            case = expand(os.path.normpath(OUTPUT_DIR + "/Methylation_Analysis/bed_files/chr{{chr_number}}/{{meth_type}}/{cases}-chr{{chr_number}}-{{meth_type}}-{{strand}}.bed"),cases=CASE),
            alu = os.path.normpath(OUTPUT_DIR + "/tmp/resources/Alu/refseq.bed_Alu_chr{chr_number}_{strand}.txt"),
            trans = os.path.normpath(OUTPUT_DIR + "/tmp/resources/Transcript/refseq.bed_Transcript_chr{chr_number}_{strand}.txt"),
            cpg = os.path.normpath(OUTPUT_DIR + "/tmp/resources/CpG/cpgi.bed_CpG_chr{chr_number}.txt")
        output:
            os.path.normpath(OUTPUT_DIR + "/Methylation_Analysis/{ref_type}/{meth_type}/cases_vs_controls/mean_table_cases_vs_controls_chr{chr_number}_{strand}.tsv")
        params:
            control = ",".join(expand(os.path.normpath(OUTPUT_DIR + "/Methylation_Analysis/bed_files/chr{{chr_number}}/{{meth_type}}/{controls}-chr{{chr_number}}-{{meth_type}}-{{strand}}.bed"),controls=CONTROL)),
            case = ",".join(expand(os.path.normpath(OUTPUT_DIR + "/Methylation_Analysis/bed_files/chr{{chr_number}}/{{meth_type}}/{cases}-chr{{chr_number}}-{{meth_type}}-{{strand}}.bed"),cases=CASE)),
            ref_folder = os.path.normpath(OUTPUT_DIR + "/tmp/resources/{ref_type}/"),
            script = os.path.normpath(PIPELINE_DIR + "/script")
        threads:1
        resources:
            mem_mb=lambda wildcards, attempt: min(10240 + 2040 * (attempt - 1 ),81920),
            time_min = (lambda wildcards, attempt: attempt * 300)
        shell:
            """
            echo {params.ref_folder}
            singularity exec --contain -B {params.script},{params.ref_folder},{OUTPUT_DIR} {SING_ENV_GENEDMR} \
            Rscript {params.script}/GeneDMRs_code{wildcards.ref_type}_final.R {params.control} {params.case} {wildcards.chr_number} {wildcards.meth_type} {wildcards.strand} {params.ref_folder} {output}
            """


rule DMR_stat_annot_graph_combination:
    input: #concat chr
        tsv = expand(os.path.normpath(OUTPUT_DIR + "/Methylation_Analysis/{{ref_type}}/{{meth_type}}/{{pair_methyl}}/mean_table_{{pair_methyl}}_chr{chr}_{strand}.tsv"), chr = CHR_NUMBER, strand = STRAND),
        ref_conversion = config["references"]["code_symbol_conversion"],
        chromSize = os.path.normpath(OUTPUT_DIR + "/tmp/resources/chromSize.txt"),
    output:
        concat = os.path.normpath(OUTPUT_DIR + "/Methylation_Analysis/{ref_type}/{meth_type}/{pair_methyl}/mean_table_{pair_methyl}_all_chr.tsv"),
        stat = os.path.normpath(OUTPUT_DIR + "/Methylation_Analysis/{ref_type}/{meth_type}/{pair_methyl}/meth_res_table_{pair_methyl}_all_chr.csv"),
        #filtered = os.path.normpath(OUTPUT_DIR + "/Methylation_Analysis/{ref_type}/{meth_type}/{pair_methyl}/filtered_table_{pair_methyl}_all_chr.tsv"),
        #annot = os.path.normpath(OUTPUT_DIR + "/Methylation_Analysis/{ref_type}/{meth_type}/{pair_methyl}/annot_table_{pair_methyl}_all_chr.tsv"),
    params:
        ref_folder = os.path.normpath(OUTPUT_DIR + "/tmp/resources/{ref_type}/"),
        script = os.path.normpath(PIPELINE_DIR + "/script"),
        ref_conversion_folder = os.path.dirname(config["references"]["code_symbol_conversion"])
    threads:1
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 10240),
        time_min = (lambda wildcards, attempt: attempt * 300)
    shell:
        """
        cat {input.tsv} > {output.concat}
        singularity exec --contain -B {params.script},{params.ref_folder},{OUTPUT_DIR},{params.ref_conversion_folder} {SING_ENV_GENEDMR} \
        Rscript {params.script}/GeneDMRs_code_all_chr_final.R {output.concat} {input.ref_conversion} {input.chromSize}
        """
        
rule DMR_stat_annot_graph_condition:
    input:
        tsv = expand(os.path.normpath(OUTPUT_DIR + "/Methylation_Analysis/{{ref_type}}/{{meth_type}}/cases_vs_controls/mean_table_cases_vs_controls_chr{chr_number}_{strand}.tsv"), chr_number = CHR_NUMBER, strand = STRAND),
        ref_conversion = config["references"]["code_symbol_conversion"],
        chromSize = os.path.normpath(OUTPUT_DIR + "/tmp/resources/chromSize.txt"),
    output:
        concat = os.path.normpath(OUTPUT_DIR + "/Methylation_Analysis/{ref_type}/{meth_type}/cases_vs_controls/mean_table_cases_vs_controls_all_chr.tsv"),
        stat = os.path.normpath(OUTPUT_DIR + "/Methylation_Analysis/{ref_type}/{meth_type}/cases_vs_controls/meth_res_table_cases_vs_controls_all_chr.csv"),
        #filtered = os.path.normpath(OUTPUT_DIR + "/Methylation_Analysis/{ref_type}/{meth_type}/cases_vs_controls/filtered_table_cases_vs_controls_all_chr.tsv"),
        #annot = os.path.normpath(OUTPUT_DIR + "/Methylation_Analysis/{ref_type}/{meth_type}/cases_vs_controls/annot_table_cases_vs_controls_all_chr.tsv"),
    threads:1
    params:
        ref_folder = os.path.normpath(OUTPUT_DIR + "/tmp/resources/{ref_type}/"),
        script = os.path.normpath(PIPELINE_DIR + "/script"),
        ref_conversion_folder = os.path.dirname(config["references"]["code_symbol_conversion"])
    resources:
        mem_mb = (lambda wildcards, attempt: attempt * 10240),
        time_min = (lambda wildcards, attempt: attempt * 300)
    shell:
        """
        cat {input.tsv} > {output.concat}
        singularity exec --contain -B {params.script},{params.ref_folder},{OUTPUT_DIR},{params.ref_conversion_folder} {SING_ENV_GENEDMR} \
        Rscript {params.script}/GeneDMRs_code_all_chr_final.R {output.concat} {input.ref_conversion} {input.chromSize}
        """


    

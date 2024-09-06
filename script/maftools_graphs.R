# 0- Set environment
library(optparse)
library(stringr)
library(maftools)

# 1- Get script parameters
option_list <- list(
  make_option("--sampleName", help="Sample name"),
  make_option("--inputMAF", help="Input MAF file absolute path"),
  make_option("--inputMAFsList", help="MAF files absolute paths, separated by a regular space. These files will be merged to get single plots for multiple samples comparison."),
  make_option("--outputDir", help="Output directory absolute path"),
  make_option("--variant_type", help="somatic or germline"),
  make_option("--genesFile", help="CSV file of genes of interest (one gene per line)."))

parser <- OptionParser(usage="Rscript %prog [options]", description = " ", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)

sample_name <- args$options$sampleName
maf_file <- args$options$inputMAF
maf_files_list <- args$options$inputMAFsList
outdir <- args$options$outputDir
genesFile <- args$options$genesFile
variant_type <- args$options$variant_type

# Check parameters
if(is.null(sample_name)) stop("--sample_name parameter must be set!")
if(is.null(maf_file)) stop("--maf_file parameter must be set!")
if(is.null(maf_file) && is.null(maf_files_list)) stop("--inputMAF or --inputMAFsList parameter must be set!")
if(is.null(outdir)) stop("--outputDir parameter must be set!")
if(is.null(variant_type)) stop("--variant_type parameter must be set!")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

print("Entered parameters are:")
print(paste("sampleName: ",sample_name))
if (!is.null(maf_file)) print(paste("inputMAF: ",maf_file))
if (!is.null(maf_files_list)) print(paste("inputMAFsList: ",maf_files_list))
print(paste("outputDir: ",outdir))
if (!is.null(genesFile)) print(paste("genesFile: ",genesFile))
print(paste("variant_type: ",variant_type))

# Create flag file
flag_file = paste0(outdir,sample_name,"_maftools_graphs_DONE.txt")
file.create(flag_file)

# 2- Load data from MAF file
if (!is.null(maf_file)) {
    if(file.exists(maf_file)) {
        print("Loading single MAF file...")
        tryCatch({
            maf_data <- read.maf(maf_file)
        }, error = function(err) {
            txt = paste("No maftools graph could be obtained as input maf file could not be read: ", err)
            write(txt, flag_file, append=TRUE)
            txt_check = paste("--> Please check that your input file was effectively annotated by SnpEff prior to using maftools and if you have non synonymous variants to look at.")
            write(txt_check, flag_file, append=TRUE)
        })
    } else {
        stop(paste(" Your maf file doesn't exist. Check your parameter: ", maf_file))
    }
} else if(!is.null(maf_files_list)){
    print("Loading list of MAF files...")
    tryCatch({
        maf_data <- merge_mafs(lapply(Sys.glob(maf_files_list), read.maf))
    }, error = function(err){
        stop(paste("There was a problem while charging the MAF file(s) which provided paths are: maf_files_list = ", maf_files_list))
    })
}

# 3- Generate plots as PDF files
if (exists("maf_data")) {
    print("Generating plots...")
    pdf(file = paste0(outdir,sample_name,"_oncoplot.pdf"))
    oncoplot(maf_data)
    dev.off()
    
    pdf(file = paste0(outdir,sample_name,"_mafSummary_plot.pdf"))
    plotmafSummary(maf = maf_data, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
    dev.off()
    
    pdf(file = paste0(outdir,sample_name,"_mafbarplot.pdf"))
    mafbarplot(maf = maf_data)
    dev.off()
    
    maf_data.titv = titv(maf = maf_data, plot = FALSE, useSyn = TRUE)
    pdf(file = paste0(outdir,sample_name,"_titv_Transitions_Transversions.pdf"))
    plotTiTv(res = maf_data.titv)
    dev.off()
    
    # 4- Generate plots for specific genes if a list of genes of interest is provided
    if (is.null(genesFile)) {
        txt="No gene_name was provided, no lollipop plot will be generated."
        print(txt)
        write(txt, flag_file, append=TRUE)
    } else {
        #Set variant_type option
        if(variant_type =="somatic") variant_type <- TRUE
        if(variant_type =="germline") variant_type <- FALSE
        
        #Read data
        genes_list <- read.table(genesFile, sep=",", header = FALSE)$V1
    
        #Read gff file from maftools R package to find associated transcripts to genes
        gff = system.file('extdata', 'protein_domains.RDs', package = 'maftools')
        gff = readRDS(file = gff)
      
        #Loop on genes list
        for (gene_name in genes_list) {
    
            #Find all transcripts for this gene
            prot = gff[HGNC %in% gene_name]
            if(dim(prot)[2] != 0){
                transcript_list <- unique(sort(prot$refseq.ID))
    
                #Loop on transcript list
                for (transcript_name in transcript_list) {
                    print(paste("Generating lollipop plots for the transcript", transcript_name, "of the gene", gene_name, "."))
                    tryCatch({
                        pdf(file = paste0(outdir,sample_name,"_",gene_name,"_",transcript_name,"_lollipopPlot_mutations.pdf"))
                        lollipopPlot(
                          maf = maf_data,
                          gene = gene_name,
                          showMutationRate = variant_type,
                          labelPos = 'all', 
                          labPosSize = 0.8,
                          refSeqID = transcript_name
                        )
                        dev.off()
                    }, error = function(err) {
                        print(err)
                        dev.off()
                        file.remove(paste0(outdir,sample_name,"_",gene_name,"_",transcript_name,"_lollipopPlot_mutations.pdf"))
                    })
                }
            }else print(paste0("Error: No transcript found for the gene ", gene_name))
        }
    }
}

#Print that the script finished in the flag file
txt = ("Finished maftools_graphs.R execution!")
write(txt, flag_file, append=TRUE)
print("Finished !")

#subsetMaf(maf = data, genes = 'TP53', mafObj = FALSE, fields ="Transcript_ID")
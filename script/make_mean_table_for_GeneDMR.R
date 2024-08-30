#script to make methylation table for differential methylation analysis depending on the regions of interest
# this script will replace GeneDMRs_codeAlu_final.R, GeneDMRs_codeCpG_final.R and GeneDMRs_codeTranscript_final.R scripts.

### Library
library(optparse)
#library(GeneDMRs)
#library(dplyr)

### Parameters
## parameters list
option_list <- list(
  make_option("--inputCase", help="Input bed file for condition."),
  make_option("--inputCtrl", help="Input bed file for control."),
  make_option("--inputRef", help="Input reference file(s). For Cpg analyse set the Cpg reference first then the Transcript reference in second, separated by comma"),
  make_option("--output_file", help="Output file"))
parser <- OptionParser(usage="Rscript %prog [options]", description = " ", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)

## get parameters
#Alu test
#inputCase <- "/mnt/beegfs/scratch/n_rabearivelo/test_pipeline/methylation/data_output/methylation_analysis/bed_files/chr22/5mCG/1033_CD14-chr22-5mCG-fwd.bed"
#inputCtrl <- "/mnt/beegfs/scratch/n_rabearivelo/test_pipeline/methylation/data_output/methylation_analysis/bed_files/chr22/5mCG/1033_CD19-chr22-5mCG-fwd.bed"
#inputRef_all <- "/mnt/beegfs/scratch/n_rabearivelo/test_pipeline/methylation/data_output/resources/Alu/refseq.bed_Alu_chr22_fwd.txt"
#output_file <- "/mnt/beegfs/scratch/n_rabearivelo/test_pipeline/methylation/data_output/methylation_analysis/Alu/5mCG/1033_CD14-1033_CD19/mean_table_1033_CD14-1033_CD19_chr22_fwd.tsv"
#CpG test
#inputCase <- "/mnt/beegfs/scratch/n_rabearivelo/test_pipeline/methylation/data_output/methylation_analysis/bed_files/chr22/5mCG/1033_CD14-chr22-5mCG-fwd.bed"
#inputCtrl <- "/mnt/beegfs/scratch/n_rabearivelo/test_pipeline/methylation/data_output/methylation_analysis/bed_files/chr22/5mCG/1033_CD19-chr22-5mCG-fwd.bed"
#inputRef_all <- "/mnt/beegfs/scratch/n_rabearivelo/test_pipeline/methylation/data_output/resources/CpG/cpgi.bed_CpG_chr22.txt,/mnt/beegfs/scratch/n_rabearivelo/test_pipeline/methylation/data_output/resources/Transcript/refseq.bed_Transcript_chr22_fwd.txt"
#output_file <- "/mnt/beegfs/scratch/n_rabearivelo/test_pipeline/methylation/data_output/methylation_analysis/CpG/5mCG/1033_CD14-1033_CD19/mean_table_1033_CD14-1033_CD19_chr22_fwd.tsv"
#Transcript test
#inputCase <- "/mnt/beegfs/scratch/n_rabearivelo/test_pipeline/methylation/data_output/methylation_analysis/bed_files/chr22/5mCG/1033_CD14-chr22-5mCG-fwd.bed"
#inputCtrl <- "/mnt/beegfs/scratch/n_rabearivelo/test_pipeline/methylation/data_output/methylation_analysis/bed_files/chr22/5mCG/1033_CD19-chr22-5mCG-fwd.bed"
#inputRef_all <- "/mnt/beegfs/scratch/n_rabearivelo/test_pipeline/methylation/data_output/resources/Transcript/refseq.bed_Transcript_chr22_fwd.txt"
#output_file <- "/mnt/beegfs/scratch/n_rabearivelo/test_pipeline/methylation/data_output/methylation_analysis/Transcript/5mCG/1033_CD14-1033_CD19/mean_table_1033_CD14-1033_CD19_chr22_fwd.tsv"


inputCase <- args$options$inputCase
inputCtrl <- args$options$inputCtrl
inputRef_all <- args$options$inputRef
output_file <- args$options$output_file

## check parameters
if(is.null(inputCase)) stop("--inputCase parameter must be set!")
if(is.null(inputCtrl)) stop("--inputCtrl parameter must be set!")
if(is.null(inputRef_all)) stop("--inputRef parameter must be set!")
if(is.null(output_file)) stop("--output_file parameter must be set!")

## get missing parameters
inputRef <- unlist(strsplit(inputRef_all,","))[1]
chr <- sub("chr","",unlist(strsplit(basename(inputRef),"_"))[3])
ref_path <- dirname(inputRef)
feature <- if(basename(ref_path) == "Alu") FALSE else TRUE
bedfile <- unlist(strsplit(basename(inputRef),"\\."))[1]
suffix <- sub(".*bed","",basename(inputRef))

### Function
read_bed_reference_file <- function(paths = paste(system.file(package = "GeneDMRs"), "/methdata", sep=""), bedfile = "refseq",
                         suffix = ".txt", feature = FALSE, featurewrite = TRUE){
    # set the paths #
    setwd(paths)
    # get regions #
    if(feature == FALSE){
        regionfile <- GeneDMRs::Regionfile_read(bedfile = bedfile, suffix = suffix)
        return(regionfile)
    }else{
        # read refseq or cpgi file #
        if(bedfile == "refseq"){
            # region file is the reference gene bed file#
            regionfile <- GeneDMRs::Regionfile_read(bedfile = bedfile, suffix = suffix)

            geneobj <- tryCatch({
                            geneobj <- genomation::readTranscriptFeatures(paste(bedfile, ".bed", suffix, sep = ""))
                        }, error = function(err) {
                            print(err)
                            print("There is probably no introns in the reference (like for MT chromosome), try without introns...")
                            
                            #inspired from https://github.com/BIMSBbioinfo/genomation/blob/master/R/readData.R "readTranscriptFeatures" function
                            library(GenomicRanges)
                            down.flank=1000
                            up.flank=1000
                            # readBed6
                            message('Reading the table...\r')
                            bed=genomation::readTableFast(paste(bedfile, ".bed", suffix, sep = ""),header=FALSE,skip="auto")                    
                            # exons
                            message('Calculating exon coordinates...\r')
                            exons    = genomation::convertBed2Exons(bed)
                            # get the locations of TSSes
                            message('Calculating TSS coordinates...\r')
                            tss=bed
                            #  + strand
                            tss[tss$V6=="+",3] = tss[tss$V6=="+",2]
                            #  - strand
                            tss[tss$V6=="-",2]=tss[tss$V6=="-",3]
                            # make Granges object
                            tssg = GRanges(seqnames=as.character(tss$V1),
                                           ranges=IRanges(start=tss$V2, end=tss$V3),
                                           strand=as.character(tss$V6),
                                           score=rep(0,nrow(tss)),
                                           name=tss$V4)
                            message('Calculating promoter coordinates...\r')
                            # get the locations of promoters
                            # + strand
                            bed[bed$V6=="+",3]=bed[bed$V6=="+",2]+down.flank
                            bed[bed$V6=="+",2]=bed[bed$V6=="+",2]-up.flank
                            #  - strand
                            bed[bed$V6=="-",2]=bed[bed$V6=="-",3]-down.flank
                            bed[bed$V6=="-",3]=bed[bed$V6=="-",3]+up.flank
                            # make the object
                            prom.df = unique(bed[,c(1,2,3,6)])
                            prom = GRanges(seqnames=as.character(prom.df$V1),
                                             ranges=IRanges(start=prom.df$V2, end=prom.df$V3),
                                             strand=as.character(prom.df$V6),
                                             score=rep(0,nrow(prom.df)),
                                             name=rep(".",nrow(prom.df)) )
                            message('Outputting the final GRangesList...\r\n')
                            geneobj <- GRangesList(exons=exons,promoters=prom,TSSes=tssg)
                            return(geneobj)
                        })
          
            # Output the defined gene features based on the .bed file #
            if(featurewrite) utils::write.table(geneobj, paste0(bedfile, ".bed", suffix,"_Genebody.txt"), col.names = F, row.names = F, quote = F)
            # check the geneobj file for promoter regions #
            geneobj <- GeneDMRs::refseqfile_check(base::as.data.frame(geneobj), regionfile)
            return(geneobj)
        }else if(bedfile == "cpgi"){
            # region file is the cpg bed file (without header) #
            cpgiobj <- genomation::readFeatureFlank(paste0(bedfile, ".bed", suffix), feature.flank.name = c("CpGisland", "Shores"))
            # Output the defined gene features based on the .bed file #
            if(featurewrite) utils::write.table(cpgiobj, paste0(bedfile, ".bed", suffix,"_Cpgifeature.txt"), col.names = F, row.names = F, quote = F)
            cpgiobj <- base::as.data.frame(cpgiobj)[,c(3,4,5,2)]
            colnames(cpgiobj) <- c("chr", "start", "end", "cpgfeature")
            # rename the "cpgi" #
            cpginum <- sum(cpgiobj$cpgfeature == "CpGisland")
            cpgiobj$cpgi <- c(paste0("cpgi", 1:cpginum), paste0("shore", 1:(nrow(cpgiobj) - cpginum)))
            # sort the file #
            cpgiobj <- dplyr::arrange(cpgiobj[,c(5,1:4)], chr, start)
            return(cpgiobj)
        }else stop("Wrong bed file names and please use 'refseq' or 'cpgi' name")
    }
}

### Main
cat(paste0("The control(s) file(s):\n", inputCtrl, "\n"))
cat(paste0("The case(s) file(s):\n", inputCase, "\n"))
cat(paste0("The reference file(s):\n", inputRef_all, "\n"))
cat(paste0("The output file:\n", output_file, "\n"))

## check if the reference is empty (for example, there is no Alu or CpG for MT chromosome in human)
if(file.size(inputRef) == 0L){
    print(paste0("Warning: reference file is empty: ", inputRef))
    print("Results will be empty too!")
    regiongeneall <- data.frame()
    
}else{

    ## open reference(s)
    print("Loading reference file")
    inputrefseqfile <- read_bed_reference_file(paths = ref_path, bedfile = bedfile, suffix = suffix, feature = feature, featurewrite = TRUE)
    if(basename(ref_path) == "CpG"){
        cpgifeaturefile <- inputrefseqfile
        inputRef <- unlist(strsplit(inputRef_all,","))[2]
        ref_path <- dirname(inputRef)
        bedfile <- unlist(strsplit(basename(inputRef),"\\."))[1]
        suffix <- sub(".*bed","",basename(inputRef))
        inputrefseqfile <- read_bed_reference_file(paths = ref_path, bedfile = bedfile, suffix = suffix, feature = feature, featurewrite = FALSE)
    } else cpgifeaturefile <- NULL
    
    ## create table containing all the sample
    print("Read the methylation files")
    inputmethfile <- GeneDMRs::Methfile_read(control_paths = inputCtrl, case_paths = inputCase)
    
    ## quality control
    print("Compute quality control")
    inputmethfile_QC <- GeneDMRs::Methfile_QC(inputmethfile)
    
    ## methylation mean
    print("Calculate the methylation mean for regions")
    date()
    regiongeneall <- GeneDMRs::Methmean_region(inputmethfile_QC = inputmethfile_QC, inputrefseqfile = inputrefseqfile, cpgifeaturefile = cpgifeaturefile, chrnum = chr, featurename = feature)
    date()
}

## save result
print("Save result")
utils::write.table(regiongeneall,output_file,sep="\t",quote=F,row.names=F)


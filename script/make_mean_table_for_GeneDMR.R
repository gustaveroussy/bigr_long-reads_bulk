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

#replace as.interger by as.numeric; needed for big interger, else value is set to NA.
bed12ToIntrons2 <- function(ref){

    #remove the genes with one exon only (they won't have any introns)
    ref=ref[ref[,10]>1,]
    ids=paste(ref[,1],ref[,2],ref[,3],ref[,4],sep="")
    ref=cbind(ref,id=ids)
    ref=unique(ref)

    # apply strsplit on columns 11 and 12 to extract exon start positions and exon sizes
    b.start.size=cbind(as.numeric(unlist(strsplit(as.character(ref$V12),"," ))),
                     as.numeric(unlist(strsplit(as.character(ref$V11),"," ))))

    # replicate rows occurs as many as its exon number
    rep.ref = ref[rep(1:nrow(ref),ref[,10]),] 

    exon.id = unlist( mapply( function(x,y) 
                                if(x=="+"){return(1:y)}
                                else{return(y:1)} ,
                                ref[,6],ref[,10]  ) )
    rep.ref$V5 = exon.id

    # calculate exon start and ends
    rep.ref$V3 = rep.ref$V2+b.start.size[,1]+b.start.size[,2] 
    rep.ref$V2 = rep.ref$V2+b.start.size[,1]
    rep.ref = rep.ref[,c(1:6,13)]

    # now try to get the exons by cbinding consecutive exons
    temp.ref        = cbind(rep.ref[1:(nrow(rep.ref)-1),],rep.ref[2:nrow(rep.ref),])
    temp.ref        = temp.ref[temp.ref[,7]==temp.ref[,14],]
    temp.ref[,2]    = temp.ref[,3]
    temp.ref[,3]    = temp.ref[,9]
    rep.ref            = temp.ref[,1:6]

    # subtract 1 from - strand exon numbers
    rep.ref[rep.ref[,6]=="-",5]=rep.ref[rep.ref[,6]=="-",5]-1 

    strands=as.character(rep.ref$V6)
    strands[strands=="."]="*"
    
    #width <- 1.0 + end - start
    #start <- rep.ref$V2+1
    #end <- rep.ref$V3
    
    df.introns <- data.frame(number_in_GRanges = rep(2,nrow(rep.ref)),
                             type = rep("introns",nrow(rep.ref)),
                             seqnames = as.character(rep.ref$V1),
                             start = rep.ref$V2+1,
                             end = rep.ref$V3,
                             width = rep.ref$V3 - rep.ref$V2,
                             strand = strands,
                             score = rep.ref$V5,
                             name = rep.ref$V4
    )

    return(df.introns)
}

#replace as.interger by as.numeric; needed for big interger, else value is set to NA.
bed12ToExons2 <- function(ref){

    ref=unique(ref)

    # apply strsplit on columns 11 and 12 to extract exon start positions and exon sizes
    b.start.size = cbind(as.numeric(unlist(strsplit(as.character(ref$V12),"," ))),
                       as.numeric(unlist(strsplit(as.character(ref$V11),"," ))))
    rep.ref = ref[rep(1:nrow(ref),ref[,10]),] # replicate rows occurs as many as its exon number

    exon.id = unlist( mapply( function(x,y)
                              if(x=="+"){return(1:y)}
                              else{return(y:1)}, 
                              ref[,6],ref[,10]))
    rep.ref$V5=exon.id

    rep.ref$V3 = rep.ref$V2+b.start.size[,1]+b.start.size[,2] # calculate exon start and ends
    rep.ref$V2 = rep.ref$V2+b.start.size[,1]

    strands = as.character(rep.ref$V6)
    strands[strands=="."]="*"
    df.exons <- data.frame(number_in_GRanges = rep(1,nrow(rep.ref)),
                             type = rep("exons",nrow(rep.ref)),
                             seqnames = as.character(rep.ref$V1),
                             start = rep.ref$V2+1,
                             end = rep.ref$V3,
                             width = rep.ref$V3 - rep.ref$V2,
                             strand = strands,
                             score = rep.ref$V5,
                             name = rep.ref$V4
    )

    return(df.exons)
}

convertBed2Introns2 <- function(bed.df = NULL){
    if(! genomation::checkBedValidity(bed.df,"intron")) stop("this is not a valid bed file")
    return(bed12ToIntrons2(bed.df))
}

convertBed2Exons2 <- function(bed.df = NULL){
    if(! genomation::checkBedValidity(bed.df,"exon")) stop("this is not a valid bed file")
    return(bed12ToExons2(bed.df))
}

readTranscriptFeatures2 <- function(location = NULL){
    #inspired from https://github.com/BIMSBbioinfo/genomation/blob/master/R/readData.R "readTranscriptFeatures" function
    library(GenomicRanges)
    down.flank=1000
    up.flank=1000
    
    # readBed6
    message('Reading the table...\r')
    bed <- genomation::readTableFast(location ,header = FALSE, skip = "auto")                    
    
    # introns
    message('Calculating intron coordinates...\r')
    introns <- tryCatch({
        introns <- convertBed2Introns2(bed)
    }, error = function(err) {
            print(err)
            print("There is probably no introns in the reference (like for MT chromosome), try without introns...")
            introns = data.frame()
            return(introns)
    })
    # exons
    message('Calculating exon coordinates...\r')
    exons <- tryCatch({
        exons <- convertBed2Exons2(bed)
    }, error = function(err) {
            print(err)
            print("There is probably no exons in the reference, try without exons...")
            exons = data.frame()
            return(exons)
    })
    
    # get the locations of TSSes
    message('Calculating TSS coordinates...\r')
    tss <- bed
    #  + strand
    tss[tss$V6=="+",3] <- tss[tss$V6=="+",2]
    #  - strand
    tss[tss$V6=="-",2] <- tss[tss$V6=="-",3]

    tssg <- data.frame(number_in_GRanges = rep(4,nrow(tss)),
                             type = rep("TSSes",nrow(tss)),
                             seqnames = as.character(tss$V1),
                             start = tss$V2,
                             end = tss$V3,
                             width = 1.0 + tss$V3 - tss$V2,
                             strand = tss$V6,
                             score = rep(0,nrow(tss)),
                             name = tss$V4
    )
    
    # get the locations of promoters
    message('Calculating promoter coordinates...\r')
    # + strand
    bed[bed$V6=="+",3]=bed[bed$V6=="+",2]+down.flank
    bed[bed$V6=="+",2]=bed[bed$V6=="+",2]-up.flank
    #  - strand
    bed[bed$V6=="-",2]=bed[bed$V6=="-",3]-down.flank
    bed[bed$V6=="-",3]=bed[bed$V6=="-",3]+up.flank

    # set unique promoters
    prom.df = unique(bed[,c(1,2,3,6)])
    prom <- data.frame(number_in_GRanges = rep(3,nrow(prom.df)),
                             type = rep("promoters",nrow(prom.df)),
                             seqnames = as.character(prom.df$V1),
                             start = prom.df$V2,
                             end = prom.df$V3,
                             width = 1.0 + prom.df$V3 - prom.df$V2,
                             strand = as.character(prom.df$V6),
                             score = rep(0,nrow(prom.df)),
                             name = rep(".",nrow(prom.df))
    )
    
    message('Outputting the final dataframe...\r\n')
    final_df <- rbind(exons, introns, prom, tssg)
    
    return(final_df)
}

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
            geneobj <- readTranscriptFeatures2(paste(bedfile, ".bed", suffix, sep = ""))
            
            # Output the defined gene features based on the .bed file #
            if(featurewrite) utils::write.table(geneobj, paste0(bedfile, ".bed", suffix,"_Genebody.txt"), col.names = F, row.names = F, quote = F)
            # check the geneobj file for promoter regions #
            geneobj <- GeneDMRs::refseqfile_check(geneobj, regionfile)

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
    if (!is.null(inputrefseqfile)){
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
    } else regiongeneall <- data.frame()  #no exons and introns
}

## save result
print("Save result")
utils::write.table(regiongeneall,output_file,sep="\t",quote=F,row.names=F)


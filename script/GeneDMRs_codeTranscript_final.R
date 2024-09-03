#script for transcrit DMR


args = commandArgs(trailingOnly=TRUE)
library(data.table)
library(GeneDMRs)
library(dplyr)
library(GenomicRanges)

control=unlist(strsplit(args[1],","))
case=unlist(strsplit(args[2],","))


chr=args[3]
meth_type=args[4]
strand=args[5]
ref_path=args[6]
output_file=args[7]
#file_path<-paste0(analysis_path,"bed_files/chr_",chr,"/",meth_type,"/")
#filelist<-list.files(path=file_path,pattern=paste0("*",strand,".bed"),full.names=F)


print("The control(s) file(s):")
print(control)
print("The case(s) file(s):")
print(case)


print("for chromosome ")
print(chr)
print("and for strand:")
print(strand)
print("ref file:")
#print(paste0(ref_path,"refseq.bed_ucsc_Transcript_chr",chr,"_grch38_",strand,".txt"))
print(paste0(ref_path,"/refseq.bed_Transcript_chr",chr,"_",strand,".txt"))


### FUNCTIONS ###

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


Bedfile_read2 <- function(paths = paste(system.file(package = "GeneDMRs"), "/methdata", sep=""), bedfile = "refseq",
                         suffix = ".txt", feature = FALSE, featurewrite = TRUE){
  # set the paths #
  setwd(paths)
  # read refseq or cpgi file #
  if(bedfile=="refseq"){
    if(feature == FALSE){
      regionfile <- Regionfile_read(bedfile, suffix)
      return(regionfile)
    }else{
      # region file is the reference gene or cpg bed file (without header) #
      genebody_name_file <- paste0(bedfile, ".bed", suffix,"_Genebody.txt")
      regionfile <- Regionfile_read(bedfile, suffix)
      geneobj <- readTranscriptFeatures2(paste(bedfile, ".bed", suffix, sep = ""))
      
      # Output the defined gene features based on refseq.bed file#
      write.table(geneobj, genebody_name_file, col.names = F, row.names = F, quote = F)
      geneobj <- read.table(genebody_name_file)
      # check the geneobj file for promoter regions #
      geneobj <- refseqfile_check(geneobj, regionfile)
      # if featurewrite == FALSE then delete geneobj file #
      if(featurewrite == FALSE){
        file.remove(genebody_name_file)
      }
      
      return(geneobj)
    }
  }else if(bedfile=="cpgi"){
    if(feature == FALSE){
      regionfile <- Regionfile_read(bedfile, suffix)
      return(regionfile)
    }else{
      # region file is the reference gene or cpg bed file (without header) #
      cpgi_name_file<-paste0(bedfile, ".bed", suffix,"_Cpgifeature.txt")
      cpgiobj <- genomation::readFeatureFlank(paste(bedfile, ".bed", suffix, sep = ""), feature.flank.name = c("CpGisland", "Shores"))
      # Output the defined gene features based on refseq.bed file #
      write.table(cpgiobj, cpgi_name_file, col.names = F, row.names = F, quote = F)
      cpgiobj <- read.table(cpgi_name_file)
      cpgiobj[,6] <- cpgiobj[,2]
      cpgiobj <- cpgiobj[,3:6]
      colnames(cpgiobj) <- c("chr", "start", "end", "cpgfeature")
      # rename the "cpgi" #
      cpginum <- sum(cpgiobj$cpgfeature == "CpGisland")
      cpginame <- c(paste("cpgi", 1:cpginum, sep = ""), paste("shore", 1:(nrow(cpgiobj) - cpginum), sep = ""))
      cpgiobj <- data.frame(cpgi = cpginame, cpgiobj)
      # sort the file #
      cpgiobj <- arrange(cpgiobj, chr, start)
      # if featurewrite == FALSE then delete cpgiobj file #
      if(featurewrite == FALSE){
        file.remove(cpgi_name_file)
      }
      return(cpgiobj)
    }
  }else{
    # inform only refseq or cpgi file can be read #
    stop("Wrong bed file names and please use refseq or cpgi name")
  }
}




### MAIN ###

#check if the reference is empty (for example, there is no Alu or CpG for MT chromosome in human)
if(file.size(paste0(ref_path, "/refseq.bed_Transcript_chr",chr,"_",strand,".txt")) == 0L){
    print(paste0("Warning: reference file is empty: ", paste0(ref_path, "/refseq.bed_Transcript_chr",chr,"_",strand,".txt")))
    print("Results will be empty too!")
    regiongeneall <- data.frame()
    
}else{

    print("Loading gtf file")
    inputrefseqfile <- Bedfile_read2(paths=ref_path,feature=T,bedfile="refseq",suffix=paste0("_Transcript_chr",chr,"_",strand,".txt"))
    

    #control<-c(combin[2,c])
    #case<-c(combin[1,c])
    
    print(paste0("Case= ", case))
    print(paste0("Control= ", control))
    
    #subDir<-paste0(analysis_path,"Transcrit/",meth_type,"/",unlist(strsplit(case,"__"))[1],"_",unlist(strsplit(control,"__"))[1],"/")
    
    #if (!file.exists(subDir)){
    #    dir.create(subDir,recursive=T)
    #}
    
    print("Read the methylation files") #create table containing all the sample
    inputmethfile <- Methfile_read(control_paths = control, case_paths = case)
    
    print("Compute quality control for the input methylation file.")
    inputmethfile_QC <- Methfile_QC(inputmethfile)
    
    print("Calculate the methylation mean for regions")
    date()
    regiongeneall <- Methmean_region(inputmethfile_QC, inputrefseqfile, chrnum = chr,featurename=T)
    date()

}

write.table(regiongeneall,output_file, sep="\t", quote=FALSE, row.names=FALSE)
        




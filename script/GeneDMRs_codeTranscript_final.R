#script for transcrit DMR


args = commandArgs(trailingOnly=TRUE)
library(data.table)
library(GeneDMRs)
library(dplyr)


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
      genebody_name_file<-paste0(bedfile, ".bed", suffix,"_Genebody.txt")
      regionfile <- Regionfile_read(bedfile, suffix)
      
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


#check if the reference is empty (for example, there is no Alu or CpG for MT chromosome in human)
if(file.size(paste0(ref_path, "/refseq.bed_Transcript_chr",chr,"_",strand,".txt")) == 0L){
    print(paste0("Warning: reference file is empty: ", paste0(ref_path, "/refseq.bed_Transcript_chr",chr,"_",strand,".txt")))
    print("Results will be empty too!")
    regiongeneall <- data.frame()
    
}else{

    print("Loading gtf file")
    inputrefseqfile<-Bedfile_read2(paths=ref_path,feature=T,bedfile="refseq",suffix=paste0("_Transcript_chr",chr,"_",strand,".txt"))
    
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
    regiongeneall<-Methmean_region(inputmethfile_QC, inputrefseqfile, chrnum = chr,featurename=T)
    date()
}

write.table(regiongeneall,output_file, sep="\t", quote=FALSE, row.names=FALSE)
        




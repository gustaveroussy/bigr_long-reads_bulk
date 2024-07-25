#script for CpG DMR with gene


args = commandArgs(trailingOnly=TRUE)
library(data.table)
library(GeneDMRs)
library(dplyr)


analysis_path=args[1]

chr=args[2]
meth_type=args[3]
strand=args[4]
ref_path=args[5]
file_path<-paste0(analysis_path,"bed_files/chr_",chr,"/",meth_type,"/")
filelist<-list.files(path=file_path,pattern=paste0("*",strand,".bed"),full.names=F)


print("The file_path is:")
print(file_path)
print("for chromosome ")
print(chr)
print("and for strand:")
print(strand)
print("ref files:")
print(paste0(ref_path,"/cpgi.bed_ucsc_cpg_chr",chr,"_grch38.txt"))
print(paste0(ref_path,"../transcripts/refseq.bed_ucsc_transcrit_chr",chr,"_grch38_",strand,".txt"))

combin<-combn(filelist,2)
print("Combination:")
print(combin)
date()

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
      geneobj <-genomation::readTranscriptFeatures(paste(bedfile, ".bed", suffix, sep = ""))
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

print("Loading gtf files")
inputcpgifeaturefile<-Bedfile_read2(paths=ref_path,feature=T,bedfile="cpgi",suffix=paste0("_ucsc_cpg_chr",chr,"_grch38.txt"),featurewrite = FALSE)
inputrefseqfile<-Bedfile_read2(paths=paste0(ref_path,"../transcripts/"),feature=T,bedfile="refseq",suffix=paste0("_ucsc_transcrit_chr",chr,"_grch38_",strand,".txt"),featurewrite = FALSE)


for (c in (1:ncol(combin))){
        control<-c(combin[2,c])
        case<-c(combin[1,c])
        
	print(paste0("Case= ", case))
	print(paste0("Control= ", control))
	
        subDir<-paste0(analysis_path,"CpG/",meth_type,"/",unlist(strsplit(case,"__"))[1],"_",unlist(strsplit(control,"__"))[1],"/")
        
        if (!file.exists(subDir)){
		dir.create(subDir,recursive=T)
	}

        print(paste0("Comparison of ",unlist(strsplit(case,"__"))[1], " against ",unlist(strsplit(control,"__"))[1]))
        inputmethfile <- Methfile_read(control_paths = paste0(file_path,control), case_paths = paste0(file_path,case))

        # quality control #
        inputmethfile_QC <- Methfile_QC(inputmethfile)


	# methylation mean #
        print("methylation mean")
        date()
        regiongeneall<-Methmean_region(inputmethfile_QC,inputrefseqfile,inputcpgifeaturefile,chrnum = chr,featurename=T)
        date()
        write.table(regiongeneall,paste0(subDir,"mean_table_chr",chr,"_",strand,"_",unlist(strsplit(case,"__"))[1],"-",unlist(strsplit(control,"__"))[1],".tsv"),sep="\t",quote=F,row.names=F)
        
    print(paste0("Combination analysis number ",c, " is finished!"))
}




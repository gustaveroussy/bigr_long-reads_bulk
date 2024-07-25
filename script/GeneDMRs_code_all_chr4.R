#!/usr/bin/env Rscript
##############################################################
### USAGE: Rscript codeDMR.R input_path chr strand #numero du chromosome SEUL
##############################################################

args = commandArgs(trailingOnly=TRUE)
library(data.table)
library(GeneDMRs)
library(dplyr)
library(pheatmap)
library(grid)
library(EnhancedVolcano)
library(ggrepel)
library(gridExtra)
library(ggplot2)


    #boxplot
#png("/mnt/beegfs/scratch/bioinfo_core/B24015_CODU_01/data_output/methylation_analysis/concat/boxplot_all_chr.png",width=1000,height=1000)
#print(Sample_boxplot(inputmethfile_QC,inputrefseqfile))
#dev.off()


path=dirname(args[1])
input_name=basename(args[1])
print(input_name)
print(path)
name=unlist(strsplit(unlist(strsplit(x=input_name,split="all_chr_|.tsv"))[2],"-"))
print(name)


print("Opening the input file and formatting it")
regiongeneall<-fread(args[1])
regiongeneall<-regiongeneall[which(regiongeneall$start!="start"),]
if (grepl(pattern="cpg", x=tolower(path))){
	regiongeneall$Methgroup1_CpGisland<-as.numeric(regiongeneall$Methgroup1_CpGisland)
	regiongeneall$Methgroup2_CpGisland<-as.numeric(regiongeneall$Methgroup2_CpGisland)
        regiongeneall$Methgroup1_Shore<-as.numeric(regiongeneall$Methgroup1_Shore)
	regiongeneall$Methgroup2_Shore<-as.numeric(regiongeneall$Methgroup2_Shore)
        regiongeneall$Readgroup1_CpGisland<-as.numeric(regiongeneall$Readgroup1_CpGisland)
        regiongeneall$Readgroup2_CpGisland<-as.numeric(regiongeneall$Readgroup2_CpGisland)
	regiongeneall$Readgroup1_Shore<-as.numeric(regiongeneall$Readgroup1_Shore)
        regiongeneall$Readgroup2_Shore<-as.numeric(regiongeneall$Readgroup2_Shore)
	regiongeneall<-na.omit(regiongeneall)
}else{
	regiongeneall$Methgroup1<-as.numeric(regiongeneall$Methgroup1)
	regiongeneall$Methgroup2<-as.numeric(regiongeneall$Methgroup2)
	regiongeneall$Readgroup1<-as.numeric(regiongeneall$Readgroup1)
	regiongeneall$Readgroup2<-as.numeric(regiongeneall$Readgroup2)
}
    #statistical test
print("Statistical test")
regiongeneall_Qvalue<-Logic_regression(regiongeneall)
write.table(regiongeneall_Qvalue,paste0(path,"/regiongeneall_stat_table_",name[1],"_vs_",name[2],".csv"),sep=",",quote=F,row.names=F)

    #significant feature filter at 10%
print("Significant filter")
regiongeneall_significant<-Significant_filter(regiongeneall_Qvalue)
write.table(regiongeneall_significant,paste0(path,"/regiongeneall_significant_",name[1],"_vs_",name[2],".csv"),sep=",",quote=F,row.names=F)

head(regiongeneall_significant)

#reverse methylation
#print("Reversing methylation")
#regiongeneall_significant$Methdiff1<-regiongeneall_significant$Methdiff1*(-1)

print("1")

############################################################## Graphic code

####pHeatmap corrected function


Heatmap_plot<-function(regiongeneall_significant,group,level,featurename = NULL,title = "Methylation level",
                         display_numbers = FALSE, number_format = "%.0f", cluster_rows = FALSE, cluster_cols = TRUE,
                                                 gaps_row = c(2,1), gaps_col = NULL){
    groupnum <- length(grep("group", colnames(regiongeneall_significant))) / 2
    methheat <- data.frame(regiongeneall_significant[, c("Methgroup1","Methgroup2")])
    #rownames(methheat) <- paste(regiongeneall_significant$transcript_name,regiongeneall_significant$feature,regiongeneall_significant$gene_name, sep = "_")
    rownames(methheat) <- regiongeneall_significant$label
    colnames(methheat) <- c(group)
    annotation_row = data.frame(chr = regiongeneall_significant$chr)
    rownames(annotation_row) <- rownames(methheat)
    annotation_col = data.frame(Group = rownames(methheat))
    rownames(annotation_col) <- rownames(methheat)
    pheatmap(methheat*100, main = paste0(title, " of ", level, "methylated region", "\n", group[2], " vs ", group[1]), display_numbers = display_numbers,
	 number_format = number_format, cluster_rows = T,cluster_cols = F,silent=T)
}

#rotate x label by 45°
draw_colnames_45 <- function (coln, gaps, ...) {
    coord = pheatmap:::find_coordinates(length(coln), gaps)
    x = coord$coord - 0.5 * coord$size
    res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 1, hjust = 1, rot = 45, gp = gpar(...))
    return(res)
}


############DMR per chromosome
geom_rrect <- function(mapping = NULL, data = NULL, # nocov start
                       stat = "identity", position = "identity",
                       radius = grid::unit(6, "pt"),
                       ...,
                       na.rm = FALSE,
                       show.legend = NA,
                       inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomRrect,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      radius = radius,
      na.rm = na.rm,
      ...
    )
  )
}
GeomRrect <- ggplot2::ggproto(
  "GeomRrect", ggplot2::Geom,
  default_aes = ggplot2::aes(
    colour = NA, fill = "grey35", size = 0.5, linetype = 1, alpha = NA
  ),
  required_aes = c("xmin", "xmax", "ymin", "ymax"),
  draw_panel = function(self, data, panel_params, coord,
                        radius = grid::unit(6, "pt")) {
    coords <- coord$transform(data, panel_params)
    lapply(1:length(coords$xmin), function(i) {
      grid::roundrectGrob(
        coords$xmin[i], coords$ymax[i],
        width = (coords$xmax[i] - coords$xmin[i]),
        height = (coords$ymax[i] - coords$ymin)[i],
        r = radius,
        default.units = "native",
        just = c("left", "top"),
        gp = grid::gpar(
          col = coords$colour[i],
          fill = alpha(coords$fill[i], coords$alpha[i]),
          lwd = coords$size[i] * .pt,
          lty = coords$linetype[i],
          lineend = "butt"
        )
      )
    }) -> gl
    grobs <- do.call(grid::gList, gl)
  },
  draw_key = ggplot2::draw_key_polygon
)

############################################################
print("2")
if(grepl(pattern="alu", x=tolower(path))){
	############ piechart of Alu type
	print("Performing Alu analysis")
	fichier<-regiongeneall_significant
	#keep only first 4 characters of Alu_name
	fichier$Alu_types<-gsub(substr(fichier$id,start=1,stop=4),pattern="_",replacement="")
	fichier$level<-ifelse(fichier$Methdiff1 < 0, "hypo", "hyper" )
	
	count<-as.data.table(table(fichier[,c(12,13)]))
	colnames(count)<-c("Alu sequence","level","value")
	
	count_hyper<-count[which(count$level=="hyper"),]%>%mutate(csum=rev(cumsum(rev(value))),
        pos=value/2+lead(csum,1),
        pos=if_else(is.na(pos),value/2,pos))
	count_hypo<-count[which(count$level=="hypo"),]%>%mutate(csum=rev(cumsum(rev(value))),
        pos=value/2+lead(csum,1),
        pos=if_else(is.na(pos),value/2,pos))
	
	count_hypo$total<-sum(count_hypo$value)
	count_hyper$total<-sum(count_hyper$value)
	
	hypo<-ggplot(count_hypo,aes(x=" ",y=value, fill =`Alu sequence`))+geom_col(width=1,color=1)+coord_polar(theta="y",start=0)+scale_fill_brewer(palette="Pastel1")+geom_label_repel(data=count_hypo,
		aes(y=pos,label=paste0(round((value/total)*100,digit=1),"%") ),size=7,nudge_x = 1,show.legend=F)+theme_void()+ggtitle(paste0("  Hypomethylated Alu regions in ",name[1]," vs ",name[2]),subtitle=paste0("   n = ",sum(count_hypo$value)))+
		theme(text = element_text(size=13))
	hyper<-ggplot(count_hyper,aes(x=" ",y=value,fill=`Alu sequence`))+geom_col(width=1,color=1)+coord_polar(theta="y",start=0)+scale_fill_brewer(palette="Pastel1")+geom_label_repel(data=count_hyper,
		aes(y=pos,label=paste0(round((value/total)*100,digit=1),"%") ),size=7,nudge_x = 1,show.legend=F)+theme_void()+ggtitle(paste0("  Hypermethylated Alu regions in ",name[1]," vs ",name[2]),subtitle=paste0("   n = ",sum(count_hyper$value)))+
		theme(text = element_text(size=13))
	
	print(paste0(path,"/piechart_significant_Alu_0.05_",name[1],"_vs_",name[2],".png"))
	png(paste0(path,"/piechart_significant_Alu_0.05_",name[1],"_vs_",name[2],".png"),width=1200,heigh=500)
	grid.arrange(hypo, hyper, ncol=2, nrow = 1)
	dev.off()
	
	###################enhanced volcano
	fichier_Qvalue<-regiongeneall_Qvalue
	fichier_Qvalue$label<-paste(gsub(fichier_Qvalue$id,pattern="_.*",replacement=""),fichier_Qvalue$start,fichier_Qvalue$end,sep="_")

	top_5<-unlist(top_n(fichier_Qvalue,-5,Qvalue1)[,12])
	
	print(paste0(path,"/Volcano_plot_Aluseq_0.05_",name[1],"_vs_",name[2],".png"))
	
	png(paste0(path,"/Volcano_plot_Aluseq_0.05_",name[1],"_vs_",name[2],".png"),width=1200,heigh=800)
	print(EnhancedVolcano(fichier_Qvalue,lab=fichier_Qvalue$label,x="Methdiff1",y="Qvalue1",selectLab=top_5,FCcutoff=0.05,pCutoff=0.05,xlab="Differential methylation Level",xlim=c(-1,1),title="Volcano plot of significant Alu sequences", subtitle="Cutoff of adjusted p-value at 0.05(5%) and differential methylation level at 0.05(5%)",legendLabels=c("NS","MethDiff","P-value","P-value and Methdiff"),legendPosition = 'right',drawConnectors = TRUE,gridlines.major = FALSE,gridlines.minor = FALSE))
	dev.off()
	
	##################pheatmap
	fichier<-regiongeneall_significant
	fichier$label<-paste(gsub(fichier$id,pattern="_.*",replacement=""),fichier$start,fichier$end,sep="_")
	hyper<-fichier[which(fichier$Methdiff1 >0),]
	hypo<-fichier[which(fichier$Methdiff1 <0),]
	
	top_30_hyper<-top_n(hyper,30,Methdiff1)
	top_30_hypo<-top_n(hypo,-30,Methdiff1)
	
	heatmap_hyper<-Heatmap_plot(top_30_hyper,c(name[2],name[1]),"hyper")
	heatmap_hypo<-Heatmap_plot(top_30_hypo,c(name[2],name[1]),"hypo")
	
	print(paste0(path,"/heatmap_top30_hypo_hyper_significant_0.05_",name[1],"_vs_",name[2],".png"))
	png(paste0(path,"/heatmap_top30_hypo_hyper_significant_0.05_",name[1],"_vs_",name[2],".png"),width=1200,heigh=800)
	grid.arrange(grobs=list(heatmap_hypo[[4]], heatmap_hyper[[4]]),ncol=2,nrow=1)
	dev.off()
	
	#dmr par chromosome
	
}


if(grepl(pattern="transcrit", x=tolower(path))){
	#Correspondance to Gene ID
	print("Transcrit analysis")
	print("Merging with corresponding gene ID")
	corr_gene<-unique(fread("/mnt/beegfs/userdata/y_mesloub/ressources/references/grch38/correspondance_table_ENST_ENSG.tsv"))
	colnames(corr_gene)<-c("gene_id","gene_name","id","transcript_name")
	
	#all_features
	regiongeneall_id<-left_join(x=regiongeneall_Qvalue,y=corr_gene,by=c("id"))
	write.table(regiongeneall_id,paste0(path,"/regiongeneall_with_SYMBOL_",name[1],"_vs_",name[2],".csv"),sep=",",quote=F,row.names=F)
	
	#significant features
	regiongeneall_sign<-regiongeneall_id[which(regiongeneall_id$Qvalue1 <=0.05),]
	write.table(regiongeneall_sign,paste0(path,"/regiongeneall_significant_0.05_",name[1],"_vs_",name[2],".csv"),sep=",",quote=F,row.names=F)
	
	################piechart
	print("piechart")
	fichier<-regiongeneall_sign
	fichier$feature_types<-gsub(fichier$feature,pattern="_.*",replacement="")
	fichier$level<-ifelse(fichier$Methdiff1 < 0,"hypo","hyper")
	
	count<-as.data.table(table(fichier[,c(16,17)]))
	colnames(count)<-c("Regulatory region","level","value")
	
	count_hyper<-count[which(count$level=="hyper"),]%>%mutate(csum=rev(cumsum(rev(value))),
        	pos=value/2+lead(csum,1),
        	pos=if_else(is.na(pos),value/2,pos))
	count_hypo<-count[which(count$level=="hypo"),]%>%mutate(csum=rev(cumsum(rev(value))),
        	pos=value/2+lead(csum,1),
        	pos=if_else(is.na(pos),value/2,pos))
	
	count_hypo$total<-sum(count_hypo$value)
	count_hyper$total<-sum(count_hyper$value)
	
	hypo<-ggplot(count_hypo, aes(x = " " , y = value, fill = `Regulatory region`)) +  geom_col(width = 1,color=1) +  coord_polar(theta = "y",start=0) + 
		scale_fill_brewer(palette = "Pastel1")+	geom_label_repel(data = count_hypo,aes(y = pos, label = paste0(round((value/total)*100,digit=1),"%") ), size = 7,
		nudge_x = 1, show.legend = F)+theme_void()+ggtitle(paste0("  Hypomethylated regions in ",name[1]," vs ",name[2]),subtitle=paste0("   n = ",sum(count_hypo$value)))+theme(text = element_text(size=13))
	
	hyper<-ggplot(count_hyper, aes(x = " " , y = value, fill = `Regulatory region`)) +  geom_col(width = 1,color=1) +  coord_polar(theta = "y",start=0) + 
		scale_fill_brewer(palette = "Pastel1")+	geom_label_repel(data = count_hyper,aes(y = pos, label = paste0(round((value/total)*100,digit=1),"%") ), size = 7, 
		nudge_x = 1, show.legend = F)+ theme_void()+ggtitle(paste0("  Hypermethylated regions in ",name[1]," vs ",name[2]),subtitle=paste0("   n = ",sum(count_hyper$value)))+theme(text = element_text(size=13))
	
	
	png(paste0(path,"/piechart_significant_features_0.05_",name[1],"_vs_",name[2],".png"),width=1200,heigh=500)
	grid.arrange(hypo, hyper, ncol=2, nrow = 1)
	dev.off()
	
	####################enhanced volcano
	fichier<-regiongeneall_id
        fichier$label<-paste(fichier$transcript_name,fichier$feature,sep="_")
	
        top_5<-unlist(top_n(fichier,-5,Qvalue1)[,16])
	
	print("volcanoplot")
        png(paste0(path,"/Volcano_plot_transcript_0.05_",name[1],"_vs_",name[2],".png"),width=1200,heigh=800)
        print(EnhancedVolcano(fichier,lab=fichier$label,x="Methdiff1",y="Qvalue1",selectLab=top_5,FCcutoff=0.05,pCutoff=0.05,xlab="Differential methylation Level",xlim=c(-1,1),title="Volcano plot of significant region", subtitle="Cutoff of adjusted p-value at 0.05(5%) and differential methylation level at 0.1(10%)",legendLabels=c("NS","MethDiff","P-value","P-value and Methdiff"),legendPosition = 'right',drawConnectors = TRUE,gridlines.major = FALSE,gridlines.minor = FALSE))
	dev.off()
		
	###################pheatmap
	print("heatmap")	
	fichier<-regiongeneall_sign
	fichier$label<-paste(fichier$transcript_name,fichier$feature,fichier$gene_name,sep="_")
        hyper<-fichier[which(fichier$Methdiff1 >0),]
        hypo<-fichier[which(fichier$Methdiff1 <0),]
	
        top_30_hyper<-top_n(hyper,30,Methdiff1)
        top_30_hypo<-top_n(hypo,-30,Methdiff1)
	
        heatmap_hyper<-Heatmap_plot(top_30_hyper,c(name[1],name[2]),"hyper")
        heatmap_hypo<-Heatmap_plot(top_30_hypo,c(name[1],name[2]),"hypo")

        png(paste0(path,"/heatmap_top30_hypo_hyper_significant_0.05_",name[1],"_vs_",name[2],".png"),width=1200,heigh=800)
        grid.arrange(grobs=list(heatmap_hypo[[4]], heatmap_hyper[[4]]),ncol=2,nrow=1)
        dev.off()
	
	#dmr par chromosome
}

if(grepl(pattern="cpg", x=tolower(path))){
	print("it is CpG analysis")
	#only on CpG for now, adding CpG Shore+shelf later
	#adding transcript and gene name
	
	print("Merging with corresponding gene ID")
        corr_gene<-unique(fread("/mnt/beegfs/userdata/y_mesloub/ressources/references/grch38/correspondance_table_ENST_ENSG.tsv"))
        colnames(corr_gene)<-c("gene_id","gene_name","id","transcript_name")

        #all_features
        regiongeneall_id<-left_join(x=regiongeneall_Qvalue,y=corr_gene,by=c("id"))
        write.table(regiongeneall_id,paste0(path,"/regiongeneall_with_SYMBOL_",name[1],"_vs_",name[2],".csv"),sep=",",quote=F,row.names=F)

        #significant features
        regiongeneall_sign<-left_join(x=regiongeneall_significant,y=corr_gene,by=c("id"))
        write.table(regiongeneall_sign,paste0(path,"/regiongeneall_significant_0.05_",name[1],"_vs_",name[2],".csv"),sep=",",quote=F,row.names=F)
	
	#piechart
	#enhanced volcano
	#label: CpG_feature_transcrit
	#heatmap
	#dmr par chromosome
}

print("Finished")

#!/usr/bin/env Rscript
##############################################################
### USAGE: Rscript GeneDMRs_code_all_chr_final.R input_mean_table ref_to_convert_ensemblID
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
library(patchwork)


## Get parameters
#input_file  <- "/mnt/beegfs/scratch/n_rabearivelo/test_pipeline/methylation/data_output2/methylation_analysis/CpG/5mCG/1033_CD14_vs_1033_CD34/mean_table_1033_CD14_vs_1033_CD34_all_chr.tsv"
#input_file  <- "/mnt/beegfs/scratch/n_rabearivelo/test_pipeline/methylation/data_output2/methylation_analysis/Transcript/5mCG/1033_CD14_vs_1033_CD34/mean_table_1033_CD14_vs_1033_CD34_all_chr.tsv"
#input_file  <- "/mnt/beegfs/scratch/n_rabearivelo/test_pipeline/methylation/data_output2/methylation_analysis/Alu/5mCG/1033_CD14_vs_1033_CD34/mean_table_1033_CD14_vs_1033_CD34_all_chr.tsv"
#reference_conversion <- "/mnt/beegfs/database/bioinfo/bigr_long-reads_bulk/REFERENCES/DMR/GRCh38/correspondance_table_ENST_ENSG.tsv"
#chromSize <- "/mnt/beegfs/scratch/n_rabearivelo/test_pipeline/methylation/data_output2/resources/chromSize.txt"
input_file  <- args[1]
reference_conversion <- args[2]
chromSize <- args[3]
Qvalue_threshold <- 0.05

## Get missing parameters
path <- dirname(input_file)
name <- unlist(strsplit(gsub("mean_table_|_all_chr.tsv", "", basename(input_file)), "_vs_"))

## Print parameters
print(paste0("input: ", input_file))
print(paste0("reference conversion file: ", reference_conversion))
print(paste0("sample or condition names: ", name))

## Load and format data
print("Opening the input file and formatting it")
input_meth_table <- fread(input_file, blank.lines.skip = TRUE) # read and remove empty lines
input_meth_table <- input_meth_table[which(input_meth_table$start!="start"),] #remove header duplicated lines
input_meth_table$start <- as.numeric(input_meth_table$start)
input_meth_table$end <- as.numeric(input_meth_table$end)
if (grepl(pattern = "cpg", x = tolower(path))){
	input_meth_table$Methgroup1_CpGisland <- as.numeric(input_meth_table$Methgroup1_CpGisland)
	input_meth_table$Methgroup2_CpGisland <- as.numeric(input_meth_table$Methgroup2_CpGisland)
    input_meth_table$Methgroup1_Shore <- as.numeric(input_meth_table$Methgroup1_Shore)
	input_meth_table$Methgroup2_Shore <- as.numeric(input_meth_table$Methgroup2_Shore)
    input_meth_table$Readgroup1_CpGisland <- as.numeric(input_meth_table$Readgroup1_CpGisland)
    input_meth_table$Readgroup2_CpGisland <- as.numeric(input_meth_table$Readgroup2_CpGisland)
	input_meth_table$Readgroup1_Shore <- as.numeric(input_meth_table$Readgroup1_Shore)
    input_meth_table$Readgroup2_Shore <- as.numeric(input_meth_table$Readgroup2_Shore)
	input_meth_table <- na.omit(input_meth_table)
}else{
	input_meth_table$Methgroup1 <- as.numeric(input_meth_table$Methgroup1)
	input_meth_table$Methgroup2 <- as.numeric(input_meth_table$Methgroup2)
	input_meth_table$Readgroup1 <- as.numeric(input_meth_table$Readgroup1)
	input_meth_table$Readgroup2 <- as.numeric(input_meth_table$Readgroup2)
}

if(dim(input_meth_table)[1] > 0){

    ## Make statistical test
    print("Statistical test")
    meth_res <- Logic_regression(input_meth_table)
    
    ## Make filtering on Qvalue
    print("Significant filter")
    #meth_res_significant<-Significant_filter(meth_res)
    meth_res_significant <- meth_res[meth_res$Qvalue1 <= Qvalue_threshold,]
    
    head(meth_res_significant)
    
    
    ## Functions to make graphs and save tables
    
    # pheatmap corrected function
    Heatmap_plot <- function(meth_res_significant, group, level, featurename = NULL, title = "Methylation level",
                             display_numbers = FALSE, number_format = "%.0f", cluster_rows = FALSE, cluster_cols = TRUE,
                             gaps_row = c(2,1), gaps_col = NULL){
        groupnum <- length(grep("group", colnames(meth_res_significant))) / 
        if("Methgroup1" %in% colnames(meth_res_significant)){
            methheat <- data.frame(meth_res_significant[, c("Methgroup1","Methgroup2")])
        }else{ 
            methheat <- data.frame(meth_res_significant[, c("Methgroup1_CpGisland","Methgroup2_CpGisland")])
        }
        #rownames(methheat) <- paste(meth_res_significant$transcript_name,meth_res_significant$feature,meth_res_significant$gene_name, sep = "_")
        rownames(methheat) <- meth_res_significant$label
        colnames(methheat) <- c(group)
        annotation_row <- data.frame(chr = meth_res_significant$chr)
        rownames(annotation_row) <- rownames(methheat)
        annotation_col <- data.frame(Group = rownames(methheat))
        rownames(annotation_col) <- rownames(methheat)
        pheatmap(methheat*100, main = paste0(title, " of ", level, "methylated region", "\n", group[2], " vs ", group[1]), display_numbers = display_numbers,
    	         number_format = number_format, cluster_rows = T,cluster_cols = F,silent=T)
    }
    
    # rotate x label by 45°
    draw_colnames_45 <- function (coln, gaps, ...) {
        coord = pheatmap:::find_coordinates(length(coln), gaps)
        x = coord$coord - 0.5 * coord$size
        res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 1, hjust = 1, rot = 45, gp = gpar(...))
        return(res)
    }
    
    #GeomRrect for DMR per chromosome
    element_roundrect <- function(fill = NULL, colour = NULL, size = NULL,
      linetype = NULL, color = NULL, r=grid::unit(0.1, "snpc"), inherit.blank = FALSE) {
      if (!is.null(color))  colour <- color
      if (!grid::is.unit(r)) r <- grid::unit(r, 'snpc')
      structure(
        list(fill = fill, 
             colour = colour, 
             size = size, 
             linetype = linetype,
             r = r,
             inherit.blank = inherit.blank),
        class = c("element_roundrect", "element_rect", "element")
      )
    }
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
    
    
    ## Make graphs and save tables
    print("Make graphs and save tables")
    
    if(grepl(pattern="/Alu/", x = path)){
        print("Alu analysis")
        
        # Save tables
        write.table(meth_res,paste0(path,"/meth_res_table_",name[1],"_vs_",name[2],"_all_chr.csv"),sep=",",quote=F,row.names=F)
        if(dim(meth_res_significant)[1]>0) {
            write.table(meth_res_significant,paste0(path,"/meth_res_significant_",name[1],"_vs_",name[2],"_all_chr.csv"),sep=",",quote=F,row.names=F)
        }else print("No differential methylation under Qvalue threshold")
        
        # Enhanced volcanoplot
    	print("Draw Volcano plot")
        meth_res_volca <- meth_res
    	#meth_res_volca$label <- paste(gsub(meth_res_volca$id,pattern="_.*",replacement=""),meth_res_volca$chr, meth_res_volca$start,meth_res_volca$end,sep="_")
    	#top_5 <- unlist(top_n(meth_res_volca,-5,Qvalue1)[,"label"])
    	png(paste0(path,"/Volcano_plot_Alu_",name[1],"_vs_",name[2],".png"),width=1200,heigh=800)
    	print(EnhancedVolcano(meth_res_volca, x = "Methdiff1", y = "Qvalue1", lab = NULL, FCcutoff = 0.05, pCutoff = Qvalue_threshold, xlab = "Differential methylation Level", ylab = bquote(~Log[10]~ 'Adusted Pvalue'), xlim = c(0,1), title = "Volcano plot of Alu regions", subtitle = paste0("Cutoff of adjusted p-value at ",Qvalue_threshold,"(",Qvalue_threshold*100,"%) and differential methylation level at 0.05(5%)"), legendLabels = c("NS","DiffMeth ","Adjusted P-value","Adjusted P-value and DiffMeth"), legendPosition = 'right', drawConnectors = TRUE, gridlines.major = FALSE, gridlines.minor = FALSE)) #, lab = meth_res_volca$label, selectLab = top_5
    	dev.off()
        	
        # MA-plot
    	print("Draw MA-plot")
        meth_res_ma <- meth_res
        #meth_res_ma$label <- paste(gsub(meth_res_ma$id,pattern="_.*",replacement=""),meth_res_ma$chr, meth_res_ma$start,meth_res_ma$end,sep="_")
        meth_res_ma$baseMean <- 1 / (10^(rowMeans(meth_res_ma[,c("Methgroup1","Methgroup2")]))) #EnhancedVolcano will internally transform the y-axis values by -log10(); here, we counteract this by creating a new column in the data for base mean.
        keyvals <- ifelse(abs(meth_res_ma$Methdiff1) > 0.05 & meth_res_ma$Qvalue1 < Qvalue_threshold, 'red2',
                   ifelse(abs(meth_res_ma$Methdiff1) > 0.05 & meth_res_ma$Qvalue1 > Qvalue_threshold, 'forestgreen',
                   ifelse(abs(meth_res_ma$Methdiff1) < 0.05 & meth_res_ma$Qvalue1 < Qvalue_threshold, 'royalblue',
                   'grey30')))
        keyvals[is.na(keyvals)] <- 'grey30'
        names(keyvals)[keyvals == 'red2'] <- 'Adjusted P-value and Differential methylation Level'
        names(keyvals)[keyvals == 'forestgreen'] <- 'Differential methylation Level'
        names(keyvals)[keyvals == 'royalblue'] <- 'Adjusted P-value'
        names(keyvals)[keyvals == 'grey30'] <- 'NS'
    	png(paste0(path,"/MA_plot_Alu_",name[1],"_vs_",name[2],".png"),width=1200,heigh=800)
    	print(EnhancedVolcano(meth_res_ma, x = "Methdiff1", y = "baseMean", lab = NULL, colCustom = keyvals, FCcutoff = 0.05, pCutoff = Qvalue_threshold, xlab = "Differential methylation Level", ylab = "Mean methylation Level", xlim = c(-1,1), ylim = c(0,1), title = 'MA plot of Alu regions', subtitle = paste0("Cutoff of adjusted p-value at ",Qvalue_threshold,"(",Qvalue_threshold*100,"%) and differential methylation level at 0.05(5%)"), legendPosition = 'right', drawConnectors = TRUE, gridlines.major = FALSE, gridlines.minor = FALSE)+ coord_flip()) #, lab = meth_res_ma$label, selectLab = top_5
    	dev.off()
    
        if(dim(meth_res_significant)[1]>0){
    
        	# Piechart
    	    print("Draw Piechart")
        	meth_res_sign_pie <- meth_res_significant
        	meth_res_sign_pie$feature_types <- gsub(substr(meth_res_sign_pie$id,start=1,stop=4),pattern="_",replacement="") #keep only first 4 characters of Alu_name #TO DO
        	meth_res_sign_pie$level <- ifelse(meth_res_sign_pie$Methdiff1 == 0, "No difference", ifelse(meth_res_sign_pie$Methdiff1 < 0, "hypo", "hyper"))
        	
        	count <- as.data.table(table(meth_res_sign_pie[,c("feature_types","level")]))
        	colnames(count) <- c("Alu regions","level","value")
        	
        	count_hyper <- count[which(count$level == "hyper"),] %>% mutate(csum = rev(cumsum(rev(value))),
                            pos = value/2+lead(csum,1),
                            pos = if_else(is.na(pos), value/2,pos))
        	count_hypo <- count[which(count$level == "hypo"),] %>% mutate(csum = rev(cumsum(rev(value))),
                            pos = value/2+lead(csum,1),
                            pos = if_else(is.na(pos), value/2,pos))
        	
        	count_hypo$total <- sum(count_hypo$value)
        	count_hyper$total <- sum(count_hyper$value)
        	
        	if(dim(count_hypo)[1] > 0 | dim(count_hyper)[1] > 0){
                if(dim(count_hypo)[1] > 0){
                    hypo <- ggplot(count_hypo, aes(x = " " , y = value, fill = `Alu regions`)) + 
                                geom_col(width = 1, color=1) +
                                coord_polar(theta = "y", start = 0) + 
                              scale_fill_brewer(palette = "Dark2") +	
                              geom_label_repel(data = count_hypo, aes(y = pos, label = paste0(round((value/total)*100,digit=1),"%")), size = 5, nudge_x = 1, show.legend = F) +
                              theme_void() +
                              ggtitle(paste0("  Hypomethylated regions in ", name[1], " vs ", name[2]),subtitle=paste0("   n = ", sum(count_hypo$value))) +
                              theme(text = element_text(size = 14))
                  } else hypo <- grid::textGrob('No hypomethylated Alu region.')
                  if(dim(count_hyper)[1] > 0){
                    hyper <- ggplot(count_hyper, aes(x = " " , y = value, fill = `Alu regions`)) +
                                geom_col(width = 1, color = 1) +
                                coord_polar(theta = "y", start = 0) +
                              scale_fill_brewer(palette = "Dark2") +
                              geom_label_repel(data = count_hyper, aes(y = pos, label = paste0(round((value/total)*100, digit = 1), "%") ), size = 5, nudge_x = 1, show.legend = F) +
                              theme_void() +
                              ggtitle(paste0("  Hypermethylated regions in ", name[1], " vs ", name[2]), subtitle = paste0("   n = ", sum(count_hyper$value))) +
                              theme(text = element_text(size = 14))
                } else hyper <- grid::textGrob('No hypermethylated Alu region.')
            	png(paste0(path,"/Piechart_significant_Alu_",Qvalue_threshold,"_",name[1],"_vs_",name[2],".png"),width=1200,heigh=500)
            	grid.arrange(hypo, hyper, ncol=2, nrow = 1)
            	dev.off()
            }
        	
        	# Heatmap
    	    print("Draw Heatmap")
        	meth_res_sign_heat <- meth_res_significant
        	meth_res_sign_heat$label <- paste(gsub(meth_res_sign_heat$id,pattern="_.*",replacement=""), meth_res_sign_heat$chr, meth_res_sign_heat$start, meth_res_sign_heat$end, sep="_")
        	hyper <- meth_res_sign_heat[which(meth_res_sign_heat$Methdiff1 > 0),]
        	hypo <- meth_res_sign_heat[which(meth_res_sign_heat$Methdiff1 < 0),]
        	
        	top_30_hyper <- top_n(hyper,30,Methdiff1)
        	top_30_hypo <- top_n(hypo,-30,Methdiff1)
        	
        	if(dim(top_30_hypo)[1] > 1 | dim(top_30_hyper)[1] > 1){
        	    if(dim(top_30_hypo)[1] > 1){
        	        heatmap_hypo <- Heatmap_plot(top_30_hypo, c(name[1], name[2]), "hypo")[[4]]
        	    } else heatmap_hypo <- grid::textGrob('No (enough) hypomethylated Alu region to plot.')
        	    if(dim(count_hyper)[1] > 1){
                    heatmap_hyper <- Heatmap_plot(top_30_hyper, c(name[1], name[2]), "hyper")[[4]]
                } else heatmap_hyper <- grid::textGrob('No (enough) hypermethylated Alu region to plot.')
            	png(paste0(path, "/Heatmap_significant_Alu_top30_hypo_hyper_", name[1],"_vs_", name[2], ".png"), width = 1200, heigh = 800)
            	grid.arrange(grobs = list(heatmap_hypo, heatmap_hyper), ncol = 2, nrow = 1)
        	    dev.off()
        	}
        	
        	# DMR par chromosome
    	    print("DMR par chromosome")	
          	meth_res_sign_DRMpc <- meth_res_significant
            meth_res_sign_DRMpc$level <- 0.25
            meth_res_sign_DRMpc$level[which(meth_res_sign_DRMpc$Methdiff1 < 0)] <- -0.25
            meth_res_sign_DRMpc$group <- "hyper"
            meth_res_sign_DRMpc$group[which(meth_res_sign_DRMpc$Methdiff1 < 0)] <- "hypo"
            meth_res_sign_DRMpc$point_shape <- 25
            meth_res_sign_DRMpc$point_shape[which(meth_res_sign_DRMpc$Methdiff1 < 0)] <- 24
              
            size_chr <- fread(chromSize)
            chr <- unique(sort(meth_res_sign_DRMpc$chr))
            p <- list()
            for (c in chr){
              	print(c)
              	assign(x = paste0("size_", c), size_chr[which(size_chr$V1 == c),])
              	p[[paste0(c)]] <- ggplot(meth_res_sign_DRMpc[which(meth_res_sign_DRMpc$chr == c),], aes(x = start , y = level, color = group)) +
              	geom_point(alpha = 0.5, aes(size = 10)) +
              	geom_segment(x = 0, y = 0, xend = as.integer(eval(parse(text = paste0("size_", c)))$V2), yend = 0, size = 10, color = "grey70") +
              	xlim(0, as.integer(max(size_chr$V2)) + 10000000) +
              	ylim(-0.50, 0.30) +
              	theme(axis.text.y = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(),  axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), legend.title = element_blank(), panel.background = element_blank(), legend.position = "none") +
              	geom_label(x = (as.integer(max(size_chr$V2)) + 10000000) , y = 0, label = paste0("Chr ", c), color = "black", size = 12)
            }
            
            png(paste0(path, "/DMR_per_chromosome_significant_Alu_", name[1], "_vs_", name[2], ".png"), width = 2000, heigh = 3000)
            grid.arrange(grobs = p, ncol = 1, nrow = length(chr))
            dev.off()
        	
        }
        
    }else if(grepl(pattern="/Transcript/", x = path)){
        print("Transcript analysis")
        	
    	# Correspondance to Gene ID
    	print("Merging with corresponding gene ID")
    	corr_gene <- unique(fread(reference_conversion))
    	colnames(corr_gene) <- c("gene_id", "gene_name", "id", "transcript_name")
    	meth_res <- left_join(x = meth_res, y = corr_gene, by = c("id"))
        meth_res_significant <- left_join(x = meth_res_significant, y = corr_gene, by = c("id"))
        
        # Save tables
        write.table(meth_res, paste0(path,"/meth_res_table_",name[1],"_vs_",name[2],"_all_chr.csv"),sep=",",quote=F,row.names=F)
        if(dim(meth_res_significant)[1]>0) {
            write.table(meth_res_significant, paste0(path,"/meth_res_significant_",name[1],"_vs_",name[2],"_all_chr.csv"),sep=",",quote=F,row.names=F)
        }else print("No differential methylation under Qvalue threshold")
    
        # Enhanced volcanoplot
    	print("Draw Volcano plot")
        meth_res_volca <- meth_res
    	#meth_res_volca$label <- paste(meth_res_volca$transcript_name,meth_res_volca$feature,sep="_")
    	#top_5 <- unlist(top_n(meth_res_volca,-5,Qvalue1)[,"label"])
    	png(paste0(path,"/Volcano_plot_Transcript_",name[1],"_vs_",name[2],".png"),width=1200,heigh=800)
    	print(EnhancedVolcano(meth_res_volca, x = "Methdiff1", y = "Qvalue1", lab = NULL, FCcutoff = 0.05, pCutoff = Qvalue_threshold, xlab = "Differential methylation Level", ylab = bquote(~Log[10]~ 'Adusted Pvalue'), xlim = c(0,1), title = "Volcano plot of Transcript regions", subtitle = paste0("Cutoff of adjusted p-value at ",Qvalue_threshold,"(",Qvalue_threshold*100,"%) and differential methylation level at 0.05(5%)"), legendLabels = c("NS","DiffMeth ","Adjusted P-value","Adjusted P-value and DiffMeth"), legendPosition = 'right', drawConnectors = TRUE, gridlines.major = FALSE, gridlines.minor = FALSE)) #, lab = meth_res_volca$label, selectLab = top_5
    	dev.off()
    
        # MA-plot
    	print("Draw MA-plot")
        meth_res_ma <- meth_res
    	#meth_res_ma$label <- paste(meth_res_ma$transcript_name,meth_res_ma$feature,sep="_")
        meth_res_ma$baseMean <- 1 / (10^(rowMeans(meth_res_ma[,c("Methgroup1","Methgroup2")]))) #EnhancedVolcano will internally transform the y-axis values by -log10(); here, we counteract this by creating a new column in the data for base mean.
        keyvals <- ifelse(abs(meth_res_ma$Methdiff1) > 0.05 & meth_res_ma$Qvalue1 < Qvalue_threshold, 'red2',
                   ifelse(abs(meth_res_ma$Methdiff1) > 0.05 & meth_res_ma$Qvalue1 > Qvalue_threshold, 'forestgreen',
                   ifelse(abs(meth_res_ma$Methdiff1) < 0.05 & meth_res_ma$Qvalue1 < Qvalue_threshold, 'royalblue',
                   'grey30')))
        keyvals[is.na(keyvals)] <- 'grey30'
        names(keyvals)[keyvals == 'red2'] <- 'Adjusted P-value and Differential methylation Level'
        names(keyvals)[keyvals == 'forestgreen'] <- 'Differential methylation Level'
        names(keyvals)[keyvals == 'royalblue'] <- 'Adjusted P-value'
        names(keyvals)[keyvals == 'grey30'] <- 'NS'
    	png(paste0(path,"/MA_plot_Transcript_",name[1],"_vs_",name[2],".png"),width=1200,heigh=800)
    	print(EnhancedVolcano(meth_res_ma, x = "Methdiff1", y = "baseMean",lab = NULL, colCustom = keyvals, FCcutoff = 0.05, pCutoff = Qvalue_threshold, xlab = "Differential methylation Level", ylab = "Mean methylation Level", xlim = c(-1,1), ylim = c(0,1), title = 'MA plot of Transcript regions', subtitle = paste0("Cutoff of adjusted p-value at ",Qvalue_threshold,"(",Qvalue_threshold*100,"%) and differential methylation level at 0.05(5%)"), legendPosition = 'right', drawConnectors = TRUE, gridlines.major = FALSE, gridlines.minor = FALSE)+ coord_flip()) #, lab = meth_res_ma$label, selectLab = top_5
    	dev.off()
    
        if(dim(meth_res_significant)[1]>0){
        
            # Piechart
    	    print("Draw Piechart")
        	meth_res_sign_pie <- meth_res_significant
        	meth_res_sign_pie$feature_types <- gsub(meth_res_sign_pie$feature, pattern="_.*", replacement="")
        	meth_res_sign_pie$level <- ifelse(meth_res_sign_pie$Methdiff1 == 0, "No difference", ifelse(meth_res_sign_pie$Methdiff1 < 0, "hypo", "hyper"))
        	
        	count <- as.data.table(table(meth_res_sign_pie[,c("feature_types","level")]))
        	colnames(count) <- c("Transcript regions", "level", "value")
        	
        	count_hyper <- count[which(count$level=="hyper"),] %>% mutate(csum=rev(cumsum(rev(value))),
                        	pos = value/2 + lead(csum, 1),
                        	pos = if_else(is.na(pos), value/2, pos))
        	count_hypo <- count[which(count$level=="hypo"),] %>% mutate(csum=rev(cumsum(rev(value))),
                        	pos = value/2 + lead(csum, 1),
                        	pos = if_else(is.na(pos), value/2, pos))
        	
        	count_hypo$total <- sum(count_hypo$value)
        	count_hyper$total <- sum(count_hyper$value)
        	
        	if(dim(count_hypo)[1] > 0 | dim(count_hyper)[1] > 0){
                if(dim(count_hypo)[1] > 0){
                    hypo <- ggplot(count_hypo, aes(x = " " , y = value, fill = `Transcript regions`)) + 
                                geom_col(width = 1, color=1) +
                                coord_polar(theta = "y", start = 0) + 
                              scale_fill_brewer(palette = "Dark2") +	
                              geom_label_repel(data = count_hypo, aes(y = pos, label = paste0(round((value/total)*100,digit=1),"%")), size = 5, nudge_x = 1, show.legend = F) +
                              theme_void() +
                              ggtitle(paste0("  Hypomethylated regions in ", name[1], " vs ", name[2]),subtitle=paste0("   n = ", sum(count_hypo$value))) +
                              theme(text = element_text(size = 14))
                  } else hypo <- grid::textGrob('No hypomethylated Transcript region.')
                  if(dim(count_hyper)[1] > 0){
                    hyper <- ggplot(count_hyper, aes(x = " " , y = value, fill = `Transcript regions`)) +
                                geom_col(width = 1, color = 1) +
                                coord_polar(theta = "y", start = 0) +
                              scale_fill_brewer(palette = "Dark2") +
                              geom_label_repel(data = count_hyper, aes(y = pos, label = paste0(round((value/total)*100, digit = 1), "%") ), size = 5, nudge_x = 1, show.legend = F) +
                              theme_void() +
                              ggtitle(paste0("  Hypermethylated regions in ", name[1], " vs ", name[2]), subtitle = paste0("   n = ", sum(count_hyper$value))) +
                              theme(text = element_text(size = 14))
                } else hyper <- grid::textGrob('No hypermethylated Transcript region.')

            	png(paste0(path,"/Piechart_significant_Transcript_",name[1],"_vs_",name[2],".png"),width=1200,heigh=500)
            	grid.arrange(hypo, hyper, ncol=2, nrow = 1)
            	dev.off()
            }
        	
        	# Heatmap
    	    print("Draw Heatmap")	
        	meth_res_sign_heat <- meth_res_significant
    	    meth_res_sign_heat$label <- paste(meth_res_sign_heat$transcript_name,meth_res_sign_heat$feature,sep="_")
            hyper <- meth_res_sign_heat[which(meth_res_sign_heat$Methdiff1 > 0),]
            hypo <- meth_res_sign_heat[which(meth_res_sign_heat$Methdiff1 < 0),]
    	
            top_30_hyper <- top_n(hyper,30,Methdiff1)
            top_30_hypo <- top_n(hypo,-30,Methdiff1)
    	
            if(dim(top_30_hypo)[1] > 1 | dim(top_30_hyper)[1] > 1){
        	    if(dim(top_30_hypo)[1] > 1){
        	        heatmap_hypo<-Heatmap_plot(top_30_hypo,c(name[1],name[2]),"hypo")[[4]]
        	    } else heatmap_hypo <- grid::textGrob('No (enough) hypomethylated Transcript region to plot.')
        	    if(dim(count_hyper)[1] > 1){
                    heatmap_hyper<-Heatmap_plot(top_30_hyper,c(name[1],name[2]),"hyper")[[4]]
                } else heatmap_hyper <- grid::textGrob('No (enough) hypermethylated Transcript region to plot.')
                png(paste0(path,"/Heatmap_significant_Transcript_top30_hypo_hyper_",name[1],"_vs_",name[2],".png"),width=1200,heigh=800)
                grid.arrange(grobs=list(heatmap_hypo, heatmap_hyper),ncol=2,nrow=1)
                dev.off()
            }
            
            # DMR par chromosome
    	    print("DMR par chromosome")	
          	meth_res_sign_DRMpc <- meth_res_significant
            meth_res_sign_DRMpc$level <- 0.25
            meth_res_sign_DRMpc$level[which(meth_res_sign_DRMpc$Methdiff1 < 0)] <- -0.25
            meth_res_sign_DRMpc$group <- "hyper"
            meth_res_sign_DRMpc$group[which(meth_res_sign_DRMpc$Methdiff1 < 0)] <- "hypo"
            meth_res_sign_DRMpc$point_shape <- 25
            meth_res_sign_DRMpc$point_shape[which(meth_res_sign_DRMpc$Methdiff1 < 0)] <- 24
              
            size_chr <- fread(chromSize)
            chr <- unique(sort(meth_res_sign_DRMpc$chr))
            p <- list()
            for (c in chr){
              	print(c)
              	assign(x = paste0("size_", c), size_chr[which(size_chr$V1 == c),])
              	p[[paste0(c)]] <- ggplot(meth_res_sign_DRMpc[which(meth_res_sign_DRMpc$chr == c),], aes(x = start , y = level, color = group)) +
              	geom_point(alpha = 0.5, aes(size = 10)) +
              	geom_segment(x = 0, y = 0, xend = as.integer(eval(parse(text = paste0("size_", c)))$V2), yend = 0, size = 10, color = "grey70") +
              	xlim(0, as.integer(max(size_chr$V2)) + 10000000) +
              	ylim(-0.50, 0.30) +
              	theme(axis.text.y = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(),  axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), legend.title = element_blank(), panel.background = element_blank(), legend.position = "none") +
              	geom_label(x = (as.integer(max(size_chr$V2)) + 10000000) , y = 0, label = paste0("Chr ", c), color = "black", size = 12)
            }
            
            png(paste0(path, "/DMR_per_chromosome_significant_Transcript_", name[1], "_vs_", name[2], ".png"), width = 2000, heigh = 3000)
            grid.arrange(grobs = p, ncol = 1, nrow = length(chr))
            dev.off()
        	
        }
        
    }else if(grepl(pattern="/CpG/", x = path)){
        print("CpG analysis")
    	#only on CpG for now, adding CpG Shore+shelf later
    	# change label for CpG_feature_transcrit
    
    	# Correspondance to Gene ID
    	print("Merging with corresponding gene ID")
    	corr_gene <- unique(fread(reference_conversion))
    	colnames(corr_gene) <- c("gene_id", "gene_name", "id", "transcript_name")
    	meth_res <- left_join(x = meth_res, y = corr_gene, by = c("id"))
        meth_res_significant <- left_join(x = meth_res_significant, y = corr_gene, by = c("id"))
        
        # Save tables
        write.table(meth_res,paste0(path, "/meth_res_table_", name[1], "_vs_", name[2], "_all_chr.csv"), sep=",", quote = F, row.names = F)
        if(dim(meth_res_significant)[1]>0) {
            write.table(meth_res_significant, paste0(path, "/meth_res_significant_", name[1], "_vs_", name[2], "_all_chr.csv"), sep=",", quote = F, row.names = F)
        }else print("No differential methylation under Qvalue threshold")
    
        # Enhanced volcanoplot
    	print("Draw Volcano plot")
        meth_res_volca <- meth_res
    	#meth_res_volca$label <- paste(meth_res_volca$transcript_name, meth_res_volca$feature, sep = "_")
    	#top_5 <- unlist(top_n(meth_res_volca,-5,Qvalue1)[,"label"])
    	png(paste0(path, "/Volcano_plot_CpG_", name[1], "_vs_", name[2], ".png"), width=1200, heigh=800)
    	print(EnhancedVolcano(meth_res_volca, x = "Methdiff1", y = "Qvalue1", lab = NULL, FCcutoff = 0.05, pCutoff = Qvalue_threshold, xlab = "Differential methylation Level", ylab = bquote(~Log[10]~ 'Adusted Pvalue'), xlim = c(0,1), title = "Volcano plot of CpG regions", subtitle = paste0("Cutoff of adjusted p-value at ",Qvalue_threshold,"(",Qvalue_threshold*100,"%) and differential methylation level at 0.05(5%)"), legendLabels = c("NS","DiffMeth ","Adjusted P-value","Adjusted P-value and DiffMeth"), legendPosition = 'right', drawConnectors = TRUE, gridlines.major = FALSE, gridlines.minor = FALSE)) #, lab = meth_res_volca$label, selectLab = top_5
    	dev.off()
    
        # MA-plot
    	print("Draw MA-plot")
        meth_res_ma <- meth_res
    	#meth_res_ma$label <- paste(meth_res_ma$transcript_name,meth_res_ma$feature,sep="_")
        meth_res_ma$baseMean <- 1 / (10^(rowMeans(meth_res_ma[,c("Methgroup1_CpGisland","Methgroup2_CpGisland")]))) #EnhancedVolcano will internally transform the y-axis values by -log10(); here, we counteract this by creating a new column in the data for base mean.
        keyvals <- ifelse(abs(meth_res_ma$Methdiff1) > 0.05 & meth_res_ma$Qvalue1 < Qvalue_threshold, 'red2',
                   ifelse(abs(meth_res_ma$Methdiff1) > 0.05 & meth_res_ma$Qvalue1 > Qvalue_threshold, 'forestgreen',
                   ifelse(abs(meth_res_ma$Methdiff1) < 0.05 & meth_res_ma$Qvalue1 < Qvalue_threshold, 'royalblue',
                   'grey30')))
        keyvals[is.na(keyvals)] <- 'grey30'
        names(keyvals)[keyvals == 'red2'] <- 'Adjusted P-value and Differential methylation Level'
        names(keyvals)[keyvals == 'forestgreen'] <- 'Differential methylation Level'
        names(keyvals)[keyvals == 'royalblue'] <- 'Adjusted P-value'
        names(keyvals)[keyvals == 'grey30'] <- 'NS'
    	png(paste0(path, "/MA_plot_CpG_", name[1], "_vs_", name[2], ".png"), width=1200, heigh=800)
    	print(EnhancedVolcano(meth_res_ma, x = "Methdiff1", y = "baseMean",lab = NULL, colCustom = keyvals, FCcutoff = 0.05, pCutoff = Qvalue_threshold, xlab = "Differential methylation Level", ylab = "Mean methylation Level", xlim = c(-1,1), ylim = c(0,1), title = 'MA plot of CpG regions', subtitle = paste0("Cutoff of adjusted p-value at ",Qvalue_threshold,"(",Qvalue_threshold*100,"%) and differential methylation level at 0.05(5%)"), legendPosition = 'right', drawConnectors = TRUE, gridlines.major = FALSE, gridlines.minor = FALSE)+ coord_flip()) #, lab = meth_res_ma$label, selectLab = top_5
    	dev.off()
    	
        if(dim(meth_res_significant)[1]>0){
        
            # Piechart
    	    print("Draw Piechart")
        	meth_res_sign_pie <- meth_res_significant
        	meth_res_sign_pie$feature_types <- gsub(meth_res_sign_pie$feature, pattern="_.*", replacement="")
        	meth_res_sign_pie$level <- ifelse(meth_res_sign_pie$Methdiff1 == 0, "No difference", ifelse(meth_res_sign_pie$Methdiff1 < 0, "hypo", "hyper"))
        	
        	count <- as.data.table(table(meth_res_sign_pie[,c("feature_types","level")]))
        	colnames(count) <- c("CpG regions", "level", "value")
        	
        	count_hyper <- count[which(count$level=="hyper"),] %>% mutate(csum=rev(cumsum(rev(value))),
                        	pos = value/2 + lead(csum, 1),
                        	pos = if_else(is.na(pos), value/2, pos))
        	count_hypo <- count[which(count$level=="hypo"),] %>% mutate(csum=rev(cumsum(rev(value))),
                        	pos = value/2 + lead(csum, 1),
                        	pos = if_else(is.na(pos), value/2, pos))
        	
        	count_hypo$total <- sum(count_hypo$value)
        	count_hyper$total <- sum(count_hyper$value)
        	
        	if(dim(count_hypo)[1] > 0 | dim(count_hyper)[1] > 0){
        	    if(dim(count_hypo)[1] > 0){
                	hypo <- ggplot(count_hypo, aes(x = " " , y = value, fill = `CpG regions`)) + 
                    	        geom_col(width = 1, color=1) +
                    	        coord_polar(theta = "y", start = 0) + 
                        		scale_fill_brewer(palette = "Dark2") +	
                        		geom_label_repel(data = count_hypo, aes(y = pos, label = paste0(round((value/total)*100,digit=1),"%")), size = 5, nudge_x = 1, show.legend = F) +
                        		theme_void() +
                        		ggtitle(paste0("  Hypomethylated regions in ", name[1], " vs ", name[2]),subtitle=paste0("   n = ", sum(count_hypo$value))) +
                        		theme(text = element_text(size = 14))
                } else hypo <- grid::textGrob('No hypomethylated CpG region.')
                if(dim(count_hyper)[1] > 0){
                	hyper <- ggplot(count_hyper, aes(x = " " , y = value, fill = `CpG regions`)) +
                	            geom_col(width = 1, color = 1) +
                	            coord_polar(theta = "y", start = 0) +
                        		scale_fill_brewer(palette = "Dark2") +
                        		geom_label_repel(data = count_hyper, aes(y = pos, label = paste0(round((value/total)*100, digit = 1), "%") ), size = 5, nudge_x = 1, show.legend = F) +
                        		theme_void() +
                        		ggtitle(paste0("  Hypermethylated regions in ", name[1], " vs ", name[2]), subtitle = paste0("   n = ", sum(count_hyper$value))) +
                        		theme(text = element_text(size = 14))
            	} else hyper <- grid::textGrob('No hypermethylated CpG region.')
            	png(paste0(path, "/Piechart_significant_CpG_", name[1], "_vs_", name[2], ".png"), width=1200, heigh=500)
            	grid.arrange(hypo, hyper, ncol=2, nrow = 1)
            	dev.off()
        	}
    
    	    # Heatmap
    	    print("Draw Heatmap")	
        	meth_res_sign_heat <- meth_res_significant
    	    meth_res_sign_heat$label <- paste(meth_res_sign_heat$transcript_name, meth_res_sign_heat$feature, sep="_")
            hyper <- meth_res_sign_heat[which(meth_res_sign_heat$Methdiff1 > 0),]
            hypo <- meth_res_sign_heat[which(meth_res_sign_heat$Methdiff1 < 0),]
            
            top_30_hyper <- top_n(hyper, 30, Methdiff1)
            top_30_hypo <- top_n(hypo, -30, Methdiff1)
            
            if(dim(top_30_hypo)[1] > 1 | dim(top_30_hyper)[1] > 1){
        	    if(dim(top_30_hypo)[1] > 1){
        	        heatmap_hypo <- Heatmap_plot(top_30_hypo, c(name[1], name[2]), "hypo")[[4]]
        	    } else heatmap_hypo <- grid::textGrob('No (enough) hypomethylated CpG region to plot.')
        	    if(dim(count_hyper)[1] > 1){
                    heatmap_hyper <- Heatmap_plot(top_30_hyper, c(name[1], name[2]), "hyper")[[4]]
                } else heatmap_hyper <- grid::textGrob('No (enough) hypermethylated CpG region to plot.')
                png(paste0(path, "/Heatmap_significant_CpG_top30_hypo_hyper_", name[1], "_vs_", name[2], ".png"), width = 1200, heigh = 800)
                grid.arrange(grobs = list(heatmap_hypo, heatmap_hyper), ncol = 2, nrow= 1)
                dev.off()
            }
    	
    	    # DMR par chromosome
    	    print("DMR par chromosome")	
        	meth_res_sign_DRMpc <- meth_res_significant
            meth_res_sign_DRMpc$level <- 0.25
            meth_res_sign_DRMpc$level[which(meth_res_sign_DRMpc$Methdiff1 < 0)] <- -0.25
            meth_res_sign_DRMpc$group <- "hyper"
            meth_res_sign_DRMpc$group[which(meth_res_sign_DRMpc$Methdiff1 < 0)] <- "hypo"
            meth_res_sign_DRMpc$point_shape <- 25
            meth_res_sign_DRMpc$point_shape[which(meth_res_sign_DRMpc$Methdiff1 < 0)] <- 24
            
            size_chr <- fread(chromSize)
            chr <- unique(sort(meth_res_sign_DRMpc$chr))
            p <- list()
            for (c in chr){
            	print(c)
            	assign(x = paste0("size_", c), size_chr[which(size_chr$V1 == c),])
            	p[[paste0(c)]] <- ggplot(meth_res_sign_DRMpc[which(meth_res_sign_DRMpc$chr == c),], aes(x = start , y = level, color = group)) +
            	geom_point(alpha = 0.5, aes(size = 10)) +
            	geom_segment(x = 0, y = 0, xend = as.integer(eval(parse(text = paste0("size_", c)))$V2), yend = 0, size = 10, color = "grey70") +
            	xlim(0, as.integer(max(size_chr$V2)) + 10000000) +
            	ylim(-0.50, 0.30) +
            	theme(axis.text.y = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(),  axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), legend.title = element_blank(), panel.background = element_blank(), legend.position = "none") +
            	geom_label(x = (as.integer(max(size_chr$V2)) + 10000000) , y = 0, label = paste0("Chr ", c), color = "black", size = 12)
            }
            
            png(paste0(path, "/DMR_per_chromosome_significant_CpG_", name[1], "_vs_", name[2], ".png"), width = 2000, heigh = 3000)
            grid.arrange(grobs = p, ncol = 1, nrow = length(chr))
            dev.off()
    
        }
    
    }
}else{
    print("No region to compare, save empty file result...")
    write.table(input_meth_table,paste0(path, "/meth_res_table_", name[1], "_vs_", name[2], ".csv"), sep=",", quote = F, row.names = F)
}

print("Finished")

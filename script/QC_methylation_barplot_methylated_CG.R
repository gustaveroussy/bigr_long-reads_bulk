#script to make QC graphs for methylation analysis

### Library
library(optparse)
library(dplyr)
library(ggplot2)

### Parameters
## parameters list
option_list <- list(
  make_option("--list_input_bed_uncomb_concat", help="List of input bed file uncombined for chrosome and strand, for all samples separated by comma."),
  make_option("--input_bed_motif_cg", help="Input bed file with all CG motifs."),
  make_option("--output_path", help="Output graphic folder."))
parser <- OptionParser(usage="Rscript %prog [options]", description = " ", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
## get parameters
list_input_bed_uncomb_concat <- args$options$list_input_bed_uncomb_concat
input_bed_motif_cg <- args$options$input_bed_motif_cg
output_path <- args$options$output_path
## check parameters
if(is.null(list_input_bed_uncomb_concat)) stop("--list_input_bed_uncomb_concat parameter must be set!")
if(is.null(input_bed_motif_cg)) stop("--input_bed_motif_cg parameter must be set!")
if(is.null(output_path)) stop("--output_path parameter must be set!")

### Set and create output path, set file_list and methylation type
path <- output_path
dir.create(path, showWarnings = FALSE, recursive = TRUE)
file_list <- unlist(stringr::str_split(list_input_bed_uncomb_concat, ","))
meth_type <- basename(dirname(file_list[1]))

### Load and format GC reference
ref <- data.table::fread(input_bed_motif_cg)
cg_motif <- nrow(ref)


concat <- data.frame()
concat_pct <- data.frame()
concat_meth <- data.frame()

for (f in file_list){
	print(f)
	
	#read data
	data <- data.table::fread(f) #V1:chromosome; V2: start; V3: end; V4:percentage of methylated reads at this position; V5: number of methylated reads; V6: number of non-methylated reads.
	
	#number of cg (position) with at least X methylation
	total_cg_min1_nb <- nrow(data)
	total_cg_min1_pct <- paste0(round((total_cg_min1_nb/cg_motif)*100,2), "%")
	data$tot <- rowSums(data[,c(5,6)])
	data_min10 <- data[which(data$tot>=10),]
	total_cg_min10_nb <- nrow(data_min10)
	total_cg_min10_pct <- paste0(round((total_cg_min10_nb/cg_motif)*100,2), "%")
	data_min20 <- data[which(data$tot>=20),]
	total_cg_min20_nb <- nrow(data_min20)
	total_cg_min20_pct <- paste0(round((total_cg_min20_nb/cg_motif)*100,2), "%")
	concat <- rbind(concat,c(cg_motif,total_cg_min1_nb,total_cg_min10_nb,total_cg_min20_nb))
	colnames(concat) <- c("Nb_CG_motifs","CG_1x_min","CG_10x_min","CG_20x_min")
	concat_pct <- rbind(concat_pct,c("100%",total_cg_min1_pct,total_cg_min10_pct,total_cg_min20_pct))
	colnames(concat_pct) <- c("Nb_CG_motifs","CG_1x_min_pc","CG_10x_min_pc","CG_20x_min_pc")

	#number of reads on cg with at least X methylation
	total_cg <- sum(data[,c(5,6)])
    meth_cg <- sum(data[,c(5)])
    total_cg_10x <- sum(data_min10$tot)
    meth_cg_10x <- sum(data_min10$V5)
	total_cg_20x <- sum(data_min20$tot)
    meth_cg_20x <- sum(data_min20$V5)
    numbers_meth <- c(total_cg, meth_cg,total_cg_10x, meth_cg_10x, total_cg_20x, meth_cg_20x)
    concat_meth <- rbind(concat_meth, numbers_meth)
    colnames(concat_meth) <- c("total_cg","meth_cg","total_cg_10x","meth_cg_10x","total_cg_20x","meth_cg_20x")
}

concat$sample <- gsub(basename(file_list), pattern = paste0("_", meth_type, "_uncomb.bed"), replacement="")
concat_pct$sample <- gsub(basename(file_list), pattern = paste0("_", meth_type, "_uncomb.bed"), replacement="")
concat_meth$sample <- gsub(basename(file_list), pattern = paste0("_", meth_type, "_uncomb.bed"), replacement="")

print(concat)
print(concat_pct)

#format concat adn concat_pct for bargraph
melting_concat <- reshape2::melt(concat,id.vars=c("sample"))
melting_concat <- melting_concat[-c(1:(length(file_list)-1)),]#remove duplicated Nb_CG_motifs lines
melting_concat$sample[1] <- "Total CG motifs\nin genome"
melting_concat_pct <- reshape2::melt(concat_pct,id.vars=c("sample"))
melting_concat_pct <- melting_concat_pct[-c(1:(length(file_list)-1)),]#remove duplicated Nb_CG_motifs lines
melting_concat_pct$sample[1] <- "Total CG motifs\nin genome"
melting_concat$value_pct <- melting_concat_pct$value

colo <- RColorBrewer::brewer.pal(n = 4, name = "Dark2")
png(paste0(path,"/barplot_methylated_CG_", meth_type, "_mqc.png"), width = 200 + (200*length(file_list)), height = 450)
ggplot(melting_concat, aes(x = factor(sample, level=c("Total CG motifs\nin genome", unique(sort(melting_concat_pct$sample[-1])))), y = value, fill = variable)) +
    geom_bar(stat = "identity", position = position_dodge(), color = "black") +
    geom_text(data = melting_concat, aes(x = sample, y = value, label = value_pct), position = position_dodge(width = 0.9), vjust = -0.25) +
    scale_fill_manual(values = c(Nb_CG_motifs = colo[1], CG_1x_min = colo[2], CG_10x_min = colo[3], CG_20x_min = colo[4]), breaks = levels(melting_concat$variable)[-1])+
    #scale_fill_brewer(palette = "Dark2") +
    ggplot2::theme_classic() +
    theme(legend.title = element_blank(),text = element_text(size=15)) +
    xlab("") +
    ylab("Number of positions") +
    ggtitle("Reads depth on CG positions")
dev.off()

## TO DO: make a graph with concat_meth


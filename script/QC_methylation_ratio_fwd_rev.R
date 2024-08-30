#script to make QC graphs for methylation analysis

### Library
library(optparse)
library(ggplot2)
library(patchwork)

### Parameters
## parameters list
option_list <- list(
  make_option("--input_bed", help="Input uncombined bed file for one chromosome."),
  make_option("--output_path", help="Output graphic folder."))
parser <- OptionParser(usage="Rscript %prog [options]", description = " ", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
## get parameters
input_bed <- args$options$input_bed
output_path <- args$options$output_path
## check parameters
if(is.null(input_bed)) stop("--input_bed parameter must be set!")
if(is.null(output_path)) stop("--output_path parameter must be set!")

### Set and create output path and input_name (for graphic name)
path <- output_path
dir.create(path, showWarnings = FALSE, recursive = TRUE)
input_name <- sub("_uncomb.bed", "", basename(input_bed))
chr_number <- gsub(".*_chr(.{1,2})_.*","\\1", input_name)

### Loading and formatiing bed file
print("Opening bed file")
bed_file <- data.table::fread(input_bed)
bed_file <- bed_file[,c(1,2,4,5,6)]

### Separate reverse and forward
print("Separating revers and forward records")
bed_file_pos <- bed_file[which(bed_file$V6=="+"),]
bed_file_neg <- bed_file[which(bed_file$V6=="-"),]

### Remove 1 before computing ratio
bed_file_neg$V2 <- bed_file_neg$V2-1

### Reconcatenating neg and pos by position 
pos_and_neg <- merge(bed_file_pos,bed_file_neg,by=c("V1","V2","V4"))
colnames(pos_and_neg) <- c("chr","position","mod","Forward","std.x","Reverse","std.y")
pos_and_neg$ratio <- pos_and_neg$Forward/pos_and_neg$Reverse

### Map number of reads per location
colnames(bed_file) <- c("chr","position","mod","Nreads","Strand")
bed_file$Strand[which(bed_file$Strand=="-")] <- "Reverse"
bed_file$Strand[which(bed_file$Strand=="+")] <- "Forward"

### Graphs
print("Graphique Nombre de read par position")
p1 <- ggplot2::ggplot(bed_file, ggplot2::aes(x = position, y = Nreads, color = Strand)) + 
        ggplot2::geom_point() + 
        ggplot2::theme_classic() +
        ggplot2::scale_color_manual(values=c("#2596be","#e12f8c")) +
        ggplot2::theme(text = element_text(size=17)) +
        ggplot2::ylab("Number of reads") + 
        ggplot2::xlab(paste0("Position on the chromosome ", chr_number)) +
        ggplot2::ggtitle(paste0("Number of reads along chromosome ", chr_number))

print("Graphique Nombre de read par position - approximation GAM")
p2 <- ggplot2::ggplot(bed_file, ggplot2::aes(x = position, y = Nreads, color = Strand)) +
        ggplot2::geom_smooth(span = 0.2) +
        ggplot2::theme_classic() +
        ggplot2::scale_color_manual(values=c("#2596be","#e12f8c")) +
        ggplot2::theme(text = element_text(size=17)) +
        ggplot2::ylab("Smoothed number of reads") + 
        ggplot2::xlab(paste0("Position on the chromosome ", chr_number)) +
        ggplot2::ggtitle(paste0("Number of reads along chromosome ", chr_number," (smoothed version)"))

### Map number ratio between forward and reverse
print("Graphique ratio forward vs reverse")
p3 <- ggplot2::ggplot(pos_and_neg, ggplot2::aes(x = position, y = ratio)) + 
        ggplot2::geom_smooth(span = 0.2, color = "#7475b4") + 
        geom_hline(yintercept = 1, linetype = "dashed", color = "#7475b4", size=2) +
        ggplot2::theme_classic() +
        ggplot2::theme(text = element_text(size=17)) +
        ggplot2::ylab("Ratio (number reads forward / number reads reverse)") + 
        ggplot2::xlab(paste0("Position on the chromosome ", chr_number)) +
        ggplot2::ggtitle(paste0("Ratio of the number of reads forward/reverse along chromosome ", chr_number))


png(paste0(path, "/per_pos_per_strand_gam_ratio_", input_name, "_mqc.png"), width = 2000, height = 1800)
print(p1 / p2 / p3)
dev.off()



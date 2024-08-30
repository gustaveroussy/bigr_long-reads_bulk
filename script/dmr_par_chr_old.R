library(data.table)
library(GeneDMRs)
library(dplyr)
library(pheatmap)
library(grid)
library(EnhancedVolcano)
library(ggrepel)
library(gridExtra)
library(ggplot2)



donnee<-fread("regiongeneall_significant_5percent_NOX4_vs_Ctrol.csv")



size_chr<-fread("/mnt/beegfs/userdata/y_mesloub/ressources/references/grch38/hg38_centroLoc_chromSize.txt")

donnee$level<-0.25
donnee$level[which(donnee$Methdiff1 < 0)]<--0.25
donnee$group<-"hyper"
donnee$group[which(donnee$Methdiff1 < 0)]<-"hypo"
donnee$point_shape<-25
donnee$point_shape[which(donnee$Methdiff1 < 0)]<-24


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
##############################################################################


p=list()

chr<-c(1:22,"X","Y")
for (c in chr){
	print(c)
	assign(x = paste0("size_", c), unlist(size_chr[which(size_chr$V1 == paste0("chr", c)),]))
	p[[c]] <- ggplot(donnee[which(donnee$chr == c),], aes(x = start , y = level, color = group)) +
	geom_point(alpha = 0.5, aes(size = 10)) +
	geom_segment(x = 0, y = 0, xend = as.integer(eval(parse(text = paste0("size_", c)))[4]), yend = 0, size = 10, color = "grey70") +
	xlim(0, as.numeric(max(size_chr$V4)) + 10000000) +
	ylim(-0.50, 0.30) +
	theme(axis.text.y = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(),  axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), legend.title = element_blank(), panel.background = element_blank(), legend.position = "none") +
	geom_label(x = (as.numeric(max(size_chr$V4)) + 10000000) , y = 0, label = paste0("Chr ", c), color = "black", size = 12)
}

png("DMR_on_chromomosome2.png", width = 2000, heigh = 3000)
grid.arrange(grobs = p, ncol = 1, nrow = length(chr))
dev.off()













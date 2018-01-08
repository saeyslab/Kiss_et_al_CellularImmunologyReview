# This should be adapted to the appropriate path on your computer -------------
# This folder should contain a folder "experiment_101718_files",              -
# as downloaded from cytobank                                                    -
# and two files parameterTypes.xlsx and sampleIDs.xlsx ------------------------

setwd("FinalExperiment/")

# Load libraries --------------------------------------------------------------

library(tidyverse)
library(xlsx)
library(flowCore)
library(flowDensity)
library(FlowSOM)
library(Rtsne)
library(pheatmap)

cytoftrans <- arcsinhTransform(transformationId="cytofTransform",
                               a=0,b=(1/5),c=0)

seed <- 1

# Directory settings ----------------------------------------------------------

fcs_dir <- "experiment_101718_files"
res_dir <- "experiment_101718_analysis"
dir.create(res_dir)
dir.create(file.path(res_dir,"agg"))

# Data ------------------------------------------------------------------------

ann <- read.table("experiment_101718_files/experiment_101718_annotations.tsv",
                  sep = "\t",
                  header = TRUE) %>% 
  dplyr::select(c(FCS.Filename,
                  Plate, 
                  Plate.Row, 
                  Plate.Column, 
                  FCS.File.Category)) %>% 
  dplyr::mutate("Sample.ID" = as.numeric(gsub("rcc", "", 
                                              gsub("[b]*.fcs", "", 
                                                   FCS.Filename)))) %>% 
  dplyr::arrange(Sample.ID) 
ann2 <- read.xlsx("sampleIDs.xlsx", "Sheet1", endRow = 79)
ann <- full_join(ann, ann2, by = "Sample.ID") %>% 
  dplyr::filter(!(Sample.ID %in% c(61, 17, 21:25)))
str(ann)

markers <- read.xlsx("parameterTypes.xlsx", 
                     "Sheet1", 
                     stringsAsFactors = FALSE) %>% 
  dplyr::select(-column) %>% 
  dplyr::mutate(Type = factor(Type))

markersToUse <- markers$name[markers$ToUse]
# Param not found: CD166
# PD-L2 = CD273, PD-L1 = CD274

# Generate an aggregate file --------------------------------------------------

set.seed(seed)
agg <- AggregateFlowFrames(file.path(fcs_dir, ann$FCS.Filename),
                           5000000,
                           writeOutput = TRUE,
                           outputFile = file.path(res_dir,
                                                  paste0("aggregate",
                                                         seed,
                                                         ".fcs")))

agg <- read.FCS(file.path(res_dir,paste0("aggregate",
                                         seed,
                                         ".fcs")))
agg <- transform(agg, transformList(colnames(agg)[12:52], cytoftrans))

 # QC plots of the Aggregate file ----------------------------------------------

png(file.path(res_dir,"agg_cellCounts.png"))
plot(table(exprs(agg)[,"File"]), ylab = "Number of cells in aggregate")
dev.off()

for(marker in markersToUse){
  prettyName <- markers %>% dplyr::filter(name == marker) %>%  pull(desc)
  png(file.path(res_dir,"agg",paste0(marker,".png")),
      height = 1000, width = 500)
  par(mar=c(5,8,4,2))
  plot(exprs(agg)[,c(marker,"File_scattered")], pch=".", col="#00000088",
       main = prettyName, 
       xlim=c(0,8),
       yaxt="n",ylab="")
  axis(2, 
       at = seq_along(ann$Sample.ID),
       labels = ann$FCS.Filename,
       las = 2)
  dev.off()
}

# FlowDensity cut-off for CD68 ------------------------------------------------

data_dir <- file.path(res_dir, "selected")
dir.create(data_dir)

cd68 <- "Nd143Di"
threshold <- flowDensity::deGate(agg, cd68, graphs = TRUE)

png(file.path(res_dir,"CD68_cutoff.png"),
    width = 3000, height = 2000)
  layout(matrix(1:72, ncol = 6, byrow = TRUE))
  par(mar = c(2,1,3,1))
  for (file in ann$FCS.Filename) {
    ff <- read.FCS(file.path(fcs_dir, file))
    ff <- transform(ff, transformList(colnames(agg)[12:52], cytoftrans))
    plot(density(exprs(ff[,cd68])),
         xlim = c(0, 8),
         main = file, ylab = "", xlab = "", yaxt = "n")
    abline(v = threshold)
    write.FCS(ff[ff[, cd68] > threshold,], file.path(data_dir, file))
  }
dev.off()

# New aggregate with only CD68 positive cells ---------------------------------

set.seed(seed)
agg <- AggregateFlowFrames(file.path(data_dir, ann$FCS.Filename),
                           5000000,
                           writeOutput = TRUE,
                           outputFile = file.path(res_dir,
                                                  paste0("aggregate",
                                                         seed,
                                                         "_CD68.fcs")))

agg <- read.FCS(file.path(res_dir,paste0("aggregate",
                                         seed,
                                         "_CD68.fcs")))

png(file.path(res_dir,"agg_cellCounts_CD68.png"))
  plot(table(exprs(agg)[,"File"]), 
       ylab = "Number of cells in CD68+ aggregate")
dev.off()

# Tsne ------------------------------------------------------------------------

set.seed(seed)
selected <- sample(seq_len(nrow(agg)), 100000)
tsne_sub <- exprs(agg)[selected, markersToUse]
rtsne <- Rtsne(tsne_sub, perplexity = 30)
saveRDS(rtsne, file = file.path(res_dir, "tsne_CD68.rds"))

# Tsne plots ------------------------------------------------------------------

rtsne <- readRDS(file.path(res_dir, "tsne_CD68.rds"))
set.seed(seed)
selected <- sample(seq_len(nrow(agg)), 100000)

colors <- c("Adhesion molecule" = "purple4",
            "Chemokine receptor" = "darkslateblue",
            "Ectoenzyme" = "darkorange2",
            "Fc & complement receptor" = "darkred",
            "Immunomodulatory molecule" = "olivedrab4",
            "Scavenger receptor" = "violetred4",
            "TLR/cytokine receptor" = "darkgreen")

markerOrder <- c("CD273", "CD274", "CD40", "Slamf7", "CD86", "HLA-ABC", 
                 "HLA-DR", "CD68", "CD169", "CD206", "CD54", "CD82",
                 "CD81", "CD163", "CD204", "CD36", "CD13", "CD38",
                 "CD16", "CD32", "CD64", "CD11b", "CD88", "CD119",
                 "CD123", "CD14", "CD71", "CD304", "CXCR4")

png(file.path(res_dir,"tSNE.png"), width = 2500, height = 1500)
layout(matrix(sort(c(seq(1, 60, by=2), rep(seq(2, 60, by=2), each=2))), 
              nrow = 5, byrow = TRUE))
par(mar = c(1,0,4,0), cex = 1.5)
for (m in markerOrder) {
  i <- which(markers$desc == m)
  plot.new()
  FlowSOM:::legendContinuous(
    colorRampPalette(c("lightgrey",
                       colors[as.character(markers$Type[i])]))(100),
    c(0,8))
  plot(rtsne$Y,
       col=colorRampPalette(c("lightgrey",
                              colors[as.character(markers$Type[i])]))(100)[
         as.numeric(cut(exprs(agg[selected,markers$name[i]]),
                        breaks = seq(0,8,length.out = 100)))],
       main=markers$desc[i],bty="n",axes=F,xlab="",ylab="",pch=19)
}
layout(1)
dev.off()

# original values of the scale:
sinh(c(0,1.52, 3.12, 4.74, 6.32, 7.92))/(1/5)

# FlowSOM ---------------------------------------------------------------------

fsom <- FlowSOM(agg,
                scale = FALSE,
                xdim = 7, ydim = 7,
                colsToUse = markersToUse,
                nClus = 15,
                seed = seed)
saveRDS(fsom, file = file.path(res_dir, "fsom_CD68.rds"))

# FlowSOM Plot ----------------------------------------------------------------

fsom <- readRDS(file.path(res_dir, "fsom_CD68.rds"))

colors <- c("Adhesion molecule" = "purple1",
            "Chemokine receptor" = "slateblue",
            "Ectoenzyme" = "orange",
            "Fc & complement receptor" = "red",
            "Immunomodulatory molecule" = "olivedrab1",
            "Scavenger receptor" = "violetred1",
            "TLR/cytokine receptor" = "green")

markers_fsom_plot <- c("CD206", "CD204",
                   "CD163", "HLA-DR")
markers_fsom_plot <- sapply(markers_fsom_plot,
                            function(x){
                              markers$name[which(markers$desc == x)]
                            })
markers_fsom_colors <- c("purple1", "deepskyblue","violetred1", "olivedrab1") 
# markers_fsom_colors <- sapply(markers_fsom_plot,
#                               function(x){
#                                 colors[as.character(markers$Type[
#                                   which(markers$name == x)])]
#                               })

pdf(file.path(res_dir,"FlowSOM.pdf"))
  metacluster_colors <-  c("#17c1a5", "#a5cde2", "#1f78b3",
                         "#b2de8b", "#359e2c", "#fbb3ad",
                         "#e2191b", "#fdbd6f", "#ff7e01",
                         "#c9b1d5", "#6b3d99", "#fffc98",
                         "#b05927", "#e62a89", "#f681bf")
  PlotStars(UpdateNodeSize(fsom$FlowSOM, reset = TRUE, maxNodeSize = 10),
            markers = markers_fsom_plot,
            colorPalette = colorRampPalette(markers_fsom_colors),
            backgroundValues = fsom$metaclustering,
            backgroundColor = metacluster_colors,
            main = "FlowSOM analysis using 29 markers")
dev.off()

png(file.path(res_dir, "tSNE_FlowSOM_mapping.png"),
    width = 1800, height = 600)
  layout(matrix(1:2,nrow=1))
  plot(rtsne$Y,
       col= colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(100)[
         fsom$FlowSOM$map$mapping[selected,1]],
       main="FlowSOM Clusters", bty="n", axes=F, xlab="", ylab="", pch=19)
 
  metacluster_colors <-  c("#17c1a5", "#a5cde2", "#1f78b3",
                           "#b2de8b", "#359e2c", "#fbb3ad",
                           "#e2191b", "#fdbd6f", "#ff7e01",
                           "#c9b1d5", "#6b3d99", "#fffc98",
                           "#b05927", "#e62a89", "#f681bf")
  plot(rtsne$Y,
       col=metacluster_colors[
         fsom$metaclustering[fsom$FlowSOM$map$mapping[selected,1]]],
       main="FlowSOM Meta clusters", bty="n", axes=F, xlab="", ylab="", pch=19)
       
dev.off()

pdf(file.path(res_dir, "FlowSOM_metaclusters_heatmap.pdf"))
  meta_medians <- apply(fsom$FlowSOM$map$medianValues[, markersToUse],
                        2, function(x){tapply(x, fsom$metaclustering, median)})
  colnames(meta_medians) <- sapply(colnames(meta_medians),
                                   function(x){
                                     markers$desc[which(markers$name == x)]
                                   })
  
  pheatmap::pheatmap(meta_medians,
                     color = colorRampPalette(
                       RColorBrewer::brewer.pal(9,"YlOrRd"))(100))
dev.off()
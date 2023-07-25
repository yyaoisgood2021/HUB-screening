library(dplyr)
library(Seurat)
library(patchwork)
library(scater)
# library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(RColorBrewer)
library(Signac)
library(ggplot2)
library(stringr)


args <- commandArgs(trailingOnly = TRUE)



data_sv_folder <- args[1]
aggr_out_folder <- args[2]




# Load the dataset from the cellranger aggr's output folder `outs/count/filtered_feature_bc_matrix`
data <- Read10X(data.dir = aggr_out_folder)       


# GEMs without cells are filtered out by min.features; and genes that are not frequently detected are removed by min.cells
SeuObj <- CreateSeuratObject(counts = data, min.cells = 200, min.features = 200)

SeuObj[["percent.mt"]] <- PercentageFeatureSet(SeuObj, pattern = "^MT-") # mitochondrial perc

head(SeuObj@meta.data, 5)
#                       orig.ident nCount_RNA nFeature_RNA percent.mt
# AAACCCAAGCCGTCGT-1 SeuratProject       5915         2350 2.60355030
# AAACCCAAGTCTAGAA-1 SeuratProject      23798         6415 2.48340197
# AAACCCACAAGTGCTT-1 SeuratProject      22056         5830 0.05894088
# AAACCCACACTAGTAC-1 SeuratProject      20520         5962 2.47076023
# AAACCCACATCTATCT-1 SeuratProject      27141         6299 2.60859954


# barcode ending with -1 indicates ess and -2 indicates 
# count how many cells are in each library and replace orig.ident based on essentiality 
new_id_list <- c(rep(c("ess"), 11052), rep(c("noness"), 7755))
SeuObj@meta.data$orig.ident <- new_id_list

tail(head(SeuObj@meta.data, 11055),5)
#                    orig.ident nCount_RNA nFeature_RNA percent.mt
# TTTGTTGTCGTCTAAG-1        ess       5615         2168   3.259127
# TTTGTTGTCTTTGCAT-1        ess       3020         1283   6.655629
# AAACCCAAGATCACTC-2     noness      18411         5325   2.368149
# AAACCCACAATAGGAT-2     noness      11659         3408   2.787546
# AAACCCACACCTGCAG-2     noness      31621         6901   2.049271


pdf(file.path(data_sv_folder, 'vlnplot.raw.pdf'), height=4, width=10)
VlnPlot(SeuObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
	group.by = 'orig.ident', ncol = 3)
dev.off()




# filter out the low quality single cells
SeuObj_remain <- subset(SeuObj, subset = (percent.mt < 5) & (nCount_RNA > 5000) & (nCount_RNA < 50000)) 
c(table(SeuObj_remain$"orig.ident"))
  #  ess noness 
  # 9333   5905 

pdf(file.path(data_sv_folder, 'vlnplot.low-qual-removed.pdf'), height=4, width=10)
VlnPlot(SeuObj_remain, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
	group.by = 'orig.ident', ncol = 3)
dev.off()

saveRDS(SeuObj_remain, file.path(vlnplot.low-qual-removed.pdf, 'hubs.high_quality.combined.s1.rds'))











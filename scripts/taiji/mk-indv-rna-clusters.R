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

seurat_obj_path <- args[1]
seurat_sv_folder_2 <- args[2]
Npcs <- args[3] %>% as.integer
Npcs_cls <- args[4] %>% as.integer



cds <- readRDS(seurat_obj_path)


# remove MT genes before analysis
avail_genes <- rownames(cds)
MT_genes <- avail_genes[grepl('^MT-', avail_genes)]
remain_genes <- avail_genes [! (avail_genes %in% MT_genes)]
sce <- cds[remain_genes, ]


# consider ess and noness datasets individually

# extact ess cells and normalized data within
sce_ess <- sce[, sce$orig.ident=='ess']
sce_ess <- NormalizeData(sce_ess, normalization.method = "LogNormalize", scale.factor = 10000)
sce_ess <- FindVariableFeatures(sce_ess, selection.method = "vst", nfeatures = 2000)

pdf(file.path(seurat_sv_folder_2,'ess-variable-genes.pdf'), height=4, width=12)
plot1 <- VariableFeaturePlot(sce_ess)
top10 <- head(VariableFeatures(sce_ess), 10)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()

sce_cur <- ScaleData(sce_ess, features = remain_genes)
sce_cur <- RunPCA(sce_cur, npcs = Npcs) # npc=50 is by default
sce_cur <- RunUMAP(sce_cur, reduction = "pca", dims = 1:Npcs_cls)
sce_cur <- FindNeighbors(sce_cur, reduction = "pca", dims = 1:Npcs_cls)
sce_cur <- FindClusters(sce_cur, resolution = c(1,2,2.5,3,1.5,0.75,0.5))

sv_file_name <- paste0('ess_rna.with_cluster_info.npcs',Npcs,'_pc',Npcs_cls,'.s2.rds')
saveRDS(sce_cur, file.path(seurat_sv_folder_2,sv_file_name))

sk_file_name <- paste0('ess_rna.cls_info.npcs',Npcs,'_pc',Npcs_cls,'.txt')
sink(file.path(seurat_sv_folder_2,sk_file_name))
for (res in c (1,2,2.5,3,1.5,0.75,0.5)) {
  print(res)
  stat_df <- data.frame(matrix(ncol = 2, nrow = 0)) 
  clo_nm <- paste0('RNA_snn_res.',res)
  colnames(stat_df) <- c('cluster_id', 'ess_num')
  for (i in unique(sce_cur@meta.data[[clo_nm]]) ){
    s <- dim(sce_cur@meta.data[sce_cur@meta.data[[clo_nm]]==i, ]) [1]
    stat_df[nrow(stat_df)+1,] <- c(i, s)
  }
  stat_df <- stat_df[order(as.numeric(stat_df$cluster_id)),]
  rownames(stat_df) <- 1:nrow(stat_df)
  print(stat_df)
}
sink()


# extact noness cells and normalized data within
sce_noness <- sce[, sce$orig.ident=='noness']
sce_noness <- NormalizeData(sce_noness, normalization.method = "LogNormalize", scale.factor = 10000)
sce_noness <- FindVariableFeatures(sce_noness, selection.method = "vst", nfeatures = 2000)

pdf(file.path(seurat_sv_folder_2,'noness-variable-genes.pdf'), height=4, width=12)
plot1 <- VariableFeaturePlot(sce_noness)
top10 <- head(VariableFeatures(sce_noness), 10)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()

sce_cur <- ScaleData(sce_noness, features = remain_genes)
sce_cur <- RunPCA(sce_cur, npcs = Npcs) # npc=50 is by default
sce_cur <- RunUMAP(sce_cur, reduction = "pca", dims = 1:Npcs_cls)
sce_cur <- FindNeighbors(sce_cur, reduction = "pca", dims = 1:Npcs_cls)
sce_cur <- FindClusters(sce_cur, resolution = c(1,2,2.5,3,1.5,0.75,0.5))

sv_file_name <- paste0('noness_rna.with_cluster_info.npcs',Npcs,'_pc',Npcs_cls,'.s2.rds')
saveRDS(sce_cur, file.path(seurat_sv_folder_2,sv_file_name))

sk_file_name <- paste0('noness_rna.cls_info.npcs',Npcs,'_pc',Npcs_cls,'.txt')
sink(file.path(seurat_sv_folder_2,sk_file_name))
for (res in c (1,2,2.5,3,1.5,0.75,0.5)) {
  print(res)
  stat_df <- data.frame(matrix(ncol = 2, nrow = 0)) 
  clo_nm <- paste0('RNA_snn_res.',res)
  colnames(stat_df) <- c('cluster_id', 'noness_num')
  for (i in unique(sce_cur@meta.data[[clo_nm]]) ){
    s <- dim(sce_cur@meta.data[sce_cur@meta.data[[clo_nm]]==i, ]) [1]
    stat_df[nrow(stat_df)+1,] <- c(i, s)
  }
  stat_df <- stat_df[order(as.numeric(stat_df$cluster_id)),]
  rownames(stat_df) <- 1:nrow(stat_df)
  print(stat_df)
}
sink()




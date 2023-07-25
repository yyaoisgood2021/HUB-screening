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
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(stringr)


##################################################################################################
# this step requires npcs and pc param, calc gene activities for atac_ds according to the specified rna_ds
analyze_gene_activity_atac <- function(atac_ds, rna_ds, sample_id, ds_save_folder, ds_save_folder_2,
                                       gene_set = 'variable',
                                       rna_expr_id = 'npcs30_pc30'
) {
  print('atac data normalization and dim reduction ...')  
  atac_ds <- RunTFIDF(atac_ds)
  atac_ds <- FindTopFeatures(atac_ds, min.cutoff = "q0")
  atac_ds <- RunSVD(atac_ds)
  atac_ds <- RunUMAP(atac_ds, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
  
  if (gene_set == 'all') {
    print('use all the genes')
    features_to_use <- rownames(rna_ds)
  } else if (gene_set == 'variable') {
    print('use top 2000 varaible genes')
    features_to_use <- VariableFeatures(rna_ds)
  } else {
    stop("use 'all' or 'variable' gene set")
  }
  
  print('calc gene activities ...')
  gene.activities <- GeneActivity(atac_ds, features = features_to_use)
  atac_ds[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
  
  # normalize gene activities
  DefaultAssay(atac_ds) <- "ACTIVITY"
  atac_ds <- NormalizeData(atac_ds)
  atac_ds <- ScaleData(atac_ds, features = rownames(atac_ds))
  
  print('saving scaled and anchor identified atac ds to ...')
  sv_path = file.path(ds_save_folder_2,
                      paste0('atac-just-before-anchor-transfer.',
                             sample_id,'.',
                             rna_expr_id,
                             '.s2.rds'))
  print(sv_path)
  saveRDS(atac_ds,sv_path)
  return (atac_ds)
}

##################################################################################################
# this step predict atac cell cls based on snn_res param
cls_pred_atac <- function(atac_ds, rna_ds, sample_id, ds_save_folder, plot_save_folder,
                          gene_set = 'variable', rna_expr_id = 'npcs30_pc30',
                          snn_res = 3, save_plots = F, save_ds = F) {
  cls_col_name <- paste0('RNA_snn_res.',snn_res)
  print(paste0('using ',cls_col_name))
  
  if (save_plots) {
    pdf(file.path(plot_save_folder, 
                  paste0('dimplot.',
                         sample_id,'.',
                         rna_expr_id,
                         '.before-anchor.pdf')), 
        height=4, width=9
    )
    p1 <- DimPlot(rna_ds, group.by = cls_col_name, label = TRUE) + NoLegend() + ggtitle("RNA")
    p2 <- DimPlot(atac_ds, group.by = "orig.ident", label = FALSE) + NoLegend() + ggtitle("ATAC")
    print(p1 + p2)
    dev.off()
  }
  if (gene_set == 'all') {
    print('use all the genes')
    features_to_use <- rownames(rna_ds)
  } else if (gene_set == 'variable') {
    print('use top 2000 varaible genes')
    features_to_use <- VariableFeatures(rna_ds)
  } else {
    stop("use 'all' or 'variable' gene set")
  }
  print('running anchor transfer ...')  
  transfer.anchors <- FindTransferAnchors(
    reference = rna_ds, 
    query = atac_ds, 
    features = features_to_use,
    reference.assay = "RNA", 
    query.assay = "ACTIVITY", 
    reduction = "cca")
  
  celltype.predictions <- TransferData(
    anchorset = transfer.anchors, 
    refdata = rna_ds@meta.data[, cls_col_name],
    weight.reduction = atac_ds[["lsi"]], 
    dims = 2:30)
  
  atac_ds <- AddMetaData(atac_ds, metadata = celltype.predictions)
  
  if (save_plots) {
    pdf(file.path(plot_save_folder, 
                  paste0('dimplot.',
                         sample_id,'.',
                         rna_expr_id,
                         '.after-anchor.pdf')), 
        height=4, width=9
    )
    
    p1 <- DimPlot(rna_ds, group.by = cls_col_name, label = TRUE) + NoLegend() + ggtitle("RNA")
    p2 <- DimPlot(atac_ds, group.by = "predicted.id", label = TRUE) + NoLegend() + ggtitle("Predicted annotation")
    print(p1 | p2)
    dev.off()
  }
  
  if (save_ds) {
    save_ds_path <- file.path(
      ds_save_folder,
      paste0('atac-after-anchor-transfer.',
             sample_id,'.',
             rna_expr_id,'.',
             snn_res,
             '.s2.rds'))
    saveRDS(atac_ds, save_ds_path)
  }
  return (atac_ds)
}


cls_pred_atac_all <- function(
    ess_atac_ds, noness_atac_ds,
    ess_rna_ds, noness_rna_ds, 
    ds_save_folder, plot_save_folder,
    rna_expr_id, res, gene_set = 'variable',
    save_plots = F, save_ds = F) {
  # N2 is ess, N7 is noness, this function will call cls_pred_atac twice to make atac_ds and pred cls

  ess_atac_ds_cls <- cls_pred_atac(
    atac_ds=ess_atac_ds, rna_ds=ess_rna_ds, sample_id='ess', 
    ds_save_folder=ds_save_folder, plot_save_folder=plot_save_folder,
    gene_set=gene_set, rna_expr_id=rna_expr_id,snn_res=res, 
    save_plots=save_plots, save_ds=save_ds
  )
  noness_atac_ds_cls <- cls_pred_atac(
    atac_ds=noness_atac_ds, rna_ds=noness_rna_ds, sample_id='noness', 
    ds_save_folder=ds_save_folder, plot_save_folder=plot_save_folder,
    gene_set=gene_set, rna_expr_id=rna_expr_id,snn_res=res, 
    save_plots=save_plots, save_ds=save_ds
  )
  result <- list(
    'ess'=ess_atac_ds_cls,
    'noness'=noness_atac_ds_cls
  )
  return (result)
}


##################################################################################################
# load scRNA data with the given npcs and pc param (already has cluster information)
args <- commandArgs(trailingOnly = TRUE)

ess_clustered_rna_seurat_obj <- args[1]
noness_clustered_rna_seurat_obj <- args[2]
atac_folder <- args[3]

Npcs <- args[4] %>% as.integer
Npcs_cls <- args[5] %>% as.integer
res <- args[6] %>% as.numeric

sv_folder <- args[7]



print(paste0(
  'using Npcs = ',
  Npcs, 
  ', Npcs_cls = ',
  Npcs_cls,
  ', res = ',
  res,
  '...'
))

sce_ess <- readRDS(ess_clustered_rna_seurat_obj)
sce_noness <- readRDS(noness_clustered_rna_seurat_obj)

##################################################################################################
# create annotation files
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
genome(annotations) <- "hg38"
new_names <- c()
for (i in seqlevels(annotations)) {
  print(paste('chr', i, sep=''))
  new_names <- append(new_names, c(paste('chr', i, sep='')))
}
new_names <- as.vector(new_names)
annotations <- renameSeqlevels(annotations, new_names)

##################################################################################################
N2.atac_remain <- readRDS(file.path(atac_folder, 'atac.ess.s1.rds'))
N7.atac_remain <- readRDS(file.path(atac_folder, 'atac.noness.s1.rds'))

##################################################################################################
# gene activity computation
rna_expr_id = paste0('npcs',Npcs,'_pc',Npcs_cls)

N2.atac_remain <- analyze_gene_activity_atac(
  atac_ds = N2.atac_remain, 
  rna_ds = sce_ess,
  sample_id = 'ess', 
  ds_save_folder = '', 
  ds_save_folder_2 = sv_folder,
  gene_set = 'variable',
  rna_expr_id = rna_expr_id)

N7.atac_remain <- analyze_gene_activity_atac(
  atac_ds = N7.atac_remain, 
  rna_ds = sce_noness,
  sample_id = 'N7', 
  ds_save_folder = '', 
  ds_save_folder_2 = sv_folder,
  gene_set = 'variable',
  rna_expr_id = rna_expr_id)

##################################################################################################
# transfer anchor and cls prediction
result <- cls_pred_atac_all(
  ess_atac_ds = N2.atac_remain, 
  noness_atac_ds = N7.atac_remain,
  ess_rna_ds = sce_ess,
  noness_rna_ds = sce_noness,
  ds_save_folder = sv_folder, 
  plot_save_folder = sv_folder,
  rna_expr_id = rna_expr_id, 
  res = res, 
  save_plots = F, save_ds = T)

ess_atac_ds_cls <- result$'ess'
noness_atac_ds_cls <- result$'noness'

print(ess_atac_ds_cls)
print(noness_atac_ds_cls)

stas_df <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(stas_df) <- c('id', 'cls', 'rna-number', 'atac-number')

ess_rna_ds <- sce_ess@meta.data
noness_rna_ds <- sce_noness@meta.data
ess_atac_ds_cls <- ess_atac_ds_cls@meta.data
noness_atac_ds_cls <- noness_atac_ds_cls@meta.data

colname_for_statdf <- paste0('RNA_snn_res.',res)
print('add cls-info to stat-df ...')

# ess
cls_list <- ess_rna_ds[[colname_for_statdf]]

for (cls_id in unique(cls_list)) {
  rna_num <- nrow(ess_rna_ds[ess_rna_ds[[colname_for_statdf]]==cls_id,])
  atac_num <- nrow(ess_atac_ds_cls[ess_atac_ds_cls$'predicted.id'==cls_id,])
  stas_df[nrow(stas_df)+1,] <- c('ess', cls_id, rna_num, atac_num)
}

# noness
cls_list <- noness_rna_ds[[colname_for_statdf]]
for (cls_id in unique(cls_list)) {
  rna_num <- nrow(noness_rna_ds[noness_rna_ds[[colname_for_statdf]]==cls_id,])
  atac_num <- nrow(noness_atac_ds_cls[noness_atac_ds_cls$'predicted.id'==cls_id,])
  stas_df[nrow(stas_df)+1,] <- c('noness', cls_id, rna_num, atac_num)
}

# write stas_df
sv_path_stas_df <- paste0('stat_df.',rna_expr_id,'.',res,'.txt')
write.table(stas_df,
            file.path(sv_folder,
                      sv_path_stas_df),
            sep='\t', row.names=F, quote=F)
print(stas_df)
print('all done!')

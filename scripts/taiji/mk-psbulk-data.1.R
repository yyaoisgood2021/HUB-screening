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

# this file loads the separate ess and noness rna object (with cls information)
# also loads separate ess and noness filtered and anchor transfered and cls predicted atac objects
# this script performs anchor transfer to predict atac cell cls

args <- commandArgs(trailingOnly = TRUE)

rna_rds_save_path <- args[1]
atac_rds_save_path <- args[2]

npcs <- args[3]
pc_for_cls <- args[4]
res <- args[5]

sv_folder_base <- args[6]

param_info <- paste0(
  'npcs',npcs,
  '_pc',pc_for_cls,
  '.',res
)

cls_col_name <- paste0('RNA_snn_res.',res)

# criteria for filtering out low confidence pseudobulk clusters
smallest_rna_community <- 100
smallest_atac_community <- 20
min_cell_for_atac_peaks <- 1

data_save_folder <- file.path(
  sv_folder_base,
  param_info)


##################################################################################################
# make folder to save data (rna bulk + atac peak.bed)
if (! file.exists(data_save_folder)){
  dir.create(data_save_folder)
} 

##################################################################################################
# already created atac object 'atac-after-anchor-transfer.{N7 or N2}.{npcs30_pc30.1.5}.s2.rds'
# load ess and noness atac object (with cls info)

ess_atac_ds <- readRDS(file.path(
  atac_rds_save_path,
  paste0(
    'atac-after-anchor-transfer.N2.',
    param_info,
    '.s2.rds'
    )))

noness_atac_ds <- readRDS(file.path(
  atac_rds_save_path,
  paste0(
    'atac-after-anchor-transfer.N7.',
    param_info,
    '.s2.rds'
  )))

DefaultAssay(ess_atac_ds) <- "peaks"
DefaultAssay(noness_atac_ds) <- "peaks"

##################################################################################################
# load rna object (separately)
# '{noness or ess}_rna.with_cluster_info.{npcs10_pc10}.s1.rds'
ess_rna_RDS_path <- paste0('ess_rna.with_cluster_info.npcs',npcs,
                       '_pc',pc_for_cls,'.s1.rds')
sce_ess <- readRDS(file.path(
  rna_rds_save_path,
  ess_rna_RDS_path))

noness_rna_RDS_path <- paste0('noness_rna.with_cluster_info.npcs',npcs,
                           '_pc',pc_for_cls,'.s1.rds')
sce_noness <- readRDS(file.path(
  rna_rds_save_path,
  noness_rna_RDS_path))


##################################################################################################
# load cls-info data 'stat_df.{npcs30_pc10.0.5}.txt'
cls_info <- read.table(
  file = file.path(
    atac_rds_save_path,
    paste0(
      'stat_df.',
      param_info,
      '.txt')
  ), header=T
)

cls_info_remained <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(cls_info_remained) <- c('cluster_id',
                                 'RNA.community_size',
                                 'ATAC.community_size')


##################################################################################################
# function for merging psbulk rna and atac data
extract_psbulk_data_for_cls_i <- function(i, rna_ds, atac_ds, cls_info_df_this_type, ds_type,
                                          smallest_rna_community,smallest_atac_community,
                                          ess_enrich_high,ess_enrich_lower,
                                          min_cell_for_atac_peaks,
                                          if_save_data=F,data_save_folder='') {
  # dstype = ess or noness
  print(paste0('processing ess type id: ', ds_type))
  print(paste0('processing cls id: ', i))
  this_cls <- cls_info_df_this_type[cls_info_df_this_type$cls==i,]
  
  if ((this_cls$'rna.number' < smallest_rna_community) | 
      (this_cls$'atac.number' < smallest_atac_community)) {
    print('this cls is too small, skip this cls')
    return('not selected')
  }
  
  print('this cls is selected')
  this_cls_info_newline <- c(paste0(ds_type,'-',i), 
      this_cls$'rna.number',
      this_cls$'atac.number'
    )
  
  # extract this cls RNA ds
  print('    proecssing RNA data ...')
  cell_ids_RNA <- rownames(rna_ds@meta.data[rna_ds@meta.data[[cls_col_name]]==i,])
  sce_this_cls <- rna_ds[,cell_ids_RNA]
  expr_table <- sce_this_cls@assays$RNA@counts %>% rowSums
  
  # save gene pscount to txt
  if (if_save_data) {
    save_path_rna_this_cls <- paste0(ds_type,'-',i,'.rna-expr-pscounts.txt')
    save_path_rna_this_cls <- file.path(data_save_folder, save_path_rna_this_cls)
    write.table(as.data.frame(expr_table), 
                save_path_rna_this_cls,
                quote=F, sep='\t', col.names=F)
  }

  # extract this cls ATAC ds, and take the subset of the peaks
  print('    proecssing ATAC data ...')
  cell_ids_atac <- rownames(atac_ds@meta.data[atac_ds@meta.data$predicted.id==i,])
  atac_this_cls <- atac_ds[,cell_ids_atac]
  # save ATAC ps-peak to bed
  atac_peaks_dt <- atac_this_cls@assays$peaks@counts %>% as.matrix
  atac_peaks_dt[atac_peaks_dt!=0] <- 1 
  # after this step, convert read counts to existence of cells (1 means there is a peak)

  atac_peaks_dt <- rowSums(atac_peaks_dt) # how many cells have this peak (fragment)
  atac_peaks_dt <- atac_peaks_dt[atac_peaks_dt>=min_cell_for_atac_peaks] %>% as.data.frame %>% rownames
  
  peak_df_to_write <- data.frame(matrix(ncol = 3, nrow = length(atac_peaks_dt)))
  for (j in (1:length(atac_peaks_dt))) {
    peak_df_to_write[j,] <- str_split(atac_peaks_dt[j], pattern='-')[[1]]
  }
  peak_df_to_write <- peak_df_to_write[peak_df_to_write$X1 %in% c(paste0('chr', 1:22),'chrX','chrY'),]
  
  # save to bed file
  if (if_save_data) {
    sv_bed_file_path <- paste0(ds_type,'-',i,'.atac-peaks.bed')
    sv_bed_file_path <- file.path(data_save_folder, sv_bed_file_path)
    write.table(peak_df_to_write, 
                sv_bed_file_path,
                quote=F, sep='\t', col.names=F, row.names=F)
  }
  
  return(this_cls_info_newline)
  }

##################################################################################################
# for each cls in the datasets, apply the filter, extract and save data
print('processing ess:')
ds_type <- 'ess'
cls_info_df_this_type <- cls_info[cls_info$id==ds_type, ]
for (i in cls_info_df_this_type$cls) {
  this_cls_info_newline <- extract_psbulk_data_for_cls_i(
    i=i, 
    rna_ds=sce_ess, 
    atac_ds=ess_atac_ds, 
    cls_info_df_this_type=cls_info_df_this_type,
    ds_type='ess',
    smallest_rna_community=smallest_rna_community, 
    smallest_atac_community=smallest_atac_community,
    min_cell_for_atac_peaks=min_cell_for_atac_peaks,
    if_save_data=T,
    data_save_folder=data_save_folder)
  if (all(this_cls_info_newline != 'not selected')) {
    cls_info_remained[nrow(cls_info_remained)+1,] <- this_cls_info_newline
  }
}

print('processing noness:')
ds_type <- 'noness'
cls_info_df_this_type <- cls_info[cls_info$id==ds_type, ]
for (i in cls_info_df_this_type$cls) {
  this_cls_info_newline <- extract_psbulk_data_for_cls_i(
    i=i, 
    rna_ds=sce_noness, 
    atac_ds=noness_atac_ds, 
    cls_info_df_this_type=cls_info_df_this_type,
    ds_type='noness',
    smallest_rna_community=smallest_rna_community, 
    smallest_atac_community=smallest_atac_community,
    min_cell_for_atac_peaks=min_cell_for_atac_peaks,
    if_save_data=T,
    data_save_folder=data_save_folder)
  if (all(this_cls_info_newline != 'not selected')) {
    cls_info_remained[nrow(cls_info_remained)+1,] <- this_cls_info_newline
  }
}

print('done creating all the datasets')
print('writing input.tsv')
# create input.tsv file based on the cls_info_remained
write_input_tsv(data_save_folder, cls_info_remained)

# create config.yml file
mk_config_file(data_save_folder)

# save cls_info_remained
write.table(cls_info_remained,
            file.path(data_save_folder, 'cls-info-remained.txt'),
            sep='\t', row.names=F, quote=F)

print('all done!')




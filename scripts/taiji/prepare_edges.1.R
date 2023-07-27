library(stringr)
library(data.table)
library(patchwork)
library(cowplot)
library(dplyr)
library(car)
library(ComplexHeatmap)
library(tidyr)
library(tidyverse)
library(ggsignif)
library(gplots)
library(circlize)
library(ComplexHeatmap)
library(cluster)
library(pheatmap)
library(RColorBrewer)
library(heatmaply)
library(factoextra)
library(ggplotify)
library(cluster)
library(ggpubr)
library(rstatix)
library(preprocessCore)
library(multcomp)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)


########################################
args <- commandArgs(trailingOnly = TRUE)

taiji_out_folder <- args[1] # taiji_results
save_folder_base <- args[2] # taiji_results_analysis


meta_df <- read.csv(
    file.path(save_folder_base,'cls-5.rep-0','meta_df.2.txt'),
    header=T, row.names=1)

network_folder_base <- file.path(taiji_out_folder, 'Network') 
network_folder_save_base <- file.path(save_folder_base, 'Network_info', 'step0.filtered_edges.all')
network_folder_save_base_2 <- file.path(save_folder_base, 'Network_info', 'step1.filtered_edges.TF')
network_folder_save_base_3 <- file.path(save_folder_base, 'Network_info', 'step2.combined_edges')


dir.create(network_folder_save_base, recursive = TRUE)
dir.create(network_folder_save_base_2, recursive = TRUE)
dir.create(network_folder_save_base_3, recursive = TRUE)


######################################## 
top_edge_perc <- 0.2
# load all edges then filter by perc, keep only top ones
dt_edge <- read.csv(file.path(
  network_folder_base, 'WT', 'edges_combined.csv'
)) %>% as.data.frame()
dt_edge <- dt_edge[order(dt_edge$weight, decreasing=T), ] 
dt_edge <- dt_edge[1:(as.integer(top_edge_perc*dim(dt_edge)[1])), ]
rownames(dt_edge) <- paste('WT.edge', (1:nrow(dt_edge)), sep='.')
dt_edge$edge_expr_type <- 'WT'
dt_edge <- dt_edge %>% mutate(percent_rank = rank(weight)/length(weight)) 
# find the perc_rank (fraction) for weight, highest=1, lowest=0 (lowest means the top `top_edge_perc` edge)
# X.START_ID    X.END_ID   weight            X.TYPE edge_expr_type
# WT.edge.1        SP1       Y_RNA 3400.211 COMBINED_REGULATE             WT
# WT.edge.2        SP1 METAZOA_SRP 2153.200 COMBINED_REGULATE             WT
# WT.edge.3       CTCF       Y_RNA 2102.746 COMBINED_REGULATE             WT
# WT.edge.4      HMGA1       Y_RNA 2053.504 COMBINED_REGULATE             WT
# WT.edge.5        SP3       Y_RNA 2022.834 COMBINED_REGULATE             WT
# percent_rank
# WT.edge.1    1.0000000
# WT.edge.2    0.9999990
# WT.edge.3    0.9999979
# WT.edge.4    0.9999969
# WT.edge.5    0.9999959

# save dt_edge
write.csv(dt_edge, file.path(network_folder_save_base, paste0('filtered_edges.','WT','.txt')), quote=F)

for (i in (meta_df$taiji_id_old)) {
  print(i)
  if (i=='WT') {
    next
  }
  dt_edge <- read.csv(file.path(network_folder_base, i, 'edges_combined.csv')) %>% as.data.frame 
  
  dt_edge <- dt_edge[order(dt_edge$weight, decreasing=T), ] 
  dt_edge <- dt_edge[1:(as.integer(top_edge_perc*dim(dt_edge)[1])), ]
  rownames(dt_edge) <- paste(i,'edge', (1:nrow(dt_edge)), sep='.')
  dt_edge$edge_expr_type <- i
  dt_edge <- dt_edge %>% mutate(percent_rank = rank(weight)/length(weight)) 
  write.csv(dt_edge, file.path(network_folder_save_base, paste0('filtered_edges.',i,'.txt')), quote=F)
}

######################################## 
# subset the edges associated with the significant TFs
tf_stats.sig_TFs <- read.csv(
  file.path(save_folder_base,'cls-5.rep-0','TF_results','tf_stats.all.2.sig_TFs.FC-2.csv'),
  header=T, row.names=1)
tf_stats.sig_TFs.down <- read.csv(
  file.path(save_folder_base,'cls-5.rep-0','TF_results','tf_stats.all.2.sig_TFs.down.FC-2.csv'),
  header=T, row.names=1)

TF_of_interest <- c(tf_stats.sig_TFs$TF, tf_stats.sig_TFs.down$TF) %>% unique
TF_of_interest %>% write.table(file=file.path(network_folder_save_base_2, 'TF_of_interest.txt'),
                               quote=F, row.names=F, col.names=F)
TF_of_interest <- read.table(file.path(network_folder_save_base_2, 'TF_of_interest.txt'))
TF_of_interest <- TF_of_interest$V1

# WT as example
dt_edge <- read.csv(file.path(
  network_folder_save_base, 
  paste0('filtered_edges.','WT','.txt')
)) %>% as.data.frame()
dt_edge <- dt_edge %>% dplyr::filter(X.START_ID %in% TF_of_interest)
write.csv(dt_edge, file.path(network_folder_save_base_2, paste0('filtered_edges.','WT','.txt')), quote=F)

for (i in (meta_df$taiji_id_old)) {
  print(i)
  if (i=='WT') {
    next
  }
  dt_edge <- read.csv(file.path(
    network_folder_save_base, 
    paste0('filtered_edges.',i,'.txt')
  )) %>% as.data.frame()
  dt_edge <- dt_edge %>% dplyr::filter(X.START_ID %in% TF_of_interest)
  write.csv(dt_edge, file.path(network_folder_save_base_2, paste0('filtered_edges.',i,'.txt')), quote=F)
}

######################################## 
# perform the next step in the next script: for each TF, extract all X.start==TF, then save to `network_folder_save_base_3`







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

args = commandArgs(trailingOnly=TRUE)
start_idx <- args[1] %>% as.integer
taiji_out_folder <- args[2]
save_folder_base <- args[3]


net_in_folder_path = file.path(save_folder_base, 'Network_info', 'step1.filtered_edges.TF')
net_out_folder_path = file.path(save_folder_base, 'Network_info', 'step2.combined_edges')


meta_df <- read.csv(
  file.path(save_folder_base,'cls-5.rep-0','meta_df.2.txt'),
  header=T, row.names=1)

meta_df

TF_of_interest <- read.table(file.path(net_in_folder_path, 'TF_of_interest.txt'))
TF_of_interest <- TF_of_interest$V1


if (start_idx <= length(TF_of_interest)) {
  TF <- TF_of_interest[start_idx]
  # list to save results
  regulatee_history <- c()
  psbulk_history <- c()
  perc_history <- c()
  
  # load all the edge files in Networks.1, then filter X.start in TFs, then save to 
  for (i in (meta_df$taiji_id_old)) {
    print(i)
    dt_edge <- read.csv(file.path(
      net_in_folder_path, 
      paste0('filtered_edges.',i,'.txt')
    )) %>% as.data.frame()
    dt_edge <- dt_edge %>% dplyr::filter(X.START_ID == TF)

    regulatee_history <- c(regulatee_history, dt_edge$X.END_ID)
    perc_history <- c(perc_history, dt_edge$percent_rank)
    psbulk_history <- c(psbulk_history, rep(i, dim(dt_edge)[1]))
    
  }
  
  result <- data.frame(
    'End' = regulatee_history,
    'perc' = perc_history,
    'psbulk' = psbulk_history
  )
  
  # save result to Networks.2
  write.csv(result, 
            file.path(net_out_folder_path, paste0('combined_edges.',TF,'.txt')), 
            quote=F)
  
}

print('All is successful!')








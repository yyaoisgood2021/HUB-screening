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
library(tibble)


########################################
args <- commandArgs(trailingOnly = TRUE)

taiji_out_folder <- args[1]
save_folder_base <- args[2]
meta_df_old_path <- args[3]


down_sample_id <- 'orig'  # {'1', '2', 'orig', 'all'}
TFs_for_pca <- 'top250' # use `all|top250` TFs to run PCA
repeat_id <- 0

N_cl <- 5 # number of kmeans clusters, need to be derived


########################################
# strategy : independently proc orig, rep1 and rep2, all the way down until sig TFs are derived, then compared commonly identified sig TFs      
down_sample_id <- 'orig'  # {'1', '2', 'orig', 'all'}


save_folder_K <- file.path(save_folder_base,
                          paste0('cls-',N_cl,'.rep-',repeat_id))

# # these folders to store edge info
# network_folder_save_base_2 <- '/stg3/data3/peiyao/HUBS/pair_can/taiji_run/20221021.separate-run/post_taiji_result/npcs30_pc30.3/top250TF/Networks.2'
# network_folder_save_base_3 <- '/stg3/data3/peiyao/HUBS/pair_can/taiji_run/20221021.separate-run/post_taiji_result/npcs30_pc30.3/top250TF/Networks.3'
# 

# top_edge_perc <- 0.2 # extract top ??% edges


# meta_df <- read.csv(
#   file.path(fig_save_folder_K,
#             'meta_df.2.txt'),
#   header=T, row.names=1)

# meta_df$km.cluster <- meta_df$km.cls
# meta_df$km.cls <- NULL
# write.csv(meta_df,
#           file.path(fig_save_folder_K, 'meta_df.2.txt'),
#           quote=F)


# meta_df


########################################
dir.create(save_folder_K, recursive = TRUE)
########################################
# load dataset
pr_all <- read.table(
  file.path(taiji_out_folder,'GeneRanks.tsv'),
)
pr_all <- pr_all[ , (pr_all %>% colnames) != 'WT.1']
pr_all %>% colnames


# preprocess and build psbulk meta data
meta_df <- read.table(meta_df_old_path,
                      header=T)

meta_df_new <- data.frame(matrix(0, nrow=0, ncol=6))
colnames(meta_df_new) <- c('taiji_id', 'type', 'seurat_cls', 'down_sample_rep',
                           'rna_cell_number', 'atac_cell_number')

for (i in (pr_all %>% colnames)) {
  if (i=='WT') {
    meta_df_new[nrow(meta_df_new)+1, ] <- c('WT', 'WT', 0, 'orig', 0, 0)
  } 
  else{
    i_split <- str_split(i, pattern='_')[[1]]
    type_i <- i_split[1]
    seurat_cls_i <-  i_split[2]
    if (type_i=='ess') {
      rep_i <- i_split[length(i_split)]
      if (rep_i!=down_sample_id) {next}
    } else {
      # noness type
      rep_i <- 'orig'}
    meta_df_this <- meta_df %>% dplyr::filter(id==type_i & cls==seurat_cls_i)
    meta_df_new[nrow(meta_df_new)+1, ] <- c(i, type_i, seurat_cls_i, rep_i, 
                            meta_df_this$rna.number%>%as.integer, 
                            meta_df_this$atac.number%>%as.integer)
  }
}

# write.csv(meta_df_new,
#           file.path(save_folder_base, 'meta_df.0.txt'),
#           quote=F)

meta_df <- meta_df_new %>% as.data.frame
mycolor <- c('#E64B35B2', '#4DBBD5B2', 'gray64') # ess, noness, WT
type_colors <- mycolor[meta_df$type %>% as.factor %>% as.numeric]
meta_df$type_color <- type_colors

write.csv(meta_df,
          file.path(save_folder_base,
                    'meta_df.0.txt'),
          quote=F)


pr_all <- pr_all[, meta_df$taiji_id] # select subset first, then run scale and PCA and so on ...
pr_all %>% dim # 1165 34

all_pr_data <- pr_all %>% as.matrix %>% t %>% scale # each row is a cell (psbulk)
all_pr_data[1:7,1:5]
#                 AC023509.3 AC138696.1        AHR       AIRE       ALX1
# WT              -1.4177775 -1.2059005 -1.2841344 -1.1774701 -1.1675198
# ess_0_sub_1     -0.6007327 -0.5658289 -0.5686747 -0.5512303 -0.5734685
# ess_0_sub_2     -0.6001719 -0.5658680 -0.6116521 -0.5296163 -0.5786832
# ess_0_sub_orig  -0.7797576 -0.7027925 -0.7575637 -0.7094926 -0.6980438
# ess_11_sub_1     3.6646364  3.7246669  3.6864395  3.6644503  3.6111540
# ess_11_sub_2     3.6003549  3.6581208  3.6978145  3.6568681  3.6591186
# ess_11_sub_orig  0.3641087  0.3985433  0.3510228  0.3813426  0.4902824


######################################## 
# pca, study if features are well reduced
# select a subset of TFs, and all the downstream analysis is based on these TFs
My_top_n_list <- list('top250'=250,
                      'all'=ncol(all_pr_data))

TopTFs <- sapply(all_pr_data %>% heatmaply::normalize(), var) %>% 
  as.data.frame %>% dplyr::rename('var' = '.') %>% 
  arrange(desc(var)) %>% 
  top_n(My_top_n_list[[TFs_for_pca]]) %>%
  row.names

input <- all_pr_data %>% as.data.frame %>% dplyr::select(TopTFs) %>% scale
input[1:6,1:6]

# PCA and performance of PCA
pca <- prcomp(input, scale = F) # perform pca only on the selected TFs
# then only look at first 4 PCs
data_reduced <- as.data.frame(pca$x) %>% dplyr::select(PC1, PC2, PC3, PC4)
data_reduced 

mycolor <- c('#E64B35B2', '#4DBBD5B2', 'gray64') # ess, noness, WT
col_to_use <- mycolor[(c(rep('WT',1), rep('ess', 15), rep('noness', 18)) 
                       %>% as.factor %>% as.numeric)]
names(col_to_use) = row.names(data_reduced)
col_to_use


#%% PCA distr of first 4 pcs
pdf(file.path(save_folder_base, 'first_four_PCs.1.pdf'))
plot(data_reduced %>% select_if(is.numeric), col=col_to_use)
# legend(1, legend=c("Ess", "Noness", 'WT'), col=c('#FF0000','#009E73','black'))
dev.off()

#%% plot explained var
pdf(file.path(save_folder_base, 'PC.var_explained.0.pdf'))
fviz_eig(pca,
         addlabels = T, 
         barcolor = "#E7B800", 
         barfill = "#E7B800", 
         linecolor = "#00AFBB", 
         choice = "variance", 
         ylim=c(0,60))
dev.off()

# save pca
saveRDS(pca, file.path(save_folder_base, "pca_dataset.rds"))

pca <- readRDS(file.path(save_folder_base, "pca_dataset.rds"))


######################################## 
# find optimal K for kmeans
pdf(file.path(save_folder_base, "fviz_nbclust.repeat.pdf"), 
    # width = 7, height = 6
)
print(fviz_nbclust(data_reduced, kmeans, method = "wss"))
print(fviz_nbclust(data_reduced, kmeans, method = "silhouette"))
gap_stat <- clusGap(data_reduced, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)
print(fviz_gap_stat(gap_stat))
dev.off()




findOptimal <- function(data_reduced, output_file){
  
  # elbow method
  fviz_nbclust(data_reduced, kmeans, method = "wss", k.max = 10)
  
  # silhouette method
  g <- fviz_nbclust(data_reduced, kmeans, method = "silhouette",k.max = 10)
  df2 <- g$data
  
  # # use pearson distance
  # res.dist <- get_dist(data_reduced, method = "pearson")
  # g <- fviz_nbclust(data_reduced, kmeans, method = "silhouette",diss=res.dist,k.max = 10)
  # df2 <- rbind(df2,g$data)
  
  # # use spearman distance
  # res.dist <- get_dist(data_reduced, method = "spearman")
  # g <- fviz_nbclust(data_reduced, kmeans, method = "silhouette",diss=res.dist,k.max = 10)
  # df2 <- rbind(df2,g$data)
  
  # # use kendall distance
  # res.dist <- get_dist(data_reduced, method = "kendall")
  # g <- fviz_nbclust(data_reduced, kmeans, method = "silhouette",diss=res.dist,k.max = 10)
  # df2 <- rbind(df2,g$data)
  
  # use manhattan distance
  res.dist <- get_dist(data_reduced, method = "manhattan")
  g <- fviz_nbclust(data_reduced, kmeans, method = "silhouette",diss=res.dist,k.max = 10)
  df2 <- rbind(df2,g$data)
  
  # visualization
  # df2$method <- c(rep("euclidean",10),rep('pearson',10),rep('spearman',10),rep('kendall',10),rep('manhattan',10))
  df2$method <- c(rep("euclidean",10),rep('manhattan',10))
  
  df2[,'clusters'] <- as.numeric(as.character(df2[,'clusters']))
  df2 <- df2[df2$clusters>=1,]
  p2 <- ggplot(df2)+aes(x=clusters,y=y,color=method,group=method)+
    geom_line()+geom_point()+
    labs(x='number of clusters K',y='average silhouette width',color='distance metrics')+
    scale_x_continuous(breaks =seq(0,10,by=5))+
    # theme(axis.title = element_text(size = 15),
    #       axis.text = element_text(size = 15),
    #       legend.title = element_text(size = 15),
    #       legend.text = element_text(size = 15)
    #       )+
    theme_light()
  # png("distanceMetrics.png",units="in", width=6, height=5, res=300)
  pdf(output_file,width=8.5,height = 5)
  print(p2)
  dev.off()
}

findOptimal(data_reduced, 
            file.path(save_folder_base, "opt_kmean_cluster.1.pdf"))



meta_df <- read.csv(file.path(
  save_folder_base, 
  'meta_df.0.txt'), 
  header=T, row.names=1)

# perform kmeans on multiple N_cl val
for (N_cl_i in 4:9) {
  set.seed(0) # cluster results are the same but still gave random cluster id
  cl <- kmeans(data_reduced, N_cl_i, iter.max = 10000)
  cluster.assignment <- cl$cluster  %>% as.data.frame 
  colnames(cluster.assignment) <- paste0('km.cls.', N_cl_i)
  cluster.assignment <- cluster.assignment %>% rownames_to_column('taiji_id')
  meta_df <- merge(meta_df, cluster.assignment, by.x='taiji_id') %>% as.data.frame
}

meta_df <- meta_df %>% arrange(km.cls.8)
meta_df

######################################## 
# now, only consider N_cl<-{}, first we have to set WT to km.cls=0
col_to_use <- paste0('km.cls.',N_cl)
WT_cls_id <- meta_df[meta_df$taiji_id=='WT', col_to_use] 
# change this id to 0, 0~WT_cls_id-1 become 1~WT_cls_id, >WT_cls_id keep unchanged
cls_new <- c()
for (cls_i in meta_df[[col_to_use]]) {
  if (cls_i<WT_cls_id) {cls_new <- c(cls_new, cls_i)}
  else if (cls_i==WT_cls_id) {cls_new <- c(cls_new, 0)}
  else {cls_new <- c(cls_new, cls_i-1)}
  
}
meta_df$km.cls <- cls_new
meta_df <- meta_df %>% arrange(km.cls)
meta_df # psbulk metadata is prepared
write.csv(data.frame(meta_df),
            file.path(save_folder_K,
                      'meta_df.1.txt'),
            quote=F)


######################################## 
# prepare km.cluster metadata
# meta_df: noness.0 has ??? atac.number and ??? rna.number
# need to prep a df like this
# df: prepare two dfs, for rna and atac
# km.cluster type `number of cells` label_ypos
# 1      noness    200  200
# 1       ess      100  200+100=300
meta_df$'km.cluster' <- meta_df$km.cls
meta_df$rna.number <- meta_df$rna_cell_number
meta_df$atac.number <- meta_df$atac_cell_number

new_taiji_id <- c()
for (psbulk_name in meta_df$taiji_id) {
  print(psbulk_name)
  if (psbulk_name=='WT') {
    new_taiji_id <- c(new_taiji_id, psbulk_name)
  } else {
    spl_values <- str_split(psbulk_name, '_')[[1]]  
    new_taiji_id <- c(new_taiji_id, paste(spl_values[1], spl_values[2], sep='_'))
  }
  
}
meta_df$taiji_id_old <- meta_df$taiji_id
meta_df$taiji_id <- new_taiji_id


write.csv(data.frame(meta_df),
          file.path(save_folder_K,
                    'meta_df.2.txt'),
          quote=F)
  
  
  
count_df_rna <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(count_df_rna) <- c('km.cluster', 'type', 'number', 'label_ypos')

count_df_atac <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(count_df_atac) <- c('km.cluster', 'type', 'number', 'label_ypos')

for (i in (meta_df$'km.cluster' %>% unique)){
  meta_this_kmcls <- meta_df[meta_df$'km.cluster'==i, c('rna.number', 'atac.number', 'type')]
  print(i)
  print(meta_this_kmcls)
  
  this_cls_rna.ess <- meta_this_kmcls[meta_this_kmcls$'type'=='ess','rna.number'] %>% as.integer %>% sum
  this_cls_rna.noness <- meta_this_kmcls[meta_this_kmcls$'type'=='noness','rna.number'] %>% as.integer %>% sum
  this_cls_atac.ess <- meta_this_kmcls[meta_this_kmcls$'type'=='ess','atac.number'] %>% as.integer %>% sum
  this_cls_atac.noness <- meta_this_kmcls[meta_this_kmcls$'type'=='noness','atac.number'] %>% as.integer %>% sum
  
  count_df_rna[nrow(count_df_rna)+1, ] <- c(
    i,
    'noness',
    this_cls_rna.noness,
    this_cls_rna.noness
  )
  count_df_rna[nrow(count_df_rna)+1, ] <- c(
    i,
    'ess',
    this_cls_rna.ess,
    this_cls_rna.ess+this_cls_rna.noness+10
  )
  
  count_df_atac[nrow(count_df_atac)+1, ] <- c(
    i,
    'noness',
    this_cls_atac.noness,
    this_cls_atac.noness
  )
  count_df_atac[nrow(count_df_atac)+1, ] <- c(
    i,
    'ess',
    this_cls_atac.ess,
    this_cls_atac.noness+this_cls_atac.ess+10
  )
  
}

count_df_rna <- count_df_rna[count_df_rna$number != 0, ]
count_df_atac <- count_df_atac[count_df_atac$number != 0, ]

count_df_rna$number <- count_df_rna$number %>% as.integer
count_df_atac$number <- count_df_atac$number %>% as.integer
count_df_rna$label_ypos <- count_df_rna$label_ypos %>% as.integer
count_df_atac$label_ypos <- count_df_atac$label_ypos %>% as.integer

my_type_color <- list('ess'='#E64B35B2','noness'='#4DBBD5B2')
count_df_rna$color <- my_type_color[count_df_rna$type] %>% as.character
count_df_atac$color <- my_type_color[count_df_atac$type] %>% as.character

print(count_df_rna)
print(count_df_atac)

write.csv(count_df_rna, 
          file.path(save_folder_K,
                    'kmcls_meta_df.rna.txt'),
          quote=F)
write.csv(count_df_atac, 
          file.path(save_folder_K,
                    'kmcls_meta_df.atac.txt'),
          quote=F)


# pdf(file.path(save_folder_K, 'rna_number_bar.pdf'))
# ggplot(data=count_df_rna, aes(x=km.cluster, y=number, fill=type)) +
#   geom_bar(stat="identity") +
#   scale_fill_manual("type", values = my_type_color) +
#   geom_text(aes(y=label_ypos, label=number), vjust=1.6,
#             color="black", size=4.5)+
#   ggtitle("") +
#   ylab('Number of cells') + xlab("km.cluster") +
#   theme(axis.title = element_text(size = 18), 
#         axis.text = element_text(size = 12),
#         legend.title = element_text(size = 15),
#         legend.text = element_text(size = 12)) 
# dev.off()

# pdf(file.path(save_folder_K, 'atac_number_bar.pdf'))
# ggplot(data=count_df_atac, aes(x=km.cluster, y=number, fill=type)) +
#   geom_bar(stat="identity") +
#   scale_fill_manual("type", values = my_type_color) +
#   geom_text(aes(y=label_ypos, label=number), vjust=1.6,
#             color="black", size=4.5)+
#   ggtitle("") +
#   ylab('Number of cells') + xlab("km.cluster") +
#   theme(axis.title = element_text(size = 18), 
#         axis.text = element_text(size = 12),
#         legend.title = element_text(size = 15),
#         legend.text = element_text(size = 12)) 
# dev.off()



# cell_counts : QC filtered cell (cells in psbulk) 
cls_info_0 <- read.table(meta_df_old_path, 
                       sep='\t', header=T)
cells_psbulk_rna_ess <- cls_info_0[cls_info_0$id=='ess', 'rna.number'] %>% sum # 9321
cells_psbulk_rna_noness <- cls_info_0[cls_info_0$id=='noness', 'rna.number'] %>% sum # 5899
cells_psbulk_atac_ess <- cls_info_0[cls_info_0$id=='ess', 'atac.number'] %>% sum # 3554
cells_psbulk_atac_noness <- cls_info_0[cls_info_0$id=='noness', 'atac.number'] %>% sum # 2559

# cell_counts :  cells in taiji input 
cells_kmcls_rna_ess <- meta_df[meta_df$type=='ess', 'rna_cell_number'] %>% sum # 5780
cells_kmcls_rna_noness <- meta_df[meta_df$type=='noness', 'rna_cell_number'] %>% sum # 4884
cells_kmcls_atac_ess <- meta_df[meta_df$type=='ess', 'atac_cell_number'] %>% sum # 2997
cells_kmcls_atac_noness <- meta_df[meta_df$type=='noness', 'atac_cell_number'] %>% sum # 2383

adjust_height_low <- 0.5
adjust_height_high <- 0.5

count_df_rna <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(count_df_rna) <- c('km.cluster', 'type', 'number', 'percentage', 'label_ypos.number', 'label_ypos.percentage')

count_df_atac <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(count_df_atac) <- c('km.cluster', 'type', 'number', 'percentage', 'label_ypos.number', 'label_ypos.percentage')

for (i in (meta_df$'km.cluster' %>% unique)){
  meta_this_kmcls <- meta_df[meta_df$'km.cluster'==i, c('rna_cell_number', 'atac_cell_number', 'type')]
  print(i)
  print(meta_this_kmcls)
  
  this_cls_rna.ess <- meta_this_kmcls[meta_this_kmcls$'type'=='ess','rna_cell_number'] %>% as.integer %>% sum
  this_cls_rna.noness <- meta_this_kmcls[meta_this_kmcls$'type'=='noness','rna_cell_number'] %>% as.integer %>% sum
  this_cls_atac.ess <- meta_this_kmcls[meta_this_kmcls$'type'=='ess','atac_cell_number'] %>% as.integer %>% sum
  this_cls_atac.noness <- meta_this_kmcls[meta_this_kmcls$'type'=='noness','atac_cell_number'] %>% as.integer %>% sum
  
  count_df_rna[nrow(count_df_rna)+1, ] <- c(
    i,
    'noness',
    this_cls_rna.noness,
    this_cls_rna.noness / cells_psbulk_rna_noness * 100,
    this_cls_rna.noness,
    this_cls_rna.noness / cells_psbulk_rna_noness * 100 - adjust_height_low
  )
  count_df_rna[nrow(count_df_rna)+1, ] <- c(
    i,
    'ess',
    this_cls_rna.ess,
    this_cls_rna.ess / cells_psbulk_rna_ess * 100,
    this_cls_rna.ess+this_cls_rna.noness+10,
    this_cls_rna.noness / cells_psbulk_rna_noness * 100 + this_cls_rna.ess / cells_psbulk_rna_ess * 100 + adjust_height_high
  )
  
  count_df_atac[nrow(count_df_atac)+1, ] <- c(
    i,
    'noness',
    this_cls_atac.noness,
    this_cls_atac.noness / cells_psbulk_atac_noness * 100,
    this_cls_atac.noness,
    this_cls_atac.noness / cells_psbulk_atac_noness * 100 - adjust_height_low
  )
  count_df_atac[nrow(count_df_atac)+1, ] <- c(
    i,
    'ess',
    this_cls_atac.ess,
    this_cls_atac.ess / cells_psbulk_atac_ess * 100,
    this_cls_atac.ess+this_cls_atac.noness+10,
    this_cls_atac.noness / cells_psbulk_atac_noness * 100 + this_cls_atac.ess / cells_psbulk_atac_ess * 100 + adjust_height_high
  )
  
}

# add `undetermined` to km.cls
count_df_rna[nrow(count_df_rna)+1, ] <- c(
  'undetermined',
  'noness',
  cells_psbulk_rna_noness - cells_kmcls_rna_noness,
  (cells_psbulk_rna_noness - cells_kmcls_rna_noness) / cells_psbulk_rna_noness * 100,
  cells_psbulk_rna_noness - cells_kmcls_rna_noness,
  (cells_psbulk_rna_noness - cells_kmcls_rna_noness) / cells_psbulk_rna_noness * 100 - adjust_height_low
)
count_df_rna[nrow(count_df_rna)+1, ] <- c(
  'undetermined',
  'ess',
  cells_psbulk_rna_ess - cells_kmcls_rna_ess,
  (cells_psbulk_rna_ess - cells_kmcls_rna_ess) / cells_psbulk_rna_ess * 100,
  (cells_psbulk_rna_noness - cells_kmcls_rna_noness)+(cells_psbulk_rna_ess - cells_kmcls_rna_ess)+10,
  (cells_psbulk_rna_noness - cells_kmcls_rna_noness) / cells_psbulk_rna_noness * 100 + (cells_psbulk_rna_ess - cells_kmcls_rna_ess) / cells_psbulk_rna_ess * 100 + adjust_height_high
)

count_df_atac[nrow(count_df_atac)+1, ] <- c(
  'undetermined',
  'noness',
  cells_psbulk_atac_noness - cells_kmcls_atac_noness,
  (cells_psbulk_atac_noness - cells_kmcls_atac_noness) / cells_psbulk_atac_noness * 100,
  cells_psbulk_atac_noness - cells_kmcls_atac_noness,
  (cells_psbulk_atac_noness - cells_kmcls_atac_noness) / cells_psbulk_atac_noness * 100 - adjust_height_low
)
count_df_atac[nrow(count_df_atac)+1, ] <- c(
  'undetermined',
  'ess',
  cells_psbulk_atac_ess - cells_kmcls_atac_ess,
  (cells_psbulk_atac_ess - cells_kmcls_atac_ess) / cells_psbulk_atac_ess * 100,
  (cells_psbulk_atac_noness - cells_kmcls_atac_noness)+(cells_psbulk_atac_ess - cells_kmcls_atac_ess)+10,
  (cells_psbulk_atac_noness - cells_kmcls_atac_noness) / cells_psbulk_atac_noness * 100 + (cells_psbulk_atac_ess - cells_kmcls_atac_ess) / cells_psbulk_atac_ess * 100 + adjust_height_high        
)

count_df_rna <- count_df_rna[count_df_rna$number != 0, ]
count_df_atac <- count_df_atac[count_df_atac$number != 0, ]

count_df_rna$number <- count_df_rna$number %>% as.integer
count_df_atac$number <- count_df_atac$number %>% as.integer

count_df_rna$label_ypos.number <- count_df_rna$label_ypos.number %>% as.integer
count_df_atac$label_ypos.number <- count_df_atac$label_ypos.number %>% as.integer

count_df_rna$percentage <- count_df_rna$percentage %>% as.numeric %>% round(digits=1)
count_df_atac$percentage <- count_df_atac$percentage %>% as.numeric %>% round(digits=1)

count_df_rna$label_ypos.percentage <- count_df_rna$label_ypos.percentage %>% as.numeric %>% round(digits=1)
count_df_atac$label_ypos.percentage <- count_df_atac$label_ypos.percentage %>% as.numeric %>% round(digits=1)


my_type_color <- list('ess'='#E64B35B2','noness'='#4DBBD5B2')
count_df_rna$color <- my_type_color[count_df_rna$type] %>% as.character
count_df_atac$color <- my_type_color[count_df_atac$type] %>% as.character


print(count_df_rna)
print(count_df_atac)

write.csv(count_df_rna,
          file.path(save_folder_K,
                    'kmcls_meta_df.percentage.rna.txt'),
          quote=F)
write.csv(count_df_atac,
          file.path(save_folder_K,
                    'kmcls_meta_df.percentage.atac.txt'),
          quote=F)


pdf(file.path(save_folder_K, 'rna_number_bar.pdf'))
ggplot(data=count_df_rna, aes(x=km.cluster, y=percentage, fill=type)) +
  geom_bar(stat="identity") +
  scale_fill_manual("type", values = my_type_color) +
  geom_text(aes(y=label_ypos.percentage, label=percentage), vjust=1.6,
            color="black", size=4.5)+
  ggtitle("") +
  ylab('Percentage of cells') + xlab("km.cluster") +
  theme(axis.title = element_text(size = 18), 
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12)) 
dev.off()

pdf(file.path(save_folder_K, 'atac_number_bar.pdf'))
ggplot(data=count_df_atac, aes(x=km.cluster, y=percentage, fill=type)) +
  geom_bar(stat="identity") +
  scale_fill_manual("type", values = my_type_color) +
  geom_text(aes(y=label_ypos.percentage, label=percentage), vjust=1.6,
            color="black", size=4.5)+
  ggtitle("") +
  ylab('Percentage of cells') + xlab("km.cluster") +
  theme(axis.title = element_text(size = 18), 
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12)) 
dev.off()


######################################## 
# pheatmap for all the data
# def function
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}


new_input_rownames <- c()
for (psbulk_name in rownames(input)) {
  print(psbulk_name)
  if (psbulk_name=='WT') {
    new_input_rownames <- c(new_input_rownames, psbulk_name)
  } else {
    spl_values <- str_split(psbulk_name, '_')[[1]]  
    new_input_rownames <- c(new_input_rownames, paste(spl_values[1], spl_values[2], sep='_'))
  }
  
}
rownames(input) <- new_input_rownames



pheatmap_func <- function(input, meta_df) {
  pheat_input <- input %>% as.data.frame %>% as.matrix # then we need to reorder based on km.cluster (WT, ...)
  pheat_input_new <- merge(pheat_input, meta_df%>%column_to_rownames('taiji_id'), by=0) %>% column_to_rownames('Row.names')
  pheat_input_new$type_number <- pheat_input_new$type %>% as.factor %>% as.numeric
  pheat_input_new$type_number <- 4 - pheat_input_new$type_number
  pheat_input_new <- pheat_input_new[order(pheat_input_new$km.cluster, pheat_input_new$type_number, pheat_input_new$seurat_cls),]
  reordered_meta_df <- pheat_input_new[,
                                       ((My_top_n_list[[TFs_for_pca]]+1) : ncol(pheat_input_new))]
  pheat_input_new <- pheat_input_new[,1:My_top_n_list[[TFs_for_pca]]] %>% t
  
  print(reordered_meta_df)
  print(pheat_input_new %>% dim) # 250 34
  
  mat_breaks <- quantile_breaks(pheat_input_new, n = 500)
  
  # make manual breaks
  reordered_meta_df$km.cluster <- paste('', reordered_meta_df$km.cluster, sep='') %>% unlist
  
  
  km_indexes = reordered_meta_df %>% pull(km.cluster) %>% table %>% as.data.frame 
  km_indexes$index = cumsum(km_indexes$Freq)
  # Make palette colors for km clusters 
  # get colors 
  my_type_color <- c('ess'='#E64B35B2','noness'='#4DBBD5B2', 'WT'='gray64')
  
  my_palette = brewer.pal(n = N_cl, name = "Set2")
  names(my_palette) = reordered_meta_df %>% pull(km.cluster) %>% unique
  
  AllColors = list(type = my_type_color, km.cluster = my_palette)
  AllColors
  
  pdf(file.path(save_folder_K, 'pheatmap.test.0.pdf'))
  a <- pheatmap::pheatmap(pheat_input_new, 
                     scale = 'none', 
                     cluster_rows = T, 
                     cluster_cols = F,
                     annotation_col = (reordered_meta_df[,c( 'km.cluster', 'type')] %>% as.data.frame),
                     annotation_names_row=F,
                     show_rownames = F,
                     clustering_distance_rows = 'euclidean',
                     clustering_method = 'average',
                     clustering_distance_cols = 'euclidean',
                     breaks = mat_breaks,
                     color = viridis(500),
                     annotation_colors = AllColors,
                     gaps_col = km_indexes$index,
                     legend = T)
  print(a)
  dev.off()
}

pheatmap_func(input, meta_df)


######################################## 
# MWU on query cluster and bg cluster
meta_df <- read.csv(file.path(save_folder_K,
                              'meta_df.2.txt'), header=T, row.names=1)
all.pageranks <- read.table(
  file.path(taiji_out_folder,'GeneRanks.tsv')
) %>% t

all.pageranks <- all.pageranks[meta_df$taiji_id_old, ]

print(meta_df)
all.pageranks[1:6,1:6]


# mwu_method_1 : one-vs-WT : background (bg) is WT-km.cls; while query is ess (1,4) or noness (2,3,5) specific km.cls
# mwu_method_4 : one-vs-all : bg is all the other km.clusters
meta_df <- meta_df %>% column_to_rownames('taiji_id_old')
all.pageranks[1:6,1:6]
query_cluster_histroy <- c()
query_cls_type_history <- c()
bg_cluster_histroy <- c()
mean_query_histroy <- c()
mean_background_histroy <- c()
mwu_pval_histroy <- c()
TFs_histroy <- c()
mwu_test_strategies_histroy <- c()
all.pageranks.norm <- all.pageranks %>% normalize # original pr score to normalize -> 0~1, then mwu test

ess_km.cls <- c(4) %>% as.character
noness_km.cls <- c(3) %>% as.character
common_km.cls <-c(1,2)
WT_km.cls <- '0'


for (TF in (all.pageranks.norm %>% colnames)){
  TF_vals = all.pageranks.norm %>% dplyr::select(TF)
  # print(TF_vals)
  for (query_km.cls in ess_km.cls) {
    # find the psbulk sample id (rownames) of a given km.cluster
    query_psbulk_samples_cls <- meta_df[meta_df$km.cluster==query_km.cls, ] %>% rownames
    TF_val_query <- TF_vals[query_psbulk_samples_cls, 1] %>% as.matrix
    query_cls_type_history <- c(query_cls_type_history, rep('ess',2))
    query_cluster_histroy <- c(query_cluster_histroy, rep(query_km.cls,2))
    mean_query_histroy <- c(mean_query_histroy, rep(mean(TF_val_query), 2))
    TFs_histroy <- c(TFs_histroy, rep(TF,2))
    
    mwu_test_strategies_histroy <- c(mwu_test_strategies_histroy, 'one-vs-WT', 'one-vs-all')
    bg_cluster_histroy <- c(bg_cluster_histroy, 'WT', 'all')
    
    # method1 one-vs-WT
    bg_samples_cls <- meta_df[meta_df$km.cluster==WT_km.cls, ] %>% rownames
    TF_val_bg <- TF_vals[bg_samples_cls, 1] %>% as.matrix
    test_query_vs_bg <- wilcox.test(TF_val_query, TF_val_bg)$p.val
    mean_background_histroy <- c(mean_background_histroy, mean(TF_val_bg))
    mwu_pval_histroy <- c(mwu_pval_histroy, test_query_vs_bg)
    # method4
    bg_samples_cls <- meta_df[meta_df$km.cluster != query_km.cls, ] %>% rownames
    TF_val_bg <- TF_vals[bg_samples_cls, 1] %>% as.matrix
    test_query_vs_bg <- wilcox.test(TF_val_query, TF_val_bg)$p.val
    mean_background_histroy <- c(mean_background_histroy, mean(TF_val_bg))
    mwu_pval_histroy <- c(mwu_pval_histroy, test_query_vs_bg)
    
  }
}

for (TF in (all.pageranks.norm %>% colnames)){
  TF_vals = all.pageranks.norm %>% dplyr::select(TF)
  # print(TF_vals)
  for (query_km.cls in noness_km.cls) {
    # find the psbulk sample id (rownames) of a given km.cluster
    query_psbulk_samples_cls <- meta_df[meta_df$km.cluster==query_km.cls, ] %>% rownames
    TF_val_query <- TF_vals[query_psbulk_samples_cls, 1] %>% as.matrix
    query_cls_type_history <- c(query_cls_type_history, rep('noness',2))
    query_cluster_histroy <- c(query_cluster_histroy, rep(query_km.cls,2))
    mean_query_histroy <- c(mean_query_histroy, rep(mean(TF_val_query), 2))
    TFs_histroy <- c(TFs_histroy, rep(TF,2))
  
    mwu_test_strategies_histroy <- c(mwu_test_strategies_histroy, 'one-vs-WT', 'one-vs-all')
    bg_cluster_histroy <- c(bg_cluster_histroy, 'WT', 'all')
  
    # method1 one-vs-WT
    bg_samples_cls <- meta_df[meta_df$km.cluster==WT_km.cls, ] %>% rownames
    TF_val_bg <- TF_vals[bg_samples_cls, 1] %>% as.matrix
    test_query_vs_bg <- wilcox.test(TF_val_query, TF_val_bg)$p.val
    mean_background_histroy <- c(mean_background_histroy, mean(TF_val_bg))
    mwu_pval_histroy <- c(mwu_pval_histroy, test_query_vs_bg)
    # method4
    bg_samples_cls <- meta_df[meta_df$km.cluster != query_km.cls, ] %>% rownames
    TF_val_bg <- TF_vals[bg_samples_cls, 1] %>% as.matrix
    test_query_vs_bg <- wilcox.test(TF_val_query, TF_val_bg)$p.val
    mean_background_histroy <- c(mean_background_histroy, mean(TF_val_bg))
    mwu_pval_histroy <- c(mwu_pval_histroy, test_query_vs_bg)
  }

}
  
for (TF in (all.pageranks.norm %>% colnames)){
  TF_vals = all.pageranks.norm %>% dplyr::select(TF)
  # print(TF_vals)
  for (query_km.cls in common_km.cls) {
    # find the psbulk sample id (rownames) of a given km.cluster
    query_psbulk_samples_cls <- meta_df[meta_df$km.cluster==query_km.cls, ] %>% rownames
    TF_val_query <- TF_vals[query_psbulk_samples_cls, 1] %>% as.matrix
    query_cls_type_history <- c(query_cls_type_history, rep('common',2))
    query_cluster_histroy <- c(query_cluster_histroy, rep(query_km.cls,2))
    mean_query_histroy <- c(mean_query_histroy, rep(mean(TF_val_query), 2))
    TFs_histroy <- c(TFs_histroy, rep(TF,2))
    
    mwu_test_strategies_histroy <- c(mwu_test_strategies_histroy, 'one-vs-WT', 'one-vs-all')
    bg_cluster_histroy <- c(bg_cluster_histroy, 'WT', 'all')
    
    # method1 one-vs-WT
    bg_samples_cls <- meta_df[meta_df$km.cluster==WT_km.cls, ] %>% rownames
    TF_val_bg <- TF_vals[bg_samples_cls, 1] %>% as.matrix
    test_query_vs_bg <- wilcox.test(TF_val_query, TF_val_bg)$p.val
    mean_background_histroy <- c(mean_background_histroy, mean(TF_val_bg))
    mwu_pval_histroy <- c(mwu_pval_histroy, test_query_vs_bg)
    # method4
    bg_samples_cls <- meta_df[meta_df$km.cluster != query_km.cls, ] %>% rownames
    TF_val_bg <- TF_vals[bg_samples_cls, 1] %>% as.matrix
    test_query_vs_bg <- wilcox.test(TF_val_query, TF_val_bg)$p.val
    mean_background_histroy <- c(mean_background_histroy, mean(TF_val_bg))
    mwu_pval_histroy <- c(mwu_pval_histroy, test_query_vs_bg)
  }
  
}


tf_stats.all_km_cls_df <- data.frame(
  'TF'=TFs_histroy,
  'query_type'=query_cls_type_history,
  'query'=query_cluster_histroy,
  'background'=bg_cluster_histroy,
  'query_mean'=mean_query_histroy,
  'background_mean'=mean_background_histroy,
  'mwu_p'=mwu_pval_histroy,
  'strategy'=mwu_test_strategies_histroy
)  
tf_stats.all_km_cls_df %>% dim  # 9320     8
tf_stats.all_km_cls_df$FC = tf_stats.all_km_cls_df$query_mean / tf_stats.all_km_cls_df$background_mean 
tf_stats.all_km_cls_df$pval.adjust = p.adjust(tf_stats.all_km_cls_df$mwu_p, method = 'BH')
tf_stats.all_km_cls_df %>% head(10)

# save raw 
dir.create(file.path(save_folder_K,'TF_results'), recursive = TRUE)
write.csv(data.frame(tf_stats.all_km_cls_df), 
          file.path(save_folder_K,
                    'TF_results',
                    'tf_stats.all.1.csv'),
          quote=F)
  


######################################## 
# extract significant TFs
tf_stats.all <- read.csv(file.path(save_folder_K,'TF_results',
                                   'tf_stats.all.1.csv'),
                         header=T, row.names=1)
# filter TF stats
FC_thresh <- 2
tf_stats.sig_TFs <- tf_stats.all %>% filter(pval.adjust < 0.05 & FC > FC_thresh & background == 'WT')
write.csv(data.frame(tf_stats.sig_TFs), 
          file.path(save_folder_K,'TF_results',
                    paste0('tf_stats.all.2.sig_TFs.FC-',FC_thresh,'.csv')),
          quote=F)

tf_stats.sig_TFs.down <- tf_stats.all %>% filter(pval.adjust < 0.05 & FC < 1/FC_thresh & background == 'WT')
write.csv(data.frame(tf_stats.sig_TFs.down), 
          file.path(save_folder_K,'TF_results',
                    paste0('tf_stats.all.2.sig_TFs.down.FC-',FC_thresh,'.csv')),
          quote=F)










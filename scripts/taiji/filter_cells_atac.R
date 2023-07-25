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
# load atac dataset from out folder, and add QC metrics and annotations
load_atac_data <- function(sample_id, data_folder, ds_save_folder, annotations,
                           genome_assembly='hg38', min.cells=20, min.features=1000) {
  print('loading counts ...')
  counts <- Read10X_h5(filename = file.path(data_folder, "filtered_peak_bc_matrix.h5"))
  print('loading metadata ...')
  metadata <- read.csv(
    file = file.path(data_folder, "singlecell.csv"),
    header = TRUE,
    row.names = 1
  )
  print('creating assay ...')
  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    genome = genome_assembly,
    fragments = file.path(data_folder, 'fragments.tsv.gz'),
    min.cells = min.cells,
    min.features = min.features
  )
  cn_from_count_data = colnames(chrom_assay@counts) 
  metadata_sub = metadata[rownames(metadata) %in% cn_from_count_data,]
  print('creating object ...')
  atac_ds <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata_sub
  )
  atac_ds@meta.data[['orig.ident']] <- sample_id
  # add annotations
  print('adding annotations ...')
  Annotation(atac_ds) <- annotations
  
  # QC of atac data
  print('QC ...')
  atac_ds <- NucleosomeSignal(object = atac_ds)
  atac_ds <- TSSEnrichment(atac_ds)
  atac_ds$pct_reads_in_peaks <- atac_ds$peak_region_fragments / atac_ds$passed_filters * 100
  # atac_ds$blacklist_ratio <- atac_ds$blacklist_region_fragments / atac_ds$peak_region_fragments # always 0
  
  # save rds, already derived all the QCs, but before filtering
  print('saving rds ...')
  saveRDS(atac_ds, 
          file.path(ds_save_folder, 
                    paste0("atac.",sample_id,".s0.rds")))
  print('done')
  return(atac_ds)
}

##################################################################################################
filter_low_quality_atac <- function(atac_ds, sample_id, ds_save_folder, plot_save_folder,
                                    peak_region_lower_bound = 3000, peak_region_higher_bound = 20000,
                                    peak_pct_lower_bound = 15, nucleo_higher_bound = 4,
                                    TSS_lower_bound = 2) {
  # remove low quality data
  atac_ds_remained <- subset(
    x = atac_ds,
    subset = peak_region_fragments > peak_region_lower_bound &
      peak_region_fragments < peak_region_higher_bound &
      pct_reads_in_peaks > peak_pct_lower_bound &
      nucleosome_signal < nucleo_higher_bound &
      TSS.enrichment > TSS_lower_bound
  )
  
  # plot QC for datasets before and after filter
  # pdf(file.path(plot_save_folder,paste0('atac.',sample_id,'.before-qc.pdf')), 
  #     height=4, width=10
  # )
  # print(VlnPlot(
  #   object = atac_ds,
  #   features = c('pct_reads_in_peaks', 'peak_region_fragments',
  #                'TSS.enrichment', 'nucleosome_signal'),
  #   pt.size = 0.1,
  #   ncol = 4
  # ))
  # dev.off()
  
  # pdf(file.path(plot_save_folder,paste0('atac.',sample_id,'.after-qc.pdf')),
  #     height=4, width=10
  # )
  # print(VlnPlot(
  #   object = atac_ds_remained,
  #   features = c('pct_reads_in_peaks', 'peak_region_fragments',
  #                'TSS.enrichment', 'nucleosome_signal'),
  #   pt.size = 0.1,
  #   ncol = 4
  # ))
  # dev.off()
  
  # save RDS after filter
  saveRDS(atac_ds_remained, 
          file.path(ds_save_folder, 
                    paste0("atac.",sample_id,".s1.rds")))
  
  print('done')
  return(atac_ds_remained)
  
}





args <- commandArgs(trailingOnly = TRUE)

ess_atac_path <- args[1]
noness_atac_path <- args[2]
sv_path <- args[3]



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

# process atac data
ess_atac_all <- load_atac_data('ess', ess_atac_path, sv_path, annotations,
                           genome_assembly='hg38', min.cells=20, min.features=1000)
noness_atac_all <- load_atac_data('noness', noness_atac_path, sv_path, annotations,
                               genome_assembly='hg38', min.cells=20, min.features=1000)  
  
# filter
ess_atac_remain <- filter_low_quality_atac(ess_atac_all, 'ess', sv_path, '',
                                    peak_region_lower_bound = 3000, peak_region_higher_bound = 20000,
                                    peak_pct_lower_bound = 15, nucleo_higher_bound = 4,
                                    TSS_lower_bound = 2)
noness_atac_remain <- filter_low_quality_atac(noness_atac_all, 'noness', sv_path, '',
                                           peak_region_lower_bound = 3000, peak_region_higher_bound = 20000,
                                           peak_pct_lower_bound = 15, nucleo_higher_bound = 4,
                                           TSS_lower_bound = 2)

print('done')

# N2.atac_remain <- readRDS('/stg3/data3/peiyao/HUBS/pair_can/ATAC-10x/results/20220926/rds/atac.N2.s1.rds')
# N7.atac_remain <- readRDS('/stg3/data3/peiyao/HUBS/pair_can/ATAC-10x/results/20220926/rds/atac.N7.s1.rds')








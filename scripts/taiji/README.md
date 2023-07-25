# Procedures for Pseudobulk Taiji analysis
1. Download single cell experiment data (scRNA-seq and scATAC-seq experiments, for both Ess and Noness libraries) from GEO [(GSE231384)](https://www.ncbi.xyz/geo/query/acc.cgi?acc=GSE231384), raw sequencing data can be obtained from the associated SRA.

2. Process scRNA-seq raw sequencing data. Use `cellranger count` to count fastq files for both Ess and Noness libraries, then use `cellranger aggr` to combine two libraries. Run all the steps with the default parameters following the [10x Genomics instructions](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) (v7.0.0). The sample csv file `cellranger_aggr_r.csv` required by the `cellranger aggr` command can be obtained from [resources](https://github.com/yyaoisgood2021/HUB-screening/tree/main/resources/taiji).
    ```bash
    # count the Ess library (named as N2 in my codes)
    cellranger count --id=N2
                     --transcriptome=/stg3/data3/peiyao/software/refdata-gex-GRCh38-2020-A \
                     --fastqs=/stg3/data3/peiyao/HUBS/pair_can/scRNA10x/data_N2 \
                     --sample=22070FL-01-02-01
    # also count the Noness library (named as N7 in my codes). change the sample id, fastqs path and sample prefix accordingly

    cellranger aggr --id=combined \
                    --csv=resources/taiji/cellranger_aggr_r.csv \
                    --nosecondary
    ```

3. Generate clusters of single cells based on scRNA-seq data using Seurat (v4.1.0) [(ref)](https://satijalab.org/seurat/). To fulfill this, you need to run `scripts/taiji/filter_cells.R` to remove low quality single cells, then run `scripts/taiji/mk-indv-rna-clusters.R`. The statistics files for the clustering results `ess_rna.cls_info.npcs30_pc30.txt` and `noness_rna.cls_info.npcs30_pc30.txt` can be found in [resources](https://github.com/yyaoisgood2021/HUB-screening/tree/main/resources/taiji).  
    ```bash
    Rscript scripts/taiji/filter_cells.R seurat_sv_folder_1 aggr_out_folder
    # seurat_sv_folder_1 is the folder to save this command's output. for this and all output_sv_folders below, make sure you have already created these folders
    # cellranger aggr's output folder should be cellranger_prefix/outs/count/filtered_feature_bc_matrix
    
    Rscript scripts/taiji/mk-indv-rna-clusters.R path/to/filtered_seurat_obj seurat_sv_folder_2 n_pc n_pc_cls 
    # path/to/filtered_seurat_obj: the saved seurat object from the last step, should be seurat_sv_folder_1/hubs.high_quality.combined.s1.rds
    # seurat_sv_folder_2: save folder for this command
    # n_pc, n_pc_cls: hyperparameters for PCA and single-cell clustering. I'm using 30, 30 in my manuscript
    ```

4. Process scATAC-seq raw sequencing data. Use `cellranger-atac count` with the default parameters following the [10x Genomics instructions](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/using/count) (v7.0.0).

    ```bash
    # count the Ess library (named as N2 in my codes)
    cellranger-atac count --id=N2 \
                          --reference=/stg3/data3/peiyao/software/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                          --fastqs=/stg3/data3/peiyao/HUBS/pair_can/ATAC-10x/data_N2
    
    # also count the Noness library (named as N7 in my codes). change the sample id and fastqs path accordingly
    ```

5. Generate pseudobulk clusters of single cells for Taiji inputs. Briefly, you need to filter out the low-quality single cells in the scATAC experiments, then integrate scATAC onto the scRNA experiments, and finally extract and sum the gene counts and ATAC fragments for the single cells from the respective pseudobulk clusters. You can achieve this by sequentially running `scripts/taiji/filter_cells_atac.R`, `scripts/taiji/mk-psbulk-data.0.R` and `scripts/taiji/mk-psbulk-data.1.R`. All hyperparameters used can be found in the scripts.

    ```bash
    Rscript scripts/taiji/filter_cells_atac.R path/to/cellranger-atac/ess/out path/to/cellranger-atac/noness/out seurat_atac_sv_folder_1
    # path/to/cellranger-atac/ess/out: outs/ folder of the cellranger-atac count command, for the Ess library. should contain filtered_peak_bc_matrix.h5
    # path/to/cellranger-atac/noness/out: out/ folder for the Noness library
    # seurat_atac_sv_folder_1: folder to save the output of this command 
    
    Rscript scripts/taiji/mk-psbulk-data.0.R path/to/ess_clustered_rna_seurat_obj path/to/noness_clustered_rna_seurat_obj \
    path/to/atac_save_folder 
    

    # path/to/ess_clustered_rna_seurat_obj: the path to the RNA-seq seurat object including clustering results for the ess library.
    # It is the output of step 3, and should be seurat_sv_folder_2/ess_rna.with_cluster_info.npcs30_pc30.s2.rds
    # path/to/noness_clustered_rna_seurat_obj: should be seurat_sv_folder_2/noness_rna.with_cluster_info.npcs30_pc30.s2.rds
    # path/to/atac_save_folder: this should be the result_save_folder of the previous step. ie, seurat_atac_sv_folder_1
    
    ```

6. Also preprae gene count file and ATAC fragment file for the K562 WT control.

7. Run Taiji following the [instruction](https://taiji-pipeline.github.io/). A sample `taiji_config.yml` file and `taiji_input.tsv` file can be found in [resources](https://github.com/yyaoisgood2021/HUB-screening/tree/main/resources/taiji).

8. After you got the Taiji results, perform K-Means analysis. Run `scripts/taiji/post_taiji_analysis.0.R` to generate K-Means clusters and to identify top transcription factors (TFs).

10. Finally, run `scripts/taiji/post_taiji_analysis.1.R` to analyze the corresponding regulatees.


# Expected results


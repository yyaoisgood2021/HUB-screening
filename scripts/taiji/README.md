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

5. Filter out the low-quality single cells in the scATAC experiments by running `scripts/taiji/filter_cells_atac.R`

   ```bash
   Rscript scripts/taiji/filter_cells_atac.R \
   path/to/cellranger-atac/ess/out \
   path/to/cellranger-atac/noness/out \
   seurat_atac_sv_folder_1
   # path/to/cellranger-atac/ess/out: outs/ folder of the cellranger-atac count command, for the Ess library. should contain filtered_peak_bc_matrix.h5
   # path/to/cellranger-atac/noness/out: out/ folder for the Noness library
   # seurat_atac_sv_folder_1: folder to save the output of this command 
   ```

6. Generate pseudobulk clusters of single cells for Taiji inputs. Briefly, you need to integrate scATAC onto the scRNA experiments, and then extract and sum the gene counts and ATAC fragments for the single cells from the respective pseudobulk clusters. You can achieve this by sequentially running `scripts/taiji/mk-psbulk-data.0.R` and then `scripts/taiji/mk-psbulk-data.1.R`. All hyperparameters used can be found in the scripts.

    ```bash
    Rscript scripts/taiji/mk-psbulk-data.0.R \
    path/to/ess_clustered_rna_seurat_obj path/to/noness_clustered_rna_seurat_obj \
    path/to/atac_save_folder n_pc n_pc_cls res \
    seurat_atac_sv_folder_2
    # path/to/ess_clustered_rna_seurat_obj: the path to the RNA-seq seurat object including the clustering results for the ess library.
    # It is the output of step 3, and should be seurat_sv_folder_2/ess_rna.with_cluster_info.npcs30_pc30.s2.rds
    # path/to/noness_clustered_rna_seurat_obj: should be seurat_sv_folder_2/noness_rna.with_cluster_info.npcs30_pc30.s2.rds
    # path/to/atac_save_folder: this should be the result_save_folder of the previous step. ie, seurat_atac_sv_folder_1
    # n_pc, n_pc_cls, res: hyperparameters used in the PCA and single-cell clustering steps, refer to the step 3. I'm using 30, 30, 3 here
    # seurat_atac_sv_folder_2: folder to save the output of this command
    ```
    After this command, you will see two seurat objects for the filtered and anchored atac-seq experiments, you will also see a statistics file called `stat_df.npcs30_pc30.3.txt`, you can check the numbers in [resources](https://github.com/yyaoisgood2021/HUB-screening/tree/main/resources/taiji).

    ```bash
    Rscript scripts/taiji/mk-psbulk-data.1.R \
    path/to/seurat/rna/folder \
    path/to/seurat/atac/folder \
    n_pc n_pc_cls res \
    taiji_data_sv_folder_base 
        
    # path/to/seurat/rna/folder: the output folder of step 3, it should be seurat_sv_folder_2 and it should contain two processed seurat files for the rna experiment: {ess|noness}_rna.with_cluster_info.npcs30_pc30.s1.rds
    # path/to/seurat/atac/folder: the output folder of the previous step, it should be seurat_atac_sv_folder_2 and it should contain two seurat files for the atac experiments: atac-after-anchor-transfer.{ess|noness}.npcs30_pc30.3.s2.rds
    # n_pc, n_pc_cls, res: hyperparameters used in the PCA and single-cell clustering steps, refer to the step 3.
    # refer to the parameters in the previous steps. I'm using 30, 30, 3 here
    # taiji_data_sv_folder_base: the base folder to save outputs of this command,
    # when using n_pc=30, n_pc_cls=30, and res=3, output are saved under taiji_data_sv_folder_base/npcs30_pc30.3
    ```




7. (optional) Also preprae gene count file and ATAC fragment file for the K562 WT control.

8. Generate the `config.yml` and `input.yml` files

9. Run Taiji following the [instruction](https://taiji-pipeline.github.io/). A sample `taiji_config.yml` file and `taiji_input.yml` file can be found in [resources](https://github.com/yyaoisgood2021/HUB-screening/tree/main/resources/taiji).

10. After you got the Taiji results, perform K-Means analysis. Run `scripts/taiji/post_taiji_analysis.0.R` to generate K-Means clusters and to identify top transcription factors (TFs).

11. Finally, run `scripts/taiji/post_taiji_analysis.1.R` to analyze the corresponding regulatees.


# Expected results


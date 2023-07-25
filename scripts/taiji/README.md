# Procedures for Pseudobulk Taiji analysis
1. Download single cell experiment data (scRNA-seq and scATAC-seq experiments, for both Ess and Noness libraries) from GEO [(GSE231384)](https://www.ncbi.xyz/geo/query/acc.cgi?acc=GSE231384), raw sequencing data can be obtained from the associated SRA.

2. Process scRNA-seq raw sequencing data. Use `cellranger count` to count fastq files for both Ess and Noness libraries, then use `cellranger aggr` to combine two libraries. Run all the steps with the default parameters following the [10x Genomics instructions](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) (v7.0.0). The sample csv file `cellranger_aggr_r.csv` required by the `cellranger aggr` command can be obtained from [resources](https://github.com/yyaoisgood2021/HUB-screening/tree/main/resources/taiji).

3. Generate clusters of single cells based on scRNA-seq data using Seurat (v4.1.0) [(ref)](https://satijalab.org/seurat/). To fulfill this, you need to run `scripts/taiji/filter_cells.R` to remove low quality single cells, then run `scripts/taiji/prep_clusters.R`.  

4. Process scATAC-seq raw sequencing data. Use `cellranger-atac count` with the default parameters following the [10x Genomics instructions](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/using/count).

5. Generate pseudobulk clusters of single cells for Taiji inputs. Briefly, you need to filter out the low-quality single cells in the scATAC experiments, then integrate scATAC onto the scRNA experiments, and finally extract and sum the gene counts and ATAC fragments for the single cells from the respective pseudobulk clusters. You can achieve this by running `scripts/taiji/mk_psbulk_data.R`. All hyper-parameters can be found in the scripts.   

6. Also preprae gene count file and ATAC fragment file for the K562 WT control.

7. Run Taiji following the [instruction](https://taiji-pipeline.github.io/). A sample `taiji_config.yml` file and `taiji_input.tsv` file can be found in [resources](https://github.com/yyaoisgood2021/HUB-screening/tree/main/resources/taiji).

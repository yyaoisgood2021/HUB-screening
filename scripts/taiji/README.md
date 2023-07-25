# Procedures for Pseudobulk Taiji analysis
1. download single cell experiment data from GEO [GSE231384](https://www.ncbi.xyz/geo/query/acc.cgi?acc=GSE231384), raw sequencing data can be obtained from the associated SRA.

2. process scRNA-seq data: using `cellranger count` and`cellranger aggr` using the default parameters following the [10x Genomics instructions](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger). The sample csv file required by cellranger aggr can be obtained from [resources](https://github.com/yyaoisgood2021/HUB-screening/tree/main/resources/taiji).

3. run `scripts/taiji/filter_cells.R` and `scripts/taiji/prep_clusters.R` 



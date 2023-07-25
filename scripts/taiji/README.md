1. download single cell experiment data from GEO [GSE231384](https://www.ncbi.xyz/geo/query/acc.cgi?acc=GSE231384), raw sequencing data can be obtained from the associated SRA.

process scRNA-seq data: using cellranger count and cellranger aggr using the default parameters following the 10x Genomics instructions. The sample csv file required by cellranger aggr can be obtained from resources

run scripts/Pseudobulk_Taiji

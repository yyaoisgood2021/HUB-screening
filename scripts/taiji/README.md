# Procedures for Pseudobulk Taiji analysis
1. Create a folder structure according to the [folder_tree.txt](https://github.com/yyaoisgood2021/HUB-screening/blob/main/folder_tree.txt). Download single cell experiment raw sequencing data to the corresponding sub-folders of `results/data`. scRNA-seq and scATAC-seq experiments, for both Ess and Noness libraries. Use SRA Toolkit from SRA (corresponding to the GEO entry [GSE231384](https://www.ncbi.xyz/geo/query/acc.cgi?acc=GSE231384)). The commands below is an example to download ess RNA-seq data:
   ```bash
   mkdir -p results/data/RNA/ess
   cd results/data/RNA/ess
   prefetch SRR24376123 
   ```
   then unzip the fastq.gz files, and repeat the codes to download fastq files for the noness library and for ATAC-seq experiments. You can also download the processed data and put them in the corresponding sub-folders of `results/proc_data`.

2. Process scRNA-seq raw sequencing data and save to sub-folders of `results/proc_data/RNA` (creating folder structrue referring to [folder_tree.txt](https://github.com/yyaoisgood2021/HUB-screening/blob/main/folder_tree.txt), same for all codes below). Use `cellranger count` to count fastq files for both Ess and Noness libraries, then use `cellranger aggr` to combine two libraries. If you download the processed files from GEO, you can put them in the corresponding sub-folders of results/proc_data/RNA and directly run `cellranger aggr`. Run all the steps with the default parameters following the [10x Genomics instructions](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) (v7.0.0). The sample csv file `cellranger_aggr_r.csv` required by the `cellranger aggr` command can be obtained from [resources](https://github.com/yyaoisgood2021/HUB-screening/tree/main/resources/taiji).
    ```bash
    # count the ess library 
    cellranger count --id=ess
                     --transcriptome=/stg3/data3/peiyao/software/refdata-gex-GRCh38-2020-A \
                     --fastqs=results/data/RNA/ess \
                     --sample=22070FL-01-02-01
    ```
    also count the Noness library. change the sample id, fastqs path and sample prefix accordingly.

    then, combine Ess and Noness libraries. 
    ```bash
    cellranger aggr --id=combined \
                    --csv=resources/taiji/cellranger_aggr_r.csv \
                    --nosecondary
    ```
3. Generate clusters of single cells based on scRNA-seq data using Seurat (v4.1.0) [(ref)](https://satijalab.org/seurat/). To fulfill this, you need to run `scripts/taiji/filter_cells.R` to remove low quality single cells, then run `scripts/taiji/mk-indv-rna-clusters.R` to cluster single cells. The statistics files for the clustering results `ess_rna.cls_info.npcs30_pc30.txt` and `noness_rna.cls_info.npcs30_pc30.txt` can be found in [resources](https://github.com/yyaoisgood2021/HUB-screening/tree/main/resources/taiji).  
    ```bash
    Rscript scripts/taiji/filter_cells.R \
    results/seurat_RNA/filtered \
    results/proc_data/RNA/combined/outs/count/filtered_feature_bc_matrix
    # results/seurat_RNA/filtered is the folder to save this command's output.
    # for this and all output_save_folders below, make sure you have already created these folders

    n_pc=30
    n_pc_cls=30
    Rscript scripts/taiji/mk-indv-rna-clusters.R \
    results/seurat_RNA/filtered/hubs.high_quality.combined.s1.rds \
    results/seurat_RNA/clustered \
    ${n_pc} ${n_pc_cls} 
    # results/seurat_RNA/filtered/hubs.high_quality.combined.s1.rds is the seurat object generated from the last step
    # results/seurat_RNA/clustered: save folder for this command
    # n_pc, n_pc_cls: hyperparameters for PCA and single-cell clustering. I'm using 30, 30 in my manuscript
    ```

4. Process scATAC-seq raw sequencing data. Use `cellranger-atac count` with the default parameters following the [10x Genomics instructions](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/using/count) (v7.0.0). If you downloaded the processed data from GEO, you can simply put them in the folder results/proc_data/ATAC, then you need to prepare an index file by running `tabix -p bed fragments.tsv.gz` in the folder results/proc_data/ATAC, then you can skip this step.

    ```bash
    cellranger-atac count --id=ess \
                          --reference=/stg3/data3/peiyao/software/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                          --fastqs=results/data/ATAC/ess
    ```
    also count the Noness library. change the sample id and fastqs path accordingly
   
5. Filter out the low-quality single cells in the scATAC experiments by running `scripts/taiji/filter_cells_atac.R`

   ```bash
   Rscript scripts/taiji/filter_cells_atac.R \
   results/proc_data/ATAC/ess/outs \
   results/proc_data/ATAC/noness/outs \
   results/seurat_ATAC/filtered
   
   # results/proc_data/ATAC/ess/outs: the outs/ folder of the cellranger-atac count command, for the Ess library.
   # It should contain important files such as filtered_peak_bc_matrix.h5
   # results/proc_data/ATAC/noness/outs: outs/ folder for the Noness library
   # results/seurat_ATAC/filtered: the folder to save the outputs of this command, the files to save the filtered ATAC objects are atac.{ess|noness}.s1.rds 
   ```

6. Generate pseudobulk clusters of single cells for Taiji inputs. Briefly, you need to integrate scATAC onto the scRNA experiments, and then extract and sum the gene counts and ATAC fragments for the single cells from the respective pseudobulk clusters. You can achieve this by sequentially running `scripts/taiji/mk-psbulk-data.0.R` and then `scripts/taiji/mk-psbulk-data.1.R`. All hyperparameters used can be found in the scripts.

    ```bash
    n_pc=30
    n_pc_cls=30
    res=3
    Rscript scripts/taiji/mk-psbulk-data.0.R \
    results/seurat_RNA/clustered/ess_rna.with_cluster_info.npcs30_pc30.s2.rds \
    results/seurat_RNA/clustered/noness_rna.with_cluster_info.npcs30_pc30.s2.rds \
    results/seurat_ATAC/filtered \
    ${n_pc} ${n_pc_cls} ${res} \ 
    results/seurat_ATAC/clustered
    
    # results/seurat_RNA/clustered/ess_rna.with_cluster_info.npcs30_pc30.s2.rds:
    # the path to the RNA-seq seurat object including the clustering results for the ess library. It is the output of step 3.
    # results/seurat_RNA/clustered/ess_rna.with_cluster_info.npcs30_pc30.s2.rds:
    # the rna-seurat path for noness library
    # results/seurat_ATAC/filtered: this should be the result_save_folder of the previous step
    # n_pc, n_pc_cls, res: hyperparameters used in the PCA and single-cell clustering steps, refer to the step 3. I'm using 30, 30, 3 here
    # results/seurat_ATAC/clustered2: the folder to save the outputs of this command
    ```
    After this command, you will see two seurat objects for the filtered and anchored atac-seq experiments, you will also see a statistics file called `stat_df.npcs30_pc30.3.txt`, you can check the numbers in [resources](https://github.com/yyaoisgood2021/HUB-screening/tree/main/resources/taiji).

   Run the following codes to prepare data:

    ```bash
    n_pc=30
    n_pc_cls=30
    res=3
    Rscript scripts/taiji/mk-psbulk-data.1.R \
    results/seurat_RNA/clustered \
    results/seurat_ATAC/clustered \
    ${n_pc} ${n_pc_cls} ${res} \ 
    results/taiji_datasets \
    results/proc_data/ATAC/ess/outs/fragments.tsv.gz \
    results/proc_data/ATAC/noness/outs/fragments.tsv.gz 
        
    # results/seurat_RNA/clustered: the output folder of step 3, it should contain two processed seurat files for the rna experiment: {ess|noness}_rna.with_cluster_info.npcs30_pc30.s2.rds
    # results/seurat_ATAC/clustered: the output folder of the previous step, and it should contain two seurat files for the atac experiments: atac-after-anchor-transfer.{ess|noness}.npcs30_pc30.3.s2.rds
    # n_pc, n_pc_cls, res: hyperparameters used in the PCA and single-cell clustering steps, refer to the step 3.
    # refer to the parameters in the previous steps. I'm using 30, 30, 3 here
    # results/taiji_datasets: the base folder to save outputs of this command, when using n_pc=30, n_pc_cls=30, and res=3, output are saved under results/taiji_datasets/npcs30_pc30.3
    # results/proc_data/ATAC/ess/outs/fragments.tsv.gz: the raw fragment file of the ATAC experiment, we will extract fragments from it
    ```
    After this command, you will see the gene count tsv file and atac fragment bed.gz file for each pseudobulk cluster. A statistics file named `cls-info-remained.txt` is also generated. You can check the numbers in [resources](https://github.com/yyaoisgood2021/HUB-screening/tree/main/resources/taiji)
   
7. Also preprae gene count tsv file and ATAC narrow-peak bed file for the K562 WT control. Download these data from ENCODE
   
   * ATAC peaks: [ENCFF976CEI](https://www.encodeproject.org/files/ENCFF976CEI/), directly download the bed narrow peak file in the GRCh38 assembly and save to `results/proc_data/ATAC/WT/ENCFF976CEI.bed`. Also download the replicate file from [ENCFF117MSK](https://www.encodeproject.org/files/ENCFF117MSK/)
   * RNA counts:
     
     i. download raw fastq data [ENCSR637VLS](https://www.encodeproject.org/experiments/ENCSR637VLS/) (single-end, all other files are paired-end), and [ENCSR000CPH](https://www.encodeproject.org/experiments/ENCSR000CPH/), [ENCSR000AEM](https://www.encodeproject.org/experiments/ENCSR000AEM/), and [ENCSR000AEO](https://www.encodeproject.org/experiments/ENCSR000AEO/)

     ii. align with [STAR](https://github.com/alexdobin/STAR) (2.7.10a) according to the manual. Use the GRCh38 assembly.
     An example script for running STAR aligner on the first replicate data of ENCSR000AEM is shown below.

     * repeat the codes below for all other replicates and for all other accessions.
     ```bash
     path_to_index_folder='add_hg38_genome_index_here'

     out_dir=results/proc_data/RNA/WT/ENCSR000AEM
     STAR --genomeDir ${path_to_index_folder} \
     --runThreadN 6 \
     --runMode alignReads \
     --readFilesIn results/data/RNA/WT/ENCSR000AEM/ENCFF001RED.fastq results/data/RNA/WT/ENCSR000AEM/ENCFF001RDZ.fastq \
     --outFileNamePrefix ${out_dir}/rep1-q30- \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattributes Standard \
     --outFilterMultimapNmax 1 \
     --outSAMmapqUnique 30 \
     --quantMode GeneCounts
     ```
     
     iii. STAR will generate a lot of files and you can see `{prefix}-q30-ReadsPerGene.out.tab`, convert gene_id (first column) to gene_short_name (gene symbol) in accordance with the outputs of 10x, use [gencode.annotation.gtf](https://www.gencodegenes.org/human/release_41.html) (v41), the example codes to implement this conversion on the first replicate data of ENCSR000AEM is shown below.
  
     * repeat the codes below for all other replicates and for all other accessions. 
     ```bash
     gencode_file_path='add_gencode_gtf_path_here' 
     prefix_in_this_example='rep1'

     work_folder=results/proc_data/RNA/WT/ENCSR000AEM
     data_to_transfer_path=${work_folder}/${prefix_in_this_example}-q30-ReadsPerGene.out.tab
     data_to_save_path=${work_folder}/${prefix_in_this_example}-q30-ReadsPerGene.name_converted.txt
     python scripts/taiji/convert_gene_short_names.py ${gencode_file_path} ${data_to_transfer_path} ${data_to_save_path}     
     ```
     
     iv. sum up all the results and save it as `results/proc_data/RNA/WT.combined/WT.rna-expr-pscounts.txt`. There's no need to normalize.

     ```bash
     python scripts/taiji/sum_gene_counts.py \
     results/taiji_datasets/npcs30_pc30.3/ess-0.rna-expr-pscounts.txt \
     results/proc_data/RNA/WT.combined/WT.rna-expr-pscounts.txt \
     results/proc_data/RNA/WT/ENCSR000AEM/rep1-q30-ReadsPerGene.name_converted.txt \
     results/proc_data/RNA/WT/ENCSR000AEM/rep2-q30-ReadsPerGene.name_converted.txt \
     {extra file 1} \
     {extra file 2} \
     {extra file N}
     # add all gene count files generated from the prevous steps

     # in this step, you need to run scripts/taiji/sum_gene_counts.py to sum over all (format changed) gene count files generated from the last steps
     # the first argument is a reference from results/taiji_datasets/npcs30_pc30.3. the purpose of adding this reference is to keep the gene set the same.
     # you can use any rna-expr-pscounts.txt file in this folder. if your selected Npcs parameters are different, you should change the folder name accordingly.
     # the second argument is the save path for the output of this command
     # then you add the paths to all gene count files, separate them with `space`
     ```

8. Prepare chromatin 3D interaction files. In this study, we used the top 10 percent of promoter-enhancer loops identified by [EpiTensor](http://wanglab.ucsd.edu/star/EpiTensor/). Run EpiTensor according to its manual and select the top 10 percent of interactions. The processed data can be found in Wang lab's server `/stg3/data1/nas-0-0.bak/share/epi_universal/human/epitensor_loop_top10p_87311.txt`

9. Prepare the `config.yml` and `input.yml` files based on the actual paths that you have. An example of `config.yml` and `input.yml` can be found in [resources](https://github.com/yyaoisgood2021/HUB-screening/blob/main/resources/taiji/). You must change the paths in the `config.yml` manually, then you can save the file with correct paths to `results/taiji_commands/config.yml`. Run the following command to generate `results/taiji_commands/input.yml`. More explanations can be found in the [scripts](https://github.com/yyaoisgood2021/HUB-screening/blob/main/scripts/taiji/prep_input_config.py).
   ```bash
   python scripts/taiji/prep_input_config.py \
   results/taiji_datasets/npcs30_pc30.3/cls-info-remained.txt \
   results/proc_data/RNA/WT.combined/WT.rna-expr-pscounts.txt \
   results/proc_data/ATAC/WT/ENCFF117MSK.bed \
   results/proc_data/ATAC/WT/ENCFF976CEI.bed \
   results/taiji_datasets/npcs30_pc30.3 \
   /stg3/data1/nas-0-0.bak/share/epi_universal/human/epitensor_loop_top10p_87311.txt \
   results/taiji_commands/input.yml
   
   # results/taiji_datasets/npcs30_pc30.3: the path to the folder that you saved the extracted rna and atac data for the sample Ess and Noness libraries. remove the "/" in the end of the string
   # /stg3/data1/nas-0-0.bak/share/epi_universal/human/epitensor_loop_top10p_87311.txt: this is the path to the Epitensor results, we use it to impute Hi-C contacts, you can copy this file to your local folder and modify this line
   # path to save the results
   ```
   
10. Run Taiji following the [instruction](https://taiji-pipeline.github.io/).
    ```bash
    ml load taiji
    progress_folder=results/taiji_progresses # this folder save the progresses so you can resume job
    mkdir -p ${progress_folder}
   
    cd ${progress_folder} # enter progress_folder and submit job
    config_path=results/taiji_commands/config.yml
    taiji run --config ${config_path} --cloud

    # I'd recommend put taiji outputs to "results/taiji_results", you can modify this in the config.yml
    ```

11. After you got the Taiji results, perform PCA and K-Means analysis. Run `scripts/taiji/post_taiji_analysis.0.R` first to generate K-Means clusters and to identify top transcription factors (TFs). Lists of the significantly-changed TFs will be generated in the folder `results/taiji_results_analysis/cls-5.rep-0/TF_results` 
    ```bash
    Rscript scripts/taiji/post_taiji_analysis.0.R \
    results/taiji_results \
    results/taiji_results_analysis \
    results/seurat_ATAC/clustered/stat_df.npcs30_pc30.3.txt 
    
    # results/taiji_results: taiji_out_folder
    # results/taiji_results_analysis: save folder base
    ```
   
12. Finally, run the following commands to derive the significantly-changed (TF -> regulatee gene) edges and the results will be generated in the folder `results/taiji_results_analysis/cls-5.rep-0/edge_results`.
    ```bash 
    Rscript scripts/taiji/prepare_edges.1.R \
    results/taiji_results \
    results/taiji_results_analysis 

    # this part of script requires mem >= 32G
    # and you may want to sbatch each job to a different node with high mem
    # so that the commands inside of the "for loop" can be run in parallel to save time. 
    for tf_id in {1..200}
    do
       Rscript scripts/taiji/weight_changed_edges.2.R ${tf_id} results/taiji_results results/taiji_results_analysis
       python scripts/taiji/find_changed_edges_ks.3.py ${tf_id} results/taiji_results_analysis
    done

    python scripts/taiji/summarize.4.py results/taiji_results_analysis
    ```


# Expected results

1. [folder_tree.txt](https://github.com/yyaoisgood2021/HUB-screening/blob/main/folder_tree.txt) displays the generated outputs in the folder architecture.

2. expected output figures
<p align="center">
  <img  height="800" src="https://github.com/yyaoisgood2021/HUB-screening/blob/main/resources/taiji/Fig5.final.new-1.png">
</p>
<p align="center">
  <img  height="291" src="https://github.com/yyaoisgood2021/HUB-screening/blob/main/resources/taiji/Fig5.final.new-4.png">
</p>

Figure 5 Single cell RNA-seq and ATAC-seq analysis on the hub deletions
(A) UMAP visualization plots of the single cells in the scRNA-seq and scATAC-seq experiments. For scRNA-seq, colors and numbers indicate the clusters identified by SNN modularity optimization (Seurat v4.1.0). For scATAC-seq, numbers indicate the clusters anchored to the corresponding clusters in the scRNA-seq analysis. Cells in each cluster were pooled into pseudobulks. As a result, 15 (Essential) and 18 (Nonessential) pseudobulk samples were used as inputs for the Taiji software.  

(B) Heatmap of the Z-scales of the PageRank scores for the TFs, arranged by sample type and K-Means clusters.

(C) Venn diagram shows the transcription factors (TFs) with significantly increased PageRank scores in the Essential (red) or Nonessential (blue) cluster. The TFs were selected using adjusted p <= 0.05 and fold change >= 2 by comparing them with the WT cluster using a Mann-Whitney U test (see Methods). 

(D) Each gray dot in the TF box (red, essential hub deletion; blue, nonessential hub deletion) represents a TF with a higher PageRank score in the hub-deleted cells listed in Figure 5C. Three TFs, BACH1 (essential-specific), JUNB (common) and ZNF787 (nonessential-specific) are representative TFs, and the outgoing edges represent regulatory interactions with their target genes with significantly increased weight (adjusted p-values <= 0.05 and percentile rank difference >= 0.5 by comparing with the WT and the cluster of the opposite essentiality type, see Methods). It is evident that essential-specific TFs regulate a significantly larger number of genes than common or nonessential-specific TFs. 



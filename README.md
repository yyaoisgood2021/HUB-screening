# HUB-screening

This is the code for the manuscript: "Distinct 3D contacts and phenotypic consequences of adjacent non-coding loci in the epigenetically quiescent regions" to implement to following sections:

### Section 1 "Constructing and training of logistic regression to facilitate candidate selection"

To run `CNN-AE-LR`, you need to prepare the sequence data for each chromosome (chr1 ~ chr22, and chrX) (hg19): download from [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/), and save files to CNN-AE-LR/data/hg19/

Then, you can run codes in scripts/CNN-AE-LR.ipynb


### Section 2 "Classifying hub essentiality with sequence and epigenetic features"

1. To run codes in this section, you need to install [juicer](https://github.com/aidenlab/juicer) (v1.6) and [bedtools](https://bedtools.readthedocs.io/en/latest/index.html).

2. Gather and format data:
   
	i. Download the supplementary table `sciadv.abi6020_table_s11.xlsx` from [Ding et al. study](https://www.science.org/doi/10.1126/sciadv.abi6020) (DOI: 10.1126/sciadv.abi6020), and put it in the designated folder `PR-LR/data`, the expected folder tree structure is displayed in [folder_tree.txt](https://github.com/yyaoisgood2021/HUB-screening/blob/main/folder_tree.txt)
	
 	ii. Download Hi-C data (`GSE63525_K562_combined_30.hic`) from [Rao's work](https://www.cell.com/fulltext/S0092-8674(14)01497-4) (DOI: 10.1016/j.cell.2014.11.021)

	iii. Download all feature data according to file [download_chip_dt.txt](https://github.com/yyaoisgood2021/HUB-screening/blob/main/resources/download_chip_dt.txt), then rename the files according to folder_tree.txt. Also download ATAC-seq narrow peak file from ENCODE [ENCFF976CEI](https://www.encodeproject.org/files/ENCFF976CEI/) (Note ENCFF330SHG and ENCFF976CEI are hg38 and you need to lift them over to hg19)

	iv. Obtain the essential genes data from [Wang et al. 2015](https://www.science.org/doi/10.1126/science.aac7041), [Morgen et al. 2016](https://www.nature.com/articles/nbt.3567), and long non-coding RNA essentiality data from [Liu et al. 2018](https://www.nature.com/articles/nbt.4283). These files have already been downloaded to [resources](https://github.com/yyaoisgood2021/HUB-screening/tree/main/resources). Identify the essential genes and lncRNAs and extract their genomic coordinates using [gencode](https://www.gencodegenes.org/human/release_19.html) (hg19), save the bed files to `PR-LR/data`.

	v. Also extract the genomic coordinates of all genes from gencode (hg19) and save to `PR-LR/data`.

	* When extracting the genomic coordinates, use the coordinates under the `gene` type in the feature column. if there're multiple gene sites, you should keep all of them.

4. then generate `{chrid}_K562_prob.5000.txt` file for each chromosome using the following bash commands

   ```bash
   juicer_path=path/to/juicer_tools.jar # change this line to the actual path to juicer_tools.jar

   sv_folder_1=PR-LR/proc_HiC/hic.dump
   sv_folder_2=PR-LR/proc_HiC/inter.prob.ref
   mkdir -p ${sv_folder_1}
   mkdir -p ${sv_folder_2}

   for chr_num in {1..22} X
   do
     chrid=chr${chr_num}
     java -jar ${juicer_path} dump observed VC PR-LR/data/GSE63525_K562_combined_30.hic ${chrid} ${chrid} BP 5000 ${sv_folder_1}/VC_observed.${chrid}.5000.txt
     java -jar ${juicer_path} dump oe VC PR-LR/data/GSE63525_K562_combined_30.hic ${chrid} ${chrid} BP 5000 ${sv_folder_1}/VC_oe.${chrid}.5000.txt
     paste ${sv_folder_1}/VC_observed.${chrid}.5000.txt ${sv_folder_1}/VC_oe.${chrid}.5000.txt | awk '{print $1 "\t" $2 "\t" $3 "\t" $6 "\t" $3 / $6}' > 
   ${sv_folder_1}/VC_combined.${chrid}.5000.txt
     python scripts/calc_pvalue.py ${sv_folder_1}/VC_combined.${chrid}.5000.txt ${sv_folder_2}/${chrid}_K562_prob.5000.txt
   done

   # the above commands first used juicer_tools.jar to dump HiC contact map to {sv_folder_1}, using VC normalization, on a resolution of 5000, on each {chrid}
   # oe represents observed/expected values
   # the commands then used calc_pvalue.py to fit a Poisson distribution to derive a p-value for each interaction on each {chrid}
   ```
   File `K562.chr1.res_5000.partial.txt` in the folder [resources](https://github.com/yyaoisgood2021/HUB-screening/tree/main/resources) displays the top 5000 lines of the expected output for chr1.

5. Generate fragment contact network (FCN) for each chromosome according to procedure:

	i. run all cells in `generate_eligible_coords.ipynb` for each chr and put the results in `PR-LR/eligible_coordinates`

 	ii. generate node_meta information using the next two steps:

 	* overlap eligible_coords with peaks of a feature using bedtools. The following example command implements overlapping of all eligible nodes on chr1 with feature `CTCF.narrow.rep-1.hg19.bed`. Repeat this code for all chromosomes, and for all features mentioned in the step 2.iii ~ 2.v (except black-list annotation)

	```bash
	chrid=chr1
	feature=CTCF.narrow.rep-1 # feature name should be ${feature}.hg19.bed
	
	save_folder=PR-LR/overlap_coords
	mkdir -p ${save_folder}
	
	bedtools intersect -a PR-LR/eligible_coordinates/eligible_coordinates.${chrid}.hg19.bed \
	-b PR-LR/data/${feature}.hg19.bed -wao > ${save_folder}/{feature}.{chrid}.bed
	```

	* run `scripts/prep_node_meta.py` to summarize all overlapping outputs, save results to `PR-LR/node_meta.0`
	```bash
 	python scripts/prep_node_meta.py PR-LR/eligible_coordinates PR-LR/overlap_coords PR-LR/node_meta.0
 	```
 
6. run all cells in `PR.ipynb` to build FCNs, and perform personalized PR. The PR results for all chromosomes will be generated in the folder `PR-LR/PR_scores/proc`. 

7. run all cells in `LR.ipynb` to train a LR classifer for hub essentiality. 

### Section 3 "Simple path analysis on the hub pairs"

run all cells in `FCN_analysis.ipynb` to repeat analysis in the manuscript

### Section 4 "Pseudobulk analysis using Taiji"

to repeat the analysis, following the [steps](https://github.com/yyaoisgood2021/HUB-screening/blob/main/scripts/taiji)


### Section 5 "Entropy analysis for chromatin accessibility"

to calculate entropy, you need to prepare `overlapped_data` according to the procedures below: 

1. prepare a `bedgraph file` for your signal (experiment)

   in the format of `bedgraph`: four-column bed file with chrid, start, end, signal_value (intensity. for ATAC-seq read count files, this column is universally 1)
    
   chr1 10000 10015 20
   
   chr1 10008 10025 27

   ...

   chr4 10234 10270 15

   to repeat results in the manuscript, the `bedgraph files` to use are pseudobulk cluster-specific atac-fragment.bed.gz files saved in the folder results/taiji_datasets/npcs30_pc30.3.

   Unzip them and add a uniform value of 1 to each line

   * Note the ATAC-seq was mapped to hg38, so we have to lift it over to hg19 at this step

   ```bash
   data_folder=results/taiji_datasets/npcs30_pc30.3
   save_folder_hg38=`enter save folder here`
   save_folder_hg19=`enter save folder here`
   
   fname="ess-0.atac-fragments.bed" # iterate through all files, without gz suffix

   mkdir -p ${save_folder_hg38}
   mkdir -p ${save_folder_hg19}
   cp ${data_folder}/${fname}.gz ${save_folder_hg38}
   gunzip ${save_folder_hg38}/${fname}.gz
   awk '{ print $1 "\t" $2 "\t" $3 "\t1"}' ${save_folder_hg38}/${fname} > ${save_folder_hg38}/${fname}graph

   liftOver ${save_folder_hg38}/${fname}graph hg38ToHg19.over.chain ${save_folder_hg19}/${fname}graph $save_folder_hg19/${fname}.unMapped
   ```
   
2. prepare a `coord file` for your genomic loci of interest

   in the format of `coord file`: three-column bed file of equal-sized (resolution) bins with no overlaps, I typically remove (1) centromere gaps, (2) black-list annotated regions, and (3) telomere regions

   chr1 10000 15000 

   chr1 15000 20000 

   ...

   chr4 15000 20000

   you can directly use results saved in the folder `PR-LR/eligible_coordinates` if you have run Section 2, just copy them to the folder `entropy/eligible_coords.hg19`

3. run `bedtools intersect` to map your signals to each genomic bin with the command:
   ```bash
   # this block shows an example of overlapping 'ess-0' sample's chr1

   chrid=chr1 # enter chrid here, iterate through all chromosomes in the folder
   psbulk_sample="ess-0" # enter sample id here according to meta df generated in the Taiji part, iterate through all samples in the folder, skip WT
   
   coord_file=entropy/eligible_coords.hg19/eligible_coordinates.${chrid}.hg19.bed 
   bedgraph_file=entropy/feature_bedgraphs.hg19/${psbulk_sample}.atac-fragments.bedgraph 
   save_folder=entropy/overlap_results.hg19
   save_fname_base="overlap"

   mkdir -p ${save_folder}
   bedtools intersect -wao -a ${coord_file} -b ${bedgraph_file} > \
   ${save_folder}/${save_fname_base}.${chrid}.${psbulk_sample}.txt # make sure all coordinates are under the same assembly hg19
   ```

   you will get a file for the overlapped results, and it should have 8 columns

   repeat the above codes for all chr and for all pseudobulk-specific ATAC-fragment bedgraphs

4. for general purpose, run `calc_entropy_genome1D.py` with the following command to compute entropy:
   ```bash
   overlap_data_path=`enter the path to previous step's output: the overlapping results`
   result_sv_path=`enter path to a .txt file`
   python scripts/calc_entropy_genome1D.py ${overlap_data_path} ${result_sv_path}
   ```

   to repeat results in the manuscript, you need to load all data from all pseudobulk clusters and normalize them, run the following codes

   ```bash
   overlap_result_save_folder="entropy/overlap_results.hg19"
   overlap_result_file_base="overlap" # modify this line or the python code accordingly
   result_save_folder="entropy/results" # path to the folder to save output of this code
   meta_df_path="results/taiji_results_analysis/cls-5.rep-0/meta_df.2.txt"
   # modify this line according to your Taiji results
   # also, you may want to change the path of the fragment files in `scripts/calc_entropy_1D.all_psbulk_cluster.py` accordingly

   python scripts/calc_entropy_1D.all_psbulk_cluster.py \
   ${overlap_result_save_folder} \
   ${result_save_folder} \
   ${meta_df_path} \
   ${overlap_result_file_base}
   ```
   repeat the codes for each chr
   






# HUB-screening

This is the code for the manuscript: "Distinct 3D contacts and phenotypic consequences of adjacent non-coding loci in the epigenetically quiescent regions" to implement to following sections:

### Section 1 "Constructing and training of logistic regression to facilitate candidate selection"

To run `CNN-AE-LR`, you need to prepare the sequence data for each chromosome (chr1 ~ chr22, and chrX) (hg19): download from [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/), and save files to CNN-AE-LR/data/hg19/

Then, you can run codes in scripts/CNN-AE-LR.ipynb


### Section 2 "Classifying hub essentiality with sequence and epigenetic features"

1. To run codes in this section, you need to install [juicer](https://github.com/aidenlab/juicer) (v1.6) and [bedtools](https://bedtools.readthedocs.io/en/latest/index.html).

2. Download the following data and put them in the designated folder `PR-LR/data`, the expected folder tree structure is displayed in [folder_tree.txt](https://github.com/yyaoisgood2021/HUB-screening/blob/main/folder_tree.txt):
   
	i. supplementary table `sciadv.abi6020_table_s11.xlsx` from [Ding et al. study](https://www.science.org/doi/10.1126/sciadv.abi6020) (DOI: 10.1126/sciadv.abi6020)
	
 	ii. Hi-C data (`GSE63525_K562_combined_30.hic`) from [Rao's work](https://www.cell.com/fulltext/S0092-8674(14)01497-4) (DOI: 10.1016/j.cell.2014.11.021)


3. then generate the `{chrid}_K562_prob.5000.txt` file for each chromosome using the following bash commands

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
  paste ${sv_folder_1}/VC_observed.${chrid}.5000.txt ${sv_folder_1}/VC_oe.${chrid}.5000.txt | awk '{print $1 "\t" $2 "\t" $3 "\t" $6 "\t" $3 / $6}' > ${sv_folder_1}/VC_combined.${chrid}.5000.txt
  python scripts/calc_pvalue.py ${sv_folder_1}/VC_combined.${chrid}.5000.txt ${sv_folder_2}/${chrid}_K562_prob.5000.txt
done

# the above commands first used juicer_tools.jar to dump HiC contact map to {sv_folder_1}, using VC normalization, on a resolution of 5000, on each {chrid}
# oe represents observed/expected values
# the commands then used calc_pvalue.py to fit a Poisson distribution to derive a p-value for each interaction on each {chrid}
```
3. then generate fragment contact network (FCN) for each chromosome according to procedure:

	i. run `generate_eligible_coords` for each chr and put the results in `coord_save_folder`

 	ii. run the following bash commands to generate node_meta, for each chromosome, for each feature

	```bash
	bedtools intersect -a path/to/eligible_coords.${chrid}.bed -b path/to/{feature}.bed_peak_file -wao > path/to/overlap_save_folder/overlap.{feature}.{chrid}.bed

 	path/to/python scripts/prep_node_meta.py path/to/coord_save_folder path/to/overlap_save_folder path/to/node_meta_save_folder
	```

 	iii. run `PR` to build FCNs, and perform PR 

	iv. run `LR` to train a LR classifer for hub essentiality (you may want to modify the codes)

### Section 3 "Simple path analysis on the hub pairs"

run `FCN_analysis` to repeat analysis in the manuscript

### Section 4 "Entropy analysis for chromatin accessibility"

to calculate entropy, you need to prepare `overlapped_data` according to the procedures below: 

1. prepare a `bedgraph file` for your signal (experiment)
format of `bedgraph file`: 
    four-column bed file with chrid, start, end, signal_value (intensity, for read count file, this column is universally 1)
    chr1 10000 10015 20
    chr1 10008 10025 27
    ...
    chr4 10234 10270 15

2. prepare a `coord file` for your genomic loci of interest
format of `coord file`: 
    three-column bed file of equal-sized (resolution) bins with no overlaps, I typically remove centromere and telomere regions
    chr1 10000 15000 
    chr1 15000 20000 
    ...
    chr4 15000 20000 

3. run bedtools intersect to map your signals to each genomic bin with the command:
    ```bash
    bedtools intersect -wao -a `coord file` -b `bedgraph file` > `overlap file save path`
    ```
	you will get a file of the overlapped results, and it should have 8 columns

4. run `calc_entropy_genome1D.py` with the following command:
    ```bash
    python scripts/calc_entropy_genome1D.py `overlap file save path` `entropy result save path`
    ```

### Section 5 "Pseudobulk analysis using Taiji"

to repeat the analysis, following the [steps](https://github.com/yyaoisgood2021/HUB-screening/blob/main/scripts/taiji)







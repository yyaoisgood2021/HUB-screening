# HUB-screening

This is the code for the manuscript: "Distinct 3D contacts and phenotypic consequences of adjacent non-coding loci in the epigenetically quiescent regions" to implement to following sections:

### Section 1 "Constructing and training of logistic regression to facilitate candidate selection"

To run `CNN-AE-LR`, you need to prepare the sequence data for each chromosome (hg19): download from [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/), and save files to resources/hg19/

Then, you can run codes in scripts/CNN-AE-LR.ipynb


### Section 2 "Classifying hub essentiality with sequence and epigenetic features"

1. you need to download the following data:
   
	i. `TableS11` from [Ding et al. study](https://www.science.org/doi/10.1126/sciadv.abi6020) (DOI: 10.1126/sciadv.abi6020)
	
 	ii. `TableS4` from this manuscript, data is copied to [resources](https://github.com/yyaoisgood2021/HUB-screening/tree/main/resources)
	
 	iii. Hi-C data (`inter_30.hic`) from [Rao's work](https://www.cell.com/fulltext/S0092-8674(14)01497-4) (DOI: 10.1016/j.cell.2014.11.021)

	iv. Download the feature data in the bed format, refer to `TableS12` and `download_chip_dt.txt` from this manuscript in [resources](https://github.com/yyaoisgood2021/HUB-screening/tree/main/resources)

2. then generate the `{chrid}_K562_prob.5000.txt` file for each chromosome using the following bash commands:

```bash
chrid=`enter chr1 ~ chrX here`

java -jar path/to/juicer_tools.jar dump observed VC path/to/inter_30.hic ${chrid} ${chrid} BP 5000 path/to/save_folder/VC_observed.${chrid}.5000.txt

java -jar path/to/juicer_tools.jar dump oe VC path/to/inter_30.hic ${chrid} ${chrid} BP 5000 path/to/save_folder/VC_oe.${chrid}.5000.txt

paste path/to/save_folder/VC_observed.${chrid}.5000.txt path/to/save_folder/VC_oe.${chrid}.5000.txt | awk '{print $1 "\t" $2 "\t" $3 "\t" $6 "\t" $3 / $6}' > path/to/save_folder/VC_combined.${chrid}.5000.txt

path/to/python scripts/calc_pvalue.py path/to/save_folder/VC_combined.${chrid}.5000.txt path/to/save_folder/${chrid}_K562_prob.5000.txt

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






# HUB-screening

This is the codes for the manuscript: "Distinct 3D contacts and phenotypic consequences of adjacent non-coding loci in the epigenetically quiescent regions" to implement to following sections:

"Constructing and training of logistic regression to facilitate candidate selection"

To run the code, you need to prepare the sequence data for each chromosome (hg19): download from 
"https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/", and save files to resource/hg19


"Classifying hub essentiality with sequence and epigenetic features"

1. you need to download the following data:
	i. TableS11 from Ding et al. study (DOI: 10.1126/sciadv.abi6020)
	ii. TableS4 from this manuscript
	iii. Hi-C data (inter_30.hic) from (DOI: 10.1016/j.cell.2014.11.021)
	iv. Download the feature data in the bed format, refer to TableS12 from this manuscript 

2. then generate the {inter_prob} file for each chromosome using the following commands:

remember to replace {chrid} with chr1 to chrX

`java -jar path/to/juicer_tools.jar dump observed VC path/to/inter_30.hic {chrid} {chrid} BP 5000 path/to/save_folder/VC_observed.{chrid}.5000.txt`

`java -jar path/to/juicer_tools.jar dump oe VC path/to/inter_30.hic {chrid} {chrid} BP 5000 path/to/save_folder/VC_oe.{chrid}.5000.txt`

`paste path/to/save_folder/VC_observed.{chrid}.5000.txt path/to/save_folder/VC_oe.{chrid}.5000.txt | awk '{print $1 "\t" $2 "\t" $3 "\t" $6 "\t" $3 / $6}' > path/to/save_folder/VC_combined.{chrid}.5000.txt`

`path/to/python script/calc_pvalue.py path/to/save_folder/VC_combined.{chrid}.5000.txt path/to/save_folder/{chrid}_K562_prob.5000.txt` 


3. then run the code to perform PageRank and LR



"Simple path analysis on the hub pairs"


"Entropy analysis for chromatin accessibility"








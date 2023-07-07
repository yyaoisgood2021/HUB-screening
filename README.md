# HUB-screening

This is the codes for the manuscript: "Distinct 3D contacts and phenotypic consequences of adjacent non-coding loci in the epigenetically quiescent regions" to implement to following sections:

### "Constructing and training of logistic regression to facilitate candidate selection"

To run the code, you need to prepare the sequence data for each chromosome (hg19): download from [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/), and save files to resources/hg19


### "Classifying hub essentiality with sequence and epigenetic features"

1. you need to download the following data:
   
	i. TableS11 from [Ding et al. study](https://www.science.org/doi/10.1126/sciadv.abi6020) (DOI: 10.1126/sciadv.abi6020)
	
 	ii. TableS4 from this manuscript, data is copied to [resources](https://github.com/yyaoisgood2021/HUB-screening/tree/main/resources)
	
 	iii. Hi-C data (inter_30.hic) from [Rao's work](https://www.cell.com/fulltext/S0092-8674(14)01497-4) (DOI: 10.1016/j.cell.2014.11.021)

	iv. Download the feature data in the bed format, refer to TableS12 from this manuscript in [resources](https://github.com/yyaoisgood2021/HUB-screening/tree/main/resources)

2. then generate the {inter_prob} file for each chromosome using the following commands:

```bash
chrid=`enter chr1 ~ chrX here`

java -jar path/to/juicer_tools.jar dump observed VC path/to/inter_30.hic ${chrid} ${chrid} BP 5000 path/to/save_folder/VC_observed.${chrid}.5000.txt

java -jar path/to/juicer_tools.jar dump oe VC path/to/inter_30.hic ${chrid} ${chrid} BP 5000 path/to/save_folder/VC_oe.${chrid}.5000.txt

paste path/to/save_folder/VC_observed.${chrid}.5000.txt path/to/save_folder/VC_oe.${chrid}.5000.txt | awk '{print $1 "\t" $2 "\t" $3 "\t" $6 "\t" $3 / $6}' > path/to/save_folder/VC_combined.${chrid}.5000.txt

path/to/python scripts/calc_pvalue.py path/to/save_folder/VC_combined.${chrid}.5000.txt path/to/save_folder/${chrid}_K562_prob.5000.txt

```
3. then generate fragment contact network (FCN) for each chromosome according to procedure:



### "Simple path analysis on the hub pairs"

run the code to build FCN




### "Entropy analysis for chromatin accessibility"








# STAR protocol for cooperative binding analysis using dSMF data

This protocol is derived from [Rao et al., 2021](https://pubmed.ncbi.nlm.nih.gov/33705711/).  

## Before you begin

Data for demo is included in this github repository, but to visualize at your
sites of interest, please download the sequencing data, and keep them in
`data_from_geo/`. Here is the list of URLs for the sequencing data. 
```
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR313/006/SRR3133326/SRR3133326_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR313/006/SRR3133326/SRR3133326_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR313/007/SRR3133327/SRR3133327_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR313/007/SRR3133327/SRR3133327_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR313/008/SRR3133328/SRR3133328_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR313/008/SRR3133328/SRR3133328_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR313/009/SRR3133329/SRR3133329_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR313/009/SRR3133329/SRR3133329_2.fastq.gz
```

## Download reference genome
The reference genome fasta file can be downloaded from [here](https://hgdownload.cse.ucsc.edu/goldenPath/dm3/bigZips/dm3.fa.gz). Place this fasta file in the `ref_genome/dm3` directory.

## To reproduce figure 5C,D

Please run the following single command. 

```
snakemake --snakefile cooperative_binding_analysis.smk plots/single_binding/suppressed_merged_demo_S2_to_example_spanning_lf_15_rf_15_extended_left_150_right_150_roi_peak_229.fp.pdf plots/single_binding/suppressed_merged_demo_S2_to_example_spanning_lf_15_rf_15_extended_left_150_right_150_roi_peak_229.methylation.pdf --configfile configs/config.yaml
```


## To reproduce figure 6C
```
snakemake  --snakefile cooperative_binding_analysis.smk plots/cobinding_bedpe/suppressed_merged_demo_S2_to_example_cobinding_lf_15_rf_15_extended_left_300_right_300_roi_peak_110_4_and_peak_110_6.fp.pdf plots/cobinding_bedpe/suppressed_merged_demo_S2_to_example_cobinding_lf_15_rf_15_extended_left_300_right_300_roi_peak_110_4_and_peak_110_6.methylation.pdf --configfile configs/config.yaml
```

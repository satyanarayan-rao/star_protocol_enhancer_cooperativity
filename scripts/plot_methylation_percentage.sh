#!/bin/bash 
# $1 : occluded edges open
# $2 : occluded edges closed 
# $3 : output png 
# head -3 occluded_edges/suppressed_merged_S2_to_mnase_peaks_open_enh_lf_15_rf_15_occluded.tsv | column -t
# peak_2105_1`SRR3133329.24391028_24391028/1_overlapping`83~163  289  24  0.5   25  260  2  25
# peak_2105_1`SRR3133329.6493044_6493044/1_overlapping`83~163    293  25  0.56  12  272  7  12
# peak_2105_1`SRR3133326.13323124_13323124/1_overlapping`99~147  294  25  0.0   0   0    0  0

cut -d'`' -f2- $1 | sort | uniq > ${1%.tmp}.tmp.tsv 
cut -d'`' -f2- $2 | sort | uniq > ${2%.tmp}.tmp.tsv 

# Rscript arguments
# args [1] : reads mapped to open enhancer
# args [2] : reads mapped to closed enhancer
# args [3] : percetage methylation field
# args [4] : bin width (5) 

Rscript scripts/plot_percentage_methylation_on_reads.R ${1%.tmp}.tmp.tsv ${2%.tmp}.tmp.tsv 4  5  ${3%.png}.pdf ${3}

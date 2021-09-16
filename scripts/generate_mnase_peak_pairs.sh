#!/bin/bash
# $1 : mnase peaks in open and closed enhancers
# $2 : dsmf reads mapped to $1
# $3 : genome_size file
# $4 : distance threshold
# $5 : span from peak
# $6 : output mnase peak pairs
# $7 : output dsmf reads mapped mnase peak pairs
# head -2 input_bed/mnase_peaks_in_open_and_closed_enhancers.bed
#chr2L	47389	47390	peak_877_1	closed	+
#chr2L	107902	107903	peak_2105_9	open	+
awk '$5=="open"' $1 | sort -k1,1  -k2,2n -k3,3n  > ${1}.peaks_in_open_enh.bed
python scripts/get_all_mnase_peak_pair_bed.py ${1}.peaks_in_open_enh.bed 30 > ${1}.all_peak_pairs.tsv

bedtools slop -b $5 -g $3 -i ${1}.all_peak_pairs.tsv > ${1}.all_peak_pairs.${5}.bed

zcat $2 | cut -f7- | sort -S4G | uniq | sort -S4G --parallel=4 -k1,1 -k2,2n -k3,3n  > ${2}.dsmf.bed 
bedtools intersect -a ${1}.all_peak_pairs.${5}.bed -b ${2}.dsmf.bed -wa -wb -f 1 -sorted > ${7%.gz} 

python scripts/count_intersected_reads_per_peak_unique.py ${7%.gz} > ${6}.dsmf_counts.tsv 

awk '{if ($(NF-1)>=15 && $NF>=15){print $0}}' ${6}.dsmf_counts.tsv > $6

cat ${7%.gz} | gzip - > $7 
# cleanup
rm ${1}.peaks_in_open_enh.bed  ${1}.all_peak_pairs.tsv   ${1}.all_peak_pairs.${5}.bed ${2}.dsmf.bed ${6}.dsmf_counts.tsv 

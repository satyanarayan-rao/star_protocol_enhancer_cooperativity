#!/bin/bash

# $1: overlapping or adjacent fragment file
# $2: regions of interest
# $3: output file
# $4: left flank
# $5: right flank

#sort -k1,1 -k2,2n -k3,3n ${2} > ${2}.tmp.bed
#sort -k1,1 -k2,2n -k3,3n ${2}  | awk '{OFS="\t"}{$2=$2-1; $3=$3+1; print $0}' | paste - ${2}.tmp.bed > ${2}.sorted.bed

# head input_bed/example_cobinding.bedpe
#chr2L	19155158	19155159	chr2L	19155265	19155266	peak_110_4_and_peak_110_6	.	.	.
# know the number of columns in bedpe
ncols=`head -1 $2 | awk -F'\t' '{print NF}'` 

awk '{print $1"\t"$2-left"\t"$6+right"\t"$0}' left=$4 right=$5 $2 | sort -k1,1 -k2,2n -k3,3n | bedtools intersect -a - -b $1 -wa -wb -f 1 -sorted | gzip - >  ${2}.intersect.bed.gz
added_cols=`echo "$ncols + 3" | bc`
zcat ${2}.intersect.bed.gz | python scripts/compress_bedpe_fields.py $added_cols |  python scripts/add_edges_to_footprints_bedpe.py 8 | gzip - > $3 

# cleanup 

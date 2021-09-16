#!/bin/bash

# $1: overlapping or adjacent fragment file
# $2: regions of interest
# $3: output file
sort -k1,1 -k2,2n -k3,3n ${2} > ${2}.tmp.bed
sort -k1,1 -k2,2n -k3,3n ${2}  | awk '{OFS="\t"}{$2=$2-1; $3=$3+1; print $0}' | paste - ${2}.tmp.bed > ${2}.sorted.bed
bedtools intersect -a ${2}.sorted.bed -b $1 -wa -wb -sorted | cut -f7-  | python scripts/add_edges_to_footprints.py | gzip - > $3 

# cleanup 

rm ${2}.sorted.bed ${2}.tmp.bed 

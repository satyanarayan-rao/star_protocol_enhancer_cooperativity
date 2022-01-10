#!/bin/bash
# $1: input list of bams "a.bam b.bam c.bam"
# $2: merged bam file
# $3: wildcard
for i in `echo $1`
do 
    echo $i 
done > ${2}.merge_list.tsv
bamtools merge -list ${2}.merge_list.tsv -out $2
# cleanup
rm ${2}.merge_list.tsv

#!/bin/bash
# $1 : read_level_methylation_status
# $2 : overlapping_or_adjacent

if [[ $OSTYPE == 'darwin'* ]]; then 
    gzcat $1 | egrep "_overlapping|_adjacent" | sort -S8G --parallel=8  -k1,1 -k2,2n -k3,3n | gzip - > $2
else 
    zcat $1 | egrep "_overlapping|_adjacent" | sort -S8G --parallel=8  -k1,1 -k2,2n -k3,3n | gzip - > $2
fi

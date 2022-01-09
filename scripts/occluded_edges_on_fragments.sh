#!/bin/bash
# $1 : fragments_covering_flank
# $2 : occluded_edges
# $3 : occluded_pkl

if [[ $OSTYPE == 'darwin'* ]]; then
    gzcat $1 |python scripts/identify_occluded_edges.py $2 $3
else
    zcat $1 |python scripts/identify_occluded_edges.py $2 $3
fi


#!/bin/bash
# $1 : fragments covering flank
# $2 : footprint dictionary

if [[ $OSTYPE == 'darwin'* ]]; then
        gzcat  $1 | python scripts/build_footprint_dict.py $2
else
        zcat  $1 | python scripts/build_footprint_dict.py $2
     
fi


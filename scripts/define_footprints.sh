#!/bin/bash
# $1 : overlapping_or_adjacent input file 
# $2 : wobble gap
# $3 : min_fp_len
# $4 : output file

if [[ $OSTYPE == 'darwin'* ]]; then
        gzcat  $1 | python scripts/call_footprints.py | python scripts/fix_wobble.py $2 | python scripts/keep_all_footprint_ge_th.py $3 | gzip - > $4 
else
        zcat  $1 | python scripts/call_footprints.py | python scripts/fix_wobble.py $2 | python scripts/keep_all_footprint_ge_th.py $3 | gzip - > $4 
     
fi


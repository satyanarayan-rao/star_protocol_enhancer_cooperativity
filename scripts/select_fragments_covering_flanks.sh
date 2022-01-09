#!bin/bash
# $1 : mapped_to_regions
# $2 : left flank
# $3 : right flank
# $4 : output: fragments covering flank
if [[ $OSTYPE == 'darwin'* ]]; then
    gzcat $1 |  python scripts/fragments_spanning_flanks.py  $2 $3 | gzip - > $4
else
    zcat $1 |  python scripts/fragments_spanning_flanks.py  $2 $3 | gzip - > $4
fi

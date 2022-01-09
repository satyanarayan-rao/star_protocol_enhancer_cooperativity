#!/bin/bash
# ${1} : fragment covering flank
# ${2} : left flank
# ${3} : right flank
# ${4} : methylation matrix
# ${5} : methylation matrix real footprint length pkl
# ${6} : methylation matrix footprint vec pkl
# ${7} : strand_agnostic_footprint_vec_pkl 
# ${8} : strand_agnostic_footprint_matrix
# ${9} : footprint_length_at_bp_resolution_pkl 
# ${10} : footprint_length_at_bp_resolution_tsv 

if [[ $OSTYPE == 'darwin'* ]]; then
    gzcat ${1} | python scripts/build_flank_methylation_matrix.py ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10}
else
    zcat ${1} | python scripts/build_flank_methylation_matrix.py ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10}
fi

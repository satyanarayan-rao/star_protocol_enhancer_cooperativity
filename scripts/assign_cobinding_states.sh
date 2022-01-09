#!/bin/bash
# $1 : extended fragments bedgz
# $2 : extended fragments verbose bedgz : have methylation and sequence information
# $3 : left flank 
# $4 : right flank 
# $5 : lextend
# $6 : rextend
# $7 : output footprint file with binding states assigned : bedgz
# $8 : output verbose file with binding states: bedgz
# $9 : output verbose file with binding states spanning 150 bp from the primary peak

if [[ $OSTYPE == 'darwin'* ]]; then
    gzcat $2  > ${2%.gz}
    gzcat $1 | python scripts/assign_cobinding_states.py ${2%.gz} ${3} ${4} ${5} ${6} ${7%.gz} ${8%.gz} ${9%.gz}
else
    zcat $2  > ${2%.gz}
    zcat $1 | python scripts/assign_cobinding_states.py ${2%.gz} ${3} ${4} ${5} ${6} ${7%.gz} ${8%.gz} ${9%.gz}
fi 

cat ${7%.gz} | gzip - > $7
cat ${8%.gz} | gzip - > $8
cat ${9%.gz} | gzip - > $9

# clean up
rm ${7%.gz} ${8%.gz} ${9%.gz}

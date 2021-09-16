#!/bin/bash
# $1 : bedgz input file
# $2 : footprint dict
# $3 : roi metadata pkl
# $4 : lflank 
# $5 : rflank 
# $6 : lextend 
# $7 : rextend
# $8 : output bedgz 
# $9 : verbose bedgz 

echo "zcat $1 | python scripts/extend_fragments_for_cobinding.py ${2} ${3} ${4} ${5} ${6} ${7} ${8%.gz} ${9%.gz}"
zcat $1 | python scripts/extend_fragments_for_cobinding.py ${2} ${3} ${4} ${5} ${6} ${7} ${8%.gz} ${9%.gz} 
cat ${8%.gz} | gzip - > $8 
cat ${9%.gz} | gzip - > $9  

# cleanup
rm ${8%.gz} ${9%.gz}

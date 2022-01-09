#!/bin/bash
# $1 : bedgz input file
# $2 : lflank 
# $3 : rflank 
# $4 : lextend 
# $5 : rextend
# $6 : output bedgz 
# $7 : verbose bedgz 

#echo "zcat $1 | python scripts/extend_fragments_for_cobinding.py ${2} ${3} ${4} ${5} ${6} ${7} ${8%.gz} ${9%.gz}"
echo "zcat $1 | python scripts/extend_fragments_for_cobinding.py ${2} ${3} ${4} ${5} ${6%.gz} ${7%.gz}" 

if [[ $OSTYPE == 'darwin'* ]]; then
    gzcat $1 | python scripts/extend_fragments_for_cobinding_bedpe.py ${2} ${3} ${4} ${5} ${6%.gz} ${7%.gz} 
else
    zcat $1 | python scripts/extend_fragments_for_cobinding_bedpe.py ${2} ${3} ${4} ${5} ${6%.gz} ${7%.gz} 
fi
cat ${6%.gz} | gzip - > $6
cat ${7%.gz} | gzip - > $7  

# cleanup
rm ${6%.gz} ${7%.gz}

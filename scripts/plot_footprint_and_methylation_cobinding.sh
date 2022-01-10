#!/bin/bash
# $1 : footprint tsv 
# $2 : methylation tsv
# $3 : mnase data
# $4 : base mnase params 
# $5 : base footprint params
# $6 : foootprint pdf
# $7 : methylation pdf   
# $8 : footprint gplt
# $9 : methylation gplt
# ${10} : lexetend
# ${11} : rextend 
# ${12} : lflank
# ${13} : rflank 
# ${14} : methylation params gplt
# ${15} : footprint matrix
# ${16} : methylation matrix
# ${17} : roi id
# ${18} : input bedpe -> to calculate the distance between peaks

cat $1 | python scripts/select_valid_states.py 9 | awk '{$1=""; print $0}' > ${15}
cat $2 | python scripts/select_valid_states.py 9 | awk '{$1=""; print $0}' > ${16}
tot_lines=`wc -l ${15} | awk '{print $1}'`


cat $1 | python scripts/select_valid_states.py 9 | awk -F '[\t#]' '{print $2}' | sort -n | uniq -c > ${1}.lc.tsv 

left_coor=`expr -${10} - 50`
label_coor=`expr -${10} - 10`
### prepare mnase data
if [[ $OSTYPE == 'darwin'* ]]; then
    gzcat $3 | grep -w "${17}" | cut -f2- | tr '\t' '\n' | awk '{print NR-1000"\t"$1}'   > ${15}.mnase.tsv
else
    zcat $3 | grep -w "${17}" | cut -f2- | tr '\t' '\n' | awk '{print NR-1000"\t"$1}'   > ${15}.mnase.tsv
fi

## Get peak pair distance
d=`grep -w ${17} ${18}| awk '{print $5 - $2 + 1}'`
d_l=`echo $d - 15 | bc`
d_r=`echo $d + 15 | bc`


############  footprint gnuplot ############### 
cat << EOT > $8 
reset
set terminal postscript size 4.5in,7.5in font 'Arial, 15'
#set terminal cairolatex  size 4.5in,7.5in font 'Arial, 15'
#set terminal pdf  size 4.5in,7.5in font 'Arial, 15'
set tmargin 1 
set bmargin 0 
set lmargin 3 
set rmargin 3 
set xra [-${10} : ${11}]  
unset xtics 
unset ytics
unset colorbox
set style rect fc lt -1 fs solid 0.1 noborder
set obj rect from -${12}, graph 0 to ${13}, graph 1   
set obj rect from ${d_l}, graph 0 to ${d_r}, graph 1   
set output '${6%.pdf}.eps'
#set output '${6}'
#set output '${15}'
set multiplot layout 12, 1
#plot '${3}' u 1:2 w l lc rgb "#31a354" lw 2  notitle
set key font 'Arial, 10' left top
plot '${15}.mnase.tsv' u 1:2 w l lc rgb "#31a354" lw 2 title 'MNase'
#plot '${15}.mnase.tsv' u 1:2 with filledcurve title 'MNase' lc rgb '#31a354' 
load '${5}'
set yra [${tot_lines} : 0 ] 
EOT

c=1
s=0
prev_val=0
echo "set arrow 20 from ${left_coor}, 0 to ${11}, 0 front lc rgb '#636363' lw 2 nohead" >> $8
while read cnt cl 
do
    state_lab=`python scripts/get_state_label_at_pair_of_sites.py ${cl}`
    s=`expr $s + $cnt`
    mid_point=`echo "($s + ${prev_val})/2" | bc`
    echo "set arrow $c from ${left_coor} , $s to ${11}, $s front lc rgb '#636363' lw 2 nohead" >> $8
    echo "set label $c 'n = $cnt (${state_lab})' at ${label_coor}, ${mid_point} right front font 'Arial, 10' " >> $8
    c=`expr $c + 1` 
    prev_val=$s
done < ${1}.lc.tsv
echo "plot '${15}' mat u (\$1 - ${10}):2:3 w image notitle" >> $8
gnuplot $8 
ps2pdf ${6%.pdf}.eps $6 #${6%.pdf}.tmp.pdf
rm ${6%.pdf}.eps 
#pdfcrop ${6%.pdf}.tmp.pdf $6 

#ps2pdf ${15} ${6}
############ methylation gnuplot ############## 

cat << EOT > $9
reset
set terminal postscript size 4.5in, 7.5in font 'Arial, 15'
#set terminal cairolatex eps size 4.5in,7.5in font 'Arial, 15'
#set terminal pdf size 4.5in,7.5in font 'Arial, 15'
set tmargin 1 
set bmargin 0 
set lmargin 3 
set rmargin 3 
set xra [-${10} : ${11}] 
unset xtics 
unset ytics
unset colorbox
set style rect fc lt -1 fs solid 0.1 noborder
set obj rect from -${12}, graph 0 to ${13}, graph 1    
set obj rect from ${d_l}, graph 0 to ${d_r}, graph 1   
set output '${7%.pdf}.eps'
#set output '${7}'
#set output '${16}'
set multiplot layout 12, 1
#plot '${3}' u 1:2 w l lc rgb "#31a354" lw 2  notitle 
set key font 'Arial, 10' left top
plot '${15}.mnase.tsv' u 1:2 w l lc rgb "#31a354" lw 2  title 'MNase'
plot '${15}.mnase.tsv' u 1:2 with filledcurve notitle 
load '${14}'
set yra [$tot_lines: 0 ]
EOT
c=1
s=0
prev_val=0

echo "set arrow 20  from ${left_coor}, 0 to ${11}, 0 front lc rgb '#636363' lw 2 nohead" >> $9
while read cnt cl 
do
    state_lab=`python scripts/get_state_label_at_pair_of_sites.py ${cl}`
    s=`expr $s + $cnt`
    mid_point=`echo "($s + ${prev_val})/2" | bc`
    echo "set arrow $c from ${left_coor}, $s to ${11}, $s front lc rgb '#636363' lw 2 nohead" >> $9
    echo "set label $c 'n = $cnt (${state_lab})' at ${label_coor}, ${mid_point} right front font 'Arial, 10'" >> $9
    c=`expr $c + 1` 
    prev_val=$s
done < ${1}.lc.tsv

echo "plot '${16}' mat u (\$1 - ${10}):2:3 w image notitle" >> $9
gnuplot $9
ps2pdf ${7%.pdf}.eps $7 #${7%.pdf}.tmp.pdf
rm ${7%.pdf}.eps 
#pdfcrop ${7%.pdf}.tmp.pdf $7
#ps2pdf ${16} ${7}


#rm ${1}.tmp ${2}.tmp

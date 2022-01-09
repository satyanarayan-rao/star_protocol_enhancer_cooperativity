#!/bin/bash

# $1 : verbos file with edges as E
# $2 : output footprint ordered_footprint vec
# $3 : output ordered methylation vector
# $4 : roi id : ex: peak_229

## Here are all the steps to get to footprint vector ## Here are all the steps to get to footprint 
roi_id=$4
if [[ $OSTYPE == 'darwin'* ]]; then
    gzcat $1 | grep ${roi_id} > ${2}.tmp.tsv
else
    zcat $1 | grep ${roi_id} > ${2}.tmp.tsv
fi
cat ${2}.tmp.tsv | awk 'NR%3==1' | grep "99~147" > ${2}.fp.99.tsv 
cat ${2}.tmp.tsv | awk 'NR%3==1' | grep "83~163" > ${2}.fp.83.tsv 
cat ${2}.tmp.tsv | awk 'NR%3==2' | grep "99~147" > ${3}.mvec.99.tsv 
cat ${2}.tmp.tsv | awk 'NR%3==2' | grep "83~163" > ${3}.mvec.83.tsv 


#python scripts/order_footprints_by_length_and_direction_left_ordered.py ${2}.fp.99.tsv ${2}.o_fp.99.tsv
#python scripts/order_footprints_by_length_and_direction_left_ordered.py ${2}.fp.83.tsv ${2}.o_fp.83.tsv  
#python scripts/order_footprints_by_length_and_direction_left_ordered.py ${3}.mvec.99.tsv ${3}.o_mvec.99.tsv
#python scripts/order_footprints_by_length_and_direction_left_ordered.py ${3}.mvec.83.tsv ${3}.o_mvec.83.tsv

python scripts/get_largest_footprint_location.py ${2}.fp.99.tsv  > ${2}.fp_to_sort.99.tsv 
python scripts/get_largest_footprint_location.py ${2}.fp.83.tsv  > ${2}.fp_to_sort.83.tsv

awk '{print $NF}' ${3}.mvec.99.tsv | paste ${2}.fp_to_sort.99.tsv - | awk -F '#' '{print $1"\t"$2}' | sort -k2,2n -k3,3g | awk '{print $1"#"$2"\t"$5"\t"$6}' > ${2}.fp_and_mvec.99.tsv

awk '{print $1"\t"$2}' ${2}.fp_and_mvec.99.tsv > ${2}.length_ordered_fp.99.tsv 
awk '{print $1"\t"$3}' ${2}.fp_and_mvec.99.tsv > ${2}.length_ordered_mvec.99.tsv

awk '{print $NF}' ${3}.mvec.83.tsv | paste ${2}.fp_to_sort.83.tsv - | awk -F '#' '{print $1"\t"$2}' | sort -k2,2n -k3,3g | awk '{print $1"#"$2"\t"$5"\t"$6}' > ${2}.fp_and_mvec.83.tsv

awk '{print $1"\t"$2}' ${2}.fp_and_mvec.83.tsv > ${2}.length_ordered_fp.83.tsv 
awk '{print $1"\t"$3}' ${2}.fp_and_mvec.83.tsv > ${2}.length_ordered_mvec.83.tsv


# Now find out what are the unique category we have by merging ${2}.fp_and_mvec.99.tsv and ${2}.fp_and_mvec.83.tsv and looking at the number after `#`

cat ${2}.fp_and_mvec.99.tsv ${2}.fp_and_mvec.83.tsv | awk -F'[#\t]' '{print $2}' | sort -n | uniq > ${2}.states.tsv 

for i in `cat ${2}.states.tsv`
do
    cat ${2}.length_ordered_fp.83.tsv | python scripts/select_by_id.py ${i} 
    cat ${2}.length_ordered_fp.99.tsv | python scripts/select_by_id.py ${i}
done > ${2}.before_numeric.tsv

for i in `cat ${2}.states.tsv`
do
    cat ${2}.length_ordered_mvec.83.tsv | python scripts/select_by_id.py ${i}
    cat ${2}.length_ordered_mvec.99.tsv | python scripts/select_by_id.py ${i}
done > ${3}.before_numeric.tsv # note the change from 2 to 3

python scripts/footprint_dot_to_digit_vec.py ${2}.before_numeric.tsv $2
python scripts/methylation_dot_to_digit_vec.py ${3}.before_numeric.tsv $3

# cleanup
rm ${2}.length_ordered_mvec.99.tsv ${2}.length_ordered_mvec.83.tsv ${2}.length_ordered_fp.83.tsv ${2}.length_ordered_fp.99.tsv ${2}.before_numeric.tsv ${3}.before_numeric.tsv ${2}.states.tsv ${2}.fp_and_mvec.83.tsv ${2}.fp_and_mvec.99.tsv ${2}.fp_to_sort.99.tsv ${2}.fp_to_sort.83.tsv ${2}.fp.99.tsv ${2}.fp.83.tsv ${3}.mvec.99.tsv ${3}.mvec.83.tsv

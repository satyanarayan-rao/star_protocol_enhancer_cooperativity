#!/bin/bash
# first collect reads from all footprint tsv files
find ./fragments_ordered_single_binding/ -name "*fp.tsv" | xargs -i cat {} | awk -F'`' '{print $2}' | awk -F'_' '{print "@"$1}' > tmp/single.tsv
find ./fragments_ordered_cobinding_bedpe/ ./fragments_ordered_cobinding/ -name "*fp.tsv" | xargs -i cat {} | awk -F'`' '{print $9}' | awk -F'_' '{print "@"$1}' > tmp/cobinding.tsv
cat tmp/single.tsv tmp/cobinding.tsv | sort | uniq  > tmp/all_selected_reads.tsv

for fq in SRR3133326 SRR3133327 SRR3133328 SRR3133329
do
    grep $fq tmp/all_selected_reads.tsv  > tmp/to_search_in_${fq}.tsv 
    echo -e "zcat ~/workplace/projects/enhancer-cooperativity_smk/raw_data/${fq}_1.fastq.gz | python \$N/select_for_fastq_read_id.py tmp/to_search_in_${fq}.tsv | gzip - > tmp/selected_reads_for_demo_data/demo_${fq}_1.fastq.gz"
    echo -e "zcat ~/workplace/projects/enhancer-cooperativity_smk/raw_data/${fq}_2.fastq.gz | python \$N/select_for_fastq_read_id.py tmp/to_search_in_${fq}.tsv | gzip - > tmp/selected_reads_for_demo_data/demo_${fq}_2.fastq.gz"
done 

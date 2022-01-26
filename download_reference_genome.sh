#!/bin/bash
cwd=`pwd`
mkdir -p ref_genome/dm3
cd ref_genome/dm3 
echo "Downloading dm3 reference genome...\n"
if [[ $OSTYPE == 'darwin'* ]]; then
    curl -JLO https://hgdownload.cse.ucsc.edu/goldenPath/dm3/bigZips/dm3.fa.gz
else
    wget https://hgdownload.cse.ucsc.edu/goldenPath/dm3/bigZips/dm3.fa.gz
fi
gunzip dm3.fa.gz
echo "Download and unzipping finished!"
cd $cwd

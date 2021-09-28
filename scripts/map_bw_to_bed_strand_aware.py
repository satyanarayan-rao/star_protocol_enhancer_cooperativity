import sys 
import os
import pyBigWig
import gzip
import math

bw_file = sys.argv[1]
bed_file = sys.argv[2]
out_file = sys.argv[3]
strand_column = int(sys.argv[5])  - 1 # -1 is for python index
bw_data = pyBigWig.open(bw_file)

bed_fp = open(bed_file, "r")
#out_fp = open(out_file, 'w')
#Only writing gzip file 
out_gz_fp = gzip.open(sys.argv[3], "wb")
out_enrichment_over_mean_fp = gzip.open (sys.argv[4], "wb") 
counter = 0
values = None
for line in bed_fp:
    line_items = line.strip().split() 
    values = bw_data.values(line_items[0], int(line_items[1]) - 1000 , int(line_items[1]) + 1000) # only from the center 
    min_val = bw_data.header()['minVal'] 
    values = [i if math.isnan(i) == False else min_val for i in values]
    enrichment_over_mean_vals = None
    total_val = sum (values)/len(values)
    if total_val != 0:
		    enrichment_over_mean_vals =[i/total_val for i in values] 
    else: 
		    enrichment_over_mean_vals =[0]*len(values) 
    print (bw_data.header()["minVal"])
    formatted_values = ["%0.3f" % el for el in values]
    formatted_enrichment_values = ["%0.4f" % el for el in enrichment_over_mean_vals]
    if line_items[strand_column] == "-": 
        formatted_values = formatted_values[::-1]
        formatted_enrichment_values = formatted_enrichment_values[::-1]
    # get the peak loci-info 
    #peak_loci = line_items[3].split("|")[0]
    #chr_info = "@".join([line_items[0], line_items[1], line_items[2], peak_loci])
    #chr_info = "@".join([line_items[0], line_items[1], line_items[2], line_items[3]])
    chr_info = "@".join(line_items)
    to_write = "\t".join(map(str, formatted_values))
    to_write_enrichment = "\t".join (map (str, formatted_enrichment_values))
    #out_fp.write(chr_info + "\t" + to_write +"\n")
    bytes_to_write = bytes(chr_info + "\t" + to_write +"\n", encoding="ascii")
    bytes_to_write_enrichment = bytes(chr_info + "\t" + to_write_enrichment +"\n", encoding="ascii")
    out_gz_fp.write(bytes_to_write)
    out_enrichment_over_mean_fp.write(bytes_to_write_enrichment) 
    counter = counter + 1 
bed_fp.close()
#out_fp.close()
out_gz_fp.close()
out_enrichment_over_mean_fp.close()

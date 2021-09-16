import os
import sys
import re
import pickle
from collections import defaultdict
from collections import OrderedDict
import numpy as np 
from length_and_loc_with_absolute import get_real_footprint_length_with_abs_start
# it takes the input file from rule : plot_footprint_length_vs_per_orange 
# input: real_footprint_length_and_per_orange 
# sample file name: actual_footprint_length_in_clusters/merged_S2_to_reordered_20clust_147_50PNE_open_v1_peak_cluster_10_lf_14_rf_14_methylation_matrix_for_site_peak_5442_nclust_3_sam_flag_83~163_footprint_length.per_orange.tsv 

inp_fp = open(sys.argv[1])
per_methylation_vec_dict = pickle.load(open(sys.argv[2], "rb"))
footprint_fp = open(sys.argv[3]) 
occluded_dict = pickle.load(open(sys.argv[4], "rb"))
footprint_dict = pickle.load(open(sys.argv[5], "rb"))
regions_metadata = pickle.load(open(sys.argv[6], "rb"))
out_fp_binding_state = open(sys.argv[7], "w")
lflank =  int(sys.argv[8])
rflank =  int(sys.argv[9]) 

label_dict = {
    "Naked-DNA" : "0",
    "TF" : "1",
    "Nuc":  "2",
    "discard" : "3"
}

fpt_dict = OrderedDict()
fpt_abs_start_dict = defaultdict(list)
fpt_length_dict = defaultdict(list)
complete_fp_len_dict = defaultdict()
for line in footprint_fp:
    d_loc = [m.start() for m in re.finditer("\t", line)]
    k = line[0:d_loc[0]]
    roi_id = k.split("`")[0]
    footprint_str = line[d_loc[0]+1:d_loc[1]]
    fpt_dict[k] = footprint_str 
    complete_footprint_vec = footprint_dict[k]["footprint"] 
    complete_footprint_chr_start = int(footprint_dict[k]["start"])
    complete_footprint_chr_end = int(footprint_dict[k]["end"])
    fp_center = regions_metadata[roi_id]["center"]
    m_vec_start = int(fp_center)  - complete_footprint_chr_start -  lflank # absolute location on read 
    m_vec_stop = int(fp_center)  - complete_footprint_chr_start +  rflank + 1  # absolute location on read
    lvec, _, _, abs_loc = get_real_footprint_length_with_abs_start(footprint_str, m_vec_start, m_vec_stop, complete_footprint_vec)
    fpt_abs_start_dict[k] = abs_loc 
    fpt_length_dict[k] = lvec
    complete_fp_len_dict[k] = len(complete_footprint_vec)
    
    
 

labels_on_reads =  defaultdict(list)
contribution_to_percentage_methylation  =  defaultdict(list) 

for line in inp_fp:  
    line_length = len(line)
    d_loc = [m.start() for m in re.finditer("\t", line)]
    percent = float(line[d_loc[0] + 1: d_loc[1]] )
    flen = int(line[d_loc[1] + 1: d_loc[2]])
    read_id = line[0:d_loc[0]] 
    per_orange_for_fp = float(line[d_loc[2]+1: line_length - 1]) 

################################# binding states #####################
# use occluded to decide naked or occluded DNA 
######################################################################
    max_of_edge = occluded_dict[read_id]["max_of_edge"]  
    pcap_total = occluded_dict[read_id]["percentage_capital"]

    if (percent <= 30) and float (per_methylation_vec_dict[read_id]["per_c"]) > 25: 
        labels_on_reads[read_id].append("Naked-DNA")
    elif (percent <= 30) and float (per_methylation_vec_dict[read_id]["per_c"]) <=25:
        labels_on_reads[read_id].append("Nuc")
   
    elif percent>30 and flen <=50 : 
        labels_on_reads[read_id].append("TF")
        contribution_to_percentage_methylation[read_id].append(per_orange_for_fp)
    else:
        labels_on_reads[read_id].append("Nuc")
        contribution_to_percentage_methylation[read_id].append(per_orange_for_fp)

for k in fpt_dict:

    if len(contribution_to_percentage_methylation[k]) == 0: # it should be naked
        to_write = k + "#" + label_dict["Naked-DNA"] + "\t" + fpt_dict[k]
        out_fp_binding_state.write(to_write + "\n")
    else: 
        max_contributor_id = np.argmax(contribution_to_percentage_methylation[k])
        abs_start_loc = fpt_abs_start_dict[k][max_contributor_id] 
        length_to_check = fpt_length_dict[k][max_contributor_id]
        lab = labels_on_reads[k][max_contributor_id] # choose whatever comes first 
        if (abs_start_loc == 0) and length_to_check <= 50:
            lab = "discard" # discard
        elif (abs_start_loc + length_to_check == complete_fp_len_dict[k]) and (length_to_check<=50):
            lab = "discard" # discard
        to_write = k + "#" + label_dict[lab] + "\t" + fpt_dict[k]
        out_fp_binding_state.write(to_write + "\n")
        
inp_fp.close()
out_fp_binding_state.close()
footprint_fp.close()

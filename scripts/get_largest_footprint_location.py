import sys
import numpy as np 
import re
from length_and_loc_with_absolute import  get_real_footprint_length_with_abs_start 

    # (m_vec, m_vec_start, m_vec_stop, complete_vec):
    # return return_list, out_vec, loc_first, abs_loc
    #"""
    #extract complete footprint length from the selected methylation vector

inp_fp = open(sys.argv[1])

for line in inp_fp:
    # example line : head -1 suppressed_merged_S2_to_all_open_and_closed_mnase_peaks_cluster_199_lf_15_rf_15_site_peak_229_3_sam_flag_83~163_extend_from_peak_left_150_right_150.cond8.footprint_ordered.tsv
    # chr2L:480305-480305^199`SRR3133328.21543131_21543131/1_overlapping`83~163#0	.................................................................................................................................................................FFFFFFMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    d_loc = [m.start() for m in re.finditer("\t", line)]
    read_id = line[0:d_loc[0]]
    fp_str = line[d_loc[0] + 1 : len(line) - 1]
    #print(read_id)
  
    m_replaced_fp = fp_str.replace("M", ".")
    lv, ov, lf, abs_loc = get_real_footprint_length_with_abs_start(m_replaced_fp, 0, len(fp_str), m_replaced_fp)
    if len(lv) >0:
        max_idx = np.argmax(lv)
        max_start = abs_loc[max_idx] - 150 
        to_write = read_id + "\t" + str(max_start) + "\t" + str(lv[max_idx]) + "\t" + fp_str
        print(to_write)
    else:
        to_write = read_id + "\t" + str(-150) + "\t" + str(0) + "\t" + fp_str
        print(to_write) 
       
    

inp_fp.close()

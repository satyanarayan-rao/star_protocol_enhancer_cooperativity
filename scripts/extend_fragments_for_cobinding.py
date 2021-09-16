import os
import sys
import pickle
import re

# zless dsmf_reads_to_peak_pairs/dsmf_reads_to_mnase_peak_pairs.bed.gz | head -2
# chr2L	107887	108067	peak_2105_9@peak_2105_1	chr2L:107902-107902	+	0	150	chr2L	107797	108072	SRR3133329.43378799_43378799/1_overlapping`83~163	.	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF......................................................................................................................................................................................FFF	...........................................................................................Z..H.Z.....H....Z....H.ZX.Z.....X..X................Z.....................................Z.....Z.....................H.Z..........................Z.....X..Z.........Z........Z...H.Z..z	ACTATAATCCATATCTTATCTAATCCTTTATTCAAAACCACAATAAATATAAAATAAAAAATCTAAATCTAAATTCACTACAATTCACTCCGAAGCGAGATAGCAACGCTAAGCGGCGAAACAGCTGCTTAAAACCCAAAAACGTAAAAATACTAAAAAAATAAAATCATTTAAAAATACCGAAAACGCTAAAAAAATTTCCCTAAAAAGCGCCCAAATAAATCTCAAAATAATAATCGAGTCAGCCGAAAAGAACCGCCTCCTCCGTTTGCGACA
# chr2L	107887	108067	peak_2105_9@peak_2105_1	chr2L:107902-107902	+	0	150	chr2L	107806	108088	SRR3133326.32099076_32099076/1_overlapping`99~147	.	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF..........................................................................................................................................................................................................	.................................................................................Z....Z.......H..Z.X....Z..Z.......X..H..............Z.....................................Z.....Z.H.....................Z.H........................Z.......XZ.........Z.H......Z.....Z..Z......H.....H.Z.H	TATATTTTGTTTGATTTTTTGTTTAAAATTATAGTGAGTATAAAATAAAGGATTTGGATTTAGATTTATTACAGTTTACTTCGGAGCGGGGTGGCAACGCTGGGCGGCGGAATAGCTGCTTAAAATCCAAAAACGTGGAGATATTGAAGAAGTAAAATTATTTGAGAGTATCGAGAACGCTAGGAAAGTTTTCCTAAAGAGCGCTTAAATAAGTTTTAAAATAATAATCGAGTTAGCCGAAGAGAATCGCTTTTTTCGTTTGCGACGATAAGCAAAAGCTCGC

footprint_dict = pickle.load(open(sys.argv[1], "rb"))
roi_dict = pickle.load(open(sys.argv[2], "rb"))
lflank = int(sys.argv[3])
rflank = int(sys.argv[4])
lextend = int(sys.argv[5])
rextend = int(sys.argv[6])

out_fp = open(sys.argv[7], "w")
out_verb_fp = open(sys.argv[8], "w")
for line in sys.stdin: 
    d_loc = [m.start() for m in re.finditer("\t", line)] 
    
    cobinding_info = line[0:d_loc[12]]
    roi_id = line[d_loc[2] + 1: d_loc[3]].split ("@")[0] 
    key_without_hash = roi_id + "`" + line[d_loc[10]+1:d_loc[11]] 
    #print (key_without_hash)
    existing_footprint = line[d_loc[12] + 1: d_loc[13]]  
    complete_footprint = footprint_dict[key_without_hash]["footprint"]
    complete_mvec = footprint_dict[key_without_hash]["mvec"]
    complete_bsseq = footprint_dict[key_without_hash]["bs_seq"]
    complete_footprint_start = int(footprint_dict[key_without_hash]["start"])
    complete_footprint_end = int(footprint_dict[key_without_hash]["end"])
    peak_center = int(roi_dict[roi_id]["center"]) 
    strand = roi_dict[roi_id]["strand"]
    m_vec_start = peak_center - complete_footprint_start -  lflank
    m_vec_stop = peak_center - complete_footprint_start  + rflank + 1
    methylation_string = complete_footprint[m_vec_start: m_vec_stop]
    ####
    # lf = -3, rf = 3 
    # m_vec_start = 26
    #  m_vec_stop = 33 
    # extend_start = m_vec_start - (lextend - lf) 
    # extend_stop = m_vec_stop + (rextend - rf)  
    # 
    # complete_footprint_start and complete_footprint_end are in closed form length of footprint = 288, `11806747 - 11806460 + 1` 
    # actual  2324252627282930313233343536 
    #          . . . . . F F F F . . . F F 
    #               -3-2-1 0 1 2 3
    # desired -5 7 
    #          . . . . . F F F F . . . F F
    #           -5-4-3-2-1 0 1 2 3 4 5 6 7 
    # 
    #### 
#lflank = int(sys.argv[3])
#rflank = int(sys.argv[4])
#lextend = int(sys.argv[5])
#rextend = int(sys.argv[6])
#    if strand == "-": 
#        lexetend = int(sys.argv[6])
#        rextend = int(sys.argv[5])
#        lflank = int(sys.argv[4])
#        rflank = int(sys.argv[3]) 
        
    desired_length = rextend - (0 - lextend) + 1 
    extend_start = peak_center - (lextend) # close interval
    extend_stop = peak_center + (rextend + 1 )  # open 
    extended_footprint = ""
    extended_mvec = ""
    extended_bsseq = ""
    #print ([extend_start, extend_stop, complete_footprint_start, complete_footprint_end])
    #sys.exit(-1)
    if (extend_start >= complete_footprint_start) and  (extend_stop <= complete_footprint_end + 1): 
        #print ([cobinding_info, "condition 1"])
        if strand == "+":
            extended_footprint = complete_footprint[extend_start - complete_footprint_start: extend_stop - complete_footprint_start]
            extended_mvec = complete_mvec[extend_start - complete_footprint_start: extend_stop - complete_footprint_start]
            extended_bsseq = complete_bsseq[extend_start - complete_footprint_start: extend_stop - complete_footprint_start]
        elif strand == "-": 
            extended_footprint = complete_footprint[extend_start - complete_footprint_start: extend_stop - complete_footprint_start][::-1]
            extended_mvec = complete_mvec[extend_start - complete_footprint_start: extend_stop - complete_footprint_start][::-1]
            extended_bsseq = complete_bsseq[extend_start - complete_footprint_start: extend_stop - complete_footprint_start][::-1]
            # recall that complete_footprint is on the watson strand always
            #print(["in", "in", key_without_hash,  extended_footprint, len(extended_footprint)])
    elif (extend_start >= complete_footprint_start) and  (extend_stop > complete_footprint_end + 1): 
        #print ([cobinding_info, "condition 2"])
        if strand == "+":
            extended_footprint = complete_footprint[extend_start - complete_footprint_start: complete_footprint_end - complete_footprint_start + 1]
            extended_mvec = complete_mvec[extend_start - complete_footprint_start: complete_footprint_end - complete_footprint_start + 1]
            extended_bsseq = complete_bsseq[extend_start - complete_footprint_start: complete_footprint_end - complete_footprint_start + 1]
            # add `M` of how much short in the end 
            extended_footprint = extended_footprint + "M"*(desired_length - len(extended_footprint))
            extended_mvec = extended_mvec + "M"*(desired_length - len(extended_mvec))
            extended_bsseq = extended_bsseq + "M"*(desired_length - len(extended_bsseq))
        elif strand == "-":
            extended_footprint = complete_footprint[extend_start - complete_footprint_start: complete_footprint_end - complete_footprint_start + 1][::-1] 
            extended_mvec = complete_mvec[extend_start - complete_footprint_start: complete_footprint_end - complete_footprint_start + 1][::-1] 
            extended_bsseq = complete_bsseq[extend_start - complete_footprint_start: complete_footprint_end - complete_footprint_start + 1][::-1] 
            #print([peak_center, extend_start, complete_footprint_start, complete_footprint_end])
            extended_footprint = "M"*(desired_length - len(extended_footprint)) + extended_footprint
            extended_mvec = "M"*(desired_length - len(extended_mvec)) + extended_mvec 
            extended_bsseq = "M"*(desired_length - len(extended_bsseq)) + extended_bsseq 
            #print(["in", "out", key_without_hash,  extended_footprint, len(extended_footprint)])
    elif (extend_start < complete_footprint_start) and (extend_stop <= complete_footprint_end + 1): 
        if strand == "+":
            extended_footprint = complete_footprint[0: extend_stop - complete_footprint_start]
            extended_mvec = complete_mvec[0: extend_stop - complete_footprint_start]
            extended_bsseq = complete_bsseq[0: extend_stop - complete_footprint_start]
            extended_footprint = "M"*(desired_length - len(extended_footprint)) + extended_footprint 
            extended_mvec = "M"*(desired_length - len(extended_mvec)) + extended_mvec 
            extended_bsseq = "M"*(desired_length - len(extended_bsseq)) + extended_bsseq
        #print(["out", "in", key_without_hash,  extended_footprint, len(extended_footprint)]) 
        elif strand == "-":
            extended_footprint = complete_footprint[0: extend_stop - complete_footprint_start][::-1]
            extended_mvec = complete_mvec[0: extend_stop - complete_footprint_start][::-1]
            extended_bsseq = complete_bsseq[0: extend_stop - complete_footprint_start][::-1]
            extended_footprint = extended_footprint + "M"*(desired_length - len(extended_footprint)) 
            extended_mvec = extended_mvec + "M"*(desired_length - len(extended_mvec)) 
            extended_bsseq = extended_bsseq + "M"*(desired_length - len(extended_bsseq)) 
            #print(["out", "in", key_without_hash,  extended_footprint, len(extended_footprint)]) 
        #print ([cobinding_info, "condition 3", extended_mvec])
        
    elif (extend_start < complete_footprint_start) and (extend_stop > complete_footprint_end + 1): 
        #print ([cobinding_info, "condition 4"])
        if strand == "+":
            to_add_in_left = "M"*(complete_footprint_start - extend_start)
            to_add_in_right = "M"*(extend_stop - (complete_footprint_end + 1) )
            extended_footprint = to_add_in_left + complete_footprint + to_add_in_right   
            extended_mvec = to_add_in_left + complete_mvec + to_add_in_right   
            extended_bsseq = to_add_in_left + complete_bsseq + to_add_in_right   
        elif strand == "-":
            to_add_in_left = "M"*(complete_footprint_start - extend_start)
            to_add_in_right = "M"*(extend_stop - (complete_footprint_end + 1) )
            extended_footprint = to_add_in_right + complete_footprint[::-1] + to_add_in_left
            extended_mvec = to_add_in_right + complete_mvec[::-1] + to_add_in_left
            extended_bsseq = to_add_in_right + complete_bsseq[::-1] + to_add_in_left
            #print(["out", "out", key_without_hash, extended_footprint, len(extended_footprint)])
    
    #print ([extend_start, extend_stop, complete_footprint_start, complete_footprint_end, extended_footprint])

    out_fp.write(cobinding_info + "\t" + extended_footprint + "\n")
    out_verb_fp.write(cobinding_info + "\t" + extended_footprint + "\n")
    out_verb_fp.write(cobinding_info + "\t" + extended_mvec + "\n")
    out_verb_fp.write(cobinding_info + "\t" + extended_bsseq + "\n")

out_fp.close()
out_verb_fp.close()
    
#    prefix = [] 
#    for left in range(lextend - lflank): 
#        if m_vec_start - left - 1 >= complete_footprint_start:
#            prefix.append(m_vec_start - complete_footprint_start - left - 1)
#        else:
#            prefix.append("M") 
#    prefix_str = "".join(prefix)[::-1] 
#    for right in range() 

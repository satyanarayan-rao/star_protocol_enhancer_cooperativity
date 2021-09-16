import os
import sys
from collections import defaultdict
import re

inp_fp = open(sys.argv[1])
cnt_dict = defaultdict(lambda : defaultdict(lambda : 0))
is_read_counted_for_peak = defaultdict(lambda : defaultdict(lambda : False))
for line in inp_fp:
    # example line : head -2  with_fp_peak_pairs.tsv 
    #chr2L	107887	108067	peak_2105_9@peak_2105_1	chr2L:107902-107902	+	0	150	chr2L	107797	108072	SRR3133329.43378799_43378799/1_overlapping`83~163	.	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF......................................................................................................................................................................................FFF	...........................................................................................Z..H.Z.....H....Z....H.ZX.Z.....X..X................Z.....................................Z.....Z.....................H.Z..........................Z.....X..Z.........Z........Z...H.Z..z	ACTATAATCCATATCTTATCTAATCCTTTATTCAAAACCACAATAAATATAAAATAAAAAATCTAAATCTAAATTCACTACAATTCACTCCGAAGCGAGATAGCAACGCTAAGCGGCGAAACAGCTGCTTAAAACCCAAAAACGTAAAAATACTAAAAAAATAAAATCATTTAAAAATACCGAAAACGCTAAAAAAATTTCCCTAAAAAGCGCCCAAATAAATCTCAAAATAATAATCGAGTCAGCCGAAAAGAACCGCCTCCTCCGTTTGCGACA
    #chr2L	107887	108067	peak_2105_9@peak_2105_1	chr2L:107902-107902	+	0	150	chr2L	107806	108088	SRR3133326.32099076_32099076/1_overlapping`99~147	.	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF..........................................................................................................................................................................................................	.................................................................................Z....Z.......H..Z.X....Z..Z.......X..H..............Z.....................................Z.....Z.H.....................Z.H........................Z.......XZ.........Z.H......Z.....Z..Z......H.....H.Z.H	TATATTTTGTTTGATTTTTTGTTTAAAATTATAGTGAGTATAAAATAAAGGATTTGGATTTAGATTTATTACAGTTTACTTCGGAGCGGGGTGGCAACGCTGGGCGGCGGAATAGCTGCTTAAAATCCAAAAACGTGGAGATATTGAAGAAGTAAAATTATTTGAGAGTATCGAGAACGCTAGGAAAGTTTTCCTAAAGAGCGCTTAAATAAGTTTTAAAATAATAATCGAGTTAGCCGAAGAGAATCGCTTTTTTCGTTTGCGACGATAAGCAAAAGCTCGC
    d_loc = [m.start() for m in re.finditer("\t", line)]
    k = line[0:d_loc[7]]
    read_id = line[d_loc[10] +1: d_loc[11]] 
    if ("83~163" in line) and (is_read_counted_for_peak[k][read_id] == False): 
        cnt_dict[k]["83~163"] +=1
        is_read_counted_for_peak[k][read_id] = True
    elif ("99~147" in line) and (is_read_counted_for_peak[k][read_id] == False):
        cnt_dict[k]["99~147"] +=1
        is_read_counted_for_peak[k][read_id] = True
        
    
for key in cnt_dict: 
    print ("{}\t{}\t{}".format(key, cnt_dict[key]["99~147"], cnt_dict[key]["83~163"])) 

inp_fp.close()

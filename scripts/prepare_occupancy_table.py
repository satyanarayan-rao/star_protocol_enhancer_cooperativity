import sys
from collections import defaultdict
import gzip
import re
import pickle

peak_to_region_dict = {}  
fp = open (sys.argv[1])
for line in fp:
    # example line: head -2 input_bed/mnase_peaks_in_open_and_closed_enhancers.bed
    # chr2L	47389	47390	peak_877_1	closed	+
    # chr2L	107902	107903	peak_2105_9	open	+
    l_items = line.strip().split()
    peak_to_region_dict[l_items[3]] = l_items[4]

fp.close()
state_dict = defaultdict(lambda: defaultdict(lambda : 0))
peak_list = [] 
peak_dict = defaultdict(lambda : False)
for line in sys.stdin:
    first_col = line[0:line.find("\t")]  
    peak_id = first_col[0:first_col.find("`")] 
    state = first_col[first_col.find("#") + 1: len(first_col)]
    state_dict[peak_id][state] +=1
    if peak_dict[peak_id] == False:
        peak_list.append(peak_id) 
        peak_dict[peak_id] = True
    
accepted_states = ['0', '1', '2'] 
state_annotations = {'0' : "Naked", '1' : "TF", '2' : "Nucleosome"}
for peak in peak_list:
    total = 0 
    for s in accepted_states:
        total += state_dict[peak][s] 
    for s in accepted_states:
        per_occup = 100 * (state_dict[peak][s]/total)  
        print ("\t".join ([peak, state_annotations[s], str (per_occup), peak_to_region_dict[peak]]))
   

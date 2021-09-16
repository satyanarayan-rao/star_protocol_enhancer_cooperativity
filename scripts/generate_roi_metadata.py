import sys
from collections import defaultdict
import gzip
import re
import pickle

roi_dict = defaultdict(dict)
out_fp = open (sys.argv[1], "wb")
for line in sys.stdin:
    # example lines: head -2 input_bed/example.bed
    #chr2L	480305	480306	peak_229	.	+
    #chr3L	616610	616611	peak_613	.	+
    l_items = line.strip().split()
    roi_dict[l_items[3]]["strand"] = l_items[5] 
    roi_dict[l_items[3]]["center"] = l_items[1]
pickle.dump(dict(roi_dict), out_fp)
out_fp.close()

import sys
from collections import defaultdict
import re

sel_dict = defaultdict(lambda : False)
for v in range (int(sys.argv[1])): 
    sel_dict[str(v)] = True 
for line in sys.stdin:
    state_id = line[0:line.find("\t")].split("#")[-1]
    if sel_dict[state_id] == True:
        print (line, end= "")

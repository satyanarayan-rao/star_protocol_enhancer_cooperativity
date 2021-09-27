import sys
from collections import defaultdict
import re

to_scan = sys.argv[1]
for line in sys.stdin:
    idx = line.strip().split("\t")[0].split("#")[-1]
    if idx == to_scan:
        print (line, end = "") 

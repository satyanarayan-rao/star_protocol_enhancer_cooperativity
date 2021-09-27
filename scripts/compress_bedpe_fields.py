import sys
from collections import defaultdict
import re
import pickle
ncols = int(sys.argv[1])

#chr2L	19155158	19155266	chr2L	19155158	19155159	chr2L	19155265	19155266	peak_110_4_and_peak_110_6	.	.	.	chr2L	19154978	19155273	SRR3133329.7675245_7675245/1_overlapping`83~163	.	....................FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF.................................FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF..................................................FFFFFFFFFFFFFFFFFFFFFFF........................................................................................	h.z...x..x....h....Z...x..x..x.........z.........................X......Z........................H...................z.................H........H....H.............Z..Z.................H.zX.................h..H..X...........X..X............Z...Z...Z...Z...H...................H..X..Z..............	ACAACAACTACATTACATCGCCTACTACTACTAATTTCCAATTTCAAATATCACATTTTATTCCAGCTCAACGTCTTTCCCAAAAATTATTAAAAAAGCCATAACAAATATATAATCATTATATTTTATATGATTGCTTTATAAGCTAAGCTTTCCCCTCTCCGCCGTTCAAGAAAAACATCTAGCAGCAATAATTATTATTATTACTGCTGCTCTTATAACAGCTGCTAAAAGGAATCGCTCGATCGCTCGATTGCCATCCCCCTTCCCCCCAAGCAGCCGCCTCACTCTCACCC
for line in sys.stdin:
    l_items = line.strip().split("\t") 
    new_output = "\t".join ([l_items[0], l_items[1], l_items[2], l_items[9], l_items[0] + ":" + l_items[4] + "-" + l_items[8], 
                              "+", "0", str(int(l_items[8]) - int (l_items[4])), "\t".join(l_items[ncols:]) ]) 
    print(new_output)

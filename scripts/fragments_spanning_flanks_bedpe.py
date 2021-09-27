import os
import sys 
import re
# input stream:  peak center bed followed by fragment bed
# sys.argv[2]: left flank 
# sys.argv[3]: right flank
# output stream: output file

# example string to process: 

############ plus strand ###############
# chr2L   11806729        11806730        chr2L:11806729-11806729^20      .       +       chr2L   11806644        11806923        SRR3133326.203452_203452/1_overlapping  .       .h..h.H..x...h.H.Z.H.........Z....x...x........X...........h......h..zxhhhh.........z.....................xh...h.....hh..x...x..xhh.......h..h.h........h....h..z..........h.................X......Z.H.......h.....Zx....H.............h........Z....Z..........z.....X................
########################################

############ minus strand ###############
# chr2L   13205596        13205597        chr2L:13205596-13205596^20      .       -       chr2L   13205465        13205732        SRR3133326.98156_98156/1_overlapping    .       .hh.hh.h.............................hhh..............Z...hhh...Z...........x....x........hhx..Z......hh............h.h.h......xz.hx...z.....z..x..h.........h.....z....h.....h.............h................hh..z.xz....z......H...Z.H...hx.....h.Z..........Z...Z...Z..xZ.
########################################

lflank = int(sys.argv[1])
rflank = int(sys.argv[2])

for line in sys.stdin:
    d_loc = [m.start() for m in re.finditer("\t", line)]
    # check the strand first
    strand = line[d_loc[4] + 1: d_loc[5]]
    if strand == "+": 
        peak_start = int(line[d_loc[0] + 1: d_loc[1]])
        read_start = int(line[d_loc[6] + 1: d_loc[7]])
        read_end = int(line[d_loc[7] + 1: d_loc[8]])
        if (peak_start - read_start >= lflank) and \
           (read_end - peak_start >= rflank): 
           print (line, end = "")
        else:
            continue
    elif strand == "-":
        peak_start = int(line[d_loc[0] + 1: d_loc[1]])
        read_start = int(line[d_loc[6] + 1: d_loc[7]])
        read_end = int(line[d_loc[7] + 1: d_loc[8]])
        if (peak_start - read_start >= rflank) and \
           (read_end - peak_start >= lflank): 
           print (line, end = "")
        else:
            continue
    elif strand == ".": # closed enhancers peak strand is not available ; just copying the + strand code here
        peak_start = int(line[d_loc[0] + 1: d_loc[1]])
        read_start = int(line[d_loc[6] + 1: d_loc[7]])
        read_end = int(line[d_loc[7] + 1: d_loc[8]])
        if (peak_start - read_start >= lflank) and \
           (read_end - peak_start >= rflank): 
           print (line, end = "")
            
        else:
            continue        
    else: 
        continue

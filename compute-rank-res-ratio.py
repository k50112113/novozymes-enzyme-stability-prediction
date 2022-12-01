import os
import sys
import numpy as np

res_range = int(sys.argv[1])

seq_id_mutation_map = {}
with open("single-mutation-map.txt","r") as fin:
    for aline in fin:
        linelist = aline.strip().split()
        seq_id_mutation_map[int(linelist[1])] = linelist[0]

mutation_msd_map = {}
with open("res-com-msd-ratio-range=%d.txt"%(res_range),"r") as fin:
    fin.readline()
    for aline in fin:
        linelist = aline.strip().split()
        mutation_msd_map[linelist[0]] = np.array([float(i) for i in linelist[1:]])

seq_id_list = [i for i in range(31390,33802+1)]
foutrescom = open("submission-res-com-range=%s.csv"%(res_range),"w")
foutrescom.write("seq_id,tm\n")
for seq_id in seq_id_list:
    mutation = seq_id_mutation_map[seq_id]
    msd_ratio = mutation_msd_map[mutation][2] 
    foutrescom.write("%d,%f\n"%(seq_id,1.0/msd_ratio))
    
foutrescom.close()

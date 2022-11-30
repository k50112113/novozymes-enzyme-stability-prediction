import os
import sys
import numpy as np
sys.path.insert(1,os.path.expanduser('~')+'/python-include/')
from ZLabPlot import ZLabPlot
import ReadMD as RM


atom_mass = np.array([14.0067,12.011,15.9994,32.06,1.00797,14.0067,12.011])

seq_id_mutation_map = {}
with open("single-mutation-map.txt","r") as fin:
    for aline in fin:
        linelist = aline.strip().split()
        seq_id_mutation_map[int(linelist[1])] = linelist[0]

mutation_natom_map = {}
with open("number_of_atoms.txt","r") as fin:
    fin.readline()
    for aline in fin:
        linelist = aline.strip().split()
        mutation_natom_map[linelist[0]] = np.array([int(i) for i in linelist[1:]])

mutation_msd_map = {}
with open("msd-1ps.txt","r") as fin:
    fin.readline()
    for aline in fin:
        linelist = aline.strip().split()
        mutation_msd_map[linelist[0]] = np.array([float(i) for i in linelist[1:]])

seq_id_list = [i for i in range(31390,33802+1)]

foutaa = open("submission-aa.csv","w")
foutaa.write("seq_id,tm\n")
foutbb = open("submission-bb.csv","w")
foutbb.write("seq_id,tm\n")
for seq_id in seq_id_list:
    mutation = seq_id_mutation_map[seq_id]
    natom_mass = mutation_natom_map[mutation]*atom_mass
    msd = mutation_msd_map[mutation]
    
    weighted_msd = msd*natom_mass

    weighted_msd_aa = np.sum(weighted_msd[:-2])/np.sum(natom_mass[:-2])
    weighted_msd_bb = np.sum(weighted_msd[-2:])/np.sum(natom_mass[-2:])

    foutaa.write("%d,%f\n"%(seq_id,1.0/weighted_msd_aa))
    foutbb.write("%d,%f\n"%(seq_id,1.0/weighted_msd_bb))
    print(seq_id, mutation, weighted_msd_aa, weighted_msd_bb)
    
foutaa.close()
foutbb.close()

import numpy as np
import os
import sys
from ZLabPlot import ZLabPlot

mass_map = {'N':14.0067,'C':12.011,'O':15.9994,'S':15.9994,'H':1.00797}
res_name_map = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

def read_xyz(xyzfile):
    #argument
    #   xyzfile (string):       xyz filename
    #return
    #   coord (numpy 3D array): coord, (n_frame, n_atom, 3)
    #   element (list):         atom type
    coord = []
    element = []
    with open(xyzfile,'r') as fin:
        i = 0
        for aline in fin:
            linelist = aline.strip().split()
            if len(linelist) > 3:
                if i == 1: element.append(linelist[0])
                coord[-1].append([float(i) for i in linelist[1:]])
            elif len(linelist) == 1:
                i += 1
                coord.append([])
    coord = np.array(coord)
    return coord, element

def recenter(coord):
    coord = np.transpose(coord, axes = (1, 2, 0))
    coord -= coord[...,:1]
    coord = np.transpose(coord, axes = (2, 0, 1))
    return coord

def compute_dr_time_avg(coord, res_indices, window):
    if res_indices == 'all': coord_res = coord
    else: coord_res = coord[:,res_indices,:]
    dr_time_avg = np.mean(np.sum((coord_res[:-window]-coord_res[window:])**2, axis = 2), axis = 0)**0.5
    return dr_time_avg

def compute_r_time_avg(coord, res_indices):
    if res_indices == 'all': coord_res = coord
    else: coord_res = coord[:,res_indices,:]
    r_time_avg = np.mean(coord_res, axis = 0)
    return r_time_avg

# prefix = sys.argv[1]
# p = os.popen("find %s -type %s -name %s"%('.','d','\"%s*\"'%(prefix)))
p = os.popen("find %s -type %s -name %s"%('.','d','\"*\"'))
folder = p.read().split('\n')
folder.pop()
del folder[0]
folder = sorted(folder)
for i in range(len(folder)): folder[i] = folder[i].replace('./', "")
folder.remove('wild')
cutoff = 100
window = 1
res_cutoff = float(sys.argv[1])

coord_wt, element_wt = read_xyz('WT/md-traj-kabsch-com.xyz')
coord_wt = coord_wt[cutoff:]
dr_wt = compute_dr_time_avg(coord_wt, 'all', window)
r_wt = compute_r_time_avg(coord_wt, 'all')

fout = open("res-com-msd-ratio-range-3d=%.2f.txt"%(res_cutoff), 'w')
fout.write("mutation %12s %12s %12s %6s\n"%('dr_mu','dr_wt','dr_mu/dr_wt', 'n_res'))

for i in range(len(folder)):
    coord_mu, element_mu = read_xyz('%s/md-traj-kabsch-com.xyz'%(folder[i]))
    seq_length = len(element_mu)
    coord_mu = coord_mu[cutoff:]
    if folder[i] == 'WT':
        print(i, folder[i])
        dr_mu_tmp = 1.0
        dr_wt_tmp = 1.0
    elif ord('Z') >= ord(folder[i][-1]) >= ord('A'):
        target_res_index = int(folder[i][1:-1]) - 1
        res_to_target_res_distance_wt = np.sum((r_wt - r_wt[target_res_index])**2, axis = 1)**0.5
        res_indices_wt = np.where(res_to_target_res_distance_wt <= res_cutoff)[0]
        r_mu = compute_r_time_avg(coord_mu, 'all')
        res_indices_mu = res_indices_wt
        print(i, folder[i], list(res_indices_mu))
    else:
        target_res_index = int(folder[i][1:]) - 1
        res_to_target_res_distance_wt = np.sum((r_wt - r_wt[target_res_index])**2, axis = 1)**0.5
        res_indices_wt = np.where(res_to_target_res_distance_wt <= res_cutoff)[0]
        res_indices_wt = np.delete(res_indices_wt, np.where(res_indices_wt == target_res_index)[0])
        res_indices_mu = np.where(res_indices_wt > target_res_index, res_indices_wt-1, res_indices_wt)
        print(i, folder[i], list(res_indices_wt), list(res_indices_mu))
    dr_mu = compute_dr_time_avg(coord_mu, res_indices_mu, window)
    dr_wt_tmp = np.mean(dr_wt[res_indices_wt])
    dr_mu_tmp = np.mean(dr_mu)
    fout.write('%8s %3.6e %3.6e %3.6e %6s\n'%(folder[i], dr_mu_tmp, dr_wt_tmp, dr_mu_tmp/dr_wt_tmp, len(res_indices_mu)))
fout.close()
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

def compute_dr(coord, res_indices, window = 1):
    #compute sqrt(<dr(t)^2>)
    if res_indices == 'all': coord_res = coord
    else: coord_res = coord[:,res_indices,:]
    dr = np.mean((coord_res[:-window]-coord_res[window:])**2, axis = (1, 2))**0.5
    return dr

def compute_dr_time_avg(coord, res_indices, window):
    if res_indices == 'all': coord_res = coord
    else: coord_res = coord[:,res_indices,:]
    dr_time_avg = np.mean(np.sum((coord_res[:-window]-coord_res[window:])**2, axis = 2), axis = 0)**0.5
    return dr_time_avg

def compute_r_r0(coord, res_indices):
    #compute |r(t)-r(0)|
    if res_indices == 'all': r = np.mean(np.sum(coord_mu[:,:,:]**2, axis = 2), axis = 1)**0.5
    else: r = np.mean(np.sum(coord_mu[:,res_indices,:]**2, axis = 2), axis = 1)**0.5
    return r

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)

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
res_range = int(sys.argv[1])

coord_wt, element_wt = read_xyz('WT/md-traj-kabsch-com.xyz')
coord_wt = coord_wt[cutoff:]
dr_wt = compute_dr_time_avg(coord_wt, 'all', window)


fout = open("res-com-msd-ratio-range=%d.txt"%(res_range), 'w')
fout.write("mutation %12s %12s %12s\n"%('dr_mu','dr_wt','dr_mu/dr_wt'))

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
        lowerbound = max(target_res_index-res_range, 0)
        upperbound = min(target_res_index+1+res_range, seq_length)
        res_indices_mu = range(lowerbound, upperbound)
        res_indices_wt = res_indices_mu
        print(i, folder[i], list(res_indices_mu))
    else:
        target_res_index = int(folder[i][1:]) - 1
        lowerbound = max(target_res_index-1-res_range, 0)
        upperbound = min(target_res_index+1+res_range, seq_length)
        res_indices_mu = range(lowerbound, upperbound)
        lowerbound = max(target_res_index-1-res_range, 0)
        upperbound = min(target_res_index+1+1+res_range, seq_length)
        res_indices_wt = range(lowerbound, upperbound)
        res_indices_wt = list(res_indices_wt)
        res_indices_wt.remove(target_res_index)
        print(i, folder[i], list(res_indices_wt), list(res_indices_mu))
    dr_mu = compute_dr_time_avg(coord_mu, res_indices_mu, window)
    dr_wt_tmp = np.mean(dr_wt[res_indices_wt])
    dr_mu_tmp = np.mean(dr_mu)
    fout.write('%8s %3.6e %3.6e %3.6e\n'%(folder[i], dr_mu_tmp, dr_wt_tmp, dr_mu_tmp/dr_wt_tmp))
fout.close()

### plotting
# prefix = sys.argv[1]
# p = os.popen("find %s -type %s -name %s"%('.','d','\"%s*\"'%(prefix)))
# folder = p.read().split('\n')
# folder.pop()
# del folder[0]
# folder = sorted(folder)
# for i in range(len(folder)): folder[i] = folder[i].replace('./', "")
# folder = ['WT'] + folder
# r = []
# dr = []
# dr_ma = []
# ma_length = 95
# target_res_index = int(prefix[1:])
# res_range = 0
# res_indices = range(target_res_index-res_range,target_res_index+res_range+1)
# res_indices = 'all'
# for a_folder in folder:
#     print(a_folder, end='\r')
#     coord_mu, element_mu = read_xyz('%s/md-traj-kabsch-com.xyz'%(a_folder))
#     coord_mu = recenter(coord_mu[cutoff:])
#     r.append(compute_r_r0(coord_mu, res_indices))
#     coord_mu = coord_mu[cutoff:]
#     dr.append(compute_dr(coord_mu, res_indices, window))
#     dr_ma.append(running_mean(dr[-1]**2, ma_length)**0.5)
# time = [(np.arange(len(dr[0])))*0.1]*len(dr)
# time_ma = [(ma_length+np.arange(len(dr_ma[0])))*0.1]*len(dr_ma)
# zp = ZLabPlot(rcmap = {'figure.figsize': (12,8)})
# zp.add_subplot()
# zp.add_data(time+time_ma, dr+dr_ma, \
# legend = folder*2, \
# lw = [0.5]*len(folder) + [2]*len(folder), \
# color = [i for i in range(len(folder))]*2, \
# xlabel = 'time (ps)', ylabel = '$\sqrt{<dr^2>} (\mathrm{\AA})$', \
# ncol = 5,
# xlim = (8, 10))
# zp.save("dr_ma", transparent = False)
# zp.clear()
# time = [np.arange(len(r[0]))*0.1]*len(r)
# zp = ZLabPlot(rcmap = {'figure.figsize': (12,8)})
# zp.add_subplot()
# zp.add_data(time, r, legend = folder, xlabel = 'time (ps)', ylabel = '$|r-r_0| (\mathrm{\AA})$', ncol = 5)
# zp.save("r", transparent = False)
# zp.clear()
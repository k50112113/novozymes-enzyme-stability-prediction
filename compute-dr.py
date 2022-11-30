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
    if res_indices == 'all': coord_res = coord
    else: coord_res = coord[:,res_indices,:]
    dr = np.mean((coord_res[:-window]-coord_res[window:])**2, axis = (1, 2))**0.5
    return dr

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)

prefix = sys.argv[1]
p = os.popen("find %s -type %s -name %s"%('.','d','\"%s*\"'%(prefix)))
folder = p.read().split('\n')
folder.pop()
del folder[0]
folder = sorted(folder)
for i in range(len(folder)): folder[i] = folder[i].replace('./', "")

folder = ['WT'] + folder[:5]
r = []
dr = []
dr_ma = []
cutoff = 100
window = 1
ma_length = 20
target_res_index = int(prefix[1:])
res_range = 10
res_indices = range(target_res_index-res_range,target_res_index+res_range+1)
# res_indices = 'all'
for a_folder in folder:
    print(a_folder, end='\r')
    coord_mu, element_mu = read_xyz('%s/md-traj-kabsch-com.xyz'%(a_folder))
    coord_mu = recenter(coord_mu[cutoff:])
    # coord_mu = coord_mu[cutoff:]
    if res_indices == 'all':
        r2_sum = np.sum(coord_mu[:,:,:]**2, axis = (1, 2))
    else:
        r2_sum = np.sum(coord_mu[:,res_indices,:]**2, axis = (1, 2))
    r.append(r2_sum**0.5)
    rstd = np.std(r2_sum[10:], ddof = 1)
    rmean = np.mean(r2_sum[10:])
    print(rmean, rstd)

    dr.append(compute_dr(coord_mu, res_indices, window))
    dr_ma.append(running_mean(dr[-1]**2, ma_length)**0.5)

time = [np.arange(len(dr[0]))*0.1]*len(dr)
zp = ZLabPlot(rcmap = {'figure.figsize': (12,8)})
zp.add_subplot()
zp.add_data(time, dr, legend = folder, xlabel = 'time (ps)', ylabel = '$\sqrt{<dr^2>} (\mathrm{\AA})$', ncol = 5)
zp.save("dr", transparent = False)
zp.clear()

time_ma = [(ma_length+np.arange(len(dr_ma[0])))*0.1]*len(dr_ma)
zp = ZLabPlot(rcmap = {'figure.figsize': (12,8)})
zp.add_subplot()
zp.add_data(time+time_ma, dr+dr_ma, \
legend = folder*2, \
lw = [0.5]*len(folder) + [2]*len(folder), \
color = [i for i in range(len(folder))]*2, \
xlabel = 'time (ps)', ylabel = '$\sqrt{<dr^2>} (\mathrm{\AA})$', \
ncol = 5)
zp.save("dr_ma", transparent = False)
zp.clear()

time = [np.arange(len(r[0]))*0.1]*len(r)
zp = ZLabPlot(rcmap = {'figure.figsize': (12,8)})
zp.add_subplot()
zp.add_data(time, r, legend = folder, xlabel = 'time (ps)', ylabel = '$|r-r_0| (\mathrm{\AA})$', ncol = 5)
zp.save("r", transparent = False)
zp.clear()
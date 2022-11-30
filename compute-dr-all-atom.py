import numpy as np
import os
import sys

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

prefix = sys.argv[1]
p = os.popen("find %s -type %s -name %s"%('.','d','\"%s*\"'%(prefix)))
# p = os.popen("find %s -type %s -name %s"%('.','d','\"*\"'))
folder = p.read().split('\n')
folder.pop()
del folder[0]
folder = sorted(folder)
for i in range(len(folder)): folder[i] = folder[i].replace('./', "")

fout = open("msd-0.1ps.txt","w")
fout.write('mutation ')
for sufix in ['N','C','O','S','H','N-bb','C-bb']: fout.write("%10s "%(sufix))
fout.write('\n')

foutstd = open("msd-0.1ps-std.txt","w")
foutstd.write('mutation ')
for sufix in ['N','C','O','S','H','N-bb','C-bb']: foutstd.write("%10s "%(sufix))
foutstd.write('\n')

foutn = open("number_of_atoms.txt","w")
foutn.write('mutation ')
for sufix in ['N','C','O','S','H','N-bb','C-bb']: foutn.write("%7s "%(sufix))
foutn.write('\n')

for a_folder in folder:
    print(a_folder)
    coord, element = read_xyz('%s/md-traj-kabsch.xyz'%(a_folder))

    indices = {\
    'N': [i for i, x in enumerate(element) if x[0] == "N"],
    'C': [i for i, x in enumerate(element) if x[0] == "C"],
    'O': [i for i, x in enumerate(element) if x[0] == "O"],
    'S': [i for i, x in enumerate(element) if x[0] == "S"],
    'H': [i for i, x in enumerate(element) if x[0] == "H"],
    'N-bb': [i for i, x in enumerate(element) if x == "N"],
    'C-bb': [i for i, x in enumerate(element) if x == "C" or x == "CA"]}

    n = [indices[i] for i in indices.keys()]
    foutn.write('%8s '%(folder[i]))
    for a in n: foutn.write('%7s '%(a))
    foutn.write('\n')

    coord = np.array(coord)[100:]
    n_frames = len(coord)
    natom = len(coord[0])

    dr_ = np.diff(coord, n = 1, axis = 0)**2

    fout.write('%8s '%(a_folder))
    foutstd.write('%8s '%(a_folder))
    for key in indices.keys():
        dr_frame = np.mean(dr_[:,indices[key],:], axis = (1, 2))
        dr_mean = np.mean(dr_frame)
        dr_std = np.std(dr_frame, ddof = 1)
        fout.write("%3.4e "%(dr_mean**0.5))
        foutstd.write("%3.4e "%(dr_std**0.5))
    fout.write('\n')
    foutstd.write('\n')

foutn.close()
fout.close()
foutstd.close()



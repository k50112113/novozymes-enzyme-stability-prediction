import numpy as np
import os
import sys

filename = sys.argv[1]

coord = []
element = []

with open(filename,'r') as fin:
    for aline in fin:
        if "step" in aline:
            coord.append([])
            natom = int(fin.readline())
        else:
            linelist = aline.strip().split()
            if len(linelist) > 3:
                coord[-1].append([float(i) for i in linelist[3:]])
                if len(coord) == 1: element.append(linelist[1])
            else: lbox = float(linelist[0])

coord = np.array(coord)
n_frames = len(coord)

fout = open(filename[:-3]+'xyz','w')
for j, a_coord in enumerate(coord):
    fout.write("%d\n\n"%(natom))
    for i in range(natom): fout.write("%s %.2f %.2f %.2f\n"%(element[i],a_coord[i][0]*10,a_coord[i][1]*10,a_coord[i][2]*10))
fout.close()





import numpy as np
import os
import sys

p = os.popen("find %s -type %s -name %s"%('./*/','f','md-traj1.xyz'))
folder = p.read().split('\n')
folder.pop()
folder = sorted(folder)
for i in range(len(folder)): folder[i] = folder[i].replace('./', "")

for filename in folder:
    print(filename, filename[:-5]+'.xyz')
    fout = open(filename[:-5]+'.xyz','w')
    with open(filename,'r') as fin:
        for aline in fin:
            linelist = aline.strip().split()
            if len(linelist) > 3:
                coord = [float(i) for i in linelist[1:]]
                element = linelist[0]
                fout.write("%s %.3f %.3f %.3f\n"%(element,coord[0]*10,coord[1]*10,coord[2]*10))
            else:
                fout.write(aline)
    fout.close()





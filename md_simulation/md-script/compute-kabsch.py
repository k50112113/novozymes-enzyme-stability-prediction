import numpy as np
import os
import sys

def kabsch(reference, target, ref_center, tar_center):
    #reference, target (N x dim)
    P = target-tar_center
    Q = reference-ref_center
    R = compute_rotation_matrix(P, Q)
    return np.transpose(np.matmul(R,np.transpose(P)))

def compute_rotation_matrix(P, Q):
    dim = P.shape[1]
    H = np.matmul(np.transpose(P),Q)
    U,S,Vh = np.linalg.svd(H)
    V = np.transpose(Vh)
    d = np.sign(np.linalg.det(np.matmul(V, np.transpose(U))))
    I = np.identity(dim)
    I[dim-1][dim-1] = d
    R = np.matmul(V, np.matmul(I,np.transpose(U)))
    return R

filename_ = sys.argv[1:]

for filename in filename_:
    print("%s running."%(filename))
    coord = []
    element = []
    with open(filename,'r') as fin:
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
    n_frames = len(coord)
    natom = len(coord[0])

    center = np.mean(coord, axis = 1)
    for i in range(n_frames):
        cc = kabsch(coord[0], coord[i], center[0], center[i])
        coord[i] = cc

    fout = open(filename[:-4]+'-kabsch.xyz','w')
    for j, a_coord in enumerate(coord):
        fout.write("%d\n\n"%(natom))
        for i in range(natom): fout.write("%s %.3f %.3f %.3f\n"%(element[i],a_coord[i][0],a_coord[i][1],a_coord[i][2]))
    fout.close()

    print("%s done."%(filename))





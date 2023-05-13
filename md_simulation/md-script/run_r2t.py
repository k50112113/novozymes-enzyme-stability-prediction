import os
import sys
sys.path.insert(1,os.path.expanduser('~')+'/zlab/python-include/')
import ReadMD as RM


p = os.popen("find %s -type %s -name %s"%('.','d','\"*\"'))
folder = p.read().split('\n')
folder.pop()
del folder[0]
folder = sorted(folder)
for i in range(len(folder)): folder[i] = folder[i].replace('./', "")
folder.pop()

folder = ['WT','wild']

foutn = open("number_of_atoms.txt","w")
foutn.write('mutation ')
for sufix in ['N','C','O','S','H','N-bb','C-bb']: foutn.write("%7s "%(sufix))
foutn.write('\n')

for i in range(len(folder)):
    print(folder[i])
    
    traj = RM.XYZ_ITR('%s/md-traj-kabsch.xyz'%(folder[i]))
    n = []
    n.append(len([i for i, x in enumerate(traj.type) if x[0] == "N"]))
    n.append(len([i for i, x in enumerate(traj.type) if x[0] == "C"]))
    n.append(len([i for i, x in enumerate(traj.type) if x[0] == "O"]))
    n.append(len([i for i, x in enumerate(traj.type) if x[0] == "S"]))
    n.append(len([i for i, x in enumerate(traj.type) if x[0] == "H"]))
    n.append(len([i for i, x in enumerate(traj.type) if x == "N"]))
    n.append(len([i for i, x in enumerate(traj.type) if x == "C" or x == "CA"]))
    foutn.write('%8s '%(folder[i]))
    for a in n: foutn.write('%7s '%(a))
    foutn.write('\n')
    
    el_map = {}
    with open("%s/md-traj-kabsch.xyz"%(folder[i]), "r") as fin:
        natom = int(fin.readline())
        fin.readline()
        for j in range(natom):
            el = fin.readline().strip().split()[0]
            if el_map.get(el[0]) is None:
                el_map[el[0]] = [el]
            elif el not in el_map[el[0]]:
                el_map[el[0]].append(el)
    el_map['N-bb'] = ['N']
    el_map['C-bb'] = ['CA','C']
    for key in el_map.keys():
        key = 'N'
        elstr = ""
        for el in el_map[key]:
            elstr += "%s "%(el)
        elstr = elstr[:-1]
        with open("r2t-%s.in"%(key),"w") as fout:
            fout.write('''-function=MeanSquaredDisplacement
-calculation_type=atom
-trajectory_file_path=%s/md-traj-kabsch.xyz
-output_file_path=%s/r2t-%s.txt
-start_frame=300
-end_frame=2000
-frame_interval=1
-dimension=3
-number_of_frames_to_average=1262
-time_scale_type=log
-trajectory_delta_time=0.005
-time_interval=1.5
-number_of_time_points=16
#-time_array_indices=
-atom_name_1=%s'''%(folder[i],folder[i],key,elstr))
        os.system('$ll_path -flagfile=r2t-%s.in'%(key))
        break
foutn.close()
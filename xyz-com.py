import numpy as np
import os
import sys

mass_map = {'N':14.0067,'C':12.011,'O':15.9994,'S':15.9994,'H':1.00797}
res_name_map = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

def get_seq_atom_map(grofile):
    #argument
    #   grofile (string):              gro filename
    #return
    #   res (string):                  protein seq
    #   res_atom_indices (list[list]): atom indices of each amino acid
    #   atom_mass (numpy 1D array):    atom mass
    res = ""              
    res_atom_indices = [] 
    atom_mass = []        
    now_res_id = -1
    with open(grofile, "r") as fin:
        fin.readline()
        natom = int(fin.readline())
        for i in range(natom):
            linelist = fin.readline().strip().split()
            res_id = int(linelist[0][:-3])
            atom_mass.append(mass_map[linelist[1][0]])
            if res_id != now_res_id:
                res += res_name_map[linelist[0][-3:]]
                res_atom_indices.append([i])
                now_res_id = res_id
            else:
                res_atom_indices[-1].append(i)
    atom_mass = np.array(atom_mass)
    return res, res_atom_indices, atom_mass

def write_xyz(xyzfile, coord, element):
    #argument
    #   xyzfile (string):       xyz filename
    #   coord (numpy 3D array): coord, (n_frame, n_atom/n_res, 3)
    #   element (list/string):  atom type/res type
    with open(xyzfile,'w') as fout:
        for i_frame in range(len(coord)):
            fout.write("%d\n\n"%(len(coord[i_frame])))
            for i_atom in range(len(coord[i_frame])):
                fout.write("%s %.3f %.3f %.3f\n"%(element[i_atom],coord[i_frame][i_atom][0],coord[i_frame][i_atom][1],coord[i_frame][i_atom][2]))

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

def compute_com(coord, res, res_atom_indices, atom_mass):
    #argument
    #   coord (numpy 3D array):        coord, (n_frame, n_atom, 3)
    #   res (string):                  protein seq
    #   res_atom_indices (list[list]): atom indices of each amino acid
    #   atom_mass (numpy 1D array):    atom mass
    #return
    #   coord_com (numpy 3D array):    com coord, (n_frame, n_res, 3)
    coord_trans = np.transpose(coord, axes = (0,2,1))
    coord_trans *= atom_mass
    coord_com = []
    for res_indices in res_atom_indices:
        res_mass = np.sum(atom_mass[res_indices])
        coord_com.append(np.sum(coord_trans[:,:,res_indices], axis = 2)/res_mass)
    coord_com = np.transpose(np.array(coord_com), axes = (1,0,2))
    return coord_com

def get_atom_type_indices_map(element):
    #argument
    #   element (list):                     atom type
    #return
    #   atom_type_indices_map (dictionary): atom indices of each element
    atom_type_indices_map = {\
    'N': [i for i, x in enumerate(element) if x[0] == "N"],
    'C': [i for i, x in enumerate(element) if x[0] == "C"],
    'O': [i for i, x in enumerate(element) if x[0] == "O"],
    'S': [i for i, x in enumerate(element) if x[0] == "S"],
    'H': [i for i, x in enumerate(element) if x[0] == "H"],
    'N-bb': [i for i, x in enumerate(element) if x == "N"],
    'C-bb': [i for i, x in enumerate(element) if x == "C" or x == "CA"]}
    return atom_type_indices_map

# p = os.popen("find %s -type %s -name %s"%('.','d','\"*\"'))
# folder = p.read().split('\n')
# folder.pop()
# del folder[0]
# folder = sorted(folder)
# for i in range(len(folder)): folder[i] = folder[i].replace('./', "")
# fout = open("msd-1ps.txt","w")
# fout.write('mutation ')
# for sufix in ['N','C','O','S','H','N-bb','C-bb']: fout.write("%10s "%(sufix))
# fout.write('\n')
# foutstd = open("msd-1ps-std.txt","w")
# foutstd.write('mutation ')
# for sufix in ['N','C','O','S','H','N-bb','C-bb']: foutstd.write("%10s "%(sufix))
# foutstd.write('\n')

folder = sys.argv[1:]
dr_list = []
for a_folder in folder:
    print(a_folder)
    coord, element = read_xyz(a_folder + '/md-traj-kabsch.xyz')
    res, res_atom_indices, atom_mass = get_seq_atom_map(a_folder + '/mutation_raw.gro')
    coord_com = compute_com(coord, res, res_atom_indices, atom_mass)
    write_xyz(a_folder + '/md-traj-kabsch-com.xyz', coord_com, res)




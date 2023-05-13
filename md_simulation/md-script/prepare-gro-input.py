import os
import sys

def create_new_workspace(pdbfile):
    workspace = pdbfile[:pdbfile.index('_')] + '/'
    if os.path.isdir(workspace) == False:
        os.system('mkdir -p %s'%(workspace))
        os.chdir('%s'%(workspace))
        return True
    return False

def back_to_root(cwd):
    os.chdir(cwd)

def make_gromacs_input(pdbfile_dir, gmx_executable, pdbfile):
    gmx_command = 'pdb2gmx'
    gro_file = 'mutation_raw'
    command = '%s %s -f %s%s -o %s.gro -water spce -ignh < ../forcefield-option.txt'%(gmx_executable, gmx_command, pdbfile_dir, pdbfile, gro_file)
    os.system(command)
    
    gmx_command = 'editconf'
    command = '%s %s -f %s.gro -o %s.gro -c -d 1.0 -bt cubic'%(gmx_executable, gmx_command, gro_file, gro_file+'_newbox')
    os.system(command)
    gro_file += '_newbox'

    gmx_command = 'solvate'
    command = '%s %s -cp %s.gro -cs spc216.gro -o %s.gro -p topol.top'%(gmx_executable, gmx_command, gro_file, gro_file+'_solvate')
    os.system(command)
    gro_file += '_solvate'
    
    gmx_command = 'grompp'
    command = '%s %s -f ../ions.mdp -c %s.gro -p topol.top -o ions.tpr'%(gmx_executable, gmx_command, gro_file)
    os.system(command)

    gmx_command = 'genion'
    gro_file += '_ions'
    command = '%s %s -s ions.tpr -o %s.gro -p topol.top -pname NA -nname CL -neutral < ../ions-option.txt'%(gmx_executable, gmx_command, gro_file)
    os.system(command)

cwd = os.getcwd()
gmx_executable = '/home/scl6/zlab/gromacs/bin/gmx_mpi'
pdbfile_dir = '/home/scl6/zlab/novozymes-enzyme-stability-prediction/3d/pdb_from_archive/'
p = os.popen("find %s -type %s -name %s"%(pdbfile_dir,'f','*pdb'))
pdbfile_list = p.read().split('\n')
pdbfile_list.pop()
pdbfile_list = sorted(pdbfile_list)
for i in range(len(pdbfile_list)): pdbfile_list[i] = pdbfile_list[i].replace(pdbfile_dir, "")

for pdbfile in pdbfile_list:
    if create_new_workspace(pdbfile):
        make_gromacs_input(pdbfile_dir, gmx_executable, pdbfile)
        back_to_root(cwd)
    else:
        print("A working space for %s already exist."%(pdbfile))
    
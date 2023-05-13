import os
import sys

def make_slurm(folder, slurm_sufix, slurm_header_command, run_command): 
    with open("run_gro_%s.slurm"%(slurm_sufix), "w") as fout:
        fout.write(slurm_header_command)
        for a_folder in folder:
            if os.path.isfile(a_folder+'/md-hot-traj.xyz') == True:
                print('Traj file for %s already exist.'%(a_folder))
                continue
            fout.write("cd %s/\n"%(a_folder))
            fout.write(run_command)
            fout.write("cd ../\n\n")


p = os.popen("find %s -type %s -name %s"%('.','d','\"*\"'))
folder = p.read().split('\n')
folder.pop()
del folder[0]
folder = sorted(folder)
for i in range(len(folder)): folder[i] = folder[i].replace('./', "")
folder.remove('wild')

# slurm_header_command = '''#!/bin/bash
# #
# #SBATCH --nodes=1
# #SBATCH --ntasks-per-node=8
# #SBATCH --cpus-per-task=8
# #SBATCH --time=100:00:00
# #SBATCH --partition=zlab
# #SBATCH --mail-user=scl6@illinois.edu
# #SBATCH --mail-type=FAIL
# #SBATCH --job-name=j-mu
# ##################################### 
# cd $SLURM_SUBMIT_DIR
# module load openmpi/4.1.0-gcc-7.2.0

# '''
# run_min_md_command = '''#echo "startrunning"
# #date
# #$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
# #mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
# echo "startrunning"
# date
# $gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
# mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
# date
# $gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
# python ../gro2xyz.py md-traj.gro
# date
# '''
# make_slurm(folder[:len(folder)//2], '0', slurm_header_command, run_min_md_command)
# make_slurm(folder[len(folder)//2:], '1', slurm_header_command, run_min_md_command)

slurm_header_command = '''#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=8
#SBATCH --time=100:00:00
#SBATCH --partition=zlab
#SBATCH --mail-user=scl6@illinois.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=j-mu
##################################### 
cd $SLURM_SUBMIT_DIR
module load openmpi/4.1.0-gcc-7.2.0

'''
run_min_md_command = '''#echo "startrunning"
echo "startrunning"
date
$gmx_path grompp -f ../md-hot.mdp -c md.gro -p topol.top -o md-hot.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md-hot
date
$gmx_path trjconv -f md-hot.trr -s md-hot.tpr -o md-hot-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-hot-traj.gro
date
'''
make_slurm(folder[:len(folder)//2], 'hot-0', slurm_header_command, run_min_md_command)
make_slurm(folder[len(folder)//2:], 'hot-1', slurm_header_command, run_min_md_command)

# slurm_header_command = '''#!/bin/bash
# #
# #SBATCH --nodes=1
# #SBATCH --ntasks-per-node=1
# #SBATCH --cpus-per-task=64
# #SBATCH --time=60:00:00
# #SBATCH --partition=zlab
# #SBATCH --mail-user=scl6@illinois.edu
# #SBATCH --mail-type=FAIL
# #SBATCH --job-name=j-mu
# ##################################### 
# cd $SLURM_SUBMIT_DIR
# module load anaconda/3
# module load openmpi/4.1.0-gcc-7.2.0

# '''
# run_trr2gro2xyz_command = '''echo "startrunning"
# date
# $gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
# python ../gro2xyz.py md-traj.gro
# date
# '''
# make_slurm(folder, 'traj-all', slurm_header_command, run_trr2gro2xyz_command)

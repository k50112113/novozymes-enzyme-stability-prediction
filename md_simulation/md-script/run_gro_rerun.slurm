#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=8
#SBATCH --time=60:00:00
#SBATCH --partition=zlab
#SBATCH --mail-user=scl6@illinois.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=j-mu
##################################### 
cd $SLURM_SUBMIT_DIR
module load openmpi/4.1.0-gcc-7.2.0

#at 1051 line of mutation_raw_newbox_solvate_ions in D197N, change to:
#   71ARG   HH11 1049   3.775   3.206   4.131

# cd D197N/
# echo "startrunning"
# date
# $gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
# mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
# echo "startrunning"
# date
# $gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
# mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
# date
# cd ../

#at 639 line of mutation_raw_newbox_solvate_ions in V127R, change to:
#   45VAL   HG21  637   4.517   3.489   4.039

# cd V127R/
# echo "startrunning"
# date
# $gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
# mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
# echo "startrunning"
# date
# $gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
# mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
# date
# cd ../


cd N63R/
# echo "startrunning"
# date
# $gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
# mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
# echo "startrunning"
# date
# $gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
cd ../
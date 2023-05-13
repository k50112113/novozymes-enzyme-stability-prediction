#!/bin/bash
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

cd A101E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A101F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A101G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A101H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A101I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A101K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A101N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A101S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A101T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A101W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A101Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A111/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A111C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A111G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A111I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A111K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A111P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A111R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A111S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A111V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A113C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A113D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A113E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A113F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A113H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A113I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A113M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A113N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A113P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A113Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A113S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A113V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A113W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A122F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A122G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A122K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A122M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A122N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A122P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A122V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A122W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A122Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A125/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A125C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A125D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A125E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A125F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A125G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A125I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A125L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A125M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A125N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A125P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A125S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A125Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A130/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A130D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A130E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A130F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A130G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A130H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A130K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A130L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A130M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A130N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A130P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A130Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A130R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A130S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A130T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A130V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A130W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A130Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A145/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A145C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A145D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A145E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A145G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A145I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A145L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A145P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A145Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A145T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A145V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A145W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A145Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A149E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A149F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A149G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A149P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A149R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A149T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A16C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A16D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A16F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A16G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A16H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A16I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A16L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A16N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A16P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A16W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A172C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A172D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A172E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A172G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A172H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A172K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A172L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A172M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A172N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A172Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A172R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A172S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A172T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A172V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A172W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A172Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A173/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A173C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A173D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A173E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A173F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A173G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A173M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A173P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A179/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A179F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A179G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A179H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A179K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A179M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A179N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A179Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A179R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A179S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A195C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A195D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A195E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A195F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A195G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A195H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A195I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A195K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A195L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A195M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A195N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A195P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A195Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A195V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A195Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A203L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A203N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A206F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A206G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A206H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A206I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A206L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A206M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A206Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A206S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A206T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A206V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A206W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A31/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A31C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A31E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A31F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A31G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A31I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A31K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A31M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A31P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A31Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A31S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A31Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A40C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A40F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A40H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A40K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A40P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A40Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A40S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A40T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A49/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A49C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A49D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A49E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A49F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A49G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A49K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A49L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A49N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A49P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A49Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A49S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A49W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A49Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A55C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A55D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A55E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A55F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A55G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A55H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A55I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A55K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A55M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A55P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A55R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A55S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A55T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A55Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A76E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A76F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A76I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A76L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A76M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A76Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A76R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A76S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A76T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A76Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A91/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A91C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A91D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A91E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A91F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A91I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A91K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A91L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A91P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A91Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A91R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A91S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A91V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A91W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D103A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D103C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D103E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D103F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D103M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D103N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D103P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D103Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D103R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D103S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D103T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D103W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D103Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D152/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D152C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D152E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D152F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D152G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D152L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D152M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D152N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D152P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D152R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D152T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D152V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D152Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D184A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D184E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D184I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D184K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D184N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D184Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D184T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D184V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D184Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D189A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D189C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D189G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D189H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D189K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D189L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D189Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D189S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D189T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D189W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D197A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D197C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D197E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D197F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D197G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D197H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D197I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D197K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D197M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D197N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D197P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D197S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D197T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D197V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D197W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D197Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D201/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D201L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D210H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D210K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D210N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D23/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D23A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D23E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D23F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D23G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D23H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D23L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D23M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D23N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D23P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D23Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D23T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D23Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D27/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D27C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D27F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D27G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D27H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D27I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D27L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D27M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D27N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D27P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D27Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D32A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D32C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D32H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D32K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D32M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D32N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D32P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D32Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D32S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D32W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D46A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D46C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D46E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D46F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D46G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D46H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D46I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D46K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D46M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D46N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D46P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D46Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D46S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D46T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D46V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D46W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D46Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D72A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D72F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D72G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D72K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D72L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D72M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D72N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D72V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D72W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D83A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D83C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D83E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D83F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D83H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D83I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D83L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D83N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D83S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D83T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd D83V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E108A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E108C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E108F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E108G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E108H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E108K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E108L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E108M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E108Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E108R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E108S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E108T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E108V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E108W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E108Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E109A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E109C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E109D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E109F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E109G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E109H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E109I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E109M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E109N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E109P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E109Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E109S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E109T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E109V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E109W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E109Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E211D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E211F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E211H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E211L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E211P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E211S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E211T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E211Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E34A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E34D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E34H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E34I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E34L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E34M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E34N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E34P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E34Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E64/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E64A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E64C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E64D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E64F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E64G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E64H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E64I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E64K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E64M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E64N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E64S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E64T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E64V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E64W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd E64Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F104I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F104M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F104N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F104P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F104W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F104Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F112D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F112G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F112H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F112L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F112M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F112P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F112Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F112R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F112S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F112V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F112W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F112Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F142A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F142C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F142E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F142G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F142M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F142N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F142R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F142S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F142V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F142W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F142Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F151A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F151C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F151D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F151E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F151G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F151I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F151K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F151M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F151N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F151P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F151Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F151S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F151V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F151W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F151Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F164/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F164A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F164C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F164D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F164E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F164H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F164I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F164L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F164M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F164N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F164P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F164S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F164T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F164V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F164W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F164Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F169/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F169D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F169E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F169I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F169L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F169N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F169T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F169V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F169W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F200H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F200N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F200W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F200Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F216C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F216I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F216L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F216M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F216Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F43A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F43I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F43L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F43M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F43T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F43Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F81C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F81T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G118E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G118F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G118I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G118L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G118M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G118N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G118P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G118R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G118S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G118V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G118W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G118Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G119A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G119D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G119E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G119F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G119I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G119K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G119L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G119M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G119N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G119P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G119R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G119S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G119T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G119V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G119W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G119Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G136A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G136C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G136D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G136E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G136F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G136I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G136K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G136L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G136M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G136N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G136P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G136Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G136R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G136S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G136T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G136V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G136W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G136Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G137/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G137A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G137C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G137D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G137E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G137H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G137I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G137K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G137L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G137M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G137N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G137P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G137Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G137R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G137S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G137T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G137V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G137W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G137Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G141/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G141A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G141C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G141D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G141I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G141L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G141M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G141N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G141P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G141Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G141R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G141T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G141V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G161/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G161A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G161D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G161H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G161N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G161R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G161W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G161Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G171A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G171D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G171E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G171H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G171K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G171M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G171Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G171R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G174A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G174C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G174E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G174F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G174H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G174I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G174M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G174N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G174P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G174Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G174R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G174T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G174V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G174W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G174Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G181A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G181C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G181D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G181F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G181H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G181I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G181L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G181M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G181N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G181P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G181S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G181V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G181W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G196A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G196D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G196F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G196K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G196L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G196M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G196N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G196P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G196Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G196R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G196S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G196T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G196V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G196W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G196Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G199A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G199C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G199D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G199E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G199F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G199H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G199I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G199K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G199L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G199N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G199Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G199R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G199S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G199T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G199Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G20/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G20D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G20E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G20F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G20H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G20I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G20L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G20M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G20N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G20P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G20Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G20R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G20S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G20T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G20V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G20Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G220A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G220E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G220H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G220N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G220P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G220T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G220V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G220W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G22A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G22C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G22D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G22E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G22F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G22I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G22K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G22L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G22P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G22Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G22R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G22T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G37A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G37C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G37D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G37E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G37F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G37H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G37I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G37K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G37L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G37N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G37P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G37Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G37R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G37S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G37T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G37W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G37Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G54A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G54C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G54E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G54F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G54H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G54I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G54K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G54P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G54S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G54T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G54W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G54Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G75F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G75K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G75M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G75N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G75P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G75Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G75R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G75S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G75Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G79D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G79H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G79N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G79P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G79T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G79V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G79W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G88/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G88A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G88C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G88E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G88F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G88H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G88I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G88M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G88P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G88S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G88T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G88W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd G88Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I123A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I123C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I123D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I123E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I123F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I123G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I123K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I123L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I123M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I123P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I123R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I123S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I123T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I123V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I123W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I166/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I166A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I166C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I166F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I166G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I166H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I166L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I166P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I166Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I166S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I166T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I166Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I194A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I194C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I194K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I194M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I194N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I194R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I194S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I194T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I194V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I194W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I29A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I29M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I29T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I29V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I50C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I50D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I50E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I50F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I50H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I50K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I50L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I50M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I50P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I50Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I50R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I50Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I89/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I89E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I89F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I89G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I89K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I89M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I89Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I89W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I89Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K105/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K105A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K105C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K105D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K105G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K105H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K105I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K105M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K105N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K105P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K105Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K105R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K105T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K105V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K105Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K147C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K147E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K147I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K147N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K147R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K155C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K155D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K155E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K155F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K155H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K155L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K155M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K155Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K155R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K155S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K155T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K155V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K155W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K155Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K160/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K160A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K160C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K160D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K160G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K160H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K160I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K160L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K160P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K160Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K160T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K160V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K168/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K168C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K168E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K168F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K168G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K168H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K168P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K168Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K168R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K168S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K168T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K168V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K178/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K178A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K178C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K178D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K178F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K178G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K178H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K178I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K178L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K178M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K178N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K178R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K178S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K178T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K178V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K178W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K18/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K185/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K185A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K185C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K185D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K185E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K185F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K185G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K185H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K185I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K185L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K185N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K185P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K185R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K185S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K185T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K185V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K185Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K18C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K18F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K18G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K18H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K18I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K18L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K18M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K18P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K18Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K18R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K18S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K18V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K18W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K190/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K190A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K190C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K190D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K190F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K190G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K190P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K190R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K190S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K190W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K190Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K192A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K192C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K192D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K192F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K192G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K192H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K192I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K192L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K192N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K192P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K192Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K192T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K192V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K192W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K204A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K204C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K204D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K204E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K204F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K204G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K204I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K204L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K204N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K204Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K204S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K204T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K204V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K204W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K204Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K212/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K212C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K212D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K212E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K212F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K212G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K212H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K212I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K212L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K212M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K212N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K212P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K212T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K212V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K212W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K212Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K215/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K215D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K215E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K215G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K215I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K215L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K215P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K215Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K215T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K215V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K215Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K221C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K221D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K221F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K221G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K221H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K221I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K221N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K221P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K221Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K221S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K221T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K221Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K30/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K30A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K30C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K30E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K30H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K30M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K30R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K30V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K30W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K30Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K36/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K36D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K36E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K36G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K36H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K36I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K36L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K36M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K36N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K36P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K36Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K36R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K36S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K36T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K36V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K36W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K36Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K53E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K53F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K53M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K53P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K53Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K53R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K53S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K53V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K53Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K65A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K65C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K65D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K65E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K65F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K65G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K65I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K65L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K65M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K65Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K65R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K65S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K65T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K65V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K65W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K65Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K67/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K67I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K67P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K67R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K67S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K67W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K78C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K78F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K78I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K78P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K78R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K78S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K78Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K82/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K82A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K82E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K82F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K82G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K82H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K82I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K82L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K82M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K82N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K82P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K82R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K82T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K82V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K82W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K82Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K86/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K86A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K86C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K86D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K86E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K86F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K86H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K86L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K86M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K86N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K86Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K86S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K86T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K86V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K86Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K90/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K90A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K90C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K90D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K90H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K90I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K90Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K90R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K90S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K90T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K90V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K95A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K95C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K95D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K95F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K95G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K95H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K95I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K95N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K95P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K95Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K95R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K95S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K95W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L116A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L116C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L116D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L116E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L116F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L116G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L116K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L116M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L116N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L116R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L116S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L116T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L116W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L116Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L124/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L124A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L124D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L124E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L124G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L124K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L124M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L124N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L124P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L124S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L124T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L124V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L124W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L124Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L129A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L129C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L129D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L129F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L129G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L129I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L129K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L129M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L129N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L129P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L129R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L129S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L129V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L129W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L129Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L139/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L139A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L139C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L139E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L139F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L139I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L139K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L139M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L139N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L139P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L139Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L139R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L139T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L139W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L139Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L17/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L17E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L17K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L17P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L17Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L17T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L17V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L17W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L17Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L180A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L180C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L180D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L180H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L180I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../


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

cd L180K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L180M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L180P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L180R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L180S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L180V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L180W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L33F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L33I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L33S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L33V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L33W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L33Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L41A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L41I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L41M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L41N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L41T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L51A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L51C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L51D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L51E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L51F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L51G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L51H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L51M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L51N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L51P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L51S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L51T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L51V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L51W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L59C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L59D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L59E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L59F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L59H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L59P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L59Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L59S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L59T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L59V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L59W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L92/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L92A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L92C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L92E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L92F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L92G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L92N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L92P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L92S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L92T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L92Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N121C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N121D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N121F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N128A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N128C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N128E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N128F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N128G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N128I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N128K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N128L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N128M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N128P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N128Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N128R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N128S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N128T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N128V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N128W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N128Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N133/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N133A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N133C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N133D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N133E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N133F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N133G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N133H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N133I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N133L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N133M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N133P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N133Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N133R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N133S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N133T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N133V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N133W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N133Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N140/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N140A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N140C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N140E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N140F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N140G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N140H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N140K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N140L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N140P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N140Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N140R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N140V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N140Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N146A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N146C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N146H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N146K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N146L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N146M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N146P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N146R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N146S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N146T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N183D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N183F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N183H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N183I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N183L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N183R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N183S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N183T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N183V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N183W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N183Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N191/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N191C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N191G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N191H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N191K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N191L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N191M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N191P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N191Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N193C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N193D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N193E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N193G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N193H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N193I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N193K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N193L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N193P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N193Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N193S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N193T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N193V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N213A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N213C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N213D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N213E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N213F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N213G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N213H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N213I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N213K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N213L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N213M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N213P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N213Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N213R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N213S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N213V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N213W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N213Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N214/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N214A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N214D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N214E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N214G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N214H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N214I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N214K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N214Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N214R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N214S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N214T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N214V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N214Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N217A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N217D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N217E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N217H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N217I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N217K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N217L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N217M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N217P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N217Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N217R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N217S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N217T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N217V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N217W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N217Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N57A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N57C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N57E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N57F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N57G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N57H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N57K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N57L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N57M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N57P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N57Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N57R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N57S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N57V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N57W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N57Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N63/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N63C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N63E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N63F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N63G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N63H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N63I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N63K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N63L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N63M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N63P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N63Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N63R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N63S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N63T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N63V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N63Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N68/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N68D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N68E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N68H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N68I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N68L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N68M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N68P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N68Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N70C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N70D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N70K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N70L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N70M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N70Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N70R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N70S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N70T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N70V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N70W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N70Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N77A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N77E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N77F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N77H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N77I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N77K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N77L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N96A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N96C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N96D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N96E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N96F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N96G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N96I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N96M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N96P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N96Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N96S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N96V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N96W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P107A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P107C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P107D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P107E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P107F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P107G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P107H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P107I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P107L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P107M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P107N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P107Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P107S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P107V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P107W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P107Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P126I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P126K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P126M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P126N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P126R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P126S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P126T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P126V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P126Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P153A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P153C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P153D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P153E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P153G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P153I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P153K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P153L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P153N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P153Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P153V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P153W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P156A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P156C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P156D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P156H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P156N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P156S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P156T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P156V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P156Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P175/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P175A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P175C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P175E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P175F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P175G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P175H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P175I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P175K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P175M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P175N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P175Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P175R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P175S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P175T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P175V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P175Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P202A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P202C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P202E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P202F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P202G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P202H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P202I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P202L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P202M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P202N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P202Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P202S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P202V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P202Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P28/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P28A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P28C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P28F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P28H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P28I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P28L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P28M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P28N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P28Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P28R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P28S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P28V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P28W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P28Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P42E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P42G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P42H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P42I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P42K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P42L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P42M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P42R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P42S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P42T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P42V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P56A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P56C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P56D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P56E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P56H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P56I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P56K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P56L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P56N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P56R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P56S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P56V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P56W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P56Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P80/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P80C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P80H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P80Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P80T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P80Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P84A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P84D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P84E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P84F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P84G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P84H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P84M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P84N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P84Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P84R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P84S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P84T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P84Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P93A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P93C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P93H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P93I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P93M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P93Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P93S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P93V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P93Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P94/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P94A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P94C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P94D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P94E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P94I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P94L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P94Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P94S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P94T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P94V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P94W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P97E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P97F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P97G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P97N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P97Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P97R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P97S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P97T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P97V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P97W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q102E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q102G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q102L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q102M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q102N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q102P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q102R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q102T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q102V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q102W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q102Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q117C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q117D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q117G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q117L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q117M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q117N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q117P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q117W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q132C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q132D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q132E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q132F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q132G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q132H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q132I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q132K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q132L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q132M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q132N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q132P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q132R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q132S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q132T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q132V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q132W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q132Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q135/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q135A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q135C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q135D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q135E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q135F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q135G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q135H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q135I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q135K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q135L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q135M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q135N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q135P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q135R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q135S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q135T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q135V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q135W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q135Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q150/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q150D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q150E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q150F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q150G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q150H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q150K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q150L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q150N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q150P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q150S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q150T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q150V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q150W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q150Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q157A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q157C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q157D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q157E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q157F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q157G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q157I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q157K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q157L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q157M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q157P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q157S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q157W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q158/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q158D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q158F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q158G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q158H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q158K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q158L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q158M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q158N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q158P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q158R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q158S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q158T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q158V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q158Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q165/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q165A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q165C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q165H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q165K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q165L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q165M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q165N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q165R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q165Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q208C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q208E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q208F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q208I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q208L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q208M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q208N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q208R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q208S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q208T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q208V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q208W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q208Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q25/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q25A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q25D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q25E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q25F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q25G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q25H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q25K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q25L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q25N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q25P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q25R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q25S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q25W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q25Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q38/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q38A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q38C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q38D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q38E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q38F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q38G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q38H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q38I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q38K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q38M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q38N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q38P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q38R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q38S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q38T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q38V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q38W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q38Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q60A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q60D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q60E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q60F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q60G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q60I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q60K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q60L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q60M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q60N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q60P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q60R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q60S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q60T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q60V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q60W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q85/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q85A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q85C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q85D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q85E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q85F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q85G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q85H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q85K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q85N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q85R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q85S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q85T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q85V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q85W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Q85Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R61A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R61C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R61D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R61E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R61F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R61G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R61H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R61I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R61K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R61L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R61M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R61P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R61S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R61T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R61V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R61W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R61Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R71C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R71D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R71E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R71F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R71G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R71H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R71K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R71L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R71M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R71N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R71S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R71T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R71V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R71W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R71Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R73D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R73E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R73N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R73P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R73Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R73S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R73V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R73W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd R73Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S100A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S100C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S100F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S100G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S100H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S100I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S100L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S100M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S100P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S100Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S100R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S100T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S100V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S100W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S106A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S106C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S106D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S106E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S106F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S106H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S106I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S106K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S106L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S106M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S106P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S106Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S106R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S106T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S106V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S106W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S106Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S114F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S114K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S114M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S114P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S114T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S114V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S114W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S115C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S115D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S115E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S115F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S115I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S115K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S115L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S115M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S115N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S115R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S115T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S115V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S115W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S131/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S131A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S131C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S131D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S131E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S131F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S131G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S131H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S131I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S131K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S131L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S131N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S131P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S131Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S131R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S131T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S131V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S131W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S134A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S134C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S134D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S134E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S134F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S134G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S134H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S134I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S134K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S134L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S134M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S134N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S134P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S134Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S134R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S134T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S134V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S134W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S134Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S144/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S144C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S144E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S144F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S144G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S144H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S144I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S144L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S144M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S144Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S144V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S144Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S154/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S154C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S154E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S154F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S154H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S154I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S154M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S154N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S154P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S154R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S154W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S182A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S182D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S182E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S182F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S182G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S182H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S182I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S182K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S182L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S182M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S182N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S182P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S182Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S182R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S182T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S182V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S182W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S182Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S186A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S186D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S186E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S186F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S186G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S186H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S186K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S186L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S186P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S186T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S186V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S186W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S21A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S21D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S21E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S21F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S21G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S21H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S21I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S21K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S21L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S21N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S21P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S21Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S21R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S24/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S24A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S24C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S24E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S24G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S24H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S24I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S24K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S24L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S24M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S24N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S24P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S24Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S24R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S24T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S24W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S24Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S26/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S26A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S26C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S26E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S26F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S26G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S26H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S26N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S26Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S26R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S26V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S26W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S26Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S39/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S39A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S39C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S39D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S39E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S39F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S39G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S39I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S39M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S39N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S39P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S39Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S39R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S39T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S39V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S39W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S69A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S69C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S69F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S69H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S69K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S69M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S69T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S69V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S69W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S74/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S74C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S74D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S74F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S74G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S74L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S74N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S74Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S74Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S98A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S98D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S98E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S98F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S98G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S98I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S98Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S98R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S98T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd S98V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T120C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T120E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T120F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T120K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T120L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T120N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T120P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T159/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T159E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T159F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T159H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T159I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T159K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T159L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T159M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T159N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T159Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T159S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T159V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T159W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T159Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T162/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T162A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T162C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T162D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T162E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T162F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T162I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T162K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T162M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T162P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T162Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T162R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T162S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T162V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T162W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T162Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T167/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T167C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T167D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T167E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T167H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T167N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T167P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T167Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T167S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T170/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T170A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T170D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T170E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T170G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T170I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T170L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T170N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T170Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T170Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T19/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T19A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T19C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T19D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T19E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T19G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T19H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T19I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T19K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T19M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T19N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T19P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T19Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T19R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T19V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T19Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T66A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T66E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T66F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T66G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T66K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T66M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T66V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T66W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V127A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V127C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V127D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V127E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V127F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V127G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V127I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V127L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V127M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V127R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V127S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V127T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V127W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V127Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V138/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V138A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V138C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V138D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V138E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V138F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V138G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V138H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V138I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V138K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V138L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V138M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V138N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V138P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V138Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V138R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V138S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V138T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V138W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V138Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V148I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V148W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V148Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V187A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V187C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V187E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V187F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V187G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V187H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V187I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V187L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V187M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V187N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V187P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V187R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V187T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V187W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V187Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V219A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V219C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V219D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V219F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V219G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V219H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V219I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V219K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V219L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V219M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V219N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V219P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V219Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V219R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V219T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V219W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V219Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V35A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V35T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V45A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V45E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V45G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V45I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V45K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V45L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V45M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V45Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V45W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V58/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V58D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V58E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V58G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V58H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V58K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V58L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V58N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V58P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V58T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V58W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V58Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V62/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V62A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V62C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V62D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V62F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V62G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V62H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V62I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V62N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V62P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V62Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V62R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V62S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V62T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd V62Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W163/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W163C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W163D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W163E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W163F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W163H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W163I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W163K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W163M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W163P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W163Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W163R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W163S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W163T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W163V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W163Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W198A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W198D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W198H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W198I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W198L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W198M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W198S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W198T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W198V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W205A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W205C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W205D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W205E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W205F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W205G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W205I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W205K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W205L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W205N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W205S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W205T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W205V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W48E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W48F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W48L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W48M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W48Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W87C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W87E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W87G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W87I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W87K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W87M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W87N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W87R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W87V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W87Y/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd W99F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd WT/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y110F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y110W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y143/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y143A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y143C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y143D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y143E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y143G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y143H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y143I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y143K/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y143L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y143M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y143N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y143P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y143Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y143T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y143V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y143W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y176/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y176A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y176C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y176E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y176F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y176I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y176N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y176P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y176S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y176V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y176W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y207A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y207C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y207E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y207F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y207G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y207H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y207I/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y207L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y207M/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y207P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y207Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y207R/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y207S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y207T/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y207V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y207W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y209C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y209E/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y209F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y209H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y209L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y209N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y209Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y209V/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y209W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y218A/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y218C/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y218D/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y218F/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y218G/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y218H/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y218L/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y218N/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y218P/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y218Q/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y218S/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y218W/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd wild/
#echo "startrunning"
#date
#$gmx_path grompp -f ../min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr
#mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm min-sd
echo "startrunning"
date
$gmx_path grompp -f ../md.mdp -c min-sd.gro -p topol.top -o md.tpr
mpirun -np 8 $gmx_path mdrun -v -ntomp 8 -deffnm md
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../


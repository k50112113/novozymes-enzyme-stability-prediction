#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=60:00:00
#SBATCH --partition=zlab
#SBATCH --mail-user=scl6@illinois.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=j-mu
##################################### 
cd $SLURM_SUBMIT_DIR
module load anaconda/3
module load openmpi/4.1.0-gcc-7.2.0

cd A101E/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A101F/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A101H/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A101I/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A101K/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A101S/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A101T/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A101W/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A101Y/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A111/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A111G/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A111I/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A111K/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A111P/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A111R/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A111S/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A111V/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A113C/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A113D/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A113E/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A113F/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A113H/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A113I/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A113M/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A113N/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A113P/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A130/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A145/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A179/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A31/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A40K/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A40P/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A40Q/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A40S/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A40T/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A49/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A49D/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A49F/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A49G/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A49K/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A49L/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A49S/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A49W/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A49Y/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A55C/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A55D/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A55E/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A55G/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A55H/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A55I/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A55K/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A55P/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd A91/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F151A/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F151C/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F151D/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F151E/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F151G/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F151I/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F151K/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F151M/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F151N/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F151P/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F151S/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F151V/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F151Y/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F164A/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F164D/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F164H/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F164I/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F164L/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd F164M/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I194W/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I29A/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I29T/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I50D/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I50E/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I50F/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I50H/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I50M/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I50P/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I50Q/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I89E/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I89F/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I89G/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I89K/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I89Q/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I89W/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd I89Y/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K105A/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K105D/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K147C/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K147E/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K147I/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K147N/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K147R/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K155C/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K155D/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K155F/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K155H/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K155L/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K155R/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K155S/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K155T/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K155V/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K155Y/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K160C/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K160D/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K160H/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K190/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K30/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K82/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd K90/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L116D/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L116K/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L116N/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L116R/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd L116S/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N121D/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N121F/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N128C/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N128E/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N128F/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N128M/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N128Q/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N128R/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N128S/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N128T/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N128V/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N128W/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N128Y/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N133/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N133C/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N133F/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N133L/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N133Q/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N133R/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N133S/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N133T/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N133V/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N133W/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N140/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N140C/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N140E/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N140F/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N140G/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N140K/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N140P/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N140R/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N146C/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N146L/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N146M/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N146S/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N146T/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N183D/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N183F/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N183H/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N183L/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N183S/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N191/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N191C/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N214/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd N63/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P107D/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P107F/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P107G/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P107H/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P107I/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P107M/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P107N/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P107Q/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P107S/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P107V/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P107Y/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P126M/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P126S/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P126T/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P126Y/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P153I/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P153K/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P153N/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P153W/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P156T/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P175E/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P175F/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P175G/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P175H/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P175N/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P175Q/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P175S/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P175T/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P28I/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P28L/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P28M/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P28N/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P28Q/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P28R/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P28S/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P28V/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P28W/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P28Y/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P42E/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P42G/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P42H/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P42I/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P42K/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P42L/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd P42M/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T162C/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T162F/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T162I/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T162M/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T162P/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T162R/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T162S/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T162V/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T162W/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T162Y/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T167C/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T167D/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T167H/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T167N/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T167Q/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd T170A/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y143L/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y143P/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y143Q/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd Y143V/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../

cd wild/long/
echo "startrunning"
date
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center < ../trjconv-option.txt
python ../gro2xyz.py md-traj.gro
date
cd ../


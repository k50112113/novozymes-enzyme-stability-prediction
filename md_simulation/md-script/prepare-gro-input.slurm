#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --partition=secondary
#SBATCH --mail-user=scl6@illinois.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=j-mut
##################################### 
cd $SLURM_SUBMIT_DIR
module load openmpi/4.1.0-gcc-7.2.0
module load anaconda/3

python prepare-gro-input.py
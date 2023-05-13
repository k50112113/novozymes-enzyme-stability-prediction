import os
import sys

p = os.popen("find %s -type %s -name %s"%('.','d','\"*\"'))
folder = p.read().split('\n')
folder.pop()
del folder[0]
folder = sorted(folder)
for i in range(len(folder)): folder[i] = folder[i].replace('./', "")

folder_batch = [37]*31 + [38]*(63-31)
ss = sum(folder_batch)
folder_batch.append(len(folder)-ss)

# folder_batch = [100]*23
# ss = sum(folder_batch)
# folder_batch.append(len(folder)-ss)

print(folder_batch)

slurm_header_command = '''#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --partition=zlab
#SBATCH --mail-user=scl6@illinois.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=j-kab
##################################### 
cd $SLURM_SUBMIT_DIR
module load openmpi/4.1.0-gcc-7.2.0
'''
with open("run_kabsch.slurm", "w") as fout:
    j_file = 0
    fout.write(slurm_header_command)
    for i_batch, batch_size in enumerate(folder_batch):
        fout.write('''python compute-kabsch.py ''')
        for x in range(batch_size):
            fout.write('%s/md-hot-traj.xyz '%(folder[j_file]))
            j_file += 1
        if i_batch < len(folder_batch) - 1: fout.write('''&''')
        fout.write(''' \n''')
# with open("run_dr.slurm", "w") as fout:
#     j_file = 0
#     fout.write(slurm_header_command)
#     for i_batch, batch_size in enumerate(folder_batch):
#         fout.write('''python xyz-dr.py ''')
#         for x in range(batch_size):
#             fout.write('%s/md-traj-kabsch.xyz '%(folder[j_file]))
#             j_file += 1
#         if i_batch < len(folder_batch) - 1: fout.write('''&''')
#         fout.write(''' \n''')
with open("run_com.slurm", "w") as fout:
    j_file = 0
    fout.write(slurm_header_command)
    for i_batch, batch_size in enumerate(folder_batch):
        fout.write('''python compute-com.py ''')
        for x in range(batch_size):
            fout.write('%s/md-hot-traj-kabsch.xyz '%(folder[j_file]))
            j_file += 1
        if i_batch < len(folder_batch) - 1: fout.write('''&''')
        fout.write(''' \n''')
    
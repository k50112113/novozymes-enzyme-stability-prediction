#module load openmpi/4.1.0-gcc-7.2.0

#make sure a correct mpirun is in your PATH
export $gmx_path=path/to/gmx_mpi or path/to/gmx
#<protein>.pdb = path/to/protein.pdb

#prepare input for minimization and MD simulation using CHARMM forcefield and SPC/E water model.
$gmx_path pdb2gmx -f <protein>.pdb -o mutation_raw.gro -water spce -ignh                                        # convert *.pdb to *gro (select 8):  <protein>.pdb                                               -> mutation_raw.gro, topol.top 
$gmx_path editconf -f mutation_raw.gro -o mutation_raw_newbox.gro -c -d 1.0 -bt cubic                           # add a box:                         mutation_raw.gro                                            -> mutation_raw_newbox.gro
$gmx_path solvate -cp mutation_raw_newbox.gro -cs spc216.gro -o mutation_raw_newbox_solvate.gro -p topol.top    # add water:                         mutation_raw_newbox.gro, spc216.gro, topol.top              -> mutation_raw_newbox_solvate.gro
$gmx_path grompp -f ions.mdp -c mutation_raw_newbox_solvate.gro -p topol.top -o ions.tpr                        # precompile a run for adding ions:  mutation_raw_newbox_solvate.gro, ions.mdp, topol.top        -> ions.tpr
$gmx_path genion -s ions.tpr -o mutation_raw_newbox_solvate_ions.gro -p topol.top -pname NA -nname CL -neutral  # add ions (select 13):              ions.tpr, topol.top                                         -> mutation_raw_newbox_solvate_ions.gro
#run minimization (relax the protein to a conformational local minimum).
$gmx_double_path grompp -f min-sd.mdp -c mutation_raw_newbox_solvate_ions.gro -p topol.top -o min-sd.tpr        # precompile a run for minimization: mutation_raw_newbox_solvate_ions.gro, min-sd.mdp, topol.top -> min-sd.tpr
mpirun -np 8 $gmx_double_path mdrun -v -ntomp 8 -deffnm min-sd                                                  # run minimization:                  min-sd.tpr                                                  -> min-sd.gro             
#run MD simulation.
$gmx_double_path grompp -f md.mdp -c min-sd.gro -p topol.top -o md.tpr                                          # precompile a run for MD:           min-sd.gro, min.mdp, topol.top                              -> md.tpr
mpirun -np 8 $gmx_double_path mdrun -v -ntomp 8 -deffnm md                                                      # run MD:                            md.tpr                                                      -> md.trr             
#extract trajectory only for protein.
$gmx_path trjconv -f md.trr -s md.tpr -o md-traj.gro -center                                                    # extract (select 1, and 1):         md.trr, md.tpr                                              -> md-traj.gro

'''
Created on 1.10.2015

@author: tohekorh
'''
import os

w       =   13

ratios  =   [12] #range(10,40)
ncores  =   8
edge    =   'ac'

path    =   '/space/tohekorh/ShearSlide/files/scripts/'

for ratio in ratios:
    script = \
'''#!/bin/bash -l
# created: Feb 6, 2015 1:46 PM
# author: tokorhon
#SBATCH -J StR{ratio}{edge}W{width}
#SBATCH --constraint="snb|hsw"
#SBATCH -o output-%j
#SBATCH -e error_f-%j
#SBATCH -p longrun
#SBATCH -n 1
#SBATCH --nodes=1
#SBATCH --cpus-per-task={ncores}
#SBATCH -t 12-00:00:00
#SBATCH --mem-per-cpu=4000

# commands to manage the batch script
#   submission command
#     sbatch [script-file]
#   status command
#     squeue -u tokorhon
#   termination command
#     scancel [jobid]

# For more information
#   man sbatch
#   more examples in Taito guide in
#   http://research.csc.fi/taito-user-guide

module purge
module load gcc/4.8.2 intelmpi/4.1.3 fftw/3.3.4
module load python-env/2.7.6

# ASE SNAPSHOT IN APPL_TAITO
export PYTHONPATH="$PYTHONPATH:$HOME/appl_taito/python-ase-3.9.0.4137"
export PATH="$PATH:$HOME/appl_taito/pytho-ase-3.9.0.4137/tools"
export ASE_TAGS=https://svn.fysik.dtu.dk/projects/ase/tags/
# END ASE

export PYTHONPATH="$PYTHONPATH:$HOME/ShearAndSlide/src/"
export PYTHONPATH="$PYTHONPATH:$HOME/appl_taito/AidPackages/src/kc"
export PYTHONPATH="$PYTHONPATH:$HOME/appl_taito/AidPackages/src/read_write"
export LAMMPS_COMMAND="$HOME/appl_taito/lammps-9Dec14/src/lmp_taitoserial"
export LAMMPS_POTENTIALS="$HOME/appl_taito/lammps-9Dec14/potentials"

# args = width, edge, ratio, ncores
srun python $HOME/appl_taito/ShearAndSlide/src/moldy/main_runStick.py {width} {edge} {ratio} {ncores}

# This script will print some usage statistics to the
# end of file: output-%j
# Use that to improve your resource request estimate
# on later jobs.
used_slurm_resources.bash
''' 
    
    
    context = {
     "ratio":ratio, 
     "edge":edge,
     "width": w,
     "ncores" : ncores
     } 
 
    if float(ratio)%2 == 0:
        
        filef   =   path + '/w=%i/r=%i/' %(w, ratio)
        if not os.path.exists(filef):
            os.makedirs(filef)
            
        with  open(filef + 'script.sh','w') as scriptfile:
            scriptfile.write(script.format(**context))
        
        
    
    
        
    
    
    
    
    
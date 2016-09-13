#!/bin/bash
#-----------------------------------------------------------------
# SLURM job script to run serial or  applications on TACC's
# Stampede system.
#
#-----------------------------------------------------------------
    
#SBATCH -J PyRosComp          # Job name
#SBATCH -o PyRosComp.%j.out    # Specify stdout output file (%j expands to jobId)
#SBATCH -e PyRosComp.%j.err    # Specify stdout output file (%j expands to jobId)
#SBATCH -p largemem           # Queue name
#SBATCH -N 1                     # Total number of nodes requested (16 cores/node)
#SBATCH -n 32                     # Total number of tasks
#SBATCH -t 2:0:0              # Run time (hh:mm:ss) - 1.5 hours
        
# Load any necessary modules
# Loading modules in the script ensures a consistent environment.
module load gcc/4.4.6
module load python

# Launch executables
./DeployPyRosetta.tacc.stampede.py --prefix=/work/02984/cwbrown/src/stampede_pyrosetta_local_compile --rosetta-source=/work/02984/cwbrown/src/stampede_pyrosetta_local_compile/rosetta_src_2016.13.58602_bundle/main/source --jobs=32
wait


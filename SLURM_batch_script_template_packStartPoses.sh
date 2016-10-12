#!/bin/bash
#-----------------------------------------------------------------
# SLURM job script to run serial or  applications on TACC's
# Stampede system.
#
#-----------------------------------------------------------------
    
#SBATCH -J 1msw_rlx           # Job name
#SBATCH -o 1msw_rlx.%j.out    # Specify stdout output file (%j expands to jobId)
#SBATCH -e 1msw_rlx.%j.err    # Specify stdout output file (%j expands to jobId)
#SBATCH -p normal           # Queue name
#SBATCH -N 4                     # Total number of nodes requested (16 cores/node)
#SBATCH -n 64                     # Total number of tasks
#SBATCH -t 8:0:0              # Run time (hh:mm:ss) - 1.5 hours
#SBATCH -A [xxx]
#SBATCH --mail-user=[xxx@xxx.xxx]
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes
        
# Load any necessary modules
# Loading modules in the script ensures a consistent environment.
module load python

# Make sure that the PyRosetta_TACC_MPI.py module is in your $PYTHONPATH!

# Launch executables

ibrun $WORK/scripts/packStartPoses.py [PATH/TO/XXX.pdb] -p -c X --database=$PYROSETTA_DB
wait

# usage: packStartPoses.py [-h] [--n_poses N_POSES] [--database DATABASE]
#                         [--scorefxn SCOREFXN]
#                         [--restrict_to_chain RESTRICT_TO_CHAIN]
#                         [--print_scores]
#                         start_pdb
#
# Performs N parallel pre-packing runs on the input pdb
#
# positional arguments:
#  start_pdb             Path to raw .pdb file
#
# optional arguments:
#   -h, --help            show this help message and exit
#   --n_poses N_POSES, -n N_POSES
#                         Number of pre-packed poses to create; defaults to number of MPI jobs available
#   --database DATABASE, -d DATABASE
#                         location of Rosetta database to use
#   --scorefxn SCOREFXN, -s SCOREFXN
#                         Name of Rosetta score function to use
#   --restrict_to_chain RESTRICT_TO_CHAIN, -c RESTRICT_TO_CHAIN
#                         Only relax residues in these chains
#   --print_scores, -p    Print pose scores to packed_pose_scores.sc



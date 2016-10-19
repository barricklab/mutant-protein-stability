#!/bin/bash
#-----------------------------------------------------------------
# SLURM job script to run serial or  applications on TACC's
# Stampede system.
#
#-----------------------------------------------------------------
    
#SBATCH -J ddG_mut           # Job name
#SBATCH -o ddG_mut.%j.out    # Specify stdout output file (%j expands to jobId)
#SBATCH -e ddG_mut.%j.err    # Specify stdout output file (%j expands to jobId)
#SBATCH -p normal           # Queue name
#SBATCH -N 16                     # Total number of nodes requested (16 cores/node)
#SBATCH -n 256                     # Total number of nodes requested (16 cores/node)
#SBATCH -t 4:0:0              # Run time (hh:mm:ss) - 1.5 hours
#SBATCH -A [xxx]
#SBATCH --mail-user=[xxx@xxx.com]
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes
        
# Load any necessary modules
# Loading modules in the script ensures a consistent environment.
module load python

# MAKE SURE THAT THE PyRosetta_TACC_MPI.py module is in your $PYTHONPATH

# Launch executables - make sure to change these for your local env

ibrun tacc_affinity $WORK/scripts/ddG_mutants.py -r 3 -k 10 ../in_1qln_start_poses_rlx/poses $CWB_PYROSETTA_DB 88A all
wait

# usage: ddG_mutants.py [-h] [--replicates REPLICATES] [--restrict_to_chain]
#                       [--dump_ref_pdb] [--dump_mut_pdb]
#                       [--dump_pdb_base DUMP_PDB_BASE]
#                       start_pose_pdbs_dir database residues mutate_to
# 
# Predict ddG (ddREU) effects of point mutations to a protein structure. Handles
# nsAAs, uses MPI
# 
# positional arguments:
#   start_pose_pdbs_dir   Directory containing starting pose .pdb files
#   database              location of Rosetta database to use
#   residues              Residue positions to mutate
#                         <PDB_residue1PDBChain[,residue2,residue3
#                         ...]|range1Start-range1End[,range2Start-
#                         range2End,...]>
#   mutate_to             AAs to mutate to (three-letter codes)
#                         <AA1[,AA2,AA3...]|'All'|'nAAs'|'nsAAs'>
#
# optional arguments:
#   -h, --help            show this help message and exit
#   --replicates REPLICATES, -r REPLICATES
#                         Replicate runs (per starting pose)
#   --restrict_to_chain, -c
#                         Only pack residues in the chains specificed in
#                         <residues>
#   --min_restrict_radius, -d
#                          Only minimize residues within the packing 
#                          radius (Default=FALSE)')
#   --min_cst_sd MIN_CST_SD, -k MIN_CST_SD
#                          Use bb coord constraints w/ this std. dev. 
#                          during minimization. Recommend 10, default
#			   is None (no constraints)
#   --dump_ref_pdb, -f
#   --dump_mut_pdb, -m
#   --dump_pdb_base DUMP_PDB_BASE, -b DUMP_PDB_BASE
#



#!/usr/bin/env python
# test run of CWB's MPI-capable PyRosetta module on TACC lonestar

import sys
import rosetta
import numpy as np
from PyRosetta_TACC_MPI import *
from mpi4py import MPI

def _main(args):

    (start_pdb, n_start_poses, n_pack, n_min, restrict_to_chain) = (None,None,None,None,None)

    if len(args) not in (4,5):
        print "usage: <start_pdb> <n_start_poses> <n_pack> <n_min> [<restrict_to_chain>]"
        sys.exit(0)
    elif len(args) == 4:
        (start_pdb, n_start_poses, n_pack, n_min) = args
    elif len(args) == 5:
        (start_pdb, n_start_poses, n_pack, n_min,restrict_to_chain) = args

    n_start_poses = int(n_start_poses)
    n_pack = int(n_pack)
    n_min = int(n_min)
    
    rosetta.mpi_init("-database /work/02984/cwbrown/Data/Rosetta/database_w_ncaas/pyrosetta_database")

    packer_job = PackSinglePoseJob(start_pdb,

def make_packmin_trial_mover(pose,packertask,minmover,scorefn,n_packing_steps=1,n_minimize_steps=1,kT=1.0):
    """
    Build a TrialMover with n_packing_steps RotamerTrialMover moves and 
    n_minimization_steps MinMover moves executed sequentially
    """

    seqmover = rosetta.SequenceMover()
        
    if n_minimize_steps > 0:
        min_repmover = rosetta.RepeatMover(minmover,n_minimize_steps)
        seqmover.add_mover(min_repmover)
       
    if n_packing_steps > 0:
        packmover = rosetta.RotamerTrialsMover(scorefn,packertask)
        pack_repmover = rosetta.RepeatMover(packmover,n_packing_steps)
        seqmover.add_mover(pack_repmover)
        
    mc = rosetta.MonteCarlo(pose,scorefn,kT)
    packmin_trial_mover = rosetta.TrialMover(seqmover,mc)
    
    return packmin_trial_mover
    
def std_dev_threshold_fn_builder(threshold,last_n=5):
    def stop(scores):
        passed = False
        if len(scores) > last_n:
            last_n_stdev = np.std(scores[-last_n:])
            passed = last_n_stdev < threshold
            print "std dev for scores %s: %f thresh: %f pass?: %s" % (str(scores[-last_n:]),last_n_stdev,  threshold,str(passed))
        return passed
    return stop


if __name__ == "__main__":
    _main(sys.argv[1:])

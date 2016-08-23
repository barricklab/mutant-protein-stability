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

    comm = MPI.COMM_WORLD

    rank = comm.Get_rank()
    size = comm.Get_size()

    scores = []
    #if rank == 0:
    #    print "packing %d poses..." % (n_start_poses)
    #    start_pose_scores = [""] * n_start_poses
    #else:
    #	start_pose_scores = None
    #
    if rank < n_start_poses:
        print " + packing %dth pose" % rank

        packed_pose = rosetta.pose_from_pdb(start_pdb)

        mm_std_sf = rosetta.create_score_function("mm_std")

        mm = rosetta.MoveMap()
        mm.set_bb(True)
        mm.set_chi(True)

        min_mover = rosetta.MinMover()
        min_mover.movemap(mm)
        min_mover.score_function(mm_std_sf)
# Note that we're now using 'dfpmin_armijo_nonmonotone' instead of 'linmin'; seems to work a lot better
        min_mover.min_type("dfpmin_armijo_nonmonotone")

        scores.append(mm_std_sf(packed_pose))
        pack_task = rosetta.standard_packer_task(packed_pose)
        pack_task.restrict_to_repacking()
        residues = []
        if restrict_to_chain:
            for r in range(1,packed_pose.n_residue()+1):
                if packed_pose.pdb_info.chain(r) == restrict_to_chain:
                    residues.append(r)
        else:
            residues = range(1,packed_pose.n_residue()+1)

        pack_residues = rosetta.utility.vector1_bool()
        for i in range(1,packed_pose.n_residue() + 1):
            pack_residues.append(i in residues)
        packer_task.restrict_to_residues(pack_residues)
        
        pack_mover = rosetta.RotamerTrialsMover(mm_std_sf,pack_task)
        seqmover = rosetta.SequenceMover()
        for i in range(0,n_pack):
            seqmover.add_mover(pack_mover)
        for i in range(0,n_min):
            seqmover.add_mover(min_mover)

        stop_fn = std_dev_threshold_fn_builder(0.1)

        while (not stop_fn(scores)):
            print "rank: %d scores: %s" % (rank,str(scores))
            #min_mover.apply(packed_pose)
            seqmover.apply(packed_pose)
            scores.append(mm_std_sf(packed_pose))
        #singlePoseColl = StartPoseCollection(start_pdb,1,rank,n_min_steps=5,n_pack_steps=0)
        #singlePoseColl.pack_all_poses()
        
        outpose = start_pdb + "_packed_rep%d.pdb" % (rank)
        packed_pose.dump_pdb(outpose)
        print "Rank %d done, final pose written to %s" % (rank,outpose)
    else:
        pass
    
    comm.barrier()
    start_pose_scores = comm.gather(scores,root=0)

    if rank==0:
        out = open("test_out.txt","w")
        for (i,sc) in enumerate(start_pose_scores):
            print >> out, "pose%d\t%s" % (i,"\t".join([str(x) for x in sc]))

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

#!/usr/bin/env python
# test run of CWB's MPI-capable PyRosetta module on TACC lonestar

import sys
import rosetta
import numpy as np
from PyRosetta_TACC_MPI import *
from mpi4py import MPI

def _main(args):

    (start_pdb, n_start_poses, n_pack, n_min, restrict_to_chain) = (None,None,None,None,None)

    if len(args) not in (5,6):
        print "usage: <start_pdb> <n_start_poses> <n_pack> <n_min> <kT> [<restrict_to_chain>]"
        sys.exit(0)
    elif len(args) == 5:
        (start_pdb, n_start_poses, n_pack, n_min, kT) = args
    elif len(args) == 6:
        (start_pdb, n_start_poses, n_pack, n_min, kT, restrict_to_chain) = args

    n_start_poses = int(n_start_poses)
    n_pack = int(n_pack)
    n_min = int(n_min)
    out_pdb = start_pdb.split("/")[-1]
    
    comm = MPI.COMM_WORLD

    rank = comm.Get_rank()
    size = comm.Get_size()

    rosetta.mpi_init(extra_options="-ignore_unrecognized_res False -database=/work/02984/cwbrown/Data/Rosetta/database_w_ncaas/pyrosetta_r96_database")

    job_div = int(np.floor(n_start_poses / size))
    job_rem = n_start_poses - (job_div * size)
    job_div_procs = size - job_rem
    n_local_jobs = None
    if rank <= job_div_procs:
        n_local_jobs = job_div
    else:
        n_local_jobs = job_div + 1

    for i in range(0,n_local_jobs):
        pack_id = (rank * job_div) + i
        packer_job = PackSinglePoseJob(start_pdb,rank,restrict_to_chain=restrict_to_chain,kT=float(kT),MCmover=False)
        packer_job.pack_pose()
        print "POSE %d in PROCESS %d COMPLETE, WRITING TO %s" % (i,rank,start_pdb + ".pack" + str(pack_id) + ".pdb")
        write_energies = open(out_pdb + ".pack" + str(pack_id) + "_energies.txt","w")
        print >> write_energies, "\t".join(str(packer_job.scores))
        write_energies.close()
        packer_job.dump_pose(out_pdb + ".pack" + str(pack_id) + ".pdb")

if __name__ == "__main__":
    _main(sys.argv[1:])

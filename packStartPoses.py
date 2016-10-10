#!/usr/bin/env python
# test run of CWB's MPI-capable PyRosetta module on TACC lonestar

import sys,os,argparse
import rosetta
import numpy as np
from PyRosetta_TACC_MPI import *
from mpi4py import MPI

def _main(args):

    parser = argparse.ArgumentParser(description="Performs N parallel pre-packing runs on the input pdb")

    parser.add_argument('start_pdb',help='Path to raw .pdb file')
    parser.add_argument('--n_poses','-n',default=None,type=int,help='Number of pre-packed poses to create')
    parser.add_argument('--database','-d',default="",help='location of Rosetta database to use')
    parser.add_argument('--scorefxn','-s',default='talaris2014',help='Name of Rosetta score function to use')
    parser.add_argument('--restrict_to_chain','-c',default=None,help='Only relax residues in these chains')
    parser.add_argument('--print_scores','-p',action="store_true",help='Print pose scores to packed_pose_scores.sc')

    parsed_args = parser.parse_args()
    #print parsed_args

    start_pdb = parsed_args.start_pdb
    n_start_poses = parsed_args.n_poses

    if parsed_args.database != "":
        database = " -database=" + parsed_args.database
 
    scorefxn = parsed_args.scorefxn
    restrict_to_chain = parsed_args.restrict_to_chain
    print_scores = parsed_args.print_scores


#    (start_pdb, n_start_poses, n_pack, n_min, restrict_to_chain) = (None,None,None,None,None)

#    if len(args) not in (5,6):
#        print "usage: <start_pdb> <n_start_poses> <n_pack> <n_min> <kT> [<restrict_to_chain>]"
#        sys.exit(0)
#    elif len(args) == 5:
#        (start_pdb, n_start_poses, n_pack, n_min, kT) = args
#    elif len(args) == 6:
#        (start_pdb, n_start_poses, n_pack, n_min, kT, restrict_to_chain) = args    

    out_pdb = start_pdb.split("/")[-1]
    
    comm = MPI.COMM_WORLD

    rank = comm.Get_rank()
    size = comm.Get_size()
    
    n_start_poses = None

    if n_start_poses:
        n_start_poses = int(n_start_poses)
    else:
        n_start_poses = size

    rosetta_options = ["-ignore_unrecognized_res False",
                       "-ex1",
                       "-ex2",
                       "-use_input_sc",
                       "-flip_HNQ",
                       "-no_optH false",
                       "-relax:constrain_relax_to_start_coords",
                       "-relax:coord_constrain_sidechains",
                       "-relax:ramp_constraints false",
                       database]

    rosetta.mpi_init(extra_options=" ".join(rosetta_options))

    n_local_jobs = None
    
    if n_start_poses > size:        
        job_div = int(np.floor(n_start_poses / size))
        job_rem = n_start_poses - (job_div * size)
        job_div_procs = size - job_rem
        if rank <= job_div_procs:
            n_local_jobs = job_div
        else:
            n_local_jobs = job_div + 1
    else:
        job_div = 1
        if rank <= n_start_poses:
            n_local_jobs = 1
        else:
            n_local_jobs = 0

    #print "n_start_poses=%d,job_div=%d,job_rem=%d,job_div_procs=%d, n_local_jobs=%d" % (n_start_poses,job_div,job_rem,job_div_procs,n_local_jobs)

    write_energies = None
    if print_scores:
        write_energies = open(out_pdb + ".pack" + str(rank) + "_energies.txt","w+")

    for i in range(0,n_local_jobs):
        pack_id = "%d-%d" % (rank,n_local_jobs)
        outfile_name = out_pdb + ".pack" + pack_id + ".pdb"
        packer_job = FastRelaxPoseJob(start_pdb,rank,scorefn=scorefxn,restrict_to_chain=restrict_to_chain)
        packer_job.pack_pose()
        print "POSE %d in PROCESS %d COMPLETE, WRITING TO %s" % (i,rank,outfile_name)
        if print_scores:
            print >> write_energies, "\t".join([outfile_name,str(packer_job.start_score),str(packer_job.final_score)])
        packer_job.dump_pose(outfile_name)
    
    if print_scores:
        write_energies.close()

if __name__ == "__main__":
    _main(sys.argv[1:])

#!/usr/bin/env python
# test run of CWB's MPI-capable PyRosetta module on TACC lonestar

import sys,os
from rosetta import *
from PyRosetta_TACC_MPI import *
from mpi4py import MPI

def _main(args):
    
    (start_pose_pdbs_dir, residues_str, AAs_str, nreplicates,restrict_to_chain) = (None,None,None,None,None)
    if len(args) not in (4,5):
        print "usage: <start_pose_pdbs_dir> <residue1[,residue2,residue3...]|range1Start-range1End[,range2Start-range2End,...]> <AA1[,AA2,AA3...]|'All'|'nAAs'|'nsAAs'> <replicate_packing_runs> <restrict_to_chain>"
        sys.exit(0)
    elif len(args) == 4:
        (start_pose_pdbs_dir, residues_str, AAs_str, nreplicates) = args
    elif len(args) == 5:
        (start_pose_pdbs_dir, residues_str, AAs_str, nreplicates,restrict_to_chain) = args

        

    nreplicates = int(nreplicates)
    AAs = None

    if AAs_str.lower() == "all":
        AAs = nAAs + nsAAs
    elif AAs_str.lower() == "naas":
        AAs = nAAs
    elif AAs_str.lower() == "nsaas":  
        AAs = nsAAs
    else:
        AAs = AAs_str.split(",")

    residues_ls = residues_str.split(",")
    residues = []
    for r in residues_ls:
        r_sp = r.split('-')
        if len(r_sp) > 1:
            r_range = range(int(r_sp[0]),int(r_sp[1])+1)
            residues.extend(r_range)
        else:
            residues.append(int(r))
    print residues
    
    comm = MPI.COMM_WORLD

    rank = comm.Get_rank()
    size = comm.Get_size()
    
    rosetta_options="-ignore_unrecognized_res False -database /work/02984/cwbrown/Data/Rosetta/database_w_ncaas/pyrosetta_database"

    start_pose_pdbs = [start_pose_pdbs_dir + "/" + x for x in os.listdir(start_pose_pdbs_dir) if len(x) > 4 and x[-4:] == '.pdb']
    
    print "Using %d replicate starting poses..." % (len(start_pose_pdbs),)

    #start_poses.pack_all_poses()
    
    #experiments = []

    #for (i,p) in enumerate(start_poses):
    	
    cur_exp = MutagenesisExperimentRunner(start_pose_pdbs,rosetta_options,comm,residues,AAs,nreplicates,restrict_to_chain)

    cur_exp.scatter_job()

    comm.barrier()
    
    cur_exp.gather_results()

    if rank == 0:
        print "==============================="
        print "All jobs complete, writing output to test_out.txt...",
        cur_exp.dump_ddg_results("ddg_out.txt")
        print "Done!"
    
if __name__=='__main__':
   _main(sys.argv[1:])

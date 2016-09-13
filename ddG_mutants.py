#!/usr/bin/env python
# test run of CWB's MPI-capable PyRosetta module on TACC lonestar

import sys,os,re
from rosetta import *
from PyRosetta_TACC_MPI import *
from mpi4py import MPI

def _main(args):
    
    (start_pose_pdbs_dir, database, residues_str, AAs_str, nreplicates,restrict_to_chain) = (None,None,None,None,None,None)
    if len(args) not in (5,6):
        print "usage: <start_pose_pdbs_dir> <database> <PDB_residue1PDBChain[,residue2,residue3...]|range1Start-range1End[,range2Start-range2End,...]> <AA1[,AA2,AA3...]|'All'|'nAAs'|'nsAAs'> <replicate_packing_runs> <restrict_to_chain>"
        sys.exit(0)
    elif len(args) == 5:
        (start_pose_pdbs_dir,database,residues_str,AAs_str, nreplicates) = args
        restrict_to_chain = True # by default
    elif len(args) == 6:
        (start_pose_pdbs_dir,database,residues_str, AAs_str, nreplicates,restrict_to_chain) = args
        if restrict_to_chain == "True":
            restrict_to_chain = True
        else:
            restrict_to_chain = False

        
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
            r_ch1 = re.search("(\d+)(\w+)",r_sp[0])
            r_ch2 = re.search("(\d+)(\w+)",r_sp[1])
            r_range_num = range(int(r_ch1.group(1)),int(r_ch2.group(1)))
            assert r_ch1.group(2) == r_ch2.group(2), "When specifying range, chains must be the same: %s" % (r)
            r_range_chain = [r_ch1.group(2),] * len(r_range_num)
            r_range = zip(r_range_num,r_range_chain)
            residues.extend(r_range)
        else:
            r_ch = re.search("(\d+)(\w+)",r)
            residues.append((int(r_ch.group(1)),r_ch.group(2)))
    print residues
    
    comm = MPI.COMM_WORLD

    rank = comm.Get_rank()
    size = comm.Get_size()
    
    rosetta_options="-ignore_unrecognized_res False -database=%s" % (database,)

    start_pose_pdbs = [start_pose_pdbs_dir + "/" + x for x in os.listdir(start_pose_pdbs_dir) if len(x) > 4 and x[-4:] == '.pdb']
    
    print "Using %d replicate starting poses..." % (len(start_pose_pdbs),)

    #start_poses.pack_all_poses()
    
    #experiments = []

    #for (i,p) in enumerate(start_poses):
    	
    cur_exp = MutagenesisExperimentRunner(start_pose_pdbs,rosetta_options,comm,residues,AAs,nreplicates,restrict_to_chain,PDB_res=True)

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

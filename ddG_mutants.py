#!/usr/bin/env python
# test run of CWB's MPI-capable PyRosetta module on TACC lonestar

import sys,os,re,argparse
from rosetta import *
from PyRosetta_TACC_MPI import *
from mpi4py import MPI


def _main(args):

    parser = argparse.ArgumentParser(description="Predict ddG (ddREU) effects of point mutations to a protein structure. Handles nsAAs, uses MPI")

    parser.add_argument('start_pose_pdbs_dir',help='Directory containing starting pose .pdb files')
    parser.add_argument('database',help='location of Rosetta database to use')
    parser.add_argument('residues',help='Residue positions to mutate\n<PDB_residue1PDBChain[,residue2,residue3...]|range1Start-range1End[,range2Start-range2End,...]>')
    parser.add_argument('mutate_to',help="AAs to mutate to (three-letter codes) <AA1[,AA2,AA3...]|'All'|'nAAs'|'nsAAs'>")
    parser.add_argument('--replicates','-r',default=10,type=int,help='Replicate runs (per starting pose)')
    parser.add_argument('--restrict_to_chain','-c',action='store_false',help='Only pack residues in the chains specificed in <residues>')
    parser.add_argument('--dump_ref_pdb','-f',action='store_true')
    parser.add_argument('--dump_mut_pdb','-m',action='store_true')
    parser.add_argument('--dump_pdb_base','-b',default="ddg_out")

    parsed_args = parser.parse_args()
    print parsed_args

    start_pose_pdbs_dir = parsed_args.start_pose_pdbs_dir
    database = parsed_args.database
    residues_str = parsed_args.residues
    AAs_str = parsed_args.mutate_to
    nreplicates = parsed_args.replicates
    restrict_to_chain = parsed_args.restrict_to_chain
    dump_ref_pdb = parsed_args.dump_ref_pdb
    dump_mut_pdb = parsed_args.dump_mut_pdb
    dump_pdb_base = parsed_args.dump_pdb_base
    
#
#    if len(args) not in (5,6):
#        print args
#        print "usage: <start_pose_pdbs_dir> <database> <PDB_residue1PDBChain[,residue2,residue3...]|range1Start-range1End[,range2Start-range2End,...]> <AA1[#,AA2,AA3...]|'All'|'nAAs'|'nsAAs'> <replicate_packing_runs> <restrict_to_chain>"
#        sys.exit(0)
#    elif len(args) == 5:
#        (start_pose_pdbs_dir,database,residues_str,AAs_str, nreplicates) = args
#        restrict_to_chain = True # by default
#    elif len(args) == 6:
#        (start_pose_pdbs_dir,database,residues_str, AAs_str, nreplicates,restrict_to_chain) = args
#        if restrict_to_chain == "True":
#            restrict_to_chain = True
#        else:
#            restrict_to_chain = False

        
#    nreplicates = int(nreplicates)
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
    	
    cur_exp = MutagenesisExperimentRunner(start_pose_pdbs,rosetta_options,comm,residues,AAs,nreplicates,restrict_to_chain,dump_ref_pdb=dump_ref_pdb,dump_mut_pdb=dump_mut_pdb,pdb_base=dump_pdb_base,PDB_res=True)

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

#!/usr/bin/env python
# test run of CWB's MPI-capable PyRosetta module on TACC lonestar

import sys,os,re,argparse
from  pyrosetta import *
from PyRosetta_TACC_MPI import *
from mpi4py import MPI


def _main(args):

    parser = argparse.ArgumentParser(description="Predict ddG (ddREU) effects of point mutations to a protein structure. Handles nsAAs, uses MPI")

    parser.add_argument('start_pose_pdbs_dir',help='Directory containing starting pose .pdb files')
    parser.add_argument('database',help='location of Rosetta database to use')
    parser.add_argument('residues',help='Residue positions to mutate\n<PDB_residue1PDBChain[,residue2,residue3...]|range1Start-range1End[,range2Start-range2End,...]>')
    parser.add_argument('mutate_to',help="AAs to mutate to (three-letter codes) <AA1[,AA2,AA3...]|'All'|'nAAs'|'nsAAs'>")
    parser.add_argument('--replicates','-r',default=10,type=int,help='Replicate runs (per starting pose)')
    parser.add_argument('--min_restrict_radius','-d',action='store_true',help='Only minimize residues within the packing radius (Default=FALSE)')
    parser.add_argument('--min_cst_sd','-k',default=None,type=float,help='Use bb coord constraints w/ this std. dev. during minimization (recommend 0.5, default=None)')
    parser.add_argument('--restrict_to_chain','-c',action='store_false',help='Only pack residues in the chains specificed in <residues>')
    parser.add_argument('--dump_ref_pdb','-f',action='store_true')
    parser.add_argument('--dump_mut_pdb','-m',action='store_true')
    parser.add_argument('--dump_pdb_base','-b',default="ddg_out")
    parser.add_argument('--constraint_file','-t',default=None,help='Add a constraint file')
    parser.add_argument('--ddg_out_file','-g',default='ddg_out.txt',help='Textfile to output ddG results')
    parser.add_argument('--verbose','-v',default=1,help="How much output. >1 may cause IO issues on Stampede for large parallel runs")

    parsed_args = parser.parse_args()
    #print parsed_args

    start_pose_pdbs_dir = parsed_args.start_pose_pdbs_dir
    database = parsed_args.database
    residues_str = parsed_args.residues
    AAs_str = parsed_args.mutate_to
    nreplicates = parsed_args.replicates
    min_cst_sd = parsed_args.min_cst_sd
    min_restrict_radius = parsed_args.min_restrict_radius
    restrict_to_chain = parsed_args.restrict_to_chain
    dump_ref_pdb = parsed_args.dump_ref_pdb
    dump_mut_pdb = parsed_args.dump_mut_pdb
    dump_pdb_base = parsed_args.dump_pdb_base
    constraint_file = parsed_args.constraint_file
    ddg_out_file = parsed_args.ddg_out_file
    verbose = parsed_args.verbose
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
    if verbose > 1:
        print residues
    
    comm = MPI.COMM_WORLD

    rank = comm.Get_rank()
    size = comm.Get_size()

    extra_options = ["mute core.pack"]
    if rank > 0:
        extra_options = ["-mute all"]

    if constraint_file:
        extra_options.append("-constraints:cst_fa_file=%s" % constraint_file)
    if len(extra_options) > 0:
        rosetta_options="-ignore_unrecognized_res False -database=%s %s" % (database," ".join(extra_options))
    else:
        rosetta_options="-ignore_unrecognized_res False -database=%s" % (database,)

    start_pose_pdbs = [start_pose_pdbs_dir + "/" + x for x in os.listdir(start_pose_pdbs_dir) if len(x) > 4 and x[-4:] == '.pdb']
    
    if rank == 0:
        print "Using %d replicate starting poses..." % (len(start_pose_pdbs),)

    #start_poses.pack_all_poses()
    
    #experiments = []

    #for (i,p) in enumerate(start_poses):
    	
    cur_exp = MutagenesisExperimentRunner(start_pose_pdbs,rosetta_options,comm,residues,AAs,nreplicates,restrict_to_chain,dump_ref_pdb=dump_ref_pdb,dump_mut_pdb=dump_mut_pdb,pdb_base=dump_pdb_base,min_cst_sd=min_cst_sd,min_restrict_radius=min_restrict_radius,PDB_res=True,constraint_file=constraint_file,verbose=verbose)

    cur_exp.scatter_job()

    comm.barrier()
    
    cur_exp.gather_results()

    if rank == 0:
        print "==============================="
        print "All jobs complete, writing output to %s..." % (ddg_out_file,),
        cur_exp.dump_ddg_results(ddg_out_file)
        print "Done!"
    
if __name__=='__main__':
   _main(sys.argv[1:])

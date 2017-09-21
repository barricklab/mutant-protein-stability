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
    parser.add_argument('center_residue',type=int,help='residue to center the sphere of residues to mutate, i.e. the primary mutation.')
    parser.add_argument('center_residue_ref_pdb',help='pdb to use for calculating CA-CA distances to define sphere of residues to mutate')
    parser.add_argument('center_residue_chain',help='pdb chain of the primary mutated residue (will only mutate residues on this chain)')
    parser.add_argument('mutate_center_res_to',help="mutant AA of the primary mutant")
    parser.add_argument('--mutate_secondary_to','-t',default='nAAs',help='set of AAs to mutate to; \'all\',\'nAAs\',\'nsAAs\', list of three-letter codes (e.g. \"AA1,AA2,...\"')
    parser.add_argument('--ddg_out_file','-o',default='ddg_out.txt',help='name of file to output final text (tab-delimited) table of results')
    parser.add_argument('--radius','-s',default=10,type=float,help="Radius defining the sphere of residues to mutate. Default=10 Angstroms")
    parser.add_argument('--replicates','-r',default=10,type=int,help='Replicate runs (per starting pose). Default=10')
    parser.add_argument('--no_chain_restriction','-c',action='store_true',help="Don't restrict to only residues w/in the same chain as center residue")
    parser.add_argument('--dump_ref_pdb','-f',action='store_true')
    parser.add_argument('--dump_mut_pdb','-m',action='store_true')
    parser.add_argument('--dump_pdb_base','-b',default="ddg_out")
    parser.add_argument('--verbose','-v',default=1,type=int,help="How much output to show. 0=no output, 1=script progress msgs, 2=1+Rosetta status msgs")

    parsed_args = parser.parse_args()
    #print parsed_args

    start_pose_pdbs_dir = parsed_args.start_pose_pdbs_dir
    database = parsed_args.database
    center_residue = parsed_args.center_residue
    center_residue_ref_pdb = parsed_args.center_residue_ref_pdb
    center_residue_chain = parsed_args.center_residue_chain
    mutate_center_res_to = parsed_args.mutate_center_res_to
    AAs_str = parsed_args.mutate_secondary_to
    radius = parsed_args.radius
    ddg_out_file = parsed_args.ddg_out_file
    nreplicates = parsed_args.replicates
    restrict_to_chain = parsed_args.no_chain_restriction
    dump_ref_pdb = parsed_args.dump_ref_pdb
    dump_mut_pdb = parsed_args.dump_mut_pdb
    dump_pdb_base = parsed_args.dump_pdb_base
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

    #residues = []
    #for r in residues_ls:
    #    r_sp = r.split('-')
    #    if len(r_sp) > 1:
    #        r_ch1 = re.search("(\d+)(\w+)",r_sp[0])
    #        r_ch2 = re.search("(\d+)(\w+)",r_sp[1])
    #        r_range_num = range(int(r_ch1.group(1)),int(r_ch2.group(1)))
    #        assert r_ch1.group(2) == r_ch2.group(2), "When specifying range, chains must be the same: %s" % (r)
    #        r_range_chain = [r_ch1.group(2),] * len(r_range_num)
    #        r_range = zip(r_range_num,r_range_chain)
    #        residues.extend(r_range)
    #    else:
    #        r_ch = re.search("(\d+)(\w+)",r)
    #        residues.append((int(r_ch.group(1)),r_ch.group(2)))
    #print residues
    
    comm = MPI.COMM_WORLD

    rank = comm.Get_rank()
    size = comm.Get_size()
    
    rosetta_mute = ""
    if verbose < 2:
        rosetta_mute = " -mute all"
    rosetta_options="-ignore_unrecognized_res False -database=%s %s" % (database,rosetta_mute)

    start_pose_pdbs = [start_pose_pdbs_dir + "/" + x for x in os.listdir(start_pose_pdbs_dir) if len(x) > 4 and x[-4:] == '.pdb']
    
    if verbose > 0:
        print "Using %d replicate starting poses..." % (len(start_pose_pdbs),)

    #start_poses.pack_all_poses()
    
    #experiments = []

    #for (i,p) in enumerate(start_poses):
    	
    cur_exp = SecondaryMutantScanExperimentRunner(start_pose_pdbs,rosetta_options,comm,center_residue,mutate_center_res_to,center_residue_chain,radius,center_residue_ref_pdb,AAs,nreplicates,restrict_to_chain,dump_ref_pdb=dump_ref_pdb,dump_mut_pdb=dump_mut_pdb,pdb_base=dump_pdb_base,PDB_res=True,verbose=verbose)

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

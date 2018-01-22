#!/usr/bin/env python
# test run of CWB's MPI-capable PyRosetta module on TACC lonestar

import sys,os,re,argparse
import itertools as it
from pyrosetta import *
from PyRosetta_TACC_MPI import *
from mpi4py import MPI


def _main(args):

    parser = argparse.ArgumentParser(description="Predict ddG (ddREU) effects of point mutations to a protein structure. Handles nsAAs, uses MPI")

    parser.add_argument('start_pose_pdbs_dir',help='Directory containing starting pose .pdb files')
    parser.add_argument('database',help='location of Rosetta database to use')
    parser.add_argument('mutant_list',help='list of mutant sequences to test, poss. with more than one AA change/mutant. format: mut1_AA1_res,mut1_AA1_chain,mut1_AA1_toAA;mut1_AA2_res,... + mut2_AA1_res,mut2_AA1_chain,mut2_AA1_toAA; mut2_AA2_res,...')
    parser.add_argument('--replicates','-r',default=10,type=int,help='Replicate runs (per starting pose)')
    parser.add_argument('--min_restrict_radius','-d',action='store_true',help='Only minimize residues within the packing radius (Default=FALSE)')
    parser.add_argument('--min_cst_sd','-k',default=None,type=float,help='Use bb coord constraints w/ this std. dev. during minimization (recommend 0.5, default=None)')
    parser.add_argument('--restrict_to_chain','-c',action='store_false',help='Only pack residues in the chains specificed in <residues>')
    parser.add_argument('--combinations','-n',default=1,type=int,help='Test all N choose 1...n combinations of the input GENOTYPES where N is the total number input as the mutant list.')
    parser.add_argument('--dump_ref_pdb','-f',action='store_true')
    parser.add_argument('--dump_mut_pdb','-m',action='store_true')
    parser.add_argument('--dump_pdb_base','-b',default="ddg_out")

    parsed_args = parser.parse_args()
    #print parsed_args

    start_pose_pdbs_dir = parsed_args.start_pose_pdbs_dir
    database = parsed_args.database
    mutant_list = parsed_args.mutant_list
    nreplicates = parsed_args.replicates
    min_cst_sd = parsed_args.min_cst_sd
    min_restrict_radius = parsed_args.min_restrict_radius
    restrict_to_chain = parsed_args.restrict_to_chain
    dump_ref_pdb = parsed_args.dump_ref_pdb
    dump_mut_pdb = parsed_args.dump_mut_pdb
    dump_pdb_base = parsed_args.dump_pdb_base
    combinations = parsed_args.combinations

#    print >> sys.stderr, "ddG_multi_mutants.py: restrict_to_chain= %s" % (str(restrict_to_chain))
        
    mutants_seqs = [[[j.strip() for j in i.split(",")] for i in [z.strip() for z in y.split(";")]] for y in [x.strip() for x in mutant_list.split("+")]]
    # make sure residue is an int
    mutants_seqs = [[(int(x[0]),x[1],x[2]) for x in y] for y in mutants_seqs]

    mutants_seqs_comb = []

    for i in range(1,combinations+1):
        mutants_seqs_comb.extend([x for x in it.combinations(mutants_seqs,i)])

    mutants_seqs = []

    for genotype in mutants_seqs_comb:
        s = []
        for mut in genotype:
            s.extend(mut)
        mutants_seqs.append(s)

    # filter any genotypes with repeat mutations
    mutants_seqs = [x for x in mutants_seqs if (len(set([y[0] for y in x])) == len([y[0] for y in x]))]

    print >> sys.stderr, "ddG_multi_mutants.py mutant_seqs = %s" % (str(mutants_seqs),)

    comm = MPI.COMM_WORLD

    rank = comm.Get_rank()
    size = comm.Get_size()
    
    rosetta_options="-ignore_unrecognized_res False -database=%s" % (database,)

    start_pose_pdbs = [start_pose_pdbs_dir + "/" + x for x in os.listdir(start_pose_pdbs_dir) if len(x) > 4 and x[-4:] == '.pdb']
    
    print "Using %d replicate starting poses..." % (len(start_pose_pdbs),)

    #start_poses.pack_all_poses()
    
    #experiments = []

    #for (i,p) in enumerate(start_poses):
    	
    cur_exp = MutantCombinationsExperimentRunner(start_pose_pdbs,rosetta_options,comm,mutants_seqs,nreplicates,restrict_to_chain=restrict_to_chain,dump_ref_pdb=dump_ref_pdb,dump_mut_pdb=dump_mut_pdb,pdb_base=dump_pdb_base,min_cst_sd=min_cst_sd,min_restrict_radius=min_restrict_radius,PDB_res=True)

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

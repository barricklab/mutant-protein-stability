#!/usr/bin/env python
import os,sys,time,random
import numpy as np
import itertools as it
#from rosetta import *
import pyrosetta
from pyrosetta.rosetta.core import conformation
from pyrosetta.rosetta.core import chemical
from pyrosetta.mpi import mpi_init

nAAs = ['ALA','ASP','LEU','ILE','VAL','GLY','SER','THR','PRO','GLU','ASN','GLN','LYS','ARG','TRP','PHE','TYR','HIS','CYS','MET']
nsAAs = ['PRK','ACK','3IY','NOY','AZF','A69','B36','NBY','PHS','PHT','PHY','SOY','MMD','MMS','MMY','MMK'] # A69 = 3-aminotyrosine, B36 = 5-hydroxytryptophan

NSAAS_PATCH = {} # Need to get rid of refs to this in code
"""
'ACK':{'cognateAA':'LYS',
                      'type':chemical.VariantType.ACETYLATION},
               'PHS':{'cognateAA':'SER',
                      'type':chemical.VariantType.PHOSPHORYLATION},
               'PHT':{'cognateAA':'THR',
                      'type':chemical.VariantType.PHOSPHORYLATION},
               'PHY':{'cognateAA':'TYR',
                      'type':chemical.VariantType.PHOSPHORYLATION}}
"""
#=========================================================================
#
# ABSTRACT CLASSES - DO NOT call directly!
#
#=========================================================================


class AbstractExperimentRunner():

    def __init__(self,start_pose_pdbs,rosetta_init_options,comm,restrict_to_chain,max_pack_rounds,min_cst_sd,min_restrict_radius,PDB_res,dump_ref_pdb,dump_mut_pdb,pdb_base,verbose,constraint_file):

        if verbose >= 2:
            print >> sys.stderr, "AbstractExperimentRunner __init__(): restrict_to_chain= %s" % (str(restrict_to_chain),)

        self.start_pose_pdbs = start_pose_pdbs
        self.rosetta_init_options = rosetta_init_options
        self.max_pack_rounds = max_pack_rounds
        self.restrict_to_chain = restrict_to_chain
        self.pdb_res = PDB_res
        self.min_cst_sd = min_cst_sd
        self.min_restrict_radius = min_restrict_radius
        self.dump_ref_pdb = dump_ref_pdb
        self.dump_mut_pdb = dump_mut_pdb
        self.pdb_base = pdb_base
        self.verbose = verbose
        self.constraint_file = constraint_file
        self.comm = None
        self.size = 1
        self.rank = 0

    def setup_MPI(self,comm):
        if not comm:
            pyrosetta.init(extra_options=self.rosetta_init_options)
        else:
            self.comm = comm
            self.size = comm.Get_size()
            self.rank = comm.Get_rank()
        # keep a list of Job objects
            mpi_init(extra_options=self.rosetta_init_options)
            
    def scatter_job(self):
        local_jobs = self.comm.scatter(self.jobs,root=0)
        for (i,job_spec) in enumerate(local_jobs):
            if self.verbose >=1:
                print "===== Process %d, running job %s [ %d / %d ]" % (self.rank,str(job_spec),i,len(local_jobs))
            try:
                ddg_job = self.packer_job_class(*job_spec,max_rounds=self.max_pack_rounds,restrict_to_chain=self.restrict_to_chain,min_restrict_radius=self.min_restrict_radius,PDB_res=self.pdb_res,min_cst_sd=self.min_cst_sd,verbose=self.verbose,constraint_file=self.constraint_file)
            
                ddg_job.run()
                self.result.append(ddg_job.get_result())
                if self.dump_ref_pdb:
                    ddg_job.dump_ref_pdb(self.pdb_base + "_".join([ddg_job.start_pose_name] + [str(x) for x in job_spec[1:]]) + "_ref.pdb")
                if self.dump_mut_pdb:
                    start_pose_idx = "p" + str(self.start_pose_pdbs.index(job_spec[0]))
                    ddg_job.dump_mut_pdb(self.pdb_base + "_".join([ddg_job.start_pose_name] + [str(x) for x in job_spec[1:]]) + "_mut.pdb")
            except RuntimeError as e:
                print >> sys.stderr, "****WARNING: RUNTIME ERROR IN PROCESS %d JOB %s; ABORTING RUN [ERROR:%s]" % (self.rank,str(job_spec),e)
                sys.exit(0)
            except PyRosettaError as e:
                print >> sys.stderr, "****WARNING: PYROSETTA ERROR IN PROCESS %d JOB %s; SKIPPING JOB [ERROR:%s]" % (self.rank,str(job_spec),e)

    def gather_results(self):
        results = self.comm.gather(self.result,root=0)
        if self.rank == 0:
            for r in results:
                self.final_results.extend(r)
    
    def run_all_jobs_serial(self,outfile='ddg_out.txt'):
        for (i,job_spec) in enumerate(self.jobs):
            if self.verbose >= 1:
                print "===== Process %d, running job %s [ %d / %d ]" % (self.rank,str(job_spec),i,len(local_jobs))
            ddg_job = MutantddGPackerJob(*job_spec,max_rounds=self.max_pack_rounds,restrict_to_chain=self.restrict_to_chain,PDB_res=self.pdb_res)
            try:
                ddg_job.run()
            except RuntimeError:
                print >> sys.stderr, "****WARNING: RUNTIME ERROR IN PROCESS %d JOB: %s %s; SKIPPING JOB" % (self.rank,str(job_spec),ddg_job.repr_str)
            self.final_results.append(ddg_job.get_result())
        self.dump_ddg_results(outfile)
    
    def dump_ddg_results(self,outfile_name,header=True):
        """
        Specification Required!
        """
        pass

class AbstractPackerJob():
    
    def __init__(self,convergence_fn=None,\
                 conv_threshold=0.1,repack_radius=10,scorefn="mm_std",\
                 mintype="dfpmin_armijo_nonmonotone",max_rounds=100,restrict_to_chain=None,verbose=1,constraint_file=None):
        self.convergence_fn = convergence_fn
        self.cnv_fn = None
        self.conv_threshold = conv_threshold
        if not convergence_fn or convergence_fn == 'stdev': # defaults to stdev for now, may add more later
            self.cnv_fn = self.std_dev_threshold_fn_builder(conv_threshold)
        self.max_rounds = max_rounds    
        self.repack_radius = repack_radius
        self.mintype = mintype
        self.restrict_to_chain = restrict_to_chain
        self.verbose = verbose
        # a PyRosetta ScoreFunction object, defaults to mm_std    
        self.scorefn = pyrosetta.create_score_function(scorefn)        
        self.constraint_file = constraint_file

    def mutate_aa(self,pose,residue,aa_name,orient_bb=True,repack_sidechain=True,clone_pose=True):
        """
        Swap w/t AA at residue number 'residue' in 'pose' with 'ncaa_name' (3-letter code)
    
        Return a new Pose object
        
        Note that this assumes the ncaa .params and .rotlib files have been permanently added to the database
        """

        mut_pose = pose
        if clone_pose:
            mut_pose = pyrosetta.Pose()
            mut_pose.assign(pose)

        res = mut_pose.residue(residue)
        ref_res_name = res.name()
        # check for disulfides and correct if needed
        if (res.name() == 'CYS:disulfide') or (res.name() == 'CYD'):
            disulfide_partner = None
            try:
                disulfide_partner = res.residue_connection_partner(
                    res.n_residue_connections())
            except AttributeError:
                disulfide_partner = res.residue_connection_partner(
                    res.n_current_residue_connections())
            temp_pose = pyrosetta.Pose()
            temp_pose.assign(mut_pose)
            # (Packing causes seg fault if current CYS residue is not
            # also converted before mutating.)
            conformation.change_cys_state(residue, 'CYS',\
                             temp_pose.conformation())
            conformation.change_cys_state(disulfide_partner, 'CYS',\
                             temp_pose.conformation())
            mut_pose = temp_pose

    
            # Get a Residue object for the desired ncAA
        mut_res_type = None
        rts = None
        chm = chemical.ChemicalManager.get_instance()
        
        try:
            rts = chm.residue_type_set("fa_standard").get()
        except AttributeError:
            # older versions of PyRosetta (e.g. the TACC stampede install of 3.4) 
            # have a slightly different call to get the residue_type_set object:
            rts = chm.residue_type_set("fa_standard")
    
            
        if aa_name in NSAAS_PATCH.keys():
            # PTMs and nsAAs using the patch system need to be treated differently
            # for TACC PyRosett install; don't know if this will work with newer versions
            # where VariantType enumerator class can be called direct
            cognate_res_type = rts.name_map( NSAAS_PATCH[aa_name]['cognateAA'] ) 
            mut_res_type = rts.get_residue_type_with_variant_added(cognate_res_type,NSAAS_PATCH[aa_name]['type'])            
        else:
            # replace the target residue with the ncAA
            #if residue == 1 or residue == mut_pose.total_residue():
            #    aa_name = aa_name + "_p"
            mut_res_type = rts.name_map(aa_name)

        if residue == 1:
            try: 
                mut_res_type = rts.get_residue_type_with_variant_added(mut_res_type,chemical.VariantType.LOWER_TERMINUS_VARIANT)                
            except AttributeError:
                mut_res_type = rts.get_residue_type_with_variant_added(mut_res_type,"LOWER_TERMINUS")
        elif residue == mut_pose.total_residue():
            try:
                mut_res_type = rts.get_residue_type_with_variant_added(mut_res_type,chemical.VariantType.UPPER_TERMINUS_VARIANT)                
            except AttributeError:
                mut_res_type = rts.get_residue_type_with_variant_added(mut_res_type,"UPPER_TERMINUS")
        mut_res = conformation.ResidueFactory.create_residue( mut_res_type )
        mut_pose.replace_residue(residue,mut_res,orient_backbone=orient_bb)

        # Highly recommended to repack the sidechain after calling Pose.replace_residue()...otherwise
        # later packing steps get caught in weird local minima :'(
        #tmp_verbose = self.verbose
        #self.verbose=2
        if repack_sidechain:
            # annoying but apparently has to be done this way?...
            repack_task = pyrosetta.rosetta.core.pack.task.TaskFactory.create_packer_task(mut_pose)
            repack_task.restrict_to_repacking()
            repack_res_list = pyrosetta.rosetta.utility.vector1_bool()
            for i in range(1,mut_pose.total_residue()+1):
                repack_res_list.append(i == residue)
            repack_task.restrict_to_residues(repack_res_list)
            sidechain_PackRotamersMover = pyrosetta.rosetta.protocols.simple_moves.PackRotamersMover(self.scorefn,repack_task)
            
            initial_mut_score = self.scorefn(mut_pose)
            if self.verbose > 1:
                print "Packer task:"
                print repack_task
                print "Repacking mutated sidechain %s %s -> %s..." % (mut_pose.pdb_info().pose2pdb(residue),ref_res_name,aa_name),
            sidechain_PackRotamersMover.apply(mut_pose)
            repacked_mut_score = self.scorefn(mut_pose)
            if self.verbose > 1:
                print "initial E = %f repacked E = %f" % (initial_mut_score,repacked_mut_score)
        #self.verbose=tmp_verbose

        return mut_pose

        
    def packer_task_repack_in_radius(self,pose,residue,radius):
        """
        Build a packer task that repacks all residues with CAlphas within distance radius
    
        Might want to remake to take centroid side-chains?  Looks like SwitchResidueTypeMover has trouble w/ nsAAs...
        """
        pack_residues = self.residue_CAs_in_radius(pose,residue,radius)
        return self.make_packer_task_with_residues(pose,pack_residues)

    def make_packer_task_with_residues(self,pose,residues=None):
        """
        Builds a packer task with the specified residues activated for repacking.
    
        Did this to avoid PackerTask.temporarily_* methods which apparently we're not supposed to use
        """
        packer_task = pyrosetta.rosetta.core.pack.task.TaskFactory.create_packer_task(pose)
        packer_task.restrict_to_repacking()
        if self.verbose > 1:
            print residues
        if residues != None:
            # Vector1 doesn't translate booleans correctly in 
            # TACC version of PyRosetta; need to build utility.vector1_bool() directly
            #pack_residues = Vector1([x in residues for x in range(1,pose.total_residue())])
            pack_residues = pyrosetta.rosetta.utility.vector1_bool()
            for i in range(1,pose.total_residue() + 1):
                if self.restrict_to_chain and not (pose.pdb_info().chain(i) in self.restrict_to_chain):
                    continue
                else:
                    #print (i,i in residues)
                    pack_residues.append(i in residues)
            packer_task.restrict_to_residues(pack_residues)
        
        return packer_task

    def residue_CAs_in_radius(self,pose,centerAA,radius):
        """
        Get a list of residues with C-alpha atoms within 'radius' distance from centerAA'a C-alpha
        """
        centerAA_CA = pose.residue(centerAA).xyz('CA')
        repack_residues = []
        if self.verbose > 2:
            print >> sys.stderr, "residue_CA_in_radius: centerAA=%d, restrict_to_chain=%s" % (centerAA,str(self.restrict_to_chain))
        for i in range(1,pose.total_residue()):
            # This is a little fudgy since restrict_to_chain could be either a list or a string depending on the calling child class
            # 
            # Note that this may break if chains are not single characters
            if self.restrict_to_chain and not (pose.pdb_info().chain(i) in self.restrict_to_chain):
                continue
            test_CA = pose.residue(i).xyz('CA')
            displacement = centerAA_CA - test_CA
            distance = displacement.norm()
            if distance <= radius:
                repack_residues.append(i)
    
        return repack_residues
    
    def make_CompoundMover(self,movers,repeats_per_round):
        
        class CompoundMover:

            def __init__(self,movers,repeats):
                self.movers = movers
                self.repeats = repeats

            def apply(self,pose):
                for (mv,reps) in zip(self.movers,self.repeats):
                    for i in range(0,reps):
                        mv.apply(pose)
        
        return CompoundMover(movers,repeats_per_round)

    def make_packmin_mover(self,pose,packertask,minmover,n_packing_steps,n_minimize_steps,kT=1.0):
        """
        Build a TrialMover with n_packing_steps RotamerTrialMover moves and 
            n_minimization_steps MinMover moves executed sequentially
        """
        
        seqmover = pyrosetta.rosetta.protocols.moves.SequenceMover()
        packmover = pyrosetta.rosetta.protocols.simple_moves.PackRotamersMover(self.scorefn,packertask)

        for i in range(0,n_packing_steps):
            #pack_repmover = rosetta.RepeatMover(packmover,n_packing_steps)
            seqmover.add_mover(packmover)
        
        for i in range(0,n_minimize_steps):
            #min_repmover = rosetta.RepeatMover(minmover,n_minimize_steps)
            seqmover.add_mover(minmover)
         

        #print >> sys.stderr, seqmover
        
        mc = pyrosetta.rosetta.protocols.moves.MonteCarlo(pose,self.scorefn,kT)
        packmin_trial_mover = pyrosetta.rosetta.protocols.moves.TrialMover(seqmover,mc)
    
        return seqmover
    
    def std_dev_threshold_fn_builder(self,threshold,last_n=5):
        def stop(scores):
            passed = False
            if len(scores) > last_n:
                last_n_stdev = np.std(scores[-last_n:])
                passed = last_n_stdev < threshold
                if self.verbose >= 2:
                    print "std dev for scores %s: %f thresh: %f pass?: %s" % (str(scores[-last_n:]),last_n_stdev,  threshold,str(passed))
            return passed
        return stop

    def build_movemap(self,pose,restrict_radius_center=None,restrict_radius=None,restrict_residues=None):
        if not restrict_radius:
            restrict_radius = self.repack_radius
        mm = pyrosetta.rosetta.core.kinematics.MoveMap()
        mm.set_bb(False)
        mm.set_chi(False)
        mm.set_jump(False)
        residues = []
        if self.restrict_to_chain:
            for i in range(1,pose.total_residue()):
                if pose.pdb_info().chain(i) in self.restrict_to_chain:
                    residues.append(i)
        else:
            residues = range(1,pose.total_residue() + 1)
        
        if restrict_radius_center:
            possible_residues = self.residue_CAs_in_radius(pose,restrict_radius_center,restrict_radius)
            restricted_residues = [r for r in possible_residues if r in resiudes]
            residues = restricted_residues
        elif restrict_residues:
            restricted_residues = [r for r in restrict_residues if r in residues]
            residues = restricted_residues

        for r in residues:
            mm.set_bb(r,True)
            mm.set_chi(r,True)
        
        return mm
        
    def make_minmover(self,mintype,movemap,pose):
        # Set up a Minimization mover using the mm_std score function
        min_mover = pyrosetta.rosetta.protocols.simple_moves.MinMover()
        min_mover.movemap(movemap)
        min_scorefn = self.scorefn.clone()
        if self.min_cst_sd:
            min_scorefn.set_weight(pyrosetta.rosetta.core.scoring.coordinate_constraint,1.0) # arbitrary weight for now
        min_mover.score_function(min_scorefn)
        min_mover.min_type(mintype)
        
        return min_mover

    def add_constraints_to_scorefxn(self,constraint_types=None,weights=None,default_weight=0.1):
        if not constraint_types:
            constraint_types = [pyrosetta.rosetta.core.scoring.constraints.atom_pair_constraint, \
                                pyrosetta.rosetta.core.scoring.constraints.angle_constraint, \
                                pyrosetta.rosetta.core.scoring.constraints.dihedral_constraint, \
                                pyrosetta.rosetta.core.scoring.constraints.coordinate_constraint, \
                                pyrosetta.rosetta.core.scoring.constraints.constant_constraint]
        if not weights:
            weights = [default_weight,] * len(constraint_types)

        for (constraint,weight) in zip(constraint_types,weights):
            if self.scorefn.get_weight(constraint)==0:
                self.scorefn.set_weight(constraint, weight)
    
    def add_constraints_to_pose(self,pose,constraint_file):
        # Still not sure how to implement this in PyRosetta4 ?
        #setup = rosetta.ConstraintSetMover()
        #setup.constraint_file(constraint_file)
        #setup.apply(pose)
        pyrosetta.rosetta.core.scoreing.constraints.add_constraints_from_cmdline_to_pose(pose)


        
class PyRosettaError(Exception):
    pass

#=========================================================================
#
# ExperimentRunner Classes
#
#=========================================================================


class MutagenesisExperimentRunner(AbstractExperimentRunner):

    def __init__(self,start_pose_pdbs,rosetta_init_options,comm,residue_list=[],AA_list=[],nreps=50,restrict_to_chain=False,max_pack_rounds=25,min_cst_sd=None,min_restrict_radius=False,PDB_res=False,dump_ref_pdb=False,dump_mut_pdb=False,pdb_base="",verbose=1,constraint_file=None):

        self.packer_job_class = MutantddGPackerJob

#def __init__(self,start_pose_pdbs,rosetta_init_options,comm,restrict_to_chain=False,max_pack_rounds=25,min_cst_sd=None,min_restrict_radius=False,PDB_res=False,dump_ref_pdb=False,dump_mut_pdb=False,pdb_base="",):

        AbstractExperimentRunner.__init__(self,start_pose_pdbs,\
                                              rosetta_init_options,\
                                              comm,\
                                              restrict_to_chain=restrict_to_chain,\
                                              max_pack_rounds=max_pack_rounds,\
                                              min_cst_sd=min_cst_sd,\
                                              min_restrict_radius=min_restrict_radius,\
                                              PDB_res=PDB_res,\
                                              dump_ref_pdb=dump_ref_pdb,\
                                              dump_mut_pdb=dump_mut_pdb,\
                                              pdb_base=pdb_base,\
                                              verbose=verbose,\
                                              constraint_file=constraint_file)

        self.setup_MPI(comm)
        self.setup_jobs(residue_list,AA_list,nreps)
        
    def setup_jobs(self,residue_list,AA_list,nreps):

        if self.rank == 0:
            self.jobs = self.build_job_list(residue_list,AA_list,nreps)
        else:
            self.jobs = None
        
        if self.rank ==0:
            self.final_results = []
        else:
            self.final_results = None
        
        self.result = []

    def build_job_list(self,residues,AAs,replicates,shuffle_jobs=True):
        # build job objects for the entire run
        jobs = []
        for start_pdb in self.start_pose_pdbs:
            for (res,chain) in residues:
                for AA in AAs:
                    for rep in np.arange(replicates):
                        jobs.append((start_pdb,res,chain,AA,rep))
        if self.verbose >=1:
            print "+++++++ Sending %d jobs to %d MPI child processes" % (len(jobs),self.size)
        if shuffle_jobs:
            random.shuffle(jobs) 
        job_div = int(np.floor(len(jobs) / self.size))
        job_rem = len(jobs) - (job_div * self.size)
        job_div_procs = self.size - job_rem

        split_jobs = [jobs[i:i+job_div] for i in range(0,job_div_procs*job_div,job_div)]
        split_jobs.extend([jobs[i:i+(job_div+1)] for i in range(job_div_procs*job_div,len(jobs),job_div+1)])
        #if self.verbose >=1:
        #    print "*********************  n_jobs: %d job_div: %d job_rem: %d" % (len(jobs),job_div,job_rem)
        #if job_rem > 0:
        #    split_jobs.append(jobs[-job_rem:])
        #    print >> sys.stderr,  "*********** split_jobs: %d" % (len(split_jobs),)
        
        return split_jobs    
    
    def dump_ddg_results(self,outfile_name,header=True):
        """
        Dump results to a .tsv file
        """
        result_fields =['start_pose_pdb',\
                        'start_pose_name',\
                        'Pose_residue',\
                        'PDB_residue',\
                        'refAA',\
                        'mutAA',\
                        'replicate',\
                        'ref_start_score',\
                        'mut_start_score',\
                        'ref_final_score',\
                        'mut_final_score',\
                        'ddG',\
                        'packing_rounds',\
                        'packing_steps',\
                        'minimize_steps',\
                        'run_time']

        outfile = open(outfile_name,"w")
        if header:
            print >> outfile, "\t".join(["jobID","MPI_rank"] + result_fields)
        for (i,job_result) in enumerate(self.final_results):
            print >> outfile, "\t".join([str(i),] + [str(self.rank),] + [str(job_result[f]) for f in result_fields])
        outfile.close()

"""
class SecondaryMutantsExperimentRunner(MutagenesisExperimentRunner):
    
    def __init__(self,start_pose_pdbs,rosetta_init_options,comm,ref_pdb,center_residue,center_residue_chain,radius,nreps=50,restrict_to_chain=False,max_pack_rounds=25,PDB_res=False,dump_ref_pdb=False,dump_mut_pdb=False,pdb_base=""):
        
"""

class MutantCombinationsExperimentRunner(AbstractExperimentRunner):
    # Set up experiment where a list of specific variants (possibly w/ more than one mutation each) are made
    #
    # Pass in mutant_list, where each list item is a mutant represented by a list of tuples of the form (residue,AA)
    #   e.g. to test a mutant with Arg234 -> Ser, Tyr397 -> 3IY we'd do [[(234,'ARG'),(397,'3IY')],...]
    def __init__(self,start_pose_pdbs,rosetta_init_options,comm,mutant_list=[],nreps=50,restrict_to_chain=False,max_pack_rounds=25,min_cst_sd=None,min_restrict_radius=False,PDB_res=False,dump_ref_pdb=False,dump_mut_pdb=False,pdb_base="",verbose=1,constraint_file=None):

        self.packer_job_class = MultiMutantPackerJob
#        AbstractExperimentRunner.__init__(self,start_pose_pdbs,rosetta_init_options,restrict_to_chain,max_pack_rounds,min_cst_sd,min_restrict_radius,PDB_res,dump_ref_pdb,dump_mut_pdb,pdb_base)
        AbstractExperimentRunner.__init__(self,start_pose_pdbs,\
                                              rosetta_init_options,\
                                              comm,\
                                              restrict_to_chain=restrict_to_chain,\
                                              max_pack_rounds=max_pack_rounds,\
                                              min_cst_sd=min_cst_sd,\
                                              min_restrict_radius=min_restrict_radius,\
                                              PDB_res=PDB_res,\
                                              dump_ref_pdb=dump_ref_pdb,\
                                              dump_mut_pdb=dump_mut_pdb,\
                                              pdb_base=pdb_base,\
                                              verbose=verbose,\
                                              constraint_file=constraint_file)
        self.setup_MPI(comm)
        self.setup_jobs(mutant_list,nreps)

    def setup_jobs(self,mutant_list,nreps):
        if self.rank == 0:
            self.jobs = self.build_job_list(mutant_list,nreps)
        else:
            self.jobs = None
        
        if self.rank ==0:
            self.final_results = []
        else:
            self.final_results = None
        
        self.result = []
        
    def build_job_list(self,mutant_list,replicates,shuffle_jobs=True):
        # build job objects for the entire run
        jobs = []
        for start_pdb in self.start_pose_pdbs:
            for mutant in mutant_list:
                for rep in np.arange(replicates):
                    jobs.append((start_pdb,mutant,rep))
        if self.verbose >= 1:
            print "+++++++ Sending %d jobs to %d MPI child processes" % (len(jobs),self.size)
        if shuffle_jobs:
            random.shuffle(jobs) 
        job_div = int(np.floor(len(jobs) / self.size))
        job_rem = len(jobs) - (job_div * self.size)
        job_div_procs = self.size - job_rem

        split_jobs = [jobs[i:i+job_div] for i in range(0,job_div_procs*job_div,job_div)]
        split_jobs.extend([jobs[i:i+(job_div+1)] for i in range(job_div_procs*job_div,len(jobs),job_div+1)])
        if self.verbose >=1:
            print >> sys.stderr, "*********************  n_jobs: %d job_div: %d job_rem: %d" % (len(jobs),job_div,job_rem)
        #if job_rem > 0:
        #    split_jobs.append(jobs[-job_rem:])
        
            print >> sys.stderr,  "*********** split_jobs: %d" % (len(split_jobs),)
            print >> sys.stderr, [len(x) for x in split_jobs]
        return split_jobs

    def dump_ddg_results(self,outfile_name,header=True):
        """
        Dump results to a .tsv file
        """
        result_fields =['start_pose_pdb',\
                        'start_pose_name',\
                        'mutation_list',\
                        'replicate',\
                        'ref_start_score',\
                        'mut_start_score',\
                        'ref_final_score',\
                        'mut_final_score',\
                        'ddG',\
                        'packing_rounds',\
                        'packing_steps',\
                        'minimize_steps',\
                        'run_time']

        outfile = open(outfile_name,"w")
        if header:
            print >> outfile, "\t".join(["jobID","MPI_rank"] + result_fields)
        for (i,job_result) in enumerate(self.final_results):
            print >> outfile, "\t".join([str(i),] + [str(self.rank),] + [str(job_result[f]) for f in result_fields])
        outfile.close()

class SecondaryMutantScanExperimentRunner(MutantCombinationsExperimentRunner):

    def __init__(self,start_pose_pdbs,rosetta_init_options,comm,center_residue,center_residue_AA,center_residue_chain,mutate_radius,center_res_ref_pose,mutate_secondary_to_AAs,nreps=50,restrict_to_chain=False,max_pack_rounds=25,min_cst_sd=None,min_restrict_radius=False,PDB_res=False,dump_ref_pdb=False,dump_mut_pdb=False,pdb_base="",verbose=1):

        # init so we can compute the residues to mutate; will re-initialize via AbstractExperimentRunner after all class-specific __init__ business is done
        mpi_init(extra_options=rosetta_init_options)

        self.center_res = None
        self.center_residue_AA = center_residue_AA
        self.mutate_radius = mutate_radius
        self.center_residue_chain = center_residue_chain
        # Set up to use center residue coordinates from a single pose for consistency for now; could change to pass this off to individual poses
        self.center_res_ref_pose = pyrosetta.pose_from_file(center_res_ref_pose)
        self.mutate_secondary_to_AAs = mutate_secondary_to_AAs
        self.verbose = verbose

        if PDB_res:
            self.center_res = self.center_res_ref_pose.pdb_info().pdb2pose(self.center_residue_chain,center_residue)
        else:
            self.center_res = center_res

        radius_list = self.residue_CAs_in_radius(self.center_res_ref_pose,self.center_residue_chain,self.center_res,self.mutate_radius)

        mutant_list = [[(self.center_res,self.center_residue_chain,self.center_residue_AA),(x,self.center_residue_chain,y)] for (x,y) in it.product(radius_list,self.mutate_secondary_to_AAs) if x != self.center_res] 
        if self.verbose >= 2:
            print >> sys.stderr, "SecondaryMutantScanExperimentRunner: radius_list=%s mutant_list=%s center_res_ref_pose=%s center_residue_chain=%s center_res=%s mutate_radius=%s" % (str(radius_list),str(mutant_list),center_res_ref_pose,str(self.center_residue_chain),str(self.center_res),str(self.mutate_radius))

        MutantCombinationsExperimentRunner.__init__(self,\
                                                    start_pose_pdbs,\
                                                    rosetta_init_options,\
                                                    comm,\
                                                    mutant_list=mutant_list,\
                                                    nreps=nreps,\
                                                    restrict_to_chain=restrict_to_chain,\
                                                    max_pack_rounds=max_pack_rounds,\
                                                    min_cst_sd=min_cst_sd,\
                                                    min_restrict_radius=min_restrict_radius,\
                                                    PDB_res=False,\
                                                    dump_ref_pdb=dump_ref_pdb,\
                                                    dump_mut_pdb=dump_mut_pdb,\
                                                    pdb_base=pdb_base,\
                                                    verbose=verbose)

    def residue_CAs_in_radius(self,center_pose,center_chain,centerAA,radius):
        """
        Get a list of residues with C-alpha atoms within 'radius' distance from centerAA'a C-alpha
        """
        centerAA_CA = center_pose.residue(centerAA).xyz('CA')
        repack_residues = []
        if self.verbose >= 2:
            print "SecondaryMutantScanExperimentRunner.residue_CA_in_radius: centerAA=%d, restrict_to_chain=%s total_res=%s" % (centerAA,str(center_chain),str(center_pose.total_residue()))
        for i in range(1,center_pose.total_residue()):
            # This is a little fudgy since restrict_to_chain could be either a list or a string depending on the calling child class
            # 
            # Note that this may break if chains are not single characters
            if not (center_pose.pdb_info().chain(i) == center_chain):
                continue
            test_CA = center_pose.residue(i).xyz('CA')
            displacement = centerAA_CA - test_CA
            distance = displacement.norm()
            if distance <= radius:
                repack_residues.append(i)
    
        return repack_residues

#        self.packer_job_class = MultiMutantPackerJob
#        AbstractExperimentRunner.__init__(self,start_pose_pdbs,rosetta_init_options,restrict_to_chain,max_pack_rounds,min_cst_sd,min_restrict_radius,PDB_res,dump_ref_pdb,dump_mut_pdb,pdb_base)
#        AbstractExperimentRunner.__init__(self,start_pose_pdbs,\
#                                              rosetta_init_options,\
#                                              comm,\
#                                              restrict_to_chain=restrict_to_chain,\
#                                              max_pack_rounds=max_pack_rounds,\
#                                              min_cst_sd=min_cst_sd,\
#                                              min_restrict_radius=min_restrict_radius,\
#                                              PDB_res=PDB_res,\
#                                              dump_ref_pdb=dump_ref_pdb,\
#                                              dump_mut_pdb=dump_mut_pdb,\
#                                              pdb_base=pdb_base)
#        self.setup_MPI(comm)
#        self.setup_jobs(mutant_list,nreps)

#=========================================================================
#
# PackerJob Classes
#
#=========================================================================    


class MutantddGPackerJob(AbstractPackerJob):

    def __init__(self,start_pose_pdb,residue,chain,AA,replicate,convergence_fn=None,\
                 conv_threshold=0.1,repack_radius=10,scorefn="mm_std",\
                 mintype="dfpmin_armijo_nonmonotone",n_pack_steps=3,n_min_steps=1,max_rounds=100,min_cst_sd=None,min_restrict_radius=False,restrict_to_chain=False,PDB_res=False,verbose=1,constraint_file=None):

        self.chain = chain
        self.restrict_packing_to_chain = None
        self.min_cst_sd = min_cst_sd
        #self.constraint_file = constraint_file
        self.min_restrict_radius = min_restrict_radius
        if restrict_to_chain:
            self.restrict_packing_to_chain = self.chain

        AbstractPackerJob.__init__(self,convergence_fn,conv_threshold,repack_radius,scorefn,mintype,max_rounds,self.restrict_packing_to_chain,verbose=verbose,constraint_file=constraint_file)
        
        # Mutant definition - residue to mutate, AA to mutate to, and replicate number
        self.start_pose_pdb = start_pose_pdb
        self.start_pose_name = os.path.split(start_pose_pdb)[1]
        start_pose_in = None
        # GAAAAAAAH Why do they not make these backwards-compatible >8-O!!!
        try:
            start_pose_in = pyrosetta.pose_from_pdb(self.start_pose_pdb)
        except AttributeError:
            start_pose_in = pyrosetta.pose_from_file(self.start_pose_pdb)
        self.start_pose = pyrosetta.Pose()
        self.start_pose.assign(start_pose_in)
        self.residue = None
        if PDB_res:
            self.residue = self.start_pose.pdb_info().pdb2pose(self.chain,residue)
        else:
            self.residue = residue
        self.AA = AA 
        self.replicate = replicate

        if self.residue <= 0:
            raise PyRosettaError("MISSING RESIDUE %s %d %d %s %s" % (self.start_pose_pdb,self.residue,residue,self.chain,self.AA))

        self.repr_str = "%s %s %s %d" % (self.start_pose.residue(self.residue).name(),self.residue,self.AA,self.replicate)

        self.n_pack_steps = n_pack_steps
        self.n_min_steps = n_min_steps

        # set up starting poses for reference and mutant
        self.ref_pose = pyrosetta.Pose()
        self.mut_pose = pyrosetta.Pose()
        self.ref_pose.assign(self.start_pose)
        self.mut_pose.assign(self.start_pose)
        
        # set up cst file constraints if specified; not sure about how this interacts w/ coord constraints?
        if self.constraint_file:
            if self.min_cst_sd:
                self.add_constraints_to_scorefxn(constraint_types=[pyrosetta.rosetta.core.scoring.constraints.atom_pair_constraint, \
                                                                   pyrosetta.rosetta.core.scoring.constraints.angle_constraint, \
                                                                   pyrosetta.rosetta.core.scoring.constraints.dihedral_constraint])
            else:
                self.add_constraints_to_scorefxn()
            self.add_constraints_to_pose(self.ref_pose,self.constraint_file)
            self.add_constraints_to_pose(self.mut_pose,self.constraint_file)

        # Add coordinate constraints if specified
        # Note for future reference that this needs to be done *AFTER* setting file constraints
        # (loading file constraints apparently deletes any existing constraints?)

        if self.min_cst_sd:
            pyrosetta.rosetta.core.scoring.constraints.add_coordinate_constraints(self.ref_pose,self.min_cst_sd,False)
            pyrosetta.rosetta.core.scoring.constraints.add_coordinate_constraints(self.mut_pose,self.min_cst_sd,False)

        
        # This is a hack to deal with chirality problems for some nsAAs when mutating from Glycine
        #
        #    Pose.replace_residue() seems to place sidechain correctly 
        #    for ALA, and from any chiral AA to the nsAAs so mutate GLY to ALA first
        #if  self.start_pose.residue(self.residue).name3() == "GLY":
        #    self.mutate_aa(self.ref_pose,self.residue,"ALA",clone_pose=False)
        #    self.mutate_aa(self.mut_pose,self.residue,"ALA",clone_pose=False)

        self.mutate_aa(self.ref_pose,self.residue,self.start_pose.residue(self.residue).name3(),clone_pose=False)
        self.mutate_aa(self.mut_pose,self.residue,self.AA,clone_pose=False)


        if self.verbose > 1:
            print "Score function and Weights for %s:" % self.repr_str
            print "Ref Pose:"
            print self.scorefn.show(self.ref_pose)
            print "Mut Pose:"
            print self.scorefn.show(self.mut_pose)

        
        # Store initial pose scores
        self.raw_start_score = self.scorefn(self.start_pose)
        self.ref_start_score = self.scorefn(self.ref_pose)
        self.mut_start_score = self.scorefn(self.mut_pose)
        self.ref_final_score = None
        self.mut_final_score = None
        
        # Set up PackerTasks and Movers 
        self.ref_packertask = self.packer_task_repack_in_radius(self.ref_pose,self.residue,self.repack_radius)
        self.mut_packertask = self.packer_task_repack_in_radius(self.mut_pose,self.residue,self.repack_radius)        

        
        self.ref_min_movemap = None
        self.mut_min_movemap = None

        if self.min_restrict_radius:
            self.ref_min_movemap = self.build_movemap(self.ref_pose,self.residue)
            self.mut_min_movemap = self.build_movemap(self.mut_pose,self.residue)
        else:
            self.ref_min_movemap = self.build_movemap(self.ref_pose)
            self.mut_min_movemap = self.build_movemap(self.mut_pose)

        self.ref_min_mover = self.make_minmover(mintype,self.ref_min_movemap,self.ref_pose)
        self.mut_min_mover = self.make_minmover(mintype,self.mut_min_movemap,self.mut_pose)
            
        self.ref_pack_mover = pyrosetta.rosetta.protocols.simple_moves.PackRotamersMover(self.scorefn,self.ref_packertask)
        self.mut_pack_mover = pyrosetta.rosetta.protocols.simple_moves.PackRotamersMover(self.scorefn,self.mut_packertask)

        self.ref_mover = self.make_CompoundMover([self.ref_pack_mover,self.ref_min_mover],[self.n_pack_steps,self.n_min_steps])
        self.mut_mover = self.make_CompoundMover([self.mut_pack_mover,self.mut_min_mover],[self.n_pack_steps,self.n_min_steps])
        
        self.ref_trialmover_scores = []
        self.mut_trialmover_scores = []
        
        # Keep track of where we are in the run
        self.started = False
        self.finished = False
        self.rnd = 0
        self.ddG = None
        self.start_time = None
        self.end_time = None      
          
    def run(self):
        # Grab the score function from w/t TrialMover object
        # Might want to check that same as mut score fxn in the future
        self.rnd = 1
        self.started = True
        self.start_time = time.time()
    
        while ((not self.cnv_fn(self.ref_trialmover_scores)) or (not self.cnv_fn(self.mut_trialmover_scores))) and (self.rnd <= self.max_rounds):
            self.mut_mover.apply(self.mut_pose)
            self.ref_mover.apply(self.ref_pose)

            self.mut_trialmover_scores.append(self.scorefn(self.mut_pose))
            self.ref_trialmover_scores.append(self.scorefn(self.ref_pose))

            if self.verbose >= 2:
                print "%s ROUND %d SCORES: w/t: %f Mut: %f" % (self.repr_str,self.rnd,self.ref_trialmover_scores[-1],self.mut_trialmover_scores[-1])
            self.rnd += 1
        
        self.end_time = time.time()
        self.run_time = self.end_time - self.start_time
        self.ref_final_score = self.ref_trialmover_scores[-1]
        self.mut_final_score = self.mut_trialmover_scores[-1]
        if self.verbose >= 1:
            print "======================="
            print "%s FINAL SCORE: w/t: %f  Mut: %f" % (self.repr_str,self.ref_final_score,self.mut_final_score)
        #print "ACCEPTED: w/t: %f  Mut: %f" % (self.ref_mover.num_accepts(),self.mut_mover.num_accepts())
        #return (wt_pose,mut_pose)
        self.ddG = self.scorefn(self.mut_pose) - self.scorefn(self.ref_pose)

    def get_result(self):
        result_dict = {'start_pose_pdb':self.start_pose_pdb,
                       'start_pose_name':self.start_pose_name,
                       'Pose_residue':self.residue,
                       'PDB_residue':self.start_pose.pdb_info().pose2pdb(self.residue),
                       'refAA':self.start_pose.residue(self.residue).name3(),
                       'mutAA':self.AA,
                       'replicate':self.replicate,
                       'ref_start_score':self.ref_start_score,
                       'mut_start_score':self.mut_start_score,
                       'ref_final_score':self.ref_final_score,
                       'mut_final_score':self.mut_final_score,
                       'ddG':self.ddG,
                       'packing_rounds':self.rnd,
                       'packing_steps':self.n_pack_steps,
                       'minimize_steps':self.n_min_steps,
                       'run_time':self.run_time}
        return result_dict

    def dump_ref_pdb(self,outfile):
        self.ref_pose.dump_pdb(outfile)

    def dump_mut_pdb(self,outfile):
        self.mut_pose.dump_pdb(outfile)

class MultiMutantPackerJob(AbstractPackerJob):

    def __init__(self,start_pose_pdb,mutation_list,replicate,convergence_fn=None,\
                 conv_threshold=0.1,repack_radius=10,scorefn="mm_std",\
                 mintype="dfpmin_armijo_nonmonotone",n_pack_steps=3,n_min_steps=1,max_rounds=100,min_cst_sd=None,min_restrict_radius=False,restrict_to_chain=False,PDB_res=False,verbose=1,constraint_file=None):

        #self.chain = chain
        self.mutation_list = mutation_list
        self.restrict_packing_to_chain = None
        self.min_cst_sd = min_cst_sd
        self.min_restrict_radius = min_restrict_radius
        if restrict_to_chain:
            self.restrict_packing_to_chain = [x[1] for x in mutation_list]

        #print >> sys.stderr, "MultiMutantPackerJob __init__(): mutation_list=%s restrict_to_chain=%s restrict_packing_to_chain=%s" % (str(self.mutation_list),str(restrict_to_chain),str(self.restrict_packing_to_chain))

        AbstractPackerJob.__init__(self,convergence_fn=convergence_fn,\
                                       conv_threshold=conv_threshold,\
                                       repack_radius=repack_radius,\
                                       scorefn=scorefn,\
                                       mintype=mintype,\
                                       max_rounds=max_rounds,\
                                       restrict_to_chain=self.restrict_packing_to_chain,\
                                       verbose=verbose,\
                                       constraint_file=constraint_file)
        
        # Mutant definition - residue to mutate, AA to mutate to, and replicate number
        self.start_pose_pdb = start_pose_pdb
        self.start_pose_name = os.path.split(start_pose_pdb)[1]
        start_pose_in = pyrosetta.pose_from_file(self.start_pose_pdb)
        self.start_pose = pyrosetta.Pose()
        self.start_pose.assign(start_pose_in)
        self.mut_str = None
        self.repr_str = None
        self.pdb_res = PDB_res

        #self.residue = None

        self.n_pack_steps = n_pack_steps
        self.n_min_steps = n_min_steps
        self.replicate = replicate
        self.ref_pose = pyrosetta.Pose()
        self.mut_pose = pyrosetta.Pose()

        self.setup_poses()

        if self.constraint_file:
            self.add_constraints_to_scorefxn()
            self.add_constraints_to_pose(self.ref_pose,self.constraint_file)
            self.add_constraints_to_pose(self.mut_pose,self.constraint_file)

        # Store initial pose scores
        self.raw_start_score = self.scorefn(self.start_pose)
        self.ref_start_score = self.scorefn(self.ref_pose)
        self.mut_start_score = self.scorefn(self.mut_pose)
        self.ref_final_score = None
        self.mut_final_score = None
        
        self.ref_min_mover = self.make_minmover(mintype,self.ref_min_movemap,self.ref_pose)
        self.mut_min_mover = self.make_minmover(mintype,self.mut_min_movemap,self.mut_pose)
            
        self.ref_pack_mover = pyrosetta.rosetta.protocols.simple_moves.PackRotamersMover(self.scorefn,self.ref_packertask)
        self.mut_pack_mover = pyrosetta.rosetta.protocols.simple_moves.PackRotamersMover(self.scorefn,self.mut_packertask)

        self.ref_mover = self.make_CompoundMover([self.ref_pack_mover,self.ref_min_mover],[self.n_pack_steps,self.n_min_steps])
        self.mut_mover = self.make_CompoundMover([self.mut_pack_mover,self.mut_min_mover],[self.n_pack_steps,self.n_min_steps])
        
        self.ref_trialmover_scores = []
        self.mut_trialmover_scores = []
        
        # Keep track of where we are in the run
        self.started = False
        self.finished = False
        self.rnd = 0
        self.ddG = None
        self.start_time = None
        self.end_time = None      

    def setup_poses(self):

        repr_str = ""

        cur_ref_pose = pyrosetta.Pose()
        cur_ref_pose.assign(self.start_pose)
        cur_mut_pose = pyrosetta.Pose()
        cur_mut_pose.assign(self.start_pose)

        pack_residues = []

        for (residue,chain,AA) in self.mutation_list:
            
            #print >> sys.stderr, "====> RESIDUE IN: %s, %d, %s, %s, PDB_res?:%s" % (self.start_pose.residue(residue).name(),residue,chain,AA,str(self.pdb_res)) 
            if self.pdb_res:
                residue = self.start_pose.pdb_info().pdb2pose(chain,residue)
            #print >> sys.stderr, "====> RESIDUE RESULT: %s, %d, %s, %s, PDB_res?:%s" % (self.start_pose.residue(residue).name(),residue,chain,AA,str(self.pdb_res)) 

            if residue <= 0:
                raise PyRosettaError("MISSING RESIDUE %s %d %s %s" % (self.start_pose_pdb,residue,chain,AA))

            repr_str += "%s_%d_%s-" % (self.start_pose.residue(residue).name(),residue,AA)
            cur_ref_pose.assign(self.mutate_aa(cur_ref_pose,residue,cur_ref_pose.residue(residue).name3()))
            cur_mut_pose.assign(self.mutate_aa(cur_mut_pose,residue,AA))

            pack_residues.extend(self.residue_CAs_in_radius(cur_ref_pose,residue,self.repack_radius))

        self.repr_str = repr_str[:-1] + ";r%d" % (self.replicate,)
        self.mut_str = repr_str[:-1]
        self.ref_pose.assign(cur_ref_pose)
        self.mut_pose.assign(cur_mut_pose)

        if self.min_cst_sd:
            pyrosetta.rosetta.core.scoring.constraints.add_coordinate_constraints(self.ref_pose,self.min_cst_sd,False)
            pyrosetta.rosetta.core.scoring.constraints.add_coordinate_constraints(self.mut_pose,self.min_cst_sd,False)

        # Set up PackerTasks and Movers 
        self.ref_packertask = self.make_packer_task_with_residues(self.ref_pose,pack_residues)
        self.mut_packertask = self.make_packer_task_with_residues(self.mut_pose,pack_residues)        
        
        self.ref_min_movemap = None
        self.mut_min_movemap = None

        if self.min_restrict_radius:
            self.ref_min_movemap = self.build_movemap(self.ref_pose,restrict_residues=pack_residues)
            self.mut_min_movemap = self.build_movemap(self.mut_pose,restrict_residues=pack_residues)
        else:
            self.ref_min_movemap = self.build_movemap(self.ref_pose)
            self.mut_min_movemap = self.build_movemap(self.mut_pose)

                                 
    def run(self):
        # Grab the score function from w/t TrialMover object
        # Might want to check that same as mut score fxn in the future
        self.rnd = 1
        self.started = True
        self.start_time = time.time()
    
        while ((not self.cnv_fn(self.ref_trialmover_scores)) or (not self.cnv_fn(self.mut_trialmover_scores))) and (self.rnd <= self.max_rounds):
            self.mut_mover.apply(self.mut_pose)
            self.ref_mover.apply(self.ref_pose)

            self.mut_trialmover_scores.append(self.scorefn(self.mut_pose))
            self.ref_trialmover_scores.append(self.scorefn(self.ref_pose))

            if self.verbose >= 2:
                print "%s ROUND %d SCORES: w/t: %f Mut: %f" % (self.repr_str,self.rnd,self.ref_trialmover_scores[-1],self.mut_trialmover_scores[-1])
            self.rnd += 1
        
        self.end_time = time.time()
        self.run_time = self.end_time - self.start_time
        self.ref_final_score = self.ref_trialmover_scores[-1]
        self.mut_final_score = self.mut_trialmover_scores[-1]
        if self.verbose >= 1:
            print "======================="
            print "%s FINAL SCORE: w/t: %f  Mut: %f" % (self.repr_str,self.ref_final_score,self.mut_final_score)
        #print "ACCEPTED: w/t: %f  Mut: %f" % (self.ref_mover.num_accepts(),self.mut_mover.num_accepts())
        #return (wt_pose,mut_pose)
        self.ddG = self.scorefn(self.mut_pose) - self.scorefn(self.ref_pose)

    def get_result(self):
        result_dict = {'start_pose_pdb':self.start_pose_pdb,
                       'start_pose_name':self.start_pose_name,
                       'mutation_list':self.mut_str, 
                       'replicate':self.replicate,
                       'ref_start_score':self.ref_start_score,
                       'mut_start_score':self.mut_start_score,
                       'ref_final_score':self.ref_final_score,
                       'mut_final_score':self.mut_final_score,
                       'ddG':self.ddG,
                       'packing_rounds':self.rnd,
                       'packing_steps':self.n_pack_steps,
                       'minimize_steps':self.n_min_steps,
                       'run_time':self.run_time}
        return result_dict

    def dump_ref_pdb(self,outfile):
        self.ref_pose.dump_pdb(outfile)

    def dump_mut_pdb(self,outfile):
        self.mut_pose.dump_pdb(outfile)

    

class PackSinglePoseJob(AbstractPackerJob):
    
    def __init__(self,in_pdb,MPI_rank,convergence_fn=None,\
                 conv_threshold=0.1,repack_radius=10,scorefn="mm_std",\
                 mintype="dfpmin_armijo_nonmonotone",n_pack_steps=1,n_min_steps=1,max_rounds=100,restrict_to_chain=None,kT=1.0,MCmover=True,verbose=1,constraint_file=None):
        
        AbstractPackerJob.__init__(self,convergence_fn,conv_threshold,repack_radius,scorefn,mintype,max_rounds,restrict_to_chain,constraint_file=constraint_file)
        
        try:
            in_pose = pyrosetta.pose_from_file(in_pdb)
        except AttributeError:
            in_pose = rosetta.pose_from_pdb(in_pdb)
        self.kT = kT
        self.pose = pyrosetta.Pose()
        self.pose.assign(in_pose)

        if self.constraint_file:
            self.add_constraints_to_scorefxn()
            self.add_constraints_to_pose(self.pose,self.constraint_file)

        self.scores = [self.scorefn(self.pose),]
        self.rnd = 0
        self.rank = MPI_rank
        self.n_pack_steps = n_pack_steps
        self.n_min_steps = n_min_steps
        self.verbose=verbose
        
        self.min_mover = self.make_minmover(self.mintype,self.pose)
        self.packer_task = self.make_packer_task_with_residues(self.pose)
        self.mover = None
        if MCmover:
            self.mover = self.make_packmin_mover(self.pose,self.packer_task,self.min_mover,self.n_pack_steps,self.n_min_steps,kT=self.kT)
        else:
            pack_mover = pyrosetta.rosetta.protocols.simple_moves.PackRotamersMover(self.scorefn,self.packer_task)
            self.mover = self.make_CompoundMover([pack_mover,self.min_mover],[self.n_pack_steps,self.n_min_steps])
            
    def pack_pose(self):
        if self.verbose >= 1:
            print "============ Packing Pose Beginning, proccess %d ===========" % (self.rank,)
            print "STARTING SCORE: %f " % (self.scores[0])
        while ((not self.cnv_fn(self.scores)) and (self.rnd <= self.max_rounds)):        
            if self.verbose >= 2:
                print "Applying Mover in process %d, round %d..." % (self.rank,self.rnd)
            #print seqmover
            #print trialmover
            self.mover.apply(self.pose)
            self.scores.append(self.scorefn(self.pose))
            if self.verbose >= 2:
                print "round %d score: %f" % (self.rnd,self.scores[-1])
            self.rnd += 1
        
        if self.verbose >= 1:
            print "============ Pose Finished, process %d ===========" % (self.rank,)
            print "FINAL SCORE: %f " % (self.scores[-1])
#
    def dump_pose(self,outfile):
        self.pose.dump_pdb(outfile)


class FastRelaxPoseJob(AbstractPackerJob):
    
    def __init__(self,in_pdb,MPI_rank,scorefn,restrict_to_chain=None,constraint_file=None):
        
        AbstractPackerJob.__init__(self,scorefn=scorefn,restrict_to_chain=restrict_to_chain,verbose=1,constraint_file=constraint_file)
        
        in_pose = pyrosetta.pose_from_file(in_pdb)

        self.pose = rosetta.Pose()
        self.pose.assign(in_pose)
        if self.constraint_file:
            self.add_constraints_to_scorefxn(constraint_types=[pyrosetta.rosetta.core.scoring.constraints.atom_pair_constraint, \
                                                               pyrosetta.rosetta.core.scoring.constraints.angle_constraint, \
                                                               pyrosetta.rosetta.core.scoring.constraints.dihedral_constraint])
            self.add_constraints_to_pose(self.pose,self.constraint_file)
        self.start_score = self.scorefn(self.pose)
        self.rank = MPI_rank
        self.rlx_mover = pyrosetta.rosetta.protocols.relax.FastRelax(self.scorefn)
        self.rlx_mover.set_movemap(self.build_movemap(self.pose))
        self.final_score = None

        
            
    def pack_pose(self):
        if self.verbose >= 1:
            print "============ Packing Pose Beginning, proccess %d ===========" % (self.rank,)
            print "STARTING SCORE: %f " % (self.scorefn(self.pose))
        self.rlx_mover.apply(self.pose)
        if self.verbose >= 1:
            print "============ Pose Finished, process %d ===========" % (self.rank,)
            print "FINAL SCORE: %f " % (self.scorefn(self.pose))
        self.final_score=self.scorefn(self.pose)

    def dump_pose(self,outfile):
        self.pose.dump_pdb(outfile)

if __name__ == "__main__":
    _main(sys.argv[1:])

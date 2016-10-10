#!/usr/bin/env python
import os,sys,time,random
import numpy as np
#from rosetta import *
import rosetta

nAAs = ['ALA','ASP','LEU','ILE','VAL','GLY','SER','THR','PRO','GLU','ASN','GLN','LYS','ARG','TRP','PHE','TYR','HIS','CYS','MET']
nsAAs = ['PRK','ACK','3IY','NOY','AZF','A69','B36','NBY'] # A69 = 3-aminotyrosine, B36 = 5-hydroxytryptophan
NSAAS_PATCH = {'NBY':{'cognateAA':'TYR',
                      'type':rosetta.VariantType.C2_AMINO_SUGAR},
               'PRK':{'cognateAA':'LYS',
                      'type':rosetta.VariantType.C3_AMINO_SUGAR},
               'ACK':{'cognateAA':'LYS',
                      'type':rosetta.VariantType.ACETYLATION}}


class MutagenesisExperimentRunner():

    def __init__(self,start_pose_pdbs,rosetta_init_options,comm,residue_list=[],AA_list=[],nreps=50,restrict_to_chain=False,max_pack_rounds=25,PDB_res=False,dump_ref_pdb=False,dump_mut_pdb=False,pdb_base=""):

        self.start_pose_pdbs = start_pose_pdbs
        self.rosetta_init_options = rosetta_init_options
        self.max_pack_rounds = max_pack_rounds
        self.restrict_to_chain = restrict_to_chain
        self.pdb_res = PDB_res
        self.dump_ref_pdb = dump_ref_pdb
        self.dump_mut_pdb = dump_mut_pdb
        self.pdb_base = pdb_base
        self.comm = None
        self.size = 1
        self.rank = 0
        if not comm:
            rosetta.init(extra_options=self.rosetta_init_options)
        else:
            self.comm = comm
            self.size = comm.Get_size()
            self.rank = comm.Get_rank()
        # keep a list of Job objects
            rosetta.mpi_init(extra_options=self.rosetta_init_options)
        
        if self.rank == 0:
            self.jobs = self.build_job_list(residue_list,AA_list,nreps)
        else:
            self.jobs = None
        
        if self.rank ==0:
            self.final_results = []
        else:
            self.final_results = None
        
        self.result = []
        # Pandas table to track results
        #self.data_table = self.init_table(residue_list,AA_list,nreps)

    def build_job_list(self,residues,AAs,replicates,shuffle_jobs=True):
        # build job objects for the entire run
        jobs = []
        for start_pdb in self.start_pose_pdbs:
            for (res,chain) in residues:
                for AA in AAs:
                    for rep in np.arange(replicates):
                        jobs.append((start_pdb,res,chain,AA,rep))
        print "+++++++ Sending %d jobs to %d MPI child processes" % (len(jobs),self.size)
        if shuffle_jobs:
            random.shuffle(jobs) 
        job_div = int(np.floor(len(jobs) / self.size))
        job_rem = len(jobs) - (job_div * self.size)
        job_div_procs = self.size - job_rem

        split_jobs = [jobs[i:i+job_div] for i in range(0,job_div_procs*job_div,job_div)]
        split_jobs.extend([jobs[i:i+(job_div+1)] for i in range(job_div_procs*job_div,len(jobs),job_div+1)])
        print >> sys.stderr, "*********************  n_jobs: %d job_div: %d job_rem: %d" % (len(jobs),job_div,job_rem)
        #if job_rem > 0:
        #    split_jobs.append(jobs[-job_rem:])
        print >> sys.stderr,  "*********** split_jobs: %d" % (len(split_jobs),)
        print >> sys.stderr, [len(x) for x in split_jobs]
        return split_jobs
    
    def scatter_job(self):
        local_jobs = self.comm.scatter(self.jobs,root=0)
        for (i,job_spec) in enumerate(local_jobs):
            print "===== Process %d, running job %s [ %d / %d ]" % (self.rank,str(job_spec),i,len(local_jobs))
            try:
                ddg_job = MutantddGPackerJob(*job_spec,max_rounds=self.max_pack_rounds,restrict_to_chain=self.restrict_to_chain,PDB_res=self.pdb_res)
            
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
            print >> sys.stderr, self.final_results
    
    def run_all_jobs_serial(self,outfile='ddg_out.txt'):
        for (i,job_spec) in enumerate(self.jobs):
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
    
class AbstractPackerJob():
    
    def __init__(self,convergence_fn=None,\
                 conv_threshold=0.1,repack_radius=10,scorefn="mm_std",\
                 mintype="dfpmin_armijo_nonmonotone",max_rounds=100,restrict_to_chain=None):
        self.convergence_fn = convergence_fn
        self.cnv_fn = None
        self.conv_threshold = conv_threshold
        if not convergence_fn or convergence_fn == 'stdev': # defaults to stdev for now, may add more later
            self.cnv_fn = self.std_dev_threshold_fn_builder(conv_threshold)
        self.max_rounds = max_rounds    
        self.repack_radius = repack_radius
        self.mintype = mintype
        self.restrict_to_chain = restrict_to_chain
        # a PyRosetta ScoreFunction object, defaults to mm_std    
        self.scorefn = rosetta.create_score_function(scorefn)        

        
    def packer_task_repack_in_radius(self,pose,residue,radius):
        """
        Build a packer task that repacks all residues with CAlphas within distance radius
    
        Might want to remake to take centroid side-chains?  Looks like SwitchResidueTypeMover has trouble w/ nsAAs...
        """
        pack_residues = self.residue_CAs_in_radius(pose,residue,radius,self.restrict_to_chain)
        return self.make_packer_task_with_residues(pose,pack_residues)

    def make_packer_task_with_residues(self,pose,residues=None):
        """
        Builds a packer task with the specified residues activated for repacking.
    
        Did this to avoid PackerTask.temporarily_* methods which apparently we're not supposed to use
        """
        packer_task = rosetta.standard_packer_task(pose) 
        packer_task.restrict_to_repacking()
       
        if residues != None:
            # Vector1 doesn't translate booleans correctly in 
            # TACC version of PyRosetta; need to build utility.vector1_bool() directly
            #pack_residues = Vector1([x in residues for x in range(1,pose.n_residue())])
            pack_residues = rosetta.utility.vector1_bool()
            for i in range(1,pose.n_residue() + 1):
                if self.restrict_to_chain and pose.pdb_info().chain(i) != self.restrict_to_chain:
                    continue
                else:
                    pack_residues.append(i in residues)
            packer_task.restrict_to_residues(pack_residues)
        
        return packer_task

    def residue_CAs_in_radius(self,pose,centerAA,radius,restrict_to_chain=None):
        """
        Get a list of residues with C-alpha atoms within 'radius' distance from centerAA'a C-alpha
        """
        centerAA_CA = pose.residue(centerAA).xyz('CA')
        repack_residues = []
    
        for i in range(1,pose.n_residue()):
            if restrict_to_chain and pose.pdb_info().chain(i) != restrict_to_chain:
                continue
            test_CA = pose.residue(i).xyz('CA')
            displacement = centerAA_CA - test_CA
            distance = displacement.norm
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
        
        seqmover = rosetta.SequenceMover()
        packmover = rosetta.PackRotamersMover(self.scorefn,packertask)

        for i in range(0,n_packing_steps):
            #pack_repmover = rosetta.RepeatMover(packmover,n_packing_steps)
            seqmover.add_mover(packmover)
        
        for i in range(0,n_minimize_steps):
            #min_repmover = rosetta.RepeatMover(minmover,n_minimize_steps)
            seqmover.add_mover(minmover)
         

        #print >> sys.stderr, seqmover
        
        mc = rosetta.MonteCarlo(pose,self.scorefn,kT)
        packmin_trial_mover = rosetta.TrialMover(seqmover,mc)
    
        return seqmover
    
    def std_dev_threshold_fn_builder(self,threshold,last_n=5):
        def stop(scores):
            passed = False
            if len(scores) > last_n:
                last_n_stdev = np.std(scores[-last_n:])
                passed = last_n_stdev < threshold
                print "std dev for scores %s: %f thresh: %f pass?: %s" % (str(scores[-last_n:]),last_n_stdev,  threshold,str(passed))
            return passed
        return stop

    def build_movemap(self,pose):
        mm = rosetta.MoveMap()
        if self.restrict_to_chain:
            chain_residues = []
            for i in range(1,pose.n_residue()):
                if self.restrict_to_chain and pose.pdb_info().chain(i) == self.restrict_to_chain:
                    chain_residues.append(i)
            chain_min = min(chain_residues)
            chain_max = max(chain_residues)
            mm.set_bb_true_range(chain_min,chain_max)
            mm.set_chi_true_range(chain_min,chain_max)
        else:
            mm.set_bb(True)
            mm.set_chi(True)
        return mm
        
    def make_minmover(self,mintype,pose):
        # Set up MoveMap.
        mm = self.build_movemap(pose)

        # Set up a Minimization mover using the mm_std score function
        min_mover = rosetta.MinMover()
        min_mover.movemap(mm)
        min_mover.score_function(self.scorefn)
        min_mover.min_type(mintype)
        
        return min_mover

class PyRosettaError(Exception):
    pass

class MutantddGPackerJob(AbstractPackerJob):

    def __init__(self,start_pose_pdb,residue,chain,AA,replicate,convergence_fn=None,\
                 conv_threshold=0.1,repack_radius=10,scorefn="mm_std",\
                 mintype="dfpmin_armijo_nonmonotone",n_pack_steps=3,n_min_steps=1,max_rounds=100,restrict_to_chain=False,PDB_res=False):

        self.chain = chain
        self.restrict_packing_to_chain = None
        if restrict_to_chain:
            self.restrict_packing_to_chain = self.chain

        AbstractPackerJob.__init__(self,convergence_fn,conv_threshold,repack_radius,scorefn,mintype,max_rounds,self.restrict_packing_to_chain)
        
        # Mutant definition - residue to mutate, AA to mutate to, and replicate number
        self.start_pose_pdb = start_pose_pdb
        self.start_pose_name = os.path.split(start_pose_pdb)[1]
        start_pose_in = None
        # GAAAAAAAH Why do they not make these backwards-compatible >8-O!!!
        try:
            start_pose_in = rosetta.pose_from_pdb(self.start_pose_pdb)
        except AttributeError:
            start_pose_in = rosetta.pose_from_file(self.start_pose_pdb)
        self.start_pose = rosetta.Pose()
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
        self.ref_pose = rosetta.Pose()
        self.mut_pose = rosetta.Pose()
        self.ref_pose.assign(self.mutate_aa(self.start_pose,self.residue,self.start_pose.residue(self.residue).name3()))
        self.mut_pose.assign(self.mutate_aa(self.start_pose,self.residue,self.AA))

        
        # Store initial pose scores
        self.raw_start_score = self.scorefn(self.start_pose)
        self.ref_start_score = self.scorefn(self.ref_pose)
        self.mut_start_score = self.scorefn(self.mut_pose)
        self.ref_final_score = None
        self.mut_final_score = None
        
        # Set up PackerTasks and Movers 
        self.ref_packertask = self.packer_task_repack_in_radius(self.ref_pose,self.residue,self.repack_radius)
        self.mut_packertask = self.packer_task_repack_in_radius(self.mut_pose,self.residue,self.repack_radius)        
        
        self.ref_min_mover = self.make_minmover(mintype,self.ref_pose)
        self.mut_min_mover = self.make_minmover(mintype,self.ref_pose)
            
        self.ref_pack_mover = rosetta.PackRotamersMover(self.scorefn,self.ref_packertask)
        self.mut_pack_mover = rosetta.PackRotamersMover(self.scorefn,self.mut_packertask)

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

        
    def mutate_aa(self,pose,residue,aa_name):
        """
        Swap w/t AA at residue number 'residue' in 'pose' with 'ncaa_name' (3-letter code)
    
        Return a new Pose object
        
        Note that this assumes the ncaa .params and .rotlib files have been permanently added to the database
        """

        mut_pose = rosetta.Pose()
        mut_pose.assign(pose)

        res = mut_pose.residue(residue)
        # check for disulfides and correct if needed
        if (res.name() == 'CYS:disulfide') or (res.name() == 'CYD'):
            disulfide_partner = None
            try:
                disulfide_partner = res.residue_connection_partner(
                    res.n_residue_connections())
            except AttributeError:
                disulfide_partner = res.residue_connection_partner(
                    res.n_current_residue_connections())
            temp_pose = rosetta.Pose()
            temp_pose.assign(mut_pose)
            # (Packing causes seg fault if current CYS residue is not
            # also converted before mutating.)
            rosetta.change_cys_state(residue, 'CYS',\
                             temp_pose.conformation())
            rosetta.change_cys_state(disulfide_partner, 'CYS',\
                             temp_pose.conformation())
            mut_pose = temp_pose

    
            # Get a Residue object for the desired ncAA
        mut_res_type = None
        rts = None
        chm = rosetta.core.chemical.ChemicalManager.get_instance()
        
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
            #if residue == 1 or residue == mut_pose.n_residue():
            #    aa_name = aa_name + "_p"
            mut_res_type = rts.name_map(aa_name)

        if residue == 1:
            try: 
                mut_res_type = rts.get_residue_type_with_variant_added(mut_res_type,rosetta.VariantType.LOWER_TERMINUS_VARIANT)                
            except AttributeError:
                mut_res_type = rts.get_residue_type_with_variant_added(mut_res_type,"LOWER_TERMINUS")
        elif residue == mut_pose.n_residue():
            try:
                mut_res_type = rts.get_residue_type_with_variant_added(mut_res_type,rosetta.VariantType.UPPER_TERMINUS_VARIANT)                
            except AttributeError:
                mut_res_type = rts.get_residue_type_with_variant_added(mut_res_type,"UPPER_TERMINUS")
        mut_res = rosetta.core.conformation.ResidueFactory.create_residue( mut_res_type )
        mut_pose.replace_residue(residue,mut_res,orient_backbone=True)
        
        return mut_pose
          
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

            print "ROUND %d SCORES: w/t: %f Mut: %f" % (self.rnd,self.ref_trialmover_scores[-1],self.mut_trialmover_scores[-1])
            self.rnd += 1
        
        self.end_time = time.time()
        self.run_time = self.end_time - self.start_time
        self.ref_final_score = self.ref_trialmover_scores[-1]
        self.mut_final_score = self.mut_trialmover_scores[-1]
        print "======================="
        print "FINAL SCORE: w/t: %f  Mut: %f" % (self.ref_final_score,self.mut_final_score)
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
        

class PackSinglePoseJob(AbstractPackerJob):
    
    def __init__(self,in_pdb,MPI_rank,convergence_fn=None,\
                 conv_threshold=0.1,repack_radius=10,scorefn="mm_std",\
                 mintype="dfpmin_armijo_nonmonotone",n_pack_steps=1,n_min_steps=1,max_rounds=100,restrict_to_chain=None,kT=1.0,MCmover=True):
        
        AbstractPackerJob.__init__(self,convergence_fn,conv_threshold,repack_radius,scorefn,mintype,max_rounds,restrict_to_chain)
        
        try:
            in_pose = rosetta.pose_from_file(in_pdb)
        except AttributeError:
            in_pose = rosetta.pose_from_pdb(in_pdb)
        self.kT = kT
        self.pose = rosetta.Pose()
        self.pose.assign(in_pose)
        self.scores = [self.scorefn(self.pose),]
        self.rnd = 0
        self.rank = MPI_rank
        self.n_pack_steps = n_pack_steps
        self.n_min_steps = n_min_steps
        
        self.min_mover = self.make_minmover(self.mintype,self.pose)
        self.packer_task = self.make_packer_task_with_residues(self.pose)
        self.mover = None
        if MCmover:
            self.mover = self.make_packmin_mover(self.pose,self.packer_task,self.min_mover,self.n_pack_steps,self.n_min_steps,kT=self.kT)
        else:
            pack_mover = rosetta.PackRotamersMover(self.scorefn,self.packer_task)
            self.mover = self.make_CompoundMover([pack_mover,self.min_mover],[self.n_pack_steps,self.n_min_steps])
            
    def pack_pose(self):
        print "============ Packing Pose Beginning, proccess %d ===========" % (self.rank,)
        print "STARTING SCORE: %f " % (self.scores[0])
        while ((not self.cnv_fn(self.scores)) and (self.rnd <= self.max_rounds)):        
            print "Applying Mover in process %d, round %d..." % (self.rank,self.rnd)
            #print seqmover
            #print trialmover
            self.mover.apply(self.pose)
            self.scores.append(self.scorefn(self.pose))
            print "round %d score: %f" % (self.rnd,self.scores[-1])
            self.rnd += 1
        
        print "============ Pose Finished, process %d ===========" % (self.rank,)
        print "FINAL SCORE: %f " % (self.scores[-1])
#
    def dump_pose(self,outfile):
        self.pose.dump_pdb(outfile)


class FastRelaxPoseJob(AbstractPackerJob):
    
    def __init__(self,in_pdb,MPI_rank,scorefn,restrict_to_chain=None):
        
        AbstractPackerJob.__init__(self,scorefn=scorefn,restrict_to_chain=restrict_to_chain)
        
        in_pose = rosetta.pose_from_file(in_pdb)

        self.pose = rosetta.Pose()
        self.pose.assign(in_pose)
        self.start_score = self.scorefn(self.pose)
        self.rank = MPI_rank
        self.rlx_mover = rosetta.FastRelax(self.scorefn)
        self.rlx_mover.set_movemap(self.build_movemap(self.pose))
        self.final_score = None

            
    def pack_pose(self):
        print "============ Packing Pose Beginning, proccess %d ===========" % (self.rank,)
        print "STARTING SCORE: %f " % (self.scorefn(self.pose))
        self.rlx_mover.apply(self.pose)
        print "============ Pose Finished, process %d ===========" % (self.rank,)
        print "FINAL SCORE: %f " % (self.scorefn(self.pose))
        self.final_score=self.scorefn(self.pose)

    def dump_pose(self,outfile):
        self.pose.dump_pdb(outfile)



def MPIJobDistributor(njobs, fun):

    comm = MPI.COMM_WORLD

    rank = comm.Get_rank()
    size = comm.Get_size()

    myjobs = []

    if rank == 0:
        jobs = range(njobs)
        jobs.extend( [None]*(size - njobs % size) )
        n = len(jobs)/size
        for i in range(size):
            queue = []  # list of jobs for individual cpu
            for j in range(n):
                queue.append(jobs[j*size+i])

            if( i == 0 ):
                myjobs = queue
            else:
                # now sending the queue to the process
                logger.info('Sending %s to node %s' % (queue, i) )
                comm.send(queue, dest=i)
    else:
        # getting decoy lists
        myjobs = comm.recv(source=0)

    logger.info('Node %s, got queue:%s' % (rank, myjobs) )

    for j in myjobs:
        if j is not None: fun(j)

if __name__ == "__main__":
    _main(sys.argv[1:])

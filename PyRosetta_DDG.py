#!/opt/apps/intel15/python/2.7.9/bin/python

######################################## 1 ########################################
# mutate_residue function taken from mutations.py sample script:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#
## @file   mutants.py
## @brief
## @author Evan H. Baugh, Johns Hopkins University
#
# mutate_residue is adapted from an original script by Sid Chaudhury
######################################## 1 ########################################



import rosetta
# from rosetta import init, pose_from_pdb, get_fa_scorefxn, standard_packer_task, Pose, MoveMap, RotamerTrialsMinMover, MinMover # TODO only import what you need to lower loading time on tacc ... does not appear to increase speed
import toolbox
import argparse
import os
import datetime
import re
import numpy
import multiprocessing
from Bio.PDB import PDBParser
from functools import partial
import sys
import shutil
import math
import random
import time

__author__ = 'ded'
"""
TODO: take .gd file input require reference to be same
TODO: blast pdb file against sequence, truncate as appropriate

"""

# Rosetta must be initialized before processing arguments to allow for global score function declaration
print "Rosetta initializing..."
rosetta.init()

parser = argparse.ArgumentParser(description="take pdb file, optimize it to given degree, read in list of mutations, optimize them to same degree")
parser.add_argument("-ref", "--reference", help="pdb file", required=True)
parser.add_argument("-o", "--output", help="Output directory for files (will be created if it doesn't exist).", default=os.curdir)
parser.add_argument("-l", "--log", help="Log file to be appended to with run information.", default="PyRosetta_DDG.py.log.txt")
parser.add_argument("-p", "--prefix", help="Output file prefix", default=str(datetime.datetime.now()).replace(" ", "_").replace(":", "/"))
parser.add_argument("-rep", "--replicates", help="how many replicate energy minimizations of the reference pdb should be attempted", default=multiprocessing.cpu_count(), type=int)
parser.add_argument("-i", "--intermediates", help="store intermediate pdb files at completion of rotamer optimization", action='store_true')
parser.add_argument("-n", "--number", help="number of optimized structures to compare. If float will take top %s of values, if integer will use that many best values, use 0 to take all optimized structures" % default.replace(r"%", r"%%"), default=1, type=float)
parser.add_argument("-s", "--score_function", help="Score function to be used", default=rosetta.get_fa_scorefxn())

parser.add_argument("--max", help="how many moves should each replicate attempt at maximum. Value of 0 will disable test.", default=1000, type=int)
parser.add_argument("--failures", help="how many times should energy optimiazations functions fail to find an improved structure before continuing", default=10, type=int)
parser.add_argument("--threshold", help="what percentage of improvement should be achieved over increment of failures. This likely relates to the minimal difference that can be observed statistically.", type=float, default=0.01)

parser.add_argument("-m", "--mutations", help="TSV file with list of mutations in coordinate system of pdb file.", default="Mutations.tsv")
parser.add_argument("-d", "--distance", help="Radius from mutated AA to be repacked", default=10.0, type=float)
parser.add_argument("-c", "--comparisons", help="number of ddg comparisons to make for each mutatation and each optimized reference file", default=100, type=int)
parser.add_argument("-v", "--verbose", action='store_true', help="Increase output")

parser.add_argument("--testing", action='store_true', help="sets constant seed for testing purposes")
args = parser.parse_args()


# rosetta re-initilaized before function read in with desired rng seeding
if args.testing:
    rosetta.init(extra_options="-constant_seed -mute basic -mute core ")  # constant seed for testing as per Refinement script. To be removed before real run
else:
    rosetta.init(extra_options="-mute basic -mute core")

# score function declared before everything else to ensure single score function used throughout
scorefxn = args.score_function  # must be declared before functions for use in functions  # TODO refactor so only scorefxn exists?


def log(string_to_print, loc_to_print=args.output.rstrip("/") + "/" + args.log):
    """
    print information to log file
    :param string_to_print: what to print
    :param loc_to_print: location to print to
    :return:
    """
    with open(loc_to_print, "a+") as output:
        print>>output, string_to_print
        if args.verbose:
            print string_to_print


# ####################################### 1 ########################################
def mutate_residue(pose, mutant_position, mutant_aa, pack_radius=args.distance, mr_scorefxn=args.score_function, proc=0):
    """
    :param pose: protein to mutate
    :param mutant_position: location to mutate
    :param mutant_aa: single letter code for new AA
    :param pack_radius: default=0.0; how far away (in Angstroms) to repack
    :param mr_scorefxn: default=get_fa_scorefxn(); scoring function to score
    Replaces the residue at  <mutant_position>  in  <pose>  with  <mutant_aa>
        and repack any residues within  <pack_radius>  Angstroms of the mutating
        residue's center (nbr_atom) using  <mr_scorefxn>
    note: <mutant_aa>  is the single letter name for the desired ResidueType

    example:
        mutate_residue(pose,30,A)
    See also:
        Pose
        PackRotamersMover
        MutateResidue
        pose_from_sequence

    :return:
    """

    # why is this here?
    # ### a MutateResidue Mover exists similar to this except it does not pack the area around the mutant residue (no pack_radius feature)
    # mutator = MutateResidue( mutant_position , mutant_aa )
    # mutator.apply( test_pose )

    if args.testing:  # need to set constant seed, but want different seeds for each process
        rosetta.init(extra_options="-constant_seed -jran %i" % proc)
    else:  # need to reinitialize to get random seed, but account for potential of multi processes to start at same time
        rosetta.init(extra_options="-use_time_as_seed -seed_offset %i" % proc)


    if not pose.is_fullatom():
        IOError('mutate_residue only works with fullatom poses')
    test_pose = rosetta.Pose()
    test_pose.assign(pose)


    # this block is accomplished by setting default value of mr_scorefxn to desired function
    # # create a standard scorefxn by default
    # if not mr_scorefxn:
    #     # mr_scorefxn = create_score_function( 'standard' )  # @ded: throws error, D070_Refinement.py script shows change to:
    #     mr_scorefxn = get_fa_scorefxn()  # appears to accomplish the same

    task = rosetta.standard_packer_task(test_pose)

    # the Vector1 of booleans (a specific object) is needed for specifying the mutation, this demonstrates another more direct method of setting PackerTask options for design
    aa_bool = rosetta.utility.vector1_bool()

    # PyRosetta uses several ways of tracking amino acids (ResidueTypes) the numbers 1-20 correspond individually to the 20 proteogenic amino acids aa_from_oneletter returns the integer representation of an amino acid from its one letter code
    # convert mutant_aa to its integer representation
    mutant_aa = rosetta.aa_from_oneletter_code(mutant_aa)

    # TODO is this the best way to cause a mutation? Specifically, creating a list of 20F and 1T bools
    # mutation is performed by using a PackerTask with only the mutant amino acid available during design to do this, construct a Vector1 of booleans indicating which amino acid (by its numerical designation, see above) to allow
    for i in range(1, 21):
        # in Python, logical expression are evaluated with priority, thus the line below appends to aa_bool the truth (True or False) of the statement i == mutant_aa
        aa_bool.append(i == mutant_aa)

    # modify the mutating residue's assignment in the PackerTask using the Vector1 of booleans across the proteogenic amino acids
    task.nonconst_residue_task(mutant_position).restrict_absent_canonical_aas(aa_bool)

    # prevent residues from packing by setting the per-residue "options" of the PackerTask
    #center = rosetta.pose.residue( mutant_position ).nbr_atom_xyz()  # throws error
    center = test_pose.residue( mutant_position ).nbr_atom_xyz()

    for i in range(1, pose.total_residue() + 1):
        # print i, "\t", center.distance_squared(test_pose.residue(i).nbr_atom_xyz()), "\t", test_pose.residue(i).name3()
        # only pack the mutating residue and any within the pack_radius
        # if not i == mutant_position or center.distance_squared(test_pose.residue(i).nbr_atom_xyz()) > pack_radius**2:  # this is the worst error i've ever had to debug NEVER EVER UNCOMMENT WITHOUT TALKING TO DAN
        if center.distance_squared(test_pose.residue(i).nbr_atom_xyz()) > pack_radius**2:  # mutated residue will always have a distance of 0,  which is less than the square of any real number
            task.nonconst_residue_task(i).prevent_repacking()
        else:
            if i == mutant_position:
                continue
            else:
                task.nonconst_residue_task(i).restrict_to_repacking()
        # else:  # this is the key block that ultimately proved the original if statement was incorrect as it only allowed the targeted residue to be repacked regardless of distance given
        #     print "THIS HAPPENED", i, mutant_position

    # apply the mutation and pack nearby residues
    before = test_pose.sequence()
    packer = rosetta.PackRotamersMover(mr_scorefxn, task)
    #print packer
    packer.apply(test_pose)  # mutation finally happens
    rotamer = mr_scorefxn(test_pose)


    # Make sure that the sequence has actually mutated, and that only 1 AA has changed
    # TODO remove testing block?
    differences = 0
    for i, b in enumerate(before):
        if b == test_pose.sequence()[i]:
            continue
        differences += 1
    assert differences == 1, "%i differences in the before and after sequence. if 0, mutation was not made, if more than 1 multiple bases mutated.\nBefore:\n%s\nAfter:\n%s" % (differences, before, test_pose.sequence())

    #print test_pose

    #attempt to free side chains
    test_pose.dump_pdb("reset_pose.pdb")
    rosetta.pose_from_pdb(test_pose, "reset_pose.pdb")


    movemap = rosetta.MoveMap()
    movemap.set_bb(True)  # change backbone confirmation
    #movemap.set_bb_true_range(1, pose.total_residue() + 1)
    movemap.set_chi(True)  # change side chain confirmation
    #movemap.set_chi_true_range(1, pose.total_residue() + 1)
    movemap.set_nu(True)
    movemap.set_branches(True)
    min_mover = rosetta.MinMover()
    min_mover.movemap(movemap)
    min_mover.score_function(mr_scorefxn)
    min_mover.min_type('dfpmin_armijo_nonmonotone')
    #print movemap
    #print min_mover
    # single move
    print min_mover
    min_mover.apply(test_pose)


    # # loop through multple moves to see if reaches same value
    # for _ in xrange(10):
    #     #rosetta.init(extra_options="-use_time_as_seed")
    #     min_mover.apply(test_pose)

    # # use small mover instead to see if deterministic
    # movemap = rosetta.MoveMap()
    # movemap.set_bb(True)
    # movemap.set_chi(True)
    # small_mover = rosetta.SmallMover(movemap, 1.0, 1)
    # small_mover.apply(test_pose)
    final = mr_scorefxn(test_pose)

    print [rotamer, final]

    return rotamer, final, test_pose

    # return test_pose

    # score = mr_scorefxn(test_pose)
    # return score

# ####################################### 1 ########################################


def minimize_energy(proc, pdb_filename=None, failures=args.failures, max_attempts=args.max, me_scorefxn=args.score_function):
    """
    :param proc: current process for output purposes
    :param pdb_filename: pdb file to be minimized
    :param failures: how many unsucessful attempts
    :param max_attempts: cap on how many times structure should be attempted to improve
    :param me_scorefxn: what scoreing function should be used
    Given a pdb file, minimize energy
    """

    # 0a. Check if structure replicate has already been completed. This is required if program times out.
    optimized_name = args.output.rstrip("/") + "/" + args.prefix + "_Optimized_Replicate_" + str(proc).zfill(len(str(args.replicates))) + ".pdb"
    if os.path.isfile(optimized_name):
        log("\t%s file already exists, therefore not re-optimizing." % optimized_name)  # messy output with multiprocessing
        optimized_structure = rosetta.Pose()
        rosetta.pose_from_pdb(optimized_structure, optimized_name)  # structure is being read back in, therefore no concern about minor differences in read/write
        return me_scorefxn(optimized_structure)

    # 0b. Initialize rosetta rng to ensure different outputs while using multiprocessing
    if args.testing:  # need to set constant seed, but want different seeds for each process
        rosetta.init(extra_options="-constant_seed -jran %i" % proc)
    else:  # need to reinitialize to get random seed, but account for potential of multi processes to start at same time
        rosetta.init(extra_options="-use_time_as_seed -seed_offset %i" % proc)

    # 0c. By Default, pdb_filename has default value of none so that functions can be defined before manipulation of reference, here, pdb is renamed to reference unless otherwise specified as could be in other instances
    if pdb_filename is None:
        pdb_filename = args.reference

    # 1. create a pose from the desired PDB file, and copies of it for the eventual output, and optimization
    pose = rosetta.Pose()
    loop_pose = rosetta.Pose()
    output_pose = rosetta.Pose()
    rosetta.pose_from_pdb(pose, pdb_filename)
    loop_pose.assign(pose)
    output_pose.assign(pose)

    # 2. Determine starting score and generate lists for loops
    initial_score = me_scorefxn(pose)  # Store starting energy mostly for testing improvements
    all_scores = [initial_score]
    # TODO implement optional testing of skipping if intermediate has already been created.

    # 3. optimize side chains only.
    log("Repacking Side Chains until < 1% of total improvement achieved in last 10 attempts. A minimum of 11 rounds of optimization will occur.")  # messy output with multiprocessing
    threshold_check = True
    pack_counter = 0
    while threshold_check:  # After much testing, rotamer optimization to be repeated untill < 1% of total improvement over last 10 attempts
        pack_counter += 1  # for bookkeeping
        if args.verbose:
            log("%s\t%s\t%s\t%s" % (proc, pack_counter, len(all_scores), all_scores))  # messy output with multiprocessing

        # Substitute "optizmized" rotamer angles or existing angle for each position in protein avoiding clashes
        task = rosetta.standard_packer_task(loop_pose)
        task.restrict_to_repacking()
        task.or_include_current(True)
        pack_rotamers_mover = rosetta.RotamerTrialsMinMover(me_scorefxn, task)
        pack_rotamers_mover.apply(loop_pose)

        all_scores.append(me_scorefxn(loop_pose))  # score and store improved score. Previous testing showed this is always an improvement
        assert all_scores[-1] < all_scores[-2], "Initial testing showed that rotamer optimization always decreased energy. This was just violated.\nCurrent score:\t%s\nPrevious score:\t%s\nTrajectory:\n" % (all_scores[-1], all_scores[-2]) + "\n".join(map(str, all_scores))
        try:  # determine if threshold still valid, will initially throw an index error as the length of all_scores is <10
            if ((all_scores[-10] - all_scores[-1]) / (all_scores[1] - all_scores[-1])) < 0.01:  # percentage calculated after first round of optimization as results in gigantic energy drop.
                if pack_counter > 11:  # len check necessary to make sure 10 improvements made after initial drop
                    threshold_check = False
        except IndexError:
            pass

    log("Repacking Side Chains complete after %i total attempts. Current score is: %f" % (pack_counter, all_scores[-1]))  # messy output with multiprocessing
    pack_counter += 1  # +1 to pack_counter required because initial score added before new scores calculated
    rotamer_score = all_scores[-1]  # last score is best, previous assert after score append ensures truth
    output_pose.assign(loop_pose)  # output best structure
    assert rotamer_score == me_scorefxn(output_pose), "Score of rotamer optimization changed somehow."
    assert len(all_scores) == pack_counter, "Score lost somewhere? %i \t %i" % (len(all_scores), pack_counter)

    if args.intermediates:  # block may be eliminated as this is being used to test importance of rotamer optimization
        intermediate_name = args.output.rstrip("/") + "/" + args.prefix + "_Rotamer_replicate_" + str(proc).zfill(len(str(args.replicates))) + ".pdb"
        if os.path.isfile(intermediate_name):
            # TODO implement warning/confirmation/assertion? that intermediate file is to be over written, or set intermediate check for file already existing
            assert not os.path.isfile(optimized_name), "%s Final file already exists, but you have regenerated an intermediate, shouldn't happen."
        output_pose.dump_pdb(intermediate_name)

    # 4a. create mover and initialize values
    log("Minimizing overall energy of backbone and side chains for maximum of %i attempts or until %i attempts fail to yield improved structure of at least %f" % (max_attempts, failures, args.threshold))  # messy output with multiprocessing
    movemap = rosetta.MoveMap()
    movemap.set_bb(True)  # change backbone confirmation
    movemap.set_chi(True)  # change side chain confirmation
    min_mover = rosetta.MinMover()
    min_mover.movemap(movemap)
    min_mover.score_function(me_scorefxn)
    min_mover.min_type('dfpmin_armijo_nonmonotone')

    output_score = rotamer_score  # reset score for book keeping
    minimization_counter = 0  # for book keeping
    threshold_check = True
    failure_count = 0
    optimal_scores = []

    # 4b. optimize backbone and rotamers

    #  TODO: determine if https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/IteratedConvergenceMover is a better method than current threshold
    while threshold_check and (failure_count < failures) and (minimization_counter <= max_attempts):
        minimization_counter += 1  # for book keeping

        min_mover.apply(loop_pose)  # attempt to minimize
        current_score = me_scorefxn(loop_pose)  # calculate once and store as separate variable as it will be used multiple times
        all_scores.append(current_score)

        if current_score < output_score:  # score has improved
            if args.verbose:
                log("Lower energy conformation with %s energy found on the %i improvement attempted after %i consecutive failed attempts to improve the energy." % (current_score, minimization_counter, failure_count))  # messy output with multiprocessing
            output_pose.assign(loop_pose)  # change output_pose to new optimal
            output_score = current_score  # change output score to
            optimal_scores.append(current_score)
            failure_count = 0  # reset failure count as new optimal found

            # determine if threshold still valid
            try:  # try less expensive than if to check range assuming usually true. differential controls number of attempts, but if threshold reached quickly only triggered limited times anyway
                threshold_percent = ((optimal_scores[-1] - optimal_scores[-failures]) / (min(all_scores[pack_counter + 1:]) - max(all_scores[pack_counter + 1:])))  # threshold based on minimization improvement alone so offset by pack_counter
                if threshold_percent < args.threshold:
                    if minimization_counter > args.failures:  # at least args.failure attempts have been made
                        threshold_check = False
            except IndexError:
                pass
        else:  # current loop_pose confirmation is not an improvement
            failure_count += 1
            loop_pose.assign(output_pose)  # reset structure to optimal. min_mover allows for slight increase in E itself

    # 4c. Optimization complete, log information about how/why it completed.
    assert [not threshold_check, (failure_count >= failures), (minimization_counter >= max_attempts)].count(True) == 1, "Multiple or unknown exiting conditions. This has not previously been encountered and is assumed to be a problem.\nAll scores:\t%s\nThreshold_check:\t%s\nfailure_count:\t%s\nminimization_counter:\t%s" % (all_scores, threshold_check, failure_count, minimization_counter)
    assert output_score == me_scorefxn(output_pose), "Score of energy minimization changed somehow"
    if not threshold_check:  # minimization stopped because of threshold value
        log("Energy minimization complete after %i total attempts because improvement over the last %i attempts has dropped to %f which is less than required value of %f\n\tCurrent score is: %f" % (minimization_counter, failures, threshold_percent, args.threshold, output_score))
    elif failure_count >= failures:  # minimization stopped because of consecutive failures to improve
        log("Energy minimization complete after %i total attempts, and %i consecutive attempts failed to improve beyond the current score: %f" % (minimization_counter, failure_count, output_score))
    elif minimization_counter >= max_attempts: # minimization stopped because of maximum attempts to continue improving
        log("Energy minimization complete after %i total attempts, the maximum allowed based on user input. Last attempts had following scores: %s\n\tCurrent score is: %f" % (minimization_counter, all_scores[-failures:], output_score))

    # 5. Finish
    # Small discrepancy in score after write to pdb file and read back in as pose. https://www.rosettacommons.org/comment/8734#comment-8734 (2nd paragraph). therefore need to write to file, and read back in and score to ensure downstream functionality.
    # optimized_name generated a start of function
    assert not os.path.isfile(optimized_name), "%s File already exists. This means other process output same file name while this process was running as first check failed. This could be a problem with multiprocessing" % optimized_name
    output_pose.dump_pdb(optimized_name)
    rosetta.pose_from_pdb(output_pose, optimized_name)
    scores = [initial_score, rotamer_score, scorefxn(output_pose)]
    if args.verbose:
        log("\nOptimization complete.\nInitial Score:\t%s\nRotamer Score:\t%s\nFinal Score:\t%s" % (scores[0], scores[1], scores[2]))
        log("\nEnergy trajectory:\n" + "\n".join(map(str, all_scores)))
    return scores[-1]  # return best score for downstream applications


# def ddg_score(proc, sref_pose, sref_packer, smut_pose, smut_packer, sminimizer, s_scorefxn=args.score_function):
#     """
#     multiprocessing attempt
#     :param proc:
#     :param sref_pose:
#     :param sref_packer:
#     :param smut_pose:
#     :param smut_packer:
#     :param sminimizer:
#     :param s_scorefxn:
#     :return:
#     """
#     if args.testing:  # need to set constant seed, but want different seeds for each process
#         rosetta.init(extra_options="-constant_seed -jran %i" % proc)
#     else:  # need to reinitialize to get random seed, but account for potential of multi processes to start at same time
#         rosetta.init(extra_options="-use_time_as_seed -seed_offset %i" % proc)
#
#     # Make sure that the sequence has actually mutated, and that only 1 AA has changed
#     sref_packer.apply(sref_pose)
#     smut_packer.apply(smut_pose)
#
#     differences = 0
#     for i, b in enumerate(sref_pose.sequence()):
#         if b == smut_pose.sequence()[i]:
#             continue
#         differences += 1
#     assert differences == 1, "%i differences in the before and after sequence. if 0, mutation was not made, if more than 1 multiple bases mutated.\nBefore:\n%s\nAfter:\n%s" % (differences, sref_pose.sequence(), smut_pose.sequence())
#
#     sminimizer(sref_pose)
#     sminimizer(smut_pose)
#
#     return s_scorefxn(sref_pose) - s_scorefxn(smut_pose)


def ddg(reference_pose, mutation_info, pack_radius=args.distance, mr_scorefxn=scorefxn):
    assert len(mutation_info) == 3 and isinstance(mutation_info, list), "mutation info must be a list with values of [postion, reference AA, mutated AA]. unclear input %s" % mutation_info

    mut_pose = rosetta.Pose()
    mut_pose.assign(reference_pose)

    mut_task = rosetta.standard_packer_task(mut_pose)
    ref_task = rosetta.standard_packer_task(reference_pose)

    # the Vector1 of booleans (a specific object) is needed for specifying the mutation, this demonstrates another more direct method of setting PackerTask options for design
    mut_aa_bool = rosetta.utility.vector1_bool()
    ref_aa_bool = rosetta.utility.vector1_bool()

    # PyRosetta uses several ways of tracking amino acids (ResidueTypes) the numbers 1-20 correspond individually to the 20 proteogenic amino acids aa_from_oneletter returns the integer representation of an amino acid from its one letter code
    # convert mutant_aa to its integer representation
    mut_aa = rosetta.aa_from_oneletter_code(mutation_info[2])
    ref_aa = rosetta.aa_from_oneletter_code(mutation_info[1])

    # TODO is this the best way to cause a mutation? Specifically, creating a list of 20F and 1T bools
    # mutation is performed by using a PackerTask with only the mutant amino acid available during design to do this, construct a Vector1 of booleans indicating which amino acid (by its numerical designation, see above) to allow
    for i in range(1, 21):
        # in Python, logical expression are evaluated with priority, thus the line below appends to aa_bool the truth (True or False) of the statement i == mutant_aa
        mut_aa_bool.append(i == mut_aa)
        ref_aa_bool.append(i == ref_aa)

    # modify the mutating residue's assignment in the PackerTask using the Vector1 of booleans across the proteogenic amino acids
    mut_task.nonconst_residue_task(mutation_info[0]).restrict_absent_canonical_aas(mut_aa_bool)
    ref_task.nonconst_residue_task(mutation_info[0]).restrict_absent_canonical_aas(ref_aa_bool)
    # prevent residues from packing by setting the per-residue "options" of the PackerTask

    #center = rosetta.pose.residue( mutant_position ).nbr_atom_xyz()  # throws error
    center = reference_pose.residue(mutation_info[0]).nbr_atom_xyz()

    for i in range(1, reference_pose.total_residue() + 1):
        # only pack the mutating residue and any within the pack_radius
        # print i, "\t", center.distance_squared(test_pose.residue(i).nbr_atom_xyz()), "\t", test_pose.residue(i).name3()  # used for testing/identifying problem with original if statment
        # if not i == mutant_position or center.distance_squared(test_pose.residue(i).nbr_atom_xyz()) > pack_radius**2:  # this is the worst error i've ever had to debug NEVER EVER UNCOMMENT WITHOUT TALKING TO DAN
        if center.distance_squared(reference_pose.residue(i).nbr_atom_xyz()) > pack_radius**2:  # mutated residue will always have a distance of 0,  which is less than the square of any real number
            mut_task.nonconst_residue_task(i).prevent_repacking()
            ref_task.nonconst_residue_task(i).prevent_repacking()
        else:
            if i == mutation_info[0]:
                continue
            else:
                mut_task.nonconst_residue_task(i).restrict_to_repacking()
                ref_task.nonconst_residue_task(i).restrict_to_repacking()
        # else:  # this is the key block that ultimately proved the original if statement was incorrect as it only allowed the targeted residue to be repacked regardless of distance given
        #     print "THIS HAPPENED", i, mutant_position

    # apply the mutation and pack nearby residues
    mut_packer = rosetta.PackRotamersMover(mr_scorefxn, mut_task)
    ref_packer = rosetta.PackRotamersMover(mr_scorefxn, ref_task)


    movemap = rosetta.MoveMap()
    movemap.set_bb(True)  # change backbone confirmation
    movemap.set_chi(True)  # change side chain confirmation
    min_mover = rosetta.MinMover()
    min_mover.movemap(movemap)
    min_mover.score_function(mr_scorefxn)
    min_mover.min_type('dfpmin_armijo_nonmonotone')

    ddg_results = []
    r_pose = rosetta.Pose()
    m_pose = rosetta.Pose()

    while len(ddg_results) < args.comparisons:
        r_pose.assign(reference_pose)
        m_pose.assign(mut_pose)

        ref_packer.apply(r_pose)
        mut_packer.apply(m_pose)

        # Make sure that the sequence has actually mutated, and that only 1 AA has changed
        differences = 0
        for i, b in enumerate(r_pose.sequence()):
            if b == m_pose.sequence()[i]:
                continue
            differences += 1
        assert differences == 1, "%i differences in the before and after sequence. if 0, mutation was not made, if more than 1 multiple bases mutated.\nBefore:\n%s\nAfter:\n%s" % (differences, r_pose.sequence(), m_pose.sequence())

        min_mover.apply(r_pose)
        min_mover.apply(m_pose)

        ddg_results.append(mr_scorefxn(r_pose) - mr_scorefxn(m_pose))

    return ddg_results

#attempt to incoproate multiprocessing
    # p = multiprocessing.Pool(multiprocessing.cpu_count())
    # fx_args = partial(ddg_score, sref_pose=reference_pose, sref_packer=ref_packer, smut_pose=mut_pose, smut_packer=mut_packer, sminimizer=min_mover)
    # print "partial constructed"
    # results = p.map(fx_args, xrange(args.comparisons))
    #
    # return results






if __name__ == '__main__':

    # 1. Verify arg values are as expected or reset based on inputs
    assert not re.search("clean", args.reference), "Please specify non-cleaned reference name, cleaned file will be created if necessary."
    assert not args.threshold > 1, "Please specify threshold value as decimal."
    if args.max == 0:
        args.max = float("inf")

    # 2. Check if output directory exists, and if not create it
    if not os.path.exists(args.output):
        os.makedirs(args.output)
        log("Specified output directory (%s) does not exist. Creating it now." % args.output)
    log("\nScript executed at: " + str(datetime.datetime.now()) + "\nWith the following command line options:")
    log("\t" + "\n\t".join(map(str, [str(arg) + "\t" + str(getattr(args, arg)) for arg in vars(args) if arg is not "score_function"])))  # log all input options except score function
    log("\t" + "\n\t".join(map(str, [str(arg) + "\t" + str(getattr(args, arg)) for arg in vars(args) if arg is "score_function"])))  # log score function, done separate due to size/formating of output

    # 3. Must start from "cleaned" pdb file, check for it and make it if it doesn't exist
    if not os.path.exists(args.reference.replace(".pdb", ".clean.pdb")):
        log("Reference file " + args.reference + " has not been previously cleaned. Cleaning now.")
        toolbox.cleanATOM(args.reference)
    args.reference = args.reference.replace(".pdb", ".clean.pdb")
    assert len([chain for chain in PDBParser().get_structure("ref", args.reference).get_chains()]) == 1, "PDB provided has multiple chains associated with it. Represents a new type of situation that may influence alignment numbers which are already bad. Please consider and consult with others as to the best way to handle this situation. At a minimum, on mutation read in the chain ID is hardcoded to 'A' which likely needs to be changed" + str([chain for chain in PDBParser().get_structure("ref", args.reference).get_chains()])

    # 4. Verify mutation list provided is correctly formatted  # TODO add blast alignment and other considerations for non-identity mutation modeling
    with open(args.mutations, "r") as f:
        assert f.readline().rstrip().lower() == "output_name\tposition\treference\tmutated", "Expected 4 column structure separated by tabs with headings of 'output_name', 'position', 'reference', and 'mutated' as first line"
        mutations = {line.rstrip().split("\t")[0]: line.rstrip().split("\t")[1:] for line in f if not line.startswith("#")}   # Read mutations into dictionary
        # mutations[output_name] = [position, reference, mutated]

    for mutation in mutations:  # verify reference AA in mutation is identical to AA in pdb file, and adjust numbering
        """
        This is a nightmare because the pose structure/sequence after cleaning is not directly reflective of the .pdb sequence.
        *NOTE* this currently assumes single chain ID, and there is an assert statement forcing this to be true.
        If the assert is to be changed, the hard coding of "A" as the chain ID in the first assert may must change as well.
        This is a terrible approximation of alignment.
        """
        verification_pose = rosetta.Pose()
        rosetta.pose_from_pdb(verification_pose, args.reference)  # all .pdb files correspond to same reference, therefor only need to check the first
        assert verification_pose.residue(verification_pose.pdb_info().pdb2pose("A", int(mutations[mutation][0]))).name1() == mutations[mutation][1], "Alignment of expected reference base and actual base failed. This Represents a case that has not been encountered yet and is a TODO. Begin development of how to better align/deal with this situation. Temporarily consider removing the following mutation from %s\n%s" % (args.mutations, mutation)  # TODO improve this. may involve simply skipping specific mutations with notice and reminder on final stats
        mutations[mutation][0] = verification_pose.pdb_info().pdb2pose("A", int(mutations[mutation][0]))  # update coordinates to be relevant to pose location, not reference/pdb locations

    log("\n%i total mutations found and all pose AA matches reference AA.\n\nBeginning reference optimization." % len(mutations))

    # TODO add checkpoint here to abort program if running on TACC headnode ... FIRST VERIFY pyrosetta can be initialized on headnode at all without angering TACC

    # 5. multiprocessing block to minimize energy
    p = multiprocessing.Pool(multiprocessing.cpu_count())
    optimized_scores = p.map(minimize_energy, xrange(args.replicates))
    optimized_scores = sorted(optimized_scores)  # final scores sorted for easier stat determination and ref assignment

    log("\nReference optimization complete with the following statistics:\nTotal Replicates:\t%i\nLowest Score:\t%f\nLowest Score Count:\t%i\nHighest Score:\t%f\nAverage Score:\t%f\nStd Dev:\t%f\nAll Scores:" % (args.replicates, min(optimized_scores), optimized_scores.count(min(optimized_scores)), max(optimized_scores), numpy.mean(optimized_scores), numpy.std(optimized_scores)))

    # 6. Determine which minimized structures to use
    # 6a. Resize optimized_scores based on user input
    if args.number == 0:
        pass  # 0 number is short hand way to request all structures be used
    elif args.number < 1:  # fraction supplied, only include minimized structures within that fraction
        optimized_scores = optimized_scores[:int(len(optimized_scores) * args.number) + 1]  # +1 included because of rounding  # TODO better use of rounding
    else:  # int supplied use that many best values
        optimized_scores = optimized_scores[:int(args.number)]

    # 6b. Determine which optimized structures are to be used for mutation comparison based on their scores
    possible_best = rosetta.Pose()
    optimized_refs_to_mutate = []

    for files in os.listdir(args.output.rstrip("/") + "/"):
        if re.search(args.prefix + "_Optimized_Replicate_", files):  # ignore files in directory that are not the output of the full optimization script
            rosetta.pose_from_pdb(possible_best, files)  # will produce slightly different values if pdb file not written/read back in
            log("\t".join(map(str, [files, scorefxn(possible_best)])))
            if scorefxn(possible_best) in optimized_scores:
                optimized_refs_to_mutate.append(files)
    log("\nThe following %i optimized reference file(s) will be mutated and scored:\n\t" % len(optimized_refs_to_mutate) + "\n\t".join(map(str, optimized_refs_to_mutate)))

    # 7. Mutation and delta delta G stuff
    log("\nMutation and Energy minimization beginning...")

    for mutation in mutations:
        all_ref_deltas = []
        log("%s ddg calculation beginning" % mutation)
        for ref in optimized_refs_to_mutate:

            optimized_pose = rosetta.Pose()
            rosetta.pose_from_pdb(optimized_pose, ref)
            assert optimized_pose.is_fullatom(), "ddg function/mutation can only work with fullatom poses. If this triggers remember this is adapted from # 1 # and might have not been necessary as an assertion. The file that caused the problem was:\n\t%s" % ref

            all_deltas = ddg(optimized_pose, mutations[mutation])
            log("\t".join(map(str, [mutation, ref, numpy.mean(all_deltas), numpy.std(all_deltas)])))
            # log("\n".join(map(str, all_deltas)))

            all_ref_deltas.extend(all_deltas)
        log("\t".join(map(str, [mutation, numpy.mean(all_ref_deltas), numpy.std(all_ref_deltas)])))






# # attempts to determine apparnt stocasticity of mutation/minimization function, is not considered graveyard yet
#         print "\t".join(map(str, ["mut", "uniq values"] + [k for k in xrange(10)]))
#         for mutation in mutations:
#             print mutation
#             print mutations[mutation]
#             x = []
#             y = []
#             for z in xrange(19, 20): # pack radius
#                 x = []
#                 y = []
#                 j = []
#                 poses = []
#                 for n in xrange(10): # true replicates
#                     rosetta.pose_from_pdb(optimized_pose, args.reference)
#                     # print scorefxn(optimized_pose)
#                     mutant_pose = rosetta.Pose()
#                     mutant_pose.assign(optimized_pose)
#                     output_pose = mutate_residue(mutant_pose, mutations[mutation][0], mutations[mutation][2], pack_radius=z)
#                     #comparison_pose = rosetta.Pose()
#                     #comparison_output = toolbox.mutate_residue(mutant_pose, mutations[mutation][0], mutations[mutation][2], pack_radius=z)
#                     #output_pose = toolbox.mutate_residue(mutant_pose, mutations[mutation][0], mutations[mutation][2], pack_radius=z**z)
#                     # print "pose"
#                     # print "output is mutant\t", output_pose == mutant_pose
#                     # print "output is reference\t", output_pose == optimized_pose
#                     # print "mutant is reference\t", mutant_pose == optimized_pose
#                     # print "seq"
#                     # print "output is mutant\t", output_pose.sequence() == mutant_pose.sequence()
#                     # print "output is reference\t", output_pose.sequence() == optimized_pose.sequence()
#                     # print "mutant is reference\t", mutant_pose.sequence() == optimized_pose.sequence()
#                     # print "scores"
#                     # print "reference\t", scorefxn(optimized_pose)
#                     # print "mutant\t", scorefxn(mutant_pose)
#                     # print "output\t", scorefxn(output_pose)
#                     if (output_pose[0] not in x) and (output_pose[1] not in y):
#
#                         x.append(output_pose[0])
#                         y.append(output_pose[1])
#                         j.append(str([output_pose[0], output_pose[1]]))
#                         poses.append(output_pose[2])
#                         if len(x) != len(j):
#                             print output_pose
#                             print x
#                             print y
#                             print j
#                             print poses
#                             assert False, "different length x,j"
#                     elif str([output_pose[0], output_pose[1]]) not in j:
#                         print output_pose
#                         print x
#                         print y
#                         print j
#                         print poses
#                         assert False, "differential output?"
#
#                     if len(x) == 2:
#                         for res_number in range(1, poses[0].total_residue() + 1):
#                             if poses[0].residue(res_number) == poses[1].residue(res_number):
#                                 print "First\n", poses[0].residue(res_number), "\nSecond\n", poses[1].residue(res_number)
#                                 sys.exit()
#
#                     #y.append(scorefxn(comparison_output))
#                     #print scorefxn(output_pose)
#                 print "\t".join(map(str, [mutation, z, len(set(x))] + x))
#                 print "\t".join(map(str, [mutation, z, len(set(y))] + y))
#                 print "\t".join(map(str, [mutation, z, len(set(j))] + j))
#             sys.exit()
# #            rosetta.MutateResidue()
#
#                 # # #x.append(mutate_residue(optimized_pose, mutations[mutation][0], mutations[mutation][2], proc=_))
#                 # # #y.append(mutate_residue(test_pose, mutations[mutation][0], mutations[mutation][2], proc=_))
#                 # # mutant_pose = rosetta.Pose()
#                 # # rosetta.pose_from_pdb(mutant_pose, ref)  # reload reference from pdb file each time
#                 # # task = rosetta.standard_packer_task(mutant_pose)
#                 # # aa_bool = rosetta.utility.vector1_bool()
#                 # # mutant_aa = rosetta.aa_from_oneletter_code(mutations[mutation][2])
#                 # # for i in range(1, 21):
#                 # #     aa_bool.append(i == mutant_aa)
#                 # # task.nonconst_residue_task(mutations[mutation][0]).restrict_absent_canonical_aas(aa_bool)
#                 # # center = mutant_pose.residue(mutations[mutation][0]).nbr_atom_xyz()
#                 # # for i in range(1, mutant_pose.total_residue() + 1):
#                 # #     if not i == mutations[mutation][0] or center.distance_squared(mutant_pose.residue(i).nbr_atom_xyz()) > 10**2:
#                 # #         task.nonconst_residue_task(i).prevent_repacking()
#                 # # packer = rosetta.PackRotamersMover(scorefxn, task)
#                 # # packer.apply(mutant_pose)
#                 # # x.append(scorefxn(mutant_pose))
#                 # # assert not mutant_pose.sequence() == optimized_pose.sequence()
#                 # rosetta.init(extra_options="-use_time_as_seed")
#                 # #rosetta.pose_from_pdb(optimized_pose, ref)
#                 # #mutant_pose = toolbox.mutate_residue(optimized_pose, mutations[mutation][0], mutations[mutation][2], 10, scorefxn)
#                 # #assert not mutant_pose.sequence() == optimized_pose.sequence()
#                 #
#                 # movemap = rosetta.MoveMap()
#                 # movemap.set_bb(True)  # change backbone confirmation
#                 # movemap.set_chi(True)  # change side chain confirmation
#                 # min_mover = rosetta.MinMover()
#                 # min_mover.movemap(movemap)
#                 # min_mover.score_function(scorefxn)
#                 # min_mover.min_type('dfpmin_armijo_nonmonotone')
#                 #
#                 # test_pose = rosetta.Pose()
#                 # rosetta.pose_from_pdb(test_pose, args.reference)
#                 #
#                 # for i in xrange(10):
#                 #     min_mover.apply(test_pose)
#                 #     x.append(scorefxn(test_pose))
#                 # #    min_mover.apply(mutant_pose)
#                 # #    x.append(scorefxn(mutant_pose))
#
#
#
#
#             #rosetta.pose_from_pdb(optimized_pose, ref)
#             #print mutation, scorefxn(optimized_pose), len(set(x)), x, len(set(y)), y
#
#
#
#             # x = []
#             # print mutation
#             # for _ in xrange(5):
#             #     mutant_pose = rosetta.Pose()
#             #     mutant_pose.assign(optimized_pose)
#             #
#             #     task = rosetta.standard_packer_task(mutant_pose)  # is this the cause of determinizim?
#             #     aa_bool = rosetta.utility.vector1_bool()
#             #     mutant_aa = rosetta.aa_from_oneletter_code(mutations[mutation][2])
#             #
#             #     for i in range(1, 21):
#             #         aa_bool.append(i == mutant_aa)
#             #
#             #     task.nonconst_residue_task(mutations[mutation][0]).restrict_absent_canonical_aas(aa_bool)
#             #
#             #     center = mutant_pose.residue(mutations[mutation][0]).nbr_atom_xyz()
#             #
#             #     for i in range(1, mutant_pose.total_residue()  + 1):
#             #         if not i == mutations[mutation][0] or center.distance_squared(mutant_pose.residue(i).nbr_atom_xyz()) > 10**2:
#             #             task.nonconst_residue_task(i).prevent_repacking()
#             #
#             #     packer = rosetta.PackRotamersMover(scorefxn, task)
#             #     packer.apply(mutant_pose)
#             #     x.append(scorefxn(mutant_pose))
#             # print x
#
#             # mutant_pose = rosetta.Pose()
#             # mutant_pose.assign(optimized_pose)
#             #
#             # task = rosetta.standard_packer_task(mutant_pose)  # is this the cause of determinizim?
#             # aa_bool = rosetta.utility.vector1_bool()
#             # mutant_aa = rosetta.aa_from_oneletter_code(mutations[mutation][2])
#             #
#             # for i in range(1, 21):
#             #     aa_bool.append(i == mutant_aa)
#             #
#             # task.nonconst_residue_task(mutations[mutation][0]).restrict_absent_canonical_aas(aa_bool)
#             #
#             # center = mutant_pose.residue(mutations[mutation][0]).nbr_atom_xyz()
#             #
#             # for i in range(1, mutant_pose.total_residue()  + 1):
#             #     if not i == mutations[mutation][0] or center.distance_squared(mutant_pose.residue(i).nbr_atom_xyz()) > 10**2:
#             #         task.nonconst_residue_task(i).prevent_repacking()
#             #
#             # packer = rosetta.PackRotamersMover(scorefxn, task)
#             # # x = []
#             # # test_pose = rosetta.Pose()
#             # # for _ in xrange(1):
#             # #     test_pose.assign(mutant_pose)
#             # #     # if args.testing:  # need to set constant seed, but want different seeds for each process
#             # #     #     rosetta.init(extra_options="-constant_seed -jran %i" % _)
#             # #     # else:  # need to reinitialize to get random seed, but account for potential of multi processes to start at same time
#             # #     #     rosetta.init(extra_options="-use_time_as_seed -seed_offset %i" % _)
#             # #
#             # #     packer.apply(test_pose)
#             # #     x.append(scorefxn(test_pose))
#             # #     print x
#             # # print len(set(x))
#             #
#             # packer.apply(mutant_pose)  # results in single movement effect with distance 0, therefore does not need to be replicated. Tested both with and without rosetta reinitialization
#
#
#             # test_pose = rosetta.Pose()
#             # movemap = rosetta.MoveMap()
#             # movemap.set_bb(True)  # change backbone confirmation
#             # movemap.set_chi(True)  # change side chain confirmation
#             # min_mover = rosetta.MinMover()
#             # min_mover.movemap(movemap)
#             # min_mover.score_function(scorefxn)
#             # min_mover.min_type('dfpmin_armijo_nonmonotone')
#
#             # for _ in xrange(10):
#             #     test_pose.assign(mutant_pose)
#             #     if args.testing:  # need to set constant seed, but want different seeds for each process
#             #         rosetta.init(extra_options="-constant_seed -jran %i" % _)
#             #     else:  # need to reinitialize to get random seed, but account for potential of multi processes to start at same time
#             #         rosetta.init(extra_options="-use_time_as_seed -seed_offset %i" % _)
#             #     min_mover.apply(test_pose)
#             #     x.append(scorefxn(test_pose))
#             #     print x
#             # print len(set(x))
#             #min_mover.apply(mutant_pose)
#             #print scorefxn(optimized_pose)
#             #print "\t".join(map(str, [mutation, scorefxn(optimized_pose), scorefxn(mutant_pose)]))
#







# Graveyard
# ######################################## 2 ########################################
# def sample_refinement(pdb_filename, kT=1.0, smallmoves=3, shearmoves=5, backbone_angle_max=7, cycles=args.moves, jobs=args.replicates, job_output=args.output.rstrip("/") + "/" + args.prefix):
#     """
#     Performs fullatom structural refinement on the input  <pdb_filename>  by
#         perturbing backbone torsion angles with a maximum perturbation of
#         <backbone_angle_max>  for  <cycles>  trials of
#         <smallmoves>  perturbations of a random residue's phi or psi and
#         <shearmoves>  perturbations of a random residue's phi and the preceding
#         residue's psi followed by gradient based backbone torsion angle
#         minimization and sidechain packing with an acceptance criteria scaled
#         by  <kT>.  <jobs>  trajectories are performed, continually exporting
#         structures to a PyMOL instance.
#         Output structures are named  <job_output>_(job#).pdb.
#     """
#
#     ######################################## 2 ########################################
# # sample_refinement function taken from D070_Refinement.py sample script with the following additional information:
# # A GENERAL EXPLANATION
# #
# # """
# # refinement.py
# #
# # This script performs simple high-resolution (fullatom) refinement on an input
# # pose. The local conformation space is sampled using small backbone torsion angle
# # perturbations followed by backbone torsion angle minimization and sidechain
# # packing. The Rosetta standard score function evaluations are used to accept or
# # reject proposed structures. Some increases in score are accepted to escape
# # local minima.
# #
# # Instructions:
# #
# # 1) ensure that your PDB file is in the current directory
# # 2) run the script:
# #     from commandline                        >python D070_Refinement.py
# #
# #     from within python/ipython              [1]: run D070_Refinement.py
# #
# # Author: Evan H. Baugh
# #     revised and motivated by Robert Schleif
# #
# # Updated by Boon Uranukul, 6/9/12
# # Simplified special constant seed initialization ~ Labonte
# #
# # References:
# #     P. Bradley, K. Misura, and D. Baker, "Toward high-resolution de novo
# #         structure prediction for small proteins," Science 309 (5742)
# #         1868-1871 (2005).
# #
# # """
# #
# # ################################################################################
# # # THE BASIC PROTOCOL, sample_refinement
# #
# # """
# # This sample script is setup for usage with
# #     commandline arguments,
# #     default running within a python interpreter,
# #     or for import within a python interpreter,
# #         (exposing the sample_refinement method)
# #
# # The method sample_refinement:
# # 1.  creates a pose from the desired PDB file
# # 2.  creates a (fullatom) reference copy of the pose
# # 3.  creates a standard fullatom ScoreFunction
# # 4.  creates a MoveMap with all backbone torsion angles free
# # 5.  sets up a SmallMover for small backbone torsion angle perturbations
# # 6.  sets up a ShearMover for small backbone torsion angle perturbations
# # 7.  sets up a MinMover for backbone torsion minimization
# # 8.  sets up a PackRotamersMover for sidechain packing
# # 9.  create a PyMOL_Mover for viewing intermediate output
# # 10. export the original structure, and scores, to PyMOL
# # 11. sets up a RepeatMover on a TrialMover of a SequenceMover
# #         -setup the TrialMover
# #             a.  create a SequenceMover with the:
# #                     >SmallMover
# #                     >ShearMover
# #                     >MinMover
# #                     >PackRotamersMover
# #                     >PyMOL_Mover
# #             b.  create a MonteCarlo object for assessing moves
# #             c.  create the TrialMover (on the SequenceMover)
# #         -create the RepeatMover (on the TrialMover)
# # 12.  creates a (Py)JobDistributor for managing multiple trajectories
# # 13.  stores the original score evaluation
# # 14.  performs the refinement protocol, for each trajectory:
# #          a. set necessary variables for the new trajectory
# #              -reload the starting pose
# #              -change the pose's PDBInfo.name, for the PyMOL_Mover
# #              -reset the MonteCarlo object
# #          b. perform the sampling and assessment using the RepeatMover
# #          c. output the (lowest scoring) decoy structure
# #              -output the decoy structure using the PyJobDistributor
# #              -export the decoy structure to PyMOL
# #              -store the decoy score
# # 15.  outputs the score evaluations
# #
# # """
# ######################################## 2 ########################################
#
#     # Steps 1-3 are typical protocol setup. Altnatively the following line could be used (if steps 1-3 commented out, but is not supported)
#     #refinement = ClassicRelax( scorefxn )
#
#
#     #  1. create a pose from the desired PDB file
#     pose = rosetta.Pose()
#     rosetta.pose_from_pdb(pose, pdb_filename)
#
#     # 2. create a reference copy of the pose in fullatom
#     starting_pose = rosetta.Pose()
#     starting_pose.assign(pose)
#
#     # 3. create a standard ScoreFunction
#     scorefxn = rosetta.get_fa_scorefxn()  # implement the desired ScoreFunction here  # TODO make function take argparse input of what ScoreFunction to use
#
#
#     #### Setup custom high-resolution refinement protocol
#     #### backbone refinement protocol
#
#     # 4. create a MoveMap, all backbone torsions free
#     movemap = rosetta.MoveMap()
#     movemap.set_bb(True)
#
#     # 5. create a SmallMover
#     """a SmallMover perturbs a random (free in the MoveMap) residue's phi or psi torsion angle for an input number of
#     times and accepts of rejects this change based on the Metropolis Criteria using the "rama" ScoreType and the
#     parameter kT"""
#     smallmover = rosetta.SmallMover(movemap, kT, smallmoves)  # set the maximum angle to backbone_angle_max, apply it smallmoves times
#
#     smallmover.angle_max(backbone_angle_max)  # angle_max is secondary structure dependent, however secondary structure has not been evaluated in this protocol, thus they are all set to the same value0 at once
#
#     # 5. ALTERNATIVE  # @DED I believe the 3 lines of code are intended as an alternative to 5
#     """use the overloaded version of the SmallMover.angle_max method if you want to use secondary structure biased moves"""
#     # smallmover.angle_max('H', backbone_angle_max)
#     # smallmover.angle_max('E', backbone_angle_max)
#     # smallmover.angle_max('L', backbone_angle_max)
#
#     # 6. create a ShearMover
#     """a ShearMover is identical to a SmallMover except that the angles perturbed are instead a random (free in the
#     MoveMap) residue's phi and the preceding residue's psi, this reduces the downstream structural change"""
#     shearmover = rosetta.ShearMover(movemap, kT, shearmoves)  # set the maximum angle to backbone_angle_max, apply it shearmoves times
#     shearmover.angle_max(backbone_angle_max)  # same angle_max restictions as SmallMover
#
#     # 6. ALTERNATIVE  # @DED like in 5, I believe the 3 lines of code are intended as an alternative to 6
#     """use the overloaded version of the SmallMover.angle_max method if you want to use secondary structure biased moves"""
#     # shearmover.angle_max('H', backbone_angle_max)
#     # shearmover.angle_max('E', backbone_angle_max)
#     # shearmover.angle_max('L', backbone_angle_max)
#
#     # 7. create a MinMover, for backbone torsion minimization
#     minmover = rosetta.MinMover()
#     minmover.movemap(movemap)
#     minmover.score_function(scorefxn)
#
#     #### sidechain refinement protocol, simple packing
#
#     # 8. setup a PackRotamersMover
#     to_pack = rosetta.standard_packer_task(starting_pose)
#     to_pack.restrict_to_repacking()    # prevents design, packing only
#     to_pack.or_include_current(True)    # considers the original sidechains
#     packmover = rosetta.PackRotamersMover(scorefxn, to_pack)
#
#     #### assess the new structure
#     # 9. create a PyMOL_Mover
#     if args.pymol:
#         pymover = rosetta.PyMOL_Mover()
#         # pymover.keep_history(True)  # load structures into sucesseive states
#         """the PyMOL_Mover slows down the protocol SIGNIFICANTLY but provides very informative displays
#         the keep_history flag (when True) tells the PyMOL_Mover to store new structures into successive states, for a single
#         trajectory, this allows you to see intermediate changes (depending on where the PyMOL_Mover is applied), when using
#         a JobDistributor or otherwise displaying multiple trajectories with a single protocol, the output can get confusing
#         to interpret, by changing the pose's PDBInfo.name the structure will load into a new PyMOL state"""
#         # pymover.update_energy(True)  # see the total score in color
#
#         # 10. export the original structure, and scores, to PyMOL
#         pymover.apply(pose)
#         scorefxn(pose)
#         pymover.send_energy(pose)
#
#     # 11. setup a RepeatMover on a TrialMover of a SequenceMover (wow!)
#     # -setup a TrialMover
#     #    a. create a SequenceMover of the previous moves
#     #### add any other moves you desire
#     combined_mover = rosetta.SequenceMover()
#     combined_mover.add_mover(smallmover)
#     combined_mover.add_mover(shearmover)
#     combined_mover.add_mover(minmover)
#     combined_mover.add_mover(packmover)
#     if args.pymol:
#         #### explore the protocol using the PyMOL_Mover, try viewing structures before they are accepted or rejected
#         combined_mover.add_mover(pymover)
#     #    b. create a MonteCarlo object to define success/failure
#     mc = rosetta.MonteCarlo(pose, scorefxn, kT)  # must reset for each trajectory!
#     # c. create the TrialMover
#     trial = rosetta.TrialMover(combined_mover, mc)
#
#     #### explore the protocol using the PyMOL_Mover, try viewing structures after acceptance/rejection, comment-out the lines below
#     # original_trial = TrialMover(combined_mover, mc)
#     # trial = SequenceMover()
#     # trial.add_mover(original_trial)
#     # trial.add_mover(pymover)
#
#     #### for each trajectory, try cycles number of applications
#
#     # -create the RepeatMover
#     refinement = rosetta.RepeatMover(trial, cycles)
#
#     # 12. create a (Py)JobDistributor
#     jd = rosetta.PyJobDistributor(job_output, jobs, scorefxn)
#     jd.native_pose = starting_pose
#
#     # 13. store the score evaluations for output
#     """printing the scores as they are produced would be difficult to read, Rosetta produces a lot of output when running"""
#     scores = [0] * (jobs + 1)
#     scores[0] = scorefxn(starting_pose)
#
#     # 14. perform the refinement protocol
#     counter = 0    # for exporting to PyMOL
#     while not jd.job_complete:
#         # a. set necessary variables for the new trajectory
#         # -reload the starting pose
#         pose.assign(starting_pose)
#         # -change the pose's PDBInfo.name, for the PyMOL_Mover
#         counter += 1
#         log("\n\nStarting replicate: %s" % counter)
#         pose.pdb_info().name(job_output + '_' + str(counter))
#         # -reset the MonteCarlo object (sets lowest_score to that of p)
#         mc.reset(pose)  # if you create a custom protocol, you may have additional variables to reset, such as kT
#
#         """if you create a custom protocol, this section will most likely change, many protocols exist as single Movers
#         or can be chained together in a sequence (see above) so you need only apply the final Mover"""
#         # b. apply the refinement protocol
#         refinement.apply(pose)
#         ####
#
#         # c. output the lowest scoring decoy structure for this trajectory
#         # -recover and output the decoy structure to a PDB file
#         mc.recover_low(pose)
#         jd.output_decoy(pose)
#         if args.pymol:  # -export the final structure to PyMOL for each trajectory
#             pose.pdb_info().name(job_output + '_' + str(counter) + '_final')
#             pymover.apply(pose)
#             pymover.send_energy(pose)    # see the total score in color
#
#         # -store the final score for this trajectory
#         scores[counter] = scorefxn(pose)
#         log(str(str(counter) + "\t" + str(scores[counter])))
#
#     # 15. output the score evaluations
#     print 'Original Score\t:\t', scores[0]
#     for i in range(1, len(scores)):    # print out the job scores
#         print job_output + '_' + str(i) + '\t:\t', scores[i]
#
#     return scores  # for other protocols
# ######################################## 2 ########################################

# # Block intending to use sample refinement
# log(str("\nEnergy Minimization to begin with %s replicates, with each job taking at most %s moves, requiring an improvement of %s percent over the course of %s moves") % (args.replicates, args.moves, str(args.threshold * 100), str(.01 * args.moves)))
# ref_e_values = sample_refinement(args.reference)
# log("\nEnergy values:")
# log("\n".join(map(str, ref_e_values)))


# def minimize_energy(pdb_filename, failures=args.failures, max_attempts=args.max):
#     """
#     :param pdb_filename: pdb file to be minimized
#     :param failures: how many unsucessful attempts
#     :param max_attempts: cap on how many times structure should be attempted to improve
#     function description goes here
#
#     Rotamer minimization is slow ... minmover is fast. how necessary is rotamer mover at all, how necessary is it for loop?
#
#     """
#     # 1. create a pose from the desired PDB file, and copies of it for the eventual output, and optimization
#     pose = rosetta.Pose()
#     loop_pose = rosetta.Pose()
#     output_pose = rosetta.Pose()
#     rosetta.pose_from_pdb(pose, pdb_filename)
#     loop_pose.assign(pose)
#     output_pose.assign(pose)
#
#     # 2. create a standard ScoreFunction and determine starting score
#     scorefxn = rosetta.get_fa_scorefxn()  # implement the desired ScoreFunction here  # TODO make function take argparse input of what ScoreFunction to use
#     initial_score = scorefxn(pose)
#     scores = [initial_score]  # need to store staring point to determine if improvements are made
#     all_scores = [initial_score]
#
#     # 3. optimize side chains only.
#     log("Repacking Side Chains until %i attempts fail to yield improved structure, " % failures)  #
#     pack_counter = 0  # for book keeping
#     threshold_check = True
#     while threshold_check and (pack_counter < max_attempts) and (len(scores) < failures + 1):  # initial testing suggests that last condition will never trigger at this stage, only hit an asymptote
#         pack_counter += 1  # for bookkeeping
#
#         # 3a. pack and minimize initial pose to remove clashes
#         task = rosetta.standard_packer_task(loop_pose)
#         task.restrict_to_repacking()
#         task.or_include_current(True)
#         pack_rotamers_mover = rosetta.RotamerTrialsMinMover(scorefxn, task)
#         pack_rotamers_mover.apply(loop_pose)
#
#         # 3b. score output, and determine if it is an improvement
#         current_score = scorefxn(loop_pose)  # store value so only calculating once
#         all_scores.append(current_score)
#         if current_score < min(scores):  # TODO verify that scores[0] would be better because it must be the minimum value or it would be reset
#             if args.testing:
#                 log("Improved side chain conformation with %s energy found on the %i improvement attempted after %i failed attempts with following energy: %s" % (current_score, pack_counter, len(scores) - 1, "\n".join(map(str, scores))))
#             output_pose.assign(loop_pose)
#             scores = [current_score]
#         else:  # current loop_pose confirmation is worse, therefore store score, but reset for next loop
#             scores.append(current_score)
#             loop_pose.assign(output_pose)
#
#         # determine if threshold still valid
#         try:  # try less expensive than if to check range assuming usually true. differential controls number of attempts, but if threshold reached quickly therefore only triggered limited number of times
#             energy_differential = sorted(all_scores)[failures + 1] - sorted(all_scores)[0]  # energy improvement over args.failures attempts
#             energy_percent_improvement = float(energy_differential) / (sorted(all_scores)[-2] - sorted(all_scores)[0])  # percent improvement over last failures attempts compared to overall improvement excluding first change (which is always the biggest)  TODO determine if first improvement is always bigger than the rest, or add check to verify
#             if energy_percent_improvement < args.threshold:
#                 threshold_check = False
#         except IndexError:  #
#             pass
#
#     else:  # TODO better control over when while loop exits and to what effect ... also move this into min_map
#         if not threshold_check:  # rotamer stopped because of threshold value
#             log("Repacking Side Chains complete after %i total attempts because improvement over the last %i attempts has dropped to %f which is less than required value of %f\n\tCurrent score is: %f" % (pack_counter, failures, energy_percent_improvement, args.threshold, scores[0]))
#         if not (len(scores) < failures + 1):  # previous attempts have not improved the score at all (DO NOT make this an elif, both can be true)
#             log("Repacking Side Chains complete after %i total attempts. Last %i attempts had following scores: %s\n\tNone better than current score %f" % (pack_counter, failures - 1, scores[1:], scores[0]))
#         if not pack_counter < max_attempts:
#             log("Repacking Side Chains complete after %i total attempts, the maximum allowed based on user input. Last attempts had following scores: %s\n\tCurrent score is: %f" % (pack_counter, all_scores[-failures:], scores[0]))
#         rotamer_score = scores[0]  # first score is the best
#         assert rotamer_score == scorefxn(output_pose), "Score of rotamer optimization changed somehow."
#         assert len(all_scores) == pack_counter + 1, "Score lost somewhere? %i \t %i" % (len(all_scores), pack_counter)  # +1 to pack_counter required because initial score added before new scores calculated
#
#     if args.intermediates:  # block may be eliminated as this is being used to test importance of rotamer optimization
#         for _ in xrange(args.replicates):
#             possible_name = args.output.rstrip("/") + "/" + args.prefix + "_Rotamer_replicate_" + str(_) + ".pdb"
#             if not os.path.isfile(possible_name):
#                 output_pose.dump_pdb(possible_name)
#
#
#     # 4a. create mover and initialize values
#     log("Minimizing overall energy of backbone and side chains until %i attempts fail to yield improved structure" % failures)
#     movemap = rosetta.MoveMap()
#     movemap.set_bb(True)  # change backbone confirmation
#     movemap.set_chi(True)  # change side chain confirmation
#     min_mover = rosetta.MinMover()
#     min_mover.movemap(movemap)
#     min_mover.score_function(scorefxn)
#     min_mover.min_type('dfpmin_armijo_nonmonotone')
#
#     scores = [rotamer_score]  # reset score for book keeping
#     minimization_counter = 0  # for book keeping
#
#     # 4b. optimize backbone and rotamers
#     threshold_check = True
#     while threshold_check and (minimization_counter < max_attempts) and (len(scores) < failures + 1):
#         minimization_counter += 1  # for book keeping
#
#         min_mover.apply(loop_pose)
#         current_score = scorefxn(loop_pose)
#         all_scores.append(current_score)
#
#         if current_score < min(scores):
#             if args.testing:
#                 log("Lower energy conformation with %s energy found on the %i improvement attempted after %i failed attempts with following energy: %s" % (current_score, minimization_counter, len(scores) - 1, "\n".join(map(str, scores))))
#             output_pose.assign(loop_pose)
#             scores = [current_score]
#         else:  # current loop_pose confirmation is worse, therefore store score, but reset to best for next loop
#             scores.append(current_score)
#             loop_pose.assign(output_pose)
#
#         # determine if threshold still valid
#         try:  # try less expensive than if to check range assuming usually true. differential controls number of attempts, but if threshold reached quickly only triggered limited times anyway
#             energy_differential = sorted(all_scores[pack_counter + 1:])[failures + 1] - sorted(all_scores[pack_counter + 1:])[0]  # energy improvement over args.failures attempts
#             energy_percent_improvement = float(energy_differential) / (sorted(all_scores[pack_counter + 1:])[-2] - sorted(all_scores[pack_counter + 1:])[0])  # percent improvement over last failures attempts compared to overall improvement excluding first change (which is always the biggest)  TODO determine if first improvement is always bigger than the rest, or add check to verify
#             if energy_percent_improvement < args.threshold:
#                 threshold_check = False
#         except IndexError:
#             pass
#
#     else:
#
#         if not threshold_check:  # minimization stopped because of threshold value
#             log("Energy minimization complete after %i total attempts because improvement over the last %i attempts has dropped to %f which is less than required value of %f\n\tCurrent score is: %f" % (minimization_counter, failures, energy_percent_improvement, args.threshold, scores[0]))
#         if not (len(scores) < failures + 1):  # previous attempts have not improved the score at all (DO NOT make this an elif, both can be true)
#             log("Energy minimization complete after %i total attempts. Last %i attempts had following scores: %s\n\tNone better than current score %f" % (minimization_counter, failures - 1, scores[1:], scores[0]))
#         if not pack_counter < max_attempts:
#             log("Energy minimization complete after %i total attempts, the maximum allowed based on user input. Last attempts had following scores: %s\n\tCurrent score is: %f" % (minimization_counter, all_scores[-failures:], scores[0]))
#
#
#         log("Energy minimization complete after %i total attempts, and %i failed attempts to improve structure further. Current score is: %i" % (minimization_counter, failures, scores[0]))
#         assert scores[0] == scorefxn(output_pose), "Score of energy minimization changed somehow"
#
#     # 5. finish
#     scores = [initial_score, rotamer_score, scores[0]]
#     log("\nOptimization complete. Initial Score:\t%s\nRotamer Score:\t%s\nFinal Score:\t%s" % (scores[0], scores[1], scores[2]))
#     if args.verbose:
#         log("\nEnergy trajectory:\n" + "\n".join(map(str, all_scores)))
#     return output_pose



    # block used for testing rotamer optimization values within minimize_energy function
    # threshold_check = True
    # while threshold_check and (pack_counter < max_attempts) and (len(scores) < failures + 1):  # initial testing suggests that last condition will never trigger at this stage, only hit an asymptote
    #     pack_counter += 1  # for bookkeeping
    #
    #     # 3a. pack and minimize initial pose to remove clashes
    #     task = rosetta.standard_packer_task(loop_pose)
    #     task.restrict_to_repacking()
    #     task.or_include_current(True)
    #     pack_rotamers_mover = rosetta.RotamerTrialsMinMover(scorefxn, task)
    #     pack_rotamers_mover.apply(loop_pose)
    #
    #     # 3b. score output, and determine if it is an improvement
    #     current_score = scorefxn(loop_pose)  # store value so only calculating once
    #     all_scores.append(current_score)
    #     if current_score < min(scores):  # TODO verify that scores[0] would be better because it must be the minimum value or it would be reset
    #         if args.testing:
    #             log("Improved side chain conformation with %s energy found on the %i improvement attempted after %i failed attempts with following energy: %s" % (current_score, pack_counter, len(scores) - 1, "\n".join(map(str, scores))))
    #         output_pose.assign(loop_pose)
    #         scores = [current_score]
    #     else:  # current loop_pose confirmation is worse, therefore store score, but reset for next loop
    #         scores.append(current_score)
    #         loop_pose.assign(output_pose)
    #
    #     # determine if threshold still valid
    #     try:  # try less expensive than if to check range assuming usually true. differential controls number of attempts, but if threshold reached quickly therefore only triggered limited number of times
    #         energy_differential = sorted(all_scores)[failures + 1] - sorted(all_scores)[0]  # energy improvement over args.failures attempts
    #         energy_percent_improvement = float(energy_differential) / (sorted(all_scores)[-2] - sorted(all_scores)[0])  # percent improvement over last failures attempts compared to overall improvement excluding first change (which is always the biggest)  TODO determine if first improvement is always bigger than the rest, or add check to verify
    #         if energy_percent_improvement < args.threshold:
    #             threshold_check = False
    #     except IndexError:
    #         pass
    #
    # else:  # TODO better control over when while loop exits and to what effect ... also move this into min_map
    #     if not threshold_check:  # rotamer stopped because of threshold value
    #         log("Repacking Side Chains complete after %i total attempts because improvement over the last %i attempts has dropped to %f which is less than required value of %f\n\tCurrent score is: %f" % (pack_counter, failures, energy_percent_improvement, args.threshold, scores[0]))
    #     if not (len(scores) < failures + 1):  # previous attempts have not improved the score at all (DO NOT make this an elif, both can be true)
    #         log("Repacking Side Chains complete after %i total attempts. Last %i attempts had following scores: %s\n\tNone better than current score %f" % (pack_counter, failures - 1, scores[1:], scores[0]))
    #     if not pack_counter < max_attempts:
    #         log("Repacking Side Chains complete after %i total attempts, the maximum allowed based on user input. Last attempts had following scores: %s\n\tCurrent score is: %f" % (pack_counter, all_scores[-failures:], scores[0]))
    #     rotamer_score = scores[0]  # first score is the best
    #     assert rotamer_score == scorefxn(output_pose), "Score of rotamer optimization changed somehow."
    #     assert len(all_scores) == pack_counter + 1, "Score lost somewhere? %i \t %i" % (len(all_scores), pack_counter)  # +1 to pack_counter required because initial score added before new scores calculated
    #

#replaced with multiprocessing
# for _ in xrange(args.replicates):
#     optimized_ref_file_name = args.output.rstrip("/") + "/" + args.prefix + "_replicate_" + str(_) + ".pdb"
#     optimized_ref = rosetta.Pose()
#
#     if not os.path.exists(optimized_ref_file_name):\t
#         optimized_ref.assign(minimize_energy())#args.reference))
#         optimized_ref.dump_pdb(optimized_ref_file_name)
#     else:
#         rosetta.pose_from_pdb(optimized_ref, optimized_ref_file_name)
#
#     optimized_score = scorefxn(optimized_ref)
#     if optimized_score < best_score:
#         best_reference.assign(optimized_ref)
#         best_score = optimized_score

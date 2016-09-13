# Compiling Rosetta and PyRosetta

## Obtaining the Rosetta Source Code

Rosetta and PyRosetta are both compiled from the Rosetta source code; in order to obtain this, you’ll need to acquire the appropriate licenses for BOTH Rosetta and PyRosetta.  Both licenses can be obtained from the Rosetta licensing page:

https://www.rosettacommons.org/software/license-and-download

Once you have your license code, follow the link to the download page and enter the login information from the email.  Click on “Latest Numbered Release” and then “Rosetta 3.6 Source”.  Source code is distributed as a single large gzipped file and is ~ 2.0Gb unzipped.  Haven’t tested this procedure on other versions of Rosetta.

Once you have the zipped source code, upload it to your $WORK/src directory on stampede.  NOTE: compiling on Lonestar 5 is possible, but requires installing an older version of the gcc compiler (GCC 4.4.7 or earlier); this is the default version on stampede.

Log in to Stampede and decompress the file (best to make a new directory to hold everything):

`$ cd $WORK/src`

`$ mkdir rosetta_local_compile`

`$ mv rosetta_src_3.6_bundle.tar.gz rosetta_local_compile`

`$ cd rosetta_local_compile`

`$ tar -xvzf rosetta_src_3.6_bundle.tar.gz`

## Compiling PyRosetta on TACC Stampede

Download these scripts from the BarrickLab github:

`barricklab/mutant-protein-stability/edit/master/PyRosetta_stampede_compile_scripts/DeployPyRosetta.tacc.stampede.py`
A modified version of the DeployPyRosetta.py script that ships with rosetta, for use on TACC Stampede

`barricklab/mutant-protein-stability/edit/master/PyRosetta_stampede_compile_scripts/run_DeployPyRosetta.tacc.stampede.sh`
A batch submission script for the above

Locate the PyRosetta.develop directory in the source file, should be located here:

`$ cd PATH/TO/ROSETTA/BUNDLE/tools/PyRosetta.develop`

Copy the scripts to the PyRosetta.develop directory

Submit the batch script using sbatch
*NOTE* This version of the batch script uses the largemem queue because the standard Stampede compute nodes run out of memory when compiling with more than ~2 cores

 `sbatch run_DeployPyRosetta.tacc.stampede.sh`
 
 This will take < 2h to complete once the job starts (largemem queue can get pretty long...)



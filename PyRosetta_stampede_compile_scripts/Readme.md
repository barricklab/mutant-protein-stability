# Compiling Rosetta and PyRosetta

## Obtaining the Rosetta Source Code

Rosetta and PyRosetta are both compiled from the Rosetta source code; in order to obtain this, you’ll need to acquire the appropriate licenses for BOTH Rosetta and PyRosetta.  Both licenses can be obtained from the Rosetta licensing page:

https://www.rosettacommons.org/software/license-and-download

Once you have your license code, follow the link to the download page and enter the login information from the email.  Click on “Latest Numbered Release” and then “Rosetta 3.6 Source”.  Source code is distributed as a single large gzipped file and is ~ 2.0Gb unzipped.  Haven’t tested this procedure on other versions of Rosetta.

Once you have the zipped source code, upload it to your $WORK/src directory on stampede.  NOTE: compiling on Lonestar 5 is possible, but requires installing an older version of the gcc compiler (GCC 4.4.7 or earlier); this is the default version on stampede.

Log in to Stampede and decompress the file (best to make a new directory to hold everything):

```
$ cd $WORK/src
$ mkdir rosetta_local_compile
$ mv rosetta_src_3.6_bundle.tar.gz rosetta_local_compile
$ cd rosetta_local_compile
$ tar -xvzf rosetta_src_3.6_bundle.tar.gz
```

## Compiling PyRosetta on TACC Stampede

Download these scripts from the BarrickLab github:

`barricklab/mutant-protein-stability/edit/master/PyRosetta_stampede_compile_scripts/DeployPyRosetta.tacc.stampede.py`
A modified version of the DeployPyRosetta.py script that ships with rosetta, for use on TACC Stampede

`barricklab/mutant-protein-stability/edit/master/PyRosetta_stampede_compile_scripts/run_DeployPyRosetta.tacc.stampede.sh`
A batch submission script for the above

Locate the PyRosetta.develop directory in the source file, should be located here:

`$ cd PATH/TO/ROSETTA/BUNDLE/tools/PyRosetta.develop`

Copy the DeployPyRosetta scripts to the PyRosetta.develop directory

*NOTE* This version of the batch script uses the largemem queue because the standard Stampede compute nodes run out of memory when compiling with more than ~2 cores

If you don't want to use `largemem`, change `#BATCH -p largemem` to `#BATCH -p normal`, the `#BATCH -t 2:0:0` to `BATCH -t 6:0:0` and the `--jobs=32` option in the script call to `--jobs=2`.  This should compile correctly, but will take more time.

Submit the batch script (`run_DeployPyRosetta.tacc.stampede.sh`) using sbatch:

`$ sbatch run_DeployPyRosetta.tacc.stampede.sh` 
 
This will take < 2h to complete once the job starts (largemem queue can get pretty long...)

Once the job is completed, you can check the `PyRosComp.[JOB_ID].out` file to make sure everything worked. 

`$ tail -n 2 PyRosComp.7609924.out`

If the compilation completed normally, the last lines should be:

```
Compiling py files...
Done!
```

All of the PyRosetta Libraries should now be compiled.

## Setting up PyRosetta on TACC Stampede

The compiled library binaries (and associated PyRosetta tools and scripts) can be found in this directory in the rosetta source bundle:

`$ cd PATH/TO/ROSETTA/BUNDLE/main/source/build/PyRosetta/linux/namespace/release`

You should be able to copy the `release` directory to run wherever you need it, but to be safe I would recommend leaving the directory as is.  Note that the database files are symlinked in the `release` directory, so you'll need to copy them separately.

You now need to run the `SetPyRosettaEnvironment.sh` script in the `release` directory:

```
$ ./SetPyRosettaEnvironment.sh
Setting PyRosetta root as: $WORK/src/stampede_pyrosetta_local_compile/rosetta_src_2016.13.58602_bundle/main/source/build/PyRosetta/linux/namespace/release
Aliasing PyRosetta Toolkit GUI to pyrosetta_toolkit
$WORK/src/stampede_pyrosetta_local_compile/rosetta_src_2016.13.58602_bundle/main/source/build/PyRosetta/linux/namespace/release
```

Finally, to be able to load the PyRosetta module from any location, you need to include the following lines in your `~/.bash_profile` and `~/.idevrc` files:

```
export LD_LIBRARY_PATH=[PATH/TO/PYROSETTA/RELEASE/DIR]release/rosetta:$LD_LIBRARY_PATH
export PYTHONPATH=[PATH/TO/PYROSETTA/RELEASE/DIR]release/:$PYTHONPATH
```

### Check that installation worked

start an `idev` session on Stampede:

`$ idev`

load the Stampede python module and fire up iPython:

```
$ module load python
$ ipython
Python 2.7.12 |Anaconda 4.1.1 (x86_64)| (default, Jul  2 2016, 17:43:17)
Type "copyright", "credits" or "license" for more information.
IPython 4.2.0 -- An enhanced Interactive Python.
?         -> Introduction and overview of IPython's features.
%quickref -> Quick reference.
help      -> Python's own help system.
object?   -> Details about 'object', use 'object??' for extra details.

In [1]: 
```

Import PyRosetta and initialize:

```
In [1]: import rosetta # Note that this can take a few minutes

In [2]: rosetta.init()
PyRosetta 2014 [Rosetta 2014 unknown:exported] retrieved from: 
(C) Copyright Rosetta Commons Member Institutions.
Created in JHU by Sergey Lyskov and PyRosetta Team.

core.init: Rosetta version  from
core.init: command: PyRosetta -ex1 -ex2aro
core.init: 'RNG device' seed mode, using '/dev/urandom', seed=-566128349 seed_offset=0 real_seed=-566128349
core.init.random: RandomGenerator:init: Normal mode, seed=-566128349 RG_type=mt19937
core.init: Resolved executable path: /opt/apps/intel16/python/2.7.11/bin/python2.7
core.init: Looking for database based on location of executable: /opt/apps/intel16/python/2.7.11/bin/../database/

In [3]: 
```

And you should be good to go!

Note that you may need to specificy a path to the database when you initialize, e.g.:

```
In [3]: rosetta.init("-database=[PATH/TO/ROSETTA/DATABASE/]database")
PyRosetta 2014 [Rosetta 2014 unknown:exported] retrieved from: 
(C) Copyright Rosetta Commons Member Institutions.
Created in JHU by Sergey Lyskov and PyRosetta Team.

core.init: Rosetta version  from
core.init: command: PyRosetta -database=[PATH/TO/ROSETTA/DATABASE/]database
core.init: 'RNG device' seed mode, using '/dev/urandom', seed=1248401244 seed_offset=0 real_seed=1248401244
core.init.random: RandomGenerator:init: Normal mode, seed=1248401244 RG_type=mt19937

In [4]: 
```

This guide is still very much a work in progress, plz contact Colin with any issues, questions or suggestions

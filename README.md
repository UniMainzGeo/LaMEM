# LaMEM 1.1.0
## Lithosphere and Mantle Evolution Model

LaMEM is a parallel 3D numerical code that can be used to model various thermomechanical 
geodynamical processes such as mantle-lithosphere interaction for rocks 
that have visco-elasto-plastic rheologies. It was developed to better understand geological 
processes, particularly related to the dynamics of the crust and  lithosphere and their 
interaction with the mantle. It can, however, also be used to solve geomechanical problems, includes (compressible) poroelasticity, can be used to compute gravity anomalies and has an (adjoint) inversion framework. The code uses a marker-in-cell approach with a staggered finite difference discretization and is build on top of PETSc such that it can run on anything from a laptop to a massively parallel machine. 

A range of (Galerkin) multigrid and iterative solvers are 
available, for both linear and non-linear rheologies, using Picard and 
quasi-Newton solvers (provided through the PETSc interface).

LaMEM has been tested on a variety of machines ranging from laptops to a massively parallel cluster with 458'752 cores.

The main developed of the current version are:

* Anton Popov         (Johannes Gutenberg University Mainz, popov@uni-mainz.de), 2011-
* Boris Kaus          (JGU Mainz, kaus@uni-mainz.de), 2011-
* Tobias Baumann      (JGU Mainz), 2011-
* Georg Reuber        (JGU Mainz), 2015-	
* Adina Puesoek       (JGU Mainz, UC San Diego), 2012-
* Naiara Fernandez    (JGU Mainz), 2011-2014
* Arthur Bauville     (JGU Mainz), 2015
* Andrea Piccolo      (JGU Mainz), 2015-
* Beatriz Montesinos  (JGU Mainz), 2015-
* Maximilian Kottwitz (JGU Mainz), 2019-
* Arne Spang          (JGU Mainz), 2019-

Older versions of LaMEM included a finite element solver as well, 
and were developed by:

* Boris J.P. Kaus (ETH Zurich, Switzerland). 2007-2011
* Dave A. May     (ETH Zurich, Switzerland). 2009-2011

Development work was supported by the European Research Council, 
with ERC Starting Grant 258830, ERC Proof of Concept Grant 713397 and ERC Consolidator Grant 771143, as well as by BMBF projects SECURE and PERMEA awarded to Boris Kaus. 

LaMEM is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License,
or (at your option) any later version.

LaMEM is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with LaMEM. If not, see <http://www.gnu.org/licenses/>.

## 1. Directory structure
LaMEM consists of the following directories:
```
/doc          -  Some documentation (remains incomplete)
/matlab       -  Various matlab files to change initial mesh and read binary output into matlab.
/scripts      -	 Various scripts, currently mainly related to Paraview VTK files.
/src          -	 LaMEM source code; compile with "make mode=opt all" and "make mode=deb all"
/tests        -	 Directory with various benchmarks. 
/utils        -	 Various non-essential files.
/bin          -  Contains binaries after compilation (/deb contains the debug version and /opt the optimized)
/input_models -  Various input models (run with ../bin/LaMEM -ParamFile *.dat). See the README file in that directory.
```

## 2. Dependencies of LaMEM

### Main dependencies
LaMEM crucially relies on:

  * PETSc 3.16.4, ideally installed with the external packages SUPERLU_DIST, MUMPS and PASTIX

and to a minor extend on:

  * Python, to run the LaMEM testing environment, and see if things work as expected. 
  * Paraview, to visualize results.
  * GIT, in order to pull the latest version of LaMEM
  * MATLAB (version not important), in order to facilitate creating more complicated input geometries
  * geomIO, to create input geometries from Inkscape (see https://geomio.bitbucket.io) 
  * Any text editor, to modify the LaMEM input files. 
  * Visual Studio Code, in case you want to develop new code (some on the development team prefer Eclipse for this)

### Dependency installation
We develop LaMEM on Linux and Mac machines, but we also have had success on Windows 10, where we recommend installing it through the (new) bash shell. In order for LaMEM to work, you'll need to install the correct version of PETSc first. PETSc is usually not backwards compatible, so it won't work with the incorrect version. Also please note that we do not always immediately update LaMEM to the latest version of PETSc, so you may have to download/install an older version.

* PETSc: 
     See http://www.mcs.anl.gov/petsc/petsc-as/documentation/installation.html
     for installation instructions. We also provide some installation instructions on how to compile 
     PETSc (and mpich, gcc, as well as PETSc) on various machines in /doc/installation. The simplest manner is sometimes to let PETSc install all additional packages (like MPICH), but that often does not result in the most efficient code. You also have to make sure that the path is set correctly in that case (for MPICH, for example). On a large scale cluster, you will have to link against the cluster libraries (for MPI, for example).

* If you want to do debugging of LaMEM as well, it is a good idea to install both a DEBUG and an OPTIMIZED version of LaMEM, in separate directories.

* Nothing else - if the correct version of PETSc is installed, LaMEM will work fine.

	
## 2. Download and build LaMEM
* Download LaMEM from BitBucket, preferable by using GIT on your system:

```
        $ git clone https://bkaus@bitbucket.org/bkaus/lamem.git ./LaMEM
``` 

* Set the environment variables in your .bashrc or .bash_profile scripts such that the LaMEM makefile knows where to look for PETSc:
```
        export PETSC_OPT=DIRECTORY_WHERE_YOU_INSTALLED_YOUR_OPTIMIZED_PETSC
        export PETSC_DEB=DIRECTORY_WHERE_YOU_INSTALLED_YOUR_DEBUG_PETSC 
```

* To build the source in the /src directory:
```
        $ make mode=opt all 
```

* Once build, you can verify that LaMEM works correctly with:
```  
        $ cd /tests
        $ make test
```

  Note that we use the PythonTestHarness framework, which is a set of Python scripts that simplify regression testing (developed by Dave May and Patrick Sanan). The first time you run the makefile, it will download the git repository and put it in the directory ```./pythontestharness```. If this fails for some reason, you can download it directly from the Dave's bitbucket repository and put it in the directory. In that case, it will ask you which batch queuing system to use, which should be ```<none>```.	

  Ideally, the tests should run all with no problems. Depending on which external direct solver packages you installed, some may fail (for example, if you did not install PASTIX). The python test harness will give some hints as where the issues came from.  

  In case you happen to have MATLAB installed on your system, additional tests will be performed in which an input setup is generated from within MATLAB. In order for this to work, you will have to let the system know where the MATLAB-binary is located by setting the following environmental variable (which can also be set in your .bashrc file):
```  
        $ which matlab
        /path_to_your_matlab_library/matlab
        $ export MATLAB=/path_to_your_matlab_library/matlab
        $ make test
```
Note that you can look at the ```tests``` directory contains subdirectories that are named: ```t?_***``` (for example ```./t1_FB1_Direct/``` which contains a test for a Falling Block setup). Within each of these directories you will find a working LaMEM input file (```*.dat```), and a python file that runs the actual tests in that directory (such as ```test_1_FB1.py```). Have a look at these files to learn more on how to run LaMEM.


## 3. Getting started
#### 3a. First simulation
  You can run your first LaMEM simulation with 

```
      $ cd /input_models/BuildInSetups
      $ ../../bin/opt/LaMEM -ParamFile FallingBlock_IterativeSolver.dat
```
  which will run a setup that has a falling block of higher density embedded in a lower density fluid for 10 timesteps.  
  Run the same in parallel with:
  
``` 
     $ mpiexec -n 4 ../../bin/opt/LaMEM -ParamFile FallingBlock_IterativeSolver.dat
```

  You can visualize the results with Paraview by opening the file ```FB_iterative.pvd```, and pushing the play button (after selecting the appropriate view, such as surface and the appropriate field such as velocity).
  
  You can change the input parameters (such as resolution) by opening the file ```FallingBlock_IterativeSolver.dat``` with a texteditor and changing the parameters.

#### 3b. Learning more
 As we do not have an extensive user-guide yet (it takes time to create one, but will come at some point..), the best way to learn LaMEM is by looking at the input files in the order that is recommended in the README files. Start with ```/BuildInSetups```, which shows various example with geometries that are specified in the LaMEM input file. 

In addition, you can also look at the [Wiki](https://bitbucket.org/bkaus/lamem/wiki/Home) page (left menu). This will be the location where we will add more extensive documentation on how to use LaMEM.

All possible input parameters in LaMEM are listed in the file ```/input_models/input/lamem_input.dat```, which is worthwhile having a look at. Note that not all of these parameters have to be set (we select useful default options in most cases). 


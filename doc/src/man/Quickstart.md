# LaMEM _ Lithosphere and Mantle Evolution Model

LaMEM is a parallel 3D numerical code that can be used to model various thermomechanical 
geodynamical processes such as mantle-lithosphere interaction for rocks 
that have visco-elasto-plastic rheologies. It was developed to better understand geological 
processes, particularly related to the dynamics of the crust and  lithosphere and their 
interaction with the mantle. It can, however, also be used to solve geomechanical problems, includes (compressible) poroelasticity, can be used to compute gravity anomalies and has an (adjoint) inversion framework. The code uses a marker-in-cell approach with a staggered finite difference discretization and is build on top of PETSc such that it can run on anything from a laptop to a massively parallel machine. 

A range of (Galerkin) multigrid and iterative solvers are available, for both linear and non-linear rheologies, using Picard and  quasi-Newton solvers (provided through the PETSc interface).

LaMEM has been tested on a variety of machines ranging from laptops to a massively parallel cluster with 458'752 cores.


## Directory structure
LaMEM consists of the following directories:
```
/bin          -  Contains binaries after compilation (/deb contains the debug version and /opt the optimized)
/doc          -  Some documentation and the current webpage
/examples     -  Various input models (run with ../bin/LaMEM -ParamFile *.dat). See the README file in that directory
/info         -	 Installation instructions and input file nomenclature
/scripts      -	 Various scripts written in Julia and Bash
/src          -	 LaMEM source code; compile with "make mode=opt all" and "make mode=deb all"
/test         -	 Directory with the testing framework
```

## 2. Download and build LaMEM
You can download pre-compiled binaries of LaMEM and directly use it (also in parallel), through the [LaMEM julia](https://github.com/JuliaGeodynamics/LaMEM.jl) package, as explained below. You can also compile LaMEM manually to get the latest update of the code.

For changing input files, logging in to remote machines etc. etc, we recommend [Visual Studio Code](https://code.visualstudio.com), along with the remoteSSH and julia plugins.

## 2.1 Using pre-compiled binaries
The recommended way to install LaMEM on your machine (windows/mac/linux) is to use the julia package manager. For this download a recent version of [julia](https://julialang.org), and start it.
```julia
julia>]
pkg> add LaMEM
pkg> test LaMEM
```
Once this is installed, you can use it from within julia with:
```julia
julia> using LaMEM
```

(Note: use the backspace to go back from the package manager to the julia REPL.)
Running a simulation can be done with
```julia
julia> ParamFile="../../examples/BuiltInSetups/FallingBlock_Multigrid.dat";
julia> run_lamem(ParamFile, 2, "-nstep_max = 1")
```
This will run the simulation on 2 processors for 1 timestep. 
You do have to make sure that the path to the input `*.dat` file is correct in your system. 
Most likely you will have to change the path, which can be done with the built-in terminal (or PowerShell) in julia:
```julia
julia>; 
shell> cd ~/LaMEM/examples/BuiltInSetups/
```

If you wish, you can also directly run the downloaded binaries from your terminal without using julia. In that case you'll need to set the correct paths to the required binaries (`LaMEM`,`mpiexec`) and required dynamic libraries, which you can show with:
```julia
julia> show_paths_LaMEM()
LaMEM & mpiexec executables path : /Users/kausb/.julia/artifacts/26630bc992637321a5e5d3c0bc66005163370db6/bin:/Users/kausb/.julia/artifacts/483cb6f025b5a8266429afcb3f4ad498c58aaaee/bin
Dynamic libraries                : /Applications/Julia-1.7.app/Contents/Resources/julia/lib/julia:
```

The [LaMEM](https://github.com/JuliaGeodynamics/LaMEM.jl) julia package has a number of other functions that may come in handy, such as reading timesteps back into julia. 

If you want test some of the LaMEM examples in this repository, either clone the repo (below) or download it (three dots at the top of this page).

**Limitations** 
Whereas the pre-build libraries are quite handy, there are some limitations:

  * On Windows the MUMPS parallel direct solver is not available. SuperLU_dist does work, so we recommend using that instead.
  * On Mac and Linux both SuperLU_dist and MUMPS available.

## 2.2 Compiling it yourself
If want, you can ofcourse also compile LaMEM yourself, which will give you the latest version of the code. On large HPC clusters, this is often necessary as you need to link PETSc to the optimized MPI implementation on that system. 
### 2.2.1 Main dependencies
LaMEM crucially relies on:

  * PETSc, ideally installed with the external packages SUPERLU_DIST, MUMPS and PASTIX.

and to a minor extend on:

  * Julia, to run the LaMEM testing environment and see if things work as expected, and to facilitate creating more complicated input geometries.
  * ParaView, to visualize results.
  * GIT, in order to pull the latest version of LaMEM
  * Any text editor, to modify the LaMEM input files. 
  
### 2.2.2 Dependency installation
We develop LaMEM on Linux and Mac machines, but we also have had success on Windows 10, where we recommend installing it through the (new) bash shell. In order for LaMEM to work, you'll need to install the correct version of PETSc first. PETSc is usually not backwards compatible, so it won't work with the incorrect version. Also please note that we do not always immediately update LaMEM to the latest version of PETSc, so you may have to download/install an older version.

* PETSc: 
     [Go here](http://www.mcs.anl.gov/petsc/petsc-as/documentation/installation.html)
     to see the installation instructions. We also provide some installation instructions on how to compile 
     PETSc (and mpich, gcc, as well as PETSc) on various machines in ```/info/installation``` . The simplest manner is sometimes to let PETSc install all additional packages (like MPICH), but that often does not result in the most efficient code. You also have to make sure that the path is set correctly in that case (for MPICH, for example). On a large scale cluster, you will have to link against the cluster libraries (for MPI, for example).

* If you want to do debugging of LaMEM as well, it is a good idea to install both a DEBUG and an OPTIMIZED version of PETSc, in separate directories.

* Nothing else - if the correct version of PETSc is installed, LaMEM will work fine.

### 2.2.3 Install LaMEM

* Download LaMEM from GitHub, preferable by using GIT on your system:

```
        $ git clone https://github.com/UniMainzGeo/LaMEM.git ./LaMEM
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
        $ cd /test
        $ make test
```
  Ideally, the tests should run all with no problems. Depending on which external direct solver packages you installed, some may fail (for example, if you did not install MUMPS). The julia test framework will give some hints as where the issues came from.  

Note that you can look at the ```tests``` directory contains subdirectories that are named: ```t?_***``` (for example ```./t01_FB1_Direct/``` which contains a test for a Falling Block setup). Within each of these directories you will find a working LaMEM input file (```*.dat```), and a Julia script that runs the actual tests in that directory (such as ```test_01_FB1.py```). Have a look at these files to learn more on how to run LaMEM.


## 3. Getting started
### 3.1 First simulation
  You can run your first LaMEM simulation with 

```
      $ cd /examples/BuiltInSetups
      $ ../../bin/opt/LaMEM -ParamFile FallingBlock_IterativeSolver.dat
```
  which will run a setup that has a falling block of higher density embedded in a lower density fluid for 10 timesteps.  
  Run the same in parallel with:
  
``` 
     $ mpiexec -n 4 ../../bin/opt/LaMEM -ParamFile FallingBlock_IterativeSolver.dat
```

  You can visualize the results with Paraview by opening the file ```FB_iterative.pvd```, and pushing the play button (after selecting the appropriate view, such as surface and the appropriate field such as velocity).
  
  You can change the input parameters (such as resolution) by opening the file ```FallingBlock_IterativeSolver.dat``` with a texteditor and changing the parameters.

### 3.2 Learning more

You can also look at the [User Guide](https://unimainzgeo.github.io/LaMEM/dev/man/Home/). This website also contains more extensive documentation on how to use LaMEM.

The best way to learn LaMEM is by looking at the input files in the order that is recommended in the README files. Start with ```/BuiltInSetups```, which shows various example with geometries that are specified in the LaMEM input file. 

All possible input parameters in LaMEM are listed in the file ```/info/options/input_file.dat```, which is worthwhile having a look at. Note that not all of these parameters have to be set (we select useful default options in most cases). 


## 4. Development team
The main developers of the current version are:

* Anton Popov (JGU Mainz, Germany) [popov@uni-mainz.de]
* Boris Kaus  (JGU Mainz, Germany) [kaus@uni-mainz.de]

These people have contributed to LaMEM over the previous years:

* Tobias Baumann
* Adina Pusok
* Arthur Bauville
* Georg Reuber
* Andrea Piccolo
* Arne Spang
* Oskar Kottwitz
* Philipp Eichheimer
* Richard Spitz
* Christian Schuler
* Naiara Fernandez
* Marine Collignon
* Beatriz Martinez
* Jianfeng Yang
* Xinxin Wang
* Patrick Sanan
* Dave May
* Garrett Apuzen-Ito
* Jana Schierjott
* Sam Howell
* Wenrong Cao

Older versions of LaMEM including a finite element solver as well, were developed by:

* Boris J.P. Kaus (ETH Zurich, Switzerland). 
* Dave A. May     (ETH Zurich, Switzerland).

Development work was supported by the European Research Council, 
with ERC Starting Grant 258830, ERC Proof of Concept Grant 713397 and ERC Consolidator Grant 771143, as well as by BMBF projects SECURE and PERMEA awarded to Boris Kaus. 

LaMEM is free software: you can redistribute it and/or modify
it under the terms of the MIT License. ```LICENSE``` file for more information.

You should have received a copy of the GNU General Public License
along with LaMEM. If not, see <http://www.gnu.org/licenses/>.
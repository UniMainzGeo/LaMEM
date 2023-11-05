# LaMEM Userguide
  
 ![Getting Started](../assets/img/LaMEM_overview.png)

  LaMEM (Lithosphere and Mantle Evolution Model) is a software package to simulate 2D/3D geological and geomechanical processes, which runs on anything from your laptop to a massively parallel machine. It takes (poro)-visco-elasto-plastic rheologies into account, and can be used to simulate anything from the collision of tectonic plates to the flow of fluids through porous rocks. 

  The purpose of this wiki is to get you started with installing and running LaMEM and give a few worked-out examples that explain how to run LaMEM on your local machine or cluster.
  
## Contents

* [1. Installation](Installation.md) 
* [2. Getting Started](GettingStarted.md)
* [3. Initial Model Setup](InitialModelSetup.md)
* [4. Examples](Examples.md)
* [5. Features](Features.md)
* [6. LaMEM Development](LaMEM_Development.md)
* [7. LaMEM Debugging](Debugging.md)


## Features
LaMEM contains a number of features, specifically tailored to simulate geological processes and complex geometries:

* 2D/3D parallel thermomechanical code for cartesian geometries

* Build from the onset to run on MPI-parallel machines; the largest we tested had 458'752 processors
  
* Support for both direct solvers and multigrid solvers

* Marker and cell approach to simulate complex geometries and large strains

* Newton solvers for nonlinear iterations 

* Multiple ways to create model geometries: 
   (1) Build-in geometrical objects,
   (2) MATLAB/Octave input files, 
   (3) [GeomIO](https://geomio.bitbucket.io) support to create 2D/3D input geometries from vector graphics,
   (4) Voxel-based input (to compute effective permeabilities of porous rocks).

* Mechanical solver for visco-elasto-plastic solvers, for both (thermo)-elastic bulk compressible and incompressible cases
  
* Nonlinear combined rock creep laws and (regularized) non-associated plasticity

* Internal free surface and sticky air approach

* Energy solver with shear heating

* Fluid pressure and Darcy solver for groundwater flow

* Phase transitions, by taking (multiple) precomputed phase diagrams into account

* Partial melting

* Simplified erosion/sedimentation algorithms

* Breakpointing/restarting options

* Adjoint formulations to perform inversions and derive scaling laws 

## Getting started

We recommend that your start with reading the [installation](Installation.md) instructions.

## Extending this userguide

The userguide consists of [Markdown](http://daringfireball.net/projects/markdown/) pages which is compiled into webpages using the julia [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) package.  
The pages are listed in the 
```
/docs
```
directory of this repository. You can extend it by adding new pages to the repository, which can be added to the side menu by modifying `make.jl`. 
It will be automatically compiled when you push
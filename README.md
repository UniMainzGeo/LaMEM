# LaMEM
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://unimainzgeo.github.io/LaMEM/dev)

LaMEM (Lithosphere and Mantle Evolution Model) is a parallel 3D numerical code that can be used to simulate various thermo-mechanical 
geodynamical processes such as mantle-lithosphere interaction for rocks 
that have visco-elasto-plastic rheologies. It was developed to better understand geological 
processes, particularly related to the dynamics of the crust and  lithosphere and their 
interaction with the mantle. It can, however, also be used to solve geomechanical problems, includes (compressible) poroelasticity, can be used to compute gravity anomalies and has an (adjoint) inversion framework. The code uses a marker-in-cell approach with a staggered finite difference discretization and is build on top of PETSc such that it can run on anything from a laptop to a massively parallel machine. 

A range of (Galerkin) multigrid and iterative solvers are 
available, for both linear and non-linear rheologies, using Picard and 
quasi-Newton solvers (provided through the PETSc interface).

LaMEM has been tested on a variety of machines ranging from laptops to a massively parallel cluster with 458'752 cores.

 ![Getting Started](./docs/src/assets/img/LaMEM_overview.png)

## Getting started

Have a look at the [documentation](https://unimainzgeo.github.io/LaMEM/dev) on how to install the code and run it. 
You can also install and run a parallel version of LaMEM with the julia package [LaMEM.jl](https://github.com/JuliaGeodynamics/LaMEM.jl). 

## Development and funding
LaMEM is an open source code that was initially developed at the Johannes-Gutenberg University in Mainz (Germany) as part of the European Research Council Grants ERC StG 258830 (MODEL), ERC PoC 713397 (SALTED) and ERC CoG 771143 (MAGMA), as well as by BMBF projects SECURE and PERMEA. Many other colleagues have contributed since (see documentation).

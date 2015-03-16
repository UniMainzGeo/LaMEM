/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2012
 **        Dave A. May [dave.may@erdw.ethz.ch]
 **        Institute of Geophysics
 **        ETH Zürich
 **        Sonneggstrasse 5
 **        CH-8092 Zürich
 **        Switzerland
 **
 **    project:    pTatin3d
 **    filename:   sub_comm.h
 **
 **
 **    pTatin3d is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, either version 3 of the License,
 **    or (at your option) any later version.
 **
 **    pTatin3d is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with pTatin3d. If not, see <http://www.gnu.org/licenses/>.
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/


#ifndef __petsc_mpi_sub_comm_h__
#define __petsc_mpi_sub_comm_h__

#include <petsc.h>

typedef struct _p_PetscMPISubComm *PetscMPISubComm;

struct _p_PetscMPISubComm {
    MPI_Comm    parent_comm;
    MPI_Comm    sub_comm;
    PetscMPIInt nranks_from_parent;
    PetscMPIInt *ranks_from_parent;
    PetscBool   parent_rank_active_in_subcomm; /* 1:true, 0:false */
};

PetscErrorCode PetscMPISubCommCreate(MPI_Comm parent_comm,PetscInt parent_reduction_factor,PetscMPISubComm *scomm);
PetscErrorCode PetscMPISubCommDestroy(PetscMPISubComm *scomm);
PetscErrorCode PetscMPISubCommGetActive(PetscMPISubComm sc,PetscBool *a);
PetscErrorCode PetscMPISubCommGetParentComm(PetscMPISubComm sc,MPI_Comm *a);
PetscErrorCode PetscMPISubCommGetComm(PetscMPISubComm sc,MPI_Comm *a);
PetscErrorCode PetscMPISubCommGetNumSubRanks(PetscMPISubComm sc,PetscMPIInt *a);
PetscErrorCode PetscMPISubCommGetActiveRanks(PetscMPISubComm sc,PetscMPIInt **a);

#endif


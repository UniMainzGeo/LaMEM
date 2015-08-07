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
 **    filename:   sub_comm.c
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


#include "mpi.h"
#include "petsc.h"
#include "sub_comm.h"

#undef __FUNCT__
#define __FUNCT__ "_PetscMPISubCommCreate"
PetscErrorCode _PetscMPISubCommCreate(PetscMPISubComm *scomm)
{
    PetscMPISubComm comm;

    PetscMalloc(sizeof(struct _p_PetscMPISubComm),&comm);
    PetscMemzero(comm,sizeof(struct _p_PetscMPISubComm));

    *scomm = comm;

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PetscMPISubCommDestroy"
PetscErrorCode PetscMPISubCommDestroy(PetscMPISubComm *scomm)
{
    PetscMPISubComm comm;
    PetscErrorCode  ierr;

    if (scomm) { comm = *scomm; }
    else {  PetscFunctionReturn(0); }

    if (comm->sub_comm)          { ierr = MPI_Comm_free(&comm->sub_comm);CHKERRQ(ierr); }
    if (comm->ranks_from_parent) { PetscFree(comm->ranks_from_parent); }
    PetscFree(comm);

    *scomm = NULL;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PetscMPISubCommGetActive"
PetscErrorCode PetscMPISubCommGetActive(PetscMPISubComm sc,PetscBool *a)
{
    *a = sc->parent_rank_active_in_subcomm;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PetscMPISubCommGetParentComm"
PetscErrorCode PetscMPISubCommGetParentComm(PetscMPISubComm sc,MPI_Comm *a)
{
    *a = sc->parent_comm;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PetscMPISubCommGetComm"
PetscErrorCode PetscMPISubCommGetComm(PetscMPISubComm sc,MPI_Comm *a)
{
    *a = sc->sub_comm;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PetscMPISubCommGetNumSubRanks"
PetscErrorCode PetscMPISubCommGetNumSubRanks(PetscMPISubComm sc,PetscMPIInt *a)
{
    *a = sc->nranks_from_parent;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PetscMPISubCommGetActiveRanks"
PetscErrorCode PetscMPISubCommGetActiveRanks(PetscMPISubComm sc,PetscMPIInt **a)
{
    *a = sc->ranks_from_parent;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PetscMPISubCommCreate_Stride"
PetscErrorCode PetscMPISubCommCreate_Stride(MPI_Comm parent_comm,PetscInt parent_reduction_factor,PetscMPISubComm *scomm)
{
    PetscMPISubComm comm;
    MPI_Comm        sub_comm;
    MPI_Group       parent_group,sub_group;
    PetscMPIInt     nproc,rank,*subranks,nsubranks,c,i;
    PetscBool       active;
    PetscErrorCode  ierr;

    if (parent_reduction_factor < 1) {
        parent_reduction_factor = 1;
        PetscPrintf(parent_comm,"Warning:PetscMPISubCommCreate: parent_reduction_factor >=1\n");
    }

    ierr = MPI_Comm_size(parent_comm,&nproc);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(parent_comm,&rank);CHKERRQ(ierr);
    nsubranks = 0;
    for (i=0; i<nproc; i++) {
        if (i%parent_reduction_factor == 0) {
            nsubranks++;
        }
    }
    PetscMalloc(sizeof(PetscMPIInt)*nsubranks,&subranks);
    c = 0;
    for (i=0; i<nproc; i++) {
        if (i%parent_reduction_factor == 0) {
            subranks[c] = i;
            c++;
        }
    }
    active = PETSC_FALSE;
    for (i=0; i<nsubranks; i++) {
        if (rank == subranks[i]) {
            active = PETSC_TRUE;
        }
    }

    ierr = _PetscMPISubCommCreate(&comm);CHKERRQ(ierr);
    ierr = MPI_Comm_group(parent_comm,&parent_group);CHKERRQ(ierr);
    if (active) {
        ierr = MPI_Group_incl(parent_group, nsubranks, subranks, &sub_group);CHKERRQ(ierr);
    } else {
        ierr = MPI_Group_excl(parent_group, nsubranks, subranks, &sub_group);CHKERRQ(ierr);
    }
    ierr = MPI_Comm_create(parent_comm, sub_group, &sub_comm);CHKERRQ(ierr);

    /*
     {
     int sr,snp;
     ierr = MPI_Comm_size(sub_comm,&snp);
     ierr = MPI_Comm_rank(sub_comm,&sr);
     printf("parent[%d of %d]: sub[%d of %d]: active = %d \n",rank,nproc,sr,snp,active);
     }
     */

    comm->parent_comm        = parent_comm;
    comm->sub_comm           = sub_comm;
    comm->nranks_from_parent = nsubranks;
    comm->ranks_from_parent  = subranks;
    comm->parent_rank_active_in_subcomm = active;

    /* We can safely free group as its been embedded inside the sub_comm */
    ierr = MPI_Group_free(&sub_group);CHKERRQ(ierr);
    *scomm = comm;

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PetscMPISubCommCreate"
PetscErrorCode PetscMPISubCommCreate(MPI_Comm parent_comm,PetscInt parent_reduction_factor,PetscMPISubComm *scomm)
{
    PetscInt       creation_type;
    PetscErrorCode ierr;

    creation_type = 1;
    PetscOptionsGetInt(NULL,"-petsc_mpi_subcomm_creation_type",&creation_type,NULL);
    if (creation_type == 1) {
        ierr = PetscMPISubCommCreate_Stride(parent_comm,parent_reduction_factor,scomm);CHKERRQ(ierr);
    } else {
        SETERRQ(parent_comm,PETSC_ERR_USER,"Invalid creation type specified");
    }
    PetscFunctionReturn(0);
}

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
 **    filename:   pc_semiredundant.c
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

#include <petsc-private/pcimpl.h>     /*I "petscpc.h" I*/
#include <petscksp.h>           /*I "petscksp.h" I*/
#include "sub_comm.h"

typedef struct {
    PetscInt nsubcomm_factor;
    PetscMPIInt nsubcomm_size;
    PetscMPISubComm subcomm;
    Mat A,B,Ared,Bred;
    Vec xred,yred;
    KSP ksp;
    IS isin;
    VecScatter scatter;
    Vec xtmp;
    PetscBool fuse_blocks;
} PC_SemiRedundant;


/* helpers for gather matrix and scattering vectors */
#undef __FUNCT__
#define __FUNCT__ "MatCreateSemiRedundant"
PetscErrorCode MatCreateSemiRedundant(Mat A,PetscMPISubComm subcomm,MatReuse reuse,Mat *_red)
{
    PetscErrorCode ierr;
    IS             isrow,iscol;
    PetscInt       start,end,nr,nc;
    PetscInt       i,j,*nnz,*onnz;
    MPI_Comm       comm;
    Mat            Alocal,*_Alocal,red;
    
    PetscFunctionBegin;
    ierr = MatGetSize(A,&nr,&nc);CHKERRQ(ierr);
    ierr = PetscObjectGetComm((PetscObject)A,&comm);CHKERRQ(ierr);
    
    ierr = ISCreateStride(comm,nc,0,1,&iscol);CHKERRQ(ierr);
    
    start = end = 0;
    if (subcomm->parent_rank_active_in_subcomm) {
        
        if (reuse == MAT_INITIAL_MATRIX) {
            ierr = MatCreate(subcomm->sub_comm,&red);CHKERRQ(ierr);
            ierr = MatSetSizes(red,PETSC_DECIDE,PETSC_DECIDE,nr,nc);CHKERRQ(ierr);
            ierr = MatSetFromOptions(red);CHKERRQ(ierr);
            ierr = MatSetUp(red);CHKERRQ(ierr);
        } else {
            red = *_red;
        }
        
        ierr = MatGetOwnershipRange(red,&start,&end);CHKERRQ(ierr);
        ierr = ISCreateStride(comm,(end-start),start,1,&isrow);CHKERRQ(ierr);
    } else {
        /* if rank not in subcomm, just fetch a local chunk of A */
        ierr = MatGetOwnershipRange(A,&start,&end);CHKERRQ(ierr);
        ierr = ISCreateStride(comm,1,start,1,&isrow);CHKERRQ(ierr);
    }
    
    ierr = MatGetSubMatrices(A,1,&isrow,&iscol,MAT_INITIAL_MATRIX,&_Alocal);CHKERRQ(ierr);
    Alocal = *_Alocal;
    
    /* insert entries */
    if (subcomm->parent_rank_active_in_subcomm) {
        PetscInt ncols,startc,endc,rowidx;
        const PetscInt *cols;
        const PetscScalar *vals;
        
        /* preallocation */
        ierr = MatGetOwnershipRange(red,&start,&end);CHKERRQ(ierr);
        ierr = MatGetOwnershipRangeColumn(red,&startc,&endc);CHKERRQ(ierr);
        
        PetscMalloc(sizeof(PetscInt)*(end-start),&nnz);
        PetscMalloc(sizeof(PetscInt)*(end-start),&onnz);
        PetscMemzero(nnz,sizeof(PetscInt)*(end-start));
        PetscMemzero(onnz,sizeof(PetscInt)*(end-start));
        
        for (i=0; i<(end-start); i++) {
            ierr = MatGetRow(Alocal,i,&ncols,&cols,NULL);CHKERRQ(ierr);
            for (j=0; j<ncols; j++) {
                if ( (cols[j] >= startc) && (cols[j] < endc) ) {
                    nnz[i]++;
                } else {
                    onnz[i]++;
                }
            }
            ierr = MatRestoreRow(Alocal,i,&ncols,&cols,NULL);CHKERRQ(ierr);
        }
        ierr = MatSeqAIJSetPreallocation(red,PETSC_DEFAULT,nnz);CHKERRQ(ierr);
        ierr = MatMPIAIJSetPreallocation(red,PETSC_DEFAULT,nnz,PETSC_DEFAULT,onnz);CHKERRQ(ierr);
        
        PetscFree(nnz);
        PetscFree(onnz);
        
        /* insert */
        ierr = MatGetOwnershipRange(red,&start,&end);CHKERRQ(ierr);
        for (i=0; i<(end-start); i++) {
            ierr = MatGetRow(Alocal,i,&ncols,&cols,&vals);CHKERRQ(ierr);
            
            rowidx = i + start;
            ierr = MatSetValues(red,1,&rowidx,ncols,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
            
            ierr = MatRestoreRow(Alocal,i,&ncols,&cols,NULL);CHKERRQ(ierr);
        }
        
        ierr = MatAssemblyBegin(red,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(red,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    }
    
    ierr = ISDestroy(&isrow);CHKERRQ(ierr);
    ierr = ISDestroy(&iscol);CHKERRQ(ierr);
    ierr = MatDestroy(&Alocal);CHKERRQ(ierr);
    
    *_red = NULL;
    if (subcomm->parent_rank_active_in_subcomm) {
        *_red = red;
    }
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatCreateSemiRedundantFuseBlocks"
PetscErrorCode MatCreateSemiRedundantFuseBlocks(Mat A,PetscMPISubComm subcomm,MatReuse reuse,Mat *_red)
{
    PetscErrorCode ierr;
    IS             isrow,iscol;
    PetscInt       start,end,nr,nc,bsize;
    PetscInt       i,j,*nnz,*onnz;
    MPI_Comm       comm,comm_sub;
    Mat            Alocal,*_Alocal,red;
    const PetscInt *ranges;
    PetscInt       fused_length,offset;
    PetscMPIInt    nsize,nsize_sub,rank,rank_sub;
    
    PetscFunctionBegin;
    ierr = MatGetSize(A,&nr,&nc);CHKERRQ(ierr);
    ierr = MatGetBlockSize(A,&bsize);CHKERRQ(ierr);
    ierr = PetscObjectGetComm((PetscObject)A,&comm);CHKERRQ(ierr);
    ierr = PetscMPISubCommGetComm(subcomm,&comm_sub);CHKERRQ(ierr);
    ierr = MatGetOwnershipRanges(A,&ranges);CHKERRQ(ierr);
    
    ierr = MPI_Comm_size(comm,&nsize);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
    
    ierr = PetscMPISubCommGetNumSubRanks(subcomm,&nsize_sub);CHKERRQ(ierr);
    if (nsize%nsize_sub != 0) {
        //printf("nsize %d nsize_sub %d \n",nsize,nsize_sub);
        SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Cannot use -pc_semiredundant_fuse_blocks when comm->size is not exactly divisible by comm_sub->size");
    }
    //PetscPrintf(PETSC_COMM_WORLD,"A->bsize = %D \n",bsize);
    
    fused_length = 0;
    //for (i=0; i<nsize; i++) {
    //	PetscPrintf(PETSC_COMM_WORLD,"[%D]: start %D --->> end %D \n",i,ranges[i],ranges[i+1]);
    //}
    
    ierr = ISCreateStride(comm,nc,0,1,&iscol);CHKERRQ(ierr);
    
    start = end = 0;
    if (subcomm->parent_rank_active_in_subcomm) {
        PetscInt f0,f1;
        
        ierr = MPI_Comm_rank(comm_sub,&rank_sub);CHKERRQ(ierr);
        
        if (reuse == MAT_INITIAL_MATRIX) {
            
            offset = nsize/nsize_sub;
            //PetscPrintf(PETSC_COMM_SELF,"[%D,%D]: offset %D\n",rank,rank_sub,offset);
            f0 = ranges[ offset * rank_sub];
            f1 = ranges[ offset * (rank_sub + 1)];
            fused_length = f1 - f0;
            //PetscPrintf(PETSC_COMM_SELF,"[%D,%D]: fused [%D -- %D] \n",rank,rank_sub,f0,f1);
            //PetscPrintf(PETSC_COMM_SELF,"[%D,%D]: fused length %D \n",rank,rank_sub,fused_length);
            
            ierr = MatCreate(subcomm->sub_comm,&red);CHKERRQ(ierr);
            //ierr = MatSetSizes(red,PETSC_DECIDE,PETSC_DECIDE,nr,nc);CHKERRQ(ierr);
            ierr = MatSetSizes(red,fused_length,fused_length,PETSC_DETERMINE,PETSC_DETERMINE);CHKERRQ(ierr);
            ierr = MatSetFromOptions(red);CHKERRQ(ierr);
        } else {
            red = *_red;
        }
        
        ierr = MatGetOwnershipRange(red,&start,&end);CHKERRQ(ierr);
        ierr = ISCreateStride(comm,(end-start),start,1,&isrow);CHKERRQ(ierr);
    } else {
        /* if rank not in subcomm, just fetch a local chunk of A */
        ierr = MatGetOwnershipRange(A,&start,&end);CHKERRQ(ierr);
        ierr = ISCreateStride(comm,1,start,1,&isrow);CHKERRQ(ierr);
    }
    
    ierr = MatGetSubMatrices(A,1,&isrow,&iscol,MAT_INITIAL_MATRIX,&_Alocal);CHKERRQ(ierr);
    Alocal = *_Alocal;
    
    /* insert entries */
    if (subcomm->parent_rank_active_in_subcomm) {
        PetscInt ncols,startc,endc,rowidx;
        const PetscInt *cols;
        const PetscScalar *vals;
        
        /* preallocation */
        ierr = MatGetOwnershipRange(red,&start,&end);CHKERRQ(ierr);
        ierr = MatGetOwnershipRangeColumn(red,&startc,&endc);CHKERRQ(ierr);
        
        PetscMalloc(sizeof(PetscInt)*(end-start),&nnz);
        PetscMalloc(sizeof(PetscInt)*(end-start),&onnz);
        PetscMemzero(nnz,sizeof(PetscInt)*(end-start));
        PetscMemzero(onnz,sizeof(PetscInt)*(end-start));
        
        for (i=0; i<(end-start); i++) {
            ierr = MatGetRow(Alocal,i,&ncols,&cols,NULL);CHKERRQ(ierr);
            for (j=0; j<ncols; j++) {
                if ( (cols[j] >= startc) && (cols[j] < endc) ) {
                    nnz[i]++;
                } else {
                    onnz[i]++;
                }
            }
            ierr = MatRestoreRow(Alocal,i,&ncols,&cols,NULL);CHKERRQ(ierr);
        }
        ierr = MatSeqAIJSetPreallocation(red,PETSC_DEFAULT,nnz);CHKERRQ(ierr);
        ierr = MatMPIAIJSetPreallocation(red,PETSC_DEFAULT,nnz,PETSC_DEFAULT,onnz);CHKERRQ(ierr);
        
        PetscFree(nnz);
        PetscFree(onnz);
        
        /* insert */
        ierr = MatGetOwnershipRange(red,&start,&end);CHKERRQ(ierr);
        for (i=0; i<(end-start); i++) {
            ierr = MatGetRow(Alocal,i,&ncols,&cols,&vals);CHKERRQ(ierr);
            
            rowidx = i + start;
            ierr = MatSetValues(red,1,&rowidx,ncols,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
            
            ierr = MatRestoreRow(Alocal,i,&ncols,&cols,NULL);CHKERRQ(ierr);
        }
        
        ierr = MatAssemblyBegin(red,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(red,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    }
    
    ierr = ISDestroy(&isrow);CHKERRQ(ierr);
    ierr = ISDestroy(&iscol);CHKERRQ(ierr);
    ierr = MatDestroy(&Alocal);CHKERRQ(ierr);
    
    *_red = NULL;
    if (subcomm->parent_rank_active_in_subcomm) {
        *_red = red;
    }
    PetscFunctionReturn(0);
}

/* implementations for SemiRedundant */
#undef __FUNCT__
#define __FUNCT__ "PCSetUp_SemiRedundant"
static PetscErrorCode PCSetUp_SemiRedundant(PC pc)
{
    PetscErrorCode   ierr;
    PC_SemiRedundant *red = (PC_SemiRedundant*)pc->data;
    MPI_Comm     comm;
    PetscMPISubComm  sc;
    MatReuse     reuse;
    PetscErrorCode (*fp_matcreate)(Mat,PetscMPISubComm,MatReuse,Mat*);
    
    PetscFunctionBegin;
    /* construction phase */
    if (!pc->setupcalled) {
        
        /* set up sub communicator */
        ierr = PetscObjectGetComm((PetscObject)pc,&comm);CHKERRQ(ierr);
        ierr = PetscMPISubCommCreate(comm,red->nsubcomm_factor,&sc);CHKERRQ(ierr);
        red->subcomm = sc;
        ierr = MPI_Comm_size(sc->sub_comm,&red->nsubcomm_size);CHKERRQ(ierr);
        
        red->Ared = NULL;
        red->Bred = NULL;
        ierr = PCGetOperators(pc,&red->A,&red->B);CHKERRQ(ierr);
        
        red->ksp = NULL;
        if (red->subcomm->parent_rank_active_in_subcomm) {
            const char     *prefix;
            
            ierr = KSPCreate(red->subcomm->sub_comm,&red->ksp);CHKERRQ(ierr);
            ierr = PetscObjectIncrementTabLevel((PetscObject)red->ksp,(PetscObject)pc,1);CHKERRQ(ierr);
            ierr = PetscLogObjectParent((PetscObject)pc,(PetscObject)red->ksp);CHKERRQ(ierr);
            
            ierr = PCGetOptionsPrefix(pc,&prefix);CHKERRQ(ierr);
            ierr = KSPSetOptionsPrefix(red->ksp,prefix);CHKERRQ(ierr);
            ierr = KSPAppendOptionsPrefix(red->ksp,"semiredundant_");CHKERRQ(ierr);
        }
    }
    
    fp_matcreate = MatCreateSemiRedundant;
    if (red->fuse_blocks) {
        fp_matcreate = MatCreateSemiRedundantFuseBlocks;
    }
    
    /* fetch redundant matrix */
    if (!pc->setupcalled) {
        
        ierr = fp_matcreate(red->A,red->subcomm,MAT_INITIAL_MATRIX,&red->Ared);CHKERRQ(ierr);
        if (red->A != red->B) {
            ierr = fp_matcreate(red->B,red->subcomm,MAT_INITIAL_MATRIX,&red->Bred);CHKERRQ(ierr);
        }
        
        if ( (red->A == red->B) && (red->Ared) ) {
            red->Bred = red->Ared;
            ierr = PetscObjectReference((PetscObject)red->Ared);CHKERRQ(ierr);
        }
        
        red->xred = NULL;
        red->yred = NULL;
        if (red->Ared) {
            PetscInt m,n;
            
            ierr = MatGetLocalSize(red->Ared,&m,&n);CHKERRQ(ierr);
            /* create xred with empty local arrays, because xdup's arrays will be placed into it */
            ierr = VecCreateMPIWithArray(red->subcomm->sub_comm,1,m,PETSC_DECIDE,NULL,&red->xred);CHKERRQ(ierr);
            
            ierr = MatGetVecs(red->Ared,NULL,&red->yred);CHKERRQ(ierr);
        }
        
    } else {
        reuse = MAT_REUSE_MATRIX;
        
        ierr = fp_matcreate(red->A,red->subcomm,reuse,&red->Ared);CHKERRQ(ierr);
        if (red->A != red->B) {
            ierr = fp_matcreate(red->B,red->subcomm,reuse,&red->Bred);CHKERRQ(ierr);
        } else {
            red->Bred = red->Ared;
        }
    }
    
    /* setup scatters */
    if (!pc->setupcalled) {
        PetscInt st,ed;
        PetscInt n,N;
        Vec      x;
        
        ierr = PetscObjectGetComm((PetscObject)pc,&comm);CHKERRQ(ierr);
        ierr = MatGetVecs(red->A,&x,NULL);CHKERRQ(ierr);
        
        if (red->xred) {
            ierr = VecGetOwnershipRange(red->xred,&st,&ed);CHKERRQ(ierr);
            ierr = ISCreateStride(comm,ed-st,st,1,&red->isin);CHKERRQ(ierr);
        } else {
            ierr = VecGetOwnershipRange(x,&st,&ed);CHKERRQ(ierr);
            ierr = ISCreateStride(comm,1,st,1,&red->isin);CHKERRQ(ierr);
        }
        
        ierr = ISGetLocalSize(red->isin,&n);CHKERRQ(ierr);
        ierr = ISGetSize(red->isin,&N);CHKERRQ(ierr);
        ierr = VecCreate(PetscObjectComm((PetscObject)red->isin),&red->xtmp);CHKERRQ(ierr);
        ierr = VecSetSizes(red->xtmp,n,N);CHKERRQ(ierr);
        ierr = VecSetType(red->xtmp,((PetscObject)x)->type_name);CHKERRQ(ierr);
        ierr = VecScatterCreate(x,red->isin,red->xtmp,NULL,&red->scatter);CHKERRQ(ierr);
        
        ierr = VecDestroy(&x);CHKERRQ(ierr);
    }
    
    /* common - no construction */
    ierr = PCGetOperators(pc,&red->A,&red->B);CHKERRQ(ierr);
    if (red->Ared) {
        //MatView(red->Ared,PETSC_VIEWER_STDOUT_(red->subcomm->sub_comm));
        ierr = KSPSetOperators(red->ksp,red->Ared,red->Bred);CHKERRQ(ierr);
        
        if (pc->setfromoptionscalled){
            ierr = KSPSetFromOptions(red->ksp);CHKERRQ(ierr);
        }
        
        ierr = KSPSetUp(red->ksp);CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCApply_SemiRedundant"
static PetscErrorCode PCApply_SemiRedundant(PC pc,Vec x,Vec y)
{
    PetscErrorCode   ierr;
    PC_SemiRedundant *red = (PC_SemiRedundant*)pc->data;
    Vec         xtmp;
    PetscInt    i,st,ed;
    VecScatter  scatter;
    PetscScalar *array;
    
    PetscFunctionBegin;
    xtmp    = red->xtmp;
    scatter = red->scatter;
    
    /* pull in vector */
    ierr = VecScatterBegin(scatter,x,xtmp,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterEnd(scatter,x,xtmp,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    
    /* solve */
    if (red->ksp) {
        
        /* define xred */
        ierr = VecGetArray(xtmp,&array);CHKERRQ(ierr);
        /* we created xred with empty local arrays, now we fill it in */
        ierr = VecPlaceArray(red->xred,(const PetscScalar*)array);CHKERRQ(ierr);
        
        ierr = KSPSolve(red->ksp,red->xred,red->yred);CHKERRQ(ierr);
        
        ierr = VecResetArray(red->xred);CHKERRQ(ierr);
        ierr = VecRestoreArray(xtmp,&array);CHKERRQ(ierr);
    }
    
    /* return vector */
    ierr = VecGetArray(xtmp,&array);CHKERRQ(ierr);
    if (red->yred) {
        PetscScalar *LA_yred;
        
        ierr = VecGetOwnershipRange(red->yred,&st,&ed);CHKERRQ(ierr);
        
        ierr = VecGetArray(red->yred,&LA_yred);CHKERRQ(ierr);
        for (i=0; i<ed-st; i++) {
            array[i] = LA_yred[i];
        }
        ierr = VecRestoreArray(red->yred,&LA_yred);CHKERRQ(ierr);
    } else {
        array[0] = 0.0;
    }
    ierr = VecRestoreArray(xtmp,&array);CHKERRQ(ierr);
    
    ierr = VecScatterBegin(scatter,xtmp,y,ADD_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
    ierr = VecScatterEnd(scatter,xtmp,y,ADD_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCReset_SemiRedundant"
static PetscErrorCode PCReset_SemiRedundant(PC pc)
{
    PetscErrorCode   ierr;
    PC_SemiRedundant *red = (PC_SemiRedundant*)pc->data;
    
    PetscFunctionBegin;
    if (red->xred) { ierr = VecDestroy(&red->xred);CHKERRQ(ierr);}
    if (red->yred) { ierr = VecDestroy(&red->yred);CHKERRQ(ierr);}
    if (red->Ared) { ierr = MatDestroy(&red->Ared);CHKERRQ(ierr);}
    if (red->Bred) { ierr = MatDestroy(&red->Bred);CHKERRQ(ierr);}
    
    ierr = ISDestroy(&red->isin);CHKERRQ(ierr);
    ierr = VecDestroy(&red->xtmp);CHKERRQ(ierr);
    ierr = VecScatterDestroy(&red->scatter);CHKERRQ(ierr);
    if (red->ksp) {  ierr = KSPReset(red->ksp);CHKERRQ(ierr);}
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCDestroy_SemiRedundant"
static PetscErrorCode PCDestroy_SemiRedundant(PC pc)
{
    PetscErrorCode   ierr;
    PC_SemiRedundant *red = (PC_SemiRedundant*)pc->data;
    
    PetscFunctionBegin;
    ierr = PCReset_SemiRedundant(pc);CHKERRQ(ierr);
    if (red->ksp) {
        ierr = KSPDestroy(&red->ksp);CHKERRQ(ierr);
    }
    ierr = PetscMPISubCommDestroy(&red->subcomm);CHKERRQ(ierr);
    ierr = PetscFree(pc->data);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCSetFromOptions_SemiRedundant"
static PetscErrorCode PCSetFromOptions_SemiRedundant(PC pc)
{
    PetscErrorCode   ierr;
    PC_SemiRedundant *red = (PC_SemiRedundant*)pc->data;
    
    PetscFunctionBegin;
    ierr = PetscOptionsHead("SemiRedundant options");CHKERRQ(ierr);
    ierr = PetscOptionsInt("-pc_semiredundant_factor","Factor to reduce parent communication size by","PCSemiRedundantSetFactor",red->nsubcomm_factor,&red->nsubcomm_factor,0);CHKERRQ(ierr);
    ierr = PetscOptionsBool("-pc_semiredundant_fuse_blocks","Fuse original matrix partitioning and preserve block size","PCSemiRedundantFuseBlocks",red->fuse_blocks,&red->fuse_blocks,0);CHKERRQ(ierr);
    ierr = PetscOptionsTail();CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PCView_SemiRedundant"
static PetscErrorCode PCView_SemiRedundant(PC pc,PetscViewer viewer)
{
    PetscErrorCode   ierr;
    PC_SemiRedundant *red = (PC_SemiRedundant*)pc->data;
    PetscBool        iascii,isstring;
    PetscViewer      subviewer;
    
    PetscFunctionBegin;
    ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&iascii);CHKERRQ(ierr);
    ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERSTRING,&isstring);CHKERRQ(ierr);
    if (iascii) {
        if (!red->subcomm) {
            ierr = PetscViewerASCIIPrintf(viewer,"  SemiRedundant: preconditioner not yet setup\n");CHKERRQ(ierr);
        } else {
            
            /* ierr = PetscViewerASCIIPrintf(viewer,"  SemiRedundant preconditioner:\n");CHKERRQ(ierr); */
            ierr = PetscViewerASCIIPrintf(viewer,"  SemiRedundant: parent comm size reduction factor = %D\n",red->nsubcomm_factor);CHKERRQ(ierr);
            ierr = PetscViewerASCIIPrintf(viewer,"  SemiRedundant: subcomm_size = %D\n",red->nsubcomm_size);CHKERRQ(ierr);
            if (!red->fuse_blocks) {
                ierr = PetscViewerASCIIPrintf(viewer,"  SemiRedundant: not preserving original matrix partition boundaries\n");CHKERRQ(ierr);
            } else {
                ierr = PetscViewerASCIIPrintf(viewer,"  SemiRedundant: preserving original matrix partition boundaries\n");CHKERRQ(ierr);
            }
            ierr = PetscViewerGetSubcomm(viewer,red->subcomm->sub_comm,&subviewer);CHKERRQ(ierr);
            ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
            if (red->subcomm->parent_rank_active_in_subcomm) {
                ierr = KSPView(red->ksp,subviewer);CHKERRQ(ierr);
            }
            ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
            ierr = PetscViewerRestoreSubcomm(viewer,red->subcomm->sub_comm,&subviewer);CHKERRQ(ierr);
        }
    } else if (isstring) { 
        ierr = PetscViewerStringSPrintf(viewer," SemiRedundant preconditioner");CHKERRQ(ierr);
    } else {
        SETERRQ1(PetscObjectComm((PetscObject)pc),PETSC_ERR_SUP,"Viewer type %s not supported for PC SemiRedundant",((PetscObject)viewer)->type_name);
    }
    PetscFunctionReturn(0);
}

EXTERN_C_BEGIN
#undef __FUNCT__  
#define __FUNCT__ "PCCreate_SemiRedundant"
PetscErrorCode PCCreate_SemiRedundant(PC pc)
{
    PetscErrorCode   ierr;
    PC_SemiRedundant *red;
    PetscMPIInt      size;
    
    PetscFunctionBegin;
    ierr = PetscNewLog(pc,&red);CHKERRQ(ierr);
    pc->data            = (void*)red; 
    
    red->nsubcomm_factor = 1;
    ierr = MPI_Comm_size(PetscObjectComm((PetscObject)pc),&size);CHKERRQ(ierr);
    red->nsubcomm_size   = size;
    red->fuse_blocks     = PETSC_FALSE;
    
    pc->ops->apply           = PCApply_SemiRedundant;
    pc->ops->applytranspose  = 0;
    pc->ops->setup           = PCSetUp_SemiRedundant;
    pc->ops->destroy         = PCDestroy_SemiRedundant;
    pc->ops->reset           = PCReset_SemiRedundant;
    pc->ops->setfromoptions  = PCSetFromOptions_SemiRedundant;
    pc->ops->view            = PCView_SemiRedundant;    
    
    PetscFunctionReturn(0);
}
EXTERN_C_END


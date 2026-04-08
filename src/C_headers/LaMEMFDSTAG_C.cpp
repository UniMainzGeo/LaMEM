#include "LaMEMFDSTAG_C.h"
#include "../LaMEM.h"
#include "../cvi.h"
#include "../JacRes.h"
#include "../fdstag.h"
#include "../advect.h"
#include "../surf.h"
#include "../bc.h"
#include "../tssolve.h"
#include "../tools.h"

extern "C" {

// MeshSeg1D functions
PetscErrorCode LaMEMFDSTAG_MeshSeg1DReadParam(void *ms, PetscScalar leng, PetscScalar gtol, 
                                              const char *dir, void *fb) {
    return MeshSeg1DReadParam((MeshSeg1D*)ms, leng, gtol, dir, (FB*)fb);
}

PetscErrorCode LaMEMFDSTAG_MeshSeg1DGenCoord(void *ms, PetscInt iseg, PetscInt nl, 
                                             PetscInt istart, PetscScalar *crd) {
    return MeshSeg1DGenCoord((MeshSeg1D*)ms, iseg, nl, istart, crd);
}

// Discret1D functions
PetscErrorCode LaMEMFDSTAG_Discret1DCreate(void *ds, PetscInt nproc, PetscInt rank, 
                                           PetscInt *nnodProc, PetscInt color, 
                                           PetscMPIInt grprev, PetscMPIInt grnext, 
                                           PetscScalar gtol) {
    return Discret1DCreate((Discret1D*)ds, nproc, rank, nnodProc, color, grprev, grnext, gtol);
}

PetscErrorCode LaMEMFDSTAG_Discret1DDestroy(void *ds) {
    return Discret1DDestroy((Discret1D*)ds);
}

PetscErrorCode LaMEMFDSTAG_Discret1DReadRestart(void *ds, FILE *fp) {
    return Discret1DReadRestart((Discret1D*)ds, fp);
}

PetscErrorCode LaMEMFDSTAG_Discret1DWriteRestart(void *ds, FILE *fp) {
    return Discret1DWriteRestart((Discret1D*)ds, fp);
}

PetscErrorCode LaMEMFDSTAG_Discret1DGetNumCells(void *ds, PetscInt **ncelProc) {
    return Discret1DGetNumCells((Discret1D*)ds, ncelProc);
}

PetscErrorCode LaMEMFDSTAG_Discret1DGenCoord(void *ds, void *ms) {
    return Discret1DGenCoord((Discret1D*)ds, (MeshSeg1D*)ms);
}

PetscErrorCode LaMEMFDSTAG_Discret1DStretch(void *ds, PetscScalar eps, PetscScalar ref) {
    return Discret1DStretch((Discret1D*)ds, eps, ref);
}

PetscErrorCode LaMEMFDSTAG_Discret1DGetColumnComm(void *ds) {
    return Discret1DGetColumnComm((Discret1D*)ds);
}

PetscErrorCode LaMEMFDSTAG_Discret1DFreeColumnComm(void *ds) {
    return Discret1DFreeColumnComm((Discret1D*)ds);
}

PetscErrorCode LaMEMFDSTAG_Discret1DGatherCoord(void *ds, PetscScalar **coord) {
    return Discret1DGatherCoord((Discret1D*)ds, coord);
}

PetscErrorCode LaMEMFDSTAG_Discret1DCheckMG(void *ds, const char *dir, PetscInt *_ncors) {
    return Discret1DCheckMG((Discret1D*)ds, dir, _ncors);
}

PetscErrorCode LaMEMFDSTAG_Discret1DgetMaxInvStep(void *ds, DM da, Vec gv, PetscInt dir, PetscScalar *_idtmax) {
    return Discret1DgetMaxInvStep((Discret1D*)ds, da, gv, dir, _idtmax);
}

PetscErrorCode LaMEMFDSTAG_Discret1DFindPoint(void *ds, PetscScalar x, PetscInt *ID) {
    return Discret1DFindPoint((Discret1D*)ds, x, *ID);
}

// DOFIndex functions
PetscErrorCode LaMEMFDSTAG_DOFIndexCreate(void *dof, DM DA_CEN, DM DA_X, DM DA_Y, DM DA_Z) {
    return DOFIndexCreate((DOFIndex*)dof, DA_CEN, DA_X, DA_Y, DA_Z);
}

PetscErrorCode LaMEMFDSTAG_DOFIndexDestroy(void *dof) {
    return DOFIndexDestroy((DOFIndex*)dof);
}

PetscErrorCode LaMEMFDSTAG_DOFIndexCompute(void *dof, PetscInt idxmod) {
    idxtype idx_type;
    switch(idxmod) {
        case 0: idx_type = IDXNONE; break;       // LAMEM_IDXNONE
        case 1: idx_type = IDXCOUPLED; break;    // LAMEM_IDXCOUPLED
        case 2: idx_type = IDXUNCOUPLED; break;  // LAMEM_IDXUNCOUPLED
        default: return PETSC_ERR_ARG_OUTOFRANGE;
    }
    return DOFIndexCompute((DOFIndex*)dof, idx_type);
}

// Main FDSTAG functions
PetscErrorCode LaMEMFDSTAG_Create(void *fs, void *fb) {
    return FDSTAGCreate((FDSTAG*)fs, (FB*)fb);
}

PetscErrorCode LaMEMFDSTAG_Destroy(void *fs) {
    return FDSTAGDestroy((FDSTAG*)fs);
}

PetscErrorCode LaMEMFDSTAG_ReadRestart(void *fs, FILE *fp) {
    return FDSTAGReadRestart((FDSTAG*)fs, fp);
}

PetscErrorCode LaMEMFDSTAG_WriteRestart(void *fs, FILE *fp) {
    return FDSTAGWriteRestart((FDSTAG*)fs, fp);
}

PetscErrorCode LaMEMFDSTAG_CreateDMDA(void *fs, PetscInt Nx, PetscInt Ny, PetscInt Nz,
                                      PetscInt Px, PetscInt Py, PetscInt Pz,
                                      PetscInt *lx, PetscInt *ly, PetscInt *lz) {
    return FDSTAGCreateDMDA((FDSTAG*)fs, Nx, Ny, Nz, Px, Py, Pz, lx, ly, lz);
}

PetscErrorCode LaMEMFDSTAG_GetNeighbProc(void *fs) {
    return FDSTAGGetNeighbProc((FDSTAG*)fs);
}

PetscErrorCode LaMEMFDSTAG_GetPointRanks(void *fs, PetscScalar *X, PetscInt *lrank, PetscMPIInt *grank) {
    return FDSTAGGetPointRanks((FDSTAG*)fs, X, lrank, grank);
}

PetscErrorCode LaMEMFDSTAG_GetAspectRatio(void *fs, PetscScalar *maxAspRat) {
    return FDSTAGGetAspectRatio((FDSTAG*)fs, maxAspRat);
}

PetscErrorCode LaMEMFDSTAG_View(void *fs) {
    return FDSTAGView((FDSTAG*)fs);
}

PetscErrorCode LaMEMFDSTAG_GetLocalBox(void *fs, PetscScalar *bx, PetscScalar *by, PetscScalar *bz,
                                       PetscScalar *ex, PetscScalar *ey, PetscScalar *ez) {
    return FDSTAGGetLocalBox((FDSTAG*)fs, bx, by, bz, ex, ey, ez);
}

PetscErrorCode LaMEMFDSTAG_GetGlobalBox(void *fs, PetscScalar *bx, PetscScalar *by, PetscScalar *bz,
                                        PetscScalar *ex, PetscScalar *ey, PetscScalar *ez) {
    return FDSTAGGetGlobalBox((FDSTAG*)fs, bx, by, bz, ex, ey, ez);
}

PetscErrorCode LaMEMFDSTAG_SaveGrid(void *fs) {
    return FDSTAGSaveGrid((FDSTAG*)fs);
}

// Wrapper functions
PetscErrorCode LaMEMFDSTAG_DMDACreate3dSetUp(MPI_Comm comm, DMBoundaryType bx, DMBoundaryType by, 
                                             DMBoundaryType bz, DMDAStencilType stencil_type,
                                             PetscInt M, PetscInt N, PetscInt P, PetscInt m, PetscInt n, PetscInt p,
                                             PetscInt dof, PetscInt s, const PetscInt lx[], const PetscInt ly[], 
                                             const PetscInt lz[], DM *da) {
    return DMDACreate3dSetUp(comm, bx, by, bz, stencil_type, M, N, P, m, n, p, dof, s, lx, ly, lz, da);
}

}
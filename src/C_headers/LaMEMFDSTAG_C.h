#ifndef LAMEMFDSTAG_C_H
#define LAMEMFDSTAG_C_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

// Index types
#define LAMEM_IDXNONE       0
#define LAMEM_IDXCOUPLED    1
#define LAMEM_IDXUNCOUPLED  2

// MeshSeg1D functions
PetscErrorCode LaMEMFDSTAG_MeshSeg1DReadParam(void *ms, PetscScalar leng, PetscScalar gtol, 
                                              const char *dir, void *fb);
PetscErrorCode LaMEMFDSTAG_MeshSeg1DGenCoord(void *ms, PetscInt iseg, PetscInt nl, 
                                             PetscInt istart, PetscScalar *crd);

// Discret1D functions
PetscErrorCode LaMEMFDSTAG_Discret1DCreate(void *ds, PetscInt nproc, PetscInt rank, 
                                           PetscInt *nnodProc, PetscInt color, 
                                           PetscMPIInt grprev, PetscMPIInt grnext, 
                                           PetscScalar gtol);
PetscErrorCode LaMEMFDSTAG_Discret1DDestroy(void *ds);
PetscErrorCode LaMEMFDSTAG_Discret1DReadRestart(void *ds, FILE *fp);
PetscErrorCode LaMEMFDSTAG_Discret1DWriteRestart(void *ds, FILE *fp);
PetscErrorCode LaMEMFDSTAG_Discret1DGetNumCells(void *ds, PetscInt **ncelProc);
PetscErrorCode LaMEMFDSTAG_Discret1DGenCoord(void *ds, void *ms);
PetscErrorCode LaMEMFDSTAG_Discret1DStretch(void *ds, PetscScalar eps, PetscScalar ref);
PetscErrorCode LaMEMFDSTAG_Discret1DGetColumnComm(void *ds);
PetscErrorCode LaMEMFDSTAG_Discret1DFreeColumnComm(void *ds);
PetscErrorCode LaMEMFDSTAG_Discret1DGatherCoord(void *ds, PetscScalar **coord);
PetscErrorCode LaMEMFDSTAG_Discret1DCheckMG(void *ds, const char *dir, PetscInt *_ncors);
PetscErrorCode LaMEMFDSTAG_Discret1DgetMaxInvStep(void *ds, DM da, Vec gv, PetscInt dir, PetscScalar *_idtmax);
PetscErrorCode LaMEMFDSTAG_Discret1DFindPoint(void *ds, PetscScalar x, PetscInt *ID);

// DOFIndex functions
PetscErrorCode LaMEMFDSTAG_DOFIndexCreate(void *dof, DM DA_CEN, DM DA_X, DM DA_Y, DM DA_Z);
PetscErrorCode LaMEMFDSTAG_DOFIndexDestroy(void *dof);
PetscErrorCode LaMEMFDSTAG_DOFIndexCompute(void *dof, PetscInt idxmod);

// Main FDSTAG functions
PetscErrorCode LaMEMFDSTAG_Create(void *fs, void *fb);
PetscErrorCode LaMEMFDSTAG_Destroy(void *fs);
PetscErrorCode LaMEMFDSTAG_ReadRestart(void *fs, FILE *fp);
PetscErrorCode LaMEMFDSTAG_WriteRestart(void *fs, FILE *fp);
PetscErrorCode LaMEMFDSTAG_CreateDMDA(void *fs, PetscInt Nx, PetscInt Ny, PetscInt Nz,
                                      PetscInt Px, PetscInt Py, PetscInt Pz,
                                      PetscInt *lx, PetscInt *ly, PetscInt *lz);
PetscErrorCode LaMEMFDSTAG_GetNeighbProc(void *fs);
PetscErrorCode LaMEMFDSTAG_GetPointRanks(void *fs, PetscScalar *X, PetscInt *lrank, PetscMPIInt *grank);
PetscErrorCode LaMEMFDSTAG_GetAspectRatio(void *fs, PetscScalar *maxAspRat);
PetscErrorCode LaMEMFDSTAG_View(void *fs);
PetscErrorCode LaMEMFDSTAG_GetLocalBox(void *fs, PetscScalar *bx, PetscScalar *by, PetscScalar *bz,
                                       PetscScalar *ex, PetscScalar *ey, PetscScalar *ez);
PetscErrorCode LaMEMFDSTAG_GetGlobalBox(void *fs, PetscScalar *bx, PetscScalar *by, PetscScalar *bz,
                                        PetscScalar *ex, PetscScalar *ey, PetscScalar *ez);
PetscErrorCode LaMEMFDSTAG_SaveGrid(void *fs);

// Wrapper functions
PetscErrorCode LaMEMFDSTAG_DMDACreate3dSetUp(MPI_Comm comm, DMBoundaryType bx, DMBoundaryType by, 
                                             DMBoundaryType bz, DMDAStencilType stencil_type,
                                             PetscInt M, PetscInt N, PetscInt P, PetscInt m, PetscInt n, PetscInt p,
                                             PetscInt dof, PetscInt s, const PetscInt lx[], const PetscInt ly[], 
                                             const PetscInt lz[], DM *da);

#ifdef __cplusplus
}
#endif

#endif
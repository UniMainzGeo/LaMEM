#ifndef LAMEMTOOLS_C_H
#define LAMEMTOOLS_C_H

#include <petsc.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

// Printing and timing functions
void LaMEMTools_PrintStart(PetscLogDouble *t_beg, const char *msg, const char *filename);
void LaMEMTools_PrintDone(PetscLogDouble t_beg);
void LaMEMTools_PrintStep(PetscInt step);

// Vector restart I/O functions
PetscErrorCode LaMEMTools_VecReadRestart(Vec x, FILE *fp);
PetscErrorCode LaMEMTools_VecWriteRestart(Vec x, FILE *fp);

// Basic statistical functions
PetscScalar LaMEMTools_getArthMean(PetscScalar *data, PetscInt n);
PetscScalar LaMEMTools_getVar(PetscScalar *data, PetscInt n);
PetscScalar LaMEMTools_getStdv(PetscScalar *data, PetscInt n);

// Memory allocation functions
PetscErrorCode LaMEMTools_makeMPIIntArray(PetscMPIInt **arr, const PetscMPIInt *init, const PetscInt n);
PetscErrorCode LaMEMTools_clearIntArray(PetscInt *arr, const PetscInt n);
PetscErrorCode LaMEMTools_makeIntArray(PetscInt **arr, const PetscInt *init, const PetscInt n);
PetscErrorCode LaMEMTools_makeScalArray(PetscScalar **arr, const PetscScalar *init, const PetscInt n);

// MPI utility functions
PetscInt LaMEMTools_ISRankZero(MPI_Comm comm);
PetscInt LaMEMTools_ISParallel(MPI_Comm comm);
PetscMPIInt LaMEMTools_getGlobalRank(PetscInt i, PetscInt j, PetscInt k, PetscInt m, PetscInt n, PetscInt p);
PetscMPIInt LaMEMTools_getGlobalRankPeriodic(PetscInt i, PetscInt j, PetscInt k,
                                             PetscInt m, PetscInt n, PetscInt p,
                                             PetscInt pi, PetscInt pj, PetscInt pk);
void LaMEMTools_getLocalRank(PetscInt *i, PetscInt *j, PetscInt *k, PetscMPIInt rank, PetscInt m, PetscInt n);

// Directory operations
PetscErrorCode LaMEMTools_DirMake(const char *name);
PetscErrorCode LaMEMTools_DirRemove(const char *name);
PetscErrorCode LaMEMTools_DirRename(const char *old_name, const char *new_name);
PetscErrorCode LaMEMTools_DirCheck(const char *name, PetscInt *exists);

// Polygon operations
void LaMEMTools_polygon_box(PetscInt *pnv, PetscScalar *vcoord, PetscScalar rtol,
                            PetscScalar *atol, PetscScalar *box);
void LaMEMTools_in_polygon(PetscInt np, PetscScalar *pcoord, PetscInt nv, PetscScalar *vcoord,
                           PetscScalar *box, PetscScalar atol, PetscInt *in);

// Mathematical utility functions
void LaMEMTools_linSpace(PetscScalar min, PetscScalar max, PetscInt N, PetscScalar *outVec);
void LaMEMTools_interpStretch(PetscScalar *Sx, PetscScalar *Sy, PetscInt numCtrlPoly,
                              PetscInt *CtrlPoly, PetscInt numPoly, PetscScalar *SxAll, PetscScalar *SyAll);
void LaMEMTools_findCenterMass(PetscScalar *coords, PetscInt nN, PetscScalar *x_cen, PetscScalar *y_cen);
void LaMEMTools_stretchPolygon(PetscScalar *coords, PetscInt nN, PetscScalar Sx, PetscScalar Sy);

// Indexing functions
PetscInt LaMEMTools_getPtrCnt(PetscInt n, PetscInt counts[], PetscInt ptr[]);
void LaMEMTools_rewindPtr(PetscInt n, PetscInt ptr[]);

// Service functions
PetscErrorCode LaMEMTools_getPhaseRatio(PetscInt n, PetscScalar *v, PetscScalar *rsum);

// Nonlinear equation solver
PetscInt LaMEMTools_solveBisect(PetscScalar a, PetscScalar b, PetscScalar tol, PetscInt maxit,
                                PetscScalar *x, PetscInt *it,
                                PetscScalar (*f)(PetscScalar, void*), void *pctx);

#ifdef __cplusplus
}
#endif

#endif
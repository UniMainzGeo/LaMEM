#include "LaMEMTools_C.h"
#include "../LaMEM.h"
#include "../tools.h"

extern "C" {

// Printing and timing functions
void LaMEMTools_PrintStart(PetscLogDouble *t_beg, const char *msg, const char *filename) {
    PrintStart(t_beg, msg, filename);
}

void LaMEMTools_PrintDone(PetscLogDouble t_beg) {
    PrintDone(t_beg);
}

void LaMEMTools_PrintStep(PetscInt step) {
    PrintStep(step);
}

// Vector restart I/O functions
PetscErrorCode LaMEMTools_VecReadRestart(Vec x, FILE *fp) {
    return VecReadRestart(x, fp);
}

PetscErrorCode LaMEMTools_VecWriteRestart(Vec x, FILE *fp) {
    return VecWriteRestart(x, fp);
}

// Basic statistical functions
PetscScalar LaMEMTools_getArthMean(PetscScalar *data, PetscInt n) {
    return getArthMean(data, n);
}

PetscScalar LaMEMTools_getVar(PetscScalar *data, PetscInt n) {
    return getVar(data, n);
}

PetscScalar LaMEMTools_getStdv(PetscScalar *data, PetscInt n) {
    return getStdv(data, n);
}

// Memory allocation functions
PetscErrorCode LaMEMTools_makeMPIIntArray(PetscMPIInt **arr, const PetscMPIInt *init, const PetscInt n) {
    return makeMPIIntArray(arr, init, n);
}

PetscErrorCode LaMEMTools_clearIntArray(PetscInt *arr, const PetscInt n) {
    return clearIntArray(arr, n);
}

PetscErrorCode LaMEMTools_makeIntArray(PetscInt **arr, const PetscInt *init, const PetscInt n) {
    return makeIntArray(arr, init, n);
}

PetscErrorCode LaMEMTools_makeScalArray(PetscScalar **arr, const PetscScalar *init, const PetscInt n) {
    return makeScalArray(arr, init, n);
}

// MPI utility functions
PetscInt LaMEMTools_ISRankZero(MPI_Comm comm) {
    return ISRankZero(comm);
}

PetscInt LaMEMTools_ISParallel(MPI_Comm comm) {
    return ISParallel(comm);
}

PetscMPIInt LaMEMTools_getGlobalRank(PetscInt i, PetscInt j, PetscInt k, PetscInt m, PetscInt n, PetscInt p) {
    return getGlobalRank(i, j, k, m, n, p);
}

PetscMPIInt LaMEMTools_getGlobalRankPeriodic(PetscInt i, PetscInt j, PetscInt k,
                                             PetscInt m, PetscInt n, PetscInt p,
                                             PetscInt pi, PetscInt pj, PetscInt pk) {
    return getGlobalRankPeriodic(i, j, k, m, n, p, pi, pj, pk);
}

void LaMEMTools_getLocalRank(PetscInt *i, PetscInt *j, PetscInt *k, PetscMPIInt rank, PetscInt m, PetscInt n) {
    getLocalRank(i, j, k, rank, m, n);
}

// Directory operations
PetscErrorCode LaMEMTools_DirMake(const char *name) {
    return DirMake(name);
}

PetscErrorCode LaMEMTools_DirRemove(const char *name) {
    return DirRemove(name);
}

PetscErrorCode LaMEMTools_DirRename(const char *old_name, const char *new_name) {
    return DirRename(old_name, new_name);
}

PetscErrorCode LaMEMTools_DirCheck(const char *name, PetscInt *exists) {
    return DirCheck(name, exists);
}

// Polygon operations
void LaMEMTools_polygon_box(PetscInt *pnv, PetscScalar *vcoord, PetscScalar rtol,
                            PetscScalar *atol, PetscScalar *box) {
    polygon_box(pnv, vcoord, rtol, atol, box);
}

void LaMEMTools_in_polygon(PetscInt np, PetscScalar *pcoord, PetscInt nv, PetscScalar *vcoord,
                           PetscScalar *box, PetscScalar atol, PetscInt *in) {
    in_polygon(np, pcoord, nv, vcoord, box, atol, in);
}

// Mathematical utility functions
void LaMEMTools_linSpace(PetscScalar min, PetscScalar max, PetscInt N, PetscScalar *outVec) {
    linSpace(min, max, N, outVec);
}

void LaMEMTools_interpStretch(PetscScalar *Sx, PetscScalar *Sy, PetscInt numCtrlPoly,
                              PetscInt *CtrlPoly, PetscInt numPoly, PetscScalar *SxAll, PetscScalar *SyAll) {
    interpStretch(Sx, Sy, numCtrlPoly, CtrlPoly, numPoly, SxAll, SyAll);
}

void LaMEMTools_findCenterMass(PetscScalar *coords, PetscInt nN, PetscScalar *x_cen, PetscScalar *y_cen) {
    findCenterMass(coords, nN, *x_cen, *y_cen);
}

void LaMEMTools_stretchPolygon(PetscScalar *coords, PetscInt nN, PetscScalar Sx, PetscScalar Sy) {
    stretchPolygon(coords, nN, Sx, Sy);
}

// Indexing functions
PetscInt LaMEMTools_getPtrCnt(PetscInt n, PetscInt counts[], PetscInt ptr[]) {
    return getPtrCnt(n, counts, ptr);
}

void LaMEMTools_rewindPtr(PetscInt n, PetscInt ptr[]) {
    rewindPtr(n, ptr);
}

// Service functions
PetscErrorCode LaMEMTools_getPhaseRatio(PetscInt n, PetscScalar *v, PetscScalar *rsum) {
    return getPhaseRatio(n, v, rsum);
}

// Nonlinear equation solver
PetscInt LaMEMTools_solveBisect(PetscScalar a, PetscScalar b, PetscScalar tol, PetscInt maxit,
                                PetscScalar *x, PetscInt *it,
                                PetscScalar (*f)(PetscScalar, void*), void *pctx) {
    return solveBisect(a, b, tol, maxit, *x, *it, f, pctx);
}

}
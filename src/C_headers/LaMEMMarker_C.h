#ifndef LAMEMMARKER_C_H
#define LAMEMMARKER_C_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

// Thin wrappers around actual marker functions from the source
PetscErrorCode LaMEMMarker_ADVMarkInit(void *actx, void *fb);
PetscErrorCode LaMEMMarker_ADVMarkInitCoord(void *actx);
PetscErrorCode LaMEMMarker_ADVMarkPerturb(void *actx);
PetscErrorCode LaMEMMarker_ADVMarkSave(void *actx);
PetscErrorCode LaMEMMarker_ADVMarkCheckMarkers(void *actx);
PetscErrorCode LaMEMMarker_ADVMarkSetTempGrad(void *actx);
PetscErrorCode LaMEMMarker_ADVMarkSetTempPhase(void *actx);
PetscErrorCode LaMEMMarker_ADVMarkSetTempFile(void *actx, void *fb);
PetscErrorCode LaMEMMarker_ADVMarkSetTempVector(void *actx);

// Specific initialization routines
PetscErrorCode LaMEMMarker_ADVMarkInitFiles(void *actx, void *fb);
PetscErrorCode LaMEMMarker_ADVMarkInitGeom(void *actx, void *fb);
PetscErrorCode LaMEMMarker_ADVMarkInitPolygons(void *actx, void *fb);
PetscErrorCode LaMEMMarker_ADVMarkReadCtrlPoly(void *fb, void *CtrlPoly, PetscInt *VolID, PetscInt *nCP);

// Phase diagram functions
PetscErrorCode LaMEMMarker_LoadPhaseDiagram(void *actx, void *phases, PetscInt i);

// Geometric primitives functions (these are typically void but we wrap them)
void LaMEMMarker_SetPhaseSphere(void *sphere, void *P);
void LaMEMMarker_SetPhaseEllipsoid(void *ellipsoid, void *P);
void LaMEMMarker_SetPhaseBox(void *box, void *P);
void LaMEMMarker_SetPhaseRidge(void *ridge, void *P);
void LaMEMMarker_SetPhaseLayer(void *layer, void *P);
void LaMEMMarker_SetPhaseHex(void *hex, void *P);
void LaMEMMarker_SetPhaseCylinder(void *cylinder, void *P);

// Temperature computation
void LaMEMMarker_ComputeTemperature(void *geom, void *P, PetscScalar *T);

// Utility functions
void LaMEMMarker_ADVMarkSecIdx(void *actx, PetscInt dir, PetscInt Islice, PetscInt *idx);
void LaMEMMarker_HexGetBoundingBox(PetscScalar *coord, PetscScalar *bounds);
PetscInt LaMEMMarker_TetPointTest(PetscScalar *coord, PetscInt *ii, PetscScalar *xp, PetscScalar tol);

#ifdef __cplusplus
}
#endif

#endif
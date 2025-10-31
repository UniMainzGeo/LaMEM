#ifndef LAMEMADJOINT_C_H
#define LAMEMADJOINT_C_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

// Main adjoint functions
PetscErrorCode LaMEMAdjoint_Main(void *IOparam);
PetscErrorCode LaMEMAdjoint_OptimizationTAO(void *tao, void *P, PetscReal *F, void *grad, void *ctx);
PetscErrorCode LaMEMAdjoint_Optimization(void *P, PetscScalar F, void *grad, void *ctx);

// Adjoint gradient object creation/destruction
PetscErrorCode LaMEMAdjoint_Create(void *aop, void *jr, void *IOparam);
PetscErrorCode LaMEMAdjoint_Destroy(void *aop, void *IOparam);

// Adjoint vectors creation/destruction
PetscErrorCode LaMEMAdjoint_VectorsCreate(void *Adjoint_vectors, void *IOparam);
PetscErrorCode LaMEMAdjoint_VectorsDestroy(void *Adjoint_vectors, void *IOparam);

// Core gradient computation functions
PetscErrorCode LaMEMAdjoint_ObjectiveAndGradientFunction(void *aop, void *jr, void *nl, void *IOparam, void *snes, void *surf);
PetscErrorCode LaMEMAdjoint_ComputeGradientsAndObjectiveFunction(void *Parameters, PetscScalar *ObjectiveValue, void *Gradient, void *IOparam);
PetscErrorCode LaMEMAdjoint_ComputeGradients(void *jr, void *aop, void *nl, void *snes, void *IOparam);

// Objective function computation
PetscErrorCode LaMEMAdjoint_ObjectiveFunction(void *aop, void *jr, void *IOparam, void *surf);

// Finite difference gradients
PetscErrorCode LaMEMAdjoint_FiniteDifferenceGradients(void *IOparam);

// Input/Output functions
PetscErrorCode LaMEMAdjoint_ReadInputSetDefaults(void *IOparam, void *Adjoint_Vectors);
PetscErrorCode LaMEMAdjoint_PrintGradientsAndObservationPoints(void *IOparam);
PetscErrorCode LaMEMAdjoint_PrintCostFunction(void *IOparam);
PetscErrorCode LaMEMAdjoint_PrintScalingLaws(void *IO_param);

// Material parameter scanning and modification
PetscErrorCode LaMEMAdjoint_ScanForMaterialParameters(void *fb, void *scaling, PetscInt *numParam, char parameterType, PetscInt *ParamID, PetscScalar *ParamValue, PetscInt *PhaseID, PetscScalar *Bounds);
PetscErrorCode LaMEMAdjoint_AddMaterialParameterToCommandLineOptions(char *name, PetscInt ID, PetscScalar val);
PetscErrorCode LaMEMAdjoint_DeleteMaterialParameterFromCommandLineOptions(char *name, PetscInt ID);
PetscErrorCode LaMEMAdjoint_CreateModifiedMaterialDatabase(void *IOparam);
PetscErrorCode LaMEMAdjoint_CopyParameterToLaMEMCommandLine(void *IOparam, PetscScalar CurVal, PetscInt j);

// Boundary conditions
PetscErrorCode LaMEMAdjoint_ApplyBCs(void *dF, void *bc);

// Projection and interpolation
PetscErrorCode LaMEMAdjoint_PointInPro(void *jr, void *aop, void *IOparam, void *surf);

// Finite difference specific functions
PetscErrorCode LaMEMAdjoint_Get_F_dFdu_Center(void *jr, void *aop, void *IOparam);
PetscErrorCode LaMEMAdjoint_GradientResetParameter(void *nl, PetscInt CurPar, PetscInt CurPhase, void *aop);
PetscErrorCode LaMEMAdjoint_FormResidualFieldFD(void *snes, void *x, void *psi, void *nl, void *aop, void *IOparam);

// Constitutive equation finite difference functions
PetscErrorCode LaMEMAdjoint_devConstEqFD(void *ctx, void *aop, void *IOparam, PetscInt ii, PetscInt jj, PetscInt k, PetscInt ik, PetscInt jk, PetscInt kk);
PetscErrorCode LaMEMAdjoint_cellConstEqFD(void *ctx, void *svCell, PetscScalar dxx, PetscScalar dyy, PetscScalar dzz, PetscScalar *sxx, PetscScalar *syy, PetscScalar *szz, PetscScalar *gres, PetscScalar *rho, void *aop, void *IOparam, PetscInt ii, PetscInt jj, PetscInt k, PetscInt ik, PetscInt jk, PetscInt kk);
PetscErrorCode LaMEMAdjoint_setUpPhaseFD(void *ctx, PetscInt ID, void *aop, void *IOparam, PetscInt ii, PetscInt jj, PetscInt k, PetscInt ik, PetscInt jk, PetscInt kk);
PetscErrorCode LaMEMAdjoint_edgeConstEqFD(void *ctx, void *svEdge, PetscScalar d, PetscScalar *s, void *aop, void *IOparam, PetscInt ii, PetscInt jj, PetscInt k, PetscInt ik, PetscInt jk, PetscInt kk);

// Utility functions
PetscErrorCode LaMEMAdjoint_Parameter_SetFDgrad_Option(PetscInt *FD_grad, char *name);
PetscErrorCode LaMEMAdjoint_swapStruct(void *A, void *B);

#ifdef __cplusplus
}
#endif

#endif
#ifndef LAMEMMELTPARAM_C_H
#define LAMEMMELTPARAM_C_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

// C wrapper for meltPar_Katz structure
typedef struct {
    PetscScalar A1, A2, A3;          // Solidus coefficients
    PetscScalar B1, B2, B3;          // Lherzolite liquidus coefficients
    PetscScalar C1, C2, C3;          // Clinopyroxene-out liquidus coefficients
    PetscScalar r1, r2;              // Cpx reaction coefficients
    PetscScalar beta1, beta2;        // Melting exponents
    PetscScalar K, gamma;            // Water depression coefficients
    PetscScalar D_water, lambda;     // Water partitioning and saturation
    PetscScalar chi1, chi2;          // Water saturation coefficients
    PetscScalar Cp, DS;              // Heat capacity and entropy of fusion
} LaMEMMeltPar_Katz;

// Thin wrappers around actual melt parameterization functions
void LaMEMMeltParam_SetMeltParamsToDefault_Katz(LaMEMMeltPar_Katz *mp);
PetscScalar LaMEMMeltParam_MPgetFReactive(PetscScalar P, PetscScalar T, PetscScalar Cf, PetscScalar M, LaMEMMeltPar_Katz *mp);
PetscScalar LaMEMMeltParam_MPgetFEquilib(PetscScalar P, PetscScalar T, PetscScalar X, PetscScalar M, LaMEMMeltPar_Katz *mp);
PetscScalar LaMEMMeltParam_MPgetTSolidus(PetscScalar P, PetscScalar X, LaMEMMeltPar_Katz *mp);
PetscScalar LaMEMMeltParam_MPgetTEquilib(PetscScalar P, PetscScalar F, PetscScalar X, PetscScalar M, LaMEMMeltPar_Katz *mp);
PetscScalar LaMEMMeltParam_MPgetFconsH(PetscScalar P, PetscScalar Ti, PetscScalar X, PetscScalar M, PetscScalar *Tf, LaMEMMeltPar_Katz *mp);

// Helper functions
PetscScalar LaMEMMeltParam_CalcDT(PetscScalar P, PetscScalar X, PetscScalar F, LaMEMMeltPar_Katz *mp);
PetscScalar LaMEMMeltParam_CalcF(PetscScalar T, PetscScalar dT, PetscScalar P, PetscScalar Fcpx, LaMEMMeltPar_Katz *mp);
PetscScalar LaMEMMeltParam_FZero(PetscScalar F, PetscScalar T, PetscScalar P, PetscScalar X, PetscScalar Fcpx, LaMEMMeltPar_Katz *mp);
PetscScalar LaMEMMeltParam_HZero(PetscScalar F, PetscScalar T, PetscScalar P, PetscScalar X, PetscScalar M, LaMEMMeltPar_Katz *mp);
PetscScalar LaMEMMeltParam_FX_bal(PetscScalar x1, PetscScalar x2, PetscScalar T, PetscScalar P, PetscScalar X, PetscScalar Fcpx, LaMEMMeltPar_Katz *mp);
PetscScalar LaMEMMeltParam_FT_bal(PetscScalar x1, PetscScalar x2, PetscScalar T, PetscScalar P, PetscScalar X, PetscScalar M, LaMEMMeltPar_Katz *mp);

// Global parameter access
void LaMEMMeltParam_SetPc(PetscScalar pc_value);
PetscScalar LaMEMMeltParam_GetPc(void);

#ifdef __cplusplus
}
#endif

#endif
/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : LaMEM
 **   License      : MIT, see LICENSE file for details
 **   Contributors : Anton Popov, Boris Kaus, see AUTHORS file for complete list
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : kaus@uni-mainz.de, popov@uni-mainz.de
 **
 **   The current framework is developed by Georg Reuber (JGU Mainz)
 **
 **   If you think it is helpful, please cite the following paper:
 **   Georg S. Reuber, Anton A. Popov, Boris J.P. Kaus, (2018)
 **   Deriving scaling laws in geodynamics using adjoint gradients,
 **   Tectonophysics, Vol. 746, p. 352-363. doi:10.1016/j.tecto.2017.07.017
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
//---------------------------------------------------------------------------
// COMPUTATION OF ADJOINT GRADIENTS
//---------------------------------------------------------------------------
#ifndef __adjoint_h__
#define __adjoint_h__
//---------------------------------------------------------------------------

struct Scaling;
struct FDSTAG;
struct FreeSurf;
struct DBMat;
struct Tensor2RN;
struct PData;
struct Material_t;
struct JacRes;
struct Controls;
struct NLSol;
struct ModParam;
struct FB;
struct BCCtx;

//---------------------------------------------------------------------------

// Structure that holds parameters for the adjoint gradient computation
struct AdjGrad
{
	PetscScalar      FD_epsilon;              // Epsilon, employed for finite difference calculation of dres/dp
	PetscScalar      CurScal;
	Vec              dF, dPardu;
	Vec              pro;
	Vec              vx, vy, vz, sty;
};

//---------------------------------------------------------------------------

// Structure that holds vectors required by TAO
struct Adjoint_Vecs
{
	Vec             Ub, Lb, grad, P;
};

//---------------------------------------------------------------------------

// Adjoint optimization driving routines
PetscErrorCode AdjointOptimisation(Vec P, PetscScalar F, Vec grad, void *ctx);
PetscErrorCode AdjointOptimisationTAO(Tao tao, Vec P, PetscReal *F, Vec grad, void *ctx);
PetscErrorCode LaMEMAdjointMain(ModParam *IOparam);

// Initialize & read input
PetscErrorCode LaMEMAdjointReadInputSetDefaults(ModParam *IOparam, Adjoint_Vecs *Adjoint_Vectors);
PetscErrorCode Adjoint_ScanForMaterialParameters(FB*, Scaling*, PetscInt*, char, PetscInt*,PetscScalar*,PetscInt*, PetscScalar*);


// Create & Destroy Adjoint_Vectors object
PetscErrorCode AdjointVectorsCreate(Adjoint_Vecs *Adjoint_vectors, ModParam *IOparam);
PetscErrorCode AdjointVectorsDestroy(Adjoint_Vecs *Adjoint_vectors, ModParam *IOparam);

// Compute the gradients for the adjoint inversion
PetscErrorCode AdjointObjectiveAndGradientFunction(AdjGrad *aop, ModParam *IOparam, SNES snes);

// Compute the gradients for the adjoint inversion
PetscErrorCode ComputeGradientsAndObjectiveFunction(Vec Parameters, PetscScalar *ObjectiveValue, Vec Gradient, ModParam *IOparam);

// Adjoint Gradients
PetscErrorCode AdjointComputeGradients(JacRes *jr, AdjGrad *aop, NLSol *nl, SNES snes, ModParam *IOparam);
PetscErrorCode Adjoint_ApplyBCs(Vec dF, BCCtx* bc);

// Cost function
PetscErrorCode AdjointObjectiveFunction(AdjGrad *aop, JacRes *jr, ModParam *IOparam, FreeSurf *surf);

// 'Brute-force' finite difference gradients
PetscErrorCode AdjointFiniteDifferenceGradients(ModParam *IOparam);
PetscErrorCode PrintGradientsAndObservationPoints(ModParam *IOparam);
PetscErrorCode PrintCostFunction(ModParam *IOparam);

// Interpolate the adjoint points and include them into the projection vector
PetscErrorCode AdjointPointInPro(JacRes *jr, AdjGrad *aop, ModParam *IOparam, FreeSurf *surf);

// PSD calculations
PetscErrorCode AdjointGet_F_dFdu_Center(JacRes *jr, AdjGrad *aop, ModParam *IOparam);

// Add or remove parameters from command-line database & update material DB
PetscErrorCode AddMaterialParameterToCommandLineOptions(char *name, PetscInt ID, PetscScalar val);
PetscErrorCode DeleteMaterialParameterFromCommandLineOptions(char *name, PetscInt ID);
PetscErrorCode CreateModifiedMaterialDatabase(ModParam *IOparam);
PetscErrorCode CopyParameterToLaMEMCommandLine(ModParam *IOparam, PetscScalar CurVal, PetscInt j);

PetscErrorCode Parameter_SetFDgrad_Option(PetscInt *FD_grad, char *name);

// Print scaling laws
PetscErrorCode PrintScalingLaws(ModParam *IO_param);

// Create & Destroy aop object
PetscErrorCode AdjointCreate(AdjGrad *aop, JacRes *jr, ModParam *IOparam);
PetscErrorCode AdjointDestroy(AdjGrad *aop, ModParam *IOparam);

// Helper functions
PetscErrorCode swapStruct(struct Material_t *A, struct Material_t *B);

//---------------------------------------------------------------------------

#endif

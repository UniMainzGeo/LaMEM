//---------------------------------------------------------------------------
//.........   LaMEM - FDSTAG CANONICAL INTERFACE CHECKING ROUTINES   ........
//---------------------------------------------------------------------------
#ifndef __check_fdstag_h__
#define __check_fdstag_h__
//---------------------------------------------------------------------------
/*
PetscErrorCode DoDarcyTests(NLCtx *nlctx, UserContext *user);

PetscErrorCode DarcyPostProcess(NLCtx *nlctx, UserContext *user);

PetscErrorCode DoMGTests(NLCtx *nlctx, PVOut *pvout);

PetscErrorCode BlockSetUpMGVelBlock(NLCtx *nlctx, MGCtx *mg);

PetscErrorCode FDSTAGOutputResTest(FDSTAG *fs, JacResCtx *jrctx, Vec f);

PetscErrorCode GetLinRes(Mat A, Vec x, Vec rhs, Vec res);

PetscErrorCode FieldSplitTest(NLCtx *nlctx, PVOut *pvout, Vec InitGuess, PetscScalar time, PetscInt itime);

PetscErrorCode MGTest(NLCtx *nlctx, PVOut *pvout, Vec InitGuess, PetscScalar time, PetscInt itime);

//---------------------------------------------------------------------------
typedef struct
{
	PetscScalar a;
	PetscScalar b;

}MaxDiff;
//---------------------------------------------------------------------------

typedef struct
{
	PetscScalar A[8];
	PetscScalar cxb, cxe;
	PetscScalar cyb, cye;
	PetscScalar czb, cze;

} BCValues;

PetscScalar InterpolateLinear3D(PetscScalar cx, PetscScalar cy, PetscScalar cz,  BCValues BC);

// Initialize velocity vectors to test strain rates.
// Selects the velocity direction x, y, z (0, 1, 2) with vectDir,
// and applies constant gradient in the gradDir direction
// between the values begVal & endVal (in the positive direction of gradDir).
// Initializes boundary ghost points in the tangential directions accordingly.

PetscErrorCode InitVelocityTest(
	JacRes      *jr,
	UserContext *usr,
	PetscInt     vectDir,
	PetscInt     gradDir,
	PetscScalar  begVal,
	PetscScalar  endVal);

PetscErrorCode JacResCtxClearVelocity(JacRes *jr, PetscInt vectDir);

PetscErrorCode StrainRateSingleComp(
	JacRes      *jr,
	UserContext *usr,
	PVOut       *pvout,
	PetscInt     vectDir,
	PetscInt     gradDir,
	PetscScalar  ttime,
	PetscInt     itime);

PetscErrorCode StrainRateInterpTest(
	JacRes      *jr,
	UserContext *usr,
	PVOut       *pvout);
*/
//---------------------------------------------------------------------------
#endif

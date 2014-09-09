//---------------------------------------------------------------------------
//.........   LaMEM - FDSTAG CANONICAL INTERFACE CHECKING ROUTINES   ........
//---------------------------------------------------------------------------
#ifndef __check_fdstag_h__
#define __check_fdstag_h__
//---------------------------------------------------------------------------

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
	FDSTAG      *fs,
	JacResCtx   *jrctx,
	UserContext *usr,
	PetscInt     vectDir,
	PetscInt     gradDir,
	PetscScalar  begVal,
	PetscScalar  endVal);

PetscErrorCode JacResCtxClearVelocity(JacResCtx *jrctx, PetscInt vectDir);

PetscErrorCode StrainRateSingleComp(
	FDSTAG      *fs,
	JacResCtx   *jrctx,
	UserContext *usr,
	PVOut       *pvout,
	PetscInt     vectDir,
	PetscInt     gradDir,
	PetscScalar  ttime,
	PetscInt     itime);

PetscErrorCode StrainRateInterpTest(
	FDSTAG      *fs,
	JacResCtx   *jrctx,
	UserContext *usr,
	PVOut       *pvout);

//---------------------------------------------------------------------------

// compute scalar residual magnitude vector
PetscErrorCode FDSTAGetResMag(FDSTAG *fs, UserContext *usr, Vec f, Vec fmag);

// check phase ratios in control volumes
//PetscErrorCode FDSTAGCheckPhaseRatios(FDSTAG *fs, JacResCtx *jrctx);

// compute sum of array elements
PetscScalar ArraySum(PetscScalar *a, PetscInt n);

// dump pahse ratios to disk
//PetscErrorCode FDSTAGDumpPhaseRatios(FDSTAG *fs, JacResCtx *jrctx);

// copy array contents
void ArrayCopy(PetscScalar *a, PetscScalar *b, PetscInt n);

// compare phase ratios with reference version
PetscErrorCode FDSTAGComparePhaseRatios(FDSTAG *fs, JacResCtx *jrctx, PetscScalar rtol, PetscScalar atol);

// access phase ratio array
PetscScalar * getPtr(PetscInt nx, PetscInt ny, PetscInt i, PetscInt j, PetscInt k, PetscInt ndof, PetscScalar *a);

// compare array contents
PetscErrorCode ArrayCompare(PetscScalar *a, PetscScalar *b, PetscInt n, PetscScalar rtol, PetscScalar atol);

PetscErrorCode FDSTAGetResidualTest(FDSTAG *fs, JacResCtx *jrctx);

PetscErrorCode InitVec(DM da, Vec v);

PetscErrorCode DumpVec(DM da, Vec v);

PetscErrorCode CountVec(DM da, Vec v, PetscScalar rtol, PetscScalar atol);

void MDiff(PetscScalar *a, PetscScalar *b, MaxDiff *md, PetscInt n);

PetscInt ArrayCount(PetscScalar *a, PetscScalar *b, PetscInt n, PetscScalar rtol, PetscScalar atol);

PetscErrorCode TestFuckingVelocities(FDSTAG *fs, JacResCtx *jrctx);

PetscErrorCode FDSTAGDumpVelocities(FDSTAG *fs, JacResCtx *jrctx);

PetscErrorCode FDSTAGCountVelocities(FDSTAG *fs, JacResCtx *jrctx, PetscScalar rtol, PetscScalar atol);

PetscErrorCode FDSTAGSetVelocities(FDSTAG *fs, JacResCtx *jrctx);

//---------------------------------------------------------------------------
#endif

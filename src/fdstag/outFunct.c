//---------------------------------------------------------------------------
//....................   FDSTAG VECTOR OUTPUT ROUTINES   ....................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "paraViewOutBin.h"
#include "outFunct.h"
#include "interpolate.h"
//---------------------------------------------------------------------------
// WARNING!
//
// ParaView symmetric tensor components ordering is: xx, yy, zz, xy, yz, xz
//
// This is diagonal (rather than row-wise) storage format !!!
//
// As usual, this isn't documented anywhere !!! Take care of this in future versions.
//---------------------------------------------------------------------------
// interpolation function header
#define COPY_FUNCTION_HEADER \
	FDSTAG      *fs; \
	Scaling     *scal; \
	PetscScalar ***buff, cf; \
	PetscInt     i, j, k, nx, ny, nz, sx, sy, sz, iter; \
	InterpFlags  iflag; \
	PetscErrorCode ierr; \
	PetscFunctionBegin; \
	fs   = outbuf->fs; \
	scal = &jr->scal; \
	iflag.update    = PETSC_FALSE; \
	iflag.use_bound = PETSC_FALSE;
//---------------------------------------------------------------------------
// access function header
#define ACCESS_FUNCTION_HEADER \
	PetscScalar cf; \
	Scaling     *scal; \
	InterpFlags  iflag; \
	PetscErrorCode ierr; \
	PetscFunctionBegin; \
	scal = &jr->scal; \
	iflag.update    = PETSC_FALSE; \
	iflag.use_bound = PETSC_FALSE;
//---------------------------------------------------------------------------
#define INTERPOLATE_COPY(da, vec, IFUNCT, FIELD, ncomp, dir) \
	ierr = DMDAGetCorners (da, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr); \
	ierr = DMDAVecGetArray(da, vec, &buff); CHKERRQ(ierr); \
	iter = 0; \
	START_STD_LOOP \
		FIELD \
	END_STD_LOOP \
	ierr = DMDAVecRestoreArray(da, vec, &buff); CHKERRQ(ierr); \
	LOCAL_TO_LOCAL(da, vec) \
	ierr = IFUNCT(fs, vec, outbuf->lbcor, iflag); CHKERRQ(ierr); \
	if(iflag.update != PETSC_TRUE) \
	{	ierr = OutBufPut3DVecComp(outbuf, ncomp, dir, cf, 0.0); CHKERRQ(ierr); }
//---------------------------------------------------------------------------
#define INTERPOLATE_ACCESS(vec, IFUNCT, ncomp, dir, shift) \
	ierr = IFUNCT(outbuf->fs, vec, outbuf->lbcor, iflag); CHKERRQ(ierr); \
	ierr = OutBufPut3DVecComp(outbuf, ncomp, dir, cf, shift); CHKERRQ(ierr);
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWritePhase"
PetscErrorCode PVOutWritePhase(JacRes *jr, OutBuf *outbuf)
{
	Material_t  *phases;
	PetscScalar *phRat, mID;
	PetscInt     jj, numPhases;

	COPY_FUNCTION_HEADER

	// macro to copy phase parameter to buffer
	#define GET_PHASE \
		phRat = jr->svCell[iter++].phRat; \
		mID = 0.0; \
		for(jj = 0; jj < numPhases; jj++) \
			mID += phRat[jj]*phases[jj].ID; \
		buff[k][j][i] = mID;

	// no scaling is necessary for the phase
	cf = scal->unit;

	// access material parameters
	phases    = jr->phases;
	numPhases = jr->numPhases;

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_PHASE, 1, 0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteDensity"
PetscErrorCode PVOutWriteDensity(JacRes *jr, OutBuf *outbuf)
{
	COPY_FUNCTION_HEADER

	// macro to copy density to buffer
	#define GET_DENSITY buff[k][j][i] = jr->svCell[iter++].svBulk.rho;

	cf = scal->density;

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_DENSITY, 1, 0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteViscosity"
PetscErrorCode PVOutWriteViscosity(JacRes *jr, OutBuf *outbuf)
{
	COPY_FUNCTION_HEADER

	// macro to copy viscosity to buffer
	#define GET_VISCOSITY buff[k][j][i] = jr->svCell[iter++].svDev.eta;

	// output viscosity logarithm in GEO-mode
	// (negative scaling requests logarithmic output)
	if(scal->utype == _GEO_) cf = -scal->viscosity;
	else                     cf =  scal->viscosity;

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_VISCOSITY, 1, 0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteVelocity"
PetscErrorCode PVOutWriteVelocity(JacRes *jr, OutBuf *outbuf)
{
	ACCESS_FUNCTION_HEADER

	cf = scal->velocity;
	iflag.use_bound = PETSC_TRUE;

	INTERPOLATE_ACCESS(jr->lvx, InterpXFaceCorner, 3, 0, 0.0)
	INTERPOLATE_ACCESS(jr->lvy, InterpYFaceCorner, 3, 1, 0.0)
	INTERPOLATE_ACCESS(jr->lvz, InterpZFaceCorner, 3, 2, 0.0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWritePressure"
PetscErrorCode PVOutWritePressure(JacRes *jr, OutBuf *outbuf)
{
	ACCESS_FUNCTION_HEADER

	cf = scal->stress;
	iflag.use_bound = PETSC_TRUE;

	INTERPOLATE_ACCESS(jr->lp, InterpCenterCorner, 1, 0, 0.0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteTemperature"
PetscErrorCode PVOutWriteTemperature(JacRes *jr, OutBuf *outbuf)
{
	ACCESS_FUNCTION_HEADER

	cf = scal->temperature;
	iflag.use_bound = PETSC_TRUE;

	INTERPOLATE_ACCESS(jr->lT, InterpCenterCorner, 1, 0, scal->Tshift)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteDevStress"
PetscErrorCode PVOutWriteDevStress(JacRes *jr, OutBuf *outbuf)
{
	// NOTE! See warning about component ordering scheme above

	COPY_FUNCTION_HEADER

	// macro to copy deviatoric stress components to buffer
	#define GET_SXX buff[k][j][i] = jr->svCell[iter++].sxx;
	#define GET_SYY buff[k][j][i] = jr->svCell[iter++].syy;
	#define GET_SZZ buff[k][j][i] = jr->svCell[iter++].szz;
	#define GET_SXY buff[k][j][i] = jr->svXYEdge[iter++].s;
	#define GET_SYZ buff[k][j][i] = jr->svYZEdge[iter++].s;
	#define GET_SXZ buff[k][j][i] = jr->svXZEdge[iter++].s;

	cf = scal->stress;

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_SXX, 6, 0)
	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_SYY, 6, 1)
	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_SZZ, 6, 2)
	INTERPOLATE_COPY(fs->DA_XY,  outbuf->lbxy,  InterpXYEdgeCorner, GET_SXY, 6, 3)
	INTERPOLATE_COPY(fs->DA_YZ,  outbuf->lbyz,  InterpYZEdgeCorner, GET_SYZ, 6, 4)
	INTERPOLATE_COPY(fs->DA_XZ,  outbuf->lbxz,  InterpXZEdgeCorner, GET_SXZ, 6, 5)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteJ2DevStress"
PetscErrorCode PVOutWriteJ2DevStress(JacRes *jr, OutBuf *outbuf)
{
	SolVarCell *svCell;
	PetscScalar s, J2;

	COPY_FUNCTION_HEADER

	// macros to copy deviatoric strain rate invariant to buffer
	#define GET_J2_STRESS_CENTER \
		svCell = &jr->svCell[iter++]; \
		s = svCell->sxx; J2  = s*s; \
		s = svCell->syy; J2 += s*s; \
		s = svCell->szz; J2 += s*s; \
		buff[k][j][i] = 0.5*J2;

	#define GET_J2_STRESS_XY_EDGE s = jr->svXYEdge[iter++].s; buff[k][j][i] = s*s;
	#define GET_J2_STRESS_YZ_EDGE s = jr->svYZEdge[iter++].s; buff[k][j][i] = s*s;
	#define GET_J2_STRESS_XZ_EDGE s = jr->svXZEdge[iter++].s; buff[k][j][i] = s*s;

	cf = scal->stress;

	iflag.update = PETSC_TRUE;

	ierr = VecSet(outbuf->lbcor, 0.0);

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_J2_STRESS_CENTER,  1, 0)
	INTERPOLATE_COPY(fs->DA_XY,  outbuf->lbxy,  InterpXYEdgeCorner, GET_J2_STRESS_XY_EDGE, 1, 0)
	INTERPOLATE_COPY(fs->DA_YZ,  outbuf->lbyz,  InterpYZEdgeCorner, GET_J2_STRESS_YZ_EDGE, 1, 0)
	INTERPOLATE_COPY(fs->DA_XZ,  outbuf->lbxz,  InterpXZEdgeCorner, GET_J2_STRESS_XZ_EDGE, 1, 0)

	// compute & store second invariant
	ierr = VecSqrtAbs(outbuf->lbcor); CHKERRQ(ierr);

	ierr = OutBufPut3DVecComp(outbuf, 1, 0, cf, 0.0); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteStrainRate"
PetscErrorCode PVOutWriteStrainRate(JacRes *jr, OutBuf *outbuf)
{
	// NOTE! See warning about component ordering scheme above

	COPY_FUNCTION_HEADER

	// macro to copy deviatoric strain rate components to buffer
	#define GET_DXX buff[k][j][i] = jr->svCell[iter++].dxx;
	#define GET_DYY buff[k][j][i] = jr->svCell[iter++].dyy;
	#define GET_DZZ buff[k][j][i] = jr->svCell[iter++].dzz;
	#define GET_DXY buff[k][j][i] = jr->svXYEdge[iter++].d;
	#define GET_DYZ buff[k][j][i] = jr->svYZEdge[iter++].d;
	#define GET_DXZ buff[k][j][i] = jr->svXZEdge[iter++].d;

	cf = scal->strain_rate;

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_DXX, 6, 0)
	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_DYY, 6, 1)
	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_DZZ, 6, 2)
	INTERPOLATE_COPY(fs->DA_XY,  outbuf->lbxy,  InterpXYEdgeCorner, GET_DXY, 6, 3)
	INTERPOLATE_COPY(fs->DA_YZ,  outbuf->lbyz,  InterpYZEdgeCorner, GET_DYZ, 6, 4)
	INTERPOLATE_COPY(fs->DA_XZ,  outbuf->lbxz,  InterpXZEdgeCorner, GET_DXZ, 6, 5)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteJ2StrainRate"
PetscErrorCode PVOutWriteJ2StrainRate(JacRes *jr, OutBuf *outbuf)
{
	SolVarCell *svCell;
	PetscScalar d, J2;

	COPY_FUNCTION_HEADER

	// macros to copy deviatoric strain rate invariant to buffer
	#define GET_J2_STRAIN_RATE_CENTER \
		svCell = &jr->svCell[iter++]; \
		d = svCell->dxx; J2  = d*d; \
		d = svCell->dyy; J2 += d*d; \
		d = svCell->dzz; J2 += d*d; \
		buff[k][j][i] = 0.5*J2;

	#define GET_J2_STRAIN_RATE_XY_EDGE d = jr->svXYEdge[iter++].d; buff[k][j][i] = d*d;
	#define GET_J2_STRAIN_RATE_YZ_EDGE d = jr->svYZEdge[iter++].d; buff[k][j][i] = d*d;
	#define GET_J2_STRAIN_RATE_XZ_EDGE d = jr->svXZEdge[iter++].d; buff[k][j][i] = d*d;

	cf = scal->strain_rate;

	iflag.update = PETSC_TRUE;

	ierr = VecSet(outbuf->lbcor, 0.0);

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_J2_STRAIN_RATE_CENTER,  1, 0)
	INTERPOLATE_COPY(fs->DA_XY,  outbuf->lbxy,  InterpXYEdgeCorner, GET_J2_STRAIN_RATE_XY_EDGE, 1, 0)
	INTERPOLATE_COPY(fs->DA_YZ,  outbuf->lbyz,  InterpYZEdgeCorner, GET_J2_STRAIN_RATE_YZ_EDGE, 1, 0)
	INTERPOLATE_COPY(fs->DA_XZ,  outbuf->lbxz,  InterpXZEdgeCorner, GET_J2_STRAIN_RATE_XZ_EDGE, 1, 0)

	// compute & store second invariant
	ierr = VecSqrtAbs(outbuf->lbcor); CHKERRQ(ierr);

	ierr = OutBufPut3DVecComp(outbuf, 1, 0, cf, 0.0); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteVolRate"
PetscErrorCode PVOutWriteVolRate(JacRes *jr, OutBuf *outbuf)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = 0; CHKERRQ(ierr);
	if(jr)     jr = NULL;
	if(outbuf) outbuf = NULL;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteVorticity"
PetscErrorCode PVOutWriteVorticity(JacRes *jr, OutBuf *outbuf)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = 0; CHKERRQ(ierr);
	if(jr)  jr = NULL;
	if(outbuf) outbuf = NULL;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteAngVelMag"
PetscErrorCode PVOutWriteAngVelMag(JacRes *jr, OutBuf *outbuf)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = 0; CHKERRQ(ierr);
	if(jr)  jr = NULL;
	if(outbuf) outbuf = NULL;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteTotStrain"
PetscErrorCode PVOutWriteTotStrain(JacRes *jr, OutBuf *outbuf)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = 0; CHKERRQ(ierr);
	if(jr)  jr = NULL;
	if(outbuf) outbuf = NULL;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWritePlastStrain"
PetscErrorCode PVOutWritePlastStrain(JacRes *jr, OutBuf *outbuf)
{
	COPY_FUNCTION_HEADER

	// macro to copy accumulated plastic strain (APS) to buffer
	#define GET_APS buff[k][j][i] = jr->svCell[iter++].svDev.APS;

	cf = scal->unit;

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_APS, 1, 0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWritePlastDissip"
PetscErrorCode PVOutWritePlastDissip(JacRes *jr, OutBuf *outbuf)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = 0; CHKERRQ(ierr);
	if(jr)  jr = NULL;
	if(outbuf) outbuf = NULL;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteTotDispl"
PetscErrorCode PVOutWriteTotDispl(JacRes *jr, OutBuf *outbuf)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = 0; CHKERRQ(ierr);
	if(jr)  jr = NULL;
	if(outbuf) outbuf = NULL;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// DEBUG VECTORS
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteMomentRes"
PetscErrorCode PVOutWriteMomentRes(JacRes *jr, OutBuf *outbuf)
{
	ACCESS_FUNCTION_HEADER

	cf = scal->volumetric_force;

	ierr = JacResCopyMomentumRes(jr, jr->gres); CHKERRQ(ierr);

	GLOBAL_TO_LOCAL(outbuf->fs->DA_X, jr->gfx, jr->lfx)
	GLOBAL_TO_LOCAL(outbuf->fs->DA_Y, jr->gfy, jr->lfy)
	GLOBAL_TO_LOCAL(outbuf->fs->DA_Z, jr->gfz, jr->lfz)

	INTERPOLATE_ACCESS(jr->lfx, InterpXFaceCorner, 3, 0, 0.0)
	INTERPOLATE_ACCESS(jr->lfy, InterpYFaceCorner, 3, 1, 0.0)
	INTERPOLATE_ACCESS(jr->lfz, InterpZFaceCorner, 3, 2, 0.0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteContRes"
PetscErrorCode PVOutWriteContRes(JacRes *jr, OutBuf *outbuf)
{
	ACCESS_FUNCTION_HEADER

	cf  = scal->strain_rate;

	ierr = JacResCopyContinuityRes(jr, jr->gres); CHKERRQ(ierr);

	GLOBAL_TO_LOCAL(outbuf->fs->DA_CEN, jr->gc, outbuf->lbcen)

	INTERPOLATE_ACCESS(jr->gc, InterpCenterCorner, 1, 0, 0.0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWritEnergRes"
PetscErrorCode PVOutWritEnergRes(JacRes *jr, OutBuf *outbuf)
{
	ACCESS_FUNCTION_HEADER

	cf = scal->dissipation_rate;

	// scatter to local vector
	GLOBAL_TO_LOCAL(outbuf->fs->DA_CEN, jr->ge, outbuf->lbcen)

	INTERPOLATE_ACCESS(jr->ge, InterpCenterCorner, 1, 0, 0.0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteDII_CEN"
PetscErrorCode PVOutWriteDII_CEN(JacRes *jr, OutBuf *outbuf)
{
	COPY_FUNCTION_HEADER

	// macros to copy effective strain rate invariant to buffer
	#define GET_DII_CENTER buff[k][j][i] = jr->svCell[iter++].svDev.DII;

	cf  = scal->strain_rate;

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_DII_CENTER,  1, 0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteDII_XY"
PetscErrorCode PVOutWriteDII_XY(JacRes *jr, OutBuf *outbuf)
{
	COPY_FUNCTION_HEADER

	// macros to copy effective strain rate invariant to buffer
	#define GET_DII_XY_EDGE buff[k][j][i] = jr->svXYEdge[iter++].svDev.DII;

	cf = scal->strain_rate;

	INTERPOLATE_COPY(fs->DA_XY, outbuf->lbxy, InterpXYEdgeCorner, GET_DII_XY_EDGE, 1, 0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteDII_XZ"
PetscErrorCode PVOutWriteDII_XZ(JacRes *jr, OutBuf *outbuf)
{
	COPY_FUNCTION_HEADER

	// macros to copy effective strain rate invariant to buffer
	#define GET_DII_XZ_EDGE buff[k][j][i] = jr->svXZEdge[iter++].svDev.DII;

	cf = scal->strain_rate;

	INTERPOLATE_COPY(fs->DA_XZ, outbuf->lbxz, InterpXZEdgeCorner, GET_DII_XZ_EDGE, 1, 0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteDII_YZ"
PetscErrorCode PVOutWriteDII_YZ(JacRes *jr, OutBuf *outbuf)
{
	COPY_FUNCTION_HEADER

	// macros to copy effective strain rate invariant to buffer
	#define GET_DII_YZ_EDGE buff[k][j][i] = jr->svYZEdge[iter++].svDev.DII;

	cf = scal->strain_rate;

	INTERPOLATE_COPY(fs->DA_YZ, outbuf->lbyz, InterpYZEdgeCorner, GET_DII_YZ_EDGE, 1, 0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

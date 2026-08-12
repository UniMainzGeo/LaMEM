/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : LaMEM
 **   License      : MIT, see LICENSE file for details
 **   Contributors : Anton Popov, Boris Kaus, see AUTHORS file for complete list
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : kaus@uni-mainz.de, popov@uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
//---------------------------------------------------------------------------
//....................   FDSTAG VECTOR OUTPUT ROUTINES   ....................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "outFunct.h"
#include "tssolve.h"
#include "scaling.h"
#include "fdstag.h"
#include "phase.h"
#include "JacRes.h"
#include "interpolate.h"
#include "paraViewOutBin.h"

//---------------------------------------------------------------------------
// WARNING!
//
// ParaView symmetric tensor components ordering is: xx, yy, zz, xy, yz, xz
//
// This is diagonal (rather than row-wise) storage format !!!
//
// ParaView TensorGlyph-Plugin requires a complete 9 component tensor, ordering is row-wise:
// xx, xy, xz, yx, yy, yz, zx, zy, zz
//
// As usual, this isn't documented anywhere !!! Take care of this in future versions.
//---------------------------------------------------------------------------
// interpolation function header
#define COPY_FUNCTION_HEADER \
    JacRes      *jr; \
    FDSTAG      *fs; \
    OutBuf      *outbuf; \
    Scaling     *scal; \
    PetscScalar ***buff, cf; \
    Vec         lbcen, lbcor; \
    PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter; \
    InterpFlags iflag; \
    PetscFunctionBeginUser; \
    outbuf = outvec->outbuf; \
    jr     = outvec->jr; \
    fs     = jr->fs; \
    scal   = jr->scal; \
    iflag.update    = 0; \
    iflag.use_bound = 0; \
    PetscCall(DMGetLocalVectorClean(fs->DA_CEN, &lbcen)); \
    PetscCall(DMGetLocalVectorClean(fs->DA_COR, &lbcor));
//---------------------------------------------------------------------------
#define COPY_FUNCTION_FOOTER \
    PetscCall(DMRestoreLocalVector(fs->DA_CEN, &lbcen)); \
    PetscCall(DMRestoreLocalVector(fs->DA_COR, &lbcor)); \
    PetscFunctionReturn(0);
//---------------------------------------------------------------------------
// access function header
#define ACCESS_FUNCTION_HEADER \
    JacRes      *jr; \
    FDSTAG      *fs; \
    OutBuf      *outbuf; \
    Scaling     *scal; \
    PetscScalar  cf; \
    InterpFlags  iflag; \
    Vec          lbcor; \
    PetscFunctionBeginUser; \
    outbuf = outvec->outbuf; \
    jr     = outvec->jr; \
    fs     = jr->fs; \
    scal   = jr->scal; \
    iflag.update    = 0; \
    iflag.use_bound = 0; \
    PetscCall(DMGetLocalVectorClean(fs->DA_COR, &lbcor));
//---------------------------------------------------------------------------
#define ACCESS_FUNCTION_FOOTER \
    PetscCall(DMRestoreLocalVector(fs->DA_COR, &lbcor)); \
    PetscFunctionReturn(0);
//---------------------------------------------------------------------------
#define COPY_TO_LOCAL_BUFFER(da, vec, FIELD) \
    PetscCall(DMDAGetCorners (da, &sx, &sy, &sz, &nx, &ny, &nz)); \
    PetscCall(DMDAVecGetArray(da, vec, &buff)); \
    iter = 0; \
    START_STD_LOOP \
        FIELD \
    END_STD_LOOP \
    PetscCall(DMDAVecRestoreArray(da, vec, &buff)); \
    LOCAL_TO_LOCAL(da, vec)
//---------------------------------------------------------------------------
#define INTERPOLATE_COPY(da, vec, IFUNCT, FIELD, ncomp, dir) \
    COPY_TO_LOCAL_BUFFER(da, vec, FIELD) \
    PetscCall(IFUNCT(fs, vec, lbcor, iflag)); \
    if(!iflag.update) { PetscCall(OutBufPut3DVecComp(outbuf, lbcor, ncomp, dir, cf, 0.0)); }
//---------------------------------------------------------------------------
#define INTERPOLATE_ACCESS(vec, IFUNCT, ncomp, dir, shift) \
    PetscCall(IFUNCT(outbuf->fs, vec, lbcor, iflag)); \
    PetscCall(OutBufPut3DVecComp(outbuf, lbcor, ncomp, dir, cf, shift));
//---------------------------------------------------------------------------
//...........  Multi-component output vector data structure .................
//---------------------------------------------------------------------------
void OutVecCreate(
    OutVec         *outvec,
    JacRes         *jr,
    OutBuf         *outbuf,
    const char     *name,
    const char     *label,
    PetscErrorCode (*OutVecWrite)(OutVec*),
    PetscInt        num,
    PetscInt       *phase_ID)
{
	PetscInt i;

	// context
	outvec->jr     = jr;
	outvec->outbuf = outbuf;

	// store name
	sprintf(outvec->name, "%s %s", name, label);

	// phase mask for phase aggregate
	if(phase_ID)
	{
		// set number of components
		outvec->ncomp = 1;

		// setup phase mask
		for(i = 0; i < num; i++)
		{
			outvec->phase_mask[phase_ID[i]] = 1;
		}
	}
	else
	{
		// set number of components
		outvec->ncomp = num;
	}

	// output function pointer
	outvec->OutVecWrite = OutVecWrite;
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWritePhase(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	Material_t  *phases;
	PetscScalar *phRat, mID;
	PetscInt     jj, numPhases;

	// macro to copy phase parameter to buffer
#define GET_PHASE \
     phRat = jr->svCell[iter++].phRat; \
     mID = 0.0; \
     for(jj = 0; jj < numPhases; jj++) \
         mID += phRat[jj]*(PetscScalar)phases[jj].visID; \
     buff[k][j][i] = mID;

	cf = scal->unit;

	// access material parameters
	phases    = jr->dbm->phases;
	numPhases = jr->dbm->numPhases;

	INTERPOLATE_COPY(fs->DA_CEN, lbcen, InterpCenterCorner, GET_PHASE, 1, 0)

	COPY_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWritePhaseAgg(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	PetscScalar *phRat, agg;
	PetscInt     jj, numPhases, *phase_mask;

	// macro to copy aggregated phase ratio to buffer
#define GET_PHASE_AGG \
     phRat = jr->svCell[iter++].phRat; \
     agg   = 0.0; \
     for(jj = 0; jj < numPhases; jj++) \
         if(phase_mask[jj]) agg += phRat[jj]; \
     buff[k][j][i] = agg;

	cf = scal->unit;

	// access material parameters
	numPhases  = jr->dbm->numPhases;
	phase_mask = outvec->phase_mask;

	INTERPOLATE_COPY(fs->DA_CEN, lbcen, InterpCenterCorner, GET_PHASE_AGG, 1, 0)

	COPY_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteDensity(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	// macro to copy density to buffer
#define GET_DENSITY buff[k][j][i] = jr->svCell[iter++].svBulk.rho;

	cf = scal->density;

	INTERPOLATE_COPY(fs->DA_CEN, lbcen, InterpCenterCorner, GET_DENSITY, 1, 0)

	COPY_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteViscTotal(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	// macro to copy viscosity to buffer
#define GET_VISC_TOTAL buff[k][j][i] = jr->svCell[iter++].svDev.eta;

	// output viscosity logarithm in GEO-mode
	// (negative scaling requests logarithmic output)
	if(scal->utype == _GEO_) cf = -scal->viscosity;
	else                     cf =  scal->viscosity;

	INTERPOLATE_COPY(fs->DA_CEN, lbcen, InterpCenterCorner, GET_VISC_TOTAL, 1, 0)

	COPY_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteViscCreep(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	// macro to copy viscosity to buffer
#define GET_VISC_CREEP buff[k][j][i] = jr->svCell[iter++].eta_cr;

	// output viscosity logarithm in GEO-mode
	// (negative scaling requests logarithmic output)
	if(scal->utype == _GEO_) cf = -scal->viscosity;
	else                     cf =  scal->viscosity;

	INTERPOLATE_COPY(fs->DA_CEN, lbcen, InterpCenterCorner, GET_VISC_CREEP, 1, 0)

	COPY_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteVelocity(OutVec* outvec)
{
	ACCESS_FUNCTION_HEADER

	Vec lvx, lvy, lvz;

	cf              = scal->velocity;
	iflag.use_bound = 1;

	// get velocity vectors
	PetscCall(JacResGetSolution(jr, jr->gsol, &lvx, &lvy, &lvz, NULL, NULL, _interp_));

	INTERPOLATE_ACCESS(lvx, InterpXFaceCorner, 3, 0, 0.0)
	INTERPOLATE_ACCESS(lvy, InterpYFaceCorner, 3, 1, 0.0)
	INTERPOLATE_ACCESS(lvz, InterpZFaceCorner, 3, 2, 0.0)

	// restore velocity vectors
	PetscCall(JacResRestoreSolution(jr, &lvx, &lvy, &lvz, NULL, NULL));

	ACCESS_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWritePressure(OutVec* outvec)
{
	PetscFunctionBeginUser;

	if(outvec->jr->ctrl.gwType != _GW_NONE_)
	{
		PetscCall(PVOutWriteTotalPress(outvec));
	}
	else
	{
		PetscCall(PVOutWriteEffPress(outvec));
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteStAngle(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	// macro to copy PSD angle to buffer
#define GET_STANGLE buff[k][j][i] = jr->svCell[iter++].svBulk.phi;

	cf = scal->unit;

	INTERPOLATE_COPY(fs->DA_CEN, lbcen, InterpCenterCorner, GET_STANGLE, 1, 0)

	COPY_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteTotalPress(OutVec* outvec)
{
	ACCESS_FUNCTION_HEADER

	Vec         lp, lpt;
	PetscScalar pShift, biot;

	biot   =  jr->ctrl.biot;
	cf     =  scal->stress;
	pShift = -cf*jr->ctrl.pShift;

	PetscCall(DMGetLocalVectorClean(fs->DA_CEN, &lpt));

	PetscCall(JacResGetSolution(jr, jr->gsol, NULL, NULL, NULL, &lp, NULL, _no_interp_));

	// compute total pressure (add pore fluid pressure)
	PetscCall(VecWAXPY(lpt, biot, jr->lp_pore, lp));

	// set corners and edges for interpolation
	PetscCall(FDSTAGSetEdgeCornerCenter(jr->fs, lpt));

	INTERPOLATE_ACCESS(lpt, InterpCenterCorner, 1, 0, pShift)

	PetscCall(JacResRestoreSolution(jr, NULL, NULL, NULL, &lp, NULL));

	PetscCall(DMRestoreLocalVector(fs->DA_CEN, &lpt));

	ACCESS_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteEffPress(OutVec* outvec)
{
	ACCESS_FUNCTION_HEADER

	Vec         lp;
	PetscScalar pShift;

	cf              = scal->stress;
	iflag.use_bound = 1;
	pShift          = -cf*jr->ctrl.pShift;

	PetscCall(JacResGetSolution(jr, jr->gsol, NULL, NULL, NULL, &lp, NULL, _interp_));

	INTERPOLATE_ACCESS(lp, InterpCenterCorner, 1, 0, pShift)

	PetscCall(JacResRestoreSolution(jr, NULL, NULL, NULL, &lp, NULL));

	ACCESS_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteOverPress(OutVec* outvec)
{
	ACCESS_FUNCTION_HEADER

	Vec lop;

	PetscScalar pShift;

	cf     =  scal->stress;
	pShift = -cf*jr->ctrl.pShift;

	PetscCall(DMGetLocalVectorClean(fs->DA_CEN, &lop));

	PetscCall(JacResGetOverPressure(jr, lop));

	INTERPOLATE_ACCESS(lop, InterpCenterCorner, 1, 0, pShift)

	PetscCall(DMRestoreLocalVector(fs->DA_CEN, &lop));

	ACCESS_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteLithoPress(OutVec* outvec)
{
	ACCESS_FUNCTION_HEADER

	cf = scal->stress;

	INTERPOLATE_ACCESS(jr->lp_lith, InterpCenterCorner, 1, 0, 0.0)

	ACCESS_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWritePorePress(OutVec* outvec)
{
	ACCESS_FUNCTION_HEADER

	cf = scal->stress;

	INTERPOLATE_ACCESS(jr->lp_pore, InterpCenterCorner, 1, 0, 0.0)

	ACCESS_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteTemperature(OutVec* outvec)
{
	ACCESS_FUNCTION_HEADER

	Vec lT;

	cf              = scal->temperature;
	iflag.use_bound = 1;

	PetscCall(JacResGetSolution(jr, jr->gsol, NULL, NULL, NULL, NULL, &lT, _interp_));

	INTERPOLATE_ACCESS(lT, InterpCenterCorner, 1, 0, scal->Tshift)

	PetscCall(JacResRestoreSolution(jr, NULL, NULL, NULL, NULL, &lT));

	ACCESS_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteConductivity(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	// macros to copy conductivity to buffer
#define GET_COND_CENTER buff[k][j][i] = jr->svCell[iter++].svBulk.cond;

	cf = scal->conductivity;

	INTERPOLATE_COPY(fs->DA_CEN, lbcen, InterpCenterCorner, GET_COND_CENTER, 1, 0)

	COPY_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteDevStress(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	// NOTE! See warning about component ordering scheme above

	SolVarEdge  *svEdge;
	SolVarCell  *svCell;
	PetscScalar  pf;
	Vec          lbxy, lbxz, lbyz;

	// get pre-factor
	if(jr->ctrl.initGuess) pf = 0.0;
	else                   pf = 2.0;

	// macro to copy deviatoric stress components to buffer
#define GET_SXX { svCell = &jr->svCell  [iter++]; buff[k][j][i] = svCell->sxx + pf*svCell->svDev.eta_st*svCell->dxx; }
#define GET_SYY { svCell = &jr->svCell  [iter++]; buff[k][j][i] = svCell->syy + pf*svCell->svDev.eta_st*svCell->dyy; }
#define GET_SZZ { svCell = &jr->svCell  [iter++]; buff[k][j][i] = svCell->szz + pf*svCell->svDev.eta_st*svCell->dzz; }
#define GET_SXY { svEdge = &jr->svXYEdge[iter++]; buff[k][j][i] = svEdge->s   + pf*svEdge->svDev.eta_st*svEdge->d;   }
#define GET_SYZ { svEdge = &jr->svYZEdge[iter++]; buff[k][j][i] = svEdge->s   + pf*svEdge->svDev.eta_st*svEdge->d;   }
#define GET_SXZ { svEdge = &jr->svXZEdge[iter++]; buff[k][j][i] = svEdge->s   + pf*svEdge->svDev.eta_st*svEdge->d;   }

	cf = scal->stress;

	PetscCall(FDSTAGGetLocalVectorEdge(fs, &lbxy, &lbxz, &lbyz));

	INTERPOLATE_COPY(fs->DA_CEN, lbcen, InterpCenterCorner, GET_SXX, 9, 0)
	INTERPOLATE_COPY(fs->DA_XY,  lbxy,  InterpXYEdgeCorner, GET_SXY, 9, 1)
	INTERPOLATE_COPY(fs->DA_XZ,  lbxz,  InterpXZEdgeCorner, GET_SXZ, 9, 2)
	INTERPOLATE_COPY(fs->DA_XY,  lbxy,  InterpXYEdgeCorner, GET_SXY, 9, 3)
	INTERPOLATE_COPY(fs->DA_CEN, lbcen, InterpCenterCorner, GET_SYY, 9, 4)
	INTERPOLATE_COPY(fs->DA_YZ,  lbyz,  InterpYZEdgeCorner, GET_SYZ, 9, 5)
	INTERPOLATE_COPY(fs->DA_XZ,  lbxz,  InterpXZEdgeCorner, GET_SXZ, 9, 6)
	INTERPOLATE_COPY(fs->DA_YZ,  lbyz,  InterpYZEdgeCorner, GET_SYZ, 9, 7)
	INTERPOLATE_COPY(fs->DA_CEN, lbcen, InterpCenterCorner, GET_SZZ, 9, 8)

	PetscCall(FDSTAGRestoreLocalVectorEdge(fs, &lbxy, &lbxz, &lbyz));

	COPY_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteJ2DevStress(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	SolVarCell  *svCell;
	SolVarEdge  *svEdge;
	PetscScalar s, J2, pf;
	Vec         lbxy, lbxz, lbyz;

	// get pre-factor
	if(jr->ctrl.initGuess) pf = 0.0;
	else                   pf = 2.0;

	// macros to copy deviatoric strain rate invariant to buffer
#define GET_J2_STRESS_CENTER \
     svCell = &jr->svCell[iter++]; \
     s = svCell->sxx + pf*svCell->svDev.eta_st*svCell->dxx; J2  = s*s; \
     s = svCell->syy + pf*svCell->svDev.eta_st*svCell->dyy; J2 += s*s; \
     s = svCell->szz + pf*svCell->svDev.eta_st*svCell->dzz; J2 += s*s; \
     buff[k][j][i] = 0.5*J2;

#define GET_J2_STRESS_XY_EDGE { svEdge = &jr->svXYEdge[iter++]; s = svEdge->s + pf*svEdge->svDev.eta_st*svEdge->d; buff[k][j][i] = s*s;}
#define GET_J2_STRESS_YZ_EDGE { svEdge = &jr->svYZEdge[iter++]; s = svEdge->s + pf*svEdge->svDev.eta_st*svEdge->d; buff[k][j][i] = s*s;}
#define GET_J2_STRESS_XZ_EDGE { svEdge = &jr->svXZEdge[iter++]; s = svEdge->s + pf*svEdge->svDev.eta_st*svEdge->d; buff[k][j][i] = s*s;}

	cf           = scal->stress;
	iflag.update = 1;

	PetscCall(FDSTAGGetLocalVectorEdge(fs, &lbxy, &lbxz, &lbyz));

	PetscCall(VecSet(lbcor, 0.0));

	INTERPOLATE_COPY(fs->DA_CEN, lbcen, InterpCenterCorner, GET_J2_STRESS_CENTER,  1, 0)
	INTERPOLATE_COPY(fs->DA_XY,  lbxy,  InterpXYEdgeCorner, GET_J2_STRESS_XY_EDGE, 1, 0)
	INTERPOLATE_COPY(fs->DA_YZ,  lbyz,  InterpYZEdgeCorner, GET_J2_STRESS_YZ_EDGE, 1, 0)
	INTERPOLATE_COPY(fs->DA_XZ,  lbxz,  InterpXZEdgeCorner, GET_J2_STRESS_XZ_EDGE, 1, 0)

	// compute & store second invariant
	PetscCall(VecSqrtAbs(lbcor));

	PetscCall(OutBufPut3DVecComp(outbuf, lbcor, 1, 0, cf, 0.0));

	PetscCall(FDSTAGRestoreLocalVectorEdge(fs, &lbxy, &lbxz, &lbyz));

	COPY_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteStrainRate(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	Vec lbxy, lbxz, lbyz;

	// NOTE! See warning about component ordering scheme above

	// macro to copy deviatoric strain rate components to buffer
#define GET_DXX buff[k][j][i] = jr->svCell[iter++].dxx;
#define GET_DYY buff[k][j][i] = jr->svCell[iter++].dyy;
#define GET_DZZ buff[k][j][i] = jr->svCell[iter++].dzz;
#define GET_DXY buff[k][j][i] = jr->svXYEdge[iter++].d;
#define GET_DYZ buff[k][j][i] = jr->svYZEdge[iter++].d;
#define GET_DXZ buff[k][j][i] = jr->svXZEdge[iter++].d;

	cf = scal->strain_rate;

	PetscCall(FDSTAGGetLocalVectorEdge(fs, &lbxy, &lbxz, &lbyz));

	INTERPOLATE_COPY(fs->DA_CEN, lbcen, InterpCenterCorner, GET_DXX, 9, 0)
	INTERPOLATE_COPY(fs->DA_XY,  lbxy,  InterpXYEdgeCorner, GET_DXY, 9, 1)
	INTERPOLATE_COPY(fs->DA_XZ,  lbxz,  InterpXZEdgeCorner, GET_DXZ, 9, 2)
	INTERPOLATE_COPY(fs->DA_XY,  lbxy,  InterpXYEdgeCorner, GET_DXY, 9, 3)
	INTERPOLATE_COPY(fs->DA_CEN, lbcen, InterpCenterCorner, GET_DYY, 9, 4)
	INTERPOLATE_COPY(fs->DA_YZ,  lbyz,  InterpYZEdgeCorner, GET_DYZ, 9, 5)
	INTERPOLATE_COPY(fs->DA_XZ,  lbxz,  InterpXZEdgeCorner, GET_DXZ, 9, 6)
	INTERPOLATE_COPY(fs->DA_YZ,  lbyz,  InterpYZEdgeCorner, GET_DYZ, 9, 7)
	INTERPOLATE_COPY(fs->DA_CEN, lbcen, InterpCenterCorner, GET_DZZ, 9, 8)

	PetscCall(FDSTAGRestoreLocalVectorEdge(fs, &lbxy, &lbxz, &lbyz));

	COPY_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteJ2StrainRate(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	Vec        lbxy, lbxz, lbyz;
	SolVarCell *svCell;
	PetscScalar d, J2;

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

	cf           = scal->strain_rate;
	iflag.update = 1;

	PetscCall(FDSTAGGetLocalVectorEdge(fs, &lbxy, &lbxz, &lbyz));

	PetscCall(VecSet(lbcor, 0.0));

	INTERPOLATE_COPY(fs->DA_CEN, lbcen, InterpCenterCorner, GET_J2_STRAIN_RATE_CENTER,  1, 0)
	INTERPOLATE_COPY(fs->DA_XY,  lbxy,  InterpXYEdgeCorner, GET_J2_STRAIN_RATE_XY_EDGE, 1, 0)
	INTERPOLATE_COPY(fs->DA_YZ,  lbyz,  InterpYZEdgeCorner, GET_J2_STRAIN_RATE_YZ_EDGE, 1, 0)
	INTERPOLATE_COPY(fs->DA_XZ,  lbxz,  InterpXZEdgeCorner, GET_J2_STRAIN_RATE_XZ_EDGE, 1, 0)

	// compute & store second invariant
	PetscCall(VecSqrtAbs(lbcor));

	PetscCall(OutBufPut3DVecComp(outbuf, lbcor, 1, 0, cf, 0.0));

	PetscCall(FDSTAGRestoreLocalVectorEdge(fs, &lbxy, &lbxz, &lbyz));

	COPY_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteFluidDensity(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	// macros to copy fluid density to buffer
#define GET_RHOPF_CENTER  buff[k][j][i] = jr->svCell[iter++].svBulk.rho_pf;

	cf = scal->density;

	INTERPOLATE_COPY(fs->DA_CEN, lbcen, InterpCenterCorner, GET_RHOPF_CENTER,  1, 0)

	COPY_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteMeltFraction(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	// macros to copy melt fraction to buffer
#define GET_MF_CENTER  buff[k][j][i] = jr->svCell[iter++].svBulk.mf;

	cf = scal->unit;

	INTERPOLATE_COPY(fs->DA_CEN, lbcen, InterpCenterCorner, GET_MF_CENTER,  1, 0)

	COPY_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteVolRate(OutVec* outvec)
{
	PetscFunctionBeginUser;

	UNUSED(outvec);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteVorticity(OutVec* outvec)
{
	PetscFunctionBeginUser;

	UNUSED(outvec);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteAngVelMag(OutVec* outvec)
{
	PetscFunctionBeginUser;

	UNUSED(outvec);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteTotStrain(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	// macro to copy accumulated total strain (ATS) to buffer
#define GET_ATS buff[k][j][i] = jr->svCell[iter++].ATS;

	cf = scal->unit;

	INTERPOLATE_COPY(fs->DA_CEN, lbcen, InterpCenterCorner, GET_ATS, 1, 0)

	COPY_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWritePlastStrain(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	// macro to copy accumulated plastic strain (APS) to buffer
#define GET_APS buff[k][j][i] = jr->svCell[iter++].svDev.APS;

	cf = scal->unit;

	INTERPOLATE_COPY(fs->DA_CEN, lbcen, InterpCenterCorner, GET_APS, 1, 0)

	COPY_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWritePlastDissip(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	Vec        lbxy, lbxz, lbyz;
	SolVarCell *svCell;
	SolVarEdge *svEdge;
	PetscScalar Hr;

	// macros to copy shear heating  to buffer
#define GET_SHEAR_HEATING_CENTER \
     svCell = &jr->svCell[iter++];  \
     Hr = svCell->svDev.Hr; \
     buff[k][j][i] = Hr;

#define GET_SHEAR_HEATING_XY_EDGE svEdge = &jr->svXYEdge[iter++]; Hr = svEdge->svDev.Hr; buff[k][j][i] = Hr;
#define GET_SHEAR_HEATING_YZ_EDGE svEdge = &jr->svYZEdge[iter++]; Hr = svEdge->svDev.Hr; buff[k][j][i] = Hr;
#define GET_SHEAR_HEATING_XZ_EDGE svEdge = &jr->svXZEdge[iter++]; Hr = svEdge->svDev.Hr; buff[k][j][i] = Hr;

	cf           = scal->dissipation_rate;
	iflag.update = 1;

	PetscCall(FDSTAGGetLocalVectorEdge(fs, &lbxy, &lbxz, &lbyz));

	PetscCall(VecSet(lbcor, 0.0));

	INTERPOLATE_COPY(fs->DA_CEN, lbcen, InterpCenterCorner, GET_SHEAR_HEATING_CENTER,  1, 0)
	INTERPOLATE_COPY(fs->DA_XY,  lbxy,  InterpXYEdgeCorner, GET_SHEAR_HEATING_XY_EDGE, 1, 0)
	INTERPOLATE_COPY(fs->DA_YZ,  lbyz,  InterpYZEdgeCorner, GET_SHEAR_HEATING_YZ_EDGE, 1, 0)
	INTERPOLATE_COPY(fs->DA_XZ,  lbxz,  InterpXZEdgeCorner, GET_SHEAR_HEATING_XZ_EDGE, 1, 0)

	PetscCall(OutBufPut3DVecComp(outbuf, lbcor, 1, 0, cf, 0.0));

	PetscCall(FDSTAGRestoreLocalVectorEdge(fs, &lbxy, &lbxz, &lbyz));

	COPY_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteTotDispl(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	cf = scal->length;

	// macros to copy displacement in cell to buffer
#define GET_DISPLX buff[k][j][i] = jr->svCell[iter++].U[0];
#define GET_DISPLY buff[k][j][i] = jr->svCell[iter++].U[1];
#define GET_DISPLZ buff[k][j][i] = jr->svCell[iter++].U[2];

	INTERPOLATE_COPY(fs->DA_CEN, lbcen, InterpCenterCorner, GET_DISPLX, 3, 0);
	INTERPOLATE_COPY(fs->DA_CEN, lbcen, InterpCenterCorner, GET_DISPLY, 3, 1);
	INTERPOLATE_COPY(fs->DA_CEN, lbcen, InterpCenterCorner, GET_DISPLZ, 3, 2);

	COPY_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteSHmax(OutVec* outvec)
{
	ACCESS_FUNCTION_HEADER

	Vec cx, cy;

	cf = scal->unit;

	PetscCall(DMGetLocalVectorClean(fs->DA_CEN, &cx));
	PetscCall(DMGetLocalVectorClean(fs->DA_CEN, &cy));

	// compute maximum horizontal compressive stress (SHmax) orientation
	PetscCall(JacResGetSHmax(jr, cx, cy));

	INTERPOLATE_ACCESS(cx, InterpCenterCorner, 3, 0, 0.0)
	INTERPOLATE_ACCESS(cy, InterpCenterCorner, 3, 1, 0.0)

	PetscCall(OutBufZero3DVecComp(outbuf, 3, 2));

	PetscCall(DMRestoreLocalVector(fs->DA_CEN, &cx));
	PetscCall(DMRestoreLocalVector(fs->DA_CEN, &cy));

	ACCESS_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteEHmax(OutVec* outvec)
{
	ACCESS_FUNCTION_HEADER

	Vec cx, cy;

	cf = scal->unit;

	PetscCall(DMGetLocalVectorClean(fs->DA_CEN, &cx));
	PetscCall(DMGetLocalVectorClean(fs->DA_CEN, &cy));

	// compute maximum horizontal extension rate (EHmax) orientation
	PetscCall(JacResGetEHmax(jr, cx, cy));

	INTERPOLATE_ACCESS(cx, InterpCenterCorner, 3, 0, 0.0)
	INTERPOLATE_ACCESS(cy, InterpCenterCorner, 3, 1, 0.0)

	PetscCall(OutBufZero3DVecComp(outbuf, 3, 2));

	PetscCall(DMRestoreLocalVector(fs->DA_CEN, &cx));
	PetscCall(DMRestoreLocalVector(fs->DA_CEN, &cy));

	ACCESS_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteYield(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	// macro to copy yield stress to buffer
#define GET_YIELD buff[k][j][i] = jr->svCell[iter++].yield;

	cf = scal->stress;

	INTERPOLATE_COPY(fs->DA_CEN, lbcen, InterpCenterCorner, GET_YIELD, 1, 0)

	COPY_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteRelDIIdif(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	// macro to copy diffusion creep relative strain rate to buffer
#define GET_DIIdif buff[k][j][i] = jr->svCell[iter++].DIIdif;

	cf = scal->unit;

	INTERPOLATE_COPY(fs->DA_CEN, lbcen, InterpCenterCorner, GET_DIIdif, 1, 0)

	COPY_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteRelDIIdis(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	// macro to copy diffusion creep relative strain rate to buffer
#define GET_DIIdis buff[k][j][i] = jr->svCell[iter++].DIIdis;

	cf = scal->unit;

	INTERPOLATE_COPY(fs->DA_CEN, lbcen, InterpCenterCorner, GET_DIIdis, 1, 0)

	COPY_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteRelDIIprl(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	// macro to copy diffusion creep relative strain rate to buffer
#define GET_DIIprl buff[k][j][i] = jr->svCell[iter++].DIIprl;

	cf = scal->unit;

	INTERPOLATE_COPY(fs->DA_CEN, lbcen, InterpCenterCorner, GET_DIIprl, 1, 0)

	COPY_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteRelDIIpl(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	// macro to copy plastic relative strain rate to buffer
#define GET_DIIpl buff[k][j][i] = jr->svCell[iter++].DIIpl;

	cf = scal->unit;

	INTERPOLATE_COPY(fs->DA_CEN, lbcen, InterpCenterCorner, GET_DIIpl, 1, 0)

	COPY_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
// DEBUG VECTORS
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteMomentRes(OutVec* outvec)
{
	ACCESS_FUNCTION_HEADER

	Vec gfx, gfy, gfz;
	Vec lfx, lfy, lfz;

	cf = scal->volumetric_force;

	// make buffer vectors
	PetscCall(FDSTAGGetGlobalVectorFace(fs, &gfx, &gfy, &gfz));
	PetscCall(FDSTAGGetLocalVectorFace (fs, &lfx, &lfy, &lfz));

	// get momentum residuals
	PetscCall(FDSTAGSplitVectors(fs, jr->gres, gfx, gfy, gfz, NULL));

	// exchange ghost values
	GLOBAL_TO_LOCAL(fs->DA_X, gfx, lfx)
	GLOBAL_TO_LOCAL(fs->DA_Y, gfy, lfy)
	GLOBAL_TO_LOCAL(fs->DA_Z, gfz, lfz)

	INTERPOLATE_ACCESS(lfx, InterpXFaceCorner, 3, 0, 0.0)
	INTERPOLATE_ACCESS(lfy, InterpYFaceCorner, 3, 1, 0.0)
	INTERPOLATE_ACCESS(lfz, InterpZFaceCorner, 3, 2, 0.0)

	PetscCall(FDSTAGRestoreGlobalVectorFace(fs, &gfx, &gfy, &gfz));
	PetscCall(FDSTAGRestoreLocalVectorFace (fs, &lfx, &lfy, &lfz));

	ACCESS_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteContRes(OutVec* outvec)
{
	ACCESS_FUNCTION_HEADER

	Vec gc, lc;

	cf = scal->strain_rate;

	PetscCall(DMGetGlobalVector    (fs->DA_CEN, &gc));
	PetscCall(DMGetLocalVectorClean(fs->DA_CEN, &lc));

	// get continuity residual
	PetscCall(FDSTAGSplitVectors(fs, jr->gres, NULL, NULL, NULL, gc));

	GLOBAL_TO_LOCAL(fs->DA_CEN, gc, lc)

	INTERPOLATE_ACCESS(lc, InterpCenterCorner, 1, 0, 0.0)

	PetscCall(DMRestoreGlobalVector(fs->DA_CEN, &gc));
	PetscCall(DMRestoreLocalVector (fs->DA_CEN, &lc));

	ACCESS_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWritEnergRes(OutVec* outvec)
{
	ACCESS_FUNCTION_HEADER

	Vec le;

	cf = scal->dissipation_rate;

	PetscCall(DMGetLocalVectorClean(fs->DA_CEN, &le));

	GLOBAL_TO_LOCAL(fs->DA_CEN, jr->ge, le)

	INTERPOLATE_ACCESS(le, InterpCenterCorner, 1, 0, 0.0)

	PetscCall(DMRestoreLocalVector(fs->DA_CEN, &le));

	ACCESS_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------
PetscErrorCode PVOutWriteVelGrad(OutVec* outvec)
{
	ACCESS_FUNCTION_HEADER

	// NOTE! See warning about component ordering scheme above

	Vec    lvx, lvy, lvz;
	Vec    dvxdx, dvxdy, dvxdz;
	Vec    dvydx, dvydy, dvydz;
	Vec    dvzdx, dvzdy, dvzdz;

	cf = scal->strain_rate;

	// get work vectors
	PetscCall(FDSTAGGetLocalVectorCenter(fs, &dvxdx, &dvydy, &dvzdz));
	PetscCall(FDSTAGGetLocalVectorEdge  (fs, &dvxdy, &dvxdz, &dvydz));
	PetscCall(FDSTAGGetLocalVectorEdge  (fs, &dvydx, &dvzdx, &dvzdy));

	// get velocity vectors
	PetscCall(JacResGetSolution(jr, jr->gsol, &lvx, &lvy, &lvz, NULL, NULL, _no_interp_));

	// compute velocity gradients
	PetscCall(JacResGetVelGrad(jr,
	                           lvx,   lvy,   lvz,
	                           dvxdx, dvxdy, dvxdz,
	                           dvydx, dvydy, dvydz,
	                           dvzdx, dvzdy, dvzdz));

	INTERPOLATE_ACCESS(dvxdx, InterpCenterCorner, 9, 0, 0.0)
	INTERPOLATE_ACCESS(dvxdy, InterpXYEdgeCorner, 9, 1, 0.0)
	INTERPOLATE_ACCESS(dvxdz, InterpXZEdgeCorner, 9, 2, 0.0)
	INTERPOLATE_ACCESS(dvydx, InterpXYEdgeCorner, 9, 3, 0.0)
	INTERPOLATE_ACCESS(dvydy, InterpCenterCorner, 9, 4, 0.0)
	INTERPOLATE_ACCESS(dvydz, InterpYZEdgeCorner, 9, 5, 0.0)
	INTERPOLATE_ACCESS(dvzdx, InterpXZEdgeCorner, 9, 6, 0.0)
	INTERPOLATE_ACCESS(dvzdy, InterpYZEdgeCorner, 9, 7, 0.0)
	INTERPOLATE_ACCESS(dvzdz, InterpCenterCorner, 9, 8, 0.0)

	// get work vectors
	PetscCall(FDSTAGRestoreLocalVectorCenter(fs, &dvxdx, &dvydy, &dvzdz));
	PetscCall(FDSTAGRestoreLocalVectorEdge  (fs, &dvxdy, &dvxdz, &dvydz));
	PetscCall(FDSTAGRestoreLocalVectorEdge  (fs, &dvydx, &dvzdx, &dvzdy));

	// restore velocity vectors
	PetscCall(JacResRestoreSolution(jr, &lvx, &lvy, &lvz, NULL, NULL));

	ACCESS_FUNCTION_FOOTER
}
//---------------------------------------------------------------------------

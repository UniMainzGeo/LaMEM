/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2011-2015, JGU Mainz, Anton Popov, Boris Kaus
 **    All rights reserved.
 **
 **    This software was developed at:
 **
 **         Institute of Geosciences
 **         Johannes-Gutenberg University, Mainz
 **         Johann-Joachim-Becherweg 21
 **         55128 Mainz, Germany
 **
 **    project:    LaMEM
 **    filename:   outFunct.c
 **
 **    LaMEM is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, version 3 of the License.
 **
 **    LaMEM is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with LaMEM. If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    Contact:
 **        Boris Kaus       [kaus@uni-mainz.de]
 **        Anton Popov      [popov@uni-mainz.de]
 **
 **
 **    Main development team:
 **         Anton Popov      [popov@uni-mainz.de]
 **         Boris Kaus       [kaus@uni-mainz.de]
 **         Tobias Baumann
 **         Adina Pusok
 **         Arthur Bauville
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
	OutBuf      *outbuf; \
	FDSTAG      *fs; \
	Scaling     *scal; \
	PetscScalar ***buff, cf; \
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter; \
	InterpFlags iflag; \
	PetscErrorCode ierr; \
	PetscFunctionBegin; \
	jr     = outvec->jr; \
	outbuf = outvec->outbuf; \
	fs     = outbuf->fs; \
	scal   = jr->scal; \
	iflag.update    = 0; \
	iflag.use_bound = 0;
//---------------------------------------------------------------------------
// access function header
#define ACCESS_FUNCTION_HEADER \
	JacRes      *jr; \
	OutBuf      *outbuf; \
	Scaling     *scal; \
	PetscScalar  cf; \
	InterpFlags  iflag; \
	PetscErrorCode ierr; \
	PetscFunctionBegin; \
	jr     = outvec->jr; \
	outbuf = outvec->outbuf; \
	scal   = jr->scal; \
	iflag.update    = 0; \
	iflag.use_bound = 0;
//---------------------------------------------------------------------------
#define COPY_TO_LOCAL_BUFFER(da, vec, FIELD) \
	ierr = DMDAGetCorners (da, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr); \
	ierr = DMDAVecGetArray(da, vec, &buff); CHKERRQ(ierr); \
	iter = 0; \
	START_STD_LOOP \
		FIELD \
	END_STD_LOOP \
	ierr = DMDAVecRestoreArray(da, vec, &buff); CHKERRQ(ierr); \
	LOCAL_TO_LOCAL(da, vec)
//---------------------------------------------------------------------------
#define INTERPOLATE_COPY(da, vec, IFUNCT, FIELD, ncomp, dir) \
	COPY_TO_LOCAL_BUFFER(da, vec, FIELD) \
	ierr = IFUNCT(fs, vec, outbuf->lbcor, iflag); CHKERRQ(ierr); \
	if(!iflag.update) { ierr = OutBufPut3DVecComp(outbuf, ncomp, dir, cf, 0.0); CHKERRQ(ierr); }
//---------------------------------------------------------------------------
#define INTERPOLATE_ACCESS(vec, IFUNCT, ncomp, dir, shift) \
	ierr = IFUNCT(outbuf->fs, vec, outbuf->lbcor, iflag); CHKERRQ(ierr); \
	ierr = OutBufPut3DVecComp(outbuf, ncomp, dir, cf, shift); CHKERRQ(ierr);
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
#undef __FUNCT__
#define __FUNCT__ "PVOutWritePhase"
PetscErrorCode PVOutWritePhase(OutVec* outvec)
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
			mID += phRat[jj]*(PetscScalar)phases[jj].visID; \
		buff[k][j][i] = mID;

	// no scaling is necessary for the phase
	cf = scal->unit;

	// access material parameters
	phases    = jr->dbm->phases;
	numPhases = jr->dbm->numPhases;

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_PHASE, 1, 0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWritePhaseAgg"
PetscErrorCode PVOutWritePhaseAgg(OutVec* outvec)
{
	PetscScalar *phRat, agg;
	PetscInt     jj, numPhases, *phase_mask;

	COPY_FUNCTION_HEADER

	// macro to copy aggregated phase ratio to buffer
	#define GET_PHASE_AGG \
		phRat = jr->svCell[iter++].phRat; \
		agg   = 0.0; \
		for(jj = 0; jj < numPhases; jj++) \
			if(phase_mask[jj]) agg += phRat[jj]; \
		buff[k][j][i] = agg;

	// no scaling is necessary for the phase
	cf = scal->unit;

	// access material parameters
	numPhases  = jr->dbm->numPhases;
	phase_mask = outvec->phase_mask;

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_PHASE_AGG, 1, 0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteDensity"
PetscErrorCode PVOutWriteDensity(OutVec* outvec)
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
#define __FUNCT__ "PVOutWriteViscTotal"
PetscErrorCode PVOutWriteViscTotal(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	// macro to copy viscosity to buffer
	#define GET_VISC_TOTAL buff[k][j][i] = jr->svCell[iter++].svDev.eta;

	// output viscosity logarithm in GEO-mode
	// (negative scaling requests logarithmic output)
	if(scal->utype == _GEO_) cf = -scal->viscosity;
	else                     cf =  scal->viscosity;

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_VISC_TOTAL, 1, 0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteViscCreep"
PetscErrorCode PVOutWriteViscCreep(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	// macro to copy viscosity to buffer
	#define GET_VISC_CREEP buff[k][j][i] = jr->svCell[iter++].eta_cr;

	// output viscosity logarithm in GEO-mode
	// (negative scaling requests logarithmic output)
	if(scal->utype == _GEO_) cf = -scal->viscosity;
	else                     cf =  scal->viscosity;

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_VISC_CREEP, 1, 0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteVelocity"
PetscErrorCode PVOutWriteVelocity(OutVec* outvec)
{
	ACCESS_FUNCTION_HEADER

	cf = scal->velocity;
	iflag.use_bound = 1;

	ierr = JacResCopyVel(jr, jr->gsol); CHKERRQ(ierr);

	INTERPOLATE_ACCESS(jr->lvx, InterpXFaceCorner, 3, 0, 0.0)
	INTERPOLATE_ACCESS(jr->lvy, InterpYFaceCorner, 3, 1, 0.0)
	INTERPOLATE_ACCESS(jr->lvz, InterpZFaceCorner, 3, 2, 0.0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWritePressure"
PetscErrorCode PVOutWritePressure(OutVec* outvec)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	if(outvec->jr->ctrl.gwType != _GW_NONE_)
	{
		ierr = PVOutWriteTotalPress(outvec); CHKERRQ(ierr);
	}
	else
	{
		ierr = PVOutWriteEffPress(outvec); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteGradient"
PetscErrorCode PVOutWriteGradient(OutVec* outvec)
{

	ACCESS_FUNCTION_HEADER

	cf = scal->unit;

	INTERPOLATE_ACCESS(jr->lgradfield, InterpCenterCorner, 1, 0, 0.0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteStAngle"
PetscErrorCode PVOutWriteStAngle(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	// macro to copy PSD angle to buffer
	#define GET_STANGLE buff[k][j][i] = jr->svCell[iter++].svBulk.phi;

	cf = scal->unit;

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_STANGLE, 1, 0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteTotalPress"
PetscErrorCode PVOutWriteTotalPress(OutVec* outvec)
{
	PetscScalar pShift, biot;

	ACCESS_FUNCTION_HEADER

	biot 	= jr->ctrl.biot;
	
	cf  	= scal->stress;

	// scale pressure shift
	pShift 	= -cf*jr->ctrl.pShift;		// minus to be consistent with output routine
	
	ierr = JacResCopyPres(jr, jr->gsol); CHKERRQ(ierr);

	// compute total pressure [add pore fluid P]
	ierr = VecWAXPY(outbuf->lbcen, biot, jr->lp_pore, jr->lp); CHKERRQ(ierr);
	
	INTERPOLATE_ACCESS(outbuf->lbcen, InterpCenterCorner, 1, 0, pShift)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteEffPress"
PetscErrorCode PVOutWriteEffPress(OutVec* outvec)
{
	PetscScalar pShift;

	ACCESS_FUNCTION_HEADER

	cf = scal->stress;
	iflag.use_bound = 1;

	// scale pressure shift
	pShift = -cf*jr->ctrl.pShift;

	INTERPOLATE_ACCESS(jr->lp, InterpCenterCorner, 1, 0, pShift)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteOverPress"
PetscErrorCode PVOutWriteOverPress(OutVec* outvec)
{
	PetscScalar pShift;
	ACCESS_FUNCTION_HEADER

	cf = scal->stress;
	
	// scale pressure shift
	pShift 	= -cf*jr->ctrl.pShift;
	ierr 	= JacResGetOverPressure(jr, outbuf->lbcen); CHKERRQ(ierr);

	INTERPOLATE_ACCESS(outbuf->lbcen, InterpCenterCorner, 1, 0, pShift)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteLithoPress"
PetscErrorCode PVOutWriteLithoPress(OutVec* outvec)
{
	ACCESS_FUNCTION_HEADER

	cf = scal->stress;

	INTERPOLATE_ACCESS(jr->lp_lith, InterpCenterCorner, 1, 0, 0.0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWritePorePress"
PetscErrorCode PVOutWritePorePress(OutVec* outvec)
{
	ACCESS_FUNCTION_HEADER

	cf = scal->stress;

	INTERPOLATE_ACCESS(jr->lp_pore, InterpCenterCorner, 1, 0, 0.0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteTemperature"
PetscErrorCode PVOutWriteTemperature(OutVec* outvec)
{
	ACCESS_FUNCTION_HEADER

	cf = scal->temperature;
	iflag.use_bound = 1;

	INTERPOLATE_ACCESS(jr->lT, InterpCenterCorner, 1, 0, scal->Tshift)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteConductivity"    // NEW
PetscErrorCode PVOutWriteConductivity(OutVec* outvec)
{

  COPY_FUNCTION_HEADER

	// macros to copy conductivity to buffer  
        #define GET_COND_CENTER buff[k][j][i] = jr->svCell[iter++].svBulk.cond;

        cf = scal->conductivity;
	
        INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_COND_CENTER, 1, 0)

        PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteDevStress"
PetscErrorCode PVOutWriteDevStress(OutVec* outvec)
{
	// NOTE! See warning about component ordering scheme above

	SolVarEdge  *svEdge;
	SolVarCell  *svCell;
	PetscScalar  pf;

	COPY_FUNCTION_HEADER

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

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_SXX, 9, 0)
	INTERPOLATE_COPY(fs->DA_XY,  outbuf->lbxy,  InterpXYEdgeCorner, GET_SXY, 9, 1)
	INTERPOLATE_COPY(fs->DA_XZ,  outbuf->lbxz,  InterpXZEdgeCorner, GET_SXZ, 9, 2)
	INTERPOLATE_COPY(fs->DA_XY,  outbuf->lbxy,  InterpXYEdgeCorner, GET_SXY, 9, 3)
	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_SYY, 9, 4)
	INTERPOLATE_COPY(fs->DA_YZ,  outbuf->lbyz,  InterpYZEdgeCorner, GET_SYZ, 9, 5)
	INTERPOLATE_COPY(fs->DA_XZ,  outbuf->lbxz,  InterpXZEdgeCorner, GET_SXZ, 9, 6)
	INTERPOLATE_COPY(fs->DA_YZ,  outbuf->lbyz,  InterpYZEdgeCorner, GET_SYZ, 9, 7)
	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_SZZ, 9, 8)
	

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteJ2DevStress"
PetscErrorCode PVOutWriteJ2DevStress(OutVec* outvec)
{
	SolVarCell  *svCell;
	SolVarEdge  *svEdge;
	PetscScalar s, J2, pf;

	COPY_FUNCTION_HEADER

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

	cf = scal->stress;

	iflag.update = 1;

	ierr = VecSet(outbuf->lbcor, 0.0); CHKERRQ(ierr);

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
PetscErrorCode PVOutWriteStrainRate(OutVec* outvec)
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

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_DXX, 9, 0)
	INTERPOLATE_COPY(fs->DA_XY,  outbuf->lbxy,  InterpXYEdgeCorner, GET_DXY, 9, 1)
	INTERPOLATE_COPY(fs->DA_XZ,  outbuf->lbxz,  InterpXZEdgeCorner, GET_DXZ, 9, 2)
	INTERPOLATE_COPY(fs->DA_XY,  outbuf->lbxy,  InterpXYEdgeCorner, GET_DXY, 9, 3)
	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_DYY, 9, 4)
	INTERPOLATE_COPY(fs->DA_YZ,  outbuf->lbyz,  InterpYZEdgeCorner, GET_DYZ, 9, 5)
	INTERPOLATE_COPY(fs->DA_XZ,  outbuf->lbxz,  InterpXZEdgeCorner, GET_DXZ, 9, 6)
	INTERPOLATE_COPY(fs->DA_YZ,  outbuf->lbyz,  InterpYZEdgeCorner, GET_DYZ, 9, 7)
	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_DZZ, 9, 8)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteJ2StrainRate"
PetscErrorCode PVOutWriteJ2StrainRate(OutVec* outvec)
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

	iflag.update = 1;

	ierr = VecSet(outbuf->lbcor, 0.0); CHKERRQ(ierr);

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
#define __FUNCT__ "PVOutWriteFluidDensity"
PetscErrorCode PVOutWriteFluidDensity(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	// macros to copy fluid density to buffer
	#define GET_RHOPF_CENTER  buff[k][j][i] = jr->svCell[iter++].svBulk.rho_pf;

	cf = scal->density;

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_RHOPF_CENTER,  1, 0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteMeltFraction"
PetscErrorCode PVOutWriteMeltFraction(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	// macros to copy melt fraction to buffer
	#define GET_MF_CENTER  buff[k][j][i] = jr->svCell[iter++].svBulk.mf;

	cf = scal->unit;

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_MF_CENTER,  1, 0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteVolRate"
PetscErrorCode PVOutWriteVolRate(OutVec* outvec)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = 0;  CHKERRQ(ierr);
	if(outvec) outvec = NULL;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteVorticity"
PetscErrorCode PVOutWriteVorticity(OutVec* outvec)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = 0;  CHKERRQ(ierr);
	if(outvec) outvec = NULL;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteAngVelMag"
PetscErrorCode PVOutWriteAngVelMag(OutVec* outvec)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = 0;  CHKERRQ(ierr);
	if(outvec) outvec = NULL;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteTotStrain"
PetscErrorCode PVOutWriteTotStrain(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	// macro to copy accumulated total strain (ATS) to buffer
	#define GET_ATS buff[k][j][i] = jr->svCell[iter++].ATS;

	cf = scal->unit;

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_ATS, 1, 0)

	PetscFunctionReturn(0);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWritePlastStrain"
PetscErrorCode PVOutWritePlastStrain(OutVec* outvec)
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
PetscErrorCode PVOutWritePlastDissip(OutVec* outvec)
{
	SolVarCell *svCell;
	SolVarEdge *svEdge;
	PetscScalar Hr;

	COPY_FUNCTION_HEADER

	// macros to copy shear heating  to buffer
	#define GET_SHEAR_HEATING_CENTER \
		svCell = &jr->svCell[iter++];  \
		Hr = svCell->svDev.Hr; \
		buff[k][j][i] = Hr;

	#define GET_SHEAR_HEATING_XY_EDGE svEdge = &jr->svXYEdge[iter++]; Hr = svEdge->svDev.Hr; buff[k][j][i] = Hr;
	#define GET_SHEAR_HEATING_YZ_EDGE svEdge = &jr->svYZEdge[iter++]; Hr = svEdge->svDev.Hr; buff[k][j][i] = Hr;
	#define GET_SHEAR_HEATING_XZ_EDGE svEdge = &jr->svXZEdge[iter++]; Hr = svEdge->svDev.Hr; buff[k][j][i] = Hr;

	cf = scal->dissipation_rate;

	iflag.update = 1;

	ierr = VecSet(outbuf->lbcor, 0.0); CHKERRQ(ierr);

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_SHEAR_HEATING_CENTER,  1, 0)
	INTERPOLATE_COPY(fs->DA_XY,  outbuf->lbxy,  InterpXYEdgeCorner, GET_SHEAR_HEATING_XY_EDGE, 1, 0)
	INTERPOLATE_COPY(fs->DA_YZ,  outbuf->lbyz,  InterpYZEdgeCorner, GET_SHEAR_HEATING_YZ_EDGE, 1, 0)
	INTERPOLATE_COPY(fs->DA_XZ,  outbuf->lbxz,  InterpXZEdgeCorner, GET_SHEAR_HEATING_XZ_EDGE, 1, 0)

	ierr = OutBufPut3DVecComp(outbuf, 1, 0, cf, 0.0); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteTotDispl"
PetscErrorCode PVOutWriteTotDispl(OutVec* outvec)
{

	COPY_FUNCTION_HEADER

	cf = scal->length;

	// macros to copy displacement in cell to buffer
	#define GET_DISPLX buff[k][j][i] = jr->svCell[iter++].U[0];
	#define GET_DISPLY buff[k][j][i] = jr->svCell[iter++].U[1];
	#define GET_DISPLZ buff[k][j][i] = jr->svCell[iter++].U[2];

	INTERPOLATE_COPY(jr->fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_DISPLX, 3, 0);
	INTERPOLATE_COPY(jr->fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_DISPLY, 3, 1);
	INTERPOLATE_COPY(jr->fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_DISPLZ, 3, 2);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteSHmax"
PetscErrorCode PVOutWriteSHmax(OutVec* outvec)
{
	ACCESS_FUNCTION_HEADER

	cf = scal->unit;

	// compute maximum horizontal compressive stress (SHmax) orientation
	ierr = JacResGetSHmax(jr); CHKERRQ(ierr);

	INTERPOLATE_ACCESS(jr->ldxx, InterpCenterCorner, 3, 0, 0.0)
	INTERPOLATE_ACCESS(jr->ldyy, InterpCenterCorner, 3, 1, 0.0)

	ierr = OutBufZero3DVecComp(outbuf, 3, 2); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteEHmax"
PetscErrorCode PVOutWriteEHmax(OutVec* outvec)
{
	ACCESS_FUNCTION_HEADER

	cf = scal->unit;

	// compute maximum horizontal extension rate (EHmax) orientation
	ierr = JacResGetEHmax(jr); CHKERRQ(ierr);

	INTERPOLATE_ACCESS(jr->ldxx, InterpCenterCorner, 3, 0, 0.0)
	INTERPOLATE_ACCESS(jr->ldyy, InterpCenterCorner, 3, 1, 0.0)

	ierr = OutBufZero3DVecComp(outbuf, 3, 2); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteYield"
PetscErrorCode PVOutWriteYield(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	// macro to copy yield stress to buffer

	#define GET_YIELD buff[k][j][i] = jr->svCell[iter++].yield;

	cf = scal->stress;

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_YIELD, 1, 0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteRelDIIdif"
PetscErrorCode PVOutWriteRelDIIdif(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	// macro to copy diffusion creep relative strain rate to buffer

	#define GET_DIIdif buff[k][j][i] = jr->svCell[iter++].DIIdif;

	cf = scal->unit;

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_DIIdif, 1, 0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteRelDIIdis"
PetscErrorCode PVOutWriteRelDIIdis(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	// macro to copy diffusion creep relative strain rate to buffer

	#define GET_DIIdis buff[k][j][i] = jr->svCell[iter++].DIIdis;

	cf = scal->unit;

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_DIIdis, 1, 0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteRelDIIprl"
PetscErrorCode PVOutWriteRelDIIprl(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	// macro to copy diffusion creep relative strain rate to buffer

	#define GET_DIIprl buff[k][j][i] = jr->svCell[iter++].DIIprl;

	cf = scal->unit;

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_DIIprl, 1, 0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteRelDIIpl"
PetscErrorCode PVOutWriteRelDIIpl(OutVec* outvec)
{
	COPY_FUNCTION_HEADER

	// macro to copy plastic relative strain rate to buffer

	#define GET_DIIpl buff[k][j][i] = jr->svCell[iter++].DIIpl;

	cf = scal->unit;

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_DIIpl, 1, 0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// DEBUG VECTORS
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteMomentRes"
PetscErrorCode PVOutWriteMomentRes(OutVec* outvec)
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
PetscErrorCode PVOutWriteContRes(OutVec* outvec)
{
	ACCESS_FUNCTION_HEADER

	cf  = scal->strain_rate;

	ierr = JacResCopyContinuityRes(jr, jr->gres); CHKERRQ(ierr);

	GLOBAL_TO_LOCAL(outbuf->fs->DA_CEN, jr->gc, outbuf->lbcen)

	INTERPOLATE_ACCESS(outbuf->lbcen, InterpCenterCorner, 1, 0, 0.0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWritEnergRes"
PetscErrorCode PVOutWritEnergRes(OutVec* outvec)
{
	FDSTAG      *fs;
	PetscScalar ***lbcen, ***ge;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;

	ACCESS_FUNCTION_HEADER

	cf = scal->dissipation_rate;

	fs = jr->fs;

	ierr = DMDAVecGetArray(fs->DA_CEN, outbuf->lbcen,  &lbcen); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(jr->DA_T,   jr->ge,         &ge);    CHKERRQ(ierr);

	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		lbcen[k][j][i] = ge[k][j][i];
	}
	END_STD_LOOP

	ierr = DMDAVecRestoreArray(fs->DA_CEN, outbuf->lbcen,  &lbcen); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(jr->DA_T,   jr->ge,         &ge);    CHKERRQ(ierr);

	LOCAL_TO_LOCAL(fs->DA_CEN, outbuf->lbcen)

	INTERPOLATE_ACCESS(outbuf->lbcen, InterpCenterCorner, 1, 0, 0.0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "&PVOutWriteVelocityGr"
PetscErrorCode PVOutWriteVelocityGr(OutVec* outvec)
{
	// NOTE! See warning about component ordering scheme above

	//COPY_FUNCTION_HEADER
	ACCESS_FUNCTION_HEADER
//	// macro to copy deviatoric strain rate components to buffer

	cf = scal->strain_rate;

	INTERPOLATE_ACCESS(jr->dvxdx, InterpCenterCorner, 9, 0,0.0)
	INTERPOLATE_ACCESS(jr->dvxdy, InterpXYEdgeCorner, 9, 1,0.0)
	INTERPOLATE_ACCESS(jr->dvxdz, InterpXZEdgeCorner, 9, 2,0.0)
	INTERPOLATE_ACCESS(jr->dvydx, InterpXYEdgeCorner, 9, 3,0.0)
	INTERPOLATE_ACCESS(jr->dvydy, InterpCenterCorner, 9, 4,0.0)
	INTERPOLATE_ACCESS(jr->dvydz, InterpYZEdgeCorner, 9, 5,0.0)
	INTERPOLATE_ACCESS(jr->dvzdx, InterpXZEdgeCorner, 9, 6,0.0)
	INTERPOLATE_ACCESS(jr->dvzdy, InterpYZEdgeCorner, 9, 7,0.0)
	INTERPOLATE_ACCESS(jr->dvzdz, InterpCenterCorner, 9, 8,0.0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

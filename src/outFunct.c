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
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "matFree.h"
#include "multigrid.h"
#include "matrix.h"
#include "lsolve.h"
#include "nlsolve.h"
#include "tools.h"
#include "interpolate.h"
#include "check_fdstag.h"
#include "paraViewOutBin.h"
#include "outFunct.h"
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
			mID += phRat[jj]*(PetscScalar)phases[jj].ID; \
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
#define __FUNCT__ "PVOutWriteViscTotal"
PetscErrorCode PVOutWriteViscTotal(JacRes *jr, OutBuf *outbuf)
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
PetscErrorCode PVOutWriteViscCreep(JacRes *jr, OutBuf *outbuf)
{
	COPY_FUNCTION_HEADER

	// macro to copy viscosity to buffer
	#define GET_VISC_CREEP buff[k][j][i] = jr->svCell[iter++].eta_creep;

	// output viscosity logarithm in GEO-mode
	// (negative scaling requests logarithmic output)
	if(scal->utype == _GEO_) cf = -scal->viscosity;
	else                     cf =  scal->viscosity;

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_VISC_CREEP, 1, 0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteViscoPlastic"
PetscErrorCode PVOutWriteViscoPlastic(JacRes *jr, OutBuf *outbuf)
{
	COPY_FUNCTION_HEADER

	// macro to copy viscosity to buffer
	#define GET_VISC_VISCOPLASTIC buff[k][j][i] = jr->svCell[iter++].eta_viscoplastic;

	// output viscosity logarithm in GEO-mode
	// (negative scaling requests logarithmic output)
	if(scal->utype == _GEO_) cf = -scal->viscosity;
	else                     cf =  scal->viscosity;

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_VISC_VISCOPLASTIC, 1, 0)

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

	ierr = JacResCopyVel(jr, jr->gsol); CHKERRQ(ierr);

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
	PetscErrorCode ierr;
	PetscFunctionBegin;

	if(jr->matLim.actPorePres == PETSC_TRUE)
	{
		ierr = PVOutWriteTotalPressure(jr, outbuf); CHKERRQ(ierr);
	}
	else
	{
		ierr = PVOutWriteEffPressure(jr, outbuf); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteTotalPressure"
PetscErrorCode PVOutWriteTotalPressure(JacRes *jr, OutBuf *outbuf)
{
	PetscScalar pShift, biot;

	ACCESS_FUNCTION_HEADER

	biot = jr->matLim.biot;
	cf  = scal->stress;

	// scale pressure shift
	pShift = cf*jr->pShift;

	ierr = JacResCopyPres(jr, jr->gsol); CHKERRQ(ierr);

	if(jr->actDarcy != PETSC_TRUE)
	{
		ierr = JacResGetPorePressure(jr); CHKERRQ(ierr);
	}
	else
	{
		ierr = JacResGetDarcyPorePressure(jr); CHKERRQ(ierr);
	}

	// compute total pressure
	ierr = VecWAXPY(outbuf->lbcen, biot, jr->lp_pore, jr->lp); CHKERRQ(ierr);

	INTERPOLATE_ACCESS(outbuf->lbcen, InterpCenterCorner, 1, 0, pShift)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteEffPressure"
PetscErrorCode PVOutWriteEffPressure(JacRes *jr, OutBuf *outbuf)
{
	PetscScalar pShift;

	ACCESS_FUNCTION_HEADER

	cf = scal->stress;
	iflag.use_bound = PETSC_TRUE;

	// scale pressure shift
	pShift = cf*jr->pShift;

	ierr = JacResCopyPres(jr, jr->gsol); CHKERRQ(ierr);

	INTERPOLATE_ACCESS(jr->lp, InterpCenterCorner, 1, 0, pShift)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteOverPressure"
PetscErrorCode PVOutWriteOverPressure(JacRes *jr, OutBuf *outbuf)
{
	PetscScalar pShift;

	ACCESS_FUNCTION_HEADER

	cf = scal->stress;

	// scale pressure shift
	pShift = cf*jr->pShift;

	ierr = JacResGetOverPressure(jr, outbuf->lbcen); CHKERRQ(ierr);

	INTERPOLATE_ACCESS(outbuf->lbcen, InterpCenterCorner, 1, 0, pShift)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteLithosPressure"
PetscErrorCode PVOutWriteLithosPressure(JacRes *jr, OutBuf *outbuf)
{
	ACCESS_FUNCTION_HEADER

	cf = scal->stress;

	ierr = JacResGetLithoStaticPressure(jr); CHKERRQ(ierr);

	INTERPOLATE_ACCESS(jr->lp_lithos, InterpCenterCorner, 1, 0, 0.0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWritePorePressure"
PetscErrorCode PVOutWritePorePressure(JacRes *jr, OutBuf *outbuf)
{
	ACCESS_FUNCTION_HEADER

	cf = scal->stress;

	if(jr->actDarcy != PETSC_TRUE)
	{
		ierr = JacResGetPorePressure(jr); CHKERRQ(ierr);
	}
	else
	{
		ierr = JacResGetDarcyPorePressure(jr); CHKERRQ(ierr);
	}

	INTERPOLATE_ACCESS(jr->lp_pore, InterpCenterCorner, 1, 0, 0.0)

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

	iflag.update = PETSC_TRUE;

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
PetscErrorCode PVOutWriteTotDispl(JacRes *jr, OutBuf *outbuf)
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
PetscErrorCode PVOutWriteSHmax(JacRes *jr, OutBuf *outbuf)
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
PetscErrorCode PVOutWriteEHmax(JacRes *jr, OutBuf *outbuf)
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
#define __FUNCT__ "PVOutWriteISA"
PetscErrorCode PVOutWriteISA(JacRes *jr, OutBuf *outbuf)
{
	ACCESS_FUNCTION_HEADER

	cf = scal->unit;

	// compute Infinite Strain Axis (ISA)
	ierr = JacResGetISA(jr); CHKERRQ(ierr);

	INTERPOLATE_ACCESS(jr->ldxx, InterpCenterCorner, 3, 0, 0.0)
	INTERPOLATE_ACCESS(jr->ldyy, InterpCenterCorner, 3, 1, 0.0)

	ierr = OutBufZero3DVecComp(outbuf, 3, 2); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteGOL"
PetscErrorCode PVOutWriteGOL(JacRes *jr, OutBuf *outbuf)
{
	ACCESS_FUNCTION_HEADER

	cf = scal->unit;

	// compute Grain Orientation Lag (GOL) parameter
	ierr = JacResGetGOL(jr); CHKERRQ(ierr);

	INTERPOLATE_ACCESS(jr->ldxx, InterpCenterCorner, 1, 0, 0.0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteYield"
PetscErrorCode PVOutWriteYield(JacRes *jr, OutBuf *outbuf)
{
	COPY_FUNCTION_HEADER

	// macro to copy density to buffer
	#define GET_YIELD buff[k][j][i] = jr->svCell[iter++].svDev.yield;

	cf = scal->stress;

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_YIELD, 1, 0)

	PetscFunctionReturn(0);
}

// From Darcy code
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteLiquidPressure"
PetscErrorCode PVOutWriteLiquidPressure(JacRes *jr, OutBuf *outbuf)
{
	ACCESS_FUNCTION_HEADER
	Vec            	s;

	ierr = DMCreateLocalVector(jr->fs->DA_CEN, &s); CHKERRQ(ierr);		// create local vector based on center DA, which has one DOF

	// Pull out the liquid vector from the solution vector
	VecStrideGather(jr->lPl,0,s,INSERT_VALUES);

	cf = scal->stress;
	iflag.use_bound = PETSC_TRUE;

	INTERPOLATE_ACCESS(s, InterpCenterCorner, 1, 0, 0.0)


	ierr = VecDestroy(&s);	CHKERRQ(ierr);
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWritePorosity"
PetscErrorCode PVOutWritePorosity(JacRes *jr, OutBuf *outbuf)
{
	// New come back to the previous way
	COPY_FUNCTION_HEADER
										//ACCESS_FUNCTION_HEADER
										//Vec            	s;
										//ierr = DMCreateLocalVector(jr->fs->DA_CEN, &s); CHKERRQ(ierr);		// create local vector based on center DA, which has one DOF

	cf = scal->unit;


	// macro to copy porosity to buffer
	#define GET_POROSITY buff[k][j][i] = jr->svCell[iter++].svBulk.Phi;

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_POROSITY, 1, 0)


										//// Pull out the porosity vector from the solution vector
										//VecStrideGather(jr->lPl,1,s,INSERT_VALUES);
										//
										//INTERPOLATE_ACCESS(s, InterpCenterCorner, 1, 0, 0.0)
										//
										//ierr = VecDestroy(&s);	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWritePermeability"
PetscErrorCode PVOutWritePermeability(JacRes *jr, OutBuf *outbuf)
{
	COPY_FUNCTION_HEADER

	// output permeability logarithm in GEO-mode
	// (negative scaling requests logarithmic output)
	if(scal->utype == _GEO_) cf =  -scal->permeability;
	else                     cf =  -scal->permeability; //!!


	// macro to copy permeability to buffer
	#define GET_PERMEABILITY buff[k][j][i] = jr->svCell[iter++].svBulk.Kphi;


	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_PERMEABILITY, 1, 0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteLiquidDensity"
PetscErrorCode PVOutWriteLiquidDensity(JacRes *jr, OutBuf *outbuf)
{
	COPY_FUNCTION_HEADER

	// macro to copy density to buffer
	#define GET_LIQUIDDENSITY buff[k][j][i] = jr->svCell[iter++].svBulk.Rhol;

	cf = scal->density;

	INTERPOLATE_COPY(fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_LIQUIDDENSITY, 1, 0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteLiquidVelocity"
PetscErrorCode PVOutWriteLiquidVelocity(JacRes *jr, OutBuf *outbuf)
{
	COPY_FUNCTION_HEADER

	cf = scal->velocity;

	// macros to copy liquid flow in cell to buffer
	#define GET_LIQVELX buff[k][j][i] = jr->svCell[iter++].svBulk.liquidvelocity[0];
	#define GET_LIQVELY buff[k][j][i] = jr->svCell[iter++].svBulk.liquidvelocity[1];
	#define GET_LIQVELZ buff[k][j][i] = jr->svCell[iter++].svBulk.liquidvelocity[2];

	INTERPOLATE_COPY(jr->fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_LIQVELX, 3, 0);
	INTERPOLATE_COPY(jr->fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_LIQVELY, 3, 1);
	INTERPOLATE_COPY(jr->fs->DA_CEN, outbuf->lbcen, InterpCenterCorner, GET_LIQVELZ, 3, 2);

	PetscFunctionReturn(0);
}


//---------------------------------------------------------------------------
// DEBUG VECTORS
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteJacTest"
PetscErrorCode PVOutWriteJacTest(JacRes *jr, OutBuf *outbuf)
{
	Vec diff;

	ACCESS_FUNCTION_HEADER

	cf = scal->unit;

	// create test vector
	ierr = VecDuplicate(jr->gsol, &diff);  CHKERRQ(ierr);
	ierr = VecSet(diff, 0.0);              CHKERRQ(ierr);

	// test closed-form Jacobian against finite difference approximation
	ierr = JacTest(jr, diff);

	// view difference
	ierr = JacResCopyMomentumRes(jr, diff); CHKERRQ(ierr);

	GLOBAL_TO_LOCAL(outbuf->fs->DA_X, jr->gfx, jr->lfx)
	GLOBAL_TO_LOCAL(outbuf->fs->DA_Y, jr->gfy, jr->lfy)
	GLOBAL_TO_LOCAL(outbuf->fs->DA_Z, jr->gfz, jr->lfz)

	INTERPOLATE_ACCESS(jr->lfx, InterpXFaceCorner, 3, 0, 0.0)
	INTERPOLATE_ACCESS(jr->lfy, InterpYFaceCorner, 3, 1, 0.0)
	INTERPOLATE_ACCESS(jr->lfz, InterpZFaceCorner, 3, 2, 0.0)

	ierr = VecDestroy(&diff); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
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

	INTERPOLATE_ACCESS(outbuf->lbcen, InterpCenterCorner, 1, 0, 0.0)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWritEnergRes"
PetscErrorCode PVOutWritEnergRes(JacRes *jr, OutBuf *outbuf)
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

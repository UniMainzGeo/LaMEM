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
 **    filename:   gravity.h
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
//.........................   GRAVITY FIELD    ..............................
//---------------------------------------------------------------------------
#ifndef __gravity_h__
#define __gravity_h__
 //-----------------------------------------------------------------------------

#define _MAX_Gravity_ 300000

struct FDSTAG;
struct JacRes;
struct FB;

//-----------------------------------------------------------------------------
/*
// Structure that holds gravity parameters - not yet used
struct gravityParams
{
	PetscInt     GetIt;
	PetscInt     SaveDebug,SaveVTK,SaveRef;
	PetscBool    UseNumerics, UseAnalytics;
	PetscInt     survey_nx, survey_ny;
	PetscScalar  survey_xs, survey_xm;
	PetscScalar  survey_ys, survey_ym;
	PetscScalar  survey_z ;
	PetscScalar  ReferenceDensity;
	PetscScalar  StdDev;
	PetscScalar  LithColDens[9],LithColDepth[8];
	PetscInt     num_intp,LithColNum;
	char         RefDatFile2load[MAX_PATH_LEN];
};

 */

//---------------------------------------------------------------------------
// survey context
struct GravitySurvey
{
	PetscInt     nx,ny, faphase;
	PetscScalar  xs,xm,ys,ym;
	PetscScalar  x,y,z,dx,dy;
	Vec          lvec_dg;
	Vec          gvec_dg;
	PetscScalar  coordx[_MAX_Gravity_],coordy[_MAX_Gravity_],gravref[_MAX_Gravity_];
	PetscMPIInt  rank;
	PetscInt     ComputeGravity, SaveGravity;
};

PetscErrorCode GRVSurveyCreate(FDSTAG *fs, GravitySurvey *survey, FB *fb);

PetscErrorCode GRVSurveyDestroy(GravitySurvey *survey);

PetscErrorCode GRVCompute(FDSTAG *fs, JacRes *jr, GravitySurvey *survey);

PetscErrorCode GetGravityEffectAnalytical(PetscScalar rho, PetscScalar *dx, PetscScalar *dy, PetscScalar *dz, PetscScalar *gsum);

PetscErrorCode GetGravityEffectAnalytical2(PetscScalar *x,PetscScalar *y, PetscScalar *z,PetscScalar *gsum);

PetscErrorCode SaveGravityField2VTK(Vec lvec_dg, PetscScalar *larray_coords,PetscScalar itime);


//---------------------------------------------------------------------------


#endif

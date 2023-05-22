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
//.........................   GRAVITY FIELD    ..............................
//---------------------------------------------------------------------------
#ifndef __gravity_h__
#define __gravity_h__
 //-----------------------------------------------------------------------------

struct FDSTAG;
struct JacRes;

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
	PetscInt     i,j,nx,ny;
	PetscInt     xs,xm,ys,ym;
	PetscInt     iter;
	PetscScalar  x,y,z,dx,dy;
	Vec          lvec_dg,lvec_dg2save,gvec_dg;
	PetscScalar *coord,*dg;
	PetscMPIInt  rank;

};

PetscErrorCode GRVSurveyCreate(GravitySurvey *survey);

PetscErrorCode GRVSurveyDestroy(GravitySurvey survey);

PetscErrorCode GRVCompute(FDSTAG *fs, JacRes *jr);


//---------------------------------------------------------------------------


#endif

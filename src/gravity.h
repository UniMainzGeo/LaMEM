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

//---------------------------------------------------------------------------
// survey context
typedef struct
{
	PetscInt     i,j,nx,ny;
	PetscInt     xs,xm,ys,ym;
	PetscInt     iter;
	PetscScalar  x,y,z,dx,dy;
	Vec          lvec_dg,lvec_dg2save,gvec_dg;
	PetscScalar *coord,*dg;
	PetscMPIInt  rank;

} GravitySurvey;

PetscErrorCode GRVSurveyCreate(UserCtx *user, GravitySurvey *survey);

PetscErrorCode GRVSurveyDestroy(GravitySurvey survey);

PetscErrorCode GRVCompute(FDSTAG *fs, UserCtx *user, JacRes *jr);


//---------------------------------------------------------------------------


#endif

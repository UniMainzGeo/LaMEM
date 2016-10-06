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
 **    filename:   gravity.c
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
#include "LaMEM.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "tools.h"
#include "gravity.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "GRVSurveyCreate"
PetscErrorCode GRVSurveyCreate(UserCtx *user, GravitySurvey *survey)
{
	if(user)   user = NULL;
	if(survey) survey = NULL;

/*
	PetscInt            n, i, j;
	PetscInt            i,n;

	PetscErrorCode ierr;
	PetscFunctionBegin;


	// create survey coordinates
	*survey.nx = user->GravityField.survey_nx;
	*survey.ny = user->GravityField.survey_ny;

	*survey.xs = user->GravityField.survey_xs;
	*survey.xm = user->GravityField.survey_xm;
	*survey.ys = user->GravityField.survey_ys;
	*survey.ym = user->GravityField.survey_ym;
	*survey.z  = user->GravityField.survey_z ;

	ierr = PetscPrintf(PETSC_COMM_WORLD,"#     gravity profile [m]: \n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"#     xs: %g, xm: %g, nx: %lld\n",*survey.xs, *survey.xm,(LLD)(*survey.nx)); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"#     ys: %g, ym: %g, ny: %lld\n",*survey.ys, *survey.ym,(LLD)(*survey.ny)); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"#     z : %g\n", *survey.z); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"#     calculation with SI-dimensional units \n"); CHKERRQ(ierr);


//	for(i=0;i<=user->GravityField.LithColNum;i++)
//	{
//		if (i<user->GravityField.LithColNum)
//		{
//			ierr = PetscPrintf(PETSC_COMM_WORLD,"#     LithColDepth[%lld]: %g \n",(LLD)i,user->GravityField.LithColDepth[i]);CHKERRQ(ierr);
//		}
//		ierr = PetscPrintf(PETSC_COMM_WORLD,"#     LithColDens[%lld] : %g \n",(LLD)i,user->GravityField.LithColDens[i]);CHKERRQ(ierr);
//	}


	// create local vector for survey
	ierr = VecCreateSeq(PETSC_COMM_SELF,(*survey.nx) * (*survey.ny), survey->lvec_dg); CHKERRQ(ierr);
	ierr = VecSet(survey->lvec_dg,0.0); CHKERRQ(ierr);

	// define 2D survey spacing
	*survey.dx = ( *survey.xm - *survey.xs )/( *survey.nx-1 );
	*survey.dy = ( *survey.ym - *survey.ys )/( *survey.ny-1 );

	// create a coordinate array
	ierr = PetscMalloc((size_t)2*((*survey.nx) * (*survey.ny))*sizeof(PetscScalar), survey->coord); CHKERRQ(ierr);
	n = 0;
	for ( i = 0; i < survey.nx; i++ )
	{
		for ( j = 0; j < (survey.ny*2); j+=2 )
		{
			n                   = i * (*survey.ny) * 2 + j;
			*survey->coord[n]   = (*survey.xs) + i * (*survey.dx);
			*survey->coord[n+1] = (*survey.ys) + j / 2 * (*survey.dy);
		}
	}

	// allocate global vector
	ierr = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE, (*survey.nx) * (*survey.ny), survey.gvec_dg); CHKERRQ(ierr);
	ierr = VecSet(*survey.gvec_dg, 0.0); CHKERRQ(ierr);
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, survey->rank); CHKERRQ(ierr);
*/
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "GRVSurveyDestroy"
PetscErrorCode GRVSurveyDestroy( GravitySurvey survey)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// destroy vector objects
	ierr = VecDestroy(&survey.gvec_dg);	CHKERRQ(ierr);
	ierr = VecDestroy(&survey.lvec_dg);	CHKERRQ(ierr);

	// free allocated memory of coordinates
	ierr = PetscFree(survey.coord); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "GRVCompute"
PetscErrorCode GRVCompute(FDSTAG *fs, UserCtx *user, JacRes *jr)
{

	if(fs)   fs = NULL;
	if(user) user = NULL;
	if(jr)   jr = NULL;

/*
	PetscInt      iter, i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar   x,y,z,dxh,dyh,dzh;
	PetscScalar   corners[8][8];
	GravitySurvey survey;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// create survey context
	ierr = GRVSurveyCreate( &user, &survey);

	// get ranges of cells in n-dir
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)



	ierr = VecGetArray(survey.lvec_dg,&survey.dg); CHKERRQ(ierr);
	for ( survey.i = 0; survey.i < survey.nx; survey.i++)
	{
		for ( survey.j=0; survey.j<survey.ny*2; survey.j+=2)
		{
			survey.x = survey.coord[i * survey.ny * 2 + survey.j];
			survey.y = survey_coord[i * survey.ny * 2 + survey.j + 1];

			iter = 0;
			START_STD_LOOP
			{
				// get density at the cell center
				rho = jrctx->svCell[iter].svBulk.rho;

				// get coordinate of cell center
				x = fs->dsx.ccoor[iter];
				y = fs->dsy.ccoor[iter];
				z = fs->dsz.ccoor[iter];

				// get mesh steps
				dxh = SIZE_CELL(i, sx, fs->dsx)/2.0;
				dyh = SIZE_CELL(j, sy, fs->dsy)/2.0;
				dzh = SIZE_CELL(k, sz, fs->dsz)/2.0;

				// compute cell volume
				dvol = 2.0*dxh * 2.0*dyh * 2.0*dzh;

				// get the coordinates at the cell corners
				corners[0][0] = x-dxh; corners[1][0] = y-dyh; corners[2][0] = z-dzh; //1
				corners[0][1] = x+dxh; corners[1][1] = y-dyh; corners[2][1] = z-dzh; //4
				corners[0][2] = x+dxh; corners[1][2] = y+dyh; corners[2][2] = z-dzh; //6
				corners[0][3] = x-dxh; corners[1][3] = y+dyh; corners[2][3] = z-dzh; //7
				corners[0][4] = x-dxh; corners[1][4] = y-dyh; corners[2][4] = z+dzh; //2
				corners[0][5] = x+dxh; corners[1][5] = y-dyh; corners[2][5] = z+dzh; //3
				corners[0][6] = x+dxh; corners[1][6] = y+dyh; corners[2][6] = z+dzh; //5
				corners[0][7] = x-dxh; corners[1][7] = y+dyh; corners[2][7] = z+dzh; //8

				// compute gravity contribution of particular cell at specified survey point
				ierr = GetGravityEffectNumerical(dvol,2,survey.x,survey.y,survey.z,cornervec,&result); CHKERRQ(ierr);
				survey.dg[survey.iter] = survey_dg[survey.iter] + rho*result;

				// reset local contribution to zero
				result              = 0.0;


				iter++;
			}
			END_STD_LOOP

			survey.iter++;
		}
	}
	ierr = VecRestoreArray(survey.lvec_dg,&survey.dg); CHKERRQ(ierr);

	// gather local surveys
	ierr = Sum_Vectors(PETSC_COMM_WORLD ,&survey.lvec_dg, survey.nx*survey.ny); CHKERRQ(ierr);
	ierr = VecSeq2VecMpi(rank,survey.lvec_dg ,&survey.gvec_dg); CHKERRQ(ierr);


	// save gravity field as binary data
	if(user->GravityField.SaveRef==1){
		ierr = PetscPrintf(PETSC_COMM_WORLD,"#  -- save reference data -- \n");							CHKERRQ(ierr);
		// --- create directory ---
		asprintf(&DirectoryName, "ReferenceData_%1.6lld",(LLD)itime);
		ierr = FDSTAGCreateOutputDirectory(DirectoryName); 												CHKERRQ(ierr);
		// --- create filename ---
		asprintf(&FileName,"%s/REF_Gravity.bin",DirectoryName);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"#     save file: %s \n",FileName);							CHKERRQ(ierr);
		// --- save binary file ---
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileName , FILE_MODE_WRITE, &view_out);			CHKERRQ(ierr);
		ierr = VecView(survey.gvec_dg,view_out); 														CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&view_out); 															CHKERRQ(ierr);
		// --- clear memory ---
		free(FileName);
		free(DirectoryName);
	}

	// save gravity field as vtk file
	if(user->GravityField.SaveVTK==1){
		// --- multiplying with gravitational constant ---
		ierr = VecDuplicate(survey.lvec_dg,&survey.lvec_dg2save);										CHKERRQ(ierr);
		ierr = VecCopy(lvec_survey_dg,survey.lvec_dg2save);												CHKERRQ(ierr);
		ierr = VecScale(survey.lvec_dg2save,G);															CHKERRQ(ierr);

		ierr = SaveGravityField2VTK(user,survey.lvec_dg2save, lvec_survey_coords,itime);				CHKERRQ(ierr);

		ierr = VecDestroy(&lvec_survey_dg2save);														CHKERRQ(ierr);
	}

	// get misfit
	if(user->Optimisation.GetIt==1){
		PetscPrintf(PETSC_COMM_WORLD,"#------------------------------------------\n");
		PetscPrintf(PETSC_COMM_WORLD,"#-- get gravity field misfit --------------\n");

		ierr = GetMisfit_GravityField(user,survey.gvec_dg);												CHKERRQ(ierr);
	}


	// destroy survey context
	ierr = GRVSurveyDestroy( &user, survey);
*/
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

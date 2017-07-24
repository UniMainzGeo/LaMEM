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
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "tools.h"
#include "gravity.h"
#include "parsing.h"
#include "phase.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "GRVSurveyCreate"
PetscErrorCode GRVSurveyCreate(FDSTAG *fs, GravitySurvey *survey, FB *fb)
{
	PetscInt           i;
	PetscInt           iter, j, k, nx, ny, nz, sx, sy, sz, rank;
	FILE              *fp;
	char              *fname;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

	ierr = getIntParam   (fb, _OPTIONAL_, "ComputeGravity",  &survey->ComputeGravity,      1, 1);   CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "SaveGravity",  &survey->SaveGravity,      1, 1);   CHKERRQ(ierr);



	ierr = getIntParam   (fb, _OPTIONAL_, "survey_ny",  &survey->ny,      1, 1000);   CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "survey_nx",  &survey->nx,      1, 1000);   CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "survey_xs",           &survey->xs,         1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "survey_xm",           &survey->xm,         1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "survey_ys",           &survey->ys,         1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "survey_ym",           &survey->ym,         1, 1.0); CHKERRQ(ierr);
	ierr = getScalarParam(fb, _OPTIONAL_, "survey_z",           &survey->z,         1, 1.0); CHKERRQ(ierr);

	// survey->faphase = faphase;

	ierr = PetscPrintf(PETSC_COMM_WORLD,"#     gravity profile: \n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"#     xs: %g, xm: %g, nx: %i\n",survey->xs, survey->xm,survey->nx); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"#     ys: %g, ym: %g, ny: %i\n",survey->ys, survey->ym,survey->ny); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"#     z : %g\n", survey->z); CHKERRQ(ierr);

	// get ranges of cells in n-dir
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)
	iter = 0;

	// Read in reference density if available

	asprintf(&fname, "Reference_Density_p%1.6i.tab", rank);
	fp=fopen(fname,"r");
	if(fp == NULL) PetscPrintf(PETSC_COMM_WORLD, "Can't open file Reference_Density.dat\n\n");
	START_STD_LOOP
	{
		fscanf(fp, "%lf\n",&survey->gravref[iter++]);
	}
	END_STD_LOOP
	fclose(fp);

	// create local vector for survey
	ierr = VecCreate(PETSC_COMM_WORLD, &survey->lvec_dg); CHKERRQ(ierr);
	// ierr = VecSetType(survey->lvec_dg, VECSEQ);
	ierr = VecSetSizes(survey->lvec_dg, (survey->nx) * (survey->ny), PETSC_DECIDE);
	ierr = VecSetType(survey->lvec_dg, VECSTANDARD);
	ierr = VecSet(survey->lvec_dg,0.0); CHKERRQ(ierr);

	// define 2D survey spacing
	survey->dx = ( survey->xm - survey->xs )/( survey->nx-1 );
	survey->dy = ( survey->ym - survey->ys )/( survey->ny-1 );

	// create a coordinate array
	ierr = PetscMalloc((size_t)2*((survey->nx))*sizeof(PetscScalar), survey->coordx); CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)2*((survey->nx))*sizeof(PetscScalar), survey->coordy); CHKERRQ(ierr);
	for ( i = 0; i < survey->nx; i++ )
	{
		survey->coordx[i] = (survey->xs) + i * (survey->dx);
	}

	for ( i = 0; i < survey->ny; i++ )
	{
		survey->coordy[i] = (survey->ys) + i * (survey->dy);
	}

	// allocate global vector
	ierr = VecDuplicate(survey->lvec_dg,&survey->gvec_dg);
	ierr = VecCopy(survey->lvec_dg,survey->gvec_dg);
	// ierr = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE, (survey->nx) * (survey->ny), &survey->gvec_dg); CHKERRQ(ierr);
	ierr = VecSet(survey->gvec_dg, 0.0); CHKERRQ(ierr);
	// ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &survey->rank); CHKERRQ(ierr);

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
	// ierr = PetscFree(survey.coordx); CHKERRQ(ierr);
	// ierr = PetscFree(survey.coordy); CHKERRQ(ierr);
	// ierr = PetscFree(survey.gravref); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "GRVCompute"
PetscErrorCode GRVCompute(FDSTAG *fs, JacRes *jr, GravitySurvey *survey)
{
	PetscInt      iter, siter, i, j, k, nx, ny, nz, sx, sy, sz, si, sj, rank, step;
	PetscScalar   x,y,z,dxh,dyh,dzh, xq[3], yq[3], zq[3], G, fac, sdx, sdy;
	PetscScalar   rho, result;
	PetscViewer   view_out;
	char         *fname;
	PetscScalar  *gravmerge, *gravpatch;
	Scaling      *scal;
	FILE         *fp;
	PetscLogDouble     cputime_start, cputime_end;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	PetscTime(&cputime_start);

	scal = jr->scal;

	if(scal->utype == _GEO_)
	{
		fac = scal->length*1e3;
	}
	else if(scal->utype == _SI_)
	{
		fac = scal->length;
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD,"Take care, gravity might not be scaled properly\n");
	}

	sdx = ((survey->xm*fac)-(survey->xs*fac))/(survey->nx-1);
	sdy = ((survey->ym*fac)-(survey->ys*fac))/(survey->ny-1);

	ierr = VecSet(survey->lvec_dg,0.0); CHKERRQ(ierr);
	ierr = VecSet(survey->gvec_dg, 0.0); CHKERRQ(ierr);

	// get ranges of cells in n-dir
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	siter = 0;
	ierr = VecGetArray(survey->lvec_dg,&gravpatch); CHKERRQ(ierr);
	for ( si = 0; si < survey->nx; si++)
	{
		for ( sj=0; sj<survey->ny; sj++)
		{
			iter = 0;
			START_STD_LOOP
			{
				// get density at the cell center
				rho = jr->svCell[iter++].svBulk.rho;

				/*// If we are at the free surface
				if(jr->svCell[iter].phRat == survey->faphase)
				{
					PetscPrintf(PETSC_COMM_WORLD,"WENT\n");
					rho = 0;
				}*/

				// get coordinate of cell center
				x = fs->dsx.ccoor[i-fs->dsx.pstart] * fac;
				y = fs->dsy.ccoor[j-fs->dsy.pstart] * fac;
				z = fs->dsz.ccoor[k-fs->dsz.pstart] * fac;

				// get mesh steps
				dxh = (SIZE_CELL(i, sx, fs->dsx)/2.0) * fac;
				dyh = (SIZE_CELL(j, sy, fs->dsy)/2.0) * fac;
				dzh = (SIZE_CELL(k, sz, fs->dsz)/2.0) * fac;

				survey->x=((survey->xs*fac)+(si)*sdx);
				survey->y=((survey->ys*fac)+(sj)*sdy);

				xq[0] = (x-dxh)-survey->x;
				xq[1] = (x+dxh)-survey->x;
				yq[0] = (y-dyh)-survey->y;
				yq[1] = (y+dyh)-survey->y;
				zq[0] = (z-dzh)-(survey->z*fac);
				zq[1] = (z+dzh)-(survey->z*fac);

				// compute gravity contribution of particular cell at specified survey point
				ierr = GetGravityEffectAnalytical(survey->gravref[iter]-rho*scal->density, xq, yq, zq, &result);

				// Debug values
				// PetscPrintf(PETSC_COMM_WORLD,"zm = %.20f\nxm = %.20f\nx = %.20f\ny = %.20f\nz = %.20f\ndzh = %.20f\nsx = %.20f\nsy = %.20f\nxq1 = %.20f\nxq2 = %.20f\nyq1 = %.20f\nyq2 = %.20f\nzq1 = %.20f\nzq2 = %.20f\nrho = %.20f\n",survey->z*fac,survey->xm*fac,x,y,z,dzh,survey->x,survey->y,xq[0],xq[1],yq[0],yq[1],zq[0],zq[1],rho*scal->density-survey->gravref[iter]);
				// PetscPrintf(PETSC_COMM_WORLD,"RESULT = %.20f\n\n",result);

				gravpatch[siter] = gravpatch[siter] + result;

				// reset local contribution to zero
				result  = 0.0;
			}
			END_STD_LOOP
			siter++;
		}
	}



	if (fs->dsz.nproc != 1 )
	{
		ierr = Discret1DGetColumnComm(&fs->dsz); CHKERRQ(ierr);
		ierr = VecGetArray(survey->gvec_dg,&gravmerge); CHKERRQ(ierr);
		ierr = MPI_Allreduce(gravpatch, gravmerge, (PetscMPIInt)(survey->nx*survey->ny), MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
		ierr = VecRestoreArray(survey->lvec_dg,&gravpatch); CHKERRQ(ierr);
		ierr = VecRestoreArray(survey->gvec_dg,&gravmerge); CHKERRQ(ierr);
	}
	else
	{
		ierr = VecRestoreArray(survey->lvec_dg,&gravpatch); CHKERRQ(ierr);
		ierr = VecCopy(survey->lvec_dg,survey->gvec_dg);
		// VecView(survey->gvec_dg,	PETSC_VIEWER_STDOUT_WORLD 	);
	}

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
	if(rank == 0)
	{
		ierr = VecGetArray(survey->gvec_dg,&gravmerge); CHKERRQ(ierr);
		asprintf(&fname, "Gravity_T%1.6i_p%1.6i.tab", jr->ts->istep, rank);
		fp = fopen(fname,"w");
		if(fp == NULL) SETERRQ1(PETSC_COMM_SELF, 1,"cannot open file %s", fname);
		free(fname);
		for(i=0;i<(survey->nx) * (survey->ny);i++)
		{
			fprintf(fp, "%lf\n",gravmerge[i]);
		}
		ierr = VecRestoreArray(survey->gvec_dg,&gravmerge); CHKERRQ(ierr);
		fclose(fp);
	}
	// MPI_Barrier(PETSC_COMM_WORLD);

	PetscTime(&cputime_end);
	PetscPrintf(PETSC_COMM_WORLD, " Gravity computation took %g (sec) \n", cputime_end - cputime_start);
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "GetGravityEffectAnalytical"
PetscErrorCode GetGravityEffectAnalytical(PetscScalar rho, PetscScalar *dx, PetscScalar *dy, PetscScalar *dz, PetscScalar *gsum)
{
	PetscInt		ind,i,j,k;
	PetscScalar		R111,R112,R121,R122,R211,R212,R221,R222,G;
	PetscScalar     g111,g112,g121,g122,g211,g212,g221,g222;
	PetscScalar     lsum=0.0;

	PetscFunctionBegin;

	// Analytical solution taken from Turcotte & Schubert (taken from chapter 12 (Taras))

	G=(6.6732e-11)*1e5;

	R111=pow(dx[0]*dx[0]+dy[0]*dy[0]+dz[0]*dz[0],0.5);
	R112=pow(dx[1]*dx[1]+dy[0]*dy[0]+dz[0]*dz[0],0.5);
	R121=pow(dx[0]*dx[0]+dy[1]*dy[1]+dz[0]*dz[0],0.5);
	R122=pow(dx[1]*dx[1]+dy[1]*dy[1]+dz[0]*dz[0],0.5);
	R211=pow(dx[0]*dx[0]+dy[0]*dy[0]+dz[1]*dz[1],0.5);
	R212=pow(dx[1]*dx[1]+dy[0]*dy[0]+dz[1]*dz[1],0.5);
	R221=pow(dx[0]*dx[0]+dy[1]*dy[1]+dz[1]*dz[1],0.5);
	R222=pow(dx[1]*dx[1]+dy[1]*dy[1]+dz[1]*dz[1],0.5);

	g111=-(dz[0]*atan((dx[0]*dy[0])/(dz[0]*R111))-dx[0]*log(R111+dy[0])-dy[0]*log(R111+dx[0]));
	g112=+(dz[0]*atan((dx[1]*dy[0])/(dz[0]*R112))-dx[1]*log(R112+dy[0])-dy[0]*log(R112+dx[1]));
	g121=+(dz[0]*atan((dx[0]*dy[1])/(dz[0]*R121))-dx[0]*log(R121+dy[1])-dy[1]*log(R121+dx[0]));
	g122=-(dz[0]*atan((dx[1]*dy[1])/(dz[0]*R122))-dx[1]*log(R122+dy[1])-dy[1]*log(R122+dx[1]));

	g211=+(dz[1]*atan((dx[0]*dy[0])/(dz[1]*R211))-dx[0]*log(R211+dy[0])-dy[0]*log(R211+dx[0]));
	g212=-(dz[1]*atan((dx[1]*dy[0])/(dz[1]*R212))-dx[1]*log(R212+dy[0])-dy[0]*log(R212+dx[1]));
	g221=-(dz[1]*atan((dx[0]*dy[1])/(dz[1]*R221))-dx[0]*log(R221+dy[1])-dy[1]*log(R221+dx[0]));
	g222=+(dz[1]*atan((dx[1]*dy[1])/(dz[1]*R222))-dx[1]*log(R222+dy[1])-dy[1]*log(R222+dx[1]));

	// Debug values
	// PetscPrintf(PETSC_COMM_WORLD,"dx1 = %.20f\ndy1 = %.20f\ndz1 = %.20f\nsq1 = %.20f\nsq2 = %.20f\nsq3 = %.20f\ng111 = %.20f\nR111 = %.20f\nR111full = %.20f\n\n",dx[0],dy[0],dz[0],pow(dx[0],2),pow(dy[0],2),pow(dz[0],2),g111,R111,pow(dx[0],2)+pow(dy[0],2)+pow(dz[0],2));

	lsum=G*rho*(g111+g112+g121+g122+g211+g212+g221+g222);

	*gsum = lsum;


/*
		PetscScalar		x,y,z,r;
		PetscInt		ind;
		PetscScalar		dummy = 0.0;

		// get gravity effect of CURRENT cell
		for (ind=0; ind<4; ind++) {
			x	=	survey_x	- cornervec[0][ind];
			y	=	survey_y	- cornervec[1][ind];
			z	=	survey_z	- cornervec[2][ind];
			r	=   pow((x*x+y*y+z*z),0.5);
			dummy = dummy - (y*log((x+r)) + x*log((y+r)) - z* atan2((x*y),(z*r)));
		}

		for (ind=4; ind<8; ind++) {
			x	=	survey_x	- cornervec[0][ind];
			y	=	survey_y	- cornervec[1][ind];
			z	=	survey_z	- cornervec[2][ind];
			r	=   pow((x*x+y*y+z*z),0.5);

			dummy = G * rho * (y*log((x+r)) + x*log((y+r)) - z* atan2((x*y),(z*r)));
		}



		*gsum = dummy;
*/

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "GetGravityEffectAnalytical2"
PetscErrorCode GetGravityEffectAnalytical2(PetscScalar *x,PetscScalar *y, PetscScalar *z,PetscScalar *gsum)
{
	// HANDLE WITH CARE; Probably not correct
	PetscInt		ind,i,j,k;
	PetscScalar		rijk=0.0,arg1,arg2,arg3;
	PetscScalar     ijk[] ={-1.0,1.0,1.0,-1.0,1.0,-1.0,-1.0,1.0};
	PetscScalar     lsum=0.0;

	PetscFunctionBegin;

	ind  = 0;
	for (i=0;i<2;i++)
	{
	    for (j=0;j<2;j++)
	    {
	        for (k=0;k<2;k++)
	        {
	            //rijk = sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]);
	            rijk = pow(x[i]*x[i]+y[j]*y[j]+z[k]*z[k],0.5); // factor 10 faster on mogon
	            arg1 = atan2( x[i]*y[j], z[k]*rijk );
	            if(arg1 < 0 )
	            {
	            	arg1 =arg1 + 2*M_PI;
	            }
	            arg2 = log(rijk +y[j]);
	            arg3 = log(rijk +x[i]);
	            lsum = lsum + ijk[ind] *( z[k]*arg1  - x[i]*arg2 -y[j]*arg3 );
	            // Error checking
	            if( PetscIsInfOrNanScalar(lsum) == 1) lsum = 0;
	            ind =ind+1;
	        }
	    }
	}

	*gsum = lsum;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

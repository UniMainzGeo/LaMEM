/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2011-2015, JGU Mainz, Anton Popov, Boris Kaus
 **    All rights reserved.
 **
 **    This sofware was developed at:
 **
 **         Institute of Geosciences
 **         Johannes-Gutenberg University, Mainz
 **         Johann-Joachim-Becherweg 21
 **         55128 Mainz, Germany
 **
 **    project:    LaMEM
 **    filename:   objFunct.c
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
//.................   COMPUTATION OF OBJECTIVE FUNCTION    ..................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "interpolate.h"
#include "surf.h"
#include "objFunct.h"
#include "Utils.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ObjFunctDestroy"
PetscErrorCode ObjFunctDestroy(ObjFunct *objf)
{
	PetscErrorCode ierr;
	PetscInt       k;
	PetscFunctionBegin;

	if(objf->CompMfit != PETSC_TRUE) PetscFunctionReturn(0);

	// free vectors
	for (k=0;k<_max_num_obs_; k++)
	{
		if (objf->otUse[k] == 1)
		{
			ierr = VecDestroy(&objf->obs[k]); CHKERRQ(ierr);
			ierr = VecDestroy(&objf->qul[k]); CHKERRQ(ierr);
		}
	}

	// free asprintf allocated char
	free(objf->infile);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ObjFunctCreate"
PetscErrorCode ObjFunctCreate(ObjFunct *objf, FreeSurf *surf)
{
	FDSTAG        *fs;
	PetscInt       sx, sy, nx, ny;
	PetscInt       gnx,gny;
	PetscScalar    header[2 + _max_num_obs_];
	PetscInt       Fsize,L,k,i,j,cnt,startf,startq,indf,indq;
	PetscScalar   *readbuff;
	PetscInt       buffsize;
	int            fd;
	PetscViewer    view_in;
	PetscScalar ***field,***qual;
	PetscBool      flg;
	PetscInt       useField;

	const char  *on[_max_num_obs_];			//static array of pointers
	const char           *velx_name="velx";
	const char           *vely_name="vely";
	const char           *velz_name="velz";
	const char           *topo_name="topo";
	const char           *boug_name="boug";
	const char           *isa_name="isa";
	const char           *shmax_name="shmax";

	PetscErrorCode ierr;
	PetscFunctionBegin;



	// compute misift?
	ierr = PetscOptionsGetBool( PETSC_NULL, "-objf_compute", &objf->CompMfit, &flg ); CHKERRQ(ierr);
	if(objf->CompMfit != PETSC_TRUE) PetscFunctionReturn(0);




	PetscPrintf(PETSC_COMM_WORLD,"# ------------------------------------------------------------------------\n");
	PetscPrintf(PETSC_COMM_WORLD,"# Objective function: \n");


	// set context
	objf->surf = surf;

	// define names for observational constraints
	on[_VELX_] = velx_name;
	on[_VELY_] = vely_name;
	on[_VELZ_] = velz_name;
	on[_TOPO_] = topo_name;
	on[_BOUG_] = boug_name;
	on[_ISA_]  = isa_name;
	on[_SHMAX_] = shmax_name;

	// read options
	ierr = ObjFunctReadFromOptions(objf, on); CHKERRQ(ierr);
	if(objf->otN != 0)
	{

		// Read file on single core (rank 0)
		if(ISRankZero(PETSC_COMM_WORLD))
		{
			// load observations file
			PetscPrintf(PETSC_COMM_WORLD,"# Load observations: %s \n", objf->infile);
			ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, objf->infile, FILE_MODE_READ, &view_in); CHKERRQ(ierr);
			ierr = PetscViewerBinaryGetDescriptor(view_in, &fd); CHKERRQ(ierr);

			// read (and ignore) the silent undocumented file header & size of file
			ierr = PetscBinaryRead(fd, &header, 2+_max_num_obs_, PETSC_SCALAR); CHKERRQ(ierr);
			Fsize = (PetscInt)(header[1]);

			for (k = 0; k <_max_num_obs_; k++)
			{
				useField = (PetscInt) header[k+2];
				if ( (useField == 0) && (objf->otUse[k] == 1) )
				{
					PetscPrintf(PETSC_COMM_WORLD,"# WARNING: Requested field [%s] is not found in input file -> Remove from objective function \n",on[k]);
					objf->otUse[k] = 0;
					objf->otN--;
					useField = 0;
				}
				else if ( (useField == 1) && (objf->otUse[k] == 0) )
				{
					PetscPrintf(PETSC_COMM_WORLD,"# WARNING: Found non-requested field [%s] in input file -> won't be part of objective function \n",on[k]);
					objf->otUse[k] = 0;
				}
				else if ( (useField == 1) && (objf->otUse[k] == 1) )
				{
					PetscPrintf(PETSC_COMM_WORLD,"# Observational constraint [%lld]: %s\n",(LLD)k,on[k] );
					objf->otUse[k] = 1;
				}
			}
		}


		// broadcast buffer to all cpus
		if (ISParallel(PETSC_COMM_WORLD))
		{
			objf->otUse[_max_num_obs_] = objf->otN;
			ierr = MPI_Bcast(objf->otUse, (PetscMPIInt)(_max_num_obs_+1), MPIU_2INT,(PetscMPIInt)0, PETSC_COMM_WORLD); CHKERRQ(ierr);
			objf->otN = objf->otUse[_max_num_obs_];
		}


		// access staggered grid layout
		fs = surf->jr->fs;

		// get global dimensions of the grid
		gnx = fs->dsx.tnods;
		gny = fs->dsy.tnods;

		// final size of buffer (number of observation types * gridsize)
		buffsize = 2 * objf->otN * gnx * gny;

		// number of observational constraints
		objf->ocN  = objf->otN * gnx * gny;

		// checks
		if (ISRankZero(PETSC_COMM_WORLD))
		{
			if (Fsize != buffsize) PetscPrintf(PETSC_COMM_WORLD,"# WARNING: Filesize does not correspond to the expected file size \n");
		}

		// allocate input buffer
		ierr = PetscMalloc((size_t)(buffsize)*sizeof(PetscScalar), &readbuff); CHKERRQ(ierr);



		if(ISRankZero(PETSC_COMM_WORLD))
		{
			// read observations into the buffer
			ierr = PetscBinaryRead(fd, readbuff, buffsize, PETSC_SCALAR); CHKERRQ(ierr);

			// destroy
			ierr = PetscViewerDestroy(&view_in); CHKERRQ(ierr);
		}

		// broadcast buffer to all cpus
		if (ISParallel(PETSC_COMM_WORLD))
		{
			ierr = MPI_Bcast(readbuff, (PetscMPIInt)buffsize, MPIU_SCALAR,(PetscMPIInt)0, PETSC_COMM_WORLD); CHKERRQ(ierr);
		}



		// get local output grid sizes
		ierr = DMDAGetCorners(surf->DA_SURF, &sx, &sy, NULL, &nx, &ny, NULL); CHKERRQ(ierr);
		L=fs->dsz.rank;

		// fill vectors with observations and respective quality weights
		cnt = 0;
		for (k=0; k<_max_num_obs_; k++)
		{

			startf = cnt*gnx*gny;
			startq = (cnt+objf->otN)*gnx*gny;

			if(objf->otUse[k] == 1)
			{
				// access buffer within local domain
				ierr = DMDAVecGetArray(surf->DA_SURF, objf->obs[k],  &field);  CHKERRQ(ierr);
				ierr = DMDAVecGetArray(surf->DA_SURF, objf->qul[k],  &qual );  CHKERRQ(ierr);

				// access buffer within local domain
				START_PLANE_LOOP
				{

					indf = startf + (i) + gnx*(j);
					indq = startq + (i) + gnx*(j);

					// get field
					field[L][j][i] = readbuff[indf];
					qual[L][j][i]  = readbuff[indq];

				}
				END_PLANE_LOOP

				// restore access
				ierr = DMDAVecRestoreArray(surf->DA_SURF, objf->obs[k],  &field);  CHKERRQ(ierr);
				ierr = DMDAVecRestoreArray(surf->DA_SURF, objf->qul[k],  &qual );  CHKERRQ(ierr);

				cnt++;
			}
		}

		// free buffer
		PetscFree(readbuff);
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD,"# WARNING: No fields requested -> Cannot create objective function \n");
		objf->errtot = 0.0;
		PetscPrintf(PETSC_COMM_WORLD,"# Total error = %g \n",objf->errtot);
		PetscFunctionReturn(0);
	}

	// compute error
	ierr = ObjFunctCompErr(objf);  CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD,"# ------------------------------------------------------------------------\n");

	PetscFunctionReturn(0);
}



//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ObjFunctReadFromOptions"
PetscErrorCode ObjFunctReadFromOptions(ObjFunct *objf, const char *on[])
{
	PetscErrorCode ierr;
	PetscBool      found, exists;
	PetscInt       k;
	char           otname [MAX_NAME_LEN];
	PetscFunctionBegin;


	// destroy asprintf

	// read filename of observation file
	asprintf(&objf->infile, "%s", "obs.bin");
	ierr = PetscOptionsGetString(PETSC_NULL,"-objf_obsfile", objf->infile, MAX_NAME_LEN, &found); CHKERRQ(ierr);
	if (!found){ PetscPrintf(PETSC_COMM_WORLD,"# WARNING: No filename given for observation file -> Use default: obs.bin \n"); }


	// number of fields to be read into the buffer
	objf->otN    = 0;



	for (k=0; k<_max_num_obs_; k++)
	{
		// initialize
		objf->otUse[k] = 0;

		// read output flags and allocate memory
		sprintf(otname,"-objf_use_%s",on[k]);
		ierr = PetscOptionsGetBool(NULL, otname, &exists, &found); CHKERRQ(ierr);
		if (found)
		{
			objf->otUse[k] = 1;
			objf->otN++;
			ierr = VecDuplicate(objf->surf->gtopo,&objf->obs[k]); CHKERRQ(ierr);
			ierr = VecDuplicate(objf->surf->gtopo,&objf->qul[k]); CHKERRQ(ierr);
			ierr = VecSet(objf->obs[k],0.0);  CHKERRQ(ierr);
			ierr = VecSet(objf->qul[k],0.0);  CHKERRQ(ierr);
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "VecErrSurf"
PetscErrorCode VecErrSurf(Vec mod, ObjFunct *objf, PetscInt field ,PetscScalar scal)
{
	PetscErrorCode    ierr;
	PetscScalar       ***lfield,***gfield;
	Vec               err;
	FreeSurf          *surf;
	FDSTAG            *fs;
	PetscInt          L,i,j,sx,sy,nx,ny;

	PetscFunctionBegin;

	fs   = objf->surf->jr->fs;
	surf = objf->surf;


	// create vector to store errors
	ierr = VecDuplicate(objf->obs[field],&err);  CHKERRQ(ierr);

	// initialize
	ierr = VecSet(err, 0.0);  CHKERRQ(ierr);
	objf->err[field] = 0.0;

	// get local output grid sizes
	ierr = DMDAGetCorners(surf->DA_SURF, &sx, &sy, NULL, &nx, &ny, NULL); CHKERRQ(ierr);
	L=fs->dsz.rank;

	ierr = VecSet(surf->vpatch, 0.0);  CHKERRQ(ierr);


	// access buffer within local domain
	ierr = DMDAVecGetArray(surf->DA_SURF, mod,          &lfield);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(surf->DA_SURF, surf->vpatch,  &gfield );  CHKERRQ(ierr);

	// access buffer within local domain
	START_PLANE_LOOP
	{
		// get field
		gfield[L][j][i] = lfield[L][j][i];
	}
	END_PLANE_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(surf->DA_SURF, mod,         &lfield);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->vpatch, &gfield );  CHKERRQ(ierr);

	// err = -scal*mod + obs
	ierr = VecWAXPY(err, -scal, surf->vpatch, objf->obs[field]);  CHKERRQ(ierr);

	// err^2
	ierr = VecPow(err, 2.0);  CHKERRQ(ierr);

	// err = err * qual
	ierr = VecPointwiseMult(err, err, objf->qul[field]);  CHKERRQ(ierr);

	// sum
	ierr = VecSum(err, &objf->err[field]);  CHKERRQ(ierr);

	//PetscSynchronizedPrintf(PETSC_COMM_SELF,"# rank[%lld,%lld,%lld] Error[%lld] = %g \n",(LLD)fs->dsx.rank,(LLD)fs->dsy.rank,(LLD)fs->dsz.rank, (LLD)field, objf->err[field]/(PetscScalar)fs->dsz.nproc);

	ierr = VecDestroy(&err);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ObjFunctCompErr"
PetscErrorCode ObjFunctCompErr(ObjFunct *objf)
{

	FreeSurf         *surf;
	PetscErrorCode    ierr;
	PetscInt          k;
	PetscScalar       velScal;
	PetscFunctionBegin;

	// surface object
	surf = objf->surf;

	// scaling factors
	velScal = surf->jr->scal.velocity;

	// compute weighted error of surface fields
	if (objf->otUse[_VELX_])	{ ierr = VecErrSurf(surf->vx,    objf, _VELX_,velScal);  CHKERRQ(ierr); }
	if (objf->otUse[_VELY_])	{ ierr = VecErrSurf(surf->vy,    objf, _VELY_,velScal);  CHKERRQ(ierr); }
	if (objf->otUse[_VELZ_])	{ ierr = VecErrSurf(surf->vz,    objf, _VELZ_,velScal);  CHKERRQ(ierr);	}
	if (objf->otUse[_TOPO_])	{ ierr = VecErrSurf(surf->ltopo, objf, _TOPO_,velScal);  CHKERRQ(ierr); }

/*
	// BOUGUER
	// ISA
	// SHMAX
*/


/*
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"residual_obs.bin", FILE_MODE_WRITE, &view_out);	CHKERRQ(ierr);
	ierr = VecView(objf->obs[_VELX_],view_out);														CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&view_out);															CHKERRQ(ierr);

	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"residual_mod.bin", FILE_MODE_WRITE, &view_out);	CHKERRQ(ierr);
	ierr = VecView(surf->vx,view_out);																CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&view_out);															CHKERRQ(ierr);

	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"err_mod.bin", FILE_MODE_WRITE, &view_out);	CHKERRQ(ierr);
	ierr = VecView(err_vx,view_out);																CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&view_out);															CHKERRQ(ierr);

	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"residual_obsvy.bin", FILE_MODE_WRITE, &view_out);	CHKERRQ(ierr);
	ierr = VecView(objf->obs[_VELY_],view_out);														CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&view_out);															CHKERRQ(ierr);

	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"residual_modvy.bin", FILE_MODE_WRITE, &view_out);	CHKERRQ(ierr);
	ierr = VecView(surf->vy,view_out);																CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&view_out);															CHKERRQ(ierr);

	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"err_modvy.bin", FILE_MODE_WRITE, &view_out);	CHKERRQ(ierr);
	ierr = VecView(err_vy,view_out);																CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&view_out);															CHKERRQ(ierr);
*/

	// total least squares error 
	objf->errtot = 0.0;
	for(k = 0; k < _max_num_obs_; k++)
	{
		if(objf->otUse[k] == 1)
		{
			objf->errtot += objf->err[k];
		}
	}

	objf->errtot = sqrt( objf->errtot / (PetscScalar) (objf->ocN * surf->jr->fs->dsz.nproc) ) ;
	PetscPrintf(PETSC_COMM_WORLD,"# Total error = %g \n",objf->errtot);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

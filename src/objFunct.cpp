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
//.................   COMPUTATION OF OBJECTIVE FUNCTION    ..................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "parsing.h"
#include "scaling.h"
#include "fdstag.h"
#include "surf.h"
#include "JacRes.h"
#include "tools.h"
#include "phase.h"
#include "objFunct.h"

//---------------------------------------------------------------------------
PetscErrorCode ObjFunctDestroy(ObjFunct *objf)
{
	
	PetscInt       k;
	PetscFunctionBeginUser;

	if(objf->CompMfit != PETSC_TRUE) PetscFunctionReturn(0);

	// free vectors
	for (k=0;k<_max_num_obs_; k++)
	{
		if (objf->otUse[k] == 1)
		{
			PetscCall(VecDestroy(&objf->obs[k]));
			PetscCall(VecDestroy(&objf->qul[k]));
		}
	}

	// free asprintf allocated char
	free(objf->infile);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode ObjFunctCreate(ObjFunct *objf, ModParam *IOparam, FreeSurf *surf, FB *fb)
{
	FDSTAG        *fs;
	PetscInt       sx, sy, nx, ny;
	PetscInt       gnx,gny;
	PetscScalar    header[2 + _max_num_obs_];
	PetscInt       Fsize=0,L,k,i,j,cnt,startf,startq,indf,indq;
	PetscScalar   *readbuff;
	PetscInt       buffsize;
	int            fd;
	PetscViewer    view_in;
	PetscScalar ***field,***qual;
	PetscInt       useField;

	const char  *on[_max_num_obs_];			//static array of pointers
	const char           *velx_name="velx";
	const char           *vely_name="vely";
	const char           *velz_name="velz";
	const char           *topo_name="topo";
	const char           *boug_name="boug";
	const char           *isa_name="isa";
	const char           *shmax_name="shmax";

	
	PetscFunctionBeginUser;

	// compute misift?
	if(IOparam == NULL) PetscFunctionReturn(0);

	// set flag
	objf->CompMfit = PETSC_TRUE;

	PetscPrintf(PETSC_COMM_WORLD," ------------------------------------------------------------------------\n");
	PetscPrintf(PETSC_COMM_WORLD," Objective function: \n");

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
	PetscCall(ObjFunctReadFromOptions(objf, on, fb));
	if(objf->otN != 0)
	{

		// Read file on single core (rank 0)
		if(ISRankZero(PETSC_COMM_WORLD))
		{
			// load observations file
			PetscPrintf(PETSC_COMM_WORLD," Load observations: %s \n", objf->infile);
			PetscCall(PetscViewerBinaryOpen(PETSC_COMM_SELF, objf->infile, FILE_MODE_READ, &view_in));
			PetscCall(PetscViewerBinaryGetDescriptor(view_in, &fd));

			// read (and ignore) the silent undocumented file header & size of file
			PetscCall(PetscBinaryRead(fd, &header, 2+_max_num_obs_, NULL, PETSC_SCALAR));
			Fsize = (PetscInt)(header[1]);

			for (k = 0; k <_max_num_obs_; k++)
			{
				useField = (PetscInt) header[k+2];
				if ( (useField == 0) && (objf->otUse[k] == 1) )
				{
					PetscPrintf(PETSC_COMM_WORLD," WARNING: Requested field [%s] is not found in input file -> Remove from objective function \n",on[k]);
					objf->otUse[k] = 0;
					objf->otN--;
					useField = 0;
				}
				else if ( (useField == 1) && (objf->otUse[k] == 0) )
				{
					PetscPrintf(PETSC_COMM_WORLD," WARNING: Found non-requested field [%s] in input file -> won't be part of objective function \n",on[k]);
					objf->otUse[k] = 0;
				}
				else if ( (useField == 1) && (objf->otUse[k] == 1) )
				{
					PetscPrintf(PETSC_COMM_WORLD," Observational constraint [%lld]: %s\n",(LLD)k,on[k] );
					objf->otUse[k] = 1;
				}
			}
		}


		// broadcast buffer to all cpus
		if (ISParallel(PETSC_COMM_WORLD))
		{
			objf->otUse[_max_num_obs_] = objf->otN;
			PetscCallMPI(MPI_Bcast(objf->otUse, (PetscMPIInt)(_max_num_obs_+1), MPIU_2INT,(PetscMPIInt)0, PETSC_COMM_WORLD));
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
			if (Fsize != buffsize) PetscPrintf(PETSC_COMM_WORLD," WARNING: Filesize does not correspond to the expected file size \n");
		}

		// allocate input buffer
		PetscCall(PetscMalloc((size_t)(buffsize)*sizeof(PetscScalar), &readbuff));



		if(ISRankZero(PETSC_COMM_WORLD))
		{
			// read observations into the buffer
			PetscCall(PetscBinaryRead(fd, readbuff, buffsize, NULL, PETSC_SCALAR));

			// destroy
			PetscCall(PetscViewerDestroy(&view_in));
		}

		// broadcast buffer to all cpus
		if (ISParallel(PETSC_COMM_WORLD))
		{
			PetscCallMPI(MPI_Bcast(readbuff, (PetscMPIInt)buffsize, MPIU_SCALAR,(PetscMPIInt)0, PETSC_COMM_WORLD));
		}



		// get local output grid sizes
		PetscCall(DMDAGetCorners(surf->DA_SURF, &sx, &sy, NULL, &nx, &ny, NULL));
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
				PetscCall(DMDAVecGetArray(surf->DA_SURF, objf->obs[k],  &field));
				PetscCall(DMDAVecGetArray(surf->DA_SURF, objf->qul[k],  &qual ));

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
				PetscCall(DMDAVecRestoreArray(surf->DA_SURF, objf->obs[k],  &field));
				PetscCall(DMDAVecRestoreArray(surf->DA_SURF, objf->qul[k],  &qual ));

				cnt++;
			}
		}

		// free buffer
		PetscCall(PetscFree(readbuff));
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD," WARNING: No fields requested -> Cannot create objective function \n");
		objf->errtot = 0.0;
		PetscPrintf(PETSC_COMM_WORLD," Total error = %g \n",objf->errtot);
		PetscFunctionReturn(0);
	}

	// compute error
	PetscCall(ObjFunctCompErr(objf));

	PetscPrintf(PETSC_COMM_WORLD," ------------------------------------------------------------------------\n");

	// transfer misfit value to IO structure
	IOparam->mfit = objf->errtot;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode ObjFunctReadFromOptions(ObjFunct *objf, const char *on[], FB *fb)
{
	
	PetscBool      found, exists;
	PetscInt       k;
	char           otname [_str_len_];
	PetscFunctionBeginUser;

	// read filename of observation file
	PetscCall(getStringParam(fb, _OPTIONAL_, "objf_obsfile", otname, "obs.bin"));

	// number of fields to be read into the buffer
	objf->otN    = 0;

	for (k=0; k<_max_num_obs_; k++)
	{
		// initialize
		objf->otUse[k] = 0;

		// read output flags and allocate memory
		sprintf(otname,"-objf_use_%s",on[k]);
		PetscCall(PetscOptionsGetBool(NULL, NULL, otname, &exists, &found));
		if (found)
		{
			objf->otUse[k] = 1;
			objf->otN++;
			PetscCall(VecDuplicate(objf->surf->gtopo,&objf->obs[k]));
			PetscCall(VecDuplicate(objf->surf->gtopo,&objf->qul[k]));
			PetscCall(VecSet(objf->obs[k],0.0));
			PetscCall(VecSet(objf->qul[k],0.0));
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode VecErrSurf(Vec mod, ObjFunct *objf, PetscInt field ,PetscScalar scal)
{
	PetscScalar       ***lfield,***gfield;
	Vec               err;
	FreeSurf          *surf;
	FDSTAG            *fs;
	PetscInt          L,i,j,sx,sy,nx,ny;

	PetscFunctionBeginUser;

	fs   = objf->surf->jr->fs;
	surf = objf->surf;

	// create vector to store errors
	PetscCall(VecDuplicate(objf->obs[field],&err));

	// initialize
	PetscCall(VecSet(err, 0.0));
	objf->err[field] = 0.0;

	// get local output grid sizes
	PetscCall(DMDAGetCorners(surf->DA_SURF, &sx, &sy, NULL, &nx, &ny, NULL));
	L=fs->dsz.rank;

	PetscCall(VecSet(surf->vpatch, 0.0));

	// access buffer within local domain
	PetscCall(DMDAVecGetArray(surf->DA_SURF, mod,          &lfield));
	PetscCall(DMDAVecGetArray(surf->DA_SURF, surf->vpatch,  &gfield ));

	// access buffer within local domain
	START_PLANE_LOOP
	{
		// get field
		gfield[L][j][i] = lfield[L][j][i];
	}
	END_PLANE_LOOP

	// restore access
	PetscCall(DMDAVecRestoreArray(surf->DA_SURF, mod,         &lfield));
	PetscCall(DMDAVecRestoreArray(surf->DA_SURF, surf->vpatch, &gfield ));

	// err = -scal*mod + obs
	PetscCall(VecWAXPY(err, -scal, surf->vpatch, objf->obs[field]));

	// err^2
	PetscCall(VecPow(err, 2.0));

	// err = err * qual
	PetscCall(VecPointwiseMult(err, err, objf->qul[field]));

	// sum
	PetscCall(VecSum(err, &objf->err[field]));

	PetscCall(VecDestroy(&err));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode ObjFunctCompErr(ObjFunct *objf)
{

	FreeSurf         *surf;
	PetscInt          k;
	PetscScalar       velScal;
	PetscFunctionBeginUser;

	// surface object
	surf = objf->surf;

	// scaling factors
	velScal = surf->jr->scal->velocity;

	// compute weighted error of surface fields
	if (objf->otUse[_VELX_])	{ PetscCall(VecErrSurf(surf->vx,    objf, _VELX_,velScal)); }
	if (objf->otUse[_VELY_])	{ PetscCall(VecErrSurf(surf->vy,    objf, _VELY_,velScal)); }
	if (objf->otUse[_VELZ_])	{ PetscCall(VecErrSurf(surf->vz,    objf, _VELZ_,velScal));	}
	if (objf->otUse[_TOPO_])	{ PetscCall(VecErrSurf(surf->ltopo, objf, _TOPO_,velScal)); }

/*
	// BOUGUER
	// ISA
	// SHMAX
*/


/*
	PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD,"residual_obs.bin", FILE_MODE_WRITE, &view_out));
	PetscCall(VecView(objf->obs[_VELX_],view_out));
	PetscCall(PetscViewerDestroy(&view_out));

	PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD,"residual_mod.bin", FILE_MODE_WRITE, &view_out));
	PetscCall(VecView(surf->vx,view_out));
	PetscCall(PetscViewerDestroy(&view_out));

	PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD,"err_mod.bin", FILE_MODE_WRITE, &view_out));
	PetscCall(VecView(err_vx,view_out));
	PetscCall(PetscViewerDestroy(&view_out));

	PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD,"residual_obsvy.bin", FILE_MODE_WRITE, &view_out));
	PetscCall(VecView(objf->obs[_VELY_],view_out));
	PetscCall(PetscViewerDestroy(&view_out));

	PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD,"residual_modvy.bin", FILE_MODE_WRITE, &view_out));
	PetscCall(VecView(surf->vy,view_out));
	PetscCall(PetscViewerDestroy(&view_out));

	PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD,"err_modvy.bin", FILE_MODE_WRITE, &view_out));
	PetscCall(VecView(err_vy,view_out));
	PetscCall(PetscViewerDestroy(&view_out));
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
	PetscPrintf(PETSC_COMM_WORLD," Total error = %g \n",objf->errtot);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

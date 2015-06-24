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
#define __FUNCT__ "ObjFunctClear"
PetscErrorCode ObjFunctClear(ObjFunct *objf)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// clear object
	ierr = PetscMemzero(objf, sizeof(ObjFunct)); CHKERRQ(ierr);

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


	PetscErrorCode ierr;
	PetscFunctionBegin;

	// compute misift?
	ierr = PetscOptionsGetBool( PETSC_NULL, "-objf_compute", &objf->CompMfit, &flg ); CHKERRQ(ierr);
	if(objf->CompMfit != PETSC_TRUE) PetscFunctionReturn(0);

	// set context
	objf->surf = surf;

	// define names for observational constraints
	objf->on[_VELX_] = "velx";
	objf->on[_VELY_] = "vely";
	objf->on[_VELZ_] = "velz";
	objf->on[_TOPO_] = "topo";
	objf->on[_BOUG_] = "boug";
	objf->on[_ISA_]  = "isa";
	objf->on[_SHMAX_] = "shmax";

	// read options
	ierr = ObjFunctReadFromOptions(objf); CHKERRQ(ierr);

	// access staggered grid layout
	fs = surf->jr->fs;

	// get global dimensions of the grid
	gnx = fs->dsx.tnods;
	gny = fs->dsy.tnods;

	// final size of buffer (number of observation types * gridsize)
	buffsize = 2 * objf->otN * gnx * gny;


	// number of observational constraints
	objf->ocN  = objf->otN * gnx * gny;

	// allocate input buffer
	ierr = PetscMalloc((size_t)(buffsize)*sizeof(PetscScalar), &readbuff); CHKERRQ(ierr);

	// Read file on single core (rank 0)
	if(!fs->dsx.rank)
	{
		// load observations file
		PetscPrintf(PETSC_COMM_WORLD," Load observations: %s \n", objf->infile);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, objf->infile, FILE_MODE_READ, &view_in); CHKERRQ(ierr);
		ierr = PetscViewerBinaryGetDescriptor(view_in, &fd); CHKERRQ(ierr);

		// read (and ignore) the silent undocumented file header & size of file 
		ierr = PetscBinaryRead(fd, &header, 2+_max_num_obs_, PETSC_SCALAR); CHKERRQ(ierr);
		Fsize = (PetscInt)(header[1]);
	
		for (k = 0; k <_max_num_obs_; k++)
		{
			objf->ocUse[k] = (PetscInt) header[k+2]; 
			if ( (objf->ocUse[k] == 0) && (objf->otUse[k] == PETSC_TRUE) )
			{
				PetscPrintf(PETSC_COMM_WORLD," WARNING: Requested field is not part in input file -> Remove from obective function \n");
				objf->otUse[k] = PETSC_FALSE;
			}
		}

		// checks
		if (Fsize != buffsize){
			// drop error message
		}


		// read observations into the buffer
		ierr = PetscBinaryRead(fd, readbuff, buffsize, PETSC_SCALAR); CHKERRQ(ierr);

		// destroy
		ierr = PetscViewerDestroy(&view_in); CHKERRQ(ierr);
	}

	// broadcast buffer to all cpus
	if (ISParallel(PETSC_COMM_WORLD))
	{
		ierr = MPI_Bcast(readbuff,buffsize, MPI_DOUBLE,(PetscMPIInt)0, PETSC_COMM_WORLD); CHKERRQ(ierr);
	}

	// get local output grid sizes
	ierr = DMDAGetCorners(fs->DA_COR, &sx, &sy, NULL, &nx, &ny, NULL); CHKERRQ(ierr);
	L=fs->dsz.rank;

	// fill vectors with observations and respective quality weights
	cnt = 0;
	for (k=0; k<_max_num_obs_; k++)
	{

		startf = cnt*gnx*gny;
		startq = (cnt+objf->otN)*gnx*gny;

		if(objf->otUse[k] == PETSC_TRUE)
		{
			// access buffer within local domain
			ierr = DMDAVecGetArray(surf->DA_SURF, objf->obs[k],  &field);  CHKERRQ(ierr);
			ierr = DMDAVecGetArray(surf->DA_SURF, objf->qul[k],  &qual );  CHKERRQ(ierr);

			// access buffer within local domain
			START_PLANE_LOOP
			{

				indf = startf + (sx+i) + gnx*(sy+j);
				indq = startq + (sx+i) + gnx*(sy+j);

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

	// compute error
	ierr = ObjFunctCompErr(objf);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}



//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ObjFunctReadFromOptions"
PetscErrorCode ObjFunctReadFromOptions(ObjFunct *objf)
{
	PetscErrorCode ierr;
	PetscBool      found;
	PetscInt       k;
	char           otname [MAX_NAME_LEN];
	PetscFunctionBegin;


	// destroy asprintf

	// read filename of observation file
	asprintf(&objf->infile, "%s", "obs.bin");
	ierr = PetscOptionsGetString(PETSC_NULL,"-objf_obsfile", objf->infile, MAX_NAME_LEN, &found); CHKERRQ(ierr);
	if (!found){ PetscPrintf(PETSC_COMM_WORLD," WARNING: No filename given for observation file -> Use default: obs.bin \n"); }


	// number of fields to be read into the buffer
	objf->otN    = 0;



	for (k=0; k<_max_num_obs_; k++)
	{
		// initialize
		objf->otUse[k] =PETSC_FALSE;

		// read output flags and allocate memory
		sprintf(otname,"-objf_use_%s",objf->on[k]);
		ierr = PetscOptionsGetBool(NULL, otname, &objf->otUse[k], &found); CHKERRQ(ierr);
		if (found)
		{
			objf->otN++;
			VecDuplicate(objf->surf->vx,&objf->obs[k]);
			VecDuplicate(objf->surf->vx,&objf->qul[k]);
		}
	}

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ObjFunctCompErr"
PetscErrorCode ObjFunctCompErr(ObjFunct *objf)
{

	FreeSurf         *surf;
	PetscErrorCode    ierr;
	PetscInt          L,i,j,k;
	PetscInt          sx, sy, nx, ny;
	PetscScalar      ***oc_vx,***oc_vy,***oc_vz,***oc_topo;
	PetscScalar      ***md_vx,***md_vy,***md_vz,***md_topo;
	PetscScalar      ***ql_vx,***ql_vy,***ql_vz,***ql_topo;
	PetscScalar      lenScal,velScal;
	PetscViewer       view_out;
	PetscFunctionBegin;

	// surface object
	surf = objf->surf;

	// scaling
	velScal = surf->jr->scal.velocity;
	lenScal = surf->jr->scal.length;




	// fill vectors with observations and respective quality weights


/*

	PetscPrintf(PETSC_COMM_WORLD," lenscal: %g , velscal: %g \n",lenScal,velScal);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"residual_obs.bin", FILE_MODE_WRITE, &view_out);	CHKERRQ(ierr);
	ierr = VecView(objf->obs[_VELX_],view_out);														CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&view_out);															CHKERRQ(ierr);

	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"residual_mod.bin", FILE_MODE_WRITE, &view_out);	CHKERRQ(ierr);
	ierr = VecView(surf->vx,view_out);														CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&view_out);															CHKERRQ(ierr);
*/


	// access buffer within local domain
	if(objf->otUse[_VELX_] == PETSC_TRUE)
	{
		ierr = DMDAVecGetArray(surf->DA_SURF, objf->obs[_VELX_],  &oc_vx);  CHKERRQ(ierr);
		ierr = DMDAVecGetArray(surf->DA_SURF, surf->vx,           &md_vx);  CHKERRQ(ierr);
		ierr = DMDAVecGetArray(surf->DA_SURF, objf->qul[_VELX_],  &ql_vx);  CHKERRQ(ierr);
	}

	if(objf->otUse[_VELY_] == PETSC_TRUE)
	{
		ierr = DMDAVecGetArray(surf->DA_SURF, objf->obs[_VELY_],  &oc_vy);  CHKERRQ(ierr);
		ierr = DMDAVecGetArray(surf->DA_SURF, surf->vy,           &md_vy);  CHKERRQ(ierr);
		ierr = DMDAVecGetArray(surf->DA_SURF, objf->qul[_VELY_],  &ql_vy);  CHKERRQ(ierr);
	}

	if(objf->otUse[_VELZ_] == PETSC_TRUE)
	{
		ierr = DMDAVecGetArray(surf->DA_SURF, objf->obs[_VELZ_],  &oc_vz);  CHKERRQ(ierr);
		ierr = DMDAVecGetArray(surf->DA_SURF, surf->vz,           &md_vz);  CHKERRQ(ierr);
		ierr = DMDAVecGetArray(surf->DA_SURF, objf->qul[_VELZ_],  &ql_vz);  CHKERRQ(ierr);
	}

	if(objf->otUse[_TOPO_] == PETSC_TRUE)
	{
		ierr = DMDAVecGetArray(surf->DA_SURF, objf->obs[_TOPO_],  &oc_topo);  CHKERRQ(ierr);
		ierr = DMDAVecGetArray(surf->DA_SURF, surf->ltopo,        &md_topo);  CHKERRQ(ierr);
		ierr = DMDAVecGetArray(surf->DA_SURF, objf->qul[_TOPO_],  &ql_topo);  CHKERRQ(ierr);
	}	

	if(objf->otUse[_BOUG_] == PETSC_TRUE)
	{
		PetscPrintf(PETSC_COMM_WORLD," WARNING: BOUGUER anomaly error is not yet implemented. \n");
	}

	if(objf->otUse[_ISA_] == PETSC_TRUE)
	{
		PetscPrintf(PETSC_COMM_WORLD," WARNING: ISA  error is not yet implemented. \n");
	}

	if(objf->otUse[_SHMAX_] == PETSC_TRUE)
	{
		PetscPrintf(PETSC_COMM_WORLD," WARNING: SHMAX  error is not yet implemented. \n");
	}

	// initialize errors
	for (k=0; k<_max_num_obs_; k++)
	{
		objf->err[k] = 0.0;
	}

	// access buffer within local domain
	ierr = DMDAGetCorners(surf->jr->fs->DA_COR, &sx, &sy, NULL, &nx, &ny, NULL); CHKERRQ(ierr);
	L = surf->jr->fs->dsz.rank;

	START_PLANE_LOOP
	{
		// compute weighted least squares
		// quality = (quality'/sigma)^2, where sigma is stdv and quality' goes from 0..1
		if(objf->otUse[_VELX_] == PETSC_TRUE)
		{
			objf->err[_VELX_] += pow((oc_vx[L][j][i] - md_vx[L][j][i]*velScal),2) * ql_vx[L][j][i];
			if(ql_vx[L][j][i] == 0.0) objf->ocN --;
		}

		if(objf->otUse[_VELY_] == PETSC_TRUE)
		{
			objf->err[_VELY_] += pow((oc_vy[L][j][i] - md_vy[L][j][i]*velScal),2) * ql_vy[L][j][i];
			if(ql_vy[L][j][i] == 0.0) objf->ocN --;
		}
		if(objf->otUse[_VELZ_] == PETSC_TRUE)
		{
			objf->err[_VELZ_] += pow((oc_vz[L][j][i] - md_vz[L][j][i]*velScal),2) * ql_vz[L][j][i];
			if(ql_vz[L][j][i] == 0.0) objf->ocN --;
		}
		if(objf->otUse[_TOPO_] == PETSC_TRUE)
		{
			objf->err[_TOPO_] += pow((oc_topo[L][j][i] - md_topo[L][j][i]*lenScal),2) * ql_topo[L][j][i];
			if(ql_topo[L][j][i] == 0.0) objf->ocN --;
		}
		
		// other fields go here
	}
	END_PLANE_LOOP

	// restore access
	if(objf->otUse[_VELX_] == PETSC_TRUE)
	{
		ierr = DMDAVecRestoreArray(surf->DA_SURF, objf->obs[_VELX_],  &oc_vx);  CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->vx,           &md_vx);  CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(surf->DA_SURF, objf->qul[_VELX_],  &ql_vx);  CHKERRQ(ierr);
	}

	if(objf->otUse[_VELY_] == PETSC_TRUE)
	{
		ierr = DMDAVecRestoreArray(surf->DA_SURF, objf->obs[_VELY_],  &oc_vy);  CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->vy,           &md_vy);  CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(surf->DA_SURF, objf->qul[_VELY_],  &ql_vy);  CHKERRQ(ierr);
	}

	if(objf->otUse[_VELZ_] == PETSC_TRUE)
	{
		ierr = DMDAVecRestoreArray(surf->DA_SURF, objf->obs[_VELZ_],  &oc_vz);  CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->vz,           &md_vz);  CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(surf->DA_SURF, objf->qul[_VELZ_],  &ql_vz);  CHKERRQ(ierr);
	}

	if(objf->otUse[_TOPO_] == PETSC_TRUE)
	{
		ierr = DMDAVecGetArray(surf->DA_SURF, objf->obs[_TOPO_],  &oc_topo);  CHKERRQ(ierr);
		ierr = DMDAVecGetArray(surf->DA_SURF, surf->ltopo,        &md_topo);  CHKERRQ(ierr);
		ierr = DMDAVecGetArray(surf->DA_SURF, objf->qul[_TOPO_],  &ql_topo);  CHKERRQ(ierr);
	}

	if(objf->otUse[_BOUG_] == PETSC_TRUE)
	{
		PetscPrintf(PETSC_COMM_WORLD," WARNING: BOUGUER anomaly error is not yet implemented. \n");
	}

	if(objf->otUse[_ISA_] == PETSC_TRUE)
	{
		PetscPrintf(PETSC_COMM_WORLD," WARNING: ISA  error is not yet implemented. \n");
	}

	if(objf->otUse[_SHMAX_] == PETSC_TRUE)
	{
		PetscPrintf(PETSC_COMM_WORLD," WARNING: SHMAX  error is not yet implemented. \n");
	}

	// MPI_Allreduce of errors
	if (ISParallel(PETSC_COMM_WORLD))
	{
		ierr = MPI_Allreduce(objf->err,objf->err,objf->otN,MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);  CHKERRQ(ierr);
	}
	// total least squares error 
	objf->errtot = 0.0;
	for (k=0; k<_max_num_obs_; k++)
	{
		if(objf->otUse[k] == PETSC_TRUE)
		{
			objf->errtot += objf->err[k];
		}
	}

	objf->errtot = objf->errtot / objf->ocN / objf->surf->jr->fs->dsz.nproc;
	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------


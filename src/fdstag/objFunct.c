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
#define __FUNCT__ "ObjFunctClean"
PetscErrorCode ObjFunctClean(ObjFunct *objf)
{
	PetscErrorCode ierr;
	PetscInt       k;
	PetscFunctionBegin;

	if(objf->CompMfit != PETSC_TRUE) PetscFunctionReturn(0);

	// free vectors
	for (k=0;k<_max_num_obs_; k++)
	{
		if (objf->otUse[k] == PETSC_TRUE)
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
		PetscPrintf(PETSC_COMM_WORLD,"# Load observations: %s \n", objf->infile);
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
				PetscPrintf(PETSC_COMM_WORLD,"# WARNING: Requested field is not part in input file -> Remove from obective function \n");
				objf->otUse[k] = PETSC_FALSE;
			}
			if (objf->ocUse[k])
				PetscPrintf(PETSC_COMM_WORLD,"# Observational constraint [%lld]: %s\n",(LLD)k,on[k] );
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

	PetscPrintf(PETSC_COMM_WORLD,"# ------------------------------------------------------------------------\n");

	PetscFunctionReturn(0);
}



//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ObjFunctReadFromOptions"
PetscErrorCode ObjFunctReadFromOptions(ObjFunct *objf, const char *on[])
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
	if (!found){ PetscPrintf(PETSC_COMM_WORLD,"# WARNING: No filename given for observation file -> Use default: obs.bin \n"); }


	// number of fields to be read into the buffer
	objf->otN    = 0;



	for (k=0; k<_max_num_obs_; k++)
	{
		// initialize
		objf->otUse[k] =PETSC_FALSE;

		// read output flags and allocate memory
		sprintf(otname,"-objf_use_%s",on[k]);
		ierr = PetscOptionsGetBool(NULL, otname, &objf->otUse[k], &found); CHKERRQ(ierr);
		if (found)
		{
			objf->otUse[k] = PETSC_TRUE;
			objf->otN++;
			ierr = VecDuplicate(objf->surf->vx,&objf->obs[k]); CHKERRQ(ierr);
			ierr = VecDuplicate(objf->surf->vx,&objf->qul[k]); CHKERRQ(ierr);
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "VecErrSurf"
PetscErrorCode VecErrSurf(Vec *err, Vec mod, ObjFunct *objf, PetscInt field ,PetscScalar scal)
{
	PetscErrorCode    ierr;

	PetscFunctionBegin;
	ierr = VecDuplicate(mod,err);  CHKERRQ(ierr);

	//err = -scal*mod + obs
	ierr = VecWAXPY(*err, -scal, mod, objf->obs[field]);  CHKERRQ(ierr);

	//err^2
	ierr = VecPow(*err, 2.0);  CHKERRQ(ierr);

	// err = err * qual
	ierr = VecPointwiseMult(*err, *err, objf->qul[field]);  CHKERRQ(ierr);

	// sum
	ierr = VecSum(*err, &objf->err[field]);  CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD,"# Error[%lld] = %g \n",(LLD)field,objf->err[field]);

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
	PetscScalar      lenScal,velScal;
	PetscViewer       view_out;
	Vec              err_vx,err_vy;
	PetscFunctionBegin;

	// surface object
	surf = objf->surf;

	// scaling factors
	velScal = surf->jr->scal.velocity;
	lenScal = surf->jr->scal.length;


	// compute weighted error of surface fields
	ierr = VecErrSurf(&err_vx, surf->vx, objf, _VELX_,velScal);  CHKERRQ(ierr);
	ierr = VecErrSurf(&err_vy, surf->vy, objf, _VELY_,velScal);  CHKERRQ(ierr);
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

	objf->errtot = sqrt( objf->errtot / (PetscScalar) (objf->ocN * objf->surf->jr->fs->dsz.nproc) ) ;
	PetscPrintf(PETSC_COMM_WORLD,"# Total error = %g \n",objf->errtot);

	// Destroy error vectors
	ierr = VecDestroy(&err_vx);  CHKERRQ(ierr);
	ierr = VecDestroy(&err_vy);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
//........................  BREAKPOINT ROUTINES   ...........................
//---------------------------------------------------------------------------
// Routines to handle restart of simulation

#include "LaMEM.h"
#include "Utils.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "multigrid.h"
#include "matrix.h"
#include "lsolve.h"
#include "nlsolve.h"
#include "interpolate.h"
#include "surf.h"
#include "advect.h"
#include "marker.h"
#include "break.h"

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BreakCheck"
PetscErrorCode BreakCheck(UserCtx *user)
{
	// check if breakpoints exist and restart new simulation is not - very useful for chain jobs

	PetscMPIInt iproc, nproc;
	FILE        *fp;
	char        *fname;
	PetscInt    res, gres;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	if(!user->restart) PetscFunctionReturn(0);

	res = 1;

	// check if directory exists
	struct stat s;
	int err = stat("./Breakpoint", &s);

	// directory missing
	if (err==-1) {
		res = 0;
	}
	// directory exists
	else
	{
		// check if breakpoint files are available
		ierr = MPI_Comm_rank(MPI_COMM_WORLD, &iproc); CHKERRQ(ierr);

		asprintf(&fname, "./Breakpoint/Breakpoint_info.%lld.out",(LLD)iproc);
		fp = fopen(fname, "r" );

		if(!fp) res = 0;
		else    fclose(fp);

		free(fname);
	}

	// check corrupted results on processors
	if(ISParallel(PETSC_COMM_WORLD))
	{
		ierr = MPI_Allreduce(&res, &gres, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
	}
	else
	{
		gres = res;
	}

	ierr = MPI_Comm_size(PETSC_COMM_WORLD, &nproc); CHKERRQ(ierr);

	if ((gres < nproc) && (gres > 0))
	{
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Number of breakpoint files does not match number of cpus!");
	}
	else if (gres == 0)
	{
		user->restart = 0;
		PetscPrintf(PETSC_COMM_WORLD,"-------------------------------------------------------------------------- \n");
		PetscPrintf(PETSC_COMM_WORLD," No breakpoints detected -> starting new simulation \n");
		PetscPrintf(PETSC_COMM_WORLD,"-------------------------------------------------------------------------- \n");
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BreakWrite"
PetscErrorCode BreakWrite(UserCtx *user, AdvCtx *actx, FreeSurf *surf, JacType jtype)
{
	// staggered grid
	FDSTAG         *fs;
	JacRes         *jr;
	FILE           *fp;
	char           *fname;
	PetscInt        n;
	PetscInt       initGuessFlag, jtypeFlag, sflatFlag;
	PetscLogDouble tstart, tend;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	PetscPrintf(PETSC_COMM_WORLD,"******************************************** \n");
	PetscPrintf(PETSC_COMM_WORLD," Writing Breakpoint files... \n");

	// start time counter
	PetscTime(&tstart);

	// initialize context
	jr = actx->jr;
	fs = actx->fs;

	// create directory
	ierr = LaMEMCreateOutputDirectory("./Breakpoint"); CHKERRQ(ierr);

	//============================================================
	//   MARKERS
	//============================================================
	// compile file name
	asprintf(&fname, "./Breakpoint/Breakpoint_mark.%lld.out",(LLD)actx->iproc);

	// open file for binary output
	fp = fopen(fname, "w" );

	// write
	fwrite(&actx->nummark, sizeof(PetscInt), 1, fp);
	fwrite(&actx->markcap, sizeof(PetscInt), 1, fp);
	fwrite( actx->markers, sizeof(Marker)*(size_t)actx->markcap, 1, fp);

	// close and free memory
	free(fname);
	fclose(fp);

	//============================================================
	//   GRID (only for background strainrate)
	//============================================================
	if(jr->bc->ExxAct == PETSC_TRUE || jr->bc->EyyAct == PETSC_TRUE)
	{
		// compile file name
		asprintf(&fname, "./Breakpoint/Breakpoint_grid.%lld.out",(LLD)actx->iproc);

		// open file for binary output
		fp = fopen(fname, "w" );

		//----------------------------------
		// Domain information
		//----------------------------------
		fwrite(&user->W      , sizeof(PetscScalar), 1, fp);
		fwrite(&user->L      , sizeof(PetscScalar), 1, fp);
		fwrite(&user->H      , sizeof(PetscScalar), 1, fp);
		fwrite(&user->x_left , sizeof(PetscScalar), 1, fp);
		fwrite(&user->y_front, sizeof(PetscScalar), 1, fp);
		fwrite(&user->z_bot  , sizeof(PetscScalar), 1, fp);

		//----------------------------------
		// FDSTAG
		//----------------------------------
		BreakWriteDiscret1D(fp, fs->dsx, fs->msx);
		BreakWriteDiscret1D(fp, fs->dsy, fs->msy);
		BreakWriteDiscret1D(fp, fs->dsz, fs->msz);

		// close and free memory
		free(fname);
		fclose(fp);
	}

	//============================================================
	//   FREE SURFACE
	//============================================================
	if (surf->UseFreeSurf==PETSC_TRUE)
	{
		// compile file name
		asprintf(&fname, "./Breakpoint/Breakpoint_surf.%lld.out",(LLD)actx->iproc);

		// open file for binary output
		fp = fopen(fname, "w" );

		// set local size of vector
		n  = fs->dsx.nnods * fs->dsy.nnods;

		// write vector
		ierr = BreakWriteVec(fp, surf->gtopo, n); CHKERRQ(ierr);

		// write flat-flag value to file
		if (surf->flat==PETSC_TRUE) sflatFlag = 1;
		else                        sflatFlag = 0;

		fwrite(&sflatFlag , sizeof(PetscInt), 1, fp);

		// close and free memory
		free(fname);
		fclose(fp);
	}

	//============================================================
	//   GSOL - Solution vectors
	//============================================================
	// compile file name
	asprintf(&fname, "./Breakpoint/Breakpoint_gsol.%lld.out",(LLD)actx->iproc);

	// open file for binary output
	fp = fopen(fname, "w" );

	//----------------------------------
	// GSOL - Stokes
	//----------------------------------
	// set local size of vector
	n  = actx->fs->dof.ln;

	// write vector
	ierr = BreakWriteVec(fp, jr->gsol, n); CHKERRQ(ierr);

	//----------------------------------
	// GT - Temperature
	//----------------------------------
	// set local size of vector
	n  = fs->dsx.ncels * fs->dsy.ncels * fs->dsz.ncels;

	// write vector
	ierr = BreakWriteVec(fp, jr->gT, n); CHKERRQ(ierr);

	// close and free memory
	free(fname);
	fclose(fp);

	//============================================================
	//   INFO
	//============================================================
	// compile file name
	asprintf(&fname, "./Breakpoint/Breakpoint_info.%lld.out",(LLD)actx->iproc);

	// open file for binary output
	fp = fopen(fname, "w" );

	//----------------------------------
	// Time Stepping
	//----------------------------------
	fwrite(&jr->ts.istep , sizeof(PetscInt   ), 1, fp);
	fwrite(&jr->ts.dt    , sizeof(PetscScalar), 1, fp);
	fwrite(&jr->ts.time  , sizeof(PetscScalar), 1, fp);

	//----------------------------------
	// Boundary Conditions
	//----------------------------------
	// pushing block center coordinates
	if (user->AddPushing)
	{
		fwrite(&user->Pushing.x_center_block , sizeof(PetscScalar), 1, fp);
		fwrite(&user->Pushing.y_center_block , sizeof(PetscScalar), 1, fp);
		fwrite(&user->Pushing.z_center_block , sizeof(PetscScalar), 1, fp);
	}

	//----------------------------------
	// Solver options
	//----------------------------------
	// store type of solver
	if      (jtype == _PICARD_) jtypeFlag = 0;
	else if (jtype == _MF_    ) jtypeFlag = 1;
	else if (jtype == _MFFD_  ) jtypeFlag = 2;
	else SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Incorrect type of jtype: %lld",(LLD)jtypeFlag);

	fwrite(&jtypeFlag , sizeof(PetscInt), 1, fp);

	// store init guess flag
	if (jr->matLim.initGuessFlg) initGuessFlag = 1;
	else                         initGuessFlag = 0;

	fwrite(&initGuessFlag , sizeof(PetscInt), 1, fp);

	//----------------------------------
	// Other
	//----------------------------------
	// store breakpoint number
	fwrite(&user->break_point_number , sizeof(PetscInt), 1, fp);

	// close and free memory
	free(fname);
	fclose(fp);

	//============================================================
	//   END WRITE BREAKPOINTS
	//============================================================

	// all other ranks should wait
	ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

	// end time counter
	PetscTime(&tend);

	PetscPrintf(PETSC_COMM_WORLD," Finished writing breakpoint files in %g s\n",tend-tstart);
	PetscPrintf(PETSC_COMM_WORLD,"******************************************** \n");
	PetscPrintf(PETSC_COMM_WORLD," \n");

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BreakRead"
PetscErrorCode BreakRead(UserCtx *user, AdvCtx *actx, JacType *jtype)
{
	JacRes      *jr;
	FDSTAG      *fs;
	FILE        *fp;
	char        *fname;
	PetscInt     n;
	JacType     j;
	PetscInt    initGuessFlag, jtypeFlag;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// set context
	jr = actx->jr;
	fs = actx->fs;

	//============================================================
	//   Solution vectors
	//============================================================
	// compile file name
	asprintf(&fname, "./Breakpoint/Breakpoint_gsol.%lld.out",(LLD)actx->iproc);

	// open file for reading
	fp = fopen(fname, "r" );

	//----------------------------------
	// GSOL
	//----------------------------------
	// set local size of vector
	n  = fs->dof.ln;

	// read vector
	ierr = BreakReadVec(fp, jr->gsol, n); CHKERRQ(ierr);

	//----------------------------------
	// GT
	//----------------------------------
	// set local size of vector
	n  = fs->dsx.ncels * fs->dsy.ncels * fs->dsz.ncels;

	// read vector
	ierr = BreakReadVec(fp, jr->gT, n); CHKERRQ(ierr);

	// close and free memory
	free(fname);
	fclose(fp);

	//============================================================
	//   INFO
	//============================================================
	// compile file name and open file
	asprintf(&fname, "./Breakpoint/Breakpoint_info.%lld.out",(LLD)actx->iproc);

	// open file for reading
	fp = fopen(fname, "r" );

	//----------------------------------
	// Time Stepping
	//----------------------------------
	fread(&jr->ts.istep , sizeof(PetscInt   ), 1, fp);
	fread(&jr->ts.dt    , sizeof(PetscScalar), 1, fp);
	fread(&jr->ts.time  , sizeof(PetscScalar), 1, fp);

	//----------------------------------
	// Boundary Conditions
	//----------------------------------
	// pushing block center coordinates
	if (user->AddPushing)
	{
		fread(&user->Pushing.x_center_block , sizeof(PetscScalar), 1, fp);
		fread(&user->Pushing.y_center_block , sizeof(PetscScalar), 1, fp);
		fread(&user->Pushing.z_center_block , sizeof(PetscScalar), 1, fp);
	}

	//----------------------------------
	// Solver options
	//----------------------------------
	// read type of solver
	fread(&jtypeFlag , sizeof(PetscInt), 1, fp);

	if      (jtypeFlag == 0) j = _PICARD_;
	else if (jtypeFlag == 1) j = _MF_    ;
	else if (jtypeFlag == 2) j = _MFFD_  ;
	else SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER, "Incorrect type of jtype: %lld",(LLD)jtypeFlag);

	*jtype = j;

	// read init guess flag
	fread(&initGuessFlag , sizeof(PetscInt), 1, fp);

	if (initGuessFlag) jr->matLim.initGuessFlg = PETSC_TRUE;
	else               jr->matLim.initGuessFlg = PETSC_FALSE;

	//----------------------------------
	// Other
	//----------------------------------
	// read breakpoint number
	fread(&user->break_point_number , sizeof(PetscInt), 1, fp);

	// close and free memory
	free(fname);
	fclose(fp);

	//============================================================
	//   END READ BREAKPOINTS
	//============================================================

	// all other ranks should wait
	ierr = MPI_Barrier(PETSC_COMM_WORLD);    CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD," \n");
	PetscPrintf(PETSC_COMM_WORLD,"RESTART from Time step = %lld \n",(LLD)jr->ts.istep);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BreakReadGrid"
PetscErrorCode BreakReadGrid(UserCtx *user, FDSTAG *fs)
{
	FILE        *fp;
	char        *fname;
	PetscMPIInt iproc;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	//----------------------------------
	// GRID information
	//----------------------------------
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &iproc); CHKERRQ(ierr);

	// compile file name - fdstag context
	asprintf(&fname, "./Breakpoint/Breakpoint_grid.%lld.out",(LLD)iproc);

	// open file for binary output
	fp = fopen(fname, "r" );

	//----------------------------------
	// Domain information
	//----------------------------------
	fread(&user->W      , sizeof(PetscScalar), 1, fp);
	fread(&user->L      , sizeof(PetscScalar), 1, fp);
	fread(&user->H      , sizeof(PetscScalar), 1, fp);
	fread(&user->x_left , sizeof(PetscScalar), 1, fp);
	fread(&user->y_front, sizeof(PetscScalar), 1, fp);
	fread(&user->z_bot  , sizeof(PetscScalar), 1, fp);

	//----------------------------------
	// FDSTAG
	//----------------------------------
	BreakReadDiscret1D(fp, fs->dsx, fs->msx);
	BreakReadDiscret1D(fp, fs->dsy, fs->msy);
	BreakReadDiscret1D(fp, fs->dsz, fs->msz);

	// close and free memory
	free(fname);
	fclose(fp);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BreakReadSurf"
PetscErrorCode BreakReadSurf(FDSTAG *fs, FreeSurf *surf)
{
	FILE        *fp;
	char        *fname;
	PetscMPIInt iproc;
	PetscInt    n;
	PetscInt    sflatFlag;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	//----------------------------------
	// FREE SURFACE
	//----------------------------------
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &iproc); CHKERRQ(ierr);

	// compile file name
	asprintf(&fname, "./Breakpoint/Breakpoint_surf.%lld.out",(LLD)iproc);

	// open file for reading
	fp = fopen(fname, "r" );

	// set local size of vector
	n  = fs->dsx.nnods * fs->dsy.nnods;

	// read vector
	ierr = BreakReadVec(fp, surf->gtopo, n); CHKERRQ(ierr);

	// read flat-flag value to file
	fread(&sflatFlag , sizeof(PetscInt), 1, fp);

	if (sflatFlag==1) surf->flat = PETSC_TRUE;
	else              surf->flat = PETSC_FALSE;

	// close and free memory
	free(fname);
	fclose(fp);

	// correct values for ltopo
	GLOBAL_TO_LOCAL(surf->DA_SURF, surf->gtopo, surf->ltopo);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BreakReadMark"
PetscErrorCode BreakReadMark(AdvCtx *actx)
{
	FILE        *fp;
	char        *fname;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	//----------------------------------
	// MARKERS
	//----------------------------------
	// compile file name
	asprintf(&fname, "./Breakpoint/Breakpoint_mark.%lld.out",(LLD)actx->iproc);

	// open file for binary output
	fp = fopen(fname, "r" );

	// read header
	fread(&actx->nummark, sizeof(PetscInt), 1, fp);
	fread(&actx->markcap, sizeof(PetscInt), 1, fp);

	// allocate memory for markers
	ierr = PetscMalloc((size_t)actx->markcap*sizeof(Marker), &actx->markers); CHKERRQ(ierr);
	ierr = PetscMemzero(actx->markers, (size_t)actx->markcap*sizeof(Marker)); CHKERRQ(ierr);

	// allocate memory for host cell numbers
	ierr = PetscMalloc((size_t)actx->markcap*sizeof(PetscInt), &actx->cellnum); CHKERRQ(ierr);
	ierr = PetscMemzero(actx->cellnum, (size_t)actx->markcap*sizeof(PetscInt)); CHKERRQ(ierr);

	// allocate memory for id marker arranging per cell
	ierr = PetscMalloc((size_t)actx->markcap*sizeof(PetscInt), &actx->markind); CHKERRQ(ierr);
	ierr = PetscMemzero(actx->markind, (size_t)actx->markcap*sizeof(PetscInt)); CHKERRQ(ierr);

	// read markers
	fread( actx->markers, sizeof(Marker)*(size_t)actx->markcap, 1, fp);

	// close and free memory
	free(fname);
	fclose(fp);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BreakWriteVec"
PetscErrorCode BreakWriteVec(FILE *fp, Vec x, PetscInt n)
{
	PetscScalar *xarray;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get vector array
	ierr = VecGetArray(x, &xarray); CHKERRQ(ierr);

	// write to file
	fwrite(xarray, sizeof(PetscScalar),(size_t)n, fp);

	// restore vector array
	ierr = VecRestoreArray(x, &xarray); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BreakReadVec"
PetscErrorCode BreakReadVec(FILE *fp, Vec x, PetscInt n)
{
	PetscScalar *xarray;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get vector array
	ierr = VecGetArray(x, &xarray); CHKERRQ(ierr);

	// read array
	fread(xarray, sizeof(PetscScalar),(size_t)n, fp);

	// restore vector array
	ierr = VecRestoreArray(x, &xarray); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
void BreakWriteDiscret1D(FILE *fp, Discret1D ds, MeshSeg1D ms)
{
	//----------------------------------
	// FDSTAG DISCRET 1D
	//----------------------------------
	fwrite(&ds.nproc  , sizeof(PetscInt),    1                 , fp);
	fwrite(&ds.rank   , sizeof(PetscInt),    1                 , fp);
	fwrite(ds.starts  , sizeof(PetscInt),    (size_t)ds.nproc+1, fp);
	fwrite(&ds.pstart , sizeof(PetscInt),    1                 , fp);
	fwrite(ds.ncoor   , sizeof(PetscScalar), (size_t)ds.nnods+1, fp);
	fwrite(ds.ccoor   , sizeof(PetscScalar), (size_t)ds.ncels+1, fp);
	fwrite(ds.nbuff   , sizeof(PetscScalar), (size_t)ds.bufsz  , fp);
	fwrite(ds.cbuff   , sizeof(PetscScalar), (size_t)ds.ncels+2, fp);
	fwrite(&ds.bufsz  , sizeof(PetscInt),    1                 , fp);
	fwrite(&ds.tnods  , sizeof(PetscInt),    1                 , fp);
	fwrite(&ds.tcels  , sizeof(PetscInt),    1                 , fp);
	fwrite(&ds.nnods  , sizeof(PetscInt),    1                 , fp);
	fwrite(&ds.ncels  , sizeof(PetscInt),    1                 , fp);
	fwrite(&ds.grprev , sizeof(PetscInt),    1                 , fp);
	fwrite(&ds.grnext , sizeof(PetscInt),    1                 , fp);
	fwrite(&ds.color  , sizeof(PetscInt),    1                 , fp);
	fwrite(&ds.h_uni  , sizeof(PetscScalar), 1                 , fp);
	fwrite(&ds.h_min  , sizeof(PetscScalar), 1                 , fp);
	fwrite(&ds.h_max  , sizeof(PetscScalar), 1                 , fp);

	//----------------------------------
	// FDSTAG STRETCH 1D
	//----------------------------------
	fwrite(&ms.nsegs  , sizeof(PetscInt),    1                 , fp);
	fwrite(ms.istart  , sizeof(PetscInt),    (size_t)ms.nsegs+1, fp);
	fwrite(ms.xstart  , sizeof(PetscScalar), (size_t)ms.nsegs+1, fp);
}
//---------------------------------------------------------------------------
void BreakReadDiscret1D(FILE *fp, Discret1D ds, MeshSeg1D ms)
{
	//----------------------------------
	// FDSTAG DISCRET 1D
	//----------------------------------
	fread(&ds.nproc  , sizeof(PetscInt),    1                 , fp);
	fread(&ds.rank   , sizeof(PetscInt),    1                 , fp);
	fread(ds.starts  , sizeof(PetscInt),    (size_t)ds.nproc+1, fp);
	fread(&ds.pstart , sizeof(PetscInt),    1                 , fp);
	fread(ds.ncoor   , sizeof(PetscScalar), (size_t)ds.nnods+1, fp);
	fread(ds.ccoor   , sizeof(PetscScalar), (size_t)ds.ncels+1, fp);
	fread(ds.nbuff   , sizeof(PetscScalar), (size_t)ds.bufsz  , fp);
	fread(ds.cbuff   , sizeof(PetscScalar), (size_t)ds.ncels+2, fp);
	fread(&ds.bufsz  , sizeof(PetscInt),    1                 , fp);
	fread(&ds.tnods  , sizeof(PetscInt),    1                 , fp);
	fread(&ds.tcels  , sizeof(PetscInt),    1                 , fp);
	fread(&ds.nnods  , sizeof(PetscInt),    1                 , fp);
	fread(&ds.ncels  , sizeof(PetscInt),    1                 , fp);
	fread(&ds.grprev , sizeof(PetscInt),    1                 , fp);
	fread(&ds.grnext , sizeof(PetscInt),    1                 , fp);
	fread(&ds.color  , sizeof(PetscInt),    1                 , fp);
	fread(&ds.h_uni  , sizeof(PetscScalar), 1                 , fp);
	fread(&ds.h_min  , sizeof(PetscScalar), 1                 , fp);
	fread(&ds.h_max  , sizeof(PetscScalar), 1                 , fp);

	//----------------------------------
	// FDSTAG STRETCH 1D
	//----------------------------------
	fread(&ms.nsegs  , sizeof(PetscInt),    1                 , fp);
	fread(ms.istart  , sizeof(PetscInt),    (size_t)ms.nsegs+1, fp);
	fread(ms.xstart  , sizeof(PetscScalar), (size_t)ms.nsegs+1, fp);
}
//---------------------------------------------------------------------------

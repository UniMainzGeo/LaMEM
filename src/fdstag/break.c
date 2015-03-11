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
#include "advect.h"
#include "marker.h"
#include "break.h"

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BreakCheck"
PetscErrorCode BreakCheck(UserCtx *user)
{
	// check if breakpoints exist and restart new simulation is not - very useful for chain jobs

	FILE        *fp;
	char        *fname;

	PetscFunctionBegin;

	// check info file name on rank 0
	if(ISRankZero(PETSC_COMM_WORLD))
	{
		asprintf(&fname, "./Breakpoint/Breakpoint_info.0.out");

		// return and restart simulation if breakpoints are not available
		fp = fopen(fname, "r" );

		if(!fp)
		{
			// new simulation
			user->restart = 0;
			PetscPrintf(PETSC_COMM_WORLD," No breakpoints detected -> starting new simulation \n");

			// free and return
			free(fname);

			PetscFunctionReturn(0);
		}

		// free and close files
		free(fname);
		fclose(fp);
	}
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BreakWrite"
PetscErrorCode BreakWrite(UserCtx *user, AdvCtx *actx, JacType jtype)
{
	// staggered grid
	FDSTAG         *fs;
	JacRes         *jr;
	FILE           *fp;
	char           *fname;
	PetscInt        n;
	PetscScalar    *gsol;
	PetscInt       initGuessFlag, jtypeFlag;
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
	if (jr->bc->bgAct)
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
	//   GSOL - Solution vector
	//============================================================
	// set local size of vector
	n  = actx->fs->dof.ln;

	// compile file name
	asprintf(&fname, "./Breakpoint/Breakpoint_gsol.%lld.out",(LLD)actx->iproc);

	// open file for binary output
	fp = fopen(fname, "w" );

	// get vector array
	ierr = VecGetArray(jr->gsol,&gsol); CHKERRQ(ierr);

	// write to file
	fwrite( gsol, sizeof(PetscScalar),(size_t)n, fp);

	// restore vector array
	ierr = VecRestoreArray(jr->gsol,&gsol); CHKERRQ(ierr);

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
	FILE        *fp;
	char        *fname;
	PetscScalar *gsol;
	PetscInt     n;
	JacType     j;
	PetscInt    initGuessFlag, jtypeFlag;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// set context
	jr = actx->jr;

	//============================================================
	//   GSOL - Solution vector
	//============================================================
	// set local size of vector
	n  = actx->fs->dof.ln;

	// compile file name
	asprintf(&fname, "./Breakpoint/Breakpoint_gsol.%lld.out",(LLD)actx->iproc);

	// open file for reading
	fp = fopen(fname, "r" );

	// get vector array
	ierr = VecGetArray(jr->gsol,&gsol); CHKERRQ(ierr);

	// read array
	fread(  gsol, sizeof(PetscScalar),(size_t)n, fp);

	// restore vector array
	ierr = VecRestoreArray(jr->gsol,&gsol); CHKERRQ(ierr);

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

	// read markers
	fread( actx->markers, sizeof(Marker)*(size_t)actx->markcap, 1, fp);

	// close and free memory
	free(fname);
	fclose(fp);

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

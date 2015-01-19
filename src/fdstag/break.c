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
#include "advect.h"
#include "marker.h"
#include "break.h"

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BreakWriteMain"
PetscErrorCode BreakWriteMain(UserCtx *user, AdvCtx *actx)
{
	// staggered grid
	FDSTAG *fs;
	JacRes *jr;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	PetscPrintf(PETSC_COMM_WORLD," Writing Breakpoint files... \n");

	// initialize context
	fs = actx->fs;
	jr = actx->jr;

	// create directory
	ierr = LaMEMCreateOutputDirectory("./Breakpoint"); CHKERRQ(ierr);

	// write markers
	ierr = BreakWriteMark(actx); CHKERRQ(ierr);

	// write FDSTAG context only for background strainrate
	if (jr->bc->bgAct) { ierr = BreakWriteGrid(user, fs, actx); CHKERRQ(ierr); }

	// write solution and other info
	ierr = BreakWriteSol(jr);              CHKERRQ(ierr);
	ierr = BreakWriteInfo(user, actx, jr); CHKERRQ(ierr);

	// all other ranks should wait
	ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD," Finished writing breakpoint files \n");
	PetscPrintf(PETSC_COMM_WORLD," \n");

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BreakReadMain"
PetscErrorCode BreakReadMain(UserCtx *user, AdvCtx *actx, JacRes *jr)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// read data from file
	ierr = BreakReadSol(jr);              CHKERRQ(ierr);
	ierr = BreakReadInfo(user, actx, jr); CHKERRQ(ierr);

	// all other ranks should wait
	ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BreakWriteMark"
PetscErrorCode BreakWriteMark(AdvCtx *actx)
{
	int          fid;
	PetscInt     imark, nprop;
	Marker      *P;
	char        *fname;
	PetscScalar *markbuf, *markptr, s_nummark;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	//----------------------------------
	// 0 - MARKERS History
	//----------------------------------
	nprop = 13;

	// create write buffer
	ierr = PetscMalloc((size_t)(nprop*actx->nummark)*sizeof(PetscScalar), &markbuf); CHKERRQ(ierr);

	// copy data from storage into buffer
	for(imark = 0, markptr = markbuf; imark < actx->nummark; imark++, markptr += nprop)
	{
		P           =              &actx->markers[imark];
		markptr[0]  =              P->X[0];
		markptr[1]  =              P->X[1];
		markptr[2]  =              P->X[2];
		markptr[3]  = (PetscScalar)P->phase;
		markptr[4]  =              P->T;
		markptr[5]  =              P->p;
		markptr[6]  =              P->APS;
		markptr[7]  =              P->S.xx;
		markptr[8]  =              P->S.xy;
		markptr[9]  =              P->S.xz;
		markptr[10] =              P->S.yy;
		markptr[11] =              P->S.yz;
		markptr[12] =              P->S.zz;
	}

	// compile file name
	asprintf(&fname, "./Breakpoint/Breakpoint_mark.%lld.out",(LLD)actx->iproc);

	// open file for binary output
	ierr = PetscBinaryOpen(fname, FILE_MODE_WRITE, &fid); CHKERRQ(ierr);

	// write binary output
	s_nummark = (PetscScalar)actx->nummark;
	ierr = PetscBinaryWrite(fid, &s_nummark, 1,                PETSC_SCALAR, PETSC_FALSE); CHKERRQ(ierr);
	ierr = PetscBinaryWrite(fid, markbuf, nprop*actx->nummark, PETSC_SCALAR, PETSC_FALSE); CHKERRQ(ierr);

	// close fid and free memory
	ierr = PetscBinaryClose(fid); CHKERRQ(ierr);
	free(fname);

	// destroy buffer
	ierr = PetscFree(markbuf); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BreakReadMark"
PetscErrorCode BreakReadMark(AdvCtx *actx)
{
	int          fid;
	PetscInt     imark, nprop, nummark;
	Marker      *P;
	char        *fname;
	PetscScalar *markbuf, *markptr, s_nummark;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	//----------------------------------
	// 0 - MARKERS History
	//----------------------------------
	// number of prop in Marker structure
	nprop = 13;

	// compile file name and open file
	asprintf(&fname, "./Breakpoint/Breakpoint_mark.%lld.out",(LLD)actx->iproc);
	ierr = PetscBinaryOpen(fname, FILE_MODE_READ, &fid); CHKERRQ(ierr);

	// read number of local of markers
	ierr = PetscBinaryRead(fid, &s_nummark, 1, PETSC_SCALAR); CHKERRQ(ierr);
	nummark = (PetscInt)s_nummark;

	// set number of markers
	actx->nummark = nummark;

	// allocate marker buffer
	ierr = PetscMalloc((size_t)(nprop*actx->nummark)*sizeof(PetscScalar), &markbuf); CHKERRQ(ierr);

	// read markers into buffer
	ierr = PetscBinaryRead(fid, markbuf, nprop*actx->nummark, PETSC_SCALAR); CHKERRQ(ierr);

	// close fid and free memory
	ierr = PetscBinaryClose(fid); CHKERRQ(ierr);
	free(fname);

	// copy buffer to marker storage
	for(imark = 0, markptr = markbuf; imark < actx->nummark; imark++, markptr += nprop)
	{
		P        = &actx->markers[imark];
		P->X[0]  = markptr[0]           ;
		P->X[1]  = markptr[1]           ;
		P->X[2]  = markptr[2]           ;
		P->phase = (PetscInt)markptr[3] ;
		P->T     = markptr[4]           ;
		P->p     = markptr[5]           ;
		P->APS   = markptr[6]           ;
		P->S.xx  = markptr[7]           ;
		P->S.xy  = markptr[8]           ;
		P->S.xz  = markptr[9]           ;
		P->S.yy  = markptr[10]          ;
		P->S.yz  = markptr[11]          ;
		P->S.zz  = markptr[12]          ;
	}

	// destroy buffer
	ierr = PetscFree(markbuf); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BreakWriteSol"
PetscErrorCode BreakWriteSol(JacRes *jr)
{
	char        *fname;
	PetscViewer viewer;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	//----------------------------------
	// 1 - SOLUTION Vector
	//----------------------------------
	// compile file name - solution
	asprintf(&fname, "./Breakpoint/Breakpoint_gsol.out");

	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,fname,FILE_MODE_WRITE,&viewer); CHKERRQ(ierr);
	ierr = VecView(jr->gsol,viewer); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

	free(fname);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BreakReadSol"
PetscErrorCode BreakReadSol(JacRes *jr)
{
	char        *fname;
	PetscViewer viewer;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	//----------------------------------
	// 1 - SOLUTION Vector
	//----------------------------------
	// compile file name - solution
	asprintf(&fname, "./Breakpoint/Breakpoint_gsol.out");

	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,fname,FILE_MODE_READ,&viewer); CHKERRQ(ierr);
	ierr = VecLoad(jr->gsol,viewer); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

	free(fname);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BreakWriteGrid"
PetscErrorCode BreakWriteGrid(UserCtx *user, FDSTAG *fs, AdvCtx *actx)
{
	int         fid;
	char        *fname;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	//----------------------------------
	// 2 - GRID information
	//----------------------------------
	// compile file name
	asprintf(&fname, "./Breakpoint/Breakpoint_grid.%lld.out",(LLD)actx->iproc);
	ierr = PetscBinaryOpen(fname, FILE_MODE_WRITE, &fid); CHKERRQ(ierr);

	//----------------------------------
	// Domain information
	//----------------------------------
	PetscBinaryWrite(fid, &user->W,        1,               PETSC_SCALAR, PETSC_FALSE);
	PetscBinaryWrite(fid, &user->L,        1,               PETSC_SCALAR, PETSC_FALSE);
	PetscBinaryWrite(fid, &user->H,        1,               PETSC_SCALAR, PETSC_FALSE);
	PetscBinaryWrite(fid, &user->x_left ,  1,               PETSC_SCALAR, PETSC_FALSE);
	PetscBinaryWrite(fid, &user->y_front,  1,               PETSC_SCALAR, PETSC_FALSE);
	PetscBinaryWrite(fid, &user->z_bot  ,  1,               PETSC_SCALAR, PETSC_FALSE);

	//----------------------------------
	// FDSTAG
	//----------------------------------
	ierr = BreakWriteDiscret1D(fid, fs->dsx); CHKERRQ(ierr);
	ierr = BreakWriteDiscret1D(fid, fs->dsy); CHKERRQ(ierr);
	ierr = BreakWriteDiscret1D(fid, fs->dsz); CHKERRQ(ierr);

	// close fid and free memory
	ierr = PetscBinaryClose(fid); CHKERRQ(ierr);
	free(fname);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BreakReadGrid"
PetscErrorCode BreakReadGrid(UserCtx *user, FDSTAG *fs)
{
	int         fid;
	char        *fname;
	PetscMPIInt iproc;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	//----------------------------------
	// 2 - GRID information
	//----------------------------------
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &iproc); CHKERRQ(ierr);

	// compile file name - fdstag context
	asprintf(&fname, "./Breakpoint/Breakpoint_grid.%lld.out",(LLD)iproc);
	ierr = PetscBinaryOpen(fname, FILE_MODE_READ, &fid); CHKERRQ(ierr);

	//----------------------------------
	// Domain information
	//----------------------------------
	PetscBinaryRead(fid, &user->W,        1,               PETSC_SCALAR);
	PetscBinaryRead(fid, &user->L,        1,               PETSC_SCALAR);
	PetscBinaryRead(fid, &user->H,        1,               PETSC_SCALAR);
	PetscBinaryRead(fid, &user->x_left ,  1,               PETSC_SCALAR);
	PetscBinaryRead(fid, &user->y_front,  1,               PETSC_SCALAR);
	PetscBinaryRead(fid, &user->z_bot  ,  1,               PETSC_SCALAR);

	//----------------------------------
	// FDSTAG
	//----------------------------------
	ierr = BreakReadDiscret1D(fid, fs->dsx); CHKERRQ(ierr);
	ierr = BreakReadDiscret1D(fid, fs->dsy); CHKERRQ(ierr);
	ierr = BreakReadDiscret1D(fid, fs->dsz); CHKERRQ(ierr);

	// close fid and free memory
	ierr = PetscBinaryClose(fid); CHKERRQ(ierr);
	free(fname);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BreakWriteDiscret1D"
PetscErrorCode BreakWriteDiscret1D(int fid, Discret1D ds)
{
	PetscFunctionBegin;

	PetscBinaryWrite(fid, &ds.nproc,  1,          PETSC_INT,    PETSC_FALSE);
	PetscBinaryWrite(fid, &ds.rank,   1,          PETSC_INT,    PETSC_FALSE);
	PetscBinaryWrite(fid, ds.starts,  ds.nproc+1, PETSC_INT,    PETSC_FALSE);
	PetscBinaryWrite(fid, &ds.pstart, 1,          PETSC_INT,    PETSC_FALSE);
	PetscBinaryWrite(fid, ds.ncoor,   ds.nnods+1, PETSC_SCALAR, PETSC_FALSE);
	PetscBinaryWrite(fid, ds.ccoor,   ds.ncels+1, PETSC_SCALAR, PETSC_FALSE);
/*	PetscBinaryWrite(fid, &ds.tnods,  1,          PETSC_INT,    PETSC_FALSE);
	PetscBinaryWrite(fid, &ds.tcels,  1,          PETSC_INT,    PETSC_FALSE);
	PetscBinaryWrite(fid, &ds.nnods,  1,          PETSC_INT,    PETSC_FALSE);
	PetscBinaryWrite(fid, &ds.ncels,  1,          PETSC_INT,    PETSC_FALSE);
	PetscBinaryWrite(fid, ds.nbuff,   ds.bufsz,   PETSC_SCALAR, PETSC_FALSE);
	PetscBinaryWrite(fid, ds.cbuff,   ds.ncels+2, PETSC_SCALAR, PETSC_FALSE);
	PetscBinaryWrite(fid, &ds.bufsz,  1,          PETSC_INT,    PETSC_FALSE);
	PetscBinaryWrite(fid, &ds.grprev, 1,          PETSC_INT,    PETSC_FALSE);
	PetscBinaryWrite(fid, &ds.grnext, 1,          PETSC_INT,    PETSC_FALSE);
	PetscBinaryWrite(fid, &ds.color,  1,          PETSC_INT,    PETSC_FALSE);
	PetscBinaryWrite(fid, &ds.h_uni,  1,          PETSC_SCALAR, PETSC_FALSE);
	PetscBinaryWrite(fid, &ds.h_min,  1,          PETSC_SCALAR, PETSC_FALSE);
	PetscBinaryWrite(fid, &ds.h_max,  1,          PETSC_SCALAR, PETSC_FALSE);*/

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BreakReadDiscret1D"
PetscErrorCode BreakReadDiscret1D(int fid, Discret1D ds)
{
	PetscFunctionBegin;

	PetscBinaryRead(fid, &ds.nproc,  1,          PETSC_INT);
	PetscBinaryRead(fid, &ds.rank,   1,          PETSC_INT);
	PetscBinaryRead(fid, ds.starts,  ds.nproc+1, PETSC_INT);
	PetscBinaryRead(fid, &ds.pstart, 1,          PETSC_INT);
	PetscBinaryRead(fid, ds.ncoor,   ds.nnods+1, PETSC_SCALAR);
	PetscBinaryRead(fid, ds.ccoor,   ds.ncels+1, PETSC_SCALAR);
/*	PetscBinaryRead(fid, &ds.tnods,  1,          PETSC_INT);
	PetscBinaryRead(fid, &ds.tcels,  1,          PETSC_INT);
	PetscBinaryRead(fid, &ds.nnods,  1,          PETSC_INT);
	PetscBinaryRead(fid, &ds.ncels,  1,          PETSC_INT);
	PetscBinaryRead(fid, ds.nbuff,   ds.bufsz,   PETSC_SCALAR);
	PetscBinaryRead(fid, ds.cbuff,   ds.ncels+2, PETSC_SCALAR);
	PetscBinaryRead(fid, &ds.bufsz,  1,          PETSC_INT);
	PetscBinaryRead(fid, &ds.grprev, 1,          PETSC_INT);
	PetscBinaryRead(fid, &ds.grnext, 1,          PETSC_INT);
	PetscBinaryRead(fid, &ds.color,  1,          PETSC_INT);
	PetscBinaryRead(fid, &ds.h_uni,  1,          PETSC_SCALAR);
	PetscBinaryRead(fid, &ds.h_min,  1,          PETSC_SCALAR);
	PetscBinaryRead(fid, &ds.h_max,  1,          PETSC_SCALAR);*/

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BreakWriteInfo"
PetscErrorCode BreakWriteInfo(UserCtx *user, AdvCtx *actx, JacRes *jr)
{
	int         fid;
	char        *fname;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	//----------------------------------
	// 3 - INFORMATION
	//----------------------------------
	// compile file name and open file for binary output
	asprintf(&fname, "./Breakpoint/Breakpoint_info.%lld.out",(LLD)actx->iproc);
	ierr = PetscBinaryOpen(fname, FILE_MODE_WRITE, &fid); CHKERRQ(ierr);

	//----------------------------------
	// Advection/Time Stepping
	//----------------------------------

	PetscBinaryWrite(fid, &jr->ts.istep,   1,               PETSC_INT,    PETSC_FALSE);
	PetscBinaryWrite(fid, &jr->ts.dt,      1,               PETSC_SCALAR, PETSC_FALSE);
	PetscBinaryWrite(fid, &jr->ts.time,    1,               PETSC_SCALAR, PETSC_FALSE);

	//----------------------------------
	// Boundary conditions related info
	//----------------------------------
	// pushing block center coordinates
	if (user->AddPushing)
	{
		PetscBinaryWrite(fid, &user->Pushing.x_center_block, 1, PETSC_SCALAR, PETSC_FALSE);
		PetscBinaryWrite(fid, &user->Pushing.y_center_block, 1, PETSC_SCALAR, PETSC_FALSE);
		PetscBinaryWrite(fid, &user->Pushing.z_center_block, 1, PETSC_SCALAR, PETSC_FALSE);
	}

	//----------------------------------
	// OTHER info
	//----------------------------------
	// store breakpoint number
	PetscBinaryWrite(fid, &user->break_point_number, 1,     PETSC_INT,    PETSC_FALSE);

	// close fid and free memory
	ierr = PetscBinaryClose(fid); CHKERRQ(ierr);
	free(fname);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BreakReadInfo"
PetscErrorCode BreakReadInfo(UserCtx *user, AdvCtx *actx, JacRes *jr)
{
	int         fid;
	char        *fname;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	//----------------------------------
	// 3 - INFORMATION
	//----------------------------------
	// compile file name and open file
	asprintf(&fname, "./Breakpoint/Breakpoint_info.%lld.out",(LLD)actx->iproc);
	ierr = PetscBinaryOpen(fname, FILE_MODE_READ, &fid); CHKERRQ(ierr);

	//----------------------------------
	// Advection/Time Stepping
	//----------------------------------
	PetscBinaryRead(fid, &jr->ts.istep,   1,               PETSC_INT   );
	PetscBinaryRead(fid, &jr->ts.dt,      1,               PETSC_SCALAR);
	PetscBinaryRead(fid, &jr->ts.time,    1,               PETSC_SCALAR);

	//----------------------------------
	// Boundary conditions related info
	//----------------------------------
	// pushing block center coordinates
	if (user->AddPushing)
	{
		PetscBinaryRead(fid, &user->Pushing.x_center_block, 1, PETSC_SCALAR);
		PetscBinaryRead(fid, &user->Pushing.y_center_block, 1, PETSC_SCALAR);
		PetscBinaryRead(fid, &user->Pushing.z_center_block, 1, PETSC_SCALAR);
	}

	//----------------------------------
	// OTHER info
	//----------------------------------
	// store breakpoint number
	PetscBinaryRead(fid, &user->break_point_number, 1,     PETSC_INT   );

	// close fid and free memory
	ierr = PetscBinaryClose(fid); CHKERRQ(ierr);
	free(fname);

	PetscPrintf(PETSC_COMM_WORLD," \n");
	PetscPrintf(PETSC_COMM_WORLD,"RESTART from Time step = %lld \n",(LLD)jr->ts.istep);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

/*
 * BFBT.cpp
 *
 *  Created on: 04.03.2020
 *      Author: daniel
 */

//---------------------------------------------------------------------------
//..........................   BFBT FUNCTIONS   .............................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "JacRes.h"
#include "phase.h"
#include "scaling.h"
#include "fdstag.h"
#include "tssolve.h"
#include "bc.h"
#include "matrix.h"
#include "multigrid.h"
#include "lsolve.h"
#include "BFBT.h"
#include "marker.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PMatBFBTCreate"
PetscErrorCode PMatBFBTCreate(PMat pm)
{
	PMatBlock      *P;
	JacRes         *jr;
	FDSTAG         *fs;
	DOFIndex       *dof;
	PetscInt        lnv, stv;
	const PetscInt *lx, *ly, *lz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// BFBT cases only
	if(pm->stype != _wBFBT_) PetscFunctionReturn(0);

	P   = (PMatBlock*)pm->data;
	jr  = pm->jr;
	fs  = jr->fs;
	dof = &fs->dof;
	lnv = dof->lnv;
	stv = dof->stv;

	// get cell center grid partitioning
	ierr = DMDAGetOwnershipRanges(fs->DA_CEN, &lx, &ly, &lz); CHKERRQ(ierr);

	// create DMDA
	ierr = DMDACreate3dSetUp(PETSC_COMM_WORLD,
		DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
		DMDA_STENCIL_STAR,
		fs->dsx.tcels, fs->dsy.tcels, fs->dsz.tcels,
		fs->dsx.nproc, fs->dsy.nproc, fs->dsz.nproc,
		1, 1, lx, ly, lz, &P->DA_P); CHKERRQ(ierr);

	// set proper interpolation type for multigrid
	ierr = DMDASetInterpolationType(P->DA_P, DMDA_Q0); CHKERRQ(ierr);

	// create matrix
	ierr = DMCreateMatrix(P->DA_P, &P->K); CHKERRQ(ierr);

	// set matrix options (development)
	ierr = MatSetOption(P->K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE); CHKERRQ(ierr);
	ierr = MatSetOption(P->K, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);   CHKERRQ(ierr);
	ierr = MatSetOption(P->K, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);       CHKERRQ(ierr);
	ierr = MatSetOption(P->K, MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_TRUE);      CHKERRQ(ierr);

	// create scaling matrix
	ierr = MatAIJCreateDiag(lnv, stv, &P->C); CHKERRQ(ierr);

	// allocate work vectors
	ierr = VecDuplicate(P->xv, &P->wv2);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PMatBFBTAssemble"
PetscErrorCode PMatBFBTAssemble(PMat pm)
{
	PMatBlock  *P;
	JacRes     *jr;
	FDSTAG     *fs;
	BCCtx      *bc;
	Vec         lvEtaCen;
	MatStencil  row[1], col[7];
	PetscInt    I1, I2, J1, J2, K1, K2;
	PetscInt    Ip1, Im1, Jp1, Jm1, Kp1, Km1;
	PetscInt    mnx, mny, mnz, mcx, mcy, mcz;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter;
	PetscScalar v[7], d[6], g[6], cfe[6], cfp[6], cfv[6];
	PetscScalar ***lEta, ***bcvx,  ***bcvy,  ***bcvz, ***bcp;
	PetscScalar bdx, fdx, bdy, fdy, bdz, fdz, dx, dy, dz;
	PetscBool	flg;
	char        pname[_str_len_];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// BFBT cases only
	if(pm->stype != _wBFBT_) PetscFunctionReturn(0);

	// access context variables
	P 	 = (PMatBlock*)pm->data;
	jr   = pm->jr;
	fs   = jr->fs;
	bc   = jr->bc;

	// initialize index bounds
	mcx = fs->dsx.tcels - 1;
	mcy = fs->dsy.tcels - 1;
	mcz = fs->dsz.tcels - 1;

	mnx  = fs->dsx.tnods - 1;
	mny  = fs->dsy.tnods - 1;
	mnz  = fs->dsz.tnods - 1;

	//===============
	// CELL VISCOSITY
	//===============

	ierr = DMGetLocalVector(fs->DA_CEN, &lvEtaCen); CHKERRQ(ierr);

	ierr = VecZeroEntries(lvEtaCen); CHKERRQ(ierr);

	ierr = DMDAVecGetArray(fs->DA_CEN, lvEtaCen, &lEta); CHKERRQ(ierr);

	iter = 0;
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD,"-------------------------------------------------------------------------- \n");
	ierr = PetscOptionsGetString(NULL, NULL, "-pcmat_BFBT_viscositySmoothing", pname, _str_len_, &flg); CHKERRQ(ierr);
	if(flg==PETSC_TRUE){
		START_STD_LOOP
		{
			lEta[k][j][i] = jr->svCell[iter++].eta_smoothed;
			if(iter%5500==0) PetscPrintf(PETSC_COMM_WORLD, "%lld \n", (LLD)lEta[k][j][i]);
		}
		END_STD_LOOP
	}else{
		START_STD_LOOP
		{
			lEta[k][j][i] = jr->svCell[iter++].svDev.eta;
			if(iter%5500==0) PetscPrintf(PETSC_COMM_WORLD, "%lld \n", (LLD)lEta[k][j][i]);
		}
		END_STD_LOOP
	}
	PetscPrintf(PETSC_COMM_WORLD,"-------------------------------------------------------------------------- \n");

	ierr = DMDAVecRestoreArray(fs->DA_CEN, lvEtaCen, &lEta); CHKERRQ(ierr);

	LOCAL_TO_LOCAL(fs->DA_CEN, lvEtaCen)

	//===============
	// SCALING MATRIX
	//===============

	ierr = DMDAVecGetArray(fs->DA_CEN, lvEtaCen, &lEta); CHKERRQ(ierr);

	// set iterator
	iter = fs->dof.stv;

	// X-points
	ierr = DMDAGetCorners(fs->DA_X, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		I1 = i;   if(I1 == mnx) I1--;
		I2 = i-1; if(I2 == -1)  I2++;

		ierr = MatSetValue(P->C, iter, iter, 1.0/sqrt((lEta[k][j][I1] + lEta[k][j][I2])/2.0), INSERT_VALUES); CHKERRQ(ierr);

		iter++;
	}
	END_STD_LOOP

	// Y-points
	ierr = DMDAGetCorners(fs->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		J1 = j;   if(J1 == mny) J1--;
		J2 = j-1; if(J2 == -1)  J2++;

		ierr = MatSetValue(P->C, iter, iter, 1.0/sqrt((lEta[k][J1][i] + lEta[k][J2][i])/2.0), INSERT_VALUES); CHKERRQ(ierr);

		iter++;
	}
	END_STD_LOOP

	// Z-points
	ierr = DMDAGetCorners(fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		K1 = k;   if(K1 == mnz) K1--;
		K2 = k-1; if(K2 == -1)  K2++;

		ierr = MatSetValue(P->C, iter, iter, 1.0/sqrt((lEta[K1][j][i] + lEta[K2][j][i])/2.0), INSERT_VALUES); CHKERRQ(ierr);

		iter++;
	}
	END_STD_LOOP

	// assemble scaling matrix
	ierr = MatAIJAssemble(P->C, bc->vNumSPC, bc->vSPCList, 1.0); CHKERRQ(ierr);

	//=======================
	// PRECONDITIONING MATRIX
	//=======================

	// clear matrix
	ierr = MatZeroEntries(P->K); CHKERRQ(ierr);

	// access boundary constraints
	ierr = DMDAVecGetArray(fs->DA_X,   bc->bcvx, &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   bc->bcvy, &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   bc->bcvz, &bcvz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, bc->bcp,  &bcp);  CHKERRQ(ierr);

	iter = 0;
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// check index bounds, get pressure TPC multipliers
		Im1 = i-1; cfp[0] = 1.0; if(Im1 < 0)   { Im1++; if(bcp[k][j][i-1] != DBL_MAX) cfp[0] = -1.0; }
		Ip1 = i+1; cfp[1] = 1.0; if(Ip1 > mcx) { Ip1--; if(bcp[k][j][i+1] != DBL_MAX) cfp[1] = -1.0; }
		Jm1 = j-1; cfp[2] = 1.0; if(Jm1 < 0)   { Jm1++; if(bcp[k][j-1][i] != DBL_MAX) cfp[2] = -1.0; }
		Jp1 = j+1; cfp[3] = 1.0; if(Jp1 > mcy) { Jp1--; if(bcp[k][j+1][i] != DBL_MAX) cfp[3] = -1.0; }
		Km1 = k-1; cfp[4] = 1.0; if(Km1 < 0)   { Km1++; if(bcp[k-1][j][i] != DBL_MAX) cfp[4] = -1.0; }
		Kp1 = k+1; cfp[5] = 1.0; if(Kp1 > mcz) { Kp1--; if(bcp[k+1][j][i] != DBL_MAX) cfp[5] = -1.0; }

		// get velocity SPC multipliers
		cfv[0] = 1.0; if(bcvx[k][j][i]   != DBL_MAX) cfv[0] = 0.0;
		cfv[1] = 1.0; if(bcvx[k][j][i+1] != DBL_MAX) cfv[1] = 0.0;
		cfv[2] = 1.0; if(bcvy[k][j][i]   != DBL_MAX) cfv[2] = 0.0;
		cfv[3] = 1.0; if(bcvy[k][j+1][i] != DBL_MAX) cfv[3] = 0.0;
		cfv[4] = 1.0; if(bcvz[k][j][i]   != DBL_MAX) cfv[4] = 0.0;
		cfv[5] = 1.0; if(bcvz[k+1][j][i] != DBL_MAX) cfv[5] = 0.0;

		// get viscosity scaling coefficients
		cfe[0] = sqrt((lEta[k][j][i] + lEta[k][j][Im1])/2.0);
		cfe[1] = sqrt((lEta[k][j][i] + lEta[k][j][Ip1])/2.0);
		cfe[2] = sqrt((lEta[k][j][i] + lEta[k][Jm1][i])/2.0);
		cfe[3] = sqrt((lEta[k][j][i] + lEta[k][Jp1][i])/2.0);
		cfe[4] = sqrt((lEta[k][j][i] + lEta[Km1][j][i])/2.0);
		cfe[5] = sqrt((lEta[k][j][i] + lEta[Kp1][j][i])/2.0);

		// get mesh steps
		dx = SIZE_CELL(i, sx, fs->dsx);
		dy = SIZE_CELL(j, sy, fs->dsy);
		dz = SIZE_CELL(k, sz, fs->dsz);

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_NODE(i, sx, fs->dsx);   fdx = SIZE_NODE(i+1, sx, fs->dsx);
		bdy = SIZE_NODE(j, sy, fs->dsy);   fdy = SIZE_NODE(j+1, sy, fs->dsy);
		bdz = SIZE_NODE(k, sz, fs->dsz);   fdz = SIZE_NODE(k+1, sz, fs->dsz);

		/*
		if(pm->BFBTtype==_unscaled_){

			g[0] =  1.0/bdx;   d[0] =  cfv[0]/dx;
			g[1] = -1.0/fdx;   d[1] = -cfv[1]/dx;
			g[2] =  1.0/bdy;   d[2] =  cfv[2]/dy;
			g[3] = -1.0/fdy;   d[3] = -cfv[3]/dy;
			g[4] =  1.0/bdz;   d[4] =  cfv[4]/dz;
			g[5] = -1.0/fdz;   d[5] = -cfv[5]/dz;

		}else{
		*/

		// get gradient vector & divergence vector scaled by viscosity and velocity SPC multipliers
		g[0] =  1.0/bdx;   d[0] =  cfv[0]/dx/cfe[0];
		g[1] = -1.0/fdx;   d[1] = -cfv[1]/dx/cfe[1];
		g[2] =  1.0/bdy;   d[2] =  cfv[2]/dy/cfe[2];
		g[3] = -1.0/fdy;   d[3] = -cfv[3]/dy/cfe[3];
		g[4] =  1.0/bdz;   d[4] =  cfv[4]/dz/cfe[4];
		g[5] = -1.0/fdz;   d[5] = -cfv[5]/dz/cfe[5];

		//}

		// compute matrix row, apply pressure TPC multipliers
		v[0] = -g[0]*d[0]*cfp[0];
		v[1] = -g[1]*d[1]*cfp[1];
		v[2] = -g[2]*d[2]*cfp[2];
		v[3] = -g[3]*d[3]*cfp[3];
		v[4] = -g[4]*d[4]*cfp[4];
		v[5] = -g[5]*d[5]*cfp[5];
		v[6] =  g[0]*d[0] + g[1]*d[1] + g[2]*d[2] + g[3]*d[3] + g[4]*d[4] + g[5]*d[5];

		// set diagonal for fully constrained cells
		if(!v[6]) v[6] = 1.0;

		// set row/column indices
		row[0].k = k;   row[0].j = j;   row[0].i = i;   row[0].c = 0;
		col[0].k = k;   col[0].j = j;   col[0].i = Im1; col[0].c = 0;
		col[1].k = k;   col[1].j = j;   col[1].i = Ip1; col[1].c = 0;
		col[2].k = k;   col[2].j = Jm1; col[2].i = i;   col[2].c = 0;
		col[3].k = k;   col[3].j = Jp1; col[3].i = i;   col[3].c = 0;
		col[4].k = Km1; col[4].j = j;   col[4].i = i;   col[4].c = 0;
		col[5].k = Kp1; col[5].j = j;   col[5].i = i;   col[5].c = 0;
		col[6].k = k;   col[6].j = j;   col[6].i = i;   col[6].c = 0;

		// set matrix coefficients
		ierr = MatSetValuesStencil(P->K, 1, row, 7, col, v, ADD_VALUES); CHKERRQ(ierr);

	}
	END_STD_LOOP

	// assemble preconditioning matrix
	ierr = MatAIJAssemble(P->K,  bc->pNumSPC, bc->pSPCList, 1.0); CHKERRQ(ierr);

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_CEN, lvEtaCen, &lEta); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_X,   bc->bcvx, &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   bc->bcvy, &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   bc->bcvz, &bcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, bc->bcp,  &bcp);  CHKERRQ(ierr);

	// restore viscosity vector
	ierr = DMRestoreLocalVector (fs->DA_CEN, &lvEtaCen); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PMatBFBTDestroy"
PetscErrorCode PMatBFBTDestroy(PMat pm)
{
	PMatBlock *P;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// BFBT cases only
	if(pm->stype != _wBFBT_) PetscFunctionReturn(0);

	P = (PMatBlock*)pm->data;

	ierr = DMDestroy (&P->DA_P); CHKERRQ(ierr);
	ierr = MatDestroy(&P->K); 	 CHKERRQ(ierr);
	ierr = MatDestroy(&P->C); 	 CHKERRQ(ierr);
	ierr = VecDestroy(&P->wv2);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesBFBTApply"
PetscErrorCode PCStokesBFBTApply(Mat JP, Vec x, Vec y)
{
	//=============
	// BLOCK w-BFBT
	//=============

	PCStokes    pc;
	PCStokesBF *bf;
	PMatBlock  *P;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	ierr = MatShellGetContext(JP, (void**)&pc); CHKERRQ(ierr);

	bf = (PCStokesBF*)pc->data;
	P  = (PMatBlock*) pc->pm->data;

	// y   = -(S^⁻1)*x
	// S⁻1 =  (K^⁻1)*B*C*A*C*B^T*(K^⁻1)
	// K   =  B*C*B^T


//	if(pm->BFBTtype==_unscaled_){
/*
		ierr = KSPSolve(bf->pksp, x, P->wp);   CHKERRQ(ierr); // wp = (K^⁻1)*x
		ierr = MatMult(P->Avp, P->wp, P->wv);  CHKERRQ(ierr); // wv = Avp*wp
		ierr = MatMult(P->Avv, P->wv, P->wv2); CHKERRQ(ierr); // wv2 = Avv * wv
		ierr = MatMult(P->Apv, P->wv2, P->wp); CHKERRQ(ierr); // wp = Apv*wv2
		ierr = KSPSolve(bf->pksp, P->wp, y);   CHKERRQ(ierr); // y = (K^⁻1)*wp
*/
	//}else if(pm->BFBTtype==scaled){

//	}else{

		ierr = KSPSolve(bf->pksp, x, P->wp);   CHKERRQ(ierr); // wp = (K^⁻1)*x

		ierr = MatMult(P->Avp, P->wp, P->wv);  CHKERRQ(ierr); // wv = Avp*wp

		ierr = MatMult(P->C, P->wv, P->wv2);   CHKERRQ(ierr); // wv2 = C*wv

		ierr = MatMult(P->Avv, P->wv2, P->wv); CHKERRQ(ierr); // wv = Avv * wv2

		ierr = MatMult(P->C, P->wv, P->wv2);   CHKERRQ(ierr); // wv2 = C*wv

		ierr = MatMult(P->Apv, P->wv2, P->wp); CHKERRQ(ierr); // wp = Apv*wv2

		ierr = KSPSolve(bf->pksp, P->wp, y);   CHKERRQ(ierr); // y = (K^⁻1)*wp

//	}



	ierr = VecScale(y, -1.0);              CHKERRQ(ierr); // y = -y

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "BFBTGaussianSmoothing"
PetscErrorCode BFBTGaussianSmoothing(JacRes *jr, PetscInt nSpheres, GeomPrim *geom)
{
	FDSTAG     *fs;
	PetscScalar	DynamicRatio, delta;
	PetscScalar etamax, etamin, radius;
	PetscScalar IndicatorValue, center[3], coords[3];
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter;
	PetscInt	iSpheres, ngeom;
	PetscScalar maximum, expfunc, vecabs;
	Discret1D  *dsx, *dsy, *dsz;
	GeomPrim   *sphere;
	PetscBool	flg;
	char        pname[_str_len_];
	SolVarCell *svCell;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// use smoothed viscosity only if option is activated
	ierr = PetscOptionsGetString(NULL, NULL, "-BFBT_viscositySmoothing", pname, _str_len_, &flg); CHKERRQ(ierr);
	if(flg == PETSC_FALSE) PetscFunctionReturn(0);

	// access context variables
	fs   = jr->fs;
	dsx  = &fs->dsx;
	dsy  = &fs->dsy;
	dsz  = &fs->dsz;

	DynamicRatio 	= 100;								// Default
	delta 			= 200;								// Default
	ierr = PetscOptionsGetScalar(NULL, NULL, "-BFBT_viscositySmoothing_DynamicRatio", &DynamicRatio, &flg); 	CHKERRQ(ierr);
	ierr = PetscOptionsGetScalar(NULL, NULL, "-BFBT_viscositySmoothing_delta", &delta, &flg); 				CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"\n");
	PetscPrintf(PETSC_COMM_WORLD,"-------------------------------------------------------------------------- \n");
	PetscPrintf(PETSC_COMM_WORLD, "   BFBT viscosity pre-smoothing  	: active \n");
	PetscPrintf(PETSC_COMM_WORLD, "   BFBT viscosity contrast     		: %lld \n", (LLD)DynamicRatio);
	PetscPrintf(PETSC_COMM_WORLD, "   BFBT smoothing exponential decay	: %lld \n", (LLD)delta);
	PetscPrintf(PETSC_COMM_WORLD,"-------------------------------------------------------------------------- \n");

	etamin 			= 1.0/sqrt(DynamicRatio);
	etamax 			= 	  sqrt(DynamicRatio);
	radius 			= 0.05;

	// Loop over central points
	iter = 0;
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);
	START_STD_LOOP
	{
		// get coordinates of actual cell center
		coords[0] = dsx->ccoor[iter];
		coords[1] = dsy->ccoor[iter];
		coords[2] = dsz->ccoor[iter];

		ngeom = 0;
		IndicatorValue = 1.0;
		for(iSpheres = 1; iSpheres < nSpheres; iSpheres++){

			GET_GEOM(sphere, geom, ngeom, _max_geom_);	// get actual sphere
			center[0]   = sphere->center[0];
			center[1]   = sphere->center[1];			// get sphere center
			center[2]   = sphere->center[2];
			radius 		= sphere->radius;				// get sphere radius

			VEC_ABS(vecabs, center, coords);
			maximum = vecabs - radius;					// 			 n
			if(maximum < 0.0){maximum = 0.0;}			// Ind(x):=product[1-exp(-delta*max(0,|center(i)-x|-radius)^2)]
			maximum = maximum * maximum;				//			i=1
			expfunc = exp(-delta * maximum);
			IndicatorValue = IndicatorValue * (1 - expfunc);
		}

		svCell = &jr->svCell[iter++];
		svCell->eta_smoothed = (etamax-etamin) * (1-IndicatorValue) + etamin; // store viscosity
		if(iter%5500==0) PetscPrintf(PETSC_COMM_WORLD, "%lld \n", svCell->eta_smoothed); 							// storing didnt work...
		if(iter%5500==0) PetscPrintf(PETSC_COMM_WORLD, "%lld \n", (etamax-etamin) * (1-IndicatorValue) + etamin); 	//
	}
	END_STD_LOOP

	PetscFunctionReturn(0);
}



//---------------------------------------------------------------------------


//---------------------------------------------------------------------------


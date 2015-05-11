//---------------------------------------------------------------------------
//......................   TEMPERATURE FUNCTIONS   ..........................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "constEq.h"
#include "Utils.h"
//---------------------------------------------------------------------------

#define SCATTER_FIELD(da, vec, FIELD) \
	ierr = DMDAGetCorners (da, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr); \
	ierr = DMDAVecGetArray(da, vec, &buff); CHKERRQ(ierr); \
	iter = 0; \
	START_STD_LOOP \
		FIELD \
	END_STD_LOOP \
	ierr = DMDAVecRestoreArray(da, vec, &buff); CHKERRQ(ierr); \
	LOCAL_TO_LOCAL(da, vec)

#define GET_KC \
	GetTempParam(numPhases, phases, jr->svCell[iter++].phRat, &kc, NULL, NULL); \
	buff[k][j][i] = kc;

#define GET_HRXY buff[k][j][i] = jr->svXYEdge[iter++].svDev.Hr;
#define GET_HRYZ buff[k][j][i] = jr->svYZEdge[iter++].svDev.Hr;
#define GET_HRXZ buff[k][j][i] = jr->svXZEdge[iter++].svDev.Hr;

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResCreateTempParam"
PetscErrorCode JacResCreateTempParam(JacRes *jr)
{
	// setup temperature parameters

	FDSTAG *fs;
	const PetscInt *lx, *ly, *lz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = jr->fs;

	// get cell center grid partitioning
	ierr = DMDAGetOwnershipRanges(fs->DA_CEN, &lx, &ly, &lz); CHKERRQ(ierr);

	// create temperature DMDA
	ierr = DMDACreate3d(PETSC_COMM_WORLD,
		DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
		DMDA_STENCIL_STAR,
		fs->dsx.tcels, fs->dsy.tcels, fs->dsz.tcels,
		fs->dsx.nproc, fs->dsy.nproc, fs->dsz.nproc,
		1, 1, lx, ly, lz, &jr->DA_T); CHKERRQ(ierr);

	// create temperature preconditioner matrix
	ierr = DMCreateMatrix(jr->DA_T, &jr->Att); CHKERRQ(ierr);

	// temperature
	ierr = DMCreateGlobalVector(fs->DA_CEN, &jr->gT); CHKERRQ(ierr);
	ierr = DMCreateLocalVector (fs->DA_CEN, &jr->lT); CHKERRQ(ierr);

	// energy residual
	ierr = DMCreateGlobalVector(fs->DA_CEN, &jr->ge);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResDestroyTempParam"
PetscErrorCode JacResDestroyTempParam(JacRes *jr)
{
	// destroy temperature parameters

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// temperature parameters
	ierr = DMDestroy (&jr->DA_T); CHKERRQ(ierr);
	ierr = MatDestroy(&jr->Att);  CHKERRQ(ierr);

	ierr = VecDestroy(&jr->gT);   CHKERRQ(ierr);
	ierr = VecDestroy(&jr->lT);   CHKERRQ(ierr);

	ierr = VecDestroy(&jr->ge);   CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResInitTemp"
PetscErrorCode JacResInitTemp(JacRes *jr)
{
	// initialize temperature from markers

	PetscScalar *T;
	FDSTAG      *fs;
	PetscInt     jj, n;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs = jr->fs;

	// copy temperatures from context storage to global vector
	ierr = VecGetArray(jr->gT, &T);  CHKERRQ(ierr);

	for(jj = 0, n = fs->nCells; jj < n; jj++) T[jj] = jr->svCell[jj].svBulk.Tn;

	ierr = VecRestoreArray(jr->gT, &T);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResCopyTemp"
PetscErrorCode JacResCopyTemp(JacRes *jr)
{
	// copy temperature from global to local vectors, enforce boundary constraints

	FDSTAG      *fs;
	BCCtx       *bc;
	PetscInt    mcx, mcy, mcz;
	PetscInt    I, J, K, fi, fj, fk;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***lT, ***bcT;
	PetscScalar *T, pmdof;
	PetscScalar *vals;
	PetscInt    num, *list;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs  =  jr->fs;
	bc  =  jr->bc;

	// initialize maximal index in all directions
	mcx = fs->dsx.tcels - 1;
	mcy = fs->dsy.tcels - 1;
	mcz = fs->dsz.tcels - 1;

	// access vector
	ierr = VecGetArray(jr->gT,  &T);   CHKERRQ(ierr);

	//=================================
	// enforce single point constraints
	//=================================

	// temperature
	num   = bc->tNumSPC;
	list  = bc->tSPCList;
	vals  = bc->tSPCVals;

	for(i = 0; i < num; i++) T[list[i]] = vals[i];

	ierr = VecRestoreArray(jr->gT, &T); CHKERRQ(ierr);

	// fill local (ghosted) version of solution vector
	GLOBAL_TO_LOCAL(fs->DA_CEN, jr->gT, jr->lT)

	// access local solution vector
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lT, &lT);  CHKERRQ(ierr);

	// access boundary constraints vector
	ierr = DMDAVecGetArray(fs->DA_CEN, bc->bcT, &bcT); CHKERRQ(ierr);

	//==============================
	// enforce two-point constraints
	//==============================

	//-----------------------------
	// central points (temperature)
	//-----------------------------
	GET_CELL_RANGE_GHOST_INT(nx, sx, fs->dsx)
	GET_CELL_RANGE_GHOST_INT(ny, sy, fs->dsy)
	GET_CELL_RANGE_GHOST_INT(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		pmdof = lT[k][j][i];

		I = i; fi = 0;
		J = j; fj = 0;
		K = k; fk = 0;

		if(i == 0)   { fi = 1; I = i-1; SET_TPC(bcT, lT, k, j, I, pmdof) }
		if(i == mcx) { fi = 1; I = i+1; SET_TPC(bcT, lT, k, j, I, pmdof) }
		if(j == 0)   { fj = 1; J = j-1; SET_TPC(bcT, lT, k, J, i, pmdof) }
		if(j == mcy) { fj = 1; J = j+1; SET_TPC(bcT, lT, k, J, i, pmdof) }
		if(k == 0)   { fk = 1; K = k-1; SET_TPC(bcT, lT, K, j, i, pmdof) }
		if(k == mcz) { fk = 1; K = k+1; SET_TPC(bcT, lT, K, j, i, pmdof) }

		if(fi*fj)    SET_EDGE_CORNER(n, lT, k, J, I, k, j, i, pmdof)
		if(fi*fk)    SET_EDGE_CORNER(n, lT, K, j, I, k, j, i, pmdof)
		if(fj*fk)    SET_EDGE_CORNER(n, lT, K, J, i, k, j, i, pmdof)
		if(fi*fj*fk) SET_EDGE_CORNER(n, lT, K, J, I, k, j, i, pmdof)
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lT, &lT);  CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fs->DA_CEN, bc->bcT, &bcT); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetTempRes"
PetscErrorCode JacResGetTempRes(JacRes *jr)
{
	// compute temperature residual vector

	FDSTAG     *fs;
	SolVarCell *svCell;
	SolVarDev  *svDev;
	SolVarBulk *svBulk;
	Material_t *phases;
	PetscInt    iter, numPhases;
	PetscInt    Ip1, Im1, Jp1, Jm1, Kp1, Km1;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, my, mz;
 	PetscScalar bkx, fkx, bky, fky, bkz, fkz;
	PetscScalar bdx, fdx, bdy, fdy, bdz, fdz;
	PetscScalar bqx, fqx, bqy, fqy, bqz, fqz;
 	PetscScalar dx, dy, dz;
	PetscScalar kc, rho, Cp, A, Tc, Tn, dt, Hr;
	PetscScalar ***ge, ***T, ***lk, ***hxy, ***hxz, ***hyz, ***buff;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access residual context variables
	fs        = jr->fs;
	numPhases = jr->numPhases; // number phases
	phases    = jr->phases;    // phase parameters
	dt        = jr->ts.dt;     // time step

	// initialize maximum cell index in all directions
	mx = fs->dsx.tcels - 1;
	my = fs->dsy.tcels - 1;
	mz = fs->dsz.tcels - 1;

	SCATTER_FIELD(fs->DA_CEN, jr->ldxx, GET_KC)
	SCATTER_FIELD(fs->DA_XY,  jr->ldxy, GET_HRXY)
	SCATTER_FIELD(fs->DA_XZ,  jr->ldxz, GET_HRYZ)
	SCATTER_FIELD(fs->DA_YZ,  jr->ldyz, GET_HRXZ)

	// access work vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ge,   &ge);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lT,   &T);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldxx, &lk);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XY,  jr->ldxy, &hxy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XZ,  jr->ldxz, &hxz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_YZ,  jr->ldyz, &hyz); CHKERRQ(ierr);

	//---------------
	// central points
	//---------------
	iter = 0;
	ierr = DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz); CHKERRQ(ierr);

	START_STD_LOOP
	{
		// access solution variables
		svCell = &jr->svCell[iter++];
		svDev  = &svCell->svDev;
		svBulk = &svCell->svBulk;

		// access
		Tc  = T[k][j][i];  // current temperature
		Tn  = svBulk->Tn;  // temperature history
		rho = svBulk->rho; // effective density

		// conductivity, heat capacity, radiogenic heat production
		GetTempParam(numPhases, phases, svCell->phRat, &kc, &Cp, &A);

		// shear heating term
		Hr = svDev->Hr +
		(hxy[k][j][i] + hxy[k][j+1][i] + hxy[k][j][i+1] + hxy[k][j+1][i+1] +
		 hxz[k][j][i] + hxz[k+1][j][i] + hxz[k][j][i+1] + hxz[k+1][j][i+1] +
		 hyz[k][j][i] + hyz[k+1][j][i] + hyz[k][j+1][i] + hyz[k+1][j+1][i])/4.0;

		// check index bounds
		Ip1 = i+1; if(Ip1 > mx) Ip1--;
		Im1 = i-1; if(Im1 < 0)  Im1++;
		Jp1 = j+1; if(Jp1 > my) Jp1--;
		Jm1 = j-1; if(Jm1 < 0)  Jm1++;
		Kp1 = k+1; if(Kp1 > mz) Kp1--;
		Km1 = k-1; if(Km1 < 0)  Km1++;

		// compute average conductivities
		bkx = (kc + lk[k][j][Im1])/2.0;      fkx = (kc + lk[k][j][Ip1])/2.0;
		bky = (kc + lk[k][Jm1][i])/2.0;      fky = (kc + lk[k][Jp1][i])/2.0;
		bkz = (kc + lk[Km1][j][i])/2.0;      fkz = (kc + lk[Kp1][j][i])/2.0;

		// get mesh steps
		bdx = SIZE_NODE(i, sx, fs->dsx);     fdx = SIZE_NODE(i+1, sx, fs->dsx);
		bdy = SIZE_NODE(j, sy, fs->dsy);     fdy = SIZE_NODE(j+1, sy, fs->dsy);
		bdz = SIZE_NODE(k, sz, fs->dsz);     fdz = SIZE_NODE(k+1, sz, fs->dsz);

		// compute heat fluxes
		bqx = bkx*(Tc - T[k][j][i-1])/bdx;   fqx = fkx*(T[k][j][i+1] - Tc)/fdx;
		bqy = bky*(Tc - T[k][j-1][i])/bdy;   fqy = fky*(T[k][j+1][i] - Tc)/fdy;
		bqz = bkz*(Tc - T[k-1][j][i])/bdz;   fqz = fkz*(T[k+1][j][i] - Tc)/fdz;

		// get mesh steps
		dx = SIZE_CELL(i, sx, fs->dsx);
		dy = SIZE_CELL(j, sy, fs->dsy);
		dz = SIZE_CELL(k, sz, fs->dsz);

		// original balance equation:

		// rho*Cp*(Tc - Tn)/dt = (fqx - bqx)/dx + (fqy - bqy)/dy + (fqz - bqz)/dz + Hr + rho*A

		// to get positive diagonal in the preconditioner matrix
		// put right hand side to the left, which gives the following:

		ge[k][j][i] = rho*Cp*(Tc - Tn)/dt - (fqx - bqx)/dx - (fqy - bqy)/dy - (fqz - bqz)/dz - Hr - rho*A;
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ge,   &ge);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lT,   &T);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldxx, &lk);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XY,  jr->ldxy, &hxy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XZ,  jr->ldxz, &hxz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_YZ,  jr->ldyz, &hyz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "JacResGetTempMat"
PetscErrorCode JacResGetTempMat(JacRes *jr)
{


	/*


	PetscErrorCode  MatSetValuesStencil(
			Mat        mat,
			PetscInt   m,
			const      MatStencil idxm[],
			PetscInt   n,
			const      MatStencil idxn[],
			const      PetscScalar v[],
			InsertMode addv)

	typedef struct {
	  PetscInt k,j,i,c;
	} MatStencil;





	//======================================================================
	// Assemble temperature equation preconditioning matrix
	//======================================================================

	FDSTAG      *fs;
	BCCtx       *bc;
	DOFIndex    *dof;
	PMatMono    *P;
	PetscInt    idx[7];
	PetscScalar v[49];
	PetscInt    iter, i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar eta, rho, IKdt, diag, pgamma, pt, dt, fssa, *grav;
	PetscScalar dx, dy, dz, bdx, fdx, bdy, fdy, bdz, fdz;
	PetscScalar ***ivx, ***ivy, ***ivz, ***ip;
	PetscScalar ***bcvx, ***bcvy, ***bcvz, ***bcp;
	PetscInt    pdofidx[7];
	PetscScalar cf[7];


	PetscScalar Cp, lambda;


	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access contexts
	jr     = pm->jr;
	fs     = jr->fs;
	bc     = jr->bc;
	dof    = &fs->dof;
	P      = (PMatMono*)pm->data;

	// get density gradient stabilization parameters
	dt   = jr->ts.dt; // time step


	// clear matrix coefficients
	ierr = MatZeroEntries(P->A); CHKERRQ(ierr);

	// access index vectors
	ierr = DMDAVecGetArray(fs->DA_X,   dof->ivx,  &ivx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   dof->ivy,  &ivy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   dof->ivz,  &ivz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, dof->ip,   &ip);   CHKERRQ(ierr);

	// access boundary constraint vectors
	ierr = DMDAVecGetArray(fs->DA_X,   bc->bcvx,  &bcvx);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   bc->bcvy,  &bcvy);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   bc->bcvz,  &bcvz);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, bc->bcp,   &bcp);   CHKERRQ(ierr);

	//---------------
	// central points
	//---------------

	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{

		// get heat capacity and conductivity

//		Cp lambda



		// get density, shear & inverse bulk viscosities
		eta  = jr->svCell[iter].svDev.eta;
		IKdt = jr->svCell[iter].svBulk.IKdt;
		rho  = jr->svCell[iter].svBulk.rho;

		iter++;

		// get mesh steps
		dx = SIZE_CELL(i, sx, fs->dsx);
		dy = SIZE_CELL(j, sy, fs->dsy);
		dz = SIZE_CELL(k, sz, fs->dsz);

		// get mesh steps for the backward and forward derivatives
		bdx = SIZE_NODE(i, sx, fs->dsx);   fdx = SIZE_NODE(i+1, sx, fs->dsx);
		bdy = SIZE_NODE(j, sy, fs->dsy);   fdy = SIZE_NODE(j+1, sy, fs->dsy);
		bdz = SIZE_NODE(k, sz, fs->dsz);   fdz = SIZE_NODE(k+1, sz, fs->dsz);

		// compute penalty term
		pt = -1.0/(pgamma*eta);

		// get pressure diagonal element (with penalty)
		diag = -IKdt + pt;

		// compute local matrix
		pm->getStiffMat(eta, diag, v, dx, dy, dz, fdx, fdy, fdz, bdx, bdy, bdz);

		// compute density gradient stabilization terms
		addDensGradStabil(fssa, v, rho, dt, grav, fdx, fdy, fdz, bdx, bdy, bdz);

		// get global indices of the points:
		// vx_(i), vx_(i+1), vy_(j), vy_(j+1), vz_(k), vz_(k+1), p
		idx[0] = (PetscInt) ivx[k][j][i];
		idx[1] = (PetscInt) ivx[k][j][i+1];
		idx[2] = (PetscInt) ivy[k][j][i];
		idx[3] = (PetscInt) ivy[k][j+1][i];
		idx[4] = (PetscInt) ivz[k][j][i];
		idx[5] = (PetscInt) ivz[k+1][j][i];
		idx[6] = (PetscInt) ip[k][j][i];

		// get boundary constraints
		pdofidx[0] = -1;   cf[0] = bcvx[k][j][i];
		pdofidx[1] = -1;   cf[1] = bcvx[k][j][i+1];
		pdofidx[2] = -1;   cf[2] = bcvy[k][j][i];
		pdofidx[3] = -1;   cf[3] = bcvy[k][j+1][i];
		pdofidx[4] = -1;   cf[4] = bcvz[k][j][i];
		pdofidx[5] = -1;   cf[5] = bcvz[k+1][j][i];
		pdofidx[6] = -1;   cf[6] = bcp[k][j][i];

		// constrain local matrix
		constrLocalMat(7, pdofidx, cf, v);

		// update global & penalty compensation matrices
		ierr = MatSetValues(P->A, 7, idx, 7, idx, v,  ADD_VALUES);    CHKERRQ(ierr);
		ierr = MatSetValue (P->M, idx[6], idx[6], pt, INSERT_VALUES); CHKERRQ(ierr);

	}
	END_STD_LOOP


	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X,   dof->ivx,  &ivx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   dof->ivy,  &ivy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   dof->ivz,  &ivz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, dof->ip,   &ip);  CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(fs->DA_X,   bc->bcvx,  &bcvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   bc->bcvy,  &bcvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   bc->bcvz,  &bcvz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, bc->bcp,   &bcp);  CHKERRQ(ierr);

	// assemble velocity-pressure matrix, remove constrained rows
	ierr = MatAIJAssemble(P->A, bc->numSPC, bc->SPCList, 1.0); CHKERRQ(ierr);
	ierr = MatAIJAssemble(P->M, bc->numSPC, bc->SPCList, 0.0); CHKERRQ(ierr);
*/
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

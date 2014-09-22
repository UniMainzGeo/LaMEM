//---------------------------------------------------------------------------
//..................   GALERKIN GEOMETRIC MULTIGRID   .......................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "bc.h"
#include "JacRes.h"
//#include "lsolve.h"
#include "matrix.h"
#include "multigrid.h"
#include "Utils.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MGCheckGrid"
PetscErrorCode MGCheckGrid(FDSTAG *fs, PetscInt *ncels)
{
	PetscInt  nx, ny, nz;
	PetscInt  lx, ly, lz;

	PetscFunctionBegin;

	// compute uniform local grid size in every direction (same on all processors)
	if(fs->dsx.tcels % fs->dsx.nproc) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "ERROR! number x-cells is not a multiple of number x-processors\n");
	if(fs->dsy.tcels % fs->dsy.nproc) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "ERROR! number y-cells is not a multiple of number y-processors\n");
	if(fs->dsz.tcels % fs->dsz.nproc) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "ERROR! number z-cells is not a multiple of number z-processors\n");

	nx = fs->dsx.tcels/fs->dsx.nproc;
	ny = fs->dsy.tcels/fs->dsy.nproc;
	nz = fs->dsz.tcels/fs->dsz.nproc;

	// get actual grid size in every direction
	lx = fs->dsx.ncels;
	ly = fs->dsy.ncels;
	lz = fs->dsz.ncels;

	// 1. local grid size must be constant in all directions
	// 2. local grid size must be power of two
	// 3. local grid size must be constant on all processors

	if(nx != ny || nx != nz)             SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "ERROR! local grid size is not constant in all directions\n");
	if(!IS_POWER_OF_TWO(nx))             SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "ERROR! Local grid size is not power of two\n");
	if(lx != nx || ly != ny || lz != nz) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "ERROR! local grid size is not constant on all processors\n");

	// output local grid size
	if(ncels) (*ncels) = nx;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MGCtxCreate"
PetscErrorCode MGCtxCreate(MGCtx *mg, FDSTAG *fs, BCCtx *bc, PC pc, idxtype idxmod)
{
	PetscInt  i, l, ncors, nlevels;
	FDSTAG   *fine, *cors;
	char      pc_type[PETSC_MAX_PATH_LEN];
	PetscBool opt_set;
	PetscInt  Nx, Ny, Nz;
	PetscInt  Px, Py, Pz;
	PetscInt  ncels;
	DOFIndex  *idfine, *idcors;

	PetscErrorCode 	 ierr;
	PetscFunctionBegin;

	// get preconditioner type
	ierr = PetscOptionsGetString(NULL, "-gmg_pc_type", pc_type, PETSC_MAX_PATH_LEN, &opt_set); CHKERRQ(ierr);

	// check whether multigrid is requested
	if(opt_set != PETSC_TRUE || strcmp(pc_type, "mg")) PetscFunctionReturn(0);

	// check multigrid mesh restrictions, store local grid size
	ierr = MGCheckGrid(fs, &ncels); CHKERRQ(ierr);

	// compute maximum number of coarsening steps (grids)
	ncors = (PetscInt)round(log((double)ncels)/M_LN2);

	// check for command line overrides
	ierr =  PetscOptionsGetInt(NULL, "-gmg_pc_mg_levels", &nlevels, &opt_set); CHKERRQ(ierr);

	if(opt_set == PETSC_TRUE)
	{	if(nlevels < ncors+1) ncors = nlevels-1;
		else SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "ERROR! too many multigrid layers.\n");
	}

	// store actual number of coarsening steps
	mg->ncors = ncors;

	// allocate storage for grids & matrices
	ierr = PetscMalloc(sizeof(FDSTAG)*(size_t)mg->ncors, &mg->mgfs); CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(Mat)   *(size_t)mg->ncors, &mg->R);    CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(Mat)   *(size_t)mg->ncors, &mg->P);    CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(BCCtx) *(size_t)mg->ncors, &mg->mgbc); CHKERRQ(ierr);

	// get number of processors
	Px = fs->dsx.nproc;
	Py = fs->dsy.nproc;
	Pz = fs->dsz.nproc;

	// get fine grid size
	Nx = fs->dsx.tcels;
	Ny = fs->dsy.tcels;
	Nz = fs->dsz.tcels;

	// setup grids & bc contexts
	for(i = 0; i < mg->ncors; i++)
	{
		// coarsen uniformly in all directions
		Nx /= 2;
		Ny /= 2;
		Nz /= 2;

		// create mesh
		ierr = FDSTAGCreate(&mg->mgfs[i], Nx+1, Ny+1, Nz+1, Px, Py, Pz); CHKERRQ(ierr);

		// create bc context
		ierr = BCCreate(&mg->mgbc[i], &mg->mgfs[i]); CHKERRQ(ierr);

		// setup bc context
		ierr = BCInit(&mg->mgbc[i], &mg->mgfs[i], idxmod); CHKERRQ(ierr);


		// check multigrid restrictions (debug)
		ierr = MGCheckGrid(&mg->mgfs[i], NULL); CHKERRQ(ierr);
	}

	// create Galerkin multigrid preconditioner
	ierr = PCSetOptionsPrefix(pc, "gmg_");        CHKERRQ(ierr);
	ierr = PCSetType(pc, PCMG);                   CHKERRQ(ierr);
	ierr = PCMGSetLevels(pc, mg->ncors+1, NULL);  CHKERRQ(ierr);
	ierr = PCMGSetType(pc, PC_MG_MULTIPLICATIVE); CHKERRQ(ierr);
	ierr = PCMGSetGalerkin(pc, PETSC_TRUE);       CHKERRQ(ierr);
	ierr = PCSetFromOptions(pc);                  CHKERRQ(ierr);

	// set fine grid
	fine = fs;

	// create & preallocate matrices
	for(i = 0, l = mg->ncors; i < mg->ncors; i++, l--)
	{
		// set coarse grid
		cors = &mg->mgfs[i];

		if(idxmod == IDXCOUPLED)   { idfine = &fine->cdof; idcors = &cors->cdof; }
		if(idxmod == IDXUNCOUPLED) { idfine = &fine->udof; idcors = &cors->udof; }

		// constant size preallocation (MAKE SURE SPACE IS ENOUGH)
		ierr = PMatCreate(idcors->numdof, idfine->numdof, 12, NULL, 4, NULL, &mg->R[i]); CHKERRQ(ierr);
		ierr = PMatCreate(idfine->numdof, idcors->numdof, 8,  NULL, 7, NULL, &mg->P[i]); CHKERRQ(ierr);

		ierr = PCMGSetRestriction  (pc, l, mg->R[i]);
		ierr = PCMGSetInterpolation(pc, l, mg->P[i]);

		// set fine grid for next step
		fine = cors;
	}

	// store fine grid contexts
	mg->fs = fs;
	mg->bc = bc;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MGCtxDestroy"
PetscErrorCode MGCtxDestroy(MGCtx *mg)
{
	PetscInt i;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	for(i = 0; i < mg->ncors; i++)
	{
		ierr = FDSTAGDestroy(&mg->mgfs[i]); CHKERRQ(ierr);
		ierr = MatDestroy   (&mg->R[i]);    CHKERRQ(ierr);
		ierr = MatDestroy   (&mg->P[i]);    CHKERRQ(ierr);
		ierr = BCDestroy    (&mg->mgbc[i]); CHKERRQ(ierr);

	}

	ierr = PetscFree(mg->mgfs); CHKERRQ(ierr);
	ierr = PetscFree(mg->R);    CHKERRQ(ierr);
	ierr = PetscFree(mg->P);    CHKERRQ(ierr);
	ierr = PetscFree(mg->mgbc); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MGCtxSetup"
PetscErrorCode MGCtxSetup(MGCtx *mg, idxtype idxmod)
{

	// Matrices are re-assembled here, just in case
	// they will be made matrix- distance- dependent.
	// Currently they depend only on boundary conditions,
	// so changing boundary condition would also require re-assembly.

	PetscInt  i;
	FDSTAG   *fine,   *cors;
	BCCtx    *bcfine, *bccors;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// set fine grid
	fine   = mg->fs;
	bcfine = mg->bc;

	// create & preallocate matrices
	for(i = 0; i < mg->ncors; i++)
	{
		// set coarse grid
		cors   = &mg->mgfs[i];
		bccors = &mg->mgbc[i];

		// assemble matrices
		ierr = SetupRestrictStep(mg->R[i], cors, fine, bccors, idxmod); CHKERRQ(ierr);
		ierr = SetupProlongStep (mg->P[i], fine, cors, bcfine, idxmod); CHKERRQ(ierr);

		// set fine grid for next step
		fine   = cors;
		bcfine = bccors;

	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MGCtxSetDiagOnLevels"
PetscErrorCode MGCtxSetDiagOnLevels(MGCtx *mg, PC pcmg)
{
	PetscInt  i, l;
	KSP       ksp;
	PC        pc;
	Mat       A;
	BCCtx    *bc;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// set dummy coarse solver
	ierr = PCMGGetCoarseSolve(pcmg, &ksp); CHKERRQ(ierr);
	ierr = KSPSetType(ksp, KSPPREONLY);    CHKERRQ(ierr);
	ierr = KSPGetPC(ksp, &pc);             CHKERRQ(ierr);
	ierr = PCSetType(pc, PCNONE);          CHKERRQ(ierr);

	// setup operators
	ierr = PCSetUp(pcmg); CHKERRQ(ierr);

	// constrain operators on all levels except the coarsest & finest
	for(i = 0, l = mg->ncors-1; i < mg->ncors-1; i++, l--)
	{
		// get boundary condition context
		bc = &mg->mgbc[i];

		// access smoother on every level
		ierr = PCMGGetSmoother(pcmg, l, &ksp); CHKERRQ(ierr);

		// get matrix
		ierr = KSPGetOperators(ksp, &A, NULL); CHKERRQ(ierr);

		// set diagonals for constrained DOF
		ierr = MatZeroRowsColumns(A, bc->numSPC, bc->SPCList, 1.0, NULL, NULL); CHKERRQ(ierr);
	}

	// constrain coarse operator
	if(mg->ncors)
	{
		// get bc context for coarse grid
		bc = &mg->mgbc[mg->ncors-1];

		// constrain coarse operators
		ierr = PCMGGetCoarseSolve(pcmg, &ksp); CHKERRQ(ierr);
		ierr = KSPGetOperators(ksp, &A, NULL); CHKERRQ(ierr);
		ierr = MatZeroRowsColumns(A, bc->numSPC, bc->SPCList, 1.0, NULL, NULL); CHKERRQ(ierr);

		// setup new solver
		ierr = KSPSetOptionsPrefix(ksp, "crs_");  CHKERRQ(ierr);
		ierr = KSPSetFromOptions(ksp);            CHKERRQ(ierr);
		ierr = KSPSetUp(ksp);                     CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "SetupRestrictStep"
PetscErrorCode SetupRestrictStep(Mat R, FDSTAG *cors, FDSTAG *fine, BCCtx *bccors, idxtype idxmod)
{
	PetscScalar v[12];
	PetscInt    idx[12];
	PetscInt    mx, my, mz;
	PetscInt    row, I, J, K, Ip1, Im1, Jp1, Jm1, Kp1, Km1;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***ivx, ***ivy, ***ivz, ***ip;
	DOFIndex    *idfine, *idcors;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	if(idxmod == IDXCOUPLED)   { idfine = &fine->cdof; idcors = &cors->cdof; }
	if(idxmod == IDXUNCOUPLED) { idfine = &fine->udof; idcors = &cors->udof; }

	// clear restriction matrix coefficients
	ierr = MatZeroEntries(R); CHKERRQ(ierr);

	// access index vectors in fine grid
	ierr = DMDAVecGetArray(fine->DA_X,   idfine->ivx, &ivx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_Y,   idfine->ivy, &ivy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_Z,   idfine->ivz, &ivz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fine->DA_CEN, idfine->ip,  &ip);  CHKERRQ(ierr);

	// get fine grid index bounds
	mx = fine->dsx.tnods;
	my = fine->dsy.tnods;
	mz = fine->dsz.tnods;

	// get global index of the first row in coarse grid
	row = idcors->istart;

	// set velocity weights
	v[0]  = 1.0/16.0;
	v[1]  = 1.0/16.0;
	v[2]  = 1.0/16.0;
	v[3]  = 1.0/16.0;
	v[4]  = 1.0/8.0;
	v[5]  = 1.0/8.0;
	v[6]  = 1.0/8.0;
	v[7]  = 1.0/8.0;
	v[8]  = 1.0/16.0;
	v[9]  = 1.0/16.0;
	v[10] = 1.0/16.0;
	v[11] = 1.0/16.0;

	//-----------------------
	// X-points (coarse grid)
	//-----------------------
	GET_NODE_RANGE(nx, sx, cors->dsx)
	GET_CELL_RANGE(ny, sy, cors->dsy)
	GET_CELL_RANGE(nz, sz, cors->dsz)

	START_STD_LOOP
	{
		// get fine grid indices
		I = 2*i;
		J = 2*j;
		K = 2*k;

		Ip1 = I+1; if(Ip1 == mx) Ip1 = I;
		Im1 = I-1; if(Im1 == -1) Im1 = I;

		idx[0]  = (PetscInt)ivx[K  ][J  ][Im1];
		idx[1]  = (PetscInt)ivx[K  ][J+1][Im1];
		idx[2]  = (PetscInt)ivx[K+1][J  ][Im1];
		idx[3]  = (PetscInt)ivx[K+1][J+1][Im1];
		idx[4]  = (PetscInt)ivx[K  ][J  ][I  ];
		idx[5]  = (PetscInt)ivx[K  ][J+1][I  ];
		idx[6]  = (PetscInt)ivx[K+1][J  ][I  ];
		idx[7]  = (PetscInt)ivx[K+1][J+1][I  ];
		idx[8]  = (PetscInt)ivx[K  ][J  ][Ip1];
		idx[9]  = (PetscInt)ivx[K  ][J+1][Ip1];
		idx[10] = (PetscInt)ivx[K+1][J  ][Ip1];
		idx[11] = (PetscInt)ivx[K+1][J+1][Ip1];

		// store full matrix row
		ierr = MatSetValues(R, 1, &row, 12, idx, v, ADD_VALUES); CHKERRQ(ierr);

		// increment row number
		row++;
	}
	END_STD_LOOP

	//-----------------------
	// Y-points (coarse grid)
	//-----------------------
	GET_CELL_RANGE(nx, sx, cors->dsx)
	GET_NODE_RANGE(ny, sy, cors->dsy)
	GET_CELL_RANGE(nz, sz, cors->dsz)

	START_STD_LOOP
	{
		// get fine grid indices
		I = 2*i;
		J = 2*j;
		K = 2*k;

		Jp1 = J+1; if(Jp1 == my) Jp1 = J;
		Jm1 = J-1; if(Jm1 == -1) Jm1 = J;

		idx[0]  = (PetscInt)ivy[K  ][Jm1][I  ];
		idx[1]  = (PetscInt)ivy[K  ][Jm1][I+1];
		idx[2]  = (PetscInt)ivy[K+1][Jm1][I  ];
		idx[3]  = (PetscInt)ivy[K+1][Jm1][I+1];
		idx[4]  = (PetscInt)ivy[K  ][J  ][I  ];
		idx[5]  = (PetscInt)ivy[K  ][J  ][I+1];
		idx[6]  = (PetscInt)ivy[K+1][J  ][I  ];
		idx[7]  = (PetscInt)ivy[K+1][J  ][I+1];
		idx[8]  = (PetscInt)ivy[K  ][Jp1][I  ];
		idx[9]  = (PetscInt)ivy[K  ][Jp1][I+1];
		idx[10] = (PetscInt)ivy[K+1][Jp1][I  ];
    	idx[11] = (PetscInt)ivy[K+1][Jp1][I+1];

		// store full matrix row
		ierr = MatSetValues(R, 1, &row, 12, idx, v, ADD_VALUES); CHKERRQ(ierr);

		// increment row number
		row++;
	}
	END_STD_LOOP

	//-----------------------
	// Z-points (coarse grid)
	//-----------------------
	GET_CELL_RANGE(nx, sx, cors->dsx)
	GET_CELL_RANGE(ny, sy, cors->dsy)
	GET_NODE_RANGE(nz, sz, cors->dsz)

	START_STD_LOOP
	{
		// get fine grid indices
		I = 2*i;
		J = 2*j;
		K = 2*k;

		Kp1 = K+1; if(Kp1 == mz) Kp1 = K;
		Km1 = K-1; if(Km1 == -1) Km1 = K;

		idx[0]  = (PetscInt)ivz[Km1][J  ][I  ];
		idx[1]  = (PetscInt)ivz[Km1][J  ][I+1];
		idx[2]  = (PetscInt)ivz[Km1][J+1][I  ];
		idx[3]  = (PetscInt)ivz[Km1][J+1][I+1];
		idx[4]  = (PetscInt)ivz[K  ][J  ][I  ];
		idx[5]  = (PetscInt)ivz[K  ][J  ][I+1];
		idx[6]  = (PetscInt)ivz[K  ][J+1][I  ];
		idx[7]  = (PetscInt)ivz[K  ][J+1][I+1];
		idx[8]  = (PetscInt)ivz[Kp1][J  ][I  ];
		idx[9]  = (PetscInt)ivz[Kp1][J  ][I+1];
		idx[10] = (PetscInt)ivz[Kp1][J+1][I  ];
    	idx[11] = (PetscInt)ivz[Kp1][J+1][I+1];

		// store full matrix row
		ierr = MatSetValues(R, 1, &row, 12, idx, v, ADD_VALUES); CHKERRQ(ierr);

		// increment row number
		row++;
	}
	END_STD_LOOP

	if(idxmod == IDXCOUPLED)
	{
		// set pressure weights
		v[0] = 1.0/8.0;
		v[1] = 1.0/8.0;
		v[2] = 1.0/8.0;
		v[3] = 1.0/8.0;
		v[4] = 1.0/8.0;
		v[5] = 1.0/8.0;
		v[6] = 1.0/8.0;
		v[7] = 1.0/8.0;

		//-----------------------
		// P-points (coarse grid)
		//-----------------------
		GET_CELL_RANGE(nx, sx, cors->dsx)
		GET_CELL_RANGE(ny, sy, cors->dsy)
		GET_CELL_RANGE(nz, sz, cors->dsz)

		START_STD_LOOP
		{
			// get fine grid indices
			I = 2*i;
			J = 2*j;
			K = 2*k;

			idx[0] = (PetscInt)ip[K  ][J  ][I  ];
			idx[1] = (PetscInt)ip[K  ][J  ][I+1];
			idx[2] = (PetscInt)ip[K  ][J+1][I  ];
			idx[3] = (PetscInt)ip[K  ][J+1][I+1];
			idx[4] = (PetscInt)ip[K+1][J  ][I  ];
			idx[5] = (PetscInt)ip[K+1][J  ][I+1];
			idx[6] = (PetscInt)ip[K+1][J+1][I  ];
			idx[7] = (PetscInt)ip[K+1][J+1][I+1];

			// store full matrix row
			ierr = MatSetValues(R, 1, &row, 8, idx, v, ADD_VALUES); CHKERRQ(ierr);

			// increment row number
			row++;

		}
		END_STD_LOOP
	}

	// restore access
	ierr = DMDAVecRestoreArray(fine->DA_X,   idfine->ivx,  &ivx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_Y,   idfine->ivy,  &ivy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_Z,   idfine->ivz,  &ivz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fine->DA_CEN, idfine->ip,   &ip);  CHKERRQ(ierr);

	// assemble restriction matrix
	ierr = PMatAssemble(R, bccors->numSPC, bccors->SPCList);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "SetupProlongStep"
PetscErrorCode SetupProlongStep(Mat P, FDSTAG *fine, FDSTAG *cors, BCCtx *bcfine, idxtype idxmod)
{
	PetscScalar v[8];
	PetscInt    idx[8];
	PetscInt    mx, my, mz;
	PetscInt    row, I, J, K, I1, J1, K1;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***ivx, ***ivy, ***ivz, ***ip;
	DOFIndex    *idcors, *idfine;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	if(idxmod == IDXCOUPLED)   { idcors = &cors->cdof; idfine = &fine->cdof; }
	if(idxmod == IDXUNCOUPLED) { idcors = &cors->udof; idfine = &fine->udof; }

	// clear prolongation matrix coefficients
	ierr = MatZeroEntries(P); CHKERRQ(ierr);

	// access index vectors in coarse grid
	ierr = DMDAVecGetArray(cors->DA_X,   idcors->ivx,  &ivx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cors->DA_Y,   idcors->ivy,  &ivy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cors->DA_Z,   idcors->ivz,  &ivz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cors->DA_CEN, idcors->ip,   &ip);  CHKERRQ(ierr);

	// get coarse grid index bounds
	mx = cors->dsx.tcels;
	my = cors->dsy.tcels;
	mz = cors->dsz.tcels;

	// get global index of the first row in fine grid
	row = idfine->istart;

	// set velocity weights
	v[0] = 9.0/32.0;
	v[1] = 3.0/32.0;
	v[2] = 3.0/32.0;
	v[3] = 1.0/32.0;
	v[4] = 9.0/32.0;
	v[5] = 3.0/32.0;
	v[6] = 3.0/32.0;
	v[7] = 1.0/32.0;

	//---------------------
	// X-points (fine grid)
	//---------------------
	GET_NODE_RANGE(nx, sx, fine->dsx)
	GET_CELL_RANGE(ny, sy, fine->dsy)
	GET_CELL_RANGE(nz, sz, fine->dsz)

	START_STD_LOOP
	{
		// get coarse grid indices
		I = i/2;
		J = j/2;
		K = k/2;

		if(i % 2) I1 = I+1; else I1 = I;

		if(j % 2) J1 = J+1; else J1 = J-1;
		if(k % 2) K1 = K+1; else K1 = K-1;

		if(J1 == -1 || J1 == my) J1 = J;
		if(K1 == -1 || K1 == mz) K1 = K;

		idx[0] = (PetscInt)ivx[K ][J ][I ];
		idx[1] = (PetscInt)ivx[K ][J1][I ];
		idx[2] = (PetscInt)ivx[K1][J ][I ];
		idx[3] = (PetscInt)ivx[K1][J1][I ];
		idx[4] = (PetscInt)ivx[K ][J ][I1];
		idx[5] = (PetscInt)ivx[K ][J1][I1];
		idx[6] = (PetscInt)ivx[K1][J ][I1];
		idx[7] = (PetscInt)ivx[K1][J1][I1];

		// store full matrix row
		ierr = MatSetValues(P, 1, &row, 8, idx, v, ADD_VALUES); CHKERRQ(ierr);

		// increment row number
		row++;
	}
	END_STD_LOOP

	//---------------------
	// Y-points (fine grid)
	//---------------------
	GET_CELL_RANGE(nx, sx, fine->dsx)
	GET_NODE_RANGE(ny, sy, fine->dsy)
	GET_CELL_RANGE(nz, sz, fine->dsz)

	START_STD_LOOP
	{
		// get coarse grid indices
		I = i/2;
		J = j/2;
		K = k/2;

		if(j % 2) J1 = J+1; else J1 = J;

		if(i % 2) I1 = I+1; else I1 = I-1;
		if(k % 2) K1 = K+1; else K1 = K-1;

		if(I1 == -1 || I1 == mx) I1 = I;
		if(K1 == -1 || K1 == mz) K1 = K;

		idx[0] = (PetscInt)ivy[K ][J ][I ];
		idx[1] = (PetscInt)ivy[K ][J ][I1];
		idx[2] = (PetscInt)ivy[K1][J ][I ];
		idx[3] = (PetscInt)ivy[K1][J ][I1];
		idx[4] = (PetscInt)ivy[K ][J1][I ];
		idx[5] = (PetscInt)ivy[K ][J1][I1];
		idx[6] = (PetscInt)ivy[K1][J1][I ];
		idx[7] = (PetscInt)ivy[K1][J1][I1];

		// store full matrix row
		ierr = MatSetValues(P, 1, &row, 8, idx, v, ADD_VALUES); CHKERRQ(ierr);

		// increment row number
		row++;
	}
	END_STD_LOOP

	//---------------------
	// Z-points (fine grid)
	//---------------------
	GET_CELL_RANGE(nx, sx, fine->dsx)
	GET_CELL_RANGE(ny, sy, fine->dsy)
	GET_NODE_RANGE(nz, sz, fine->dsz)

	START_STD_LOOP
	{
		// get coarse grid indices
		I = i/2;
		J = j/2;
		K = k/2;

		if(k % 2) K1 = K+1; else K1 = K;

		if(i % 2) I1 = I+1; else I1 = I-1;
		if(j % 2) J1 = J+1; else J1 = J-1;

		if(I1 == -1 || I1 == mx) I1 = I;
		if(J1 == -1 || J1 == my) J1 = J;

		idx[0] = (PetscInt)ivz[K ][J ][I ];
		idx[1] = (PetscInt)ivz[K ][J ][I1];
		idx[2] = (PetscInt)ivz[K ][J1][I ];
		idx[3] = (PetscInt)ivz[K ][J1][I1];
		idx[4] = (PetscInt)ivz[K1][J ][I ];
		idx[5] = (PetscInt)ivz[K1][J ][I1];
		idx[6] = (PetscInt)ivz[K1][J1][I ];
		idx[7] = (PetscInt)ivz[K1][J1][I1];

		// store full matrix row
		ierr = MatSetValues(P, 1, &row, 8, idx, v, ADD_VALUES); CHKERRQ(ierr);

		// increment row number
		row++;
	}
	END_STD_LOOP

	if(idxmod == IDXCOUPLED)
	{
		// set pressure weights (direct injection)
		v[0] = 1.0;

		//---------------------
		// P-points (fine grid)
		//---------------------
		GET_CELL_RANGE(nx, sx, fine->dsx)
		GET_CELL_RANGE(ny, sy, fine->dsy)
		GET_CELL_RANGE(nz, sz, fine->dsz)

		START_STD_LOOP
		{
			// get coarse grid indices
			I = i/2;
			J = j/2;
			K = k/2;

			idx[0] = (PetscInt)ip[K][J][I];

			// store full matrix row
			ierr = MatSetValues(P, 1, &row, 1, idx, v, ADD_VALUES); CHKERRQ(ierr);

			// increment row number
			row++;

		}
		END_STD_LOOP

	}
	// restore access
	ierr = DMDAVecRestoreArray(cors->DA_X,   idcors->ivx,  &ivx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(cors->DA_Y,   idcors->ivy,  &ivy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(cors->DA_Z,   idcors->ivz,  &ivz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(cors->DA_CEN, idcors->ip,   &ip);  CHKERRQ(ierr);

	// assemble prolongation matrix
	ierr = PMatAssemble(P, bcfine->numSPC, bcfine->SPCList);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

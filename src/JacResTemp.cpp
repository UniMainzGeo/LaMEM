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
//......................   TEMPERATURE FUNCTIONS   ..........................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "JacRes.h"
#include "phase.h"
#include "scaling.h"
#include "fdstag.h"
#include "tssolve.h"
#include "bc.h"
#include "matrix.h"
#include "surf.h"
#include "dike.h"
#include "heatzone.h"

//---------------------------------------------------------------------------

#define SCATTER_FIELD(da, vec, lT, FIELD)				\
	PetscCall(DMDAGetCorners (da, &sx, &sy, &sz, &nx, &ny, &nz)); \
	PetscCall(DMDAVecGetArray(da, vec, &buff)); \
	iter = 0; \
	START_STD_LOOP \
		FIELD \
	END_STD_LOOP \
	PetscCall(DMDAVecRestoreArray(da, vec, &buff)); \
	LOCAL_TO_LOCAL(da, vec)

#define GET_KC \
  PetscCall(JacResGetTempParam(jr, jr->svCell[iter++].phRat, &kc, NULL, NULL, lT[k][j][i], COORD_CELL(j,sy,fs->dsy), COORD_CELL(i,sx,fs->dsx), COORD_CELL(k,sz,fs->dsz), j-sy, 1.0)); \
  buff[k][j][i] = kc;   // added one NULL because of the new variables that are passed

#define GET_HRXY buff[k][j][i] = jr->svXYEdge[iter++].svDev.Hr;
#define GET_HRXZ buff[k][j][i] = jr->svXZEdge[iter++].svDev.Hr;
#define GET_HRYZ buff[k][j][i] = jr->svYZEdge[iter++].svDev.Hr;

//---------------------------------------------------------------------------
// Temperature parameters functions
//---------------------------------------------------------------------------
PetscErrorCode JacResGetTempParam(
	JacRes      *jr,
	PetscScalar *phRat,
	PetscScalar *k_,      // conductivity
	PetscScalar *rho_Cp_, // volumetric heat capacity
	PetscScalar *rho_A_,  // volumetric radiogenic heat   
	PetscScalar Tc,       // temperature of cell
    PetscScalar y_c,      // center of cell in y-direction
    PetscScalar x_c,      // center of cell in x-direction
    PetscScalar z_c,      // center of cell in z-direction
    PetscInt J,           // coordinate of cell
    PetscScalar sxx_eff_ave_cell) // lithospheric sxx

{
	// compute effective energy parameters in the cell

	PetscInt    i, numPhases, AirPhase;
	Material_t  *phases, *M;
	Controls    ctrl;
	PetscScalar cf, k, rho, rho_Cp, rho_A, density, nu_k, T_Nu; 

	//PetscErrorCode ierr;

	PetscFunctionBeginUser;

	// initialize
	k         = 0.0;
	rho_Cp    = 0.0;
	rho_A     = 0.0;
	nu_k      = 0.0;
	T_Nu	  = 0.0;
	
	numPhases = jr->dbm->numPhases;
	phases    = jr->dbm->phases;
	density   = jr->scal->density;
	AirPhase  = jr->surf->AirPhase;

	// access the control which contains switch for T-dep conductivity
	ctrl      = jr->ctrl;  

	// average all phases
	for(i = 0; i < numPhases; i++)
	{
		M       = &phases[i];
		cf      =  phRat[i];
		rho     =  M->rho;

		// override air phase density
		if(AirPhase != -1 && i == AirPhase)
		{
			rho = 1.0/density;
		}

		k      +=  cf*M->k;
		rho_Cp +=  cf*M->Cp*rho;
		rho_A  +=  cf*M->A*rho;

		// Temperature-dependent conductivity: phase-dependent nusselt number
		if(ctrl.useTk)
		{
		    if(! M->nu_k)
		      {
			// set Nusselt number = 1 if not defined 
			M->nu_k = 1.0;
		      }
		    nu_k +=  cf*M->nu_k;
		    T_Nu +=  cf*M->T_Nu;
		}
		
	}

	// switch and temperature condition to use T-dep conductivity
	if (ctrl.useTk && Tc <= T_Nu) 
	{
	    k = k*nu_k;
	}

	if (ctrl.actDike && ctrl.dikeHeat)
	{
		PetscCall(Dike_k_heatsource(jr, phases, Tc, phRat, k, rho_A, y_c, J, sxx_eff_ave_cell));
	}

	if (ctrl.actHeatZone)
	{
		PetscCall(GetHeatZoneSource(jr, phases, Tc, phRat, rho_A, y_c, x_c, z_c, J, sxx_eff_ave_cell));
	}

	// store
	if(k_)      (*k_)      = k;
	if(rho_Cp_) (*rho_Cp_) = rho_Cp;
	if(rho_A_)  (*rho_A_)  = rho_A;
	
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode JacResCheckTempParam(JacRes *jr)
{
	// check whether thermal material parameters are properly defined

    Material_t  *phases, *M;
	PetscInt    i, numPhases, AirPhase;

	PetscFunctionBeginUser;

	// temperature diffusion cases only
	if(!jr->ctrl.actTemp) PetscFunctionReturn(0);

	// initialize
	numPhases = jr->dbm->numPhases;
	phases    = jr->dbm->phases;
	AirPhase  = jr->surf->AirPhase;

	// check all phases
	for(i = 0; i < numPhases; i++)
	{
		M = &phases[i];

		// check density of the rock phases
		if((AirPhase != -1 && i != AirPhase) || AirPhase == -1)
		{
			if(M->rho == 0.0) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Define density of phase %lld\n", (LLD)i);
		}
			if(M->k   == 0.0) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Define conductivity of phase %lld\n", (LLD)i);
			if(M->Cp  == 0.0) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Define heat capacity of phase %lld\n", (LLD)i);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode JacResCreateTempParam(JacRes *jr)
{
	// setup temperature parameters

	FDSTAG *fs;
	const PetscInt *lx, *ly, *lz;

	PetscFunctionBeginUser;

	fs = jr->fs;

	// create local temperature vector using box-stencil central DMDA
	PetscCall(DMCreateLocalVector(fs->DA_CEN, &jr->lT));

	// temperature diffusion cases only
	if(!jr->ctrl.actTemp) PetscFunctionReturn(0);

	// get cell center grid partitioning
	PetscCall(DMDAGetOwnershipRanges(fs->DA_CEN, &lx, &ly, &lz));

	// create temperature DMDA
	PetscCall(DMDACreate3dSetUp(PETSC_COMM_WORLD,
		DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
		DMDA_STENCIL_STAR,
		fs->dsx.tcels, fs->dsy.tcels, fs->dsz.tcels,
		fs->dsx.nproc, fs->dsy.nproc, fs->dsz.nproc,
		1, 1, lx, ly, lz, &jr->DA_T));

	// create temperature preconditioner matrix
	PetscCall(DMCreateMatrix(jr->DA_T, &jr->Att));

	// set matrix options (development)
	PetscCall(MatSetOption(jr->Att, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));
	PetscCall(MatSetOption(jr->Att, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE));
	PetscCall(MatSetOption(jr->Att, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE));
	PetscCall(MatSetOption(jr->Att, MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_TRUE));

	// temperature solution vector
	PetscCall(DMCreateGlobalVector(jr->DA_T, &jr->dT));

	// energy residual
	PetscCall(DMCreateGlobalVector(jr->DA_T, &jr->ge));

	// create temperature diffusion solver
	PetscCall(KSPCreate(PETSC_COMM_WORLD, &jr->tksp));
	PetscCall(KSPSetOptionsPrefix(jr->tksp,"ts_"));
	PetscCall(KSPSetFromOptions(jr->tksp));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode JacResDestroyTempParam(JacRes *jr)
{
	// destroy temperature parameters

	PetscFunctionBeginUser;

	PetscCall(VecDestroy(&jr->lT));

	// temperature diffusion cases only
	if(!jr->ctrl.actTemp) PetscFunctionReturn(0);

	// temperature parameters
	PetscCall(DMDestroy (&jr->DA_T));
	PetscCall(MatDestroy(&jr->Att));

	PetscCall(VecDestroy(&jr->dT));

	PetscCall(VecDestroy(&jr->ge));

	PetscCall(KSPDestroy(&jr->tksp));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode JacResInitTemp(JacRes *jr)
{
	// initialize temperature from markers

	FDSTAG      *fs;
	BCCtx       *bc;
	PetscScalar ***lT, ***bcT, T;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter;

	PetscFunctionBeginUser;

	// access context
	fs = jr->fs;
	bc = jr->bc;

	PetscCall(VecZeroEntries(jr->lT));

	PetscCall(DMDAVecGetArray(fs->DA_CEN, jr->lT,  &lT));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, bc->bcT, &bcT));

	iter = 0;

	PetscCall(DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		T = bcT[k][j][i];

		if(T == DBL_MAX) T = jr->svCell[iter].svBulk.Tn;

		lT[k][j][i] = T;

		iter++;
	}
	END_STD_LOOP

	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, jr->lT,  &lT));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, bc->bcT, &bcT));

	// apply two-point constraints
	PetscCall(JacResApplyTempBC(jr));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode JacResUpdateTemp(JacRes *jr)
{
	// correct temperature for diffusion (Newton update)

	FDSTAG      *fs;
	PetscScalar ***lT, ***dT;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;

	PetscFunctionBeginUser;

	fs = jr->fs;

	PetscCall(DMDAVecGetArray(fs->DA_CEN, jr->lT, &lT));
	PetscCall(DMDAVecGetArray(jr->DA_T,   jr->dT, &dT));

	PetscCall(DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		lT[k][j][i] -= dT[k][j][i];
	}
	END_STD_LOOP

	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, jr->lT, &lT));
	PetscCall(DMDAVecRestoreArray(jr->DA_T,   jr->dT, &dT));

	// apply two-point constraints
	PetscCall(JacResApplyTempBC(jr));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode JacResApplyTempBC(JacRes *jr)
{
	// apply temperature two-point constraints

	FDSTAG      *fs;
	BCCtx       *bc;
	PetscScalar pmdof;
	PetscScalar ***lT, ***bcT;
	PetscInt    mcx, mcy, mcz;
	PetscInt    I, J, K, fi, fj, fk;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	fs  =  jr->fs;
	bc  =  jr->bc;

	// initialize maximal index in all directions
	mcx = fs->dsx.tcels - 1;
	mcy = fs->dsy.tcels - 1;
	mcz = fs->dsz.tcels - 1;

	// exchange internal ghost points
	LOCAL_TO_LOCAL(fs->DA_CEN, jr->lT)

	// access local solution & boundary constraints
	PetscCall(DMDAVecGetArray(fs->DA_CEN, jr->lT,  &lT));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, bc->bcT, &bcT));

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

        // Set BC @ edges and corners
        if(fi && fj )           SET_EDGE_CORNER(n, lT, k, J, I, k, j, i, pmdof)
        if(fi && fk )           SET_EDGE_CORNER(n, lT, K, j, I, k, j, i, pmdof)
		if(fj && fk )           SET_EDGE_CORNER(n, lT, K, J, i, k, j, i, pmdof)
		if(fi && fj && fk )     SET_EDGE_CORNER(n, lT, K, J, I, k, j, i, pmdof)

        /* 
            Note: a special case occurs for 2D setups, in which nel_y==1
            In that we need to split the setting of edges & corners in two parts 
            to ensure that both front & back side are accounted for.
        */
       	J = j; fj = 0;  if(j == 0)   { fj = 1; J = j-1; }
        if(fi && fj )           SET_EDGE_CORNER(n, lT, k, J, I, k, j, i, pmdof)
        if(fj && fk )           SET_EDGE_CORNER(n, lT, K, J, i, k, j, i, pmdof)
		if(fi && fj && fk   )   SET_EDGE_CORNER(n, lT, K, J, I, k, j, i, pmdof)

	}
	END_STD_LOOP

	// restore access
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, jr->lT,  &lT));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, bc->bcT, &bcT));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode JacResGetTempRes(JacRes *jr, PetscScalar dt)
{
	// compute temperature residual vector
	// STEADY STATE solution is activated by setting time step to zero

	FDSTAG     *fs;
	BCCtx      *bc;
	SolVarCell *svCell;
    SolVarDev  *svDev;
	SolVarBulk *svBulk;
	Dike       *dike; // *djking
	Discret1D  *dsz; // *djking
	Controls    ctrl;
	PetscInt    iter, num, *list;
	PetscInt    Ip1, Im1, Jp1, Jm1, Kp1, Km1;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, my, mz;
	PetscInt    nD, L; // *djking
	PetscScalar ***gsxx_eff_ave, sxx_eff_ave_cell; // *djking
 	PetscScalar bkx, fkx, bky, fky, bkz, fkz;
	PetscScalar bdx, fdx, bdy, fdy, bdz, fdz;
	PetscScalar bqx, fqx, bqy, fqy, bqz, fqz;
	PetscScalar bdpdx, bdpdy, bdpdz, fdpdx, fdpdy, fdpdz;
 	PetscScalar dx, dy, dz;
	PetscScalar invdt, kc, rho_Cp, rho_A, Tc, Pc, Tn, Hr, Ha, cond;
	PetscScalar ***ge, ***lT, ***lk, ***hxy, ***hxz, ***hyz, ***buff, *e,***P;
	PetscScalar ***vx,***vy,***vz;
	PetscScalar y_c, x_c, z_c;

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access residual context variables
	fs    = jr->fs;
	bc    = jr->bc;
	num   = bc->tNumSPC;
	list  = bc->tSPCList;

	// establishing z rank for varM dike heating *djking
	dsz = &fs->dsz;
	L   =  (PetscInt)dsz->rank;

	// access controls
	ctrl = jr->ctrl;

	// initialize maximum cell index in all directions
	mx = fs->dsx.tcels - 1;
	my = fs->dsy.tcels - 1;
	mz = fs->dsz.tcels - 1;

	// compute inverse time step
	if(dt) invdt = 1.0/dt;
	else   invdt = 0.0;

	PetscCall(DMDAVecGetArray(fs->DA_CEN, jr->lT,   &lT));

	SCATTER_FIELD(fs->DA_CEN, jr->ldxx, lT, GET_KC)
	SCATTER_FIELD(fs->DA_XY,  jr->ldxy, lT, GET_HRXY)
	SCATTER_FIELD(fs->DA_XZ,  jr->ldxz, lT, GET_HRXZ)
	SCATTER_FIELD(fs->DA_YZ,  jr->ldyz, lT, GET_HRYZ)

	// access work vectors
	PetscCall(DMDAVecGetArray(jr->DA_T,   jr->ge,   &ge));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, jr->ldxx, &lk));
	PetscCall(DMDAVecGetArray(fs->DA_XY,  jr->ldxy, &hxy));
	PetscCall(DMDAVecGetArray(fs->DA_XZ,  jr->ldxz, &hxz));
	PetscCall(DMDAVecGetArray(fs->DA_YZ,  jr->ldyz, &hyz));
	PetscCall(DMDAVecGetArray(fs->DA_X,   jr->lvx,  &vx) );
	PetscCall(DMDAVecGetArray(fs->DA_Y,   jr->lvy,  &vy) );
	PetscCall(DMDAVecGetArray(fs->DA_Z,   jr->lvz,  &vz) );
	PetscCall(DMDAVecGetArray(fs->DA_CEN, jr->lp_lith, &P ));


	//---------------
	// central points
	//---------------
	iter = 0;
	PetscCall(DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz));

	// *djking
	if (jr->ctrl.actDike)
	{
		nD = 0; // sets dike number to 0 for calculation of sxx_eff_ave across entire domain
		dike = jr->dbdike->matDike + nD;
		PetscCall(DMDAVecGetArray(jr->DA_CELL_2D, dike->sxx_eff_ave, &gsxx_eff_ave)); // *revisit (can we disconnect from individual dike?)
	}

	START_STD_LOOP
	{
		// access solution variables
		svCell = &jr->svCell[iter++];
		svDev  = &svCell->svDev;
		svBulk = &svCell->svBulk;
		
		// access
		Tc  = lT[k][j][i]; // current temperature
		Tn  = svBulk->Tn;  // temperature history
		Pc  = P[k][j][i] ; // Current Pressure

		y_c = COORD_CELL(j, sy, fs->dsy);
		x_c = COORD_CELL(i, sx, fs->dsx);
		z_c = COORD_CELL(k, sz, fs->dsz);

		// conductivity, heat capacity, radiogenic heat production
		if (jr->ctrl.actDike)
		{
			sxx_eff_ave_cell = gsxx_eff_ave[L][j][i];

			PetscCall(JacResGetTempParam(jr, svCell->phRat, &kc, &rho_Cp, &rho_A, Tc, y_c, x_c, z_c, j-sy, sxx_eff_ave_cell));
		}
		else
		{
			PetscCall(JacResGetTempParam(jr, svCell->phRat, &kc, &rho_Cp, &rho_A, Tc, y_c, x_c, z_c, j-sy, 1.0));
		}
		

		// shear heating term (effective)
		Hr = svDev->Hr +
		(hxy[k][j][i] + hxy[k][j+1][i] + hxy[k][j][i+1] + hxy[k][j+1][i+1] +
		 hxz[k][j][i] + hxz[k+1][j][i] + hxz[k][j][i+1] + hxz[k+1][j][i+1] +
		 hyz[k][j][i] + hyz[k+1][j][i] + hyz[k][j+1][i] + hyz[k+1][j+1][i])/4.0;
		Hr = Hr * jr->ctrl.shearHeatEff;

		// check index bounds
		Im1 = i-1; if(Im1 < 0)  Im1++;
		Ip1 = i+1; if(Ip1 > mx) Ip1--;
		Jm1 = j-1; if(Jm1 < 0)  Jm1++;
		Jp1 = j+1; if(Jp1 > my) Jp1--;
		Km1 = k-1; if(Km1 < 0)  Km1++;
		Kp1 = k+1; if(Kp1 > mz) Kp1--;

		// to output as a paraview-field
		cond = kc;
		svBulk->cond = cond;
		
		// compute average conductivities
		bkx = (kc + lk[k][j][Im1])/2.0;      fkx = (kc + lk[k][j][Ip1])/2.0;
		bky = (kc + lk[k][Jm1][i])/2.0;      fky = (kc + lk[k][Jp1][i])/2.0;
		bkz = (kc + lk[Km1][j][i])/2.0;      fkz = (kc + lk[Kp1][j][i])/2.0;

		// get mesh steps
		bdx = SIZE_NODE(i, sx, fs->dsx);     fdx = SIZE_NODE(i+1, sx, fs->dsx);
		bdy = SIZE_NODE(j, sy, fs->dsy);     fdy = SIZE_NODE(j+1, sy, fs->dsy);
		bdz = SIZE_NODE(k, sz, fs->dsz);     fdz = SIZE_NODE(k+1, sz, fs->dsz);

		// compute heat fluxes
		bqx = bkx*(Tc - lT[k][j][i-1])/bdx;   fqx = fkx*(lT[k][j][i+1] - Tc)/fdx;
		bqy = bky*(Tc - lT[k][j-1][i])/bdy;   fqy = fky*(lT[k][j+1][i] - Tc)/fdy;
		bqz = bkz*(Tc - lT[k-1][j][i])/bdz;   fqz = fkz*(lT[k+1][j][i] - Tc)/fdz;

		// Compute the pressure gradient
		if(jr->ctrl.initGuess == 0)
		{
			bdpdx = ((Pc - P[k][j][i-1])/bdx)*vx[k][j][i];        fdpdx = ((P[k][j][i+1] - Pc)/fdx)*vx[k][j][i+1];
			bdpdy = ((Pc - P[k][j-1][i])/bdy)*vy[k][j][i];        fdpdy = ((P[k][j+1][i] - Pc)/fdy)*vy[k][j+1][i];
			bdpdz = ((Pc - P[k-1][j][i])/bdz)*vz[k][j][i];        fdpdz = ((P[k+1][j][i] - Pc)/fdz)*vz[k+1][j][i];

			// Adiabatic Heat term
			Ha = jr->ctrl.AdiabHeat*(Tc*svBulk->alpha*((bdpdx+fdpdx)*0.5+(bdpdy+fdpdy)*0.5+(bdpdz+fdpdz)*0.5));
		}
		else
		{
			Ha = 0.0;
		}

		svBulk->Ha = Ha;

		Ha = jr->ctrl.AdiabHeat*Ha;


		// get mesh steps
		dx = SIZE_CELL(i, sx, fs->dsx);
		dy = SIZE_CELL(j, sy, fs->dsy);
		dz = SIZE_CELL(k, sz, fs->dsz);

		// original balance equation:

		// rho*Cp*(Tc - Tn)/dt = (fqx - bqx)/dx + (fqy - bqy)/dy + (fqz - bqz)/dz + Hr + rho*A

		// to get positive diagonal in the preconditioner matrix
		// put right hand side to the left, which gives the following:

		ge[k][j][i] = rho_Cp*(invdt*(Tc - Tn)) - (fqx - bqx)/dx - (fqy - bqy)/dy - (fqz - bqz)/dz - Hr - rho_A - Ha;
	}
	END_STD_LOOP

	// restore access
	PetscCall(DMDAVecRestoreArray(jr->DA_T,   jr->ge,   &ge));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, jr->lT,   &lT));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, jr->ldxx, &lk));
	PetscCall(DMDAVecRestoreArray(fs->DA_XY,  jr->ldxy, &hxy));
	PetscCall(DMDAVecRestoreArray(fs->DA_XZ,  jr->ldxz, &hxz));
	PetscCall(DMDAVecRestoreArray(fs->DA_YZ,  jr->ldyz, &hyz));
	PetscCall(DMDAVecRestoreArray(fs->DA_X,   jr->lvx,     &vx) );
	PetscCall(DMDAVecRestoreArray(fs->DA_Y,   jr->lvy,     &vy) );
	PetscCall(DMDAVecRestoreArray(fs->DA_Z,   jr->lvz,     &vz) );
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, jr->lp_lith, &P)  );


	// impose primary temperature constraints
	PetscCall(VecGetArray(jr->ge, &e));

	for(i = 0; i < num; i++) e[list[i]] = 0.0;

	PetscCall(VecRestoreArray(jr->ge, &e));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode JacResGetTempMat(JacRes *jr, PetscScalar dt)
{
	// assemble temperature preconditioner matrix
	// STEADY STATE solution is activated by setting time step to zero
	// COMPLETE SINGLE-POINT CONSTRIANT IMLEMENTATION !!!

	FDSTAG     *fs;
	BCCtx      *bc;
	SolVarCell *svCell;
	SolVarBulk *svBulk; // NEW
	PetscInt    iter, num, *list;
	PetscInt    Ip1, Im1, Jp1, Jm1, Kp1, Km1;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, my, mz;
	PetscScalar bkx, fkx, bky, fky, bkz, fkz;
	PetscScalar bdx, fdx, bdy, fdy, bdz, fdz;
 	PetscScalar dx, dy, dz;
	PetscScalar v[7], cf[6], kc, rho_Cp, invdt, Tc, cond;
	MatStencil  row[1], col[7];
	PetscScalar ***lk, ***bcT, ***buff, ***lT;
	PetscScalar y_c, x_c, z_c;
	
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// access residual context variables
	fs   = jr->fs;
	bc   = jr->bc;
	num  = bc->tNumSPC;
	list = bc->tSPCList;

	// compute inverse time step
	if(dt) invdt = 1.0/dt;
	else   invdt = 0.0;

	// initialize maximum cell index in all directions
	mx = fs->dsx.tcels - 1;
	my = fs->dsy.tcels - 1;
	mz = fs->dsz.tcels - 1;

	PetscCall(DMDAVecGetArray(fs->DA_CEN, jr->lT,   &lT));

	SCATTER_FIELD(fs->DA_CEN, jr->ldxx, lT, GET_KC)

	// clear matrix coefficients
	PetscCall(MatZeroEntries(jr->Att));

	// access work vectors
	PetscCall(DMDAVecGetArray(fs->DA_CEN, jr->ldxx, &lk));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, bc->bcT,  &bcT));
	
	//---------------
	// central points
	//---------------
	iter = 0;
	PetscCall(DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		// access solution variables
		svCell = &jr->svCell[iter++];
		svBulk = &svCell->svBulk;

		y_c = COORD_CELL(j, sy, fs->dsy);
		x_c = COORD_CELL(i, sx, fs->dsx);
		z_c = COORD_CELL(k, sz, fs->dsz);
		
		// access
		Tc  = lT[k][j][i]; // current temperature
		
		// conductivity, heat capacity
		PetscCall(JacResGetTempParam(jr, svCell->phRat, &kc, &rho_Cp, NULL, Tc, y_c, x_c, z_c, j-sy, 1.0));

		// check index bounds and TPC multipliers
		Im1 = i-1; cf[0] = 1.0; if(Im1 < 0)  { Im1++; if(bcT[k][j][i-1] != DBL_MAX) cf[0] = -1.0; }
		Ip1 = i+1; cf[1] = 1.0; if(Ip1 > mx) { Ip1--; if(bcT[k][j][i+1] != DBL_MAX) cf[1] = -1.0; }
		Jm1 = j-1; cf[2] = 1.0; if(Jm1 < 0)  { Jm1++; if(bcT[k][j-1][i] != DBL_MAX) cf[2] = -1.0; }
		Jp1 = j+1; cf[3] = 1.0; if(Jp1 > my) { Jp1--; if(bcT[k][j+1][i] != DBL_MAX) cf[3] = -1.0; }
		Km1 = k-1; cf[4] = 1.0; if(Km1 < 0)  { Km1++; if(bcT[k-1][j][i] != DBL_MAX) cf[4] = -1.0; }
		Kp1 = k+1; cf[5] = 1.0; if(Kp1 > mz) { Kp1--; if(bcT[k+1][j][i] != DBL_MAX) cf[5] = -1.0; }

        // to output as a paraview-field
		cond = kc;
		svBulk->cond = cond;
		
 		// compute average conductivities
		bkx = (kc + lk[k][j][Im1])/2.0;      fkx = (kc + lk[k][j][Ip1])/2.0;
		bky = (kc + lk[k][Jm1][i])/2.0;      fky = (kc + lk[k][Jp1][i])/2.0;
		bkz = (kc + lk[Km1][j][i])/2.0;      fkz = (kc + lk[Kp1][j][i])/2.0;

		// get mesh steps
		bdx = SIZE_NODE(i, sx, fs->dsx);     fdx = SIZE_NODE(i+1, sx, fs->dsx);
		bdy = SIZE_NODE(j, sy, fs->dsy);     fdy = SIZE_NODE(j+1, sy, fs->dsy);
		bdz = SIZE_NODE(k, sz, fs->dsz);     fdz = SIZE_NODE(k+1, sz, fs->dsz);

		// get mesh steps
		dx = SIZE_CELL(i, sx, fs->dsx);
		dy = SIZE_CELL(j, sy, fs->dsy);
		dz = SIZE_CELL(k, sz, fs->dsz);

		// set row/column indices
		row[0].k = k;   row[0].j = j;   row[0].i = i;   row[0].c = 0;
		col[0].k = k;   col[0].j = j;   col[0].i = Im1; col[0].c = 0;
		col[1].k = k;   col[1].j = j;   col[1].i = Ip1; col[1].c = 0;
		col[2].k = k;   col[2].j = Jm1; col[2].i = i;   col[2].c = 0;
		col[3].k = k;   col[3].j = Jp1; col[3].i = i;   col[3].c = 0;
		col[4].k = Km1; col[4].j = j;   col[4].i = i;   col[4].c = 0;
		col[5].k = Kp1; col[5].j = j;   col[5].i = i;   col[5].c = 0;
		col[6].k = k;   col[6].j = j;   col[6].i = i;   col[6].c = 0;

		// set values including TPC multipliers
		v[0] = -bkx/bdx/dx*cf[0];
		v[1] = -fkx/fdx/dx*cf[1];
		v[2] = -bky/bdy/dy*cf[2];
		v[3] = -fky/fdy/dy*cf[3];
		v[4] = -bkz/bdz/dz*cf[4];
		v[5] = -fkz/fdz/dz*cf[5];
		v[6] =  invdt*rho_Cp
		+       (bkx/bdx + fkx/fdx)/dx
		+       (bky/bdy + fky/fdy)/dy
		+       (bkz/bdz + fkz/fdz)/dz;

		// set matrix coefficients
		PetscCall(MatSetValuesStencil(jr->Att, 1, row, 7, col, v, ADD_VALUES));

		// NOTE! since only TPC are active, no SPC modification is necessary
	}
	END_STD_LOOP

	// restore access
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, jr->ldxx, &lk));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, bc->bcT, &bcT));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, jr->lT,   &lT));

	// assemble temperature matrix
	PetscCall(MatAIJAssemble(jr->Att, num, list, 1.0));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
/*
Diffusion term expansion

		bqx = bkx*(Tc - T[k][j][i-1])/bdx;   fqx = fkx*(T[k][j][i+1] - Tc)/fdx;
		bqy = bky*(Tc - T[k][j-1][i])/bdy;   fqy = fky*(T[k][j+1][i] - Tc)/fdy;
		bqz = bkz*(Tc - T[k-1][j][i])/bdz;   fqz = fkz*(T[k+1][j][i] - Tc)/fdz;

[1]
		-(fqx - bqx)/dx
		-(fqy - bqy)/dy
		-(fqz - bqz)/dz
[2]
		(bqx - fqx)/dx
		(bqy - fqy)/dy
		(bqz - fqz)/dz
[3]
		(bkx*(Tc - T[k][j][i-1])/bdx - fkx*(T[k][j][i+1] - Tc)/fdx)/dx
		(bky*(Tc - T[k][j-1][i])/bdy - fky*(T[k][j+1][i] - Tc)/fdy)/dy
		(bkz*(Tc - T[k-1][j][i])/bdz - fkz*(T[k+1][j][i] - Tc)/fdz)/dz
[4]
		(bkx*(Tc - T[k][j][i-1])/bdx + fkx*(Tc - T[k][j][i+1])/fdx)/dx
		(bky*(Tc - T[k][j-1][i])/bdy + fky*(Tc - T[k][j+1][i])/fdy)/dy
		(bkz*(Tc - T[k-1][j][i])/bdz + fkz*(Tc - T[k+1][j][i])/fdz)/dz
[5]
		bkx/bdx/dx*(Tc - T[k][j][i-1]) + fkx/fdx/dx*(Tc - T[k][j][i+1])
		bky/bdy/dy*(Tc - T[k][j-1][i]) + fky/fdy/dy*(Tc - T[k][j+1][i])
		bkz/bdz/dz*(Tc - T[k-1][j][i]) + fkz/fdz/dz*(Tc - T[k+1][j][i])
[6]
		(bkx/bdx/dx + fkx/fdx/dx)*Tc - bkx/bdx/dx*T[k][j][i-1] - fkx/fdx/dx*T[k][j][i+1]
		(bky/bdy/dy + fky/fdy/dy)*Tc - bky/bdy/dy*T[k][j-1][i] - fky/fdy/dy*T[k][j+1][i]
		(bkz/bdz/dz + fkz/fdz/dz)*Tc - bkz/bdz/dz*T[k-1][j][i] - fkz/fdz/dz*T[k+1][j][i]
[7]

		(bkx/bdx + fkx/fdx)/dx*Tc - bkx/bdx/dx*T[k][j][i-1] - fkx/fdx/dx*T[k][j][i+1]
		(bky/bdy + fky/fdy)/dy*Tc - bky/bdy/dy*T[k][j-1][i] - fky/fdy/dy*T[k][j+1][i]
		(bkz/bdz + fkz/fdz)/dz*Tc - bkz/bdz/dz*T[k-1][j][i] - fkz/fdz/dz*T[k+1][j][i]
*/
//---------------------------------------------------------------------------

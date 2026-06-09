/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : LaMEM
 **   License      : MIT, see LICENSE file for details
 **   Contributors : Anton Popov, Boris Kaus, see AUTHORS file for complete list
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : kaus@uni-mainz.de, popov@uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
#include "LaMEM.h"
#include "scaling.h"
#include "bc.h"
#include "JacRes.h"
#include "fdstag.h"
#include "surf.h"
#include "phase.h"
#include "tools.h"
#include "Tensor.h"
#include "parsing.h"
//---------------------------------------------------------------------------
PetscErrorCode JacResGetSHmax(JacRes *jr)
{
	// compute maximum horizontal compressive stress (SHmax) orientation

	FDSTAG      *fs;
	SolVarCell  *svCell;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter;
	PetscScalar v1[3], v2[3], sxx, syy, sxy, s1, s2;
	PetscScalar ***dx, ***dy, ***lsxy;

	
	PetscFunctionBeginUser;

	// access context
	fs = jr->fs;

	// setup shear stress vector
	PetscCall(DMDAVecGetArray(fs->DA_XY, jr->ldxy, &lsxy));

	PetscCall(DMDAGetCorners(fs->DA_XY, &sx, &sy, &sz, &nx, &ny, &nz));

	iter = 0;

	START_STD_LOOP
	{
		lsxy[k][j][i] = jr->svXYEdge[iter++].s;
	}
	END_STD_LOOP

	PetscCall(DMDAVecRestoreArray(fs->DA_XY, jr->ldxy, &lsxy));

	LOCAL_TO_LOCAL(fs->DA_XY, jr->ldxy);

	// get SHmax
	PetscCall(DMDAVecGetArray(fs->DA_CEN, jr->ldxx, &dx));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, jr->ldyy, &dy));
	PetscCall(DMDAVecGetArray(fs->DA_XY,  jr->ldxy, &lsxy));

	PetscCall(DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz));

	iter = 0;

	START_STD_LOOP
	{
		svCell = &jr->svCell[iter++];
		sxx    = svCell->sxx;
		syy    = svCell->syy;
		sxy    = (lsxy[k][j][i] + lsxy[k][j][i+1] + lsxy[k][j+1][i] + lsxy[k][j+1][i+1])/4.0;

		// maximum compressive stress orientation is the eigenvector of the SMALLEST eigenvalue
		// (stress is negative in compression)
		PetscCall(Tensor2RS2DSpectral(sxx, syy, sxy, &s1, &s2, v1, v2, 1e-12));

		// get common sense
		if(v2[0] < 0.0 || (v2[0] == 0.0 && v2[1] < 0.0))
		{
			v2[0] = -v2[0];
			v2[1] = -v2[1];
		}

		// store direction vector for output
		dx[k][j][i] = v2[0];
		dy[k][j][i] = v2[1];
	}
	END_STD_LOOP

	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, jr->ldxx, &dx));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, jr->ldyy, &dy));
	PetscCall(DMDAVecRestoreArray(fs->DA_XY,  jr->ldxy, &lsxy));

	LOCAL_TO_LOCAL(fs->DA_CEN, jr->ldxx);
	LOCAL_TO_LOCAL(fs->DA_CEN, jr->ldyy);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode JacResGetEHmax(JacRes *jr)
{

	// compute maximum horizontal extension rate (EHmax) orientation

	FDSTAG      *fs;
	SolVarCell  *svCell;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter;
	PetscScalar v1[3], v2[3], dxx, dyy, dxy, d1, d2;
	PetscScalar ***dx, ***dy, ***ldxy;

	
	PetscFunctionBeginUser;

	// access context
	fs = jr->fs;

	// setup shear strain rate vector
	PetscCall(DMDAVecGetArray(fs->DA_XY, jr->ldxy, &ldxy));

	PetscCall(DMDAGetCorners(fs->DA_XY, &sx, &sy, &sz, &nx, &ny, &nz));

	iter = 0;

	START_STD_LOOP
	{
		ldxy[k][j][i] = jr->svXYEdge[iter++].d;
	}
	END_STD_LOOP

	PetscCall(DMDAVecRestoreArray(fs->DA_XY, jr->ldxy, &ldxy));

	LOCAL_TO_LOCAL(fs->DA_XY, jr->ldxy);

	// get EHmax
	PetscCall(DMDAVecGetArray(fs->DA_CEN, jr->ldxx, &dx));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, jr->ldyy, &dy));
	PetscCall(DMDAVecGetArray(fs->DA_XY,  jr->ldxy, &ldxy));

	PetscCall(DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz));

	iter = 0;

	START_STD_LOOP
	{
		svCell = &jr->svCell[iter++];
		dxx    = svCell->dxx;
		dyy    = svCell->dyy;
		dxy    = (ldxy[k][j][i] + ldxy[k][j][i+1] + ldxy[k][j+1][i] + ldxy[k][j+1][i+1])/4.0;

		// maximum extension rate orientation is the eigenvector of the LARGEST eigenvalue
		// (strain rate is positive in extension)
		PetscCall(Tensor2RS2DSpectral(dxx, dyy, dxy, &d1, &d2, v1, v2, 1e-12));

		// get common sense
		if(v1[0] < 0.0 || (v1[0] == 0.0 && v1[1] < 0.0))
		{
			v1[0] = -v1[0];
			v1[1] = -v1[1];
		}

		// store direction vector for output
		dx[k][j][i] = v1[0];
		dy[k][j][i] = v1[1];
	}
	END_STD_LOOP

	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, jr->ldxx, &dx));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, jr->ldyy, &dy));
	PetscCall(DMDAVecRestoreArray(fs->DA_XY,  jr->ldxy, &ldxy));

	LOCAL_TO_LOCAL(fs->DA_CEN, jr->ldxx);
	LOCAL_TO_LOCAL(fs->DA_CEN, jr->ldyy);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode JacResGetOverPressure(JacRes *jr, Vec lop)
{
	// compute overpressure

	FDSTAG      *fs;
	PetscScalar ***op, ***p, ***p_lith;
	PetscInt    i, j, k, sx, sy, sz, nx, ny, nz;

	
	PetscFunctionBeginUser;

	// access context
	fs  =  jr->fs;

	// get local grid sizes
	PetscCall(DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz));

	// access pressure vectors
	PetscCall(VecZeroEntries(lop));

	PetscCall(DMDAVecGetArray(fs->DA_CEN, lop,         &op));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, jr->lp,      &p));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, jr->lp_lith, &p_lith));

	START_STD_LOOP
	{
		// overpressure = dynamic - lithostatic
		op[k][j][i] = p[k][j][i] - p_lith[k][j][i];
	}
	END_STD_LOOP

	// restore buffer and pressure vectors
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, lop,         &op));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, jr->lp,      &p));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, jr->lp_lith, &p_lith));

	// fill ghost points
	LOCAL_TO_LOCAL(fs->DA_CEN, lop)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode JacResGetLithoStaticPressure(JacRes *jr)
{
	// compute lithostatic pressure

	Vec         vbuff;
	FDSTAG      *fs;
	Discret1D   *dsz;
	MPI_Request srequest, rrequest;
	PetscScalar ***lp, ***ibuff, *lbuff, dz, dp, g, rho;
	PetscInt    i, j, k, sx, sy, sz, nx, ny, nz, iter, L;

	
	PetscFunctionBeginUser;

	// access context
	fs  =  jr->fs;
	dsz = &fs->dsz;
	L   =  (PetscInt)dsz->rank;
	g   =   PetscAbsScalar(jr->ctrl.grav[2]);

	// initialize
	PetscCall(VecZeroEntries(jr->lp_lith));

	// get local grid sizes
	PetscCall(DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz));

	// get integration/communication buffer
	PetscCall(DMGetGlobalVector(jr->DA_CELL_2D, &vbuff));

	PetscCall(VecZeroEntries(vbuff));

	// open index buffer for computation
	PetscCall(DMDAVecGetArray(jr->DA_CELL_2D, vbuff, &ibuff));

	// open linear buffer for send/receive
	PetscCall(VecGetArray(vbuff, &lbuff));

	// access lithostatic pressure
	PetscCall(DMDAVecGetArray(fs->DA_CEN, jr->lp_lith, &lp));

	// start receiving integral from top domain (next)
	if(dsz->nproc != 1 && dsz->grnext != -1)
	{
		PetscCallMPI(MPI_Irecv(lbuff, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grnext, 0, PETSC_COMM_WORLD, &rrequest));
	}

	// copy density
	iter = 0;

	START_STD_LOOP
	{
		lp[k][j][i] = jr->svCell[iter++].svBulk.rho;
	}
	END_STD_LOOP

	// finish receiving
	if(dsz->nproc != 1 && dsz->grnext != -1)
	{
		PetscCallMPI(MPI_Wait(&rrequest, MPI_STATUSES_IGNORE));
	}

	// compute local integral from top to bottom
	for(k = sz + nz - 1; k >= sz; k--)
	{
		START_PLANE_LOOP
		{
			// get density, cell size, pressure increment
			rho = lp[k][j][i];
			dz  = SIZE_CELL(k, sz, (*dsz));
			dp  = rho*g*dz;

			// store  lithostatic pressure
			lp[k][j][i] = ibuff[L][j][i] + dp/2.0;

			// update lithostatic pressure integral
			ibuff[L][j][i] += dp;
		}
		END_PLANE_LOOP
	}

	// send integral to bottom domain (previous)
	if(dsz->nproc != 1 && dsz->grprev != -1)
	{
		PetscCallMPI(MPI_Isend(lbuff, (PetscMPIInt)(nx*ny), MPIU_SCALAR, dsz->grprev, 0, PETSC_COMM_WORLD, &srequest));

		PetscCallMPI(MPI_Wait(&srequest, MPI_STATUSES_IGNORE));
	}

	// restore buffer and pressure vectors
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, jr->lp_lith, &lp));

	PetscCall(DMDAVecRestoreArray(jr->DA_CELL_2D, vbuff, &ibuff));

	PetscCall(VecRestoreArray(vbuff, &lbuff));

	PetscCall(DMRestoreGlobalVector(jr->DA_CELL_2D, &vbuff));

	// fill ghost points
	LOCAL_TO_LOCAL(fs->DA_CEN, jr->lp_lith)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode JacResGetPorePressure(JacRes *jr)
{
	// compute pore pressure
	FDSTAG      *fs;
	Controls    *ctrl;
	Material_t  *phases, *mat;
	PetscScalar ***lp_pore, ***lp_lith, *phRat;
	PetscScalar ztop, g, gwLevel=0.0, rho_fluid, depth, p_hydro, rp_cv, rp;
	PetscInt    numPhases, i, j, k, iter, iphase, sx, sy, sz, nx, ny, nz;

	
	PetscFunctionBeginUser;

	// initialize
	PetscCall(VecZeroEntries(jr->lp_pore));

	// return if not activated
	if(jr->ctrl.gwType == _GW_NONE_) PetscFunctionReturn(0);

	// access context
	fs        =  jr->fs;
	phases    =  jr->dbm->phases;
	numPhases =  jr->dbm->numPhases;
	ctrl      = &jr->ctrl;
	rho_fluid =  ctrl->rho_fluid;
	g         =  PetscAbsScalar(ctrl->grav[2]);

	// get top boundary coordinate
	PetscCall(FDSTAGGetGlobalBox(fs, NULL, NULL, NULL, NULL, NULL, &ztop));

	// set ground water level
	if     (ctrl->gwType == _GW_TOP_)   gwLevel = ztop;
	else if(ctrl->gwType == _GW_SURF_)  gwLevel = jr->surf->avg_topo;
	else if(ctrl->gwType == _GW_LEVEL_) gwLevel = ctrl->gwLevel;

	// get local grid sizes
	PetscCall(DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz));

	// access vectors
	PetscCall(DMDAVecGetArray(fs->DA_CEN, jr->lp_pore, &lp_pore));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, jr->lp_lith, &lp_lith));

	iter = 0;
	START_STD_LOOP
	{
		// access phase ratio array
		phRat = jr->svCell[iter++].phRat;

		// compute depth of the current control volume
		depth = gwLevel - COORD_CELL(k, sz, fs->dsz);
		if(depth < 0.0) depth = 0.0;				// we don't want these calculations in the 'air'

		// Evaluate pore pressure ratio in control volume
		rp_cv = 0.0;
		// scan all phases
		for(iphase = 0; iphase < numPhases; iphase++)
		{
			// update present phases only
			if(phRat[iphase])
			{
				// get reference to material parameters table
				mat = &phases[iphase];

				// get and check pore pressure ratio of each phase
				if(mat->rp<0.0)      mat->rp = 0.0;
				else if(mat->rp>1.0) mat->rp = 1.0;
				rp = mat->rp;

				// compute average pore pressure ratio
				rp_cv +=  phRat[iphase] * rp;
			}
		}

		// hydrostatic pressure (based on the water column)
		p_hydro = rho_fluid * g * PetscAbsScalar(depth);

		// compute the pore pressure as product of lithostatic pressure and porepressure ratio of the control volume
		lp_pore[k][j][i] =  p_hydro + rp_cv * (lp_lith[k][j][i]-p_hydro);
	}
	END_STD_LOOP

	// restore buffer and pressure vectors
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, jr->lp_pore, &lp_pore));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, jr->lp_lith, &lp_lith));

	// fill ghost points
	LOCAL_TO_LOCAL(fs->DA_CEN, jr->lp_pore)

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode JacResGetPermea(JacRes *jr, PetscInt bgPhase, PetscInt step, char *outfile)
{
	FILE        *db;
	FDSTAG      *fs;
	BCCtx       *bc;
	Material_t  *phases;
	Scaling     *scal;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***vz, nZFace, lvel, gvel, dp, eta, ks, bz, ez;
	char        path[_str_len_];

	
	PetscFunctionBeginUser;

	// check activation
	if(!jr->ctrl.getPermea || !step) PetscFunctionReturn(0);

	// access context variables
	fs     = jr->fs;
	bc     = jr->bc;
	phases = jr->dbm->phases;
	scal   = jr->scal;

	// get grid bounds in z-direction
	PetscCall(FDSTAGGetGlobalBox(fs, NULL, NULL, &bz, NULL, NULL, &ez));

	// get total number of z-faces
	nZFace = (PetscScalar)(fs->dsx.tcels*fs->dsy.tcels*fs->dsz.tnods);

	// get fluid viscosity (background phase is fluid)
	eta = 1.0/(2.0*phases[bgPhase].Bd);

	// get pressure gradient
	dp = (bc->pbot - bc->ptop)/(ez-bz);

	// get local grid sizes
	PetscCall(DMDAGetCorners(fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz));

	// access z-velocity vector
	PetscCall(DMDAVecGetArray(fs->DA_Z, jr->lvz, &vz));

	// compute volume average absolute velocity
	lvel = 0.0;

	//---------
	// Z-points
	//---------

	START_STD_LOOP
	{
		lvel += PetscAbsScalar(vz[k][j][i]);
	}
	END_STD_LOOP

	// restore access
	PetscCall(DMDAVecRestoreArray(fs->DA_Z, jr->lvz, &vz));

	// compute global sum
	if(ISParallel(PETSC_COMM_WORLD))
	{
		PetscCallMPI(MPI_Allreduce(&lvel, &gvel, 1, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD));
	}
	else
	{
		gvel = lvel;
	}

	// normalize
	gvel /= nZFace;

	// compute permeability
	ks = gvel*eta/dp;

	ks = PetscAbsScalar(ks);

	// output to the file
	if(ISRankZero(PETSC_COMM_WORLD))
	{

		memset(path, 0, _str_len_);
		strcpy(path, outfile);
		strcat(path, ".darcy.dat");

		db = fopen(path, "wb");

		fprintf(db,"# ==============================================\n");
		fprintf(db,"# EFFECTIVE PERMEABILITY CONSTANT: %E %s \n ", ks*scal->area_si, scal->lbl_area_si);
		fprintf(db,"# ==============================================\n");

		fclose(db);

		// Also output this to the screen (for regression testing, for example)
		PetscPrintf(PETSC_COMM_WORLD,"\n");
		PetscPrintf(PETSC_COMM_WORLD,"==========================================================================\n");
		PetscPrintf(PETSC_COMM_WORLD,"EFFECTIVE PERMEABILITY CONSTANT: %E %s\n", ks*scal->area_si, scal->lbl_area_si);
		PetscPrintf(PETSC_COMM_WORLD,"==========================================================================\n");
		
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode JacResGetVelGrad(JacRes *jr,
		Vec dvxdx, Vec dvxdy, Vec dvxdz,
		Vec dvydx, Vec dvydy, Vec dvydz,
		Vec dvzdx, Vec dvzdy, Vec dvzdz)
{
	// compute velocity gradients for output

	FDSTAG     *fs;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar dx, dy, dz;
	PetscScalar ***vx,  ***vy,  ***vz;
	PetscScalar ***dxx, ***dxy, ***dxz;
	PetscScalar ***dyx, ***dyy, ***dyz;
	PetscScalar ***dzx, ***dzy, ***dzz;

	
	PetscFunctionBeginUser;

	fs = jr->fs;

	// access velocity and gradients
	PetscCall(DMDAVecGetArray(fs->DA_X, jr->lvx, &vx));
	PetscCall(DMDAVecGetArray(fs->DA_Y, jr->lvy, &vy));
	PetscCall(DMDAVecGetArray(fs->DA_Z, jr->lvz, &vz));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, dvxdx, &dxx));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, dvydy, &dyy));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, dvzdz, &dzz));
	PetscCall(DMDAVecGetArray(fs->DA_XY,  dvxdy, &dxy));
	PetscCall(DMDAVecGetArray(fs->DA_XY,  dvydx, &dyx));
	PetscCall(DMDAVecGetArray(fs->DA_XZ,  dvxdz, &dxz));
	PetscCall(DMDAVecGetArray(fs->DA_XZ,  dvzdx, &dzx));
	PetscCall(DMDAVecGetArray(fs->DA_YZ,  dvydz, &dyz));
	PetscCall(DMDAVecGetArray(fs->DA_YZ,  dvzdy, &dzy));

	//-------------------------------
	// central points (dxx, dyy, dzz)
	//-------------------------------
	PetscCall(DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		// get mesh steps
		dx = SIZE_CELL(i, sx, fs->dsx);
		dy = SIZE_CELL(j, sy, fs->dsy);
		dz = SIZE_CELL(k, sz, fs->dsz);

		// compute velocity gradients
		dxx[k][j][i] = (vx[k][j][i+1] - vx[k][j][i])/dx;
		dyy[k][j][i] = (vy[k][j+1][i] - vy[k][j][i])/dy;
		dzz[k][j][i] = (vz[k+1][j][i] - vz[k][j][i])/dz;
	}
	END_STD_LOOP

	//-------------------------------
	// xy edge points (dxy)
	//-------------------------------
	PetscCall(DMDAGetCorners(fs->DA_XY, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		// get mesh steps
		dx = SIZE_NODE(i, sx, fs->dsx);
		dy = SIZE_NODE(j, sy, fs->dsy);

		// compute velocity gradients
		dxy[k][j][i] = (vx[k][j][i] - vx[k][j-1][i])/dy;
		dyx[k][j][i] = (vy[k][j][i] - vy[k][j][i-1])/dx;
	}
	END_STD_LOOP

	//-------------------------------
	// xz edge points (dxz)
	//-------------------------------
	PetscCall(DMDAGetCorners(fs->DA_XZ, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		// get mesh steps
		dx = SIZE_NODE(i, sx, fs->dsx);
		dz = SIZE_NODE(k, sz, fs->dsz);

		// compute velocity gradients
		dxz[k][j][i] = (vx[k][j][i] - vx[k-1][j][i])/dz;
		dzx[k][j][i] = (vz[k][j][i] - vz[k][j][i-1])/dx;
	}
	END_STD_LOOP

	//-------------------------------
	// yz edge points (dyz)
	//-------------------------------
	PetscCall(DMDAGetCorners(fs->DA_YZ, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		// get mesh steps
		dy = SIZE_NODE(j, sy, fs->dsy);
		dz = SIZE_NODE(k, sz, fs->dsz);

		// compute velocity gradients
		dyz[k][j][i] = (vy[k][j][i] - vy[k-1][j][i])/dz;
		dzy[k][j][i] = (vz[k][j][i] - vz[k][j-1][i])/dy;
	}
	END_STD_LOOP

	// restore access
	PetscCall(DMDAVecRestoreArray(fs->DA_X, jr->lvx, &vx));
	PetscCall(DMDAVecRestoreArray(fs->DA_Y, jr->lvy, &vy));
	PetscCall(DMDAVecRestoreArray(fs->DA_Z, jr->lvz, &vz));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, dvxdx, &dxx));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, dvydy, &dyy));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, dvzdz, &dzz));
	PetscCall(DMDAVecRestoreArray(fs->DA_XY,  dvxdy, &dxy));
	PetscCall(DMDAVecRestoreArray(fs->DA_XY,  dvydx, &dyx));
	PetscCall(DMDAVecRestoreArray(fs->DA_XZ,  dvxdz, &dxz));
	PetscCall(DMDAVecRestoreArray(fs->DA_XZ,  dvzdx, &dzx));
	PetscCall(DMDAVecRestoreArray(fs->DA_YZ,  dvydz, &dyz));
	PetscCall(DMDAVecRestoreArray(fs->DA_YZ,  dvzdy, &dzy));

	// communicate boundary strain-rate values
	LOCAL_TO_LOCAL(fs->DA_CEN, dvxdx);
	LOCAL_TO_LOCAL(fs->DA_CEN, dvydy);
	LOCAL_TO_LOCAL(fs->DA_CEN, dvzdz);
	LOCAL_TO_LOCAL(fs->DA_XY,  dvxdy);
	LOCAL_TO_LOCAL(fs->DA_XY,  dvydx);
	LOCAL_TO_LOCAL(fs->DA_XZ,  dvxdz);
	LOCAL_TO_LOCAL(fs->DA_XZ,  dvzdx);
	LOCAL_TO_LOCAL(fs->DA_YZ,  dvydz);
	LOCAL_TO_LOCAL(fs->DA_YZ,  dvzdy);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

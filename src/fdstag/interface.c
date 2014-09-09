//---------------------------------------------------------------------------
//..............   LaMEM - FDSTAG CANONICAL INTERFACE ROUTINES   ............
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "bc.h"
#include "JacRes.h"
#include "interface.h"
#include "Utils.h"
//---------------------------------------------------------------------------
/*
#undef __FUNCT__
#define __FUNCT__ "FDSTAGSetCoordDMDA"
PetscErrorCode FDSTAGSetCoordDMDA(FDSTAG *fs, DM da)
{
	// copy local coordinates from DMDA object

	DM             cda;
	Vec            gc;
	DMDACoor3d     ***coors;
	PetscScalar    *crdx, *crdy, *crdz;
	PetscInt       i, j, k, m, n, p, sm, sn, sp;
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get grid sizes
	GET_NODE_RANGE(m, sm, fs->dsx);
	GET_NODE_RANGE(n, sn, fs->dsy);
	GET_NODE_RANGE(p, sp, fs->dsz);

	// open access to local coordinate vector
	ierr = DMGetCoordinates(da, &gc);      CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(da, &cda);    CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda, gc, &coors); CHKERRQ(ierr);

	// create local copies
	ierr = makeScalArray(&crdx, 0, m); CHKERRQ(ierr);
	ierr = makeScalArray(&crdy, 0, n); CHKERRQ(ierr);
	ierr = makeScalArray(&crdz, 0, p); CHKERRQ(ierr);

	// copy local coordinates of the nodes to 1D vectors
	for(i = sm; i < sm+m; i++) crdx[i-sm] = coors[sp][sn][i ].x;
	for(j = sn; j < sn+n; j++) crdy[j-sn] = coors[sp][j ][sm].y;
	for(k = sp; k < sp+p; k++) crdz[k-sp] = coors[k ][sn][sm].z;

	// close access
	ierr = DMDAVecRestoreArray(cda, gc, &coors);  CHKERRQ(ierr);

	// setup all necessary local coordinates for FDSTAG
//	ierr = Discret1DSetCoord(&fs->dsx, crdx); CHKERRQ(ierr);
//	ierr = Discret1DSetCoord(&fs->dsy, crdy); CHKERRQ(ierr);
//	ierr = Discret1DSetCoord(&fs->dsz, crdz); CHKERRQ(ierr);

	// free arrays
	ierr = PetscFree(crdx); CHKERRQ(ierr);
	ierr = PetscFree(crdy); CHKERRQ(ierr);
	ierr = PetscFree(crdz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
*/
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAGInitMaterialProps"
PetscErrorCode FDSTAGInitMaterialProps(JacResCtx *jrctx, UserContext *usr)
{
	// initialize material properties in the FDSTAG data structures

	PetscInt    i, viscLaw;
	PetscScalar eta, D, n;
	Material_t *phases;
	PhaseProps *PhaseProperties;
	MatParLim  *matLim;

	PetscFunctionBegin;

	// access phase properties in user context variables
	phases          = jrctx->phases;
	PhaseProperties = &usr->PhaseProperties;
	matLim          = &jrctx->matLim;

	// clear phase properties array
	memset(phases, 0, sizeof(Material_t)*(size_t)jrctx->numPhases);

	// copy phase parameters
	for(i = 0; i < jrctx->numPhases; i++)
	{
		// density and power law parameters
		phases[i].ID     = i;
		phases[i].rho    = PhaseProperties->rho[i];
		phases[i].frSoft = NULL;
		phases[i].chSoft = NULL;

		// get viscosity law
		viscLaw = PhaseProperties->ViscosityLaw[i];

		if(viscLaw == 1)
		{
			// set diffusion creep constant
			phases[i].Bd = 1.0/(2.0*PhaseProperties->mu[i]);
		}
		else if(viscLaw == 2)
		{
			// power-law creep parameters
			eta = PhaseProperties->mu[i];
			D   = PhaseProperties->Powerlaw_e0[i];
			n   = PhaseProperties->n_exponent[i];

			// convert power-law creep parameters to dislocation creep parameters
			phases[i].Bn = pow(2.0*eta, -n)*pow(D, 1-n);
			phases[i].n  = n;

			// store reference strain rate
			matLim->DII_ref = D;
		}
		else
		{	// unsupported viscosity law
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "ERROR! Unsupported viscosity law used\n");
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
/*
#undef __FUNCT__
#define __FUNCT__ "FDSTAGInitPhaseRatios"
PetscErrorCode FDSTAGInitPhaseRatios(FDSTAG *fs, JacResCtx *jrctx, UserContext *usr)
{
	// initialize phase ratios in the FDSTAG data structures

	PetscScalar *phRat;
	FDSTAG_MaterialProperties *fdstag;
	PetscScalar ***Center_PhaseProportions  [max_num_phases];
	PetscScalar ***XYPoints_PhaseProportions[max_num_phases];
	PetscScalar ***XZPoints_PhaseProportions[max_num_phases];
	PetscScalar ***YZPoints_PhaseProportions[max_num_phases];
	PetscInt i, j, k, nx, ny, nz, sx, sy, sz, iter, jj, numPhases;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access solution variables
	fdstag    = &usr->FDSTAG;
	numPhases = jrctx->numPhases;

	// gain access to LaMEM storage scheme
	for(jj = 0; jj < numPhases; jj++)
	{
		ierr = DMDAVecGetArray(fdstag->DA_CENTER,    fdstag->Center_PhaseProportions  [jj], &Center_PhaseProportions  [jj]); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(fdstag->DA_XY_POINTS, fdstag->XYPoints_PhaseProportions[jj], &XYPoints_PhaseProportions[jj]); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(fdstag->DA_XZ_POINTS, fdstag->XZPoints_PhaseProportions[jj], &XZPoints_PhaseProportions[jj]); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(fdstag->DA_YZ_POINTS, fdstag->YZPoints_PhaseProportions[jj], &YZPoints_PhaseProportions[jj]); CHKERRQ(ierr);
	}

	//---------------
	// central points
	//---------------
	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		phRat = jrctx->svCell[iter++].svDev.phRat;

		for(jj = 0; jj < numPhases; jj++)
			phRat[jj] = Center_PhaseProportions[jj][k][j][i];
	}
	END_STD_LOOP

	//---------------
	// xy edge points
	//---------------
	iter = 0;
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		phRat = jrctx->svXYEdge[iter++].svDev.phRat;

		for(jj = 0; jj < numPhases; jj++)
			phRat[jj] = XYPoints_PhaseProportions[jj][k][j][i];
	}
	END_STD_LOOP

	//---------------
	// xz edge points
	//---------------
	iter = 0;
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		phRat = jrctx->svXZEdge[iter++].svDev.phRat;

		for(jj = 0; jj < numPhases; jj++)
			phRat[jj] = XZPoints_PhaseProportions[jj][k][j][i];
	}
	END_STD_LOOP

	//---------------
	// yz edge points
	//---------------
	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		phRat = jrctx->svYZEdge[iter++].svDev.phRat;

		for(jj = 0; jj < numPhases; jj++)
			phRat[jj] = YZPoints_PhaseProportions[jj][k][j][i];
	}
	END_STD_LOOP

	// restore access
	for(jj = 0; jj < jrctx->numPhases; jj++)
	{
		ierr = DMDAVecRestoreArray(fdstag->DA_CENTER,    fdstag->Center_PhaseProportions  [jj], &Center_PhaseProportions  [jj]); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(fdstag->DA_XY_POINTS, fdstag->XYPoints_PhaseProportions[jj], &XYPoints_PhaseProportions[jj]); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(fdstag->DA_XZ_POINTS, fdstag->XZPoints_PhaseProportions[jj], &XZPoints_PhaseProportions[jj]); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(fdstag->DA_YZ_POINTS, fdstag->YZPoints_PhaseProportions[jj], &YZPoints_PhaseProportions[jj]); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
*/
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAGProjectEffVisc"
PetscErrorCode FDSTAGProjectEffVisc(FDSTAG *fs, JacResCtx *jrctx, UserContext *usr)
{

	// project effective viscosity values to LaMEM data structures

	FDSTAG_MaterialProperties *fdstag;
	PetscScalar ***Center_EffectiveViscosity;
	PetscScalar ***XYPoints_EffectiveViscosity;
	PetscScalar ***XZPoints_EffectiveViscosity;
	PetscScalar ***YZPoints_EffectiveViscosity;
	PetscInt i, j, k, nx, ny, nz, sx, sy, sz, iter;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access solution variables
	fdstag = &usr->FDSTAG;

	// gain access to LaMEM storage scheme
	ierr = DMDAVecGetArray(fdstag->DA_CENTER,    fdstag->Center_EffectiveViscosity,   &Center_EffectiveViscosity);   CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fdstag->DA_XY_POINTS, fdstag->XYPoints_EffectiveViscosity, &XYPoints_EffectiveViscosity); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fdstag->DA_XZ_POINTS, fdstag->XZPoints_EffectiveViscosity, &XZPoints_EffectiveViscosity); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fdstag->DA_YZ_POINTS, fdstag->YZPoints_EffectiveViscosity, &YZPoints_EffectiveViscosity); CHKERRQ(ierr);

	//---------------
	// central points
	//---------------
	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		Center_EffectiveViscosity[k][j][i] = jrctx->svCell[iter++].svDev.eta;
	}
	END_STD_LOOP

	//---------------
	// xy edge points
	//---------------
	iter = 0;
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		XYPoints_EffectiveViscosity[k][j][i] = jrctx->svXYEdge[iter++].svDev.eta;
	}
	END_STD_LOOP

	//---------------
	// xz edge points
	//---------------
	iter = 0;
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		XZPoints_EffectiveViscosity[k][j][i] = jrctx->svXZEdge[iter++].svDev.eta;
	}
	END_STD_LOOP

	//---------------
	// yz edge points
	//---------------
	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		YZPoints_EffectiveViscosity[k][j][i] = jrctx->svYZEdge[iter++].svDev.eta;
	}
	END_STD_LOOP

	ierr = DMDAVecRestoreArray(fdstag->DA_CENTER,    fdstag->Center_EffectiveViscosity,   &Center_EffectiveViscosity);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fdstag->DA_XY_POINTS, fdstag->XYPoints_EffectiveViscosity, &XYPoints_EffectiveViscosity); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fdstag->DA_XZ_POINTS, fdstag->XZPoints_EffectiveViscosity, &XZPoints_EffectiveViscosity); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fdstag->DA_YZ_POINTS, fdstag->YZPoints_EffectiveViscosity, &YZPoints_EffectiveViscosity); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAGCopyResidual"
PetscErrorCode FDSTAGCopyResidual(FDSTAG *fs, JacResCtx *jrctx, UserContext *usr, Vec f, Vec g)
{
	// copy residual from FDSTAG to LaMEM collocated layout
	// NOTE: this is completely local operation!
	// WARNING: function assumes FREE-SLIP BOX!

	Field       ***f_dest;
	PetscScalar ***fx_src,  ***fy_src,  ***fz_src, ***gc_src, ***g_dest;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, my, mz;

	PetscErrorCode 	ierr;
	PetscFunctionBegin;

	// initialize maximal node index in all directions
	mx = fs->dsx.tnods - 1;
	my = fs->dsy.tnods - 1;
	mz = fs->dsz.tnods - 1;

	// zero out LaMEM residual vectors
	ierr = VecZeroEntries(f); CHKERRQ(ierr);
	ierr = VecZeroEntries(g); CHKERRQ(ierr);

	// access LaMEM residual vectors (destination)
	ierr = DMDAVecGetArray(usr->DA_Vel,  f, &f_dest); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(usr->DA_Pres, g, &g_dest); CHKERRQ(ierr);

	// access FDSTAG residual vectors (source)
	ierr = DMDAVecGetArray(fs ->DA_X,   jrctx->gfx, &fx_src); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs ->DA_Y,   jrctx->gfy, &fy_src); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs ->DA_Z,   jrctx->gfz, &fz_src); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs ->DA_CEN, jrctx->gc,  &gc_src); CHKERRQ(ierr);

	//---------
	// X points
	//---------
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		if(i != 0 && i != mx) f_dest[k][j][i].Vx = fx_src[k][j][i];
		else                  fx_src[k][j][i]    = 0.0;
	}
	END_STD_LOOP

	//---------
	// Y points
	//---------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		if(j != 0 && j != my) f_dest[k][j][i].Vy = fy_src[k][j][i];
		else                  fy_src[k][j][i]    = 0.0;

	}
	END_STD_LOOP

	//---------
	// Z points
	//---------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		if(k != 0 && k != mz) f_dest[k][j][i].Vz = fz_src[k][j][i];
		else                  fz_src[k][j][i]    = 0.0;

	}
	END_STD_LOOP

	//-------------------------------------
	// central points (divergence residaul)
	//-------------------------------------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		g_dest[k][j][i] = gc_src[k][j][i];
	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(usr->DA_Vel,  f,          &f_dest); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(usr->DA_Pres, g,          &g_dest); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs ->DA_X,    jrctx->gfx, &fx_src); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs ->DA_Y,    jrctx->gfy, &fy_src); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs ->DA_Z,    jrctx->gfz, &fz_src); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs ->DA_CEN,  jrctx->gc,  &gc_src); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "FDSTAGCopySolution"
PetscErrorCode FDSTAGCopySolution(FDSTAG *fs, JacResCtx *jrctx, UserContext *usr, Vec v, Vec p)
{
	// copy solution from LaMEM collocated layout to FDSTAG & set boundary values
	// WARNING: function assumes FREE-SLIP BOX!

	Field       ***v_src;
	PetscScalar ***p_src;
	PetscScalar ***vx_dest,  ***vy_dest,  ***vz_dest, ***p_dest;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, my, mz;

	PetscErrorCode 	ierr;
	PetscFunctionBegin;

	// initialize maximal node index in all directions
	mx = fs->dsx.tnods - 1;
	my = fs->dsy.tnods - 1;
	mz = fs->dsz.tnods - 1;

	// zero out global FDSTAG solution vectors
	ierr = VecZeroEntries(jrctx->gvx); CHKERRQ(ierr);
	ierr = VecZeroEntries(jrctx->gvy); CHKERRQ(ierr);
	ierr = VecZeroEntries(jrctx->gvz); CHKERRQ(ierr);
	ierr = VecZeroEntries(jrctx->gp);  CHKERRQ(ierr);

	// access LaMEM solution vectors (source)
	ierr = DMDAVecGetArray(usr->DA_Vel,  v, &v_src); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(usr->DA_Pres, p, &p_src); CHKERRQ(ierr);

	// access global FDSTAG solution vectors (destination)
	ierr = DMDAVecGetArray(fs->DA_X,   jrctx->gvx, &vx_dest); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Y,   jrctx->gvy, &vy_dest); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_Z,   jrctx->gvz, &vz_dest); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jrctx->gp,  &p_dest);  CHKERRQ(ierr);

	//---------
	// X points
	//---------
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		if(i != 0 && i != mx) vx_dest[k][j][i] = v_src[k][j][i].Vx;
	}
	END_STD_LOOP

	//---------
	// Y points
	//---------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		if(j != 0 && j != my) vy_dest[k][j][i] = v_src[k][j][i].Vy;
	}
	END_STD_LOOP

	//---------
	// Z points
	//---------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		if(k != 0 && k != mz) vz_dest[k][j][i] = v_src[k][j][i].Vz;
	}
	END_STD_LOOP

	//--------------------------
	// central points (pressure)
	//--------------------------
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		p_dest[k][j][i] = p_src[k][j][i];
	}
	END_STD_LOOP

	// restore access LaMEM solution vectors (source)
	ierr = DMDAVecRestoreArray(usr->DA_Vel,  v, &v_src); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(usr->DA_Pres, p, &p_src); CHKERRQ(ierr);

	// restore access global FDSTAG solution vectors
	ierr = DMDAVecRestoreArray(fs->DA_X,   jrctx->gvx, &vx_dest); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Y,   jrctx->gvy, &vy_dest); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_Z,   jrctx->gvz, &vz_dest); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jrctx->gp,  &p_dest);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

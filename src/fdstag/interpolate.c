//---------------------------------------------------------------------------
//.............   FDSTAG OUTPUT VECTOR INTERPOLATION ROUTINES   .............
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "interpolate.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "InterpXFaceCorner"
PetscErrorCode InterpXFaceCorner(FDSTAG *fs, Vec XFace, Vec Corner, InterpFlags iflag)
{
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, my, mz;
	PetscScalar cf, ***lXFace, ***lCorner, A1, A2, A3, A4, B1, B2, E1, E2;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access vectors
	ierr = DMDAVecGetArray(fs->DA_X,   XFace,  &lXFace);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_COR, Corner, &lCorner); CHKERRQ(ierr);

	// set index boundaries in Y & Z directions
	my = fs->dsy.tnods - 1;
	mz = fs->dsz.tnods - 1;

	// interpolate x-face vector to corners
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// access basic source values
		A1 = lXFace[k-1][j-1][i];
		A2 = lXFace[k-1][j  ][i];
		A3 = lXFace[k  ][j-1][i];
		A4 = lXFace[k  ][j  ][i];

		if(iflag.use_bound != PETSC_TRUE)
		{
			// set ghost values on boundaries if not defined
			if(j == 0)	{ A1 = A2; A3 = A4; }
			if(j == my) { A2 = A1; A4 = A3; }
			if(k == 0)	{ A1 = A3; A2 = A4; }
			if(k == mz) { A3 = A1; A4 = A2; }
		}
		else
		{
			// interpolate given ghost values to boundary edges
			if(j == 0  && k == 0)   A1 = 0.5*(A2 + A3);
			if(j == my && k == 0)   A2 = 0.5*(A1 + A4);
			if(j == 0  && k == mz)  A3 = 0.5*(A1 + A4);
			if(j == my && k == mz)  A4 = 0.5*(A2 + A3);
		}

		// get weight coefficients
		E1 = WEIGHT_NODE(j, sy, fs->dsy); B1 = 1.0 - E1;
		E2 = WEIGHT_NODE(k, sz, fs->dsz); B2 = 1.0 - E2;

		// interpolate in Y-Z plane
		cf = A1*B1*B2 + A2*E1*B2 + A3*B1*E2 + A4*E1*E2;

		// store
		if(iflag.update != PETSC_TRUE) lCorner[k][j][i]  = cf;
		else                           lCorner[k][j][i] += cf;

	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_X,   XFace,  &lXFace);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_COR, Corner, &lCorner);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "InterpYFaceCorner"
PetscErrorCode InterpYFaceCorner(FDSTAG *fs, Vec YFace, Vec Corner, InterpFlags iflag)
{
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, mz;
	PetscScalar cf, ***lYFace, ***lCorner, A1, A2, A3, A4, B1, B2, E1, E2;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access vectors
	ierr = DMDAVecGetArray(fs->DA_Y,   YFace,  &lYFace);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_COR, Corner, &lCorner); CHKERRQ(ierr);

	// set index boundaries in X & Z directions
	mx = fs->dsx.tnods - 1;
	mz = fs->dsz.tnods - 1;

	// interpolate y-face vector to corners
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// access basic source values
		A1 = lYFace[k-1][j][i-1];
		A2 = lYFace[k-1][j][i  ];
		A3 = lYFace[k  ][j][i-1];
		A4 = lYFace[k  ][j][i  ];

		if(iflag.use_bound != PETSC_TRUE)
		{
			// set ghost values on boundaries if not defined
			if(i == 0)	{ A1 = A2; A3 = A4; }
			if(i == mx) { A2 = A1; A4 = A3; }
			if(k == 0)	{ A1 = A3; A2 = A4; }
			if(k == mz) { A3 = A1; A4 = A2; }
		}
		else
		{
			// interpolate given ghost values to boundary edges
			if(i == 0  && k == 0)   A1 = 0.5*(A2 + A3);
			if(i == mx && k == 0)   A2 = 0.5*(A1 + A4);
			if(i == 0  && k == mz)  A3 = 0.5*(A1 + A4);
			if(i == mx && k == mz)  A4 = 0.5*(A2 + A3);
		}

		// get weight coefficients
		E1 = WEIGHT_NODE(i, sx, fs->dsx); B1 = 1.0 - E1;
		E2 = WEIGHT_NODE(k, sz, fs->dsz); B2 = 1.0 - E2;

		// interpolate in X-Z plane
		cf = A1*B1*B2 + A2*E1*B2 + A3*B1*E2 + A4*E1*E2;

		// store
		if(iflag.update != PETSC_TRUE) lCorner[k][j][i]  = cf;
		else                           lCorner[k][j][i] += cf;

	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_Y,   YFace,  &lYFace);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_COR, Corner, &lCorner);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "InterpZFaceCorner"
PetscErrorCode InterpZFaceCorner(FDSTAG *fs, Vec ZFace, Vec Corner, InterpFlags iflag)
{
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, my;
	PetscScalar cf, ***lZFace, ***lCorner, A1, A2, A3, A4, B1, B2, E1, E2;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access vectors
	ierr = DMDAVecGetArray(fs->DA_Z,   ZFace,  &lZFace);  CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_COR, Corner, &lCorner); CHKERRQ(ierr);

	// set index boundaries in X & Y directions
	mx = fs->dsx.tnods - 1;
	my = fs->dsy.tnods - 1;

	// interpolate z-face vector to corners
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// access basic source values
		A1 = lZFace[k][j-1][i-1];
		A2 = lZFace[k][j-1][i  ];
		A3 = lZFace[k][j  ][i-1];
		A4 = lZFace[k][j  ][i  ];

		if(iflag.use_bound != PETSC_TRUE)
		{
			// set ghost values on boundaries if not defined
			if(i == 0)	{ A1 = A2; A3 = A4; }
			if(i == mx) { A2 = A1; A4 = A3; }
			if(j == 0)	{ A1 = A3; A2 = A4; }
			if(j == my) { A3 = A1; A4 = A2; }
		}
		else
		{
			// interpolate given ghost values to boundary edges
			if(i == 0  && j == 0)   A1 = 0.5*(A2 + A3);
			if(i == mx && j == 0)   A2 = 0.5*(A1 + A4);
			if(i == 0  && j == my)  A3 = 0.5*(A1 + A4);
			if(i == mx && j == my)  A4 = 0.5*(A2 + A3);
		}

		// get weight coefficients
		E1 = WEIGHT_NODE(i, sx, fs->dsx); B1 = 1.0 - E1;
		E2 = WEIGHT_NODE(j, sy, fs->dsy); B2 = 1.0 - E2;

		// interpolate in X-Y plane
		cf = A1*B1*B2 + A2*E1*B2 + A3*B1*E2 + A4*E1*E2;

		// store
		if(iflag.update != PETSC_TRUE) lCorner[k][j][i]  = cf;
		else                           lCorner[k][j][i] += cf;

	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_Z,   ZFace,  &lZFace);   CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_COR, Corner, &lCorner);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "InterpCenterCorner"
PetscErrorCode InterpCenterCorner(FDSTAG *fs, Vec Center, Vec Corner, InterpFlags iflag)
{
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, my, mz, I1, I2, J1, J2, K1, K2;
	PetscScalar cf, ***lCenter, ***lCorner, A1, A2, A3, A4, A5, A6, A7, A8, B1, B2, B3, E1, E2, E3;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access vectors
	ierr = DMDAVecGetArray(fs->DA_CEN, Center, &lCenter); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_COR, Corner, &lCorner); CHKERRQ(ierr);

	// set index boundaries in all directions
	mx = fs->dsx.tnods - 1;
	my = fs->dsy.tnods - 1;
	mz = fs->dsz.tnods - 1;

	// interpolate center vector to corners
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// initialize indices
		I1 = i;
		I2 = i-1;
		J1 = j;
		J2 = j-1;
		K1 = k;
		K2 = k-1;

		if(iflag.use_bound != PETSC_TRUE)
		{
			// check index bounds if ghost points are undefined
			if(I1 == mx) I1--;
			if(I2 == -1) I2++;
			if(J1 == my) J1--;
			if(J2 == -1) J2++;
			if(K1 == mz) K1--;
			if(K2 == -1) K2++;
		}

		// access basic source values
		A1 = lCenter[K2][J2][I2];
		A2 = lCenter[K2][J2][I1];
		A3 = lCenter[K2][J1][I2];
		A4 = lCenter[K2][J1][I1];
		A5 = lCenter[K1][J2][I2];
		A6 = lCenter[K1][J2][I1];
		A7 = lCenter[K1][J1][I2];
		A8 = lCenter[K1][J1][I1];

		// get weight coefficients
		E1 = WEIGHT_NODE(i, sx, fs->dsx); B1 = 1.0 - E1;
		E2 = WEIGHT_NODE(j, sy, fs->dsy); B2 = 1.0 - E2;
		E3 = WEIGHT_NODE(k, sz, fs->dsz); B3 = 1.0 - E3;

		// interpolate in 3D cube
		cf = A1*B1*B2*B3 + A2*E1*B2*B3 + A3*B1*E2*B3 + A4*E1*E2*B3
		+    A5*B1*B2*E3 + A6*E1*B2*E3 + A7*B1*E2*E3 + A8*E1*E2*E3;

		// store
		if(iflag.update != PETSC_TRUE) lCorner[k][j][i]  = cf;
		else                           lCorner[k][j][i] += cf;

	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_CEN, Center, &lCenter);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_COR, Corner, &lCorner);  CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "InterpXYEdgeCorner"
PetscErrorCode InterpXYEdgeCorner(FDSTAG *fs, Vec XYEdge, Vec Corner, InterpFlags iflag)
{
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mz, K1, K2;
	PetscScalar cf, ***lXYEdge, ***lCorner, A1, A2, B1, E1;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access vectors
	ierr = DMDAVecGetArray(fs->DA_XY,  XYEdge, &lXYEdge); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_COR, Corner, &lCorner); CHKERRQ(ierr);

	// set index boundaries in Z direction
	mz = fs->dsz.tnods - 1;

	// interpolate xy-edge vector to corners
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// set index bounds
		K1 = k;   if(K1 == mz) K1--;
		K2 = k-1; if(K2 == -1) K2++;

		// access basic source values
		A1 = lXYEdge[K2][j][i];
		A2 = lXYEdge[K1][j][i];

		// get weight coefficients
		E1 = WEIGHT_NODE(k, sz, fs->dsz); B1 = 1.0 - E1;

		// interpolate along Z-edge
		cf = A1*B1 + A2*E1;

		// store
		if(iflag.update != PETSC_TRUE) lCorner[k][j][i]  = cf;
		else                           lCorner[k][j][i] += cf;

	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_XY,  XYEdge, &lXYEdge);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_COR, Corner, &lCorner);  CHKERRQ(ierr);


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "InterpXZEdgeCorner"
PetscErrorCode InterpXZEdgeCorner(FDSTAG *fs, Vec XZEdge, Vec Corner, InterpFlags iflag)
{
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, my, J1, J2;
	PetscScalar cf, ***lXZEdge, ***lCorner, A1, A2, B1, E1;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access vectors
	ierr = DMDAVecGetArray(fs->DA_XZ,  XZEdge, &lXZEdge); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_COR, Corner, &lCorner); CHKERRQ(ierr);

	// set index boundaries in Y direction
	my = fs->dsy.tnods - 1;

	// interpolate xz-edge vector to corners
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// set index bounds
		J1 = j;   if(J1 == my) J1--;
		J2 = j-1; if(J2 == -1) J2++;

		// access basic source values
		A1 = lXZEdge[k][J2][i];
		A2 = lXZEdge[k][J1][i];

		// get weight coefficients
		E1 = WEIGHT_NODE(j, sy, fs->dsy); B1 = 1.0 - E1;

		// interpolate along Y-edge
		cf = A1*B1 + A2*E1;

		// store
		if(iflag.update != PETSC_TRUE) lCorner[k][j][i]  = cf;
		else                           lCorner[k][j][i] += cf;

	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_XZ,  XZEdge, &lXZEdge);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_COR, Corner, &lCorner);  CHKERRQ(ierr);


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "InterpYZEdgeCorner"
PetscErrorCode InterpYZEdgeCorner(FDSTAG *fs, Vec YZEdge, Vec Corner, InterpFlags iflag)
{
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, mx, I1, I2;
	PetscScalar cf, ***lYZEdge, ***lCorner, A1, A2, B1, E1;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access vectors
	ierr = DMDAVecGetArray(fs->DA_YZ,  YZEdge, &lYZEdge); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_COR, Corner, &lCorner); CHKERRQ(ierr);

	// set index boundaries in X direction
	mx = fs->dsx.tnods - 1;

	// interpolate yz-edge vector to corners
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// set index bounds
		I1 = i;   if(I1 == mx) I1--;
		I2 = i-1; if(I2 == -1) I2++;

		// access basic source values
		A1 = lYZEdge[k][j][I2];
		A2 = lYZEdge[k][j][I1];

		// get weight coefficients
		E1 = WEIGHT_NODE(i, sx, fs->dsx); B1 = 1.0 - E1;

		// interpolate along X-edge
		cf = A1*B1 + A2*E1;

		// store
		if(iflag.update != PETSC_TRUE) lCorner[k][j][i]  = cf;
		else                           lCorner[k][j][i] += cf;

	}
	END_STD_LOOP

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_YZ,  YZEdge, &lYZEdge);  CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_COR, Corner, &lCorner);  CHKERRQ(ierr);


	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

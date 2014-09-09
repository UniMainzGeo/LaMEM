/*
LaMEM
Lithosphere and Mantle Evolution Model

A 3D numerical code that can be used to model geodynamical processes such as
mantle-lithosphere interaction.

Originally developed by:

Boris J.P. Kaus (Uni-Mainz, Germany). 2007-
kaus@uni-mainz.de

Further-development:

Dave May (ETH Zurich, Switzerland). 2009-
dave.may@erdw.ethz.ch

Tobias Baumann (Uni-Mainz, Germany). 2011-
baumann@uni-mainz.de

Anton Popov (Uni-Mainz, Germany). 2011-
popov@uni-mainz.de

ApplyBoundaryConditions.c, contains the following subroutines:

DefineInternalBC		-	Computes the nodal indices, where internal BC's are applied, from coordinates (defined via command line) 

$Id:$
$Date:$

This version is compatible with PETSc version 3.3-p2
The code cannot be distributed, used for commercial purposes or handed on,
without the explicit agreement of Boris Kaus.
*/


#include "LaMEM.h"
#include "ApplyBoundaryConditions.h"

/*==========================================================================================================*/
void DefineInternalBC( UserContext *user )
{
	PetscScalar    	dy;
	PetscInt		i, j, k, m, n, p, mstart, nstart, pstart, nx, ny, nz;
	Vec 			gc, global;
	DMDACoor3d		***coors;
	DM				cda;
	
	DMGetCoordinateDM(user->DA_Vel,&cda);
	DMGetCoordinatesLocal(user->DA_Vel,&gc);
	DMDAVecGetArray(cda,gc,&coors);
	DMDAGetCorners(cda,&mstart,&nstart,&pstart,&m,&n,&p);
	DMDAGetInfo(user->DA_Vel,0,&nx,&ny,&nz,0,0,0,0,0,0,0,0,0);
	dy = user->L/((double) (ny-1));
	for (i=mstart; i<mstart+m; i++) {
		for (j=nstart; j<nstart+n; j++) {
			for (k=pstart; k<pstart+p; k++) {
				if ( (coors[k][j][i].y-(user->internalBC_coord/user->Characteristic.Length)) <= dy ) {
					user->internalBC = j;
					if (user->internalBC % 2 == 1){
						user->internalBC = user->internalBC+1;
					}
					user->internalBC_node = user->internalBC;
					user->internalBC_backel = user->internalBC_node/2 - 1;
					user->internalBC_frontel = user->internalBC_node/2;			
					/*goto endnestedloop;*/
				}
			}
		}
	}

	/*endnestedloop:*/

	PetscPrintf(PETSC_COMM_WORLD," dy: %g \n",dy);
	PetscPrintf(PETSC_COMM_WORLD," internalBC_coord: %lld \n",(LLD)(user->internalBC_coord));
	PetscPrintf(PETSC_COMM_WORLD," internalBC: %lld \n",(LLD)(user->internalBC));
	PetscPrintf(PETSC_COMM_WORLD," internalBC_node: %lld \n",(LLD)(user->internalBC_node));
	PetscPrintf(PETSC_COMM_WORLD," internalBC_backel: %lld \n",(LLD)(user->internalBC_backel));
	PetscPrintf(PETSC_COMM_WORLD," internalBC_frontel: %lld \n",(LLD)(user->internalBC_frontel));
	PetscPrintf(PETSC_COMM_WORLD," zdepth_BC_node: %lld \n",(LLD)(user->zdepth_BC_node));
	PetscPrintf(PETSC_COMM_WORLD," zdepth_BC_el: %lld \n",(LLD)(user->zdepth_BC_el));

	DMDAVecRestoreArray(cda,gc,&coors);
	DMGetCoordinates(user->DA_Vel,&global);
	DMLocalToGlobalBegin(cda,gc,INSERT_VALUES,global);
	DMLocalToGlobalEnd  (cda,gc,INSERT_VALUES,global);
	
	//VecDestroy(gc);
	//VecDestroy(global);
	//DMDestroy(cda);

}
/*==========================================================================================================*/

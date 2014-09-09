/*
LaMEM
Lithosphere and Mantle Evolution Model

A 3D numerical code that can be used to model geodynamical processes such as
mantle-lithosphere interaction.

Originally developed by:

Boris J.P. Kaus (Uni-Mainz, Germany). 2007-
kaus@uni-mainz.de

Further-development:

Dave May (ETH Zurich, Switzerland). 2009-2011
dave.may@erdw.ethz.ch

Tobias Baumann (Uni-Mainz, Germany). 2011-
baumann@uni-mainz.de

Anton Popov (Uni-Mainz, Germany). 2011-
popov@uni-mainz.de

Assembly_FDSTAG.c contains routines required for the staggered FD stokes solver:

ComputeStiffnessMatrixes_FDSTAG 	- computes the Staggered FD stiffness matrix

$Id:$
$Date:$

This version is compatible with PETSc version 3.3-p2
The code cannot be distributed, used for commercial purposes or handed on,
without the explicit agreement of Boris Kaus.
*/

//#define LAMEM_DEBUG_FDSTAG


#include "LaMEM.h"
#include "Assembly_FDSTAG.h"
#include "Elements.h"
#include "Addition_FDSTAG.h"

/*==========================================================================================================
 * DAGetGlobalIndex computes the global index from i,j,k coordinates and a dof
 *
 * PetscInt i,j,k,d;  (i,j,k) = GLOBAL node indices, d = dof index
 *
 * Created by dmay on 11/10/10.
 * Copyright 2010 Geophysical Fluid Dynamics. All rights reserved.
 */
#undef __FUNCT__
#define __FUNCT__ "DAGetGlobalIndex"
PetscErrorCode DAGetGlobalIndex(DM da,PetscInt i,PetscInt j,PetscInt k,PetscInt d,PetscInt *gidx)
{
	const PetscInt      *indices;
	PetscInt 		    dof;
	PetscInt 		    gnx,gny,gnz, M,N,P; 	/* ghost nodes in (x,y,z), global nodes in x,y,z */
	PetscInt 		    si,sj,sk; 				/* start indices */
	PetscInt		    global_idx;
	ISLocalToGlobalMapping ltogm;

	PetscErrorCode 	ierr;
	PetscFunctionBegin;

	ierr = DMDAGetGhostCorners(da, &si,&sj,&sk, &gnx,&gny,&gnz );	CHKERRQ(ierr);
	ierr = DMDAGetInfo(da,0, &M,&N,&P, 0,0,0, &dof,0,0,0,0,0 );		CHKERRQ(ierr);

	/* check valid local range */
	if(i<0) 		SETERRQ3(PETSC_COMM_WORLD, PETSC_ERR_USER,"i = %lld outside valid global range [%lld <= i <= %lld] \n",(LLD)i,0LL,(LLD)(M-1) );
	if(i>=M) 		SETERRQ3(PETSC_COMM_WORLD, PETSC_ERR_USER,"i = %lld outside valid global range [%lld <= i <= %lld] \n",(LLD)i,0LL,(LLD)(M-1) );

	if(j<0) 		SETERRQ3(PETSC_COMM_WORLD, PETSC_ERR_USER,"j = %lld outside valid global range [%lld <= j <= %lld] \n",(LLD)j,0LL,(LLD)(N-1) );
	if(j>=N) 		SETERRQ3(PETSC_COMM_WORLD, PETSC_ERR_USER,"j = %lld outside valid global range [%lld <= j <= %lld] \n",(LLD)j,0LL,(LLD)(N-1) );

	if(k<0) 		SETERRQ3(PETSC_COMM_WORLD, PETSC_ERR_USER,"k = %lld outside valid global range [%lld <= k <= %lld] \n",(LLD)k,0LL,(LLD)(P-1) );
	if(k>=P) 		SETERRQ3(PETSC_COMM_WORLD, PETSC_ERR_USER,"k = %lld outside valid global range [%lld <= k <= %lld] \n",(LLD)k,0LL,(LLD)(P-1) );

	/* check valid local range */
	if(i<si) 		SETERRQ3(PETSC_COMM_WORLD, PETSC_ERR_USER,"i = %lld outside valid local range [%lld <= i <= %lld] \n",(LLD)i,(LLD)si,(LLD)(si+gnx-1) );
	if(i>=si+gnx) 	SETERRQ3(PETSC_COMM_WORLD, PETSC_ERR_USER,"i = %lld outside valid local range [%lld <= i <= %lld] \n",(LLD)i,(LLD)si,(LLD)(si+gnx-1) );

	if(j<sj) 		SETERRQ3(PETSC_COMM_WORLD, PETSC_ERR_USER,"j = %lld outside valid local range [%lld <= j <= %lld] \n",(LLD)j,(LLD)sj,(LLD)(sj+gny-1) );
	if(j>=sj+gny) 	SETERRQ3(PETSC_COMM_WORLD, PETSC_ERR_USER,"j = %lld outside valid local range [%lld <= j <= %lld] \n",(LLD)j,(LLD)sj,(LLD)(sj+gny-1) );

	if(k<sk) 		SETERRQ3(PETSC_COMM_WORLD, PETSC_ERR_USER,"k = %lld outside valid local range [%lld <= k <= %lld] \n",(LLD)k,(LLD)sk,(LLD)(sk+gnz-1) );
	if(k>=sk+gnz) 	SETERRQ3(PETSC_COMM_WORLD, PETSC_ERR_USER,"k = %lld outside valid local range [%lld <= k <= %lld] \n",(LLD)k,(LLD)sk,(LLD)(sk+gnz-1) );

//	ierr = DMDAGetGlobalIndices(da, &nidx, &indices); CHKERRQ(ierr);

	ierr = DMGetLocalToGlobalMapping(da, &ltogm);             CHKERRQ(ierr);
	ierr = ISLocalToGlobalMappingGetIndices(ltogm, &indices); CHKERRQ(ierr);

	global_idx 	= indices[ dof*((i-si) + (j-sj)*gnx + (k-sk)*gnx*gny) + d ];

	*gidx = global_idx;

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/

#undef __FUNCT__
#define __FUNCT__ "FDSTAG_GetStencilValues"
PetscErrorCode FDSTAG_GetStencilValues(DM dav,DM dap,PetscInt vsL,MatStencil vs[],PetscInt psL,MatStencil ps[],PetscInt vgidx[],PetscInt pgidx[])
{
	PetscInt       n;
	PetscErrorCode ierr;

	PetscFunctionBegin;
	for (n=0; n<vsL; n++) {
		vgidx[n] = -1;
		if ( (vs[n].i>=0) && (vs[n].j>=0) && (vs[n].k>=0) ) {
			ierr = DAGetGlobalIndex(dav,vs[n].i,vs[n].j,vs[n].k,vs[n].c,&vgidx[n]);CHKERRQ(ierr);
		}
	}
	for (n=0; n<psL; n++) {
		pgidx[n] = -1;
		if ( (ps[n].i>=0) && (ps[n].j>=0) && (ps[n].k>=0) ) {
			ierr = DAGetGlobalIndex(dap,ps[n].i,ps[n].j,ps[n].k,ps[n].c,&pgidx[n]);CHKERRQ(ierr);
		}
	}
	PetscFunctionReturn(0);
}

/*==========================================================================================================*/
/* Compute the stiffness matrix for a given level with a staggered-grid finite difference formulation using
 * the 'classical' viscous approximation to the Jacobian.
 * */
#undef __FUNCT__
#define __FUNCT__ "ComputeStiffnessMatrixes_FDSTAG"
PetscErrorCode ComputeStiffnessMatrixes_FDSTAG(
		DM da, DM da_pres,
		Mat VV_MAT, Mat PP_MAT, Mat PV_MAT, Mat VP_MAT,
		UserContext *user, PetscScalar dt, Vec rhs_add_BC,
		Vec rhs_add_BC_loc, Vec pv_rhs_push, Vec ViscosityScaling, Mat approx_S )
{
	PetscMPIInt     rank;
	PetscBool		eliminate_stickyair_from_system, subtract_volumetric_strainrates;
	PetscInt        nel_x, nel_y, nel_z, nnode_x, nnode_y, nnode_z, number_coefficients;
	PetscInt        i,j,k,xm,ym,zm,xs,ys,zs;
	PetscInt        iel_x, iel_y, iel_z;
	PetscInt        xmp,ymp,zmp,xsp,ysp,zsp,zs_FreeSurface;
	PetscInt 		*rowid_DirichletBC_array, *rowid_addPoints_array, numDirichlet, numAddpoints;
	PetscInt        AirPhase;
	const PetscInt  *P_globalindices, *Vel_globalindices;
	PetscInt        FreeSurfaceCells[7];
	PetscScalar     dx_P[2], dy_P[2], dz_P[2], dx_vec[7], dy_vec[7], dz_vec[7], Eta_Center[7], Eta_XY[4], Eta_XZ[4], Eta_YZ[4], dRho_dxdydz[3],z_FreeSurface;
	PetscScalar     ***viscosity_center, ***viscosity_XY,  ***viscosity_XZ,  ***viscosity_YZ, ***density_center, ***PhaseProportionsAir_Center, ***LocalSurfaceTopography;
	PetscScalar 	z_center, FreeSurface_Fraction;
	PetscScalar 	ScalingParameterDiagonal;
	DM              cda, cda_SurfaceTopo;
	Vec             gc, Viscosity_Center_local, Viscosity_XY_local, Viscosity_XZ_local, Viscosity_YZ_local, Density_Center_local, PhaseProportionsAir_Center_Vec;
	Vec 			LocalSurfaceTopography_vec, gc_SurfaceTopo;
	DMDACoor3d        ***coords;
	Field			 	    ***rhs_for_BoundaryConditions;
	PetscLogDouble  cputime_start, cputime_end;
	MatStencil      row_pp[1];   //  indices for pp matrix
	PetscScalar     v_pp[1];     //	values  for pp matrix
	MatStencil      col_pv[6];                                      //  indices for pv matrix
	PetscScalar     v_pv[6];                                        //	values  for pv matrix
	MatStencil      col_vp[2];                                                       //  indices for vp matrix
	PetscScalar     v_vp[2];                                                         //	values  for vp matrix
	MatStencil      row_vv[1],col_vv[32];                      //  indices for vv matrix
	PetscScalar     v_vv[32];                                  //	values  for vv matrix
	PetscInt        rowidx, vgidx[32],pgidx[2];
	PetscLogDouble  t0,t1;
	PetscErrorCode  ierr;
	DMDACoor3d		***coors_SurfaceTopo;
	Mat				PV_MAT_push;
	PetscInt		m,n;
	ISLocalToGlobalMapping P_ltogm, Vel_ltogm;

	PetscFunctionBegin;



	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	// Create matrix for pushing BC
	ierr = MatCreate(PETSC_COMM_WORLD,&PV_MAT_push);CHKERRQ(ierr);
	ierr = MatGetSize(PV_MAT,&m,&n);CHKERRQ(ierr);
	ierr = MatSetSizes(PV_MAT_push,PETSC_DECIDE, PETSC_DECIDE,m,n); CHKERRQ(ierr);
	ierr = MatSetType(PV_MAT_push,	MATAIJ); CHKERRQ(ierr);
	ierr = MatSeqAIJSetPreallocation(PV_MAT_push,6,PETSC_NULL); CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(PV_MAT_push,6,PETSC_NULL,6,PETSC_NULL); CHKERRQ(ierr); // preallocation
	ierr = MatZeroEntries(PV_MAT_push);   CHKERRQ(ierr);

	ierr = DMGetCoordinateDM(da,&cda); CHKERRQ(ierr);	//coordinates
	ierr = DMGetCoordinatesLocal(da,&gc); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,gc,&coords); CHKERRQ(ierr);

	ierr = DMGetLocalVector(da,&rhs_add_BC_loc); CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da, rhs_add_BC, INSERT_VALUES, rhs_add_BC_loc); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(da,rhs_add_BC, INSERT_VALUES, rhs_add_BC_loc); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da,rhs_add_BC_loc,&rhs_for_BoundaryConditions); CHKERRQ(ierr);

	/* print some info */
	ierr = DMDAGetInfo(user->DA_Processors, 0, &nel_x, &nel_y, &nel_z,	0,0,0,0,0,0,0,0,0);	CHKERRQ(ierr); 	// # of elements in all directions
	ierr = DMDAGetInfo(da, 0, &nnode_x, &nnode_y, &nnode_z,	0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);  // # of nodes 	 in all directions
	PetscPrintf(PETSC_COMM_WORLD,"#  Forming FD stiffness matrix mx,my,mz=(%lld, %lld, %lld) ...  \n",(LLD)nnode_x,(LLD)nnode_y,(LLD)nnode_z);

	PetscTime(&cputime_start);
	/* Corners */
	ierr = DMDAGetCorners(user->DA_Processors,&xsp,&ysp,&zsp,&xmp,&ymp,&zmp); CHKERRQ(ierr);
	ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm); CHKERRQ(ierr);

	/* Obtain mapping from current proc to global numbering for both velocity and pressure */
//	ierr = DMDAGetGlobalIndices(da_pres,&np,&P_globalindices); CHKERRQ(ierr);
//	ierr = DMDAGetGlobalIndices(da,     &nvel,&Vel_globalindices); CHKERRQ(ierr);

	ierr = DMGetLocalToGlobalMapping(da_pres, &P_ltogm);   CHKERRQ(ierr);
	ierr = DMGetLocalToGlobalMapping(da,      &Vel_ltogm); CHKERRQ(ierr);

	ierr = ISLocalToGlobalMappingGetIndices(P_ltogm,   &P_globalindices);   CHKERRQ(ierr);
	ierr = ISLocalToGlobalMappingGetIndices(Vel_ltogm, &Vel_globalindices); CHKERRQ(ierr);

	/* Get viscosity @ center points, including ghost node values */
	ierr = DMGetLocalVector     (user->FDSTAG.DA_CENTER,&Viscosity_Center_local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin (user->FDSTAG.DA_CENTER,user->FDSTAG.Center_EffectiveViscosity,INSERT_VALUES,Viscosity_Center_local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd   (user->FDSTAG.DA_CENTER,user->FDSTAG.Center_EffectiveViscosity,INSERT_VALUES,Viscosity_Center_local); CHKERRQ(ierr);

	/* Get density @ center points, including ghost node values */
	ierr = DMGetLocalVector     (user->FDSTAG.DA_CENTER,&Density_Center_local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin (user->FDSTAG.DA_CENTER,user->FDSTAG.Center_Density,INSERT_VALUES,Density_Center_local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd   (user->FDSTAG.DA_CENTER,user->FDSTAG.Center_Density,INSERT_VALUES,Density_Center_local); CHKERRQ(ierr);


	/* Get viscosity @ corner points, including ghost node values */
	ierr = DMGetLocalVector     (user->FDSTAG.DA_XY_POINTS,&Viscosity_XY_local);	CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin (user->FDSTAG.DA_XY_POINTS,user->FDSTAG.XYPoints_EffectiveViscosity,INSERT_VALUES,Viscosity_XY_local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd   (user->FDSTAG.DA_XY_POINTS,user->FDSTAG.XYPoints_EffectiveViscosity,INSERT_VALUES,Viscosity_XY_local); CHKERRQ(ierr);


	ierr = DMGetLocalVector     (user->FDSTAG.DA_XZ_POINTS,&Viscosity_XZ_local);	CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin (user->FDSTAG.DA_XZ_POINTS,user->FDSTAG.XZPoints_EffectiveViscosity,INSERT_VALUES,Viscosity_XZ_local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd   (user->FDSTAG.DA_XZ_POINTS,user->FDSTAG.XZPoints_EffectiveViscosity,INSERT_VALUES,Viscosity_XZ_local); CHKERRQ(ierr);

	ierr = DMGetLocalVector     (user->FDSTAG.DA_YZ_POINTS,&Viscosity_YZ_local);	CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin (user->FDSTAG.DA_YZ_POINTS,user->FDSTAG.YZPoints_EffectiveViscosity,INSERT_VALUES,Viscosity_YZ_local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd   (user->FDSTAG.DA_YZ_POINTS,user->FDSTAG.YZPoints_EffectiveViscosity,INSERT_VALUES,Viscosity_YZ_local); CHKERRQ(ierr);

	/* Get arrays for viscosity @ center, XY, XZ & YZ points */
	ierr = DMDAVecGetArray(user->FDSTAG.DA_CENTER,		Viscosity_Center_local, &viscosity_center ); CHKERRQ(ierr); // Viscosity @ center
	ierr = DMDAVecGetArray(user->FDSTAG.DA_CENTER,		Density_Center_local,   &density_center ); CHKERRQ(ierr); 	// Density @ center
	ierr = DMDAVecGetArray(user->FDSTAG.DA_XY_POINTS, 	Viscosity_XY_local,     &viscosity_XY     ); CHKERRQ(ierr); // Viscosity @ Sxy points
	ierr = DMDAVecGetArray(user->FDSTAG.DA_XZ_POINTS, 	Viscosity_XZ_local,     &viscosity_XZ     ); CHKERRQ(ierr); // Viscosity @ Sxz points
	ierr = DMDAVecGetArray(user->FDSTAG.DA_YZ_POINTS, 	Viscosity_YZ_local,     &viscosity_YZ     ); CHKERRQ(ierr); // Viscosity @ Syz points

	/* In case we use an internal free surface, we need this: */
	if (user->ErosionParameters.UseInternalFreeSurface==1){
		AirPhase 	=	user->ErosionParameters.StickyAirPhase;
	}
	else{
		AirPhase	=	0;
	}
	ierr = DMGetLocalVector     (user->FDSTAG.DA_CENTER,&PhaseProportionsAir_Center_Vec);	CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin (user->FDSTAG.DA_CENTER,	user->FDSTAG.Center_PhaseProportions[AirPhase],INSERT_VALUES,PhaseProportionsAir_Center_Vec); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd   (user->FDSTAG.DA_CENTER,	user->FDSTAG.Center_PhaseProportions[AirPhase],INSERT_VALUES,PhaseProportionsAir_Center_Vec); CHKERRQ(ierr);
	ierr = DMDAVecGetArray		(user->FDSTAG.DA_CENTER,  	PhaseProportionsAir_Center_Vec, 	&PhaseProportionsAir_Center	 ); CHKERRQ(ierr);

	/* Get free surface information [to correct stencil for internal free surface] */
	ierr = DMGetCoordinateDM(user->DA_SurfaceTopography,						&cda_SurfaceTopo		); 	CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(user->DA_SurfaceTopography,				&gc_SurfaceTopo			); 	CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda_SurfaceTopo,gc_SurfaceTopo,						&coors_SurfaceTopo		); 	CHKERRQ(ierr);

	ierr = DMGetLocalVector(user->DA_SurfaceTopography,	&LocalSurfaceTopography_vec); CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(user->DA_SurfaceTopography, user->SurfaceTopography, INSERT_VALUES, LocalSurfaceTopography_vec); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(user->DA_SurfaceTopography,   user->SurfaceTopography, INSERT_VALUES, LocalSurfaceTopography_vec); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(user->DA_SurfaceTopography,LocalSurfaceTopography_vec, 	&LocalSurfaceTopography	);	CHKERRQ(ierr);
	ierr = DMDAGetCorners(user->DA_SurfaceTopography,PETSC_NULL,PETSC_NULL,&zs_FreeSurface,PETSC_NULL,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);


	/* Do we eliminate sticky air or not? */
	eliminate_stickyair_from_system = PETSC_FALSE;
	ierr = PetscOptionsGetBool( PETSC_NULL, "-eliminate_stickyair_from_system", &eliminate_stickyair_from_system,PETSC_NULL ); CHKERRQ(ierr);

	/* Do we extract the volumetric component of strain rates or not (Exx_dev = Exx-1/3*(Exx+Eyy+Ezz), or Exx=Exx)	 */
	subtract_volumetric_strainrates = PETSC_FALSE;
	ierr = PetscOptionsGetBool( PETSC_NULL, "-subtract_volumetric_strainrates", &subtract_volumetric_strainrates,PETSC_NULL ); CHKERRQ(ierr);
	if (subtract_volumetric_strainrates){
		number_coefficients = 32;			// we have more equations as we subtract the volumetric components
	}
	else{
		number_coefficients = 20;
	}



	PetscTime(&t0);
	/* Make a loop over all local elements, construct the element stiffness matrix */
	for (iel_z=zsp; iel_z<zsp+zmp; iel_z++) {
		for (iel_y=ysp; iel_y<ysp+ymp; iel_y++) {
			for (iel_x=xsp; iel_x<xsp+xmp; iel_x++) {
				i = iel_x; j = iel_y; k = iel_z;		// Initialize variables

				/* Extract coordinates and spacing of the local element & neighboring elements
				 * (required in order to form the FD discretization)
				 */
				ierr = FDSTAG_ComputeSpacing(coords, iel_x,iel_y,iel_z, user, dx_vec, dy_vec, dz_vec, dx_P, dy_P, dz_P, &z_center); CHKERRQ(ierr);


				/* Extract material parameters */
				ierr = FDSTAG_ExtractMaterialParameters(viscosity_center, viscosity_XY, viscosity_YZ, viscosity_XZ, density_center, PhaseProportionsAir_Center, LocalSurfaceTopography,
														iel_x,iel_y,iel_z,zs_FreeSurface,dx_P, dy_P, dz_P, user, Eta_Center, Eta_XY, Eta_YZ, Eta_XZ, dRho_dxdydz, FreeSurfaceCells, &z_FreeSurface); CHKERRQ(ierr);

				// debugging
				/*
					Eta_Center[0]=Eta_Center[1]=Eta_Center[2]=Eta_Center[3]=Eta_Center[4]=Eta_Center[5]=Eta_Center[6]=1;
				//Eta_XZ[0]=Eta_XZ[1]=Eta_XZ[2]=Eta_XZ[3]=1;
				*/
				//Eta_XY[0]=Eta_XY[1]=Eta_XY[2]=Eta_XY[3]=1;
				//Eta_YZ[0]=Eta_YZ[1]=Eta_YZ[2]=Eta_YZ[3]=1;

				//PetscPrintf(PETSC_COMM_WORLD,"#   FDSTAG i,j,k=[%lld,%lld,%lld], Eta_XY=[%f,%f,%f,%f], Eta_YZ=[%f,%f,%f,%f] \n",	i,j,k,Eta_XY[0],Eta_XY[1],Eta_XY[2],Eta_XY[3],Eta_YZ[0],Eta_YZ[1],Eta_YZ[2],Eta_YZ[3]);



				/* If we have an internal free surface, compute the height of that relatve to the center pressure node */
				FreeSurface_Fraction=1.0;
				if ((FreeSurfaceCells[1]==0) && (FreeSurfaceCells[6]==1) && (eliminate_stickyair_from_system)){
					// Current cell is below FS & cell above is above FS
					FreeSurface_Fraction = (z_FreeSurface-z_center)/dz_P[1];

			//		PetscPrintf(PETSC_COMM_WORLD,"#   FDSTAG i,j,k=[%lld,%lld,%lld], z_center=%f, z_FS=%f, FreeSurface_Fraction=%f dz=%f FreeSurfaceCells=[%i,%i,%i,%i,%i,%i,%i] \n",	i,j,k,z_center, z_FreeSurface,FreeSurface_Fraction,dz_P[1],FreeSurfaceCells[0],FreeSurfaceCells[1],FreeSurfaceCells[2],FreeSurfaceCells[3],FreeSurfaceCells[4],FreeSurfaceCells[5],FreeSurfaceCells[6]);

				}


				/* */
				//PetscPrintf(PETSC_COMM_SELF,"#   FDSTAG i,j,k=[%lld,%lld,%lld] dRho_dxdydz=[%1.10e,%1.10e,%1.10e] z_FreeSurface=%f dx_P=[%f,%f]\n", (LLD)i,(LLD)j,(LLD)k, dRho_dxdydz[0]*user->Gravity*dt,dRho_dxdydz[1]*user->Gravity*dt,dRho_dxdydz[2]*user->Gravity*dt,z_FreeSurface,dx_P[0],dx_P[1]);


				/* The nomenclature for the global stiffness matrixes is:
				 *
				 * | VV_MAT VP_MAT | | V |    |Rhs |
				 * |        	   | |	 |  = |    |
				 * | PV_MAT PP_MAT | | P |    | 0  |
				 *
				 */

				/* (1) Create FD stencil for the incompressibility equation
				 *
				 * 		(Vx(i+1)-Vx(i))/dx + (Vy(j+1)-Vy(j))/dy + (Vz(k+1)-Vz(k))/dz = 1/gamma*P(i,j,k)
				 *
				 *	BC's are not imposed on this equation, so we don't have to worry about it.
				 */
				if ( (i < (nnode_x-1)) && (j < (nnode_y-1)) && (k < (nnode_z-1)) ) {
					ierr = FDSTAG_FDStencil_Incompressibility(i,j,k,dx_vec, dy_vec, dz_vec, row_pp,v_pp,col_pv, v_pv); CHKERRQ(ierr);

					/* Find global indices */
					ierr = DAGetGlobalIndex(da_pres,row_pp[0].i,row_pp[0].j,row_pp[0].k,row_pp[0].c,&rowidx); CHKERRQ(ierr);
					ierr = FDSTAG_GetStencilValues(da,da_pres,6,col_pv,0,PETSC_NULL,vgidx,pgidx); CHKERRQ(ierr);

					/* Add pressure BC's for the free surface */
					if ((eliminate_stickyair_from_system) && (FreeSurfaceCells[1]==1)){
						v_pv[0] = v_pv[1] = v_pv[2] = v_pv[3] = v_pv[4] = v_pv[5] = 0;		// add a constant pressure BC here
					}


					ierr = MatSetValues(PV_MAT,1,&rowidx,6,vgidx,v_pv,ADD_VALUES); CHKERRQ(ierr);

					// Pushing BC - create copy matrix for PV_MAT
					if (user->AddPushing == 1){
						ierr = MatSetValues(PV_MAT_push,1,&rowidx,6,vgidx,v_pv,ADD_VALUES); CHKERRQ(ierr);
					}

					/* Create the approximate Schur matrix, with 1/eta on the diagonal */
					v_pp[0] = -1.0/Eta_Center[1];
					ierr = MatSetValue(approx_S,rowidx,rowidx,v_pp[0],INSERT_VALUES); CHKERRQ(ierr);

					/* Store vector with viscosity scaling */
					ierr = VecSetValue(ViscosityScaling,rowidx,1/Eta_Center[1],INSERT_VALUES); CHKERRQ(ierr);


					/* Create the PP matrix, with 0 on the diagonal */
					if (user->StokesSolver==1){
						v_pp[0] = -1;		// powell-hesteness scales it later
					}
					else {
						v_pp[0] = 0;		// as the system is incompressible
					}

					/* Add pressure BC's for the free surface */
					if ((eliminate_stickyair_from_system) && (FreeSurfaceCells[1]==1)){
						v_pp[0] = 1.0;		// add a constant pressure BC for sticky-air cells
					}
					ierr = MatSetValue(PP_MAT,rowidx,rowidx,v_pp[0],INSERT_VALUES); CHKERRQ(ierr);

				}

				/* (2)  Create FD stencil for the first force balance equation:
				 * -(P(i+1)-P(i))/dx + (Txx(i+1)-Txx(i))/dx + (Txy(j+1)-Txy(j))/dy + (Txz(k+1)-Txz(k))/dz = 0
				 *
				 *	This equation is located @ Vx nodes.
				 */
				if ((i > 0) && (i < (nnode_x-1)) ) {
					/* Compute FD stencil (w/out BCs) */
					ierr = FDSTAG_FDStencil_ForceBalance1(i,j,k,dx_vec,dy_vec, dz_vec, dx_P,dy_P,dz_P, Eta_Center, Eta_XY, Eta_XZ, row_vv, col_vv, v_vv, col_vp,  v_vp); CHKERRQ(ierr);

					/* Correct stencil for BCs (front/back/lower/upper) */
					ierr = FDSTAG_SetBCs_FDStencil_ForceBalance1(j,k,nel_y,nel_z, v_vv, user); CHKERRQ(ierr);

					/* Correct stencil for internal free surface */
					ierr = DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx); CHKERRQ(ierr);
					ierr = FDSTAG_GetStencilValues(da,da_pres,number_coefficients,col_vv,2,col_vp,vgidx,pgidx); CHKERRQ(ierr);

					ierr = MatSetValues(VV_MAT,1,&rowidx,number_coefficients,vgidx,v_vv,ADD_VALUES); CHKERRQ(ierr);
					ierr = MatSetValues(VP_MAT,1,&rowidx,2, pgidx,v_vp,ADD_VALUES); CHKERRQ(ierr);

				}

				/* (3)  Create FD stencil for the second force balance equation:
				 * -(P(j+1)-P(j))/dy + (Txy(i+1)-Txy(i))/dx + (Tyy(j+1)-Tyy(j))/dy + (Tyz(k+1)-Tyz(k))/dz = 0
				 *
				 * This equation is located @ Vy nodes.
				 */
				if ((j > 0) && (j < (nnode_y-1))) {
					FDSTAG_FDStencil_ForceBalance2(i,j,k,dx_vec,dy_vec, dz_vec, dx_P,dy_P,dz_P, Eta_Center, Eta_XY, Eta_YZ,	row_vv, col_vv, v_vv, col_vp,  v_vp);

					/* Correct stencil for BCs (left/right/lower/upper) */
					ierr = FDSTAG_SetBCs_FDStencil_ForceBalance2(i,k,nel_x,nel_z, v_vv, user); CHKERRQ(ierr);

					/* Global indices */
					ierr = DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx); CHKERRQ(ierr);
					ierr = FDSTAG_GetStencilValues(da,da_pres,number_coefficients,col_vv,2,col_vp,vgidx,pgidx); CHKERRQ(ierr);

					ierr = MatSetValues(VV_MAT,1,&rowidx,number_coefficients,vgidx,v_vv,ADD_VALUES); CHKERRQ(ierr);
					ierr = MatSetValues(VP_MAT,1,&rowidx,2, pgidx,v_vp,ADD_VALUES); CHKERRQ(ierr);

				}

				/* (4)  Create FD stencil for the third force balance equation:
				 * -(P(k)-P(k-1))/dz + (Txz(i+1)-Txz(i))/dx + (Syz(j+1)-Syz(j))/dy + (Tzz(k+1)-Tzz(k))/dz = 0
				 *
				 * This equation is located @ Vz nodes
				 */
				if ( (k > 0) && (k < (nnode_z-1)) ) {

					FreeSurface_Fraction=1.0;
					FDSTAG_FDStencil_ForceBalance3(i,j,k,dx_vec,dy_vec, dz_vec, dx_P,dy_P,dz_P, Eta_Center, Eta_XZ, Eta_YZ, row_vv, col_vv, v_vv, col_vp,  v_vp, FreeSurface_Fraction);

					/* Correct stencil for BCs (left/right/front/back) */
					ierr = FDSTAG_SetBCs_FDStencil_ForceBalance3(i,j,nel_x,nel_y, v_vv, v_vp, FreeSurfaceCells, user, eliminate_stickyair_from_system); CHKERRQ(ierr);

					/* Find global indices */
					ierr = DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx); CHKERRQ(ierr);
					ierr = FDSTAG_GetStencilValues(da,da_pres,number_coefficients,col_vv,2,col_vp,vgidx,pgidx); CHKERRQ(ierr);

					ierr = MatSetValues(VV_MAT,1,&rowidx,number_coefficients,vgidx,v_vv,ADD_VALUES); CHKERRQ(ierr);
					ierr = MatSetValues(VP_MAT,1,&rowidx,2, pgidx,v_vp,ADD_VALUES); CHKERRQ(ierr);

					/* Add FSSA algorithm as described in Kaus et al (2010) and Popov and Sobolev (2008). The algorithm described in Duretz et al., G-cubed (2011) introduces asymmetry terms */
					/* There is still the discussion whether the correction is negative or positive. Here it is: sigma_ij,j + FSSA*dRho/dz*g_i = rho*g_i with g_i = [0 0 -g] */
					v_vv[0] = user->FSSA*dRho_dxdydz[2]*user->Gravity*dt;
					ierr = MatSetValue(VV_MAT,rowidx,rowidx,v_vv[0],ADD_VALUES); CHKERRQ(ierr);		// add FSSA dRho_dz (Vz) contribution

				}

			}
		}
	}
	PetscTime(&t1);
	PetscPrintf(PETSC_COMM_WORLD,"#  FDSTAG: Time to insert stencil into operators %1.4e \n", t1-t0);





#if 1
		ScalingParameterDiagonal 		=	 1.0;		// The value set on the diagonal; it should be on

		// Allocate array that will contain indices of the boundary conditions
		PetscMalloc((size_t)(2*xmp*ymp + 2*ymp*zmp + 2*xmp*zmp)*sizeof(PetscInt),&rowid_DirichletBC_array);

		/* Set Dirichlet BCs for velocity normal to the side boundaries
		 *
		 * Note:  we do set "dummy" values here as PETSC requires the full matrix to be assembled in order to symmetrize it.
		 * Yet, these values will be overwritten at a later stage once we symmetrize the matrix
		 * */
		numDirichlet=0;
		PetscTime(&t0);
		for (k=zs; k<zs+zm; k++) {
			for (j=ys; j<ys+ym; j++) {
				for (i=xs; i<xs+xm; i++) {

					if ( ((i == 0) || (i == nnode_x-1)) && (j < nnode_y-1) && (k < nnode_z-1) ) {
						/* Left or right boundary [Vx=constant] */
						row_vv[0].i = i; row_vv[0].j = j; row_vv[0].k = k; row_vv[0].c = 0;
						v_vv[0] = ScalingParameterDiagonal;
						ierr = DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx); CHKERRQ(ierr);
						ierr = MatSetValue(VV_MAT,rowidx,rowidx,v_vv[0],ADD_VALUES); CHKERRQ(ierr);

						rowid_DirichletBC_array[numDirichlet]	= 	rowidx;
						numDirichlet 							=	numDirichlet + 1;

#ifdef  LAMEM_DEBUG_FDSTAG
					//	ierr = StencilToLocalNumbering(VV_MAT,1,row_vv, row_vv_index); CHKERRQ(ierr);
					//	PetscPrintf(PETSC_COMM_WORLD,"#   FDSTAG BC Assembly:i,j,k=[%lld,%lld,%lld] c=%lld, row=%lld\n", (LLD)i,(LLD)j,(LLD)k, (LLD)(row_vv[0].c), (LLD)row_vv_index[0]);
						PetscPrintf(PETSC_COMM_WORLD,"#   FDSTAG BC Assembly Vx:i,j,k=[%lld,%lld,%lld] c=%lld, row=%lld\n", (LLD)i,(LLD)j,(LLD)k, (LLD)(row_vv[0].c), (LLD)rowidx);
#endif
					}

					if ( ((j == 0) || (j == nnode_y-1)) && (i < nnode_x-1) && (k < nnode_z-1)) {
						/* Front or back boundary [Vy=constant]*/
						row_vv[0].i = i; row_vv[0].j = j; row_vv[0].k = k; row_vv[0].c = 1;
						v_vv[0] = ScalingParameterDiagonal;
						ierr = DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx); CHKERRQ(ierr);
						ierr = MatSetValue(VV_MAT,rowidx,rowidx,v_vv[0],ADD_VALUES); CHKERRQ(ierr);
						//ierr = MatZeroRowsColumns(VV_MAT,1,&rowidx,v_vv[0],PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);

						rowid_DirichletBC_array[numDirichlet]	= 	rowidx;
						numDirichlet 							=	numDirichlet + 1;

	#ifdef  LAMEM_DEBUG_FDSTAG
					//	ierr = StencilToLocalNumbering(VV_MAT,1,row_vv, row_vv_index); CHKERRQ(ierr);
					//	PetscPrintf(PETSC_COMM_WORLD,"#   FDSTAG BC Assembly:i,j,k=[%lld,%lld,%lld] c=%lld, row=%lld\n", (LLD)i,(LLD)j,(LLD)k, (LLD)(row_vv[0].c), (LLD)row_vv_index[0]);
						PetscPrintf(PETSC_COMM_WORLD,"#   FDSTAG BC Assembly Vy:i,j,k=[%lld,%lld,%lld] c=%lld, row=%lld\n", (LLD)i,(LLD)j,(LLD)k, (LLD)(row_vv[0].c), (LLD)rowidx);

	#endif
					}

					if ( (k == nnode_z-1) && (i < nnode_x-1) && (j < nnode_y-1) ) {
						/* Upper boundary */

						if (user->BC.UpperBound == 1) {
							//  Free slip: Vz=0
							row_vv[0].i = i; row_vv[0].j = j; row_vv[0].k = k; row_vv[0].c = 2;
							v_vv[0] = ScalingParameterDiagonal;
							ierr = DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx); CHKERRQ(ierr);
							ierr = MatSetValue(VV_MAT,rowidx,rowidx,v_vv[0],ADD_VALUES); CHKERRQ(ierr);
							//ierr = MatZeroRowsColumns(VV_MAT,1,&rowidx,v_vv[0],PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);

							rowid_DirichletBC_array[numDirichlet]	= 	rowidx;
							numDirichlet 							=	numDirichlet + 1;

						} else if (user->BC.UpperBound ==0 ) {
							// stress-free; Szz=0

							if ( (j == 0) || (j == nnode_y-1) || (i == 0) || (i == nnode_x-1) ) {
								// nodes at the side boundaries (not sure whether this is necessary)
								row_vv[0].i = i; row_vv[0].j = j; row_vv[0].k = k; row_vv[0].c = 2;
								v_vv[0] = ScalingParameterDiagonal;
								//					ierr = MatSetValuesStencil(VV_MAT,1,row_vv,1,row_vv,v_vv,ADD_VALUES);		CHKERRQ(ierr);
								ierr = DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx); CHKERRQ(ierr);
								ierr = MatSetValue(VV_MAT,rowidx,rowidx,v_vv[0],ADD_VALUES); CHKERRQ(ierr);

								rowid_DirichletBC_array[numDirichlet]	= 	rowidx;
								numDirichlet 							=	numDirichlet + 1;

							} else {
								// center nodes

								FreeSurface_Fraction=1.0;
								FDSTAG_FDStencil_ForceBalance3(i,j,k,dx_vec,dy_vec, dz_vec, dx_P,dy_P,dz_P, Eta_Center, Eta_XZ, Eta_YZ, row_vv, col_vv, v_vv, col_vp,  v_vp, FreeSurface_Fraction);

								/* Correct stencil for Stress-free upper BC */
								// delete coeffs by hand
								v_vp[1] 	= 0;				//	P(iz+1/2)=0
								v_vv[0] 	= 0;				// 	tau_zz(iz+1/2)=0
								v_vv[1] 	= 0;
								col_vv[0] = col_vv[1];		// 	Put 'fake' coefficients as the routines below otherwise complain that it is outside the domain
								col_vp[1] = col_vp[0];

								/* Global indices */
								ierr = DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx); CHKERRQ(ierr);
								ierr = FDSTAG_GetStencilValues(da,da_pres,20,col_vv,2,col_vp,vgidx,pgidx); CHKERRQ(ierr);

								ierr = MatSetValues(VV_MAT,1,&rowidx,20,vgidx,v_vv,ADD_VALUES); CHKERRQ(ierr);
								ierr = MatSetValues(VP_MAT,1,&rowidx,2, pgidx,v_vp,ADD_VALUES); CHKERRQ(ierr);

								rowid_DirichletBC_array[numDirichlet]	= 	rowidx;
								numDirichlet 							=	numDirichlet + 1;

							}


						}
					} else if ( (k == 0) && (i < nnode_x-1) && (j < nnode_y-1) ) {
						/* Lower or upper boundary [Vz=constant]*/
						row_vv[0].i = i; row_vv[0].j = j; row_vv[0].k = k; row_vv[0].c = 2;
						v_vv[0] 								= 	ScalingParameterDiagonal;
						ierr 									= 	DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx); CHKERRQ(ierr);
						ierr 									=	MatSetValue(VV_MAT,rowidx,rowidx,v_vv[0],ADD_VALUES); CHKERRQ(ierr);
						rowid_DirichletBC_array[numDirichlet]	= 	rowidx;
						numDirichlet 							=	numDirichlet + 1;


	#ifdef  LAMEM_DEBUG_FDSTAG
					//	ierr = StencilToLocalNumbering(VV_MAT,1,row_vv, row_vv_index); CHKERRQ(ierr);
					//	PetscPrintf(PETSC_COMM_WORLD,"#   FDSTAG BC Assembly:i,j,k=[%lld,%lld,%lld] c=%lld, row=%lld\n", (LLD)i,(LLD)j,(LLD)k, (LLD)(row_vv[0].c), (LLD)row_vv_index[0]);
						PetscPrintf(PETSC_COMM_WORLD,"#   FDSTAG BC Assembly Vz:i,j,k=[%lld,%lld,%lld] c=%lld, row=%lld\n", (LLD)i,(LLD)j,(LLD)k, (LLD)(row_vv[0].c), (LLD)rowidx);

	#endif

					}

				}
			}
		}


#if 1
	/* The FDSTAG formulation, in combination with the DA, results in some velocity points that are not used in the discretization.
	 * Yet, the matrix should not contain any rows that are zero.
	 * Therefore, we have to set the diagonal of the matrix to 1 for these velocity points */
	ierr = DMDAGetInfo(da, 0, &nnode_x, &nnode_y,	&nnode_z,	0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);  // # of nodes 	 in all directions
	ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm); CHKERRQ(ierr);

	// Allocate array that will contain indices of the boundary conditions
	PetscMalloc((size_t)(2*xm*ym + 2*ym*zm + 2*xm*zm)*sizeof(PetscInt),&rowid_addPoints_array);

	numAddpoints = 0;

	/* Set additional Vx points */
	for (j=ys; j<ys+ym; j++) {
		for	(i=xs; i<xs+xm; i++) {
			if (zs+zm == nnode_z) { // the local processor borders this boundary
				/* only if the proc borders the upper boundary */
				k = nnode_z-1;
				row_vv[0].i = i; row_vv[0].j = j; row_vv[0].k = k; row_vv[0].c = 0;
				v_vv[0] = ScalingParameterDiagonal;

				ierr = DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx); CHKERRQ(ierr);
				ierr = MatSetValue(VV_MAT,rowidx,rowidx,v_vv[0],ADD_VALUES); CHKERRQ(ierr);
				//ierr = MatZeroRowsColumns(VV_MAT,1,&rowidx,v_vv[0],PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);

				rowid_addPoints_array[numAddpoints]	= 	rowidx;
				numAddpoints 						=	numAddpoints + 1;

			}
		}
	}

	for (k=zs; k<zs+zm; k++) {
		for (i=xs; i<xs+xm; i++) {
			if (ys+ym == nnode_y) {
				/* only if the proc borders the back boundary */
				j = nnode_y-1;
				row_vv[0].i = i; row_vv[0].j = j; row_vv[0].k = k; row_vv[0].c = 0;
				v_vv[0] = ScalingParameterDiagonal;

				ierr = DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx); CHKERRQ(ierr);
				ierr = MatSetValue(VV_MAT,rowidx,rowidx,v_vv[0],ADD_VALUES); CHKERRQ(ierr);
				//ierr = MatZeroRowsColumns(VV_MAT,1,&rowidx,v_vv[0],PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);

				rowid_addPoints_array[numAddpoints]	= 	rowidx;
				numAddpoints 						=	numAddpoints + 1;
			}
		}
	}

	/* Set additional Vy points */
	for (j=ys; j<ys+ym; j++) {
		for (i=xs; i<xs+xm; i++) {
			if (zs+zm == nnode_z) {
				/* only if the proc borders the upper boundary */
				k = nnode_z-1;
				row_vv[0].i = i; row_vv[0].j = j; row_vv[0].k = k; row_vv[0].c = 1;
				v_vv[0] = 1.0;

				ierr = DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx); CHKERRQ(ierr);
				ierr = MatSetValue(VV_MAT,rowidx,rowidx,v_vv[0],ADD_VALUES); CHKERRQ(ierr);
				//ierr = MatZeroRowsColumns(VV_MAT,1,rowidx,v_vv[0],PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);

				rowid_addPoints_array[numAddpoints]	= 	rowidx;
				numAddpoints 						=	numAddpoints + 1;
			}
		}
	}

	for (k=zs; k<zs+zm; k++) {
		for (j=ys; j<ys+ym; j++) {
			if (xs+xm == nnode_x) {
				/* only if the proc borders the right boundary */
				i = nnode_x-1;
				row_vv[0].i = i; row_vv[0].j = j; row_vv[0].k = k; row_vv[0].c = 1;
				v_vv[0] = ScalingParameterDiagonal;

				ierr = DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx); CHKERRQ(ierr);
				ierr = MatSetValue(VV_MAT,rowidx,rowidx,v_vv[0],ADD_VALUES); CHKERRQ(ierr);
				//ierr = MatZeroRowsColumns(VV_MAT,1,&rowidx,v_vv[0],PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);

				rowid_addPoints_array[numAddpoints]	= 	rowidx;
				numAddpoints 						=	numAddpoints + 1;
			}
		}
	}

	/* Set additional Vz points */
	for (k=zs; k<zs+zm; k++) {
		for (i=xs; i<xs+xm; i++) {
			if (ys+ym == nnode_y) {
				/* only if the proc borders the back boundary */
				j = nnode_y-1;
				row_vv[0].i = i; row_vv[0].j = j; row_vv[0].k = k; row_vv[0].c = 2;
				v_vv[0] = ScalingParameterDiagonal;

				ierr = DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx); CHKERRQ(ierr);
				ierr = MatSetValue(VV_MAT,rowidx,rowidx,v_vv[0],ADD_VALUES); CHKERRQ(ierr);
				//ierr = MatZeroRowsColumns(VV_MAT,1,&rowidx,v_vv[0],PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);

				rowid_addPoints_array[numAddpoints]	= 	rowidx;
				numAddpoints 						=	numAddpoints + 1;
			}

		}
	}

	for (k=zs; k<zs+zm; k++) {
		for (j=ys; j<ys+ym; j++) {
			if (xs+xm == nnode_x) {
				/* only if the proc borders the right boundary */
				i = nnode_x-1;
				row_vv[0].i = i; row_vv[0].j = j; row_vv[0].k = k; row_vv[0].c = 2;
				v_vv[0] = ScalingParameterDiagonal;

				ierr = DAGetGlobalIndex(da,row_vv[0].i,row_vv[0].j,row_vv[0].k,row_vv[0].c,&rowidx);CHKERRQ(ierr);
				ierr = MatSetValue(VV_MAT,rowidx,rowidx,v_vv[0],ADD_VALUES); CHKERRQ(ierr);
				//ierr = MatZeroRowsColumns(VV_MAT,1,rowidx,v_vv[0],PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);

				rowid_addPoints_array[numAddpoints]	= 	rowidx;
				numAddpoints 						=	numAddpoints + 1;
			}
		}
	}
	PetscTime(&t1);
	PetscPrintf(PETSC_COMM_WORLD,"#  FDSTAG: Time to insert bc's into operators %1.4e \n", t1-t0);

#endif

#endif // setting Dirichlet BC

	ierr = DMDAVecRestoreArray(user->FDSTAG.DA_CENTER,      Viscosity_Center_local, &viscosity_center); 	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(user->FDSTAG.DA_CENTER,      Density_Center_local,   &density_center); 	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(user->FDSTAG.DA_XY_POINTS,   Viscosity_XY_local,     &viscosity_XY); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(user->FDSTAG.DA_XZ_POINTS,   Viscosity_XZ_local,     &viscosity_XZ); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(user->FDSTAG.DA_YZ_POINTS,   Viscosity_YZ_local,     &viscosity_YZ); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(user->FDSTAG.DA_CENTER,  		PhaseProportionsAir_Center_Vec, 	&PhaseProportionsAir_Center	 ); CHKERRQ(ierr);


	ierr = DMRestoreLocalVector(user->FDSTAG.DA_CENTER,		&Viscosity_Center_local);	CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(user->FDSTAG.DA_CENTER,		&Density_Center_local);	CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(user->FDSTAG.DA_XY_POINTS	,&Viscosity_XY_local); CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(user->FDSTAG.DA_XZ_POINTS,	&Viscosity_XZ_local); CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(user->FDSTAG.DA_YZ_POINTS,	&Viscosity_YZ_local); CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(user->FDSTAG.DA_CENTER,		&PhaseProportionsAir_Center_Vec);	CHKERRQ(ierr);

	// Get information and modify PV_MAT for pushing BC - this needs to be done before the MatAssembly routine for PV_MAT
	ierr = MatAssemblyBegin(PV_MAT_push,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(PV_MAT_push,MAT_FINAL_ASSEMBLY); 	 CHKERRQ(ierr);

	ierr = MatAssemblyBegin(PV_MAT,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(PV_MAT,MAT_FLUSH_ASSEMBLY);   CHKERRQ(ierr);

	ierr = AddPushingToModelPV_MAT(user,PV_MAT,PV_MAT_push,pv_rhs_push); CHKERRQ(ierr);
	ierr = MatDestroy(&PV_MAT_push); CHKERRQ(ierr);

	/* Finalize the matrixes */
	PetscTime(&t0);
	ierr = MatAssemblyBegin(VV_MAT,	MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);	// Assemble global velocity stiffness matrix
	ierr = MatAssemblyEnd  (VV_MAT,	MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);  // Assemble global velocity stiffness matrix

	ierr = MatAssemblyBegin(PP_MAT,	MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);	// Assemble global pressure mass matrix
	ierr = MatAssemblyEnd  (PP_MAT,	MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);  // Assemble global pressure mass matrix

	ierr = MatAssemblyBegin(VP_MAT,	MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);	// Assemble global VP matrix
	ierr = MatAssemblyEnd  (VP_MAT,	MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);  // Assemble global VP matrix

	ierr = MatAssemblyBegin(PV_MAT,	MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);	// Assemble global PV matrix
	ierr = MatAssemblyEnd  (PV_MAT,	MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);  // Assemble global PV matrix
	PetscTime(&t1);
	PetscPrintf(PETSC_COMM_WORLD,"#  FDSTAG: Time to assemble operators %1.4e \n", t1-t0);

	ierr = DMRestoreLocalVector(da,&rhs_add_BC_loc); CHKERRQ(ierr);

	ierr = MatAssemblyBegin(approx_S,	MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);	// Assemble schur preconditioning matrix
	ierr = MatAssemblyEnd  (approx_S, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);  // Assemble schur preconditioning matrix

	ierr = VecAssemblyBegin(ViscosityScaling); CHKERRQ(ierr);	//Assemble ViscosityScaling
	ierr = VecAssemblyEnd  (ViscosityScaling); CHKERRQ(ierr);  //Assemble ViscosityScaling

	//coordinates
	ierr = DMDAVecRestoreArray(cda,gc,&coords); CHKERRQ(ierr);
	//ierr = DMDestroy(cda); CHKERRQ(ierr);
	//ierr = VecDestroy(gc); CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(da,rhs_add_BC_loc,&rhs_for_BoundaryConditions); CHKERRQ(ierr);

	/* Cleaning up */
	ierr 	= 	DMDAVecRestoreArray(cda_SurfaceTopo,gc_SurfaceTopo,	&coors_SurfaceTopo); 							CHKERRQ(ierr);
	ierr 	= 	DMDAVecRestoreArray(user->DA_SurfaceTopography,LocalSurfaceTopography_vec, 	&LocalSurfaceTopography	);	CHKERRQ(ierr);
	ierr 	= 	DMRestoreLocalVector(user->DA_SurfaceTopography,	&LocalSurfaceTopography_vec); CHKERRQ(ierr);




	/* 				Make the VV and VP matrix symmetric
	 *
	 * 	Doing this requires us to know on which rows we have set BC's; we have stored that above for
	 * 		rowid_DirichletBC_array - Dirichlet BC's  (essentially normal velocities) to the boundary
	 * 		rowid_addPoints_array 	- Additional points added since we use a DMDA with 3 dof
	 *
	 * */
	IS				is;
	Vec				DirichletVelocityValues;

	// Create temporary vector in which we will set the correct boundary velocities
	ierr = VecDuplicate(rhs_add_BC, &DirichletVelocityValues); 								CHKERRQ(ierr);
	ierr = VecZeroEntries(DirichletVelocityValues);				 							CHKERRQ(ierr);
	ierr = VecZeroEntries(rhs_add_BC);				 										CHKERRQ(ierr);
	ierr = SetBoundaryConditionsRHS_FDSTAG(user->DA_Vel, user, DirichletVelocityValues, 1); CHKERRQ(ierr);

	// Symmetrize Dirichlet BC's - in case we use non-zero velocities at the boundaries, the rhs-contributions will be added to rhs_add_BC
	ierr = ISCreateGeneral(PETSC_COMM_WORLD,numDirichlet,rowid_DirichletBC_array,PETSC_COPY_VALUES,&is);  CHKERRQ(ierr);
	ierr = MatZeroRowsColumnsIS(VV_MAT,is,1.0,DirichletVelocityValues,rhs_add_BC);  		CHKERRQ(ierr); 				// make VV matrix symmetric, and add additional terms to RHS
	ierr = MatZeroRowsIS(VP_MAT,is,0.0,PETSC_NULL,PETSC_NULL);  							CHKERRQ(ierr); 				// delete corresponding rows in VP matrix
	ierr = ISDestroy(&is); 																	CHKERRQ(ierr);

	// Symmetrize additional Points (as they are only added since we use a DMDA to create a grid for the FDSTAG discretization)
	ierr = ISCreateGeneral(PETSC_COMM_WORLD,numAddpoints,rowid_addPoints_array,PETSC_COPY_VALUES,&is);  CHKERRQ(ierr);
	ierr = MatZeroRowsColumnsIS(VV_MAT,is,1.0,PETSC_NULL,PETSC_NULL);  						CHKERRQ(ierr);
	ierr = MatZeroRowsIS(VP_MAT,is,0.0,PETSC_NULL,PETSC_NULL);  							CHKERRQ(ierr);							// delete corresponding rows in VP matrix
	ierr = ISDestroy(&is); 																	CHKERRQ(ierr);

	/* Clean up */
	ierr = PetscFree(rowid_DirichletBC_array); 												CHKERRQ(ierr);
	ierr = PetscFree(rowid_addPoints_array);												CHKERRQ(ierr);
	ierr = VecDestroy(&DirichletVelocityValues);											CHKERRQ(ierr);







	PetscTime(&cputime_end);

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Computes the staggered FD stencil for the incompressibility equation  */
#undef __FUNCT__
#define __FUNCT__ "FDSTAG_FDStencil_Incompressibility"
PetscErrorCode FDSTAG_FDStencil_Incompressibility(PetscInt i, PetscInt j, PetscInt k,
		PetscScalar dx_vec[7], PetscScalar dy_vec[7], PetscScalar dz_vec[7],
		MatStencil  row_pp[], 	PetscScalar v_pp[],
		MatStencil 	col_pv[], 	PetscScalar v_pv[])
{
	PetscScalar		dx,dy,dz;
	PetscScalar		v;
	MatStencil		row, col;

	PetscFunctionBegin;

	// Equation  1/gamma*P(k,j,i)
	row.i = i; row.j = j; row.k = k; row.c = 0;		v  = 1.0/1e3;		row_pp[0] 	= 	row;	v_pp[0] = v;

	dx 			=	dx_vec[1];
	dy 			=	dy_vec[1];
	dz 			=	dz_vec[1];

	// (Vx(i+1)-Vx(i))/dx
	col.i = 	i+1; col.j = j; col.k = k; col.c = 0.0;	v =  1.0/dx;	col_pv[0] = col;	v_pv[0] = v;	// 	Vx(i+1)/dx
	col.i = 	i  ; col.j = j; col.k = k; col.c = 0.0;	v = -1.0/dx; 	col_pv[1] = col;	v_pv[1] = v; 	// 	-Vx(i)/dx

	// (Vy(j+1)-Vy(j))/dy
	col.i = 	i; col.j = j+1; col.k = k; col.c = 1.0;	v =  1.0/dy;	col_pv[2] = col;	v_pv[2] = v;	//   Vy(i+1)/dy
	col.i = 	i; col.j = j  ; col.k = k; col.c = 1.0;	v = -1.0/dy;	col_pv[3] = col;	v_pv[3] = v;	//  -Vy(i)/dy


	// (Vz(k+1)-Vz(k))/dz
	col.i = 	i; col.j = j; col.k = k+1; col.c = 2.0;	v =  1.0/dz;	col_pv[4] = col;	v_pv[4] = v;	//   Vz(i+1)/dz
	col.i = 	i; col.j = j; col.k = k  ; col.c = 2.0;	v = -1.0/dz;	col_pv[5] = col;	v_pv[5] = v;	//  -Vz(i)/dz

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Computes the staggered FD stencil for the 1th force balance equation, which is located at Vx nodes.
 *
 */
#undef __FUNCT__
#define __FUNCT__ "FDSTAG_FDStencil_ForceBalance1"
PetscErrorCode FDSTAG_FDStencil_ForceBalance1(PetscInt i, PetscInt j, PetscInt k,
		PetscScalar dx_vec[7], 		PetscScalar dy_vec[7], 		PetscScalar dz_vec[7],
		PetscScalar   dx_P[2], 		PetscScalar	dy_P[2], 		PetscScalar   dz_P[2],
		PetscScalar	Eta_Center[7], 	PetscScalar Eta_XY[4],		PetscScalar Eta_XZ[4],
		MatStencil  row_vv[], 		MatStencil  col_vv[], 		PetscScalar	v_vv[],
		MatStencil  col_vp[],  		PetscScalar v_vp[])
{
	PetscScalar		dx_P_west,dx_west,dx_cen, dy_P_south, dy_P_north, dy_cen, dz_cen, dz_P_top, dz_P_bottom;
	PetscScalar		mu_west, mu_cen, mu_northwest, mu_southwest, mu_bottomwest, mu_topwest;
	PetscScalar		v;
	MatStencil		row, col;

	PetscFunctionBegin;

	/* Extract info needed to construct the FD stencil */
	dx_P_west 	=	dx_P[0];
	dx_west 	=	dx_vec[0];
	dx_cen 		=	dx_vec[1];
	dy_cen 		=	dy_vec[1];
	dz_cen 		=	dz_vec[1];

	dy_P_north 	=	dy_P[1];
	dy_P_south 	=	dy_P[0];

	dz_P_top	= 	dz_P[1];
	dz_P_bottom	= 	dz_P[0];


	mu_west 	=	Eta_Center[0];
	mu_cen 		=	Eta_Center[1];
	mu_northwest=   Eta_XY[2];
	mu_southwest=   Eta_XY[0];

	mu_topwest 		=	Eta_XZ[1];
	mu_bottomwest 	=	Eta_XZ[0];



	/* Equation number [1th Force balance is located @ Vx nodes] */
	row.i = i;   row.j = j; row.k = k; row.c = 0;	row_vv[0] = row;		// 1th force balance equation is located @ Vx(i) nodes

	// -[ P(i) - P(i-1) ]/dx_P_west		- Pressure equations
	col.i = i-1; col.j = j; col.k = k; col.c = 0; v= 1/dx_P_west;	col_vp[0] = col;	v_vp[0] = v; 	//  P(i-1)/dx
	col.i = i  ; col.j = j; col.k = k; col.c = 0; v=-1/dx_P_west;	col_vp[1] = col;	v_vp[1] = v;		// -P(i  )/dx

	// [Txx(i)-Txx(i-1)]/dx
	col.i = 	i+1; col.j = j; col.k = k; col.c = 0;	v =    2.0*mu_cen/(dx_cen*dx_P_west);								col_vv[0] = col; v_vv[0] = v;			//  Vx(i+1,j,k)
	col.i = 	i  ; col.j = j; col.k = k; col.c = 0;	v =  (-2.0*mu_cen/dx_cen)/dx_P_west;								col_vv[1] = col; v_vv[1] = v;			//  Vx(i  ,j,k)
	col.i = 	i  ; col.j = j; col.k = k; col.c = 0;	v =  (-2.0*mu_west/dx_west)/dx_P_west;								col_vv[2] = col; v_vv[2] = v;			//  Vx(i  ,j,k)
	col.i = 	i-1; col.j = j; col.k = k; col.c = 0;	v =    2.0*mu_west/(dx_west*dx_P_west);								col_vv[3] = col; v_vv[3] = v;			//  Vx(i-1,j,k)

	// [Sxy(j)-Sxy(j-1)]/dy
	// Sxy(j)/dy
	col.i = 	i; col.j = j+1; col.k = k; col.c = 0;	v =  mu_northwest/(dy_P_north*dy_cen);								col_vv[4] = col; v_vv[4] = v;			//  Vx(i,j+1)
	col.i = 	i; col.j = j;   col.k = k; col.c = 0;	v =  (-mu_northwest/dy_P_north)/dy_cen;								col_vv[5] = col; v_vv[5] = v;			//  Vx(i,j  )
	col.i =   i  ; col.j = j+1; col.k = k; col.c = 1;	v =  mu_northwest/(dx_P_west*dy_cen);								col_vv[6] = col; v_vv[6] = v;			//  Vy(i  ,j+1)
	col.i =	  i-1; col.j = j+1; col.k = k; col.c = 1;	v = -mu_northwest/(dx_P_west*dy_cen);								col_vv[7] = col; v_vv[7] = v;			//  Vy(i-1,j+1)

	// -Sxy(j-1)/dy
	col.i = 	i; col.j = j;   col.k = k; col.c = 0;	v =  (-mu_southwest/dy_P_south)/dy_cen;								col_vv[8] = col; v_vv[8] = v;			//  Vx(i,j  )
	col.i = 	i; col.j = j-1; col.k = k; col.c = 0;	v =  mu_southwest/(dy_P_south*dy_cen);								col_vv[9] = col; v_vv[9] = v;			//  Vx(i,j-1)
	col.i =   i  ; col.j = j  ; col.k = k; col.c = 1;	v = -mu_southwest/(dx_P_west*dy_cen);								col_vv[10]= col; v_vv[10] = v;			//  Vy(i  ,j  )
	col.i =   i-1; col.j = j  ; col.k = k; col.c = 1;	v =  mu_southwest/(dx_P_west*dy_cen);								col_vv[11]= col; v_vv[11]= v;			//  Vy(i-1,j  )


	// [Sxz(k)-Sxz(k-1)]/dz
	// Sxz(k)/dz
	col.i = 	i; col.j = j; col.k = k+1; col.c = 0;	v =  mu_topwest/(dz_P_top*dz_cen);									col_vv[12] = col; v_vv[12] = v;	//  Vx(i,k+1)
	col.i = 	i; col.j = j; col.k = k;   col.c = 0;	v =  (    -mu_topwest/dz_P_top  )/dz_cen;							col_vv[13] = col; v_vv[13] = v;	//  Vx(i,k  )
	col.i =   i  ; col.j = j; col.k = k+1; col.c = 2;	v =  mu_topwest/(dx_P_west*dz_cen);									col_vv[14] = col; v_vv[14] = v;	//  Vz(i  ,k+1)
	col.i =	  i-1; col.j = j; col.k = k+1; col.c = 2;	v = -mu_topwest/(dx_P_west*dz_cen);									col_vv[15] = col; v_vv[15] = v;	//  Vz(i-1,k+1)

	// -Sxz(k-1)/dz
	col.i = 	i; col.j = j; col.k = k;   col.c = 0;	v =  (-mu_bottomwest/dz_P_bottom)/dz_cen;							col_vv[16] = col; v_vv[16] = v;	//  Vx(i,k  )
	col.i = 	i; col.j = j; col.k = k-1; col.c = 0;	v =  mu_bottomwest/(dz_P_bottom*dz_cen);							col_vv[17] = col; v_vv[17] = v;	//  Vx(i,k-1)
	col.i =   i  ; col.j = j; col.k = k;   col.c = 2;	v = -mu_bottomwest/(dx_P_west*dz_cen);								col_vv[18] = col; v_vv[18] = v;	//  Vz(i  ,k  )
	col.i =   i-1; col.j = j; col.k = k;   col.c = 2;	v =  mu_bottomwest/(dx_P_west*dz_cen);								col_vv[19] = col; v_vv[19] = v;	//  Vz(i-1,k  )



	// add volumetric component to Txx derivatives, which is
	// (Txx_vol(i)-Txx_vol(i-1))/dx, where:
	// Txx_vol(i)   = -2/3*mu_eff *( (Vx(i+1)-Vx(i  ))/dx + (Vy(j+1)-Vy(j))/dy + (Vz(k+1)-Vz(k))/dz )
	// Txx_vol(i-1) = -2/3*mu_eff *( (Vx(i  )-Vx(i-1))/dx + (Vy(j+1)-Vy(j))/dy + (Vz(k+1)-Vz(k))/dz )

	// Txx_vol(i)/dx:
	col.i = 	i+1; col.j = j;   col.k = k  ; col.c = 0;	v =  -2.0/3.0*mu_cen/(dx_cen*dx_P_west);							col_vv[20] = col; v_vv[20] = v;			//  Vx(i+1,j  ,k)
	col.i = 	i  ; col.j = j;   col.k = k  ; col.c = 0;	v = ( 2.0/3.0*mu_cen/dx_cen)/dx_P_west;								col_vv[21] = col; v_vv[21] = v;			//  Vx(i  ,j  ,k)
	col.i = 	i  ; col.j = j+1; col.k = k  ; col.c = 1;	v =  -2.0/3.0*mu_cen/(dx_cen*dx_P_west);							col_vv[22] = col; v_vv[22] = v;			//  Vy(i  ,j+1,k)
	col.i = 	i  ; col.j = j;   col.k = k  ; col.c = 1;	v = ( 2.0/3.0*mu_cen/dx_cen)/dx_P_west;								col_vv[23] = col; v_vv[23] = v;			//  Vy(i  ,j  ,k)
	col.i = 	i  ; col.j = j;   col.k = k+1; col.c = 2;	v =  -2.0/3.0*mu_cen/(dx_cen*dx_P_west);							col_vv[24] = col; v_vv[24] = v;			//  Vz(i  ,j  ,k+1)
	col.i = 	i  ; col.j = j;   col.k = k  ; col.c = 2;	v = ( 2.0/3.0*mu_cen/dx_cen)/dx_P_west;								col_vv[25] = col; v_vv[25] = v;			//  Vz(i  ,j  ,k)

	// -Txx_vol(i-1)/dx:
	col.i = 	i  ; col.j = j;   col.k = k  ; col.c = 0;	v =   2.0/3.0*mu_west/(dx_cen*dx_P_west);							col_vv[26] = col; v_vv[26] = v;			//  Vx(i  ,j  ,k)
	col.i = 	i-1; col.j = j;   col.k = k  ; col.c = 0;	v = (-2.0/3.0*mu_west/dx_cen)/dx_P_west;							col_vv[27] = col; v_vv[27] = v;			//  Vx(i-1,j  ,k)
	col.i = 	i-1; col.j = j+1; col.k = k  ; col.c = 1;	v =   2.0/3.0*mu_west/(dx_cen*dx_P_west);							col_vv[28] = col; v_vv[28] = v;			//  Vy(i-1,j+1,k)
	col.i = 	i-1; col.j = j;   col.k = k  ; col.c = 1;	v = (-2.0/3.0*mu_west/dx_cen)/dx_P_west;							col_vv[29] = col; v_vv[29] = v;			//  Vy(i-1,j  ,k)
	col.i = 	i-1; col.j = j;   col.k = k+1; col.c = 2;	v =   2.0/3.0*mu_west/(dx_cen*dx_P_west);							col_vv[30] = col; v_vv[30] = v;			//  Vz(i-1,j  ,k+1)
	col.i = 	i-1; col.j = j;   col.k = k  ; col.c = 2;	v = (-2.0/3.0*mu_west/dx_cen)/dx_P_west;							col_vv[31] = col; v_vv[31] = v;			//  Vz(i-1,j  ,k)


	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Set BC's inside the stencil ('ghost node approach') of the 1th force balance equation  */
#undef __FUNCT__
#define __FUNCT__ "FDSTAG_SetBCs_FDStencil_ForceBalance1"
PetscErrorCode FDSTAG_SetBCs_FDStencil_ForceBalance1(PetscInt j, PetscInt k, PetscInt nel_y, PetscInt nel_z,
		PetscScalar v_vv[], UserContext *user)
{

	/* Lower/Upper BC */
	if (k==0){
		if (user->BC.LowerBound==1){
			// Free slip  - Sxz=0
			v_vv[16]=0; v_vv[17]=0; v_vv[18]=0; v_vv[19]=0;
		}
		else if (user->BC.LowerBound==2){
			// no slip [see explanation in FDSTAG_SetBCs_FDStencil_ForceBalance2]
			v_vv[16] = 2*v_vv[16];
			v_vv[17] = 0;
		}
		else{
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER,"BC.LowerBound=%lld is not implemented for FDSTAG! \n",(LLD)(user->BC.LowerBound));
		}

	}
	else if (k==nel_z-1){
		if (user->BC.UpperBound==1){
			// Free slip - Sxz=0
			v_vv[12]=0; v_vv[13]=0; v_vv[14]=0; v_vv[15]=0;
		}
		else if (user->BC.UpperBound==0){
			// Stress free - Sxz=0  and dVz/dz=0
			v_vv[12]=0; v_vv[13]=0; v_vv[14]=0; v_vv[15]=0;
		}
		else{
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER,"BC.UpperBound=%lld is not implemented for FDSTAG! \n",(LLD)(user->BC.UpperBound));
		}
	}

	/* Front/Back BC */
	if (j==0){
		if (user->BC.FrontBound==1){
			// Free slip  - Sxy=0
			v_vv[8]=0; v_vv[9]=0; v_vv[10]=0; v_vv[11]=0;
		}
		else{
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER,"BC.FrontBound=%lld is not implemented for FDSTAG! \n",(LLD)(user->BC.FrontBound));
		}
	}
	if (j==nel_y-1){
		if (user->BC.BackBound==1){
			// Free slip  - Sxy=0
			v_vv[4]=0; v_vv[5]=0; v_vv[6]=0; v_vv[7]=0;
		}
		else{
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER,"BC.BackBound=%lld is not implemented for FDSTAG! \n",(LLD)(user->BC.BackBound));
		}
	}

#if 0	// Check if we have a 'sticky-air' cell in the current cellor in the cell above.
	if (FreeSurfaceCells[1]==1){
		// current cell is sticky air
		// Set it to:  Szz(k-1) = 0

		v_vp[1]=0; // P(k)
	//	v_vp[0]=0;
	//	v_vv[0]=v_vv[1]=v_vv[2]=v_vv[3]=0;   //Szz=0 in center of current cell
		v_vv[0]=v_vv[1]=0;	//Szz(k)=0

		// stresses are zero too
		v_vv[4]=v_vv[5]=v_vv[6]=v_vv[7]=0;   //Sxz=0 in current cell
		v_vv[8]=v_vv[9]=v_vv[10]=v_vv[11]=0;  //Sxz=0 in current cell
		v_vv[12]=v_vv[13]=v_vv[14]=v_vv[15]=0;  //Syz=0 in current cell
		v_vv[16]=v_vv[17]=v_vv[18]=v_vv[18]=0;  //Syz=0 in current cell

	}
#endif

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/



/*==========================================================================================================*/
/* Computes the staggered FD stencil for the 2nd force balance equation, which is located at Vy nodes.
 *
 */
#undef __FUNCT__
#define __FUNCT__ "FDSTAG_FDStencil_ForceBalance2"
PetscErrorCode FDSTAG_FDStencil_ForceBalance2(PetscInt i, PetscInt j, PetscInt k,
		PetscScalar dx_vec[7], 		PetscScalar dy_vec[7], 		PetscScalar dz_vec[7],
		PetscScalar   dx_P[2], 		PetscScalar	dy_P[2], 		PetscScalar   dz_P[2],
		PetscScalar	Eta_Center[7], 	PetscScalar Eta_XY[4],		PetscScalar Eta_YZ[4],
		MatStencil  row_vv[], 		MatStencil  col_vv[], 		PetscScalar	v_vv[],
		MatStencil  col_vp[],  		PetscScalar v_vp[])
{
	PetscScalar		dx_P_west, dy_south, dy_P_south, dy_cen, dx_cen, dx_P_east, dz_cen, dz_P_top, dz_P_bottom;
	PetscScalar		mu_south, mu_cen, mu_southeast, mu_southwest, mu_topsouth, mu_bottomsouth;
	PetscScalar		v;
	MatStencil		row, col;

	PetscFunctionBegin;

	/* Extract info needed to construct the FD stencil */
	dx_cen 			=	dx_vec[1];
	dy_cen 			=	dy_vec[1];
	dz_cen 			=	dz_vec[1];
	dy_south 		=	dy_vec[3];

	dz_P_top		=	dz_P[1];
	dz_P_bottom		=	dz_P[0];
	dy_P_south 		=	dy_P[0];
	dx_P_west 		=	dx_P[0];
	dx_P_east   	=   dx_P[1];

	mu_south 		=	Eta_Center[3];
	mu_cen 			=	Eta_Center[1];
	mu_southeast	=	Eta_XY[1];
	mu_southwest	=   Eta_XY[0];

	mu_topsouth		=	Eta_YZ[1];
	mu_bottomsouth	=	Eta_YZ[0];


	/* Equation number [2nd Force balance is located @ Vy nodes] */
	row.i = i;   row.j = j; row.k = k; row.c = 1;	row_vv[0] = row;		// 2nd force balance equation is located @ Vy(i,j,k) nodes


	// -[ P(j) - P(j-1) ]/dy_P_south		- Pressure equations
	col.i = i; col.j = j-1; col.k = k; col.c = 0; v= 1/dy_P_south;	col_vp[0] = col;	v_vp[0] = v; 	//  P(j-1)/dy
	col.i = i; col.j = j  ; col.k = k; col.c = 0; v=-1/dy_P_south;	col_vp[1] = col;	v_vp[1] = v;	// -P(j  )/dy

	// [Syy(j)-Syy(j-1)]/dy
	col.i = 	i; col.j = j+1; col.k = k; col.c = 1;	v =    2.0*mu_cen/(dy_cen*dy_P_south);								col_vv[0] = col; v_vv[0] = v;			//  Vy(i,j+1,k)
	col.i = 	i; col.j = j  ; col.k = k; col.c = 1;	v =  (-2.0*mu_cen/dy_cen)/dy_P_south;								col_vv[1] = col; v_vv[1] = v;			//  Vy(i,j  ,k)
	col.i = 	i; col.j = j  ; col.k = k; col.c = 1;	v =  (-2.0*mu_south/dy_south)/dy_P_south;							col_vv[2] = col; v_vv[2] = v;			//  Vy(i,j  ,k)
	col.i = 	i; col.j = j-1; col.k = k; col.c = 1;	v =    2.0*mu_south/(dy_south*dy_P_south);							col_vv[3] = col; v_vv[3] = v;			//  Vy(i,j-1,k)

	// [Sxy(i+1)-Sxy(i)]/dx
	// Sxy(i+1)/dx
	col.i =   i+1; col.j = j;   col.k = k; col.c = 1;	v =  mu_southeast/(dx_P_east*dx_cen);								col_vv[4] = col; v_vv[4] = v;			//  Vy(i+1,j)
	col.i = 	i; col.j = j;   col.k = k; col.c = 1;	v =  (-mu_southeast/dx_P_east)/dx_cen;								col_vv[5] = col; v_vv[5] = v;			//  Vy(i  ,j  )
	col.i =   i+1; col.j = j-1; col.k = k; col.c = 0;	v = -mu_southeast/(dy_P_south*dx_cen);								col_vv[6] = col; v_vv[6] = v;			//  Vx(i+1,j-1)
	col.i =   i+1; col.j = j  ; col.k = k; col.c = 0;	v =  mu_southeast/(dy_P_south*dx_cen);								col_vv[7] = col;  v_vv[7]= v;			//  Vx(i+1,j  )

	// -Sxy(i)/dx
	col.i = 	i; col.j = j;   col.k = k; col.c = 1;	v =  (-mu_southwest/dx_P_west)/dx_cen;								col_vv[8] = col; v_vv[8] = v;			//  Vy(i  ,j  )
	col.i =   i-1; col.j = j;   col.k = k; col.c = 1;	v =  mu_southwest/(dx_P_west*dx_cen);								col_vv[9] = col; v_vv[9] = v;			//  Vy(i-1,j)
	col.i =	  i  ; col.j = j-1; col.k = k; col.c = 0;	v =  mu_southwest/(dy_P_south*dx_cen);								col_vv[10] = col; v_vv[10] = v;			//  Vx(i  ,j-1)
	col.i =   i  ; col.j = j  ; col.k = k; col.c = 0;	v = -mu_southwest/(dy_P_south*dx_cen);								col_vv[11]= col; v_vv[11]= v;			//  Vx(i  ,j  )

	// [Syz(k+1)-Syz(k)]/dz
	// Syz(k+1)/dz
	col.i = i; col.j = j; col.k = k+1; col.c = 1;	v =  mu_topsouth/(dz_P_top*dz_cen);										col_vv[12] = col; v_vv[12] = v;	//  Vy(i,j,k+1)
	col.i = i; col.j = j; col.k = k;   col.c = 1;	v =  (-mu_topsouth/dz_P_top)/dz_cen;									col_vv[13] = col; v_vv[13] = v;	//  Vy(i,j,k  )
	col.i = i; col.j = j  ; col.k = k+1; col.c = 2;	v =  mu_topsouth/(dy_P_south*dz_cen);									col_vv[14] = col; v_vv[14] = v;	//  Vz(i  ,j,k+1)
	col.i =	i; col.j = j-1; col.k = k+1; col.c = 2;	v = -mu_topsouth/(dy_P_south*dz_cen);									col_vv[15] = col; v_vv[15] = v;	//  Vz(i-1,j,k+1)

	// -Syz(k)/dz
	col.i = i; col.j = j; col.k = k;   col.c = 1;	v =  (-mu_bottomsouth/dz_P_bottom)/dz_cen;								col_vv[16] = col; v_vv[16] = v;	//  Vy(i,j,k  )
	col.i = i; col.j = j; col.k = k-1; col.c = 1;	v =  mu_bottomsouth/(dz_P_bottom*dz_cen);								col_vv[17] = col; v_vv[17] = v;	//  Vy(i,j,k-1)
	col.i = i; col.j = j  ; col.k = k;   col.c = 2;	v = -mu_bottomsouth/(dy_P_south*dz_cen);								col_vv[18] = col; v_vv[18] = v;	//  Vz(i  ,j,k  )
	col.i = i; col.j = j-1; col.k = k;   col.c = 2;	v =  mu_bottomsouth/(dy_P_south*dz_cen);								col_vv[19] = col; v_vv[19] = v;	//  Vz(i-1,j,k  )



	// add volumetric component to Tyy derivatives, which is

	// Tyy_vol(j)/dy:
	col.i = 	i+1; col.j = j;   col.k = k  ; col.c = 0;	v =  -2.0/3.0*mu_cen/(dy_cen*dy_P_south);							col_vv[20] = col; v_vv[20] = v;			//  Vx(i+1,j  ,k)
	col.i = 	i  ; col.j = j;   col.k = k  ; col.c = 0;	v = ( 2.0/3.0*mu_cen/dy_cen)/dy_P_south;							col_vv[21] = col; v_vv[21] = v;			//  Vx(i  ,j  ,k)
	col.i = 	i  ; col.j = j+1; col.k = k  ; col.c = 1;	v =  -2.0/3.0*mu_cen/(dy_cen*dy_P_south);							col_vv[22] = col; v_vv[22] = v;			//  Vy(i  ,j+1,k)
	col.i = 	i  ; col.j = j;   col.k = k  ; col.c = 1;	v = ( 2.0/3.0*mu_cen/dy_cen)/dy_P_south;							col_vv[23] = col; v_vv[23] = v;			//  Vy(i  ,j  ,k)
	col.i = 	i  ; col.j = j;   col.k = k+1; col.c = 2;	v =  -2.0/3.0*mu_cen/(dy_cen*dy_P_south);							col_vv[24] = col; v_vv[24] = v;			//  Vz(i  ,j  ,k+1)
	col.i = 	i  ; col.j = j;   col.k = k  ; col.c = 2;	v = ( 2.0/3.0*mu_cen/dy_cen)/dy_P_south;							col_vv[25] = col; v_vv[25] = v;			//  Vz(i  ,j  ,k)

	// -Tyy_vol(j-1)/dy:
	col.i = 	i+1; col.j = j-1;   col.k = k  ; col.c = 0;	v =   2.0/3.0*mu_south/(dy_cen*dy_P_south);							col_vv[26] = col; v_vv[26] = v;
	col.i = 	i  ; col.j = j-1;   col.k = k  ; col.c = 0;	v = (-2.0/3.0*mu_south/dy_cen)/dy_P_south;							col_vv[27] = col; v_vv[27] = v;
	col.i = 	i  ; col.j = j  ;   col.k = k  ; col.c = 1;	v =   2.0/3.0*mu_south/(dy_cen*dy_P_south);							col_vv[28] = col; v_vv[28] = v;
	col.i = 	i  ; col.j = j-1;   col.k = k  ; col.c = 1;	v = (-2.0/3.0*mu_south/dy_cen)/dy_P_south;							col_vv[29] = col; v_vv[29] = v;
	col.i = 	i  ; col.j = j-1;   col.k = k+1; col.c = 2;	v =   2.0/3.0*mu_south/(dy_cen*dy_P_south);							col_vv[30] = col; v_vv[30] = v;
	col.i = 	i  ; col.j = j-1;   col.k = k  ; col.c = 2;	v = (-2.0/3.0*mu_south/dy_cen)/dy_P_south;							col_vv[31] = col; v_vv[31] = v;


	PetscFunctionReturn(0);
}
/*==========================================================================================================*/

/*==========================================================================================================*/
/* Set BC's inside the stencil ('ghost node approach') of the 2nd force balance equation  */
#undef __FUNCT__
#define __FUNCT__ "FDSTAG_SetBCs_FDStencil_ForceBalance2"
PetscErrorCode FDSTAG_SetBCs_FDStencil_ForceBalance2(PetscInt i, PetscInt k, PetscInt nel_x, PetscInt nel_z,
		PetscScalar v_vv[], UserContext *user)
{
	/* Lower/Upper BC */
	if (k==0){
		if (user->BC.LowerBound==1){
			// Free slip  - Syz=0
			v_vv[16]=0; v_vv[17]=0; v_vv[18]=0; v_vv[19]=0;
		}
		else if (user->BC.LowerBound==2){
			// No slip - Vx=Vy=0 at z=z_bot,  which implies that Vy(i,j,k-1) = -Vy(i,j,k)   (as this will precisely give 0 at k-1/2 which is where the boundary is)
			//
			// this changes Sxy(k)/dz= 1/dz*( [Vy(k)- Vy(k-1)]/dz + [Vz(i)-Vz(i-1)]/dx  )
			//  becomes: 	Sxy(k)/dz= 1/dz*( [2*Vy(k)    )  ]/dz + [Vz(i)-Vz(i-1)]/dx  )
			v_vv[16] = v_vv[16]*2.0;
			v_vv[17] = 0;

		}
		else{
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER,"BC.LowerBound=%lld is not implemented for FDSTAG! \n",(LLD)(user->BC.LowerBound));
		}

	}
	else if (k==nel_z-1){
		if (user->BC.UpperBound==1){
			// Free slip - Syz=0
			v_vv[12]=0; v_vv[13]=0; v_vv[14]=0; v_vv[15]=0;
		}
		else if (user->BC.UpperBound==0){
			// Free surface - Syz=0
			v_vv[12]=0; v_vv[13]=0; v_vv[14]=0; v_vv[15]=0;
		}
		else{
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER,"BC.UpperBound=%lld is not implemented for FDSTAG! \n",(LLD)(user->BC.UpperBound));
		}
	}

	/* Left/Right BC */
	if (i==0){
		if (user->BC.LeftBound==1){
			// Free slip  - Sxy=0
			v_vv[8]=0; v_vv[9]=0; v_vv[10]=0; v_vv[11]=0;
		}
		else{
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER,"BC.LeftBound=%lld is not implemented for FDSTAG! \n",(LLD)(user->BC.FrontBound));
		}
	}
	if (i==nel_x-1){
		if (user->BC.RightBound==1){
			// Free slip  - Sxy=0
			v_vv[4]=0; v_vv[5]=0; v_vv[6]=0; v_vv[7]=0;
		}
		else{
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER,"BC.RightBound=%lld is not implemented for FDSTAG! \n",(LLD)(user->BC.BackBound));
		}
	}


#if 0	// Check if we have a 'sticky-air' cell in the current cellor in the cell above.
	if (FreeSurfaceCells[1]==1){
		// current cell is sticky air
		// Set it to:  Szz(k-1) = 0

		v_vp[1]=0; // P(k)
	//	v_vp[0]=0;
	//	v_vv[0]=v_vv[1]=v_vv[2]=v_vv[3]=0;   //Szz=0 in center of current cell
		v_vv[0]=v_vv[1]=0;	//Szz(k)=0

		// stresses are zero too
		v_vv[4]=v_vv[5]=v_vv[6]=v_vv[7]=0;   //Sxz=0 in current cell
		v_vv[8]=v_vv[9]=v_vv[10]=v_vv[11]=0;  //Sxz=0 in current cell
		v_vv[12]=v_vv[13]=v_vv[14]=v_vv[15]=0;  //Syz=0 in current cell
		v_vv[16]=v_vv[17]=v_vv[18]=v_vv[18]=0;  //Syz=0 in current cell

	}
#endif

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/

/*==========================================================================================================*/
/* Computes the staggered FD stencil for the 3rd force balance equation, which is located at Vz nodes.
 *
 */
#undef __FUNCT__
#define __FUNCT__ "FDSTAG_FDStencil_ForceBalance3"
PetscErrorCode FDSTAG_FDStencil_ForceBalance3(PetscInt i, PetscInt j, PetscInt k,
		PetscScalar dx_vec[7], 		PetscScalar dy_vec[7], 		PetscScalar dz_vec[7],
		PetscScalar   dx_P[2], 		PetscScalar	dy_P[2], 		PetscScalar   dz_P[2],
		PetscScalar	Eta_Center[7], 	PetscScalar Eta_XZ[4],		PetscScalar Eta_YZ[4],
		MatStencil  row_vv[], 		MatStencil  col_vv[], 		PetscScalar	v_vv[],
		MatStencil  col_vp[],  		PetscScalar v_vp[], 		PetscScalar	FreeSurface_Fraction)
{
	PetscScalar		dz_P_bottom, dz_cen, dz_bottom, dx_cen, dy_cen, dx_P_east, dx_P_west, dy_P_south, dy_P_north;
	PetscScalar		mu_cen, mu_bottom, mu_bottomeast, mu_bottomwest, mu_bottomnorth, mu_bottomsouth;
	PetscScalar		v;
	MatStencil		row, col;

	PetscFunctionBegin;

	/* Extract info needed to construct the FD stencil */
	dz_P_bottom		=	dz_P[0];
	dz_cen 			=	dz_vec[1];
	dz_bottom 		=	dz_vec[5];
	dx_cen 			=	dx_vec[1];
	dy_cen 			=	dy_vec[1];
	dx_P_west 		=	dx_P[0];
	dx_P_east 		=	dx_P[1];
	dy_P_south 		=	dy_P[0];
	dy_P_north 		=	dy_P[1];

	mu_cen 			=	Eta_Center[1];
	mu_bottom 		=	Eta_Center[5];
	mu_bottomwest	=   Eta_XZ[0];
	mu_bottomeast	=	Eta_XZ[2];
	mu_bottomsouth	=	Eta_YZ[0];
	mu_bottomnorth	=	Eta_YZ[2];



	/* Equation number [3rd Force balance is located @ Vz nodes] */
	row.i = i;   row.j = j; row.k = k; row.c = 2;	row_vv[0] = row;		// 3rd force balance equation is located @ Vz(i,j,k) nodes


	// -[ P(k) - P(k-1) ]/dz_P_bottom		- Pressure equations
	col.i = i; col.j = j; col.k = k-1; col.c = 0; v= 1/(dz_P_bottom*FreeSurface_Fraction);	col_vp[0] = col;	v_vp[0] = v; 	//  P(k-1)/dz
	col.i = i; col.j = j; col.k = k  ; col.c = 0; v=-1/(dz_P_bottom*FreeSurface_Fraction);	col_vp[1] = col;	v_vp[1] = v;	// -P(k  )/dz


	// [Szz(k)-Szz(k-1)]/dz
	col.i = 	i; col.j = j; col.k = k+1; col.c = 2;	v =    2.0*mu_cen/(dz_cen*dz_P_bottom);								col_vv[0] = col; v_vv[0] = v;			//  Vz(i,j,k+1)
	col.i = 	i; col.j = j; col.k = k  ; col.c = 2;	v =  (-2.0*mu_cen/dz_cen)/dz_P_bottom;								col_vv[1] = col; v_vv[1] = v;			//  Vz(i,j,k)
	col.i = 	i; col.j = j; col.k = k  ; col.c = 2;	v =  (-2.0*mu_bottom/dz_bottom)/dz_P_bottom;						col_vv[2] = col; v_vv[2] = v;			//  Vz(i,j,k)
	col.i = 	i; col.j = j; col.k = k-1; col.c = 2;	v =    2.0*mu_bottom/(dz_bottom*dz_P_bottom);						col_vv[3] = col; v_vv[3] = v;			//  Vz(i,j,k-1)


	// [Sxz(i+1)-Sxz(i)]/dx
	// Sxz(i+1)/dx
	col.i = 	i+1; col.j = j; col.k = k; col.c = 2;	v =  mu_bottomeast/(dx_P_east*dx_cen);								col_vv[4] = col; v_vv[4] = v;			//  Vz(i+1,j,k)
	col.i = 	i; 	 col.j = j; col.k = k; col.c = 2;	v =  (-mu_bottomeast/dx_P_east)/dx_cen;								col_vv[5] = col; v_vv[5] = v;			//  Vz(i  ,j,k)
	col.i =   i+1; col.j = j; col.k = k  ; col.c = 0;	v = mu_bottomeast/(dz_P_bottom*dx_cen);								col_vv[6] = col; v_vv[6] = v;			//  Vx(i+1,j, k  )
	col.i =	  i+1; col.j = j; col.k = k-1; col.c = 0;	v =-mu_bottomeast/(dz_P_bottom*dx_cen);								col_vv[7] = col; v_vv[7] = v;			//  Vx(i+1,j, k-1)

	// -Sxz(i)/dx
	col.i = 	i; 	 col.j = j; col.k = k; col.c = 2;	v =  (-mu_bottomwest/dx_P_west)/dx_cen;								col_vv[8] = col; v_vv[8] = v;			//  Vz(i  ,j,k)
	col.i = 	i-1; col.j = j; col.k = k; col.c = 2;	v =  mu_bottomwest/(dx_P_west*dx_cen);								col_vv[9] = col; v_vv[9] = v;			//  Vz(i-1,j,k)
	col.i =   i  ; col.j = j; col.k = k  ; col.c = 0;	v =-mu_bottomwest/(dz_P_bottom*dx_cen);								col_vv[10] = col; v_vv[10] = v;			//  Vx(i  ,j, k  )
	col.i =   i  ; col.j = j; col.k = k-1; col.c = 0;	v = mu_bottomwest/(dz_P_bottom*dx_cen);								col_vv[11]= col; v_vv[11]= v;			//  Vx(i  ,j, k-1)

	// [Syz(j+1)-Syz(j)]/dy
	// Syz(j+1)/dy
	col.i = i; col.j = j+1; col.k = k; col.c = 2;	v =  mu_bottomnorth/(dy_P_north*dy_cen);								col_vv[12] = col; v_vv[12] = v;	//  Vz(i,j+1,k)
	col.i = i; col.j = j  ; col.k = k; col.c = 2;	v =  (-mu_bottomnorth/dy_P_north)/dy_cen;								col_vv[13] = col; v_vv[13] = v;	//  Vz(i,j  ,k)
	col.i = i; col.j = j+1; col.k = k  ; col.c = 1;	v =  mu_bottomnorth/(dz_P_bottom*dy_cen);								col_vv[14] = col; v_vv[14] = v;	//  Vy(i,j+1,k  )
	col.i =	i; col.j = j+1; col.k = k-1; col.c = 1;	v = -mu_bottomnorth/(dz_P_bottom*dy_cen);								col_vv[15] = col; v_vv[15] = v;	//  Vy(i,j+1,k-1)

	// -Syz(j)/dy
	col.i = i; col.j = j  ; col.k = k; col.c = 2;	v =  (-mu_bottomsouth/dy_P_south)/dy_cen;								col_vv[16] = col; v_vv[16] = v;	//  Vz(i,j  ,k)
	col.i = i; col.j = j-1; col.k = k; col.c = 2;	v =  mu_bottomsouth/(dy_P_south*dy_cen);								col_vv[17] = col; v_vv[17] = v;	//  Vz(i,j-1,k)
	col.i = i; col.j = j;   col.k = k  ; col.c = 1;	v = -mu_bottomsouth/(dz_P_bottom*dy_cen);								col_vv[18] = col; v_vv[18] = v;	//  Vy(i,j  ,k  )
	col.i = i; col.j = j;   col.k = k-1; col.c = 1;	v =  mu_bottomsouth/(dz_P_bottom*dy_cen);								col_vv[19] = col; v_vv[19] = v;	//  Vy(i,j  ,k-1)


	// add volumetric component to Tzz derivatives, which is

	// Tzz_vol(j)/dz:
	col.i = 	i+1; col.j = j;   col.k = k  ; col.c = 0;	v =  -2.0/3.0*mu_cen/(dz_cen*dz_P_bottom);							col_vv[20] = col; v_vv[20] = v;			//  Vx(i+1,j  ,k)
	col.i = 	i  ; col.j = j;   col.k = k  ; col.c = 0;	v = ( 2.0/3.0*mu_cen/dz_cen)/dz_P_bottom;							col_vv[21] = col; v_vv[21] = v;			//  Vx(i  ,j  ,k)
	col.i = 	i  ; col.j = j+1; col.k = k  ; col.c = 1;	v =  -2.0/3.0*mu_cen/(dz_cen*dz_P_bottom);							col_vv[22] = col; v_vv[22] = v;			//  Vy(i  ,j+1,k)
	col.i = 	i  ; col.j = j;   col.k = k  ; col.c = 1;	v = ( 2.0/3.0*mu_cen/dz_cen)/dz_P_bottom;							col_vv[23] = col; v_vv[23] = v;			//  Vy(i  ,j  ,k)
	col.i = 	i  ; col.j = j;   col.k = k+1; col.c = 2;	v =  -2.0/3.0*mu_cen/(dz_cen*dz_P_bottom);							col_vv[24] = col; v_vv[24] = v;			//  Vz(i  ,j  ,k+1)
	col.i = 	i  ; col.j = j;   col.k = k  ; col.c = 2;	v = ( 2.0/3.0*mu_cen/dz_cen)/dz_P_bottom;							col_vv[25] = col; v_vv[25] = v;			//  Vz(i  ,j  ,k)

	// -Tzz_vol(k-1)/dz:
	col.i = 	i+1; col.j = j  ;   col.k = k-1; col.c = 0;	v =   2.0/3.0*mu_bottom/(dz_cen*dz_P_bottom);							col_vv[26] = col; v_vv[26] = v;
	col.i = 	i  ; col.j = j  ;   col.k = k-1; col.c = 0;	v = (-2.0/3.0*mu_bottom/dz_cen)/dz_P_bottom;							col_vv[27] = col; v_vv[27] = v;
	col.i = 	i  ; col.j = j+1;   col.k = k-1; col.c = 1;	v =   2.0/3.0*mu_bottom/(dz_cen*dz_P_bottom);							col_vv[28] = col; v_vv[28] = v;
	col.i = 	i  ; col.j = j  ;   col.k = k-1; col.c = 1;	v = (-2.0/3.0*mu_bottom/dz_cen)/dz_P_bottom;							col_vv[29] = col; v_vv[29] = v;
	col.i = 	i  ; col.j = j  ;   col.k = k  ; col.c = 2;	v =   2.0/3.0*mu_bottom/(dz_cen*dz_P_bottom);							col_vv[30] = col; v_vv[30] = v;
	col.i = 	i  ; col.j = j  ;   col.k = k-1; col.c = 2;	v = (-2.0/3.0*mu_bottom/dz_cen)/dz_P_bottom;							col_vv[31] = col; v_vv[31] = v;



	PetscFunctionReturn(0);
}
/*==========================================================================================================*/

/*==========================================================================================================*/
/* Set BC's inside the stencil ('ghost node approach') of the 3rd force balance equation  */
#undef __FUNCT__
#define __FUNCT__ "FDSTAG_SetBCs_FDStencil_ForceBalance3"
PetscErrorCode FDSTAG_SetBCs_FDStencil_ForceBalance3(PetscInt i, PetscInt j, PetscInt nel_x, PetscInt nel_y,
		PetscScalar v_vv[], PetscScalar v_vp[], PetscInt FreeSurfaceCells[], UserContext *user, PetscBool eliminate_stickyair_from_system)
{
	/* Front/Back BC */
	if (j==0){
		if (user->BC.FrontBound==1){
			// Free slip  - Syz=0
			v_vv[16]=0; v_vv[17]=0; v_vv[18]=0; v_vv[19]=0;
		}
		else{
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER,"BC.FrontBound=%lld is not implemented for FDSTAG! \n",(LLD)(user->BC.FrontBound));
		}

	}
	else if (j==nel_y-1){
		if (user->BC.BackBound==1){
			// Free slip - Syz=0
			v_vv[12]=0; v_vv[13]=0; v_vv[14]=0; v_vv[15]=0;
		}
		else{
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER,"BC.BackBound=%lld is not implemented for FDSTAG! \n",(LLD)(user->BC.BackBound));
		}
	}

	/* Left/Right BC */
	if (i==0){
		if (user->BC.LeftBound==1){
			// Free slip  - Sxz=0
			v_vv[8]=0; v_vv[9]=0; v_vv[10]=0; v_vv[11]=0;
		}
		else{
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER,"BC.LeftBound=%lld is not implemented for FDSTAG! \n",(LLD)(user->BC.LeftBound));
		}
	}
	if (i==nel_x-1){
		if (user->BC.RightBound==1){
			// Free slip  - Sxz=0
			v_vv[4]=0; v_vv[5]=0; v_vv[6]=0; v_vv[7]=0;
		}
		else{
			SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER,"BC.RightBound=%lld is not implemented for FDSTAG! \n",(LLD)(user->BC.RightBound));
		}
	}

	// Check if we have a 'sticky-air' cell in the current cellor in the cell above.
	if ((eliminate_stickyair_from_system) && (FreeSurfaceCells[1]==1)){
		// current cell is sticky air
		// Set it to:  Szz(k-1) = 0

		// Set Szz(k-1)=0
		v_vp[1]=0; // P(k)
		v_vv[0]=v_vv[1]=0;	//Szz(k)=0

		// stresses are zero too
		v_vv[4]=v_vv[5]=v_vv[6]=v_vv[7]=0;   //Sxz=0 in current cell
		v_vv[8]=v_vv[9]=v_vv[10]=v_vv[11]=0;  //Sxz=0 in current cell
		v_vv[12]=v_vv[13]=v_vv[14]=v_vv[15]=0;  //Syz=0 in current cell

		// 
		v_vv[16]=v_vv[17]=v_vv[18]=v_vv[19]=0;  //Syz=0 in current cell

	}

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Compute the right-hand-side vector */
#undef __FUNCT__
#define __FUNCT__ "ComputeRHS_FDSTAG"
PetscErrorCode ComputeRHS_FDSTAG(DM da, Vec b, UserContext *user)
{
	Vec			 			bloc;
	Field					***rhs;
	PetscErrorCode 			ierr;
	PetscScalar				***density_Vz, rho;
	Vec						Density_local;
	PetscInt       			i,j,k;
	PetscInt 				xsp,ysp,zsp,xmp,ymp,zmp;
	PetscInt		 		iel_x, iel_y, iel_z;
	PetscScalar	 			Gravity_z;

	PetscFunctionBegin;

	Gravity_z 		=	 cos(user->GravityAngle/180*M_PI)*user->Gravity;

	/* Get rhs vector */
	ierr = DMGetLocalVector(da,&bloc); 													CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin	(da, 	b, 		INSERT_VALUES, bloc); 							CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd	(da, 	b, 		INSERT_VALUES, bloc); 								CHKERRQ(ierr);
	ierr = DMDAVecGetArray		(da,	bloc,	&rhs);	CHKERRQ(ierr);


	/* Extract density @ Vz points */
	ierr = DMGetLocalVector		(user->FDSTAG.DA_CENTER,&Density_local);	CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin	(user->FDSTAG.DA_CENTER,user->FDSTAG.Center_Density,INSERT_VALUES,Density_local);	CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd  	(user->FDSTAG.DA_CENTER,user->FDSTAG.Center_Density,INSERT_VALUES,Density_local);	CHKERRQ(ierr);
	ierr = DMDAVecGetArray(user->FDSTAG.DA_CENTER,   	Density_local, 		&density_Vz	 ); CHKERRQ(ierr); // Density @ Vz points

	/* Compute RHS for each element; ignoring BC's */
	ierr = DMDAGetCorners(user->DA_Processors,&xsp,&ysp,&zsp,&xmp,&ymp,&zmp);CHKERRQ(ierr);
	for (iel_z=zsp; iel_z<zsp+zmp; iel_z++){
		for (iel_y=ysp; iel_y<ysp+ymp; iel_y++){
			for(iel_x=xsp; iel_x<xsp+xmp; iel_x++){
				i = iel_x;	j = iel_y; k = iel_z;		// Initialize variables

				/* Developer note:  in case we use a 'sticky-air' setup, we must ensure that rho=0 within the sticky air
				 * A check for this is done in Utils.c
				 * If this is not done correctly, setting a stress BC (Szz=0) within the air results in spurious (large!) velocities within the air
				 *
				 */

				if (iel_z>0){

					/* Density at Vz points is averaged from density at center points */
					rho = (density_Vz[k-1][j][i]+density_Vz[k][j][i])/2;

				}
				else{
					rho = density_Vz[k][j][i];
				}

				/* Equation number */
				rhs[k][j][i].Vz = rho*Gravity_z;


			}
		}
	}

	/* restore density arrays */
	ierr = DMDAVecRestoreArray(user->FDSTAG.DA_CENTER,   	Density_local, 		&density_Vz	 ); CHKERRQ(ierr); // density @ Vz points
	ierr = DMRestoreLocalVector(user->FDSTAG.DA_CENTER,&Density_local);	CHKERRQ(ierr);



	/* restore rhs vectors */
	ierr = DMDAVecRestoreArray(da,bloc,&rhs);	CHKERRQ(ierr);
	// update including ghost points
	ierr = DMLocalToGlobalBegin(da, bloc, ADD_VALUES, b);	 CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd  (da, bloc, ADD_VALUES, b);    CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(da,&bloc); CHKERRQ(ierr);

	SetBoundaryConditionsRHS_FDSTAG(da, user, b, 0);


	PetscFunctionReturn(0);
}
/*==========================================================================================================*/

/*==========================================================================================================*/
/* Set's the boundary conditions in the RHS vector */
#undef __FUNCT__
#define __FUNCT__ "SetBoundaryConditionsRHS_FDSTAG"
PetscErrorCode SetBoundaryConditionsRHS_FDSTAG( DM da, UserContext *user, Vec b, PetscInt SetNonZeroValuesFlag)
{
	PetscErrorCode 		ierr;
	PetscInt			i,j,k,xs,xm,ys,ym,zs,zm, mx,my,mz;
	PetscScalar 		factor;
	DM					cda;
	Field				***rhs;
	Vec					gc;
	DMDACoor3d			***coords;

	ierr = DMDAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0);	CHKERRQ(ierr);  	// # of nodes in all directions
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);	CHKERRQ(ierr);

	ierr = DMGetCoordinateDM(da,&cda); CHKERRQ(ierr);	//coordinates
	ierr = DMGetCoordinatesLocal(da,&gc); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,gc,&coords);			CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, b,	   &rhs);			CHKERRQ(ierr);


	/* If we symmetrize BC's at the element-level, we should NOT any values to the RHS here, since they
	* were already added to the BC vector*/
	// not yet implemented for FDSTAG

	factor = 1.0;

	for (k=zs; k<zs+zm; k++){
		for (j=ys; j<ys+ym; j++){
			for(i=xs; i<xs+xm; i++){
				if (i==0  ){	/* Left boundary */
					if 		(user->BC.LeftBound==1){ // free slip w. specified BG strainrate
						if (SetNonZeroValuesFlag==1){
							rhs[k][j][i].Vx  = -user->BC.Exx*coords[k  ][j  ][i  ].x*factor;
						}
						else{
							rhs[k][j][i].Vx  = 0.0;
						}

					}
					else if (user->BC.LeftBound==2){ // no slip; Vx=0,Vy=0,Vz=0
						rhs[k][j][i].Vx  = 0.0;	rhs[k][j][i].Vy  = 0.0;	  	rhs[k][j][i].Vz  = 0.0;
					}
					else{
						SETERRQ1(PETSC_COMM_WORLD, 1,"Unknown left boundary condition BC.LeftBound=%lld",(LLD)(user->BC.LeftBound));
					}

				}
				if (i==mx-1  ){	/* Right boundary */
					if 		((user->BC.RightBound==1) || (user->BC.RightBound==7)){ // free slip w. specified BG strainrate, growth rate
						if (SetNonZeroValuesFlag==1){
							rhs[k][j][i].Vx  = -user->BC.Exx*coords[k  ][j  ][i  ].x*factor;
						}
						else {
							rhs[k][j][i].Vx  = 0.0;
						}
					}
					else if (user->BC.RightBound==2){ // no slip; Vx=0,Vy=0,Vz=0
						rhs[k][j][i].Vx  = 0.0;	rhs[k][j][i].Vy  = 0.0;	  	rhs[k][j][i].Vz  = 0.0;
					}
					else if (user->BC.RightBound==3){
						// periodic - taken care of already

					}
					else{
						SETERRQ1(PETSC_COMM_WORLD, 1,"Unknown right boundary condition BC.RightBound=%lld",(LLD)(user->BC.RightBound));
					}
				}
				if (j==0){					/* Front boundary */
					if 		(user->BC.FrontBound==1){ // free slip w. specified BG strainrate
						if (SetNonZeroValuesFlag==1){
							rhs[k][j][i].Vy  = -user->BC.Eyy*coords[k  ][j  ][i  ].y*factor;
						}
						else{
							rhs[k][j][i].Vy  = 0.0;
						}

					}
					else if (user->BC.FrontBound==2){ // no slip
						rhs[k][j][i].Vx  = 0.0;  	rhs[k][j][i].Vy  = 0.0;	  	rhs[k][j][i].Vz  = 0.0;
					}
					else{
						SETERRQ1(PETSC_COMM_WORLD, 1,"Unknown front boundary condition BC.FrontBound=%lld",(LLD)(user->BC.FrontBound));
					}

				}
				if (j==my-1){				/* Back boundary */
					if 		((user->BC.BackBound==1) || (user->BC.BackBound==7)){ // free slip w. specified BG strainrate, growth-rate
						if (SetNonZeroValuesFlag==1){
							rhs[k][j][i].Vy  = -user->BC.Eyy*coords[k  ][j  ][i  ].y*factor;
						}
						else {
							rhs[k][j][i].Vy  = 0.0;
						}
					}
					else if (user->BC.BackBound==2){ // no slip
						rhs[k][j][i].Vx  = 0.0;	rhs[k][j][i].Vy  = 0.0;	 	rhs[k][j][i].Vz  = 0.0;
					}
					else{
						SETERRQ1(PETSC_COMM_WORLD, 1,"Unknown back boundary condition BC.BackBound=%lld",(LLD)(user->BC.BackBound));
					}
				}
				if (k==0 	){	/* Bottom boundary */
					if 		((user->BC.LowerBound==1) || (user->BC.LowerBound==6) || (user->BC.LowerBound==7)){ // free slip, growth rate, thin-sheet
						if (SetNonZeroValuesFlag==1){
							rhs[k][j][i].Vz  = (user->BC.Exx+user->BC.Eyy)*coords[k  ][j  ][i  ].z*factor;
						}
						else {
							rhs[k][j][i].Vz  = 0.0;
						}
					}
					else if (user->BC.LowerBound==2){ // no slip
						rhs[k][j][i].Vx  = 0.0;	rhs[k][j][i].Vy  = 0.0;	 	rhs[k][j][i].Vz  = 0.0;
					}
					else{
						SETERRQ1(PETSC_COMM_WORLD, 1,"Unknown lower boundary condition BC.LowerBound=%lld",(LLD)(user->BC.LowerBound));
					}
				}
				if (k==mz-1	){	/* Top boundary */
					if 		(user->BC.UpperBound==1){ // free slip
						if (SetNonZeroValuesFlag==1){
							rhs[k][j][i].Vz  = (user->BC.Exx+user->BC.Eyy)*coords[k  ][j  ][i  ].z*factor;
						}
						else {
							rhs[k][j][i].Vz  = 0.0;
						}

					}
					else if (user->BC.UpperBound==2){ // no slip
						rhs[k][j][i].Vx  = 0.0;	rhs[k][j][i].Vy  = 0.0;	 	rhs[k][j][i].Vz  = 0.0;
					}
					else if (user->BC.UpperBound==0){ // no stress
						rhs[k][j][i].Vz  = 0.0;
					}
					else{
						SETERRQ1(PETSC_COMM_WORLD, 1,"Unknown upper boundary condition BC.UpperBound=%lld",(LLD)(user->BC.UpperBound));
					}

				}

			}
		}
	}
	ierr = DMDAVecRestoreArray(da ,b,	   &rhs);	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(cda,gc,&coords);	CHKERRQ(ierr);

//	DMDestroy( cda );
//	VecDestroy( gc );


	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Set an initial velocity solution vector (for debugging purposes - check if the force balance equations are set in the stiffness matrix)
 *
 * BK, November 2010
 */
#undef __FUNCT__
#define __FUNCT__ "FDSTAG_SetInitialSolutionVector"
PetscErrorCode FDSTAG_SetInitialSolutionVector(UserContext *user, Vec sol_vel, Vec Pressure)
{
	PetscErrorCode 	ierr;
	PetscInt		i,j,k, xsp,ysp,zsp,xmp,ymp,zmp, test_equation;
	PetscScalar		TwoPi, W,L,H;
	PetscScalar		X_Vx, Y_Vx, Z_Vx;
	PetscScalar		X_Vy, Y_Vy, Z_Vy;
	PetscScalar		X_Vz, Y_Vz, Z_Vz;
	PetscScalar		X_P,  Y_P,  Z_P;
	PetscScalar		***pressure;
	Field 			***velocity;
	DM			 	cda;
	Vec			 	gc;
	DMDACoor3d		***coords;

	PetscFunctionBegin;

	TwoPi	= 2.0*M_PI;
	TwoPi=1;
	W 		=	user->W;
	L 		=	user->L;
	H 		=	user->H;


	/* Obtain coordinates */
	ierr = DMGetCoordinateDM(user->DA_Vel,&cda); 		CHKERRQ(ierr);	//coordinates
	ierr = DMGetCoordinatesLocal(user->DA_Vel,&gc); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,gc,&coords); 	CHKERRQ(ierr);

	/* Extract solution */
	ierr = DMDAVecGetArray(user->DA_Vel,   sol_vel, &velocity); 	CHKERRQ(ierr);
	ierr = DMDAVecGetArray(user->DA_Pres,	Pressure,&pressure); 	CHKERRQ(ierr);


	/* Make a loop over all local elements, construct the element stiffness matrix */
	test_equation = 0; 		// 0 - just numbers, 1-incompressibility, 2-1th force balance (isoviscous), 3-2nd FB (isoviscous), 4-3rd FB (isoviscous); 5-1th FB exp(x) viscosity

	// Vx
	ierr = DMDAGetCorners(user->DA_Vel,&xsp,&ysp,&zsp,&xmp,&ymp,&zmp);CHKERRQ(ierr);
	for (k=zsp; k<zsp+zmp-1; k++){
		for (j=ysp; j<ysp+ymp-1; j++){
			for(i=xsp; i<xsp+xmp; i++){

				/* Extract the coordinates for Vx, Vy & Vz (as they are shifted we have to reconstruct them here) */
				X_Vx =  ( coords[k  ][j  ][i  ].x + coords[k+1][j  ][i  ].x + coords[k  ][j+1][i  ].x + coords[k+1][j+1][i  ].x )/4.0;
				Y_Vx =  ( coords[k  ][j  ][i  ].y + coords[k+1][j  ][i  ].y + coords[k  ][j+1][i  ].y + coords[k+1][j+1][i  ].y )/4.0;
				Z_Vx =  ( coords[k  ][j  ][i  ].z + coords[k+1][j  ][i  ].z + coords[k  ][j+1][i  ].z + coords[k+1][j+1][i  ].z )/4.0;

				//===========================================================================
				if 		(test_equation==0){
					velocity[k][j][i].Vx = 1;
				}
				else if (test_equation==1){
					// Test incompressibility equation - See MAPLE derivation
					velocity[k][j][i].Vx = cos(TwoPi/W*X_Vx)*cos(TwoPi/L*Y_Vx)*cos(TwoPi/H*Z_Vx);
				}
				else if (test_equation==2){
					// Test 1th force balance equation - See MAPLE derivation
					velocity[k][j][i].Vx = cos(TwoPi/W*X_Vx)*cos(TwoPi/L*Y_Vx)*cos(TwoPi/H*Z_Vx);
				}
				else if (test_equation==3){
					// Test 2nd force balance equation - See MAPLE derivation
					velocity[k][j][i].Vx = cos(TwoPi/W*X_Vx)*cos(TwoPi/L*Y_Vx)*cos(TwoPi/H*Z_Vx);
				}
				else if (test_equation==4){
					// Test 3rd force balance equation - See MAPLE derivation
					velocity[k][j][i].Vx = cos(TwoPi/W*X_Vx)*cos(TwoPi/L*Y_Vx)*cos(TwoPi/H*Z_Vx);
				}
				else if (test_equation==5){
					// Test 3rd force balance equation - See MAPLE derivation
					velocity[k][j][i].Vx = cos(TwoPi/W*X_Vx)*cos(TwoPi/L*Y_Vx)*cos(TwoPi/H*Z_Vx);
				}
				//===========================================================================

			}
		}
	}

	// Vy
	for (k=zsp; k<zsp+zmp-1; k++){
		for (j=ysp; j<ysp+ymp; j++){
			for(i=xsp; i<xsp+xmp-1; i++){

				/* Extract the coordinates for Vx, Vy & Vz (as they are shifted we have to reconstruct them here) */
				X_Vy =  ( coords[k  ][j  ][i  ].x + coords[k+1][j  ][i  ].x + coords[k  ][j  ][i+1].x + coords[k+1][j  ][i+1].x )/4.0;
				Y_Vy =  ( coords[k  ][j  ][i  ].y + coords[k+1][j  ][i  ].y + coords[k  ][j  ][i+1].y + coords[k+1][j  ][i+1].y )/4.0;
				Z_Vy =  ( coords[k  ][j  ][i  ].z + coords[k+1][j  ][i  ].z + coords[k  ][j  ][i+1].z + coords[k+1][j  ][i+1].z )/4.0;

				//===========================================================================
				if (test_equation==0){
					velocity[k][j][i].Vy = 2;
				}
				else if (test_equation==1){
					velocity[k][j][i].Vy = sin(TwoPi/W*X_Vy)*sin(TwoPi/L*Y_Vy)*cos(TwoPi/H*Z_Vy)/2.0;
				}
				else if (test_equation==2){
					velocity[k][j][i].Vy = -sin(TwoPi/W*X_Vy)*sin(TwoPi/L*Y_Vy)*cos(TwoPi/H*Z_Vy);
				}
				else if (test_equation==3){
					velocity[k][j][i].Vy = -sin(TwoPi/W*X_Vy)*sin(TwoPi/L*Y_Vy)*cos(TwoPi/H*Z_Vy);
				}
				else if (test_equation==4){
					velocity[k][j][i].Vy = -sin(TwoPi/W*X_Vy)*sin(TwoPi/L*Y_Vy)*cos(TwoPi/H*Z_Vy);
				}
				else if (test_equation==5){
					velocity[k][j][i].Vy = -sin(TwoPi/W*X_Vy)*sin(TwoPi/L*Y_Vy)*cos(TwoPi/H*Z_Vy);
				}
				//===========================================================================


//				PetscPrintf(PETSC_COMM_WORLD,"  coord = [%g,%g,%g] - [%g,%g,%g] - [%g,%g,%g] \n",X_Vx, Y_Vx, Z_Vx, X_Vy, Y_Vy, Z_Vy, X_Vz, Y_Vz, Z_Vz);


			}
		}
	}

	// Vz
	for (k=zsp; k<zsp+zmp; k++){
		for (j=ysp; j<ysp+ymp-1; j++){
			for(i=xsp; i<xsp+xmp-1; i++){

				/* Extract the coordinates for Vx, Vy & Vz (as they are shifted we have to reconstruct them here) */
				X_Vz =  ( coords[k  ][j  ][i  ].x + coords[k  ][j+1][i  ].x + coords[k  ][j  ][i+1].x + coords[k  ][j+1][i+1].x )/4.0;
				Y_Vz =  ( coords[k  ][j  ][i  ].y + coords[k  ][j+1][i  ].y + coords[k  ][j  ][i+1].y + coords[k  ][j+1][i+1].y )/4.0;
				Z_Vz =  ( coords[k  ][j  ][i  ].z + coords[k  ][j+1][i  ].z + coords[k  ][j  ][i+1].z + coords[k  ][j+1][i+1].z )/4.0;

				//===========================================================================
				if (test_equation==0){
					velocity[k][j][i].Vz = 3; //sin(TwoPi/W*X_Vz)*cos(TwoPi/L*Y_Vz)*sin(TwoPi/H*Z_Vz)/2.0;
				}
				else if (test_equation==1){
					velocity[k][j][i].Vz = sin(TwoPi/W*X_Vz)*cos(TwoPi/L*Y_Vz)*sin(TwoPi/H*Z_Vz)/2.0;
				}
				else if (test_equation==2){
					velocity[k][j][i].Vz = -sin(TwoPi/W*X_Vz)*cos(TwoPi/L*Y_Vz)*sin(TwoPi/H*Z_Vz);
				}
				else if (test_equation==3){
					velocity[k][j][i].Vz = -sin(TwoPi/W*X_Vz)*cos(TwoPi/L*Y_Vz)*sin(TwoPi/H*Z_Vz);
				}
				else if (test_equation==4){
					velocity[k][j][i].Vz = -sin(TwoPi/W*X_Vz)*cos(TwoPi/L*Y_Vz)*sin(TwoPi/H*Z_Vz);
				}
				else if (test_equation==5){
					velocity[k][j][i].Vz = -sin(TwoPi/W*X_Vz)*cos(TwoPi/L*Y_Vz)*sin(TwoPi/H*Z_Vz);
				}
				//===========================================================================


//				PetscPrintf(PETSC_COMM_WORLD,"  coord = [%g,%g,%g] - [%g,%g,%g] - [%g,%g,%g] \n",X_Vx, Y_Vx, Z_Vx, X_Vy, Y_Vy, Z_Vy, X_Vz, Y_Vz, Z_Vz);


			}
		}
	}


	// Pressure
	ierr = DMDAGetCorners(user->DA_Processors,&xsp,&ysp,&zsp,&xmp,&ymp,&zmp);CHKERRQ(ierr);
	for (k=zsp; k<zsp+zmp; k++){
		for (j=ysp; j<ysp+ymp; j++){
			for(i=xsp; i<xsp+xmp; i++){

				/* Extract the coordinates for Vx, Vy & Vz (as they are shifted we have to reconstruct them here) */
				X_P  =  ( coords[k  ][j  ][i  ].x + coords[k  ][j+1][i  ].x + coords[k  ][j  ][i+1].x + coords[k  ][j+1][i+1].x +
						  coords[k+1][j  ][i  ].x + coords[k+1][j+1][i  ].x + coords[k+1][j  ][i+1].x + coords[k+1][j+1][i+1].x )/8.0;

				Y_P  =  ( coords[k  ][j  ][i  ].y + coords[k  ][j+1][i  ].y + coords[k  ][j  ][i+1].y + coords[k  ][j+1][i+1].y +
					      coords[k+1][j  ][i  ].y + coords[k+1][j+1][i  ].y + coords[k+1][j  ][i+1].y + coords[k+1][j+1][i+1].y )/8.0;
				Z_P  =  ( coords[k  ][j  ][i  ].z + coords[k  ][j+1][i  ].z + coords[k  ][j  ][i+1].z + coords[k  ][j+1][i+1].z +
						  coords[k+1][j  ][i  ].z + coords[k+1][j+1][i  ].z + coords[k+1][j  ][i+1].z + coords[k+1][j+1][i+1].z )/8.0;

				//===========================================================================
				if (test_equation==0){
					pressure[k][j][i] = 4;
				}
				else if (test_equation==2){
					pressure[k][j][i] =  -6.0*sin(TwoPi/W*X_P)*cos(TwoPi/L*Y_P)*cos(TwoPi/H*Z_P);
				}
				else if (test_equation==3){
					pressure[k][j][i] =  -6.0*sin(TwoPi/W*X_P)*cos(TwoPi/L*Y_P)*cos(TwoPi/H*Z_P);
				}
				else if (test_equation==4){
					pressure[k][j][i] =  -6.0*sin(TwoPi/W*X_P)*cos(TwoPi/L*Y_P)*cos(TwoPi/H*Z_P);
				}
				else if (test_equation==5){

					/* For 1th force balance equation: */
					//pressure[k][j][i] =  -4.0*exp(X_P)*sin(TwoPi/W*X_P)*cos(TwoPi/L*Y_P)*cos(TwoPi/H*Z_P) - 2.0*exp(X_P)*cos(TwoPi/W*X_P)*cos(TwoPi/L*Y_P)*cos(TwoPi/H*Z_P);	// exp(x)
					//pressure[k][j][i] =  -6.0*exp(Y_P)*sin(TwoPi/W*X_P)*cos(TwoPi/L*Y_P)*cos(TwoPi/H*Z_P) - 2.0*exp(Y_P)*sin(TwoPi/W*X_P)*sin(TwoPi/L*Y_P)*cos(TwoPi/H*Z_P);	// exp(y)
					pressure[k][j][i] =  -6.0*exp(Z_P)*sin(TwoPi/W*X_P)*cos(TwoPi/L*Y_P)*cos(TwoPi/H*Z_P) - 2.0*exp(Z_P)*sin(TwoPi/W*X_P)*cos(TwoPi/L*Y_P)*sin(TwoPi/H*Z_P);	// exp(z)
//					pressure[k][j][i] =  -4.0*exp(Z_P)*sin(TwoPi/W*X_P)*cos(TwoPi/L*Y_P)*cos(TwoPi/H*Z_P) ;	// exp(z)


					/* For 2nd force balance equation: */
					//pressure[k][j][i] =  2.0*exp(X_P)*cos(TwoPi/W*X_P)*cos(TwoPi/L*Y_P)*cos(TwoPi/H*Z_P) - 6.0*exp(X_P)*sin(TwoPi/W*X_P)*cos(TwoPi/L*Y_P)*cos(TwoPi/H*Z_P);	// exp(x)
					//pressure[k][j][i] =  -4.0*exp(Y_P)*sin(TwoPi/W*X_P)*cos(TwoPi/L*Y_P)*cos(TwoPi/H*Z_P) + 2.0*exp(Y_P)*sin(TwoPi/W*X_P)*sin(TwoPi/L*Y_P)*cos(TwoPi/H*Z_P);	// exp(y)
					//pressure[k][j][i] =  -6.0*exp(Z_P)*sin(TwoPi/W*X_P)*cos(TwoPi/L*Y_P)*cos(TwoPi/H*Z_P) - 2.0*exp(Z_P)*sin(TwoPi/W*X_P)*cos(TwoPi/L*Y_P)*sin(TwoPi/H*Z_P);	// exp(z)

					/* For 3rd force balance equation: */
					//pressure[k][j][i] =  -6.0*exp(X_P)*sin(TwoPi/W*X_P)*cos(TwoPi/L*Y_P)*cos(TwoPi/H*Z_P) + 2.0*exp(X_P)*cos(TwoPi/W*X_P)*cos(TwoPi/L*Y_P)*cos(TwoPi/H*Z_P);	// exp(x)
					//pressure[k][j][i] =  -6.0*exp(Y_P)*sin(TwoPi/W*X_P)*cos(TwoPi/L*Y_P)*cos(TwoPi/H*Z_P) - 2.0*exp(Y_P)*sin(TwoPi/W*X_P)*sin(TwoPi/L*Y_P)*cos(TwoPi/H*Z_P);	// exp(y)
					//pressure[k][j][i] =  -4.0*exp(Z_P)*sin(TwoPi/W*X_P)*cos(TwoPi/L*Y_P)*cos(TwoPi/H*Z_P) + 2.0*exp(Z_P)*sin(TwoPi/W*X_P)*cos(TwoPi/L*Y_P)*sin(TwoPi/H*Z_P);	// exp(z)

				}
				//===========================================================================

			}
		}
	}

	ierr = DMDAVecRestoreArray(cda,gc,&coords); 							CHKERRQ(ierr);
	//ierr = DMDestroy(cda); 												CHKERRQ(ierr);	//coordinates
	//ierr = VecDestroy(gc); 												CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(user->DA_Vel,	sol_vel	,&velocity); 	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(user->DA_Pres,	Pressure,&pressure); 	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Interpolates velocities to corner points as they are vertex-centered in the FDSTAG formulation
 * (although the DMDA is constructed using a ghost-point approach)
 *
 * BK, December 2010
 */
#undef __FUNCT__
#define __FUNCT__ "FDSTAG_InterpolateVelocityPressureToCornerpoints"
PetscErrorCode FDSTAG_InterpolateVelocityPressureToCornerpoints(UserContext *user, Vec sol_vel)
{
	PetscFunctionBegin;

	Vec 			sol_vel_local;
	Field 			***velocity, ***velocity_local;
	PetscErrorCode	ierr;
	PetscInt		xs,ys,zs,xm,ym,zm,i,j,k;

	/* Copy velocity vector to local vector  */
	ierr = 	DMGetLocalVector(user->DA_Vel,&sol_vel_local);	CHKERRQ(ierr);
	ierr = 	DMGlobalToLocalBegin(user->DA_Vel,sol_vel,INSERT_VALUES,sol_vel_local);	CHKERRQ(ierr);
	ierr = 	DMGlobalToLocalEnd  (user->DA_Vel,sol_vel,INSERT_VALUES,sol_vel_local);	CHKERRQ(ierr);


	/* Extract velocity arrays (global and local with ghost points) */
	ierr = 	DMDAVecGetArray(user->DA_Vel, 	sol_vel_local, 	&velocity_local); 	CHKERRQ(ierr);
	ierr = 	DMDAVecGetArray(user->DA_Vel,  sol_vel, 		&velocity); 		CHKERRQ(ierr);

	ierr = DMDAGetCorners(user->DA_Vel,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);


	//PetscPrintf(PETSC_COMM_WORLD,"FDSTAG_InterpolateVelocityPressureToCornerpoints: [xs,ys,zs],[xm,ym,zm]=[%lld,%lld,%lld]-[%lld,%lld,%lld] nnode=[%lld,%lld,%lld]\n",
	//		(LLD)xs,(LLD)ys,(LLD)zs,(LLD)xm,(LLD)ym,(LLD)zm,(LLD)(user->finest_nnode_x),(LLD)(user->finest_nnode_y),(LLD)(user->finest_nnode_z));


	for (k=zs; k<zs+zm; k++){
		for (j=ys; j<ys+ym; j++){
			for(i=xs; i<xs+xm; i++){


				/* Vx nodes*/
				if ((j>0) && (k>0)  && (j<user->finest_nnode_y-1) && (k<user->finest_nnode_z-1)){
					velocity[k][j][i].Vx = 	(	velocity_local[k-1][j-1][i  ].Vx + velocity_local[k-1][j  ][i  ].Vx
											+ 	velocity_local[k  ][j-1][i  ].Vx + velocity_local[k  ][j  ][i  ].Vx)/4.0;
				}
				else if ( ((j==0) ) && (k>0)   && (k<user->finest_nnode_z-1)){
					velocity[k][j][i].Vx = 	(	velocity_local[k-1][j  ][i  ].Vx + velocity_local[k  ][j  ][i  ].Vx )/2.0;
				}
				else if  ( ((j==user->finest_nnode_y-1) ) && (k>0)   && (k<user->finest_nnode_z-1)){  // added nodes because of DA
					velocity[k][j][i].Vx = 	(	velocity_local[k-1][j-1][i  ].Vx + velocity_local[k  ][j-1][i  ].Vx )/2.0;
				}
				else if ( ((k==0) ) && (j>0)   && (j<user->finest_nnode_y-1)){
					velocity[k][j][i].Vx = 	(	velocity_local[k][j-1][i  ].Vx + velocity_local[k  ][j  ][i  ].Vx )/2.0;
				}
				else if ( ((k==user->finest_nnode_z-1) ) && (j>0)   && (j<user->finest_nnode_y-1)){
					velocity[k][j][i].Vx = 	(	velocity_local[k-1][j-1][i  ].Vx + velocity_local[k-1][j  ][i  ].Vx )/2.0;
				}
				else if  ( ((k==0) ) && (j==0)   ){
					velocity[k][j][i].Vx = 	velocity_local[k  ][j  ][i  ].Vx;
				}
				else if  ( ((k==0) ) && (j==user->finest_nnode_y-1)   ){
					velocity[k][j][i].Vx = 	velocity_local[k  ][j-1][i  ].Vx;
				}
				else if  ( ((k==user->finest_nnode_z-1) ) && (j==0)   ){
					velocity[k][j][i].Vx = 	velocity_local[k-1][j  ][i  ].Vx;
				}
				else if  ( ((k==user->finest_nnode_z-1) ) && (j==user->finest_nnode_y-1)   ){
					velocity[k][j][i].Vx = 	velocity_local[k-1][j-1][i  ].Vx;
				}



				/* Vy nodes*/
				if ((i>0) && (k>0)  && (i<user->finest_nnode_x-1) && (k<user->finest_nnode_z-1)){
					velocity[k][j][i].Vy = 	(	velocity_local[k-1][j  ][i-1].Vy + velocity_local[k-1][j  ][i  ].Vy
											+ 	velocity_local[k  ][j  ][i-1].Vy + velocity_local[k  ][j  ][i  ].Vy)/4.0;
				}
				else if ( ((i==0) ) && (k>0)   && (k<user->finest_nnode_z-1)){
					velocity[k][j][i].Vy = 	(	velocity_local[k-1][j  ][i  ].Vy + velocity_local[k  ][j  ][i  ].Vy )/2.0;
				}
				else if ( ((i==user->finest_nnode_x-1) ) && (k>0)   && (k<user->finest_nnode_z-1)){
					velocity[k][j][i].Vy = 	(	velocity_local[k-1][j  ][i-1].Vy + velocity_local[k  ][j  ][i-1].Vy )/2.0;
				}
				else if ( ((k==0) ) && (i>0)   && (i<user->finest_nnode_x-1)){
					velocity[k][j][i].Vy = 	(	velocity_local[k  ][j  ][i-1].Vy + velocity_local[k  ][j  ][i  ].Vy )/2.0;
				}
				else if ( ((k==user->finest_nnode_z-1) ) && (i>0)   && (i<user->finest_nnode_x-1)){
					velocity[k][j][i].Vy = 	(	velocity_local[k-1][j  ][i-1].Vy + velocity_local[k-1][j  ][i  ].Vy )/2.0;
				}
				else if  ( ((k==0) ) && (i==0)   ){
					velocity[k][j][i].Vy = 	velocity_local[k  ][j  ][i  ].Vy;
				}
				else if  ( ((k==0) ) && (i==user->finest_nnode_x-1)   ){
					velocity[k][j][i].Vy = 	velocity_local[k  ][j  ][i -1].Vy;
				}
				else if  ( ((k==user->finest_nnode_z-1) ) && (i==0)   ){
					velocity[k][j][i].Vy = 	velocity_local[k-1][j  ][i  ].Vy;
				}
				else if  ( ((k==user->finest_nnode_z-1) ) && (i==user->finest_nnode_x-1)   ){
					velocity[k][j][i].Vy = 	velocity_local[k-1][j  ][i-1].Vy;
				}



				/* Vz nodes*/
				if ((i>0) && (j>0)  && (i<user->finest_nnode_x-1) && (j<user->finest_nnode_y-1)){
					velocity[k][j][i].Vz = 	(	velocity_local[k  ][j-1][i-1].Vz + velocity_local[k  ][j  ][i-1].Vz
											+ 	velocity_local[k  ][j-1][i  ].Vz + velocity_local[k  ][j  ][i  ].Vz)/4.0;
				}
				else if ( ((j==0) ) && (i>0)   && (i<user->finest_nnode_x-1)){
					velocity[k][j][i].Vz = 	(	velocity_local[k ][j  ][i-1].Vz + velocity_local[k  ][j  ][i  ].Vz )/2.0;
				}
				else if ( ((j==user->finest_nnode_y-1) ) && (i>0)   && (i<user->finest_nnode_x-1)){
					velocity[k][j][i].Vz = 	(	velocity_local[k ][j-1][i-1].Vz + velocity_local[k  ][j-1][i  ].Vz )/2.0;
				}
				else if ( ((i==0) ) && (j>0)   && (j<user->finest_nnode_y-1)){
					velocity[k][j][i].Vz = 	(	velocity_local[k ][j-1][i  ].Vz + velocity_local[k  ][j  ][i  ].Vz )/2.0;
				}
				else if ( ((i==user->finest_nnode_x-1) ) && (j>0)   && (j<user->finest_nnode_y-1)){
					velocity[k][j][i].Vz = 	(	velocity_local[k ][j-1][i-1].Vz + velocity_local[k  ][j  ][i-1].Vz )/2.0;
				}
				else if  ( ((i==0) ) && (j==0)   ){
					velocity[k][j][i].Vz = 	velocity_local[k  ][j  ][i  ].Vz;
				}
				else if  ( ((i==0) ) && (j==user->finest_nnode_y-1)   ){
					velocity[k][j][i].Vz = 	velocity_local[k  ][j-1][i ].Vz;
				}
				else if  ( ((i==user->finest_nnode_x-1) ) && (j==0)   ){
					velocity[k][j][i].Vz = 	velocity_local[k  ][j  ][i-1].Vz;
				}
				else if  ( ((i==user->finest_nnode_x-1) ) && (j==user->finest_nnode_y-1)   ){
					velocity[k][j][i].Vz = 	velocity_local[k  ][j-1][i-1].Vz;
				}
			}
		}
	}

	ierr = 	DMDAVecRestoreArray(user->DA_Vel, sol_vel, &velocity); CHKERRQ(ierr);
	ierr = 	DMDAVecRestoreArray(user->DA_Vel, sol_vel_local, &velocity_local); CHKERRQ(ierr);
	ierr =	DMRestoreLocalVector(user->DA_Vel,	&sol_vel_local);	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Sets a viscosity field on center & corner nodes to test the FDSTAG scheme.
 *
 * BK, November 2010
 */
#undef __FUNCT__
#define __FUNCT__ "FDSTAG_SetInitialViscosityFields"
PetscErrorCode FDSTAG_SetInitialViscosityFields(UserContext *user)
{
	PetscFunctionBegin;
	DM				cda;
	PetscErrorCode	ierr;
	PetscInt		i,j,k,zsp,zmp,ysp,ymp,xsp,xmp,xs,ys,zs,xm,ym,zm;
	PetscInt		ViscosityStructure;
	PetscScalar		X,Y,Z;
	Vec			 	gc;
	DMDACoor3d		***coords;
	PetscScalar 	***viscosity_center;
	PetscScalar		***viscosity_XY,***viscosity_XZ,***viscosity_YZ;



	ViscosityStructure = 1;	//0-isoviscous (1); 1-exp(X)

	/* Obtain coordinates */
	ierr = DMGetCoordinateDM(user->DA_Vel,&cda); 		CHKERRQ(ierr);	//coordinates
	ierr = DMGetCoordinatesLocal(user->DA_Vel,&gc); 	CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,gc,&coords); 	CHKERRQ(ierr);

	/* Set viscosity @ cell centers */
	ierr = DMDAVecGetArray(user->FDSTAG.DA_CENTER,   			   user->FDSTAG.Center_EffectiveViscosity, 	&viscosity_center	 ); CHKERRQ(ierr); // Viscosity @ center
	ierr = DMDAGetCorners(user->FDSTAG.DA_CENTER,&xsp,&ysp,&zsp,&xmp,&ymp,&zmp); CHKERRQ(ierr);
	for (k=zsp; k<zsp+zmp; k++){
		for (j=ysp; j<ysp+ymp; j++){
			for(i=xsp; i<xsp+xmp; i++){

				/* Extract the coordinates for Vx, Vy & Vz (as they are shifted we have to reconstruct them here) */
				X  =  ( coords[k  ][j  ][i  ].x + coords[k  ][j+1][i  ].x + coords[k  ][j  ][i+1].x + coords[k  ][j+1][i+1].x +
						  coords[k+1][j  ][i  ].x + coords[k+1][j+1][i  ].x + coords[k+1][j  ][i+1].x + coords[k+1][j+1][i+1].x )/8.0;

				Y  =  ( coords[k  ][j  ][i  ].y + coords[k  ][j+1][i  ].y + coords[k  ][j  ][i+1].y + coords[k  ][j+1][i+1].y +
						  coords[k+1][j  ][i  ].y + coords[k+1][j+1][i  ].y + coords[k+1][j  ][i+1].y + coords[k+1][j+1][i+1].y )/8.0;
				Z  =  ( coords[k  ][j  ][i  ].z + coords[k  ][j+1][i  ].z + coords[k  ][j  ][i+1].z + coords[k  ][j+1][i+1].z +
						  coords[k+1][j  ][i  ].z + coords[k+1][j+1][i  ].z + coords[k+1][j  ][i+1].z + coords[k+1][j+1][i+1].z )/8.0;

				//===========================================================================
				if (ViscosityStructure==0){
					viscosity_center[k][j][i] 	= 1.0;
				}
				else if (ViscosityStructure==1){
					viscosity_center[k][j][i] 	= exp(X);
					viscosity_center[k][j][i] 	= exp(Y);
					viscosity_center[k][j][i] 	= exp(Z);

				}
				//===========================================================================

			}
		}
	}
	ierr = DMDAVecRestoreArray(user->FDSTAG.DA_CENTER,   				user->FDSTAG.Center_EffectiveViscosity, 		&viscosity_center); 				CHKERRQ(ierr);


	/* Set viscosity @ Sxy points */
	ierr = DMDAVecGetArray(user->FDSTAG.DA_XY_POINTS,   user->FDSTAG.XYPoints_EffectiveViscosity, 	&viscosity_XY		 ); CHKERRQ(ierr); // Viscosity @ Sxy coordinates
	ierr = DMDAGetCorners(user->FDSTAG.DA_XY_POINTS,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);		// we loop over corner points
	for (k=zs; k<zs+zm; k++){
		for (j=ys; j<ys+ym; j++){
			for (i=xs; i<xs+xm; i++){

				/* coordinate of Sxy point */
				X  =    coords[k  ][j  ][i  ].x ;
				Y  =    coords[k  ][j  ][i  ].y ;
				Z  =  ( coords[k  ][j  ][i  ].z + coords[k+1][j  ][i  ].z )/2.0;

				//===========================================================================
				if (ViscosityStructure==0){
					viscosity_XY[k][j][i] 	= 1.0;
				}
				else if (ViscosityStructure==1){
					viscosity_XY[k][j][i] 	= exp(X);
					viscosity_XY[k][j][i] 	= exp(Y);
					viscosity_XY[k][j][i] 	= exp(Z);
				}
				//===========================================================================
			}
		}
	}
	ierr = DMDAVecRestoreArray(user->FDSTAG.DA_XY_POINTS,    	user->FDSTAG.XYPoints_EffectiveViscosity, 		&viscosity_XY		); CHKERRQ(ierr);

	/* Set viscosity @ Sxz points */
	ierr = DMDAVecGetArray(user->FDSTAG.DA_XZ_POINTS,   user->FDSTAG.XZPoints_EffectiveViscosity, 	&viscosity_XZ		 ); CHKERRQ(ierr); // Viscosity @ Sxz coordinates
	ierr = DMDAGetCorners(user->FDSTAG.DA_XZ_POINTS,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);		// we loop over corner points
	PetscPrintf(PETSC_COMM_WORLD,"[xs,ys,zs,xm,ym,zm]=[%lld,%lld,%lld,%lld,%lld,%lld]  \n", (LLD)xs,(LLD)ys,(LLD)zs,(LLD)xm,(LLD)ym,(LLD)zm);
	for (k=zs; k<zs+zm; k++){
		for (j=ys; j<ys+ym; j++){
			for (i=xs; i<xs+xm; i++){

				/* coordinate of Sxy point */
				X  =    coords[k  ][j  ][i  ].x ;
				Z  =    coords[k  ][j  ][i  ].z ;
				Y  =  ( coords[k  ][j  ][i  ].y + coords[k  ][j+1][i  ].y )/2.0;

				//===========================================================================
				if (ViscosityStructure==0){
					viscosity_XZ[k][j][i] 	= 1.0;
				}
				else if (ViscosityStructure==1){
					viscosity_XZ[k][j][i] 	= exp(X);
					viscosity_XZ[k][j][i] 	= exp(Y);
					viscosity_XZ[k][j][i] 	= exp(Z);
				}
				//===========================================================================
			}
		}
	}
	ierr = DMDAVecRestoreArray(user->FDSTAG.DA_XZ_POINTS,   user->FDSTAG.XZPoints_EffectiveViscosity, 		&viscosity_XZ		); CHKERRQ(ierr);


	/* Set viscosity @ Syz points */
	ierr = DMDAVecGetArray(user->FDSTAG.DA_YZ_POINTS,   user->FDSTAG.YZPoints_EffectiveViscosity, 	&viscosity_YZ		 ); CHKERRQ(ierr); // Viscosity @ Syz coordinates
	ierr = DMDAGetCorners(user->FDSTAG.DA_YZ_POINTS,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);		// we loop over corner points
	for (k=zs; k<zs+zm; k++){
		for (j=ys; j<ys+ym; j++){
			for (i=xs; i<xs+xm; i++){

				/* coordinate of Sxy point */
				Y  =    coords[k  ][j  ][i  ].y ;
				Z  =    coords[k  ][j  ][i  ].z ;
				X  =  ( coords[k  ][j  ][i  ].x + coords[k  ][j  ][i+1].x )/2.0;

				//===========================================================================
				if (ViscosityStructure==0){
					viscosity_YZ[k][j][i] 	= 1.0;
				}
				else if (ViscosityStructure==1){
					//viscosity_YZ[k][j][i] 	= exp(X);
					//viscosity_YZ[k][j][i] 	= exp(Y);
					viscosity_YZ[k][j][i] 	= exp(Z);
				}
				//===========================================================================
			}
		}
	}
	ierr = DMDAVecRestoreArray(user->FDSTAG.DA_YZ_POINTS,   user->FDSTAG.YZPoints_EffectiveViscosity, 		&viscosity_YZ		); CHKERRQ(ierr);



	ierr = DMDAVecRestoreArray(cda,gc,&coords); 							CHKERRQ(ierr);
	//ierr = DMDestroy(cda); 												CHKERRQ(ierr);	//coordinates
	//ierr = VecDestroy(gc); 												CHKERRQ(ierr);




	/* Dump viscosity top disk (just to check) */
	{
		PetscViewer	vv;

		/* Write center viscosity and density */
		ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "Viscosity_Center_FDSTAG.vtk", &vv);CHKERRQ(ierr);
		ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
		ierr = PetscObjectSetName( (PetscObject)user->FDSTAG.Center_EffectiveViscosity, "Center_Viscosity" );CHKERRQ(ierr);
		ierr = PetscObjectSetName( (PetscObject)user->FDSTAG.Center_Density, "Center_Density" );CHKERRQ(ierr);
		ierr = DMView(user->FDSTAG.DA_CENTER, vv);CHKERRQ(ierr);
		ierr = VecView(user->FDSTAG.Center_EffectiveViscosity , vv);CHKERRQ(ierr);
		ierr = VecView(user->FDSTAG.Center_Density, vv);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&vv);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"FDSTAG DEBUG: wrote center viscosity & density info to  Viscosity_Center_FDSTAG.vtk \n");

		/* Write viscosity at Sxy points */
		ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "Viscosity_XY_FDSTAG.vtk", &vv);CHKERRQ(ierr);
		ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
		ierr = PetscObjectSetName( (PetscObject)user->FDSTAG.XYPoints_EffectiveViscosity, "Viscosity_XY" );CHKERRQ(ierr);
		ierr = DMView(user->FDSTAG.DA_XY_POINTS, vv);CHKERRQ(ierr);
		ierr = VecView(user->FDSTAG.XYPoints_EffectiveViscosity, vv);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&vv);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"FDSTAG DEBUG: wrote corner viscosity info to  Viscosity_XY_FDSTAG.vtk \n");

		/* Write viscosity at Sxz points */
		ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "Viscosity_XZ_FDSTAG.vtk", &vv);CHKERRQ(ierr);
		ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
		ierr = PetscObjectSetName( (PetscObject)user->FDSTAG.XZPoints_EffectiveViscosity, "XZ_Viscosity" );CHKERRQ(ierr);
		ierr = DMView(user->FDSTAG.DA_XZ_POINTS, vv);CHKERRQ(ierr);
		ierr = VecView(user->FDSTAG.XZPoints_EffectiveViscosity, vv);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&vv);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"FDSTAG DEBUG: wrote corner viscosity info to  Viscosity_XZ_FDSTAG.vtk \n");

		/* Write viscosity at Syz points */
		ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "Viscosity_YZ_FDSTAG.vtk", &vv);CHKERRQ(ierr);
		ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
		ierr = PetscObjectSetName( (PetscObject)user->FDSTAG.YZPoints_EffectiveViscosity, "YZ_Viscosity" );CHKERRQ(ierr);
		ierr = DMView(user->FDSTAG.DA_YZ_POINTS, vv);CHKERRQ(ierr);
		ierr = VecView(user->FDSTAG.XZPoints_EffectiveViscosity, vv);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&vv);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"FDSTAG DEBUG: wrote corner viscosity info to  Viscosity_YZ_FDSTAG.vtk \n");

	}



PetscFunctionReturn(0);
}
/*==========================================================================================================*/

/*==========================================================================================================*/
/* Compute the temperature stiffness matrix using FDSTAG
 *
 * rho*Cp*(dT/dt) = grad(k*div(T)) + Q
 *
 */
#undef __FUNCT__
#define __FUNCT__ "ComputeStiffnessMatrixRHSTemperature_FDSTAG"
PetscErrorCode ComputeStiffnessMatrixRHSTemperature_FDSTAG(DM da_temp, DM da, Mat T_MAT, UserContext *user,
		Vec Temp, Vec rhs_Temp_local, Vec rhs_Temp, PetscScalar dt )
{

	PetscErrorCode ierr;
	PetscInt			nnode_x, nnode_y, nnode_z, ix,iy,iz,xs,ys,zs,xm,ym,zm;
	PetscScalar			dx_1,dx1,dy_1,dy1,dz_1,dz1;
	PetscScalar 		kx_1, kx1, ky_1, ky1, kz_1, kz1, rho_cp, Q, v[13], dx_c, dy_c, dz_c;
	PetscScalar			***Density_Nodes, ***Cp_Nodes, ***Q_Nodes, ***k_Nodes, ***rhs, ***T_old;
	DM			 		cda;
	Vec			 		gc, Density_Nodes_local, Cp_Nodes_local, k_Nodes_local, Q_Nodes_local;
	DMDACoor3d		 	***coords;
	MatStencil     		row[1], col[13];


	ierr = DMGetCoordinateDM(da,&cda); 		CHKERRQ(ierr);	//coordinates
	ierr = DMGetCoordinatesLocal(da,&gc); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,gc,&coords); 	CHKERRQ(ierr);

	ierr = MatZeroEntries(T_MAT); CHKERRQ(ierr);

	/* print some info */
	ierr = DMDAGetInfo(da_temp,						0,&nnode_x,  &nnode_y,	&nnode_z,	0,0,0,0,0,0,0,0,0);		CHKERRQ(ierr);  // # of nodes 	 in all directions
	PetscPrintf(PETSC_COMM_WORLD,"#  Forming FD stiffness matrix for energy equation with mx,my,mz=(%lld, %lld, %lld) with Top/Bottom T=[%g,%g]...  \n",(LLD)nnode_x,(LLD)nnode_y,(LLD)nnode_z,  user->Temp_top,  user->Temp_bottom);


	/* Get density @ nodal points, including ghost node values */
	ierr = DMGetLocalVector(user->FDSTAG.DA_CORNER,&Density_Nodes_local);	CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(user->FDSTAG.DA_CORNER,user->FDSTAG.Corner_Density,INSERT_VALUES,Density_Nodes_local);	CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd  (user->FDSTAG.DA_CORNER,user->FDSTAG.Corner_Density,INSERT_VALUES,Density_Nodes_local);	CHKERRQ(ierr);

	ierr = DMGetLocalVector(user->FDSTAG.DA_CORNER,&Cp_Nodes_local);	CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(user->FDSTAG.DA_CORNER,user->FDSTAG.Corner_HeatCapacity,INSERT_VALUES,Cp_Nodes_local);	CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd  (user->FDSTAG.DA_CORNER,user->FDSTAG.Corner_HeatCapacity,INSERT_VALUES,Cp_Nodes_local);	CHKERRQ(ierr);

	ierr = DMGetLocalVector(user->FDSTAG.DA_CORNER,&Q_Nodes_local);	CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(user->FDSTAG.DA_CORNER,user->FDSTAG.Corner_RadioactiveHeat,INSERT_VALUES,Q_Nodes_local);	CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd  (user->FDSTAG.DA_CORNER,user->FDSTAG.Corner_RadioactiveHeat,INSERT_VALUES,Q_Nodes_local);	CHKERRQ(ierr);

	ierr = DMGetLocalVector(user->FDSTAG.DA_CORNER,&k_Nodes_local);	CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(user->FDSTAG.DA_CORNER,user->FDSTAG.Corner_Conductivity,INSERT_VALUES,k_Nodes_local);	CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd  (user->FDSTAG.DA_CORNER,user->FDSTAG.Corner_Conductivity,INSERT_VALUES,k_Nodes_local);	CHKERRQ(ierr);

	/* Get arrays for material properties @ nodes */
	ierr = DMDAVecGetArray(user->FDSTAG.DA_CORNER,    Density_Nodes_local, 			&Density_Nodes	); CHKERRQ(ierr); // density
	ierr = DMDAVecGetArray(user->FDSTAG.DA_CORNER,    Cp_Nodes_local, 					&Cp_Nodes		); CHKERRQ(ierr); // heat capacity
	ierr = DMDAVecGetArray(user->FDSTAG.DA_CORNER,    k_Nodes_local, 					&k_Nodes		); CHKERRQ(ierr); // conductivity
	ierr = DMDAVecGetArray(user->FDSTAG.DA_CORNER,    Q_Nodes_local, 					&Q_Nodes		); CHKERRQ(ierr); // radioactive heat production

	ierr = DMDAVecGetArray(da_temp,    				Temp, 							&T_old		); CHKERRQ(ierr); // old temperature

	/* Rhs */
	VecSet(rhs_Temp,0.0);  	VecSet(rhs_Temp_local,0.0);
	DMGlobalToLocalBegin(da_temp, rhs_Temp, INSERT_VALUES, rhs_Temp_local);
	DMGlobalToLocalEnd(da_temp,   rhs_Temp, INSERT_VALUES, rhs_Temp_local);
	ierr = DMDAVecGetArray(da_temp, rhs_Temp_local,		 &rhs);		CHKERRQ(ierr);


	/* Corners */
	ierr = DMDAGetCorners(da,      			      &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);

	/* Make a loop over all local nodes, construct the element stiffness matrix */
	for (iz=zs; iz<zs+zm; iz++){
		for (iy=ys; iy<ys+ym; iy++){
			for(ix=xs; ix<xs+xm; ix++){

				if ((ix>0) && (ix<user->finest_nnode_x-1) && (iy>0) && (iy<user->finest_nnode_y-1) && (iz>0) && (iz<user->finest_nnode_z-1)){
					/* Central nodes */

					/* Compute spacing */
					dx_1 	= coords[iz  ][iy  ][ix  ].x-coords[iz  ][iy  ][ix-1].x;
					dx1  	= coords[iz  ][iy  ][ix+1].x-coords[iz  ][iy  ][ix  ].x;
					dx_c 	= (dx_1 + dx1)/2;

					dy_1 	= coords[iz  ][iy  ][ix  ].y-coords[iz  ][iy-1][ix  ].y;
					dy1 	= coords[iz  ][iy+1][ix  ].y-coords[iz  ][iy  ][ix  ].y;
					dy_c 	= (dy_1 + dy1)/2;

					dz_1 	= coords[iz  ][iy  ][ix  ].z-coords[iz-1][iy  ][ix  ].z;
					dz1 	= coords[iz+1][iy  ][ix  ].z-coords[iz  ][iy  ][ix  ].z;
					dz_c 	= (dz_1 + dz1)/2;

					/* Material properties */
					rho_cp 	=	 Density_Nodes[iz  ][iy  ][ix  ]*Cp_Nodes[iz  ][iy  ][ix  ];  	// rho*cp at (ix,iy,iz)
					Q 		=	 Q_Nodes[iz  ][iy  ][ix  ];  									// Q at (ix    ,iy    ,iz    )
					kx1 	=	 (k_Nodes[iz  ][iy  ][ix+1] + k_Nodes[iz  ][iy  ][ix ])/2.0;	// k at (ix+1/2,iy    ,iz    )
					kx_1 	=	 (k_Nodes[iz  ][iy  ][ix-1] + k_Nodes[iz  ][iy  ][ix ])/2.0;	// k at (ix-1/2,iy    ,iz    )
					ky1 	=	 (k_Nodes[iz  ][iy+1][ix  ] + k_Nodes[iz  ][iy  ][ix ])/2.0;	// k at (ix    ,iy+1/2,iz    )
					ky_1 	=	 (k_Nodes[iz  ][iy-1][ix  ] + k_Nodes[iz  ][iy  ][ix ])/2.0;	// k at (ix    ,iy-1/2,iz    )
					kz1 	=	 (k_Nodes[iz+1][iy  ][ix  ] + k_Nodes[iz  ][iy  ][ix ])/2.0;	// k at (ix    ,iy    ,iz+1/2)
					kz_1 	=	 (k_Nodes[iz-1][iy  ][ix  ] + k_Nodes[iz  ][iy  ][ix ])/2.0;	// k at (ix    ,iy    ,iz-1/2)

					/* Form stencil */
					row[0].i = ix; row[0].j = iy; row[0].k = iz; row[0].c = 0;		// equation number

					// d/dx(k*dT/dx)
					col[0].i = ix+1; col[0].j = iy; col[0].k = iz;		v[0] =  kx1/dx1/dx_c;
					col[1].i = ix  ; col[1].j = iy; col[1].k = iz;		v[1] = -kx1/dx1/dx_c;

					col[2].i = ix  ; col[2].j = iy; col[2].k = iz;		v[2] =  -kx_1/dx_1/dx_c;
					col[3].i = ix-1; col[3].j = iy; col[3].k = iz;		v[3] =   kx_1/dx_1/dx_c;

					// d/dy(k*dT/dy)
					col[4].i = ix  ; col[4].j = iy+1; col[4].k = iz;	v[4] =  ky1/dy1/dy_c;
					col[5].i = ix  ; col[5].j = iy  ; col[5].k = iz;	v[5] = -ky1/dy1/dy_c;

					col[6].i = ix  ; col[6].j = iy  ; col[6].k = iz;	v[6] =  -ky_1/dy_1/dy_c;
					col[7].i = ix  ; col[7].j = iy-1; col[7].k = iz;	v[7] =   ky_1/dy_1/dy_c;

					// d/dz(k*dT/dz)
					col[8].i = ix  ; col[8].j = iy  ; col[8].k = iz+1;	v[8] =  kz1/dz1/dz_c;
					col[9].i = ix  ; col[9].j = iy  ; col[9].k = iz  ;	v[9] = -kz1/dz1/dz_c;

					col[10].i= ix  ; col[10].j= iy  ; col[10].k= iz  ;	v[10]=  -kz_1/dz_1/dz_c;
					col[11].i= ix  ; col[11].j= iy  ; col[11].k= iz-1;	v[11]=   kz_1/dz_1/dz_c;

					// -rho*Cp/dt*T
					col[12].i= ix  ; col[12].j= iy  ; col[12].k= iz  ;	v[12]= -rho_cp/dt;


					// Set values in matrix
					ierr = MatSetValuesStencil(T_MAT,1,row,13,col,v,ADD_VALUES);		CHKERRQ(ierr);

					// Set values in RHS
					rhs[iz  ][iy  ][ix  ] = -rho_cp/dt*T_old[iz  ][iy  ][ix  ] - Q;

				}

				/* Set BC's. Note: we do it in a primitive, first order, manner here (easier to code). Could be changed in future */
				if ((ix==0) && (iz>0) && (iz<user->finest_nnode_z-1) && (iy>0) && (iy<user->finest_nnode_y-1)){
					// left boundary
					row[0].i = ix; row[0].j = iy; row[0].k = iz; row[0].c = 0;		// equation number

					// dT/dx=0
					col[0].i = ix+1; col[0].j = iy; col[0].k = iz;		v[0] = 1;
					col[1].i = ix  ; col[1].j = iy; col[1].k = iz;		v[1] = -1;

					// Set values in matrix
					ierr = MatSetValuesStencil(T_MAT,1,row,2,col,v,ADD_VALUES);		CHKERRQ(ierr);

					// Set values in RHS
					rhs[iz  ][iy  ][ix  ] = 0;

				}
				else if ((ix==user->finest_nnode_x-1) && (iz>0) && (iz<user->finest_nnode_z-1) && (iy>0) && (iy<user->finest_nnode_y-1)){
					// right boundary
					row[0].i = ix; row[0].j = iy; row[0].k = iz; row[0].c = 0;		// equation number

					// dT/dx=0
					col[0].i = ix-1; col[0].j = iy; col[0].k = iz;		v[0] = -1;
					col[1].i = ix  ; col[1].j = iy; col[1].k = iz;		v[1] =  1;

					// Set values in matrix
					ierr = MatSetValuesStencil(T_MAT,1,row,2,col,v,ADD_VALUES);		CHKERRQ(ierr);

					// Set values in RHS
					rhs[iz  ][iy  ][ix  ] = 0;
				}
				else if ((iy==0) && (iz>0) && (iz<user->finest_nnode_z-1) ){
					// front boundary
					row[0].i = ix; row[0].j = iy; row[0].k = iz; row[0].c = 0;		// equation number

					// dT/dy=0
					col[0].i = ix  ; col[0].j = iy+1; col[0].k = iz;		v[0] = 1;
					col[1].i = ix  ; col[1].j = iy; col[1].k = iz;		v[1] = -1;

					// Set values in matrix
					ierr = MatSetValuesStencil(T_MAT,1,row,2,col,v,ADD_VALUES);		CHKERRQ(ierr);

					// Set values in RHS
					rhs[iz  ][iy  ][ix  ] = 0;

				}
				else if ((iy==user->finest_nnode_y-1) && (iz>0) && (iz<user->finest_nnode_z-1)){
					// back boundary
					row[0].i = ix; row[0].j = iy; row[0].k = iz; row[0].c = 0;		// equation number

					// dT/dx=0
					col[0].i = ix  ; col[0].j = iy-1; col[0].k = iz;	v[0] = -1;
					col[1].i = ix  ; col[1].j = iy; col[1].k = iz;		v[1] =  1;

					// Set values in matrix
					ierr = MatSetValuesStencil(T_MAT,1,row,2,col,v,ADD_VALUES);		CHKERRQ(ierr);

					// Set values in RHS
					rhs[iz  ][iy  ][ix  ] = 0;
				}

				if (iz==0) {
					// Bottom boundary
					row[0].i = ix; row[0].j = iy; row[0].k = iz; row[0].c = 0;		// equation number
					col[0].i = ix; col[0].j = iy; col[0].k = iz;     v[0] = 1;

					// Set values in matrix
					ierr = MatSetValuesStencil(T_MAT,1,row,1,col,v,ADD_VALUES);		CHKERRQ(ierr);

					// Set values in RHS
					rhs[iz  ][iy  ][ix  ] = user->Temp_bottom;
				}
				else if (iz==user->finest_nnode_z-1) {
					// Bottom boundary
					row[0].i = ix; row[0].j = iy; row[0].k = iz; row[0].c = 0;		// equation number
					col[0].i = ix; col[0].j = iy; col[0].k = iz;     v[0] = 1;

					// Set values in matrix
					ierr = MatSetValuesStencil(T_MAT,1,row,1,col,v,ADD_VALUES);		CHKERRQ(ierr);

					// Set values in RHS
					rhs[iz  ][iy  ][ix  ] = user->Temp_top;
				}


			}
		}
	}


	//coordinates
	ierr = DMDAVecRestoreArray(cda,gc,&coords); CHKERRQ(ierr);
	//ierr = DMDestroy(cda); CHKERRQ(ierr);
	//ierr = VecDestroy(gc); CHKERRQ(ierr);



	/* Finalize the matrixes */
	ierr = MatAssemblyBegin(T_MAT,	MAT_FINAL_ASSEMBLY);	CHKERRQ(ierr);	// Assemble global temperature stiffness matrix
	ierr = MatAssemblyEnd(T_MAT,	MAT_FINAL_ASSEMBLY);	CHKERRQ(ierr);  // Assemble global temperature stiffness matrix


	//
	ierr = DMDAVecRestoreArray(user->FDSTAG.DA_CORNER,    	Density_Nodes_local, 	&Density_Nodes	); CHKERRQ(ierr); 	// density
	ierr = DMDAVecRestoreArray(user->FDSTAG.DA_CORNER,    	Cp_Nodes_local, 		&Cp_Nodes		); CHKERRQ(ierr); 	// heat capacity
	ierr = DMDAVecRestoreArray(user->FDSTAG.DA_CORNER,    	k_Nodes_local, 			&k_Nodes		); CHKERRQ(ierr); 	// conductivity
	ierr = DMDAVecRestoreArray(user->FDSTAG.DA_CORNER,    	Q_Nodes_local, 			&Q_Nodes		); CHKERRQ(ierr); 	// radioactive heat generation
	ierr = DMDAVecRestoreArray(da_temp,					rhs_Temp_local,			&rhs			); CHKERRQ(ierr);	// right-hand-side vector
	ierr = DMDAVecRestoreArray(da_temp,    				Temp, 					&T_old			); CHKERRQ(ierr); 	// old temperature


	ierr = DMLocalToGlobalBegin(da_temp,rhs_Temp_local, ADD_VALUES, rhs_Temp);				CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da_temp,rhs_Temp_local, ADD_VALUES, rhs_Temp);					CHKERRQ(ierr);



	ierr = DMRestoreLocalVector(user->FDSTAG.DA_CORNER,&Density_Nodes_local);	CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(user->FDSTAG.DA_CORNER,&Cp_Nodes_local);			CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(user->FDSTAG.DA_CORNER,&Q_Nodes_local);			CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(user->FDSTAG.DA_CORNER,&k_Nodes_local);			CHKERRQ(ierr);



	PetscFunctionReturn(0);

}
/*==========================================================================================================*/




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

Utils.c, contains the following subroutines:

LaMEMMod							-
Mod									-	Simple routine to compute result = Mod(x,y)
GetGlobalIndex						-	Get global row number, given i,j and k
StencilToLocalNumbering				-	Takes a stencil & creates an array with local numbering
ReadInputFile						-	Read input file if required
InitializeCode						-	Set default code parameters and read input file, if required
ComputeGlobalProperties				-	Compute global properties that are interesting to record
DAUpdateCoordinates3d_only_local	-	Update coordinates in DA, but not the ghost coordinates
DAUpdateCoordinates3d_with_ghosts	-	Update coordinates in DA, with ghost coordinates
PetscErrorCode DASetGhostedCoordinates  -       Updates the ghosted nodes

$Id$
$Date$

This version is compatible with PETSc version 3.3-p2
The code cannot be distributed, used for commercial purposes or handed on,
without the explicit agreement of Boris Kaus.
*/


#include "LaMEM.h"
#include "Utils.h"
#include "Parsing.h"
#include "ParaViewOutput.h"
#include "Material.h"
#include "NonDimensionalisation.h"

#define __DACompareStructures(var,v1,v2)PetscSynchronizedPrintf( comm, "[%lld]: DA->%s are different { %lld , %lld } \n", rank, (var), (LLD)((v1)),(LLD)((v2)) )

#undef __FUNCT__
#define __FUNCT__ "DACompareStructures"
PetscErrorCode DACompareStructures( DM da1, DM da2, PetscBool *flg )
{
	PetscErrorCode ierr;
	PetscInt si1,sj1,sk1 , si2,sj2,sk2;
	PetscInt mx1,my1,mz1 , mx2,my2,mz2;
	PetscInt M1,N1,P1 , M2,N2,P2;
	PetscInt cx1,cy1,cz1 , cx2,cy2,cz2;
	PetscInt dim1 , dim2;
	PetscInt sw1 , sw2;
	DMBoundaryType wrap1[3] , wrap2[3];
	DMDAStencilType st1 , st2;
	MPI_Comm comm;
	PetscMPIInt rank;


	*flg = PETSC_TRUE;


	ierr = DMDAGetInfo( da1, &dim1, &M1,&N1,&P1, &cx1,&cy1,&cz1, 0, &sw1, &wrap1[0],&wrap1[1],&wrap1[2], &st1 ); CHKERRQ(ierr);
	ierr = DMDAGetInfo( da2, &dim2, &M2,&N2,&P2, &cx2,&cy2,&cz2, 0, &sw2, &wrap2[0],&wrap2[1],&wrap2[2], &st2 ); CHKERRQ(ierr);
	
	ierr = DMDAGetCorners( da1, &si1,&sj1,&sk1 , &mx1,&my1,&mz1 ); CHKERRQ(ierr);
	ierr = DMDAGetCorners( da2, &si2,&sj2,&sk2 , &mx2,&my2,&mz2 ); CHKERRQ(ierr);


	PetscObjectGetComm( (PetscObject)da1, &comm );
	MPI_Comm_rank( comm, &rank );

	if( dim1 != dim2 ) {	__DACompareStructures( "dim",dim1,dim2 );		*flg = PETSC_FALSE;	}

	if( M1 != M2 ) {	__DACompareStructures( "M",M1,M2 );		*flg = PETSC_FALSE;	}
	if( N1 != N2 ) {	__DACompareStructures( "N",N1,N2 );		*flg = PETSC_FALSE;	}
	if( P1 != P2 ) {	__DACompareStructures( "P",P1,P2 );		*flg = PETSC_FALSE;	}

	if( cx1 != cx2 ) {	__DACompareStructures( "px",cx1,cx2 );		*flg = PETSC_FALSE;	}
	if( cy1 != cy2 ) {	__DACompareStructures( "py",cy1,cy2 );		*flg = PETSC_FALSE;	}
	if( cz1 != cz2 ) {	__DACompareStructures( "pz",cz1,cz2 );		*flg = PETSC_FALSE;	}

	if( sw1 != sw2 ) {		__DACompareStructures( "stencil_width",sw1,sw2 );		*flg = PETSC_FALSE;	}
	if( wrap1[0] != wrap2[0] ) {	__DACompareStructures( "wrapX",wrap1[0],wrap2[0] );			*flg = PETSC_FALSE;	}
	if( wrap1[1] != wrap2[1] ) {	__DACompareStructures( "wrapY",wrap1[1],wrap2[1] );			*flg = PETSC_FALSE;	}
	if( wrap1[2] != wrap2[2] ) {	__DACompareStructures( "wrapZ",wrap1[2],wrap2[2] );			*flg = PETSC_FALSE;	}
	if( st1 != st2 ) {		__DACompareStructures( "stencil_type",st1,st2 );		*flg = PETSC_FALSE;	}

	if( si1 != si2 ) {	__DACompareStructures( "si",si1,si2 );		*flg = PETSC_FALSE;	}
	if( sj1 != sj2 ) {	__DACompareStructures( "sj",sj1,sj2 );		*flg = PETSC_FALSE;	}
	if( sk1 != sk2 ) {	__DACompareStructures( "sk",sk1,sk2 );		*flg = PETSC_FALSE;	}

	if( mx1 != mx2 ) {	__DACompareStructures( "mx",mx1,mx2 );		*flg = PETSC_FALSE;	}
	if( my1 != my2 ) {	__DACompareStructures( "my",my1,my2 );		*flg = PETSC_FALSE;	}
	if( mz1 != mz2 ) {	__DACompareStructures( "mz",mz1,mz2 );		*flg = PETSC_FALSE;	}

	if( *flg == PETSC_FALSE ) {
		PetscSynchronizedFlush(comm, PETSC_STDOUT);
	}

	PetscFunctionReturn(0);
}


/*==========================================================================================================*/
/*
    Code for manipulating distributed regular arrays in parallel.

    DAGetProcessorSubset_VerticalDirection - Returns a communicator consisting only of the
    processors in a DMDA that own a particular global x, y grid point
    (corresponding to a line in a 3D grid).

    Collective on DA

    Input Parameters:
 +  da - the distributed array
 .  dir - Cartesian direction, either DA_X, DA_Y, or DA_Z
 -  gp - global grid point number in this direction

    Output Parameters:
 .  comm - new communicator

    Level: advanced

 	Notes:
	This routine is modified from the PETSC routine DAGetProcessorSubset, which
	returns communicators to planes.

	All processors that share the DMDA must call this with the same gp value

   	 This routine is particularly useful to compute boundary conditions
     or other application-specific calculations that require manipulating
     sets of data throughout a logical plane of grid points.

 */
//#include private/daimpl.h

#undef __FUNCT__
#define __FUNCT__ "DAGetProcessorSubset_VerticalDirection"
PetscErrorCode  DAGetProcessorSubset_VerticalDirection(DM da,PetscInt *NumProcs_Z, MPI_Comm *comm)
{
	PetscInt       m,n,p, rank_x, rank_y, rank_z, rank_column;
	PetscMPIInt    size, rank;

	// DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);
	MPI_Comm_size(((PetscObject)da)->comm,&size);
	MPI_Comm_rank(((PetscObject)da)->comm,&rank);

	DMDAGetInfo(da, 0, 0, 0, 0, &m, &n, &p, 0, 0, 0, 0, 0, 0);	// # of processors in every direction

	/* PETSC orders the DA's in a logical manner, which we use here
	 *
	 * First find in which x/y collumn we have to search in the vertical direction
	 * and get the number of the collumn
	 */
	DMDAGetProcessorRank(da, &rank_x, &rank_y, &rank_z, &rank_column);


	/* Split communicator if we have more than 1 PROC in z-direction*/
	 MPI_Comm_split(PETSC_COMM_WORLD,rank_column,rank,comm);

	 *NumProcs_Z = p;


	return(0);
}
/*==========================================================================================================*/


/*==========================================================================================================
 * Locates the processor rank in a 3D DA
 *
 * */
/*
#undef __FUNCT__
#define __FUNCT__ "ProcessorRank_Within_3DDA"
PetscErrorCode ProcessorRank_Within_3DDA(DM da, PetscInt *rank_x, PetscInt *rank_y, PetscInt *rank_z, PetscInt *rank_collumn )
{
	PetscMPIInt	rank;
	PetscInt	m,n,p, ix,iy,iz, irank;

	MPI_Comm_rank(((PetscObject)da)->comm,&rank);

	DMDAGetInfo(da,0,0,0,0,&m,&n,&p,0,0,0,0,0,0);	// # of processors in every direction

	// PETSC orders the DA's in a logical manner, which we use here

	// First find in which x/y collumn we have to search in the vertical direction

	irank = 0;
	for (iz=0; iz<p; iz++){
		for (iy=0; iy<n; iy++){
			for (ix=0; ix<m; ix++){
				if (irank==rank){
					*rank_x = ix;
					*rank_y = iy;
					*rank_z = iz;
				}
				irank=irank+1;
			}
		}

		*rank_collumn = (*rank_x) + (*rank_y)*m;

	}

	PetscFunctionReturn(0);

}
*/

/*==========================================================================================================*/
#undef __FUNCT__
#define __FUNCT__ "LaMEMMod"
PetscErrorCode LaMEMMod( PetscInt x, PetscInt y, PetscInt *result )
{
	PetscFunctionBegin;
	PetscInt     i;

	i   = (x / y);
	x   = x - i * y;
	while (x > y) x -= y;
	*result = x;

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Simple routine to compute result = Mod(x,y) */
PetscInt Mod( PetscInt x,PetscInt y )
{
	PetscInt     i, result;

	i   = x/y; 	x = x- i*y;
	while (x > y) x -= y;
	return result = x;
}
/*==========================================================================================================*/

/* Point in polygon routine, used here to test whether a point is inside a triangle.
 * Modified from PNPOLY by W. Randolph Franklin
 *
 * http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
 *
 * modified to output PetscBool and to use Petsc definitions
 */
PetscBool pnpoly(PetscInt nvert, PetscScalar *vertx, PetscScalar *verty, PetscScalar testx, PetscScalar testy)
{
	PetscInt 		i, j;
	PetscBool		c = PETSC_FALSE;

	for (i = 0, j = nvert-1; i < nvert; j = i++) {
		if ( ((verty[i]>testy) != (verty[j]>testy)) &&
				(testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )

			c = PETSC_NEGATE(c);
	}
	return c;
}



/*==========================================================================================================*/
/* Interpolate a point within a quadratic linear element [used to interpolate the internal free surface]
 * */
PetscScalar InterpolateWithin2DLinearElement(PetscScalar *p_x,PetscScalar *p_y,PetscScalar p_z[4], PetscScalar x, PetscScalar y)
{
	PetscErrorCode ierr;
	PetscScalar	eta, xsi,z_interp, N[4];
	PetscInt 	i;

	z_interp=p_z[0];

	/* Find local coordinates within element */
	ierr = NaturalCoords_Linear2D_Element( &eta, &xsi, x, y, p_x, p_y); CHKERRQ(ierr);

	/* Compute linear shape function */
	N[0] 			=	0.25*(1-eta)*(1-xsi);
	N[1] 			=	0.25*(1+eta)*(1-xsi);
	N[2] 			=	0.25*(1+eta)*(1+xsi);
	N[3] 			=	0.25*(1-eta)*(1+xsi);

	/* Interpolate elevation using shape function */
	z_interp 		=	0;
	for (i=0; i<4; i++){
		z_interp	+= N[i]*p_z[i];
	}

	//PetscPrintf(PETSC_COMM_WORLD,"Reached InterpolateWithin2DLinearElement, with eta=%f, xsi=%f z_interp=%f\n",eta,xsi,z_interp);

	return z_interp;
}
/*==========================================================================================================*/

/*==========================================================================================================*/
//[eta, xsi, N, dNdui_deta, dNdui_dxsi] = NaturalCoordsQuad4(x_real, z_real, eta0, xsi0, ECOORD_x, ECOORD_z)

#undef __FUNCT__
#define __FUNCT__ "NaturalCoords_Linear2D_Element"
PetscErrorCode  NaturalCoords_Linear2D_Element( PetscScalar *eta_output, PetscScalar *xsi_output,PetscScalar x_real,
		PetscScalar y_real, PetscScalar ECOORD_x[4], PetscScalar ECOORD_y[4])
{
	PetscScalar	eta0=0, xsi0=0, Error=100, N[4], dNdui_deta[4], dNdui_dxsi[4], eta, xsi;
	PetscScalar	Jx[2], Jy[2], detJ, invdetJ, invJx[2],invJy[2];
	PetscScalar	Point_realX,Point_realY, diff_x, diff_y, dcoord_x, dcoord_y, factor;
	PetscInt	niter, i;

	eta         =   eta0;
	xsi         =   xsi0;
	niter       =   1;
	while ((Error>5e-6) && (niter<30)) { // can't improve error much more, since single-precision data


		/*
	 	 Compute shape function and derivative @ the current point
		 */
		N[0] 			=	0.25*(1-eta)*(1-xsi);
		N[1] 			=	0.25*(1+eta)*(1-xsi);
		N[2] 			=	0.25*(1+eta)*(1+xsi);
		N[3] 			=	0.25*(1-eta)*(1+xsi);

		dNdui_deta[0] 	=   -0.25+0.25*xsi;
		dNdui_deta[1] 	=    0.25-0.25*xsi;
		dNdui_deta[2] 	=    0.25+0.25*xsi;
		dNdui_deta[3] 	=   -0.25-0.25*xsi;

		dNdui_dxsi[0] 	=   -0.25+0.25*eta;
		dNdui_dxsi[1] 	=   -0.25-0.25*eta;
		dNdui_dxsi[2] 	=    0.25+0.25*eta;
		dNdui_dxsi[3] 	=    0.25-0.25*eta;

		/* Compute jacobian */
		Jx[0]=0; 	Jx[1]=0;
		Jy[0]=0; 	Jy[1]=0;
		for (i=0; i<4; i++){
			Jx[0] += ECOORD_x[i]*dNdui_deta[i];
			Jx[1] += ECOORD_x[i]*dNdui_dxsi[i];
			Jy[0] += ECOORD_y[i]*dNdui_deta[i];
			Jy[1] += ECOORD_y[i]*dNdui_dxsi[i];
		}
		detJ 	= Jx[0]*Jy[1] - Jx[1]*Jy[0];
		invdetJ	=	1.0/detJ;

		/* Inverse jacobian */
		invJx[0]  = +Jy[1]*invdetJ;
		invJx[1]  = -Jy[0]*invdetJ;
		invJy[0]  = -Jx[1]*invdetJ;
		invJy[1]  = +Jx[0]*invdetJ;


		/* Compute real coordinates from natural coordinates */
		Point_realX = 0;	Point_realY = 0;
		for (i=0; i<4; i++){
			Point_realX	+= N[i]*ECOORD_x[i];
			Point_realY += N[i]*ECOORD_y[i];
		}

		/* error */
		diff_x = x_real - Point_realX;
		diff_y = y_real - Point_realY;

		/* compute update */
		dcoord_x	=	diff_x*invJx[0] +   diff_y*invJy[0];
		dcoord_y	=	diff_x*invJx[1] +   diff_y*invJy[1];


		/* add update */
		factor 		=	1.0;
		eta        =    eta + factor*dcoord_x;
		xsi        =    xsi + factor*dcoord_y;

		Error 		=	PetscAbsScalar(diff_x + diff_y);

		//	PetscPrintf(PETSC_COMM_WORLD,"Reached NaturalCoords_Linear2D_Element, iter=%i x,y=[%f,%f] & pointx,y=[%f,%f] with eta=%f, xsi=%f \n",niter, x_real, y_real,Point_realX,Point_realY, eta,xsi);

		niter = niter+1;

	}
	*eta_output = eta;
	*xsi_output = xsi;

	PetscFunctionReturn(0);
}

/*==========================================================================================================*/
/* This is completely stolen from the PETSC routine MatSetValuesStencil 									*/
#undef __FUNCT__
#define __FUNCT__ "StencilToLocalNumbering"
PetscErrorCode  StencilToLocalNumbering(Mat mat,PetscInt m,const MatStencil idxm[], PetscInt jdxm[])
{
	PetscInt        j,i,dim = mat->stencil.dim,*dims = mat->stencil.dims+1,tmp;
	PetscInt       *starts = mat->stencil.starts;
	const PetscInt *dxm = (const PetscInt*)idxm;
	PetscInt 	    sdim = dim - (1 - (PetscInt)mat->stencil.noc);

	if (!m ) return(0); // no values to insert

	if (m > 128) SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_SUP,"Can only set 128 rows at a time; trying to set %lld",(LLD)m);

	for (i=0; i<m; i++) {
		for (j=0; j<3-sdim; j++) dxm++;
		tmp = *dxm++ - starts[0];
		for (j=0; j<dim-1; j++) {
			if ((*dxm++ - starts[j+1]) < 0 || tmp < 0) tmp = PETSC_MIN_INT;
			else                                       tmp = tmp*dims[j] + *(dxm-1) - starts[j+1];
		}
		if (mat->stencil.noc) dxm++;

		jdxm[i] = tmp;
	}

	PetscFunctionReturn(0);
}

/*==========================================================================================================
 * Set initial material properties
 */
#undef __FUNCT__
#define __FUNCT__ "LaMEMInitializeMaterialProperties"
PetscErrorCode  LaMEMInitializeMaterialProperties( UserContext *user )
{
	PetscInt 	i;

	for( i=0; i<user->num_phases; i++ ){

		if (user->DimensionalUnits==1){
			user->PhaseProperties.ViscosityLaw[i]			=	1;		// constant viscosity
			user->PhaseProperties.mu[i]  					=	1e20;
			user->PhaseProperties.n_exponent[i]  			=	1.0;
			user->PhaseProperties.FrankKamenetskii[i]		=	0.0;
			user->PhaseProperties.Powerlaw_e0[i]  			=	1.0;
			user->PhaseProperties.A[i]						=	1e-20;
			user->PhaseProperties.E[i]						=	1.0;

			user->PhaseProperties.DensityLaw[i]				=	1;		// T-dependent density
			user->PhaseProperties.rho[i]  					=	2800;
			user->PhaseProperties.Density_T0[i]				=	273;	// Kelvin (temperature at which rho=rho0, in eq. rho=rho0*(1-alpha*(T-T0)))
			user->PhaseProperties.ThermalExpansivity[i]		=	0;		// no coupling mechanics - thermics


			user->PhaseProperties.PlasticityLaw[i]					=	0;		// none
			user->PhaseProperties.Cohesion[i]  						=	1e100;	// effectively switches off plasticity
			user->PhaseProperties.CohesionAfterWeakening[i]			=	1e100;	// effectively switches off plasticity
			user->PhaseProperties.Weakening_PlasticStrain_Begin[i]	=	0;		//
			user->PhaseProperties.Weakening_PlasticStrain_End[i]	=	0;		//
			user->PhaseProperties.FrictionAngle[i] 					=	0;		// effectively switches off plasticity
			user->PhaseProperties.FrictionAngleAfterWeakening[i] 	=	0;		// effectively switches off plasticity

			user->PhaseProperties.ElasticShearModule[i]  	=	1e100;  // will make the model effectively viscous

			user->PhaseProperties.T_Conductivity[i] 		=	3;		// reasonable value
			user->PhaseProperties.HeatCapacity[i]			=	1050;	// reasonable value
			user->PhaseProperties.RadioactiveHeat[i] 		=	0;
		}
		else{
			user->PhaseProperties.ViscosityLaw[i]			=	1;		// constant viscosity
			user->PhaseProperties.mu[i]  					=	1;
			user->PhaseProperties.n_exponent[i]  			=	1.0;
			user->PhaseProperties.FrankKamenetskii[i]		=	0.0;
			user->PhaseProperties.Powerlaw_e0[i]  			=	1.0;
			user->PhaseProperties.A[i]						=	1.0;
			user->PhaseProperties.E[i]						=	1.0;

			user->PhaseProperties.DensityLaw[i]				=	1;		// T-dependent density
			user->PhaseProperties.rho[i]  					=	1;
			user->PhaseProperties.Density_T0[i]				=	0;	// Kelvin (temperature at which rho=rho0, in eq. rho=rho0*(1-alpha*(T-T0)))
			user->PhaseProperties.ThermalExpansivity[i]		=	0;		// no coupling mechanics - thermics


			user->PhaseProperties.PlasticityLaw[i]					=	0;		// none
			user->PhaseProperties.Cohesion[i]  						=	1e100;	// effectively switches off plasticity
			user->PhaseProperties.CohesionAfterWeakening[i]			=	1e100;	// effectively switches off plasticity
			user->PhaseProperties.Weakening_PlasticStrain_Begin[i]	=	0;		//
			user->PhaseProperties.Weakening_PlasticStrain_End[i]	=	0;		//
			user->PhaseProperties.FrictionAngle[i] 					=	0;		// effectively switches off plasticity
			user->PhaseProperties.FrictionAngleAfterWeakening[i] 	=	0;		// effectively switches off plasticity

			user->PhaseProperties.ElasticShearModule[i]  	=	1e100;  // will make the model effectively viscous

			user->PhaseProperties.T_Conductivity[i] 		=	1;		// reasonable value
			user->PhaseProperties.HeatCapacity[i]			=	1;	// reasonable value
			user->PhaseProperties.RadioactiveHeat[i] 		=	0;

		}


	}

	if (user->DimensionalUnits==1){
		// Add a higher density and viscosity to phase 1
		user->PhaseProperties.rho[1] = 3000;
		user->PhaseProperties.mu[1]  = 1e25;
	}
	else {
		// Add a higher density and viscosity to phase 1
		user->PhaseProperties.rho[1] = 2;
		user->PhaseProperties.mu[1]  = 1e3;
	}


	PetscFunctionReturn(0);
}

/*==========================================================================================================
 * Parse the input file
 */
#undef __FUNCT__
#define __FUNCT__ "LaMEMReadInputFile"
PetscErrorCode   LaMEMReadInputFile( UserContext *user )
{
	FILE *fp;
	PetscInt found;
	double d_values[1000], data;
	PetscInt i_values[1000];
	PetscInt nv, iphase,i;
	const PetscInt max_vals = 1000;
	PetscInt found_data;
	char setup_name[PETSC_MAX_PATH_LEN];


	PetscFunctionBegin;

	fp = fopen( user->ParamFile, "r" );
	if( fp == NULL ) {
		PetscPrintf( PETSC_COMM_WORLD, "LaMEMReadInputFile: Cannot open input file %s \n", user->ParamFile );
		MPI_Abort(PETSC_COMM_WORLD,1);
	}

	parse_GetInt( fp, "nnode_x", &user->nnode_x, &found );
	parse_GetInt( fp, "nnode_y", &user->nnode_y, &found );
	parse_GetInt( fp, "nnode_z", &user->nnode_z, &found );

	/* read the # of elements @ coarsest level ; this override any specification of #nodes */
	parse_GetInt( fp, "nel_x",   &user->nel_x, &found );
	parse_GetInt( fp, "nel_y",   &user->nel_y, &found );
	parse_GetInt( fp, "nel_z",   &user->nel_z, &found );

	/* Characteristic values, used to non-dimensionalize parameters ------------------------------------------------------------- */
	parse_GetInt( fp,    "DimensionalUnits", &user->DimensionalUnits, &found );
	parse_GetDouble( fp, "Characteristic.Length", &data, &found );			if (user->DimensionalUnits==1){		user->Characteristic.Length       = data;	}// read data
	parse_GetDouble( fp, "Characteristic.Viscosity", &data, &found );		if (user->DimensionalUnits==1){		user->Characteristic.Viscosity    = data;	}// read data
	parse_GetDouble( fp, "Characteristic.Temperature", &data, &found );		if (user->DimensionalUnits==1){		user->Characteristic.Temperature  = data; 	}// read data
	parse_GetDouble( fp, "Characteristic.Stress", &data, &found );			if (user->DimensionalUnits==1){		user->Characteristic.Stress  	  = data;	}// read data
	/* ------------------------------------------------------------------------------------------------------------------------- */

	parse_GetDouble( fp, "L", &user->L, &found );
	parse_GetDouble( fp, "W", &user->W, &found );
	parse_GetDouble( fp, "H", &user->H, &found );
	parse_GetDouble( fp, "x_left", &user->x_left, &found );
	parse_GetDouble( fp, "y_front", &user->y_front, &found );
	parse_GetDouble( fp, "z_bot", &user->z_bot, &found );

	// read model setup
	parse_GetString(fp, "msetup", setup_name, PETSC_MAX_PATH_LEN-1, &found);
	if(found)
	{
		if     (!strcmp(setup_name, "parallel"))   user->msetup = PARALLEL;
		else if(!strcmp(setup_name, "redundant"))  user->msetup = REDUNDANT;
		else if(!strcmp(setup_name, "diapir"))     user->msetup = DIAPIR;
		else if(!strcmp(setup_name, "block"))      user->msetup = BLOCK;
		else if(!strcmp(setup_name, "subduction")) user->msetup = SUBDUCTION;
		else if(!strcmp(setup_name, "folding"))    user->msetup = FOLDING;
		else if(!strcmp(setup_name, "detachment")) user->msetup = DETACHMENT;
		else if(!strcmp(setup_name, "slab"))       user->msetup = SLAB;
		else if(!strcmp(setup_name, "spheres"))    user->msetup = SPHERES;
		else if(!strcmp(setup_name, "bands"))      user->msetup = BANDS;
		else SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER,"#ERROR! Incorrect model setup: %s \n", setup_name);
	}

	parse_GetInt( fp,    "Setup.Model", &user->Setup.Model, &found );
	parse_GetDouble( fp, "ampl2D", &user->ampl2D, &found );
	parse_GetDouble( fp, "ampl3D", &user->ampl3D, &found );
	parse_GetDouble( fp, "amplNoise", &user->amplNoise, &found );
	parse_GetDouble( fp, "Hinterface", &user->Hinterface, &found );

	parse_GetString( fp, "OutputFile", user->OutputFile, PETSC_MAX_PATH_LEN-1, &found );
	parse_GetInt( fp,    "save_timesteps", &user->save_timesteps, &found );
	parse_GetInt( fp,    "time_end", &user->time_end, &found );
	parse_GetInt( fp,    "time_end_temp", &user->time_end_temp, &found );
	parse_GetDouble( fp, "CFL", &user->CFL, &found );
	parse_GetDouble( fp, "dt_max", &user->dt_max, &found );
	parse_GetDouble( fp, "dt_temp", &user->dt_temp, &found );

	parse_GetDouble( fp, "BC.Eyy", &user->BC.Eyy, &found );
	parse_GetDouble( fp, "BC.Exx", &user->BC.Exx, &found );

	parse_GetInt( fp, "internalBC_coord", &user->internalBC_coord, &found );
	parse_GetInt( fp, "internalBC_frontel", &user->internalBC_frontel, &found );
	parse_GetInt( fp, "internalBC_backel", &user->internalBC_backel, &found );
	parse_GetInt( fp, "internalBC_node", &user->internalBC_node, &found );
	parse_GetInt( fp, "zdepth_BC_el", &user->zdepth_BC_el, &found );
	parse_GetInt( fp, "zdepth_BC_node", &user->zdepth_BC_node, &found );
	parse_GetInt( fp, "BC.InternalBound", &user->BC.InternalBound, &found );
	parse_GetInt( fp, "BC.LeftBound", &user->BC.LeftBound, &found );
	parse_GetInt( fp, "BC.RightBound", &user->BC.RightBound, &found );
	parse_GetInt( fp, "BC.FrontBound", &user->BC.FrontBound, &found );
	parse_GetInt( fp, "BC.BackBound", &user->BC.BackBound, &found );
	parse_GetInt( fp, "BC.LowerBound", &user->BC.LowerBound, &found );
	parse_GetInt( fp, "BC.UpperBound", &user->BC.UpperBound, &found );
	parse_GetDouble( fp, "Temp_top", &user->Temp_top, &found );
	parse_GetDouble( fp, "Temp_bottom", &user->Temp_bottom, &found );
	parse_GetInt( fp, "temp_initialize", &user->temp_initialize, &found );

	parse_GetInt( fp, "UseInternalFreeSurface", 		&user->ErosionParameters.UseInternalFreeSurface, 		&found );
	parse_GetInt( fp, "StickyAirPhase", 				&user->ErosionParameters.StickyAirPhase, 				&found );
	parse_GetDouble( fp, "FSSA", 	 					&user->FSSA, 											&found );	// FSSA parameter
	parse_GetDouble( fp, "InitialFreeSurfaceHeight", 	&user->ErosionParameters.InitialFreeSurfaceHeight, 		&found );

	// -- parameters related to the erosion model employed
	parse_GetInt( fp, "ErosionModel", 				&user->ErosionParameters.ErosionModel, 		&found );	// which erosion model do we employ?	[0=default=none]

	// parameters in case we use the FD_Erosion model
	parse_GetInt   ( fp, "FE_ErosionCode.ResolutionFactorX", 		&user->ErosionParameters.FE_ErosionCode.ResolutionFactorX, 		&found );
	parse_GetInt   ( fp, "FE_ErosionCode.ResolutionFactorY", 		&user->ErosionParameters.FE_ErosionCode.ResolutionFactorY, 		&found );
	parse_GetDouble( fp, "FE_ErosionCode.InitialRandomNoise_m",  	&user->ErosionParameters.FE_ErosionCode.InitialRandomNoise_m, 	&found );
	parse_GetDouble( fp, "FE_ErosionCode.InitialUpliftedSide_m",	&user->ErosionParameters.FE_ErosionCode.InitialUpliftedSide_m, 	&found );
	parse_GetDouble( fp, "FE_ErosionCode.rain_m_year", 				&user->ErosionParameters.FE_ErosionCode.rain_m_year, 			&found );
	parse_GetDouble( fp, "FE_ErosionCode.k0", 						&user->ErosionParameters.FE_ErosionCode.k0, 					&found ); // k = k0 + c*q^n [q=fluvial discharge,all the others are parameters specified here]
	parse_GetDouble( fp, "FE_ErosionCode.c", 						&user->ErosionParameters.FE_ErosionCode.c, 						&found );
	parse_GetDouble( fp, "FE_ErosionCode.n", 						&user->ErosionParameters.FE_ErosionCode.n, 						&found );
	parse_GetDouble( fp, "FE_ErosionCode.dt", 						&user->ErosionParameters.FE_ErosionCode.dt, 					&found ); // in years
	user->ErosionParameters.FE_ErosionCode.dt = user->ErosionParameters.FE_ErosionCode.dt*3600*24*365.25;	// in seconds
    parse_GetInt   ( fp, "FE_ErosionCode.fill_lake",                &user->ErosionParameters.FE_ErosionCode.fill_lake,              &found );
    parse_GetInt   ( fp, "FE_ErosionCode.BC",                       &user->ErosionParameters.FE_ErosionCode.BC,                     &found );
    parse_GetInt   ( fp, "FE_ErosionCode.mode_river",               &user->ErosionParameters.FE_ErosionCode.mode_river,             &found );
    parse_GetInt   ( fp, "FE_ErosionCode.nbre_river",               &user->ErosionParameters.FE_ErosionCode.nbre_river,             &found );
    parse_GetDouble( fp, "FE_ErosionCode.rain_river_year", 			&user->ErosionParameters.FE_ErosionCode.rain_river_year, 		&found );


	parse_GetInt( fp, "SedimentationModel", 				&user->ErosionParameters.SedimentationModel, 		&found );	// which sedimentation model do we employ?
	parse_GetDouble( fp, "SedimentationRate_cmYr", 	 	&user->ErosionParameters.SedimentationRate_cmYr, 		&found );
	parse_GetInt( fp, 		"PhaseFirstSedimentedLayer", 		&user->ErosionParameters.PhaseFirstSedimentedLayer, 	&found );
	parse_GetInt( fp, 		"PhaseLastSedimentedLayer", 		&user->ErosionParameters.PhaseLastSedimentedLayer, 		&found );
	parse_GetDouble( fp, 	"SedimentLayerThicknessYears", 		&user->ErosionParameters.SedimentLayerThicknessYears, 	&found );

	parse_GetInt( fp,    "GridAdvectionMethod", &user->GridAdvectionMethod, &found );
	parse_GetDouble( fp, "FactorSurfaceLayer", &user->FactorSurfaceLayer, &found );
	parse_GetInt( fp,    "num_subdt", &user->num_subdt, &found );

	// linear solver options
	parse_GetInt( fp,    "StokesSolver", &user->StokesSolver, &found );
	parse_GetInt( fp,    "VelocitySolver", &user->VelocitySolver, &found );

	// nonlinear solver options
	parse_GetDouble( fp, "NonlinearIterationsAccuracy", &user->NonlinearIterationsAccuracy, &found );
	parse_GetInt( fp, "MaxNonlinearIterations", &user->MaxNonlinearIterations, &found );

	/* Particle related variables */
	parse_GetInt( fp,    "ParticleInput", &user->ParticleInput, &found );
	parse_GetInt( fp,    "LoadInitialParticlesFromDisc", &user->LoadInitialParticlesFromDisc, &found );
	parse_GetString( fp, "ParticleFilename", user->ParticleFilename, PETSC_MAX_PATH_LEN-1, &found );
	parse_GetString( fp, "LoadInitialParticlesDirectory", user->LoadInitialParticlesDirectory, PETSC_MAX_PATH_LEN-1, &found );
	if (!found){
		sprintf(user->LoadInitialParticlesDirectory, "InitialParticles");
	}
	parse_GetString( fp, "SaveInitialParticlesDirectory", user->SaveInitialParticlesDirectory, PETSC_MAX_PATH_LEN-1, &found );
	if (!found){
		sprintf(user->SaveInitialParticlesDirectory, "InitialParticles");
	}
	parse_GetInt( fp,    "SaveParticles", &user->SaveParticles, &found );

	parse_GetInt( fp,    "NumPartX", &user->NumPartX, &found );
	parse_GetInt( fp,    "NumPartY", &user->NumPartY, &found );
	parse_GetInt( fp,    "NumPartZ", &user->NumPartZ, &found );

	parse_GetInt( fp,    "InitialMeshFromFile", &user->InitialMeshFromFile, &found );
	parse_GetString( fp, "InitialMeshFileName", user->InitialMeshFileName, PETSC_MAX_PATH_LEN-1, &found );
	parse_GetInt( fp,    "InitialMantleLevel", &user->InitialMantleLevel, &found );

	parse_GetInt( fp,    "InitialErosionSurfaceFromFile", &user->InitialErosionSurfaceFromFile, &found );

	/* Read material properties - This manner of setting material properties will be disabled in the near future, and
	 * replaced with a more general routine
	 */

	// no need for that anymore (will be counted automatically)
//	parse_GetInt( fp, "num_phases", &user->num_phases, &found );

	parse_GetDouble( fp, "LowerViscosityCutoff", &user->LowerViscosityCutoff, &found );
	parse_GetDouble( fp, "UpperViscosityCutoff", &user->UpperViscosityCutoff, &found );

	parse_GetDouble( fp, "Gravity", &user->Gravity, &found );
	parse_GetInt( fp,    "PlasticityModel", &user->PlasticityModel, &found );
	parse_GetDouble( fp, "Xi", &user->Xi, &found );
	parse_GetDouble( fp, "GasConstant", &user->GasConstant, &found );

	// --- SurfaceVelocity related input parameters ---
	parse_GetInt( fp,    "get_SurfVelField", 	&user->SurfVelField.GetIt, 		&found );
	parse_GetInt( fp,    "SurfVelField_SaveRef", &user->SurfVelField.SaveRef,	&found );
	parse_GetString( fp, "SurfVelField_RefDatFile", user->SurfVelField.RefDatFile2load, PETSC_MAX_PATH_LEN-1, &found );
	parse_GetDouble( fp, "SurfVelField_VxStdDev", &user->SurfVelField.VxStdDev, &found );
	parse_GetDouble( fp, "SurfVelField_VyStdDev", &user->SurfVelField.VyStdDev, &found );
	parse_GetDouble( fp, "SurfVelField_VzStdDev", &user->SurfVelField.VzStdDev, &found );
	// --- Optimization related input parameters ---
	parse_GetInt( fp,    "get_Misfit", 				&user->Optimisation.GetIt, 		&found );

	// --- Gravity related input parameters ---
	parse_GetInt( fp,    "get_GravityField", 		&user->GravityField.GetIt, 		&found );
	parse_GetInt( fp,    "GravityField_SaveRef", 	&user->GravityField.SaveRef, 	&found );
	parse_GetInt( fp,    "GravityField_SaveDebug", 	&user->GravityField.SaveDebug, 	&found );
	parse_GetInt( fp,    "GravityField_SaveVTK", 	&user->GravityField.SaveVTK, 	&found );
	parse_GetInt( fp,    "GravityField_survey_nx", &user->GravityField.survey_nx	, &found );
	parse_GetInt( fp,    "GravityField_survey_ny", &user->GravityField.survey_ny	, &found );
	parse_GetDouble( fp, "GravityField_survey_xs", &user->GravityField.survey_xs	, &found );
	parse_GetDouble( fp, "GravityField_survey_xm", &user->GravityField.survey_xm	, &found );
	parse_GetDouble( fp, "GravityField_survey_ys", &user->GravityField.survey_ys	, &found );
	parse_GetDouble( fp, "GravityField_survey_ym", &user->GravityField.survey_ym	, &found );
	parse_GetDouble( fp, "GravityField_survey_z"	, &user->GravityField.survey_z 	, &found );
	parse_GetDouble( fp, "GravityField_ReferenceDensity", &user->GravityField.ReferenceDensity, &found );
	parse_GetInt( fp,    "GravityField_num_intp"	, &user->GravityField.num_intp		, &found );
	parse_GetString( fp, "GravityField_RefDatFile", user->GravityField.RefDatFile2load, PETSC_MAX_PATH_LEN-1, &found );
	parse_GetDouble( fp, "GravityField_StdDev", &user->GravityField.StdDev, &found );
	parse_GetInt( fp,    "GravityField_LithColNum", &user->GravityField.LithColNum, &found );
	parse_GetDoubleArray(fp,"GravityField_LithColDepth",&nv,d_values, &found );
    if (found!=0){
        for( i=0; i<user->GravityField.LithColNum; i++ )
        {
            user->GravityField.LithColDepth[i] = d_values[i] ;
            PetscPrintf(PETSC_COMM_WORLD,"# LithColDepth[%lld] = %g \n",(LLD) i,user->GravityField.LithColDepth[i]);
		}
	}

	parse_GetDoubleArray( fp,    "GravityField_LithColDens",&nv,d_values, &found );
    if (found!=0){
        for( i=0; i<user->GravityField.LithColNum+1; i++ )
        {
            user->GravityField.LithColDens[i] = d_values[i] ;
            PetscPrintf(PETSC_COMM_WORLD,"# LithColDens[%lld] = %g \n",(LLD) i,user->GravityField.LithColDens[i]);
		}
	}

	// --- Isostasy related input parameters ---
	parse_GetInt( fp,    "get_AiryIsostasy", &user->Isostasy.GetIt, &found );
	parse_GetInt( fp,    "Isostasy_SaveRef", &user->Isostasy.SaveRef, &found );
	parse_GetInt( fp, "Isostasy_ref_xi", &user->Isostasy.ref_xi, &found );
	parse_GetInt( fp, "Isostasy_ref_yi", &user->Isostasy.ref_yi, &found );
	parse_GetDouble( fp, "Isostasy_corr_topo", &user->Isostasy.corr_topo,  &found );
	parse_GetDouble( fp, "Isostasy_TisoStdDev", &user->GravityField.StdDev, &found );
	parse_GetDouble( fp, "Isostasy_RefRho", &user->Isostasy.ref_rho, &found );
	parse_GetString( fp, "Isostasy_RefDatFile", user->Isostasy.RefDatFile2load, PETSC_MAX_PATH_LEN-1, &found );


	parse_GetInt( fp,    "save_breakpoints", &user->save_breakpoints, &found );

	// --- Output related input parameters ---
	parse_GetInt( fp,    "Output.velocity", &user->Output.velocity, &found );
	parse_GetInt( fp,    "Output.temperature", &user->Output.temperature, &found );
	parse_GetInt( fp,    "Output.surface_topography", &user->Output.surface_topography, &found );
	parse_GetInt( fp,    "Output.bottom_topography", &user->Output.bottom_topography, &found );
	parse_GetInt( fp,    "Output.quadrature", &user->Output.quadrature, &found );

	// --- Pushing related input parameters ---
	parse_GetInt( fp,    "AddPushing", &user->AddPushing, &found );
	parse_GetInt( fp,    "Pushing.num_changes", &user->Pushing.num_changes, &found );
	parse_GetInt( fp,    "Pushing.reset_pushing_coord", &user->Pushing.reset_pushing_coord, &found );
	parse_GetDoubleArray( fp,    "Pushing.time",  &nv, d_values, &found );
	for( i=0; i<user->Pushing.num_changes+1; i++ ) {
		user->Pushing.time[i] = d_values[i] ;
		//PetscPrintf(PETSC_COMM_WORLD,"# time stage = %g \n",user->Pushing.time[i]);
	}
	parse_GetDoubleArray( fp, "Pushing.V_push", &nv, d_values, &found );
	for( i=0; i<user->Pushing.num_changes; i++ ) {
		user->Pushing.V_push[i] = d_values[i] ;
		//PetscPrintf(PETSC_COMM_WORLD,"# V_push stage = %g \n",user->Pushing.V_push[i]);
	}
	parse_GetDoubleArray( fp, "Pushing.omega",&nv, d_values, &found  );
	for( i=0; i<user->Pushing.num_changes; i++ ) {
		user->Pushing.omega[i] = d_values[i] ;
		//PetscPrintf(PETSC_COMM_WORLD,"# omega = %g \n",user->Pushing.omega[i]);
	}
	parse_GetIntArray( fp, "Pushing.coord_advect",&nv, i_values, &found);
	for( i=0; i<user->Pushing.num_changes; i++ ) {
		user->Pushing.coord_advect[i] = i_values[i] ;
		//PetscPrintf(PETSC_COMM_WORLD,"# coord_advect = %d \n",user->Pushing.coord_advect[i]);
	}
	parse_GetIntArray( fp, "Pushing.dir",&nv, i_values, &found);
	for( i=0; i<user->Pushing.num_changes; i++ ) {
		user->Pushing.dir[i] = i_values[i] ;
		//PetscPrintf(PETSC_COMM_WORLD,"# direction = %d \n",user->Pushing.dir[i]);
	}
	parse_GetDouble( fp, "Pushing.L_block", &user->Pushing.L_block, &found );
	parse_GetDouble( fp, "Pushing.W_block", &user->Pushing.W_block, &found );
	parse_GetDouble( fp, "Pushing.H_block", &user->Pushing.H_block, &found );
	parse_GetDouble( fp, "Pushing.x_center_block", &user->Pushing.x_center_block, &found );
	parse_GetDouble( fp, "Pushing.y_center_block", &user->Pushing.y_center_block, &found );
	parse_GetDouble( fp, "Pushing.z_center_block", &user->Pushing.z_center_block, &found );

	/* -------------------------------------------------------------------------------------------------------------------------
	 * Read phase transitions related information -
	 * The input structure of this is also likely to change but at a later stage
	 */
	parse_GetInt( fp,    "num_phase_transitions", &user->num_phase_transitions, &found );

	parse_GetIntAllInstances( fp,    "TransitionType", &nv, i_values, max_vals, &found );	//	if( nv != user->num_phase_transitions ) {  PetscPrintf( PETSC_COMM_WORLD, "# Warning: num_phase_transitions is inconsistent\n");  }
	for( iphase=0; iphase<user->num_phase_transitions; iphase++ ) {
		user->PhaseTransitions[iphase].TransitionType = i_values[iphase] ;
	}

	parse_GetIntAllInstances( fp,    "TransitionBelow", &nv, i_values, max_vals, &found );	//	if( nv != user->num_phase_transitions ) {  PetscPrintf( PETSC_COMM_WORLD, "# Warning: num_phase_transitions is inconsistent\n");  }
	for( iphase=0; iphase<user->num_phase_transitions; iphase++ ) {
		user->PhaseTransitions[iphase].TransitionBelow = i_values[iphase] ;
	}

	parse_GetDoubleAllInstances( fp, "TransitionDepth", &nv, d_values, max_vals, &found );	//	if( nv != user->num_phase_transitions ) {  PetscPrintf( PETSC_COMM_WORLD, "# Warning: num_phase_transitions is inconsistent\n");  }
	for( iphase=0; iphase<user->num_phase_transitions; iphase++ ) {
		user->PhaseTransitions[iphase].TransitionDepth = d_values[iphase] ;
	}

	parse_GetDoubleAllInstances( fp, "TransitionP0", &nv, d_values, max_vals, &found );		//	if( nv != user->num_phase_transitions ) {  PetscPrintf( PETSC_COMM_WORLD, "# Warning: num_phase_transitions is inconsistent\n");  }
	for( iphase=0; iphase<user->num_phase_transitions; iphase++ ) {
		user->PhaseTransitions[iphase].TransitionP0 = d_values[iphase] ;
	}

	parse_GetDoubleAllInstances( fp, "TransitionAlpha", &nv, d_values, max_vals, &found );	//	if( nv != user->num_phase_transitions ) {  PetscPrintf( PETSC_COMM_WORLD, "# Warning: num_phase_transitions is inconsistent\n");  }
	for( iphase=0; iphase<user->num_phase_transitions; iphase++ ) {
		user->PhaseTransitions[iphase].TransitionAlpha = d_values[iphase] ;
	}

	parse_GetIntAllInstances( fp,    "InitialPhase", &nv, i_values, max_vals, &found );		//	if( nv != user->num_phase_transitions ) {  PetscPrintf( PETSC_COMM_WORLD, "# Warning: num_phase_transitions is inconsistent\n");  }
	for( iphase=0; iphase<user->num_phase_transitions; iphase++ ) {
		user->PhaseTransitions[iphase].InitialPhase = i_values[iphase] ;
	}

	parse_GetIntAllInstances( fp,    "TransformedPhase", &nv, i_values, max_vals, &found );	//	if( nv != user->num_phase_transitions ) {  PetscPrintf( PETSC_COMM_WORLD, "# Warning: num_phase_transitions is inconsistent\n");  }
	for( iphase=0; iphase<user->num_phase_transitions; iphase++ ) {
		user->PhaseTransitions[iphase].TransformedPhase = i_values[iphase] ;
	}
	/* ------------------------------------------------------------------------------------------------------------------------- */

	fclose( fp );

	/* -------------------------------------------------------------------------------------------------------------------------
	 * Read material parameters from the input file (new, more general, input format for material properties)
	 *
	 * At this stage, the new data is read from the input file and added to the phase properties array.
	 *
	 */

	MaterialCreate( &user->PhaseMaterialProperties );
	MaterialReadFromFile( user->PhaseMaterialProperties, 	user->ParamFile, &ConstantViscosityReadFromFile,			&found_data	);   	/* Constant viscosity params */
	MaterialReadFromFile( user->PhaseMaterialProperties, 	user->ParamFile, &PowerLawViscosityReadFromFile,			&found_data	);   	/* Power-law viscosity params */
	MaterialReadFromFile( user->PhaseMaterialProperties, 	user->ParamFile, &TempDepViscosityReadFromFile,				&found_data	);   	/* Temp-dep viscosity params */
	MaterialReadFromFile( user->PhaseMaterialProperties, 	user->ParamFile, &TempDepNoPowerLawViscosityReadFromFile,	&found_data	);   	/* Temp-dep viscosity params without power law */
	MaterialReadFromFile( user->PhaseMaterialProperties, 	user->ParamFile, &ConstantElasticityReadFromFile,			&found_data	);  	/* Elasticity params */
	MaterialReadFromFile( user->PhaseMaterialProperties, 	user->ParamFile, &TemperatureDependentDensityReadFromFile, 	&found_data );  	/* Density params */
	MaterialReadFromFile( user->PhaseMaterialProperties, 	user->ParamFile, &ConvectionDensityReadFromFile, 			&found_data );  	/* Density params */
	MaterialReadFromFile( user->PhaseMaterialProperties, 	user->ParamFile, &ArtificialTemperatureDependentDensityReadFromFile, 	&found_data );  	/* Density params */

	MaterialReadFromFile( user->PhaseMaterialProperties, 	user->ParamFile, &ConstantEnergyReadFromFile, 				&found_data );  	/* Energy equation params */
	MaterialReadFromFile( user->PhaseMaterialProperties, 	user->ParamFile, &ConstantPlasticityReadFromFile, 			&found_data );  	/* Plasticity params */

// DUPLICATION; ALREADY CALLED THIS ROUTINE!!!
//	PetscOptionsReadFromFile(user->ParamFile, &found_data, 1);

	/* Print an overview of Material parameters read from file */
	PetscPrintf(PETSC_COMM_WORLD,"# --------------------------------------------------------------------------\n");
	PetscPrintf(PETSC_COMM_WORLD,"# Phase material parameters read from %s: \n",user->ParamFile);

	{
		PetscInt      *check_phase;
		PetscBool      empty_phase;
		PetscErrorCode ierr;

		PetscInt n_phases,n_attrs, n_types;
		Phase *phases, p;
		Attribute *attrs, a;
		AttributeType *types, t;
		PetscInt PP,AA,TT;
		PetscInt attr_match, type_match;
		PetscInt p_idx, a_idx, t_idx;

		MaterialGetAllPhases( user->PhaseMaterialProperties, &n_phases, &phases );

		// check total number of phases found
		if(n_phases > max_num_phases)
		{
			SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many phases in the input file! Actual: %lld, Maximum: %lld\n", (LLD)n_phases, (LLD)max_num_phases);
		}

		// initialize phase checking array
		ierr = makeIntArray(&check_phase, NULL, n_phases); CHKERRQ(ierr);

		for(PP = 0; PP < n_phases; PP++) check_phase[PP] = 0;

		// store total number of phases
		user->num_phases = n_phases;

		for( PP=0; PP<n_phases; PP++ )
		{
			p = phases[PP];

			// check phase numbering
			if(p->phase_number > n_phases-1)
			{
				SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_USER, "Phase numbering out of bound! Phase ID: %lld, Maximum ID: %lld", (LLD)p->phase_number, (LLD)(n_phases-1));
			}
			check_phase[p->phase_number] = 1;

			p_idx = p->P_id;

			PhaseGetAllAttributes( p, &n_attrs, &attrs );

			for( AA=0; AA<n_attrs; AA++ ) {
				a = attrs[AA];
				a_idx = a->A_id;
				AttributeGetAllTypes( a, &n_types, &types );

				for( TT=0; TT<n_types; TT++ ) {
					t = types[TT];
					t_idx = t->T_id;

					AttributeCompare( a, "VISCOSITY", &attr_match );
					if( attr_match == _TRUE ) {
						// constant
						AttributeTypeCompare( t, "constant", &type_match );
						if( type_match == _TRUE ) {
							double eta0;

							MaterialGetConstantViscosityParams( user->PhaseMaterialProperties, p_idx,a_idx,t_idx, &eta0 );
							PetscPrintf(PETSC_COMM_WORLD,"#    %s (id=%lld), VISCOSITY, constant [%lld,%lld,%lld] = %e \n", p->name, (LLD)(p->phase_number), (LLD)p_idx,(LLD)a_idx,(LLD)t_idx, eta0 );

							// Add to PhaseProperties
							user->PhaseProperties.ViscosityLaw[p->phase_number] = 1; 		// constant viscosity
							user->PhaseProperties.mu[p->phase_number] 			= eta0;



						}
						AttributeTypeCompare( t, "power_law", &type_match );
						if( type_match == _TRUE ) {
							double eta0;
							double N_exp;
							double e0;

							MaterialGetPowerLawViscosityParams( user->PhaseMaterialProperties, p_idx,a_idx,t_idx, &eta0, &N_exp, &e0 );
							PetscPrintf(PETSC_COMM_WORLD,"#    %s (id=%lld), VISCOSITY, power-law [%lld,%lld,%lld] = %e, %e %e \n", p->name, (LLD)(p->phase_number), (LLD)p_idx,(LLD)a_idx,(LLD)t_idx, eta0, N_exp, e0 );


							// Add to PhaseProperties
							user->PhaseProperties.ViscosityLaw[p->phase_number] = 2; 		// powerlaw viscosity, given by eta=eta0*(e2nd/e0)^(1/n-1)
							user->PhaseProperties.mu[p->phase_number] 			= eta0;
							user->PhaseProperties.n_exponent[p->phase_number] 	= N_exp;
							user->PhaseProperties.Powerlaw_e0[p->phase_number] 	= e0;		//	e0 in (e2nd/e0)



						}
						AttributeTypeCompare( t, "temp_dep", &type_match );
						if( type_match == _TRUE ) {
							double PreExpFactor;
							double N_exp;
							double e0;
							double ActivationEnergy;

							MaterialGetTempDepViscosityParams( user->PhaseMaterialProperties, p_idx,a_idx,t_idx, &PreExpFactor, &N_exp, &e0, &ActivationEnergy );
							PetscPrintf(PETSC_COMM_WORLD,"#    %s (id=%lld), VISCOSITY, temp-dep [%lld,%lld,%lld] = %e, %e %e %e \n", p->name, (LLD)(p->phase_number), (LLD)p_idx,(LLD)a_idx,(LLD)t_idx, PreExpFactor, N_exp, e0, ActivationEnergy );


							// Add to PhaseProperties
							user->PhaseProperties.ViscosityLaw[p->phase_number] = 4;
							user->PhaseProperties.A[p->phase_number] 			= PreExpFactor;
							user->PhaseProperties.n_exponent[p->phase_number] 	= N_exp;
							user->PhaseProperties.Powerlaw_e0[p->phase_number] 	= e0;		//	e0 in (e2nd/e0)
							user->PhaseProperties.E[p->phase_number]			= ActivationEnergy;


						}
						AttributeTypeCompare( t, "tempdep_nopowerlaw", &type_match );
						if( type_match == _TRUE ) {
							double PreExpFactor;
							double N_exp;
							double e0;
							double ActivationEnergy;

							MaterialGetTempDepNoPowerLawViscosityParams( user->PhaseMaterialProperties, p_idx,a_idx,t_idx, &PreExpFactor, &N_exp, &e0, &ActivationEnergy );
							PetscPrintf(PETSC_COMM_WORLD,"#    %s (id=%lld), VISCOSITY, tempdep_nopowerlaw [%lld,%lld,%lld] = %e, %e %e %e \n", p->name, (LLD)(p->phase_number), (LLD)p_idx,(LLD)a_idx,(LLD)t_idx, PreExpFactor, N_exp, e0, ActivationEnergy );


							// Add to PhaseProperties
							user->PhaseProperties.ViscosityLaw[p->phase_number] = 5;
							user->PhaseProperties.A[p->phase_number] 			= PreExpFactor;
							user->PhaseProperties.n_exponent[p->phase_number] 	= N_exp;
							user->PhaseProperties.Powerlaw_e0[p->phase_number] 	= e0;		//	e0 in (e2nd/e0)
							user->PhaseProperties.E[p->phase_number]			= ActivationEnergy;


						}
						// others
					}
					AttributeCompare( a, "ELASTICITY", &attr_match );
					if( attr_match == _TRUE ) {
						// constant
						AttributeTypeCompare( t, "constant", &type_match );
						if( type_match == _TRUE ) {
							double shear, bulk;

							MaterialGetConstantElasticityParams( user->PhaseMaterialProperties, p_idx,a_idx,t_idx, &shear, &bulk );
							PetscPrintf(PETSC_COMM_WORLD,"#    %s (id=%lld), ELASTICITY, constant [%lld,%lld,%lld] = %e, %e \n", p->name, (LLD)(p->phase_number), (LLD)p_idx,(LLD)a_idx,(LLD)t_idx, shear,bulk );

							user->PhaseProperties.ElasticShearModule[p->phase_number] 		= shear; 	// add to PhaseProperties

						}
						// others

					}
					AttributeCompare( a, "DENSITY", &attr_match );
					if( attr_match == _TRUE ) {
						// temp dep.
						AttributeTypeCompare( t, "temperature_dependent", &type_match );
						if( type_match == _TRUE ) {
							double rho0, alpha, T0;

							MaterialGetTemperatureDependentDensityParams( user->PhaseMaterialProperties, p_idx,a_idx,t_idx, &rho0, &alpha, &T0 );
							PetscPrintf(PETSC_COMM_WORLD,"#    %s (id=%lld), DENSITY, temperature-dependent [%lld,%lld,%lld] = %e, %e, %e \n", p->name, (LLD)(p->phase_number), (LLD)p_idx,(LLD)a_idx,(LLD)t_idx, rho0,alpha, T0 );


							user->PhaseProperties.DensityLaw[p->phase_number] 	= 1; 		// T-dependent density
							user->PhaseProperties.rho[p->phase_number] 			= rho0; 	// add to PhaseProperties
							user->PhaseProperties.ThermalExpansivity[p->phase_number] 		= alpha; 	// add to PhaseProperties
							user->PhaseProperties.Density_T0[p->phase_number] 	= T0; 		// add to PhaseProperties

						}
						// convection setup
						AttributeTypeCompare( t, "convection", &type_match );
						if( type_match == _TRUE ) {
							double rho0, Ra, T0;

							MaterialGetConvectionDensityParams( user->PhaseMaterialProperties, p_idx,a_idx,t_idx, &rho0, &Ra, &T0 );
							PetscPrintf(PETSC_COMM_WORLD,"#    %s (id=%lld), DENSITY, convection [%lld,%lld,%lld] = %e, %e, %e \n", p->name, (LLD)(p->phase_number), (LLD)p_idx,(LLD)a_idx,(LLD)t_idx, rho0,Ra, T0 );


							user->PhaseProperties.DensityLaw[p->phase_number] 	= 2; 				// Density used in Ra-convection simulations where the rhs of the force balance 2 is Ra*T*g instead of rho*g
							user->PhaseProperties.rho[p->phase_number] 			= rho0; 			// add to PhaseProperties
							user->PhaseProperties.Ra[p->phase_number] 			= Ra; 				// add to PhaseProperties
							user->PhaseProperties.Density_T0[p->phase_number] 	= T0; 				// add to PhaseProperties
						}
						// artificial temp dep.
						AttributeTypeCompare( t, "artificial_temperature_dependent", &type_match );
						if( type_match == _TRUE ) {
							double rho0;

							MaterialGetArtificialTemperatureDependentDensityParams( user->PhaseMaterialProperties, p_idx,a_idx,t_idx, &rho0);
							PetscPrintf(PETSC_COMM_WORLD,"#    %s (id=%lld), DENSITY, artificial-temperature-dependent [%lld,%lld,%lld] = %e \n", p->name, (LLD)(p->phase_number), (LLD)p_idx,(LLD)a_idx,(LLD)t_idx, rho0);

							user->PhaseProperties.DensityLaw[p->phase_number] 	= 3; 		// artificial T-dependent density
							user->PhaseProperties.rho[p->phase_number] 			= rho0; 	// add to PhaseProperties
						}

						// others
					}
					AttributeCompare( a, "ENERGY", &attr_match );
					if( attr_match == _TRUE ) {
						// constant
						AttributeTypeCompare( t, "constant", &type_match );
						if( type_match == _TRUE ) {
							double  ThermalConductivity, HeatCapacity, RadioactiveHeat, ShearHeating;

							MaterialGetConstantEnergyParams( user->PhaseMaterialProperties, p_idx,a_idx,t_idx, &ThermalConductivity, &HeatCapacity, &RadioactiveHeat, &ShearHeating);
							PetscPrintf(PETSC_COMM_WORLD,"#    %s (id=%lld), ENERGY, constant [%lld,%lld,%lld] = %e %e %e %e \n", p->name, (LLD)(p->phase_number), (LLD)p_idx,(LLD)a_idx,(LLD)t_idx, ThermalConductivity, HeatCapacity, RadioactiveHeat, ShearHeating);


							user->PhaseProperties.T_Conductivity[p->phase_number] 		= ThermalConductivity; 	// add to PhaseProperties
							user->PhaseProperties.HeatCapacity[p->phase_number] 		= HeatCapacity; 		// add to PhaseProperties
							user->PhaseProperties.RadioactiveHeat[p->phase_number] 		= RadioactiveHeat; 		// add to PhaseProperties


						}
						// others
					}
					AttributeCompare( a, "PLASTICITY", &attr_match );
					if( attr_match == _TRUE ) {
						// constant
						AttributeTypeCompare( t, "DruckerPrager", &type_match );
						if( type_match == _TRUE ) {
							double  Cohesion, FrictionAngle;
							double Weakening_PlasticStrain_Begin,Weakening_PlasticStrain_End;
							double  CohesionAfterWeakening, FrictionAngleAfterWeakening;

							MaterialGetConstantPlasticityParams( user->PhaseMaterialProperties, p_idx,a_idx,t_idx, &Cohesion, &FrictionAngle,
									&Weakening_PlasticStrain_Begin, &Weakening_PlasticStrain_End,
									&CohesionAfterWeakening, &FrictionAngleAfterWeakening);
							PetscPrintf(PETSC_COMM_WORLD,"#    %s (id=%lld), PLASTICITY, constant [%lld,%lld,%lld] = %e %e %e %e %e %e\n", p->name, (LLD)(p->phase_number), (LLD)p_idx,(LLD)a_idx,(LLD)t_idx,
									Cohesion, FrictionAngle, Weakening_PlasticStrain_Begin, Weakening_PlasticStrain_End, CohesionAfterWeakening, FrictionAngleAfterWeakening);

							user->PhaseProperties.PlasticityLaw[p->phase_number] 	= 1; 					// Drucker-Prager
							user->PhaseProperties.Cohesion[p->phase_number] 		= Cohesion; 			// add to PhaseProperties
							user->PhaseProperties.FrictionAngle[p->phase_number] 	= FrictionAngle; 		// add to PhaseProperties

							user->PhaseProperties.Weakening_PlasticStrain_Begin[p->phase_number] = Weakening_PlasticStrain_Begin; 			// add to PhaseProperties
							user->PhaseProperties.Weakening_PlasticStrain_End[p->phase_number] 	 = Weakening_PlasticStrain_End; 		  	// add to PhaseProperties
							user->PhaseProperties.CohesionAfterWeakening[p->phase_number] 	 	 = CohesionAfterWeakening; 		  						// add to PhaseProperties
							user->PhaseProperties.FrictionAngleAfterWeakening[p->phase_number] 	 = FrictionAngleAfterWeakening; 		  	// add to PhaseProperties

						}
						// others
					}
				}
			}
		}

		// check empty phases
		empty_phase = PETSC_FALSE;

		for(PP = 0; PP < n_phases; PP++)
		{
			if(!check_phase[PP])
			{
				PetscPrintf(PETSC_COMM_WORLD, "Phase %lld is not initialized\n", (LLD)PP);

				empty_phase = PETSC_TRUE;
			}
		}

		ierr = PetscFree(check_phase); CHKERRQ(ierr);

		if(empty_phase == PETSC_TRUE)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Incorrect phase numbering");
		}

	}
	/* ------------------------------------------------------------------------------------------------------------------------- */
	// destroy attribute database right after use (it's not used anywhere else anymore, because data is copied to internal structures)
	MaterialDestroy(user->PhaseMaterialProperties);

	PetscFunctionReturn(0);
}

/*==========================================================================================================*/
/* Set default code parameters and read input file, if required */
#undef __FUNCT__
#define __FUNCT__ "InitializeCode"
PetscErrorCode InitializeCode( UserContext *user )
{
	PetscMPIInt    size;
	PetscErrorCode ierr;
	PetscInt       nx,ny,nz, mod, i, nel_array[3], nel_input_max;
    PetscInt       n_int;
	PetscScalar    dx, dy, dz, SecYear;
	PetscBool      found,flg;
	char          *all_options;
	char           setup_name[PETSC_MAX_PATH_LEN];

	// NEW OPTIONS WILL BE ADDED *** BEFORE *** ALREADY SPECIFIED  
	PetscOptionsGetAll( &all_options ); /* copy all command line args */

	SecYear 			=	3600*24*365.25;		// seconds per year

	/* Set default values */
	user->W 		  	= 		1.0;
	user->L 			= 		1.0;
	user->H 		 	= 		1.0;
	user->y_front 		= 		0.0;
	user->ampl2D 	  	= 		0.0;
	user->ampl3D  		= 		1e-2;
	user->amplNoise 	= 		0.0;
	user->Hinterface 	=		0.5;			// average interface height for diapir setup [0-1]
	user->mumax			= 		1.0;

	// set default grid sizes
	user->nel_x  		=       8;
	user->nel_y  		=       8;
	user->nel_z			=       8;
	user->refinex     	=     	2;
	user->refiney   	=     	2;
	user->refinez  		=  		2;

	user->time_start 	=  		0.0;
	user->time_end 		= 		1.0;
	user->time_end_temp	= 		1;
	user->save_timesteps=		1;
	user->time		  	=   	0.0;
	user->CFL 		 	= 		0.5;
	user->Temp_top  		= 	0.0;
	user->Temp_bottom 		=	1.0;
	user->temp_initialize	=	0.0;
	user->DimensionalUnits 	= 	0;
	user->Setup.Model 		=	2; // 0-Diapir, 1-Single Layer Fold, 2-Falling block, 3-Particle file (defined below), 4-multilayer detachment fold, 6-subduction w/sticky air
	user->Gravity   		= 	1.0;
	user->GasConstant  		= 	1.0;
	user->BC.Vy_front		=	0;
	user->BC.Vy_back		=	0;
	user->BC.Vz_bot 		=	0;
	user->BC.Vz_top   		=	0;
	user->BC.Vx_left		=	0;
	user->BC.Vx_right		=	0;
	user->BC.Exx			=	0;
	user->BC.Eyy		 	=	0;
	sprintf(user->OutputFile, 		"FallingBlock3D");
	sprintf(user->ParticleFilename, "ParticlesInput.dat");
	sprintf(user->LoadInitialParticlesDirectory, "InitialParticles");

	// FDSTAG Canonical Model Setup
	user->msetup = BLOCK;

	user->MatlabOutputFiles	=	1;		// write MATLAB output
	user->VTKOutputFiles	=	1;		// write VTK output
	user->AVDPhaseViewer  	= 	0;
	user->SavePartitioning  = 	PETSC_FALSE;		// write partitioning
	user->x_left 	  		= 	-0.5*(user->W)*0.0;
	user->y_front	  		= 	-0.5*(user->L)*0.0;
	user->z_bot 			= 	-0.5*(user->H)*0.0;

	// linear solver settings
	user->SkipStokesSolver     = PETSC_FALSE;
	user->StokesSolver         = 1; // 1 - Powell-Hesteness iterations; 2 - Schur Complement Reduction; 3 - Fully Coupled Solver; 4 - MatVec Test;
	user->VelocitySolver       = 1; // 0 - User defined; 1 - Direct (MUMPS); 2 - Galerkin geometric multigrid; 3 - Fieldsplit + Algebraic Multigrid (ML), 4 - GCR & HYPRE preconditioner on full velocity block
	user->VelocityTest         = PETSC_FALSE; // Request to perform single velocity solve for test purposes & quit
	user->ScaleSystem          = PETSC_FALSE; // Request to scale linear system before solution
	user->use_fdstag_canonical = PETSC_FALSE; // request native staggered grid discretization

	// nonlinear solver settings
	user->MaxNonlinearIterations      =	25;
	user->NonlinearIterationsAccuracy =	1e-3;

	user->internalBC_coord			= 0.0;
	user->internalBC_frontel		= 0.0;
	user->internalBC_backel			= 0.0;
	user->internalBC_node		    = 0.0;
	user->zdepth_BC_el			    = 0.0;
	user->zdepth_BC_node		    = 0.0;
	user->BC.InternalBound			=   0;	// 0-free surface, 1-free slip 					, 2-no-slip
	user->BC.UpperBound				=   1;	// 0-free surface, 1-free slip 					, 2-no-slip
	user->BC.LowerBound				=	1;	// 0-free surface, 1-free slip					, 2-no-slip
	user->BC.LeftBound				=	1;	// 0-free surface, 1-free slip w. BG strainrate	, 2-no slip
	user->BC.RightBound				=	1;	// 0-free surface, 1-free slip w. BG strainrate	, 2-no-slip
	user->BC.FrontBound				=	1;	// 0-free surface, 1-free slip w. BG strainrate	, 2-no-slip
	user->BC.BackBound				=	1;	// 0-free surface, 1-free slip w. BG strainrate	, 2-no-slip
	user->NumParticlesToStartInjection = 0;  	// how many particles should an element have before we start injecting particles?
	user->ParticleInjectionPhase       = 0;  	// which is the phase of the injected particle?
	user->ParticleInput = 1;					// 0-do not use particles to track phases; 1-do use particles to track phases
	user->LoadInitialParticlesFromDisc = 0;		// Read the initial particles from disc
	user->remesh = 0;
	user->CriticalDiagonalRatio = 0.55;			// save range is [0.4-1.0]

	// Check if we are performing a benchmark with a build-in analytical benchmark
	ierr 	 = PetscOptionsGetBool( PETSC_NULL, "-Benchmark_SolCx",   		  		&user->AnalyticalBenchmark, PETSC_NULL); CHKERRQ(ierr);
	ierr 	 = PetscOptionsGetBool( PETSC_NULL, "-Benchmark_FallingBlock",   		&user->AnalyticalBenchmark, PETSC_NULL); CHKERRQ(ierr);
	ierr 	 = PetscOptionsGetBool( PETSC_NULL, "-Benchmark_ArcTanFallingBlock",   &user->AnalyticalBenchmark, PETSC_NULL); CHKERRQ(ierr);
	ierr 	 = PetscOptionsGetBool( PETSC_NULL, "-Benchmark_VerticalDensityCollumn",   &user->AnalyticalBenchmark, PETSC_NULL); CHKERRQ(ierr);

	user->GridAdvectionMethod  = 0; 	 // 0-Fully Eulerian, 1-Fully Lagrangian, 2-ALE with remeshing @ surface layer
	user->EulerianAfterTimestep=-1;		//	if >0, you can specify after which timestep to switch to eulerian mode
	user->FactorSurfaceLayer   = 0.2; 	 // how thick is the surface layer?
	user->num_subdt			   =  10;    // How many sub-timestep iterations if an ALE mode is selected?
	user->SaveParticles 	   =  0;	 // Save particles or not?
	user->num_phase_transitions= 0;
	user->InitialMeshFromFile  = 0;		 // In case you want to read an initial mesh from file
	user->InitialErosionSurfaceFromFile = 0; //
	user->InitialMantleLevel   = 10;     // for cases in which a lithosphere is modeled with hand-set crustal thickness
	user->dt_max               = 1e6;    // maximum timestep
	user->dt_temp              = 1e6;    // timestep to initialize temperature
	user->dt				   = 1e-3;	 // initial timestep
	user->Vx_Front             = 8.0;      // x-velocity at front boundary if BC==3 @ this boundary
	user->Vx_Back              = 0.0;      // x-velocity at back boundary if BC==3 @ this boundary
	user->Vy_Front             = 4.0;      // y-velocity at front boundary if BC==4 @ this boundary
	user->Vy_Back              = 0.0;      // y-velocity at back boundary if BC==4 @ this boundary
	user->Vy_Partx             = 0.5;    // Apply velocity condition at part of the x-domain, if BC=4
	user->Vy_Partz             = 0.0;      // Apply velocity condition at part of the x-domain, if BC=4

	user->MaximumSurfaceAngle  = 80.0;	 // maximum surface angle allowed
	user->MuMeanMethod         = 1;      // 0-compute mat. props @ integration points, 1-Arith. element-average, 2-Geom. element-average, 3-Harm. element-average
	user->NumPartX			   = 2;      // Number of particles per cell in x-direction
	user->NumPartY			   = 2;      // Number of particles per cell in y-direction
	user->NumPartZ			   = 2;      // Number of particles per cell in z-direction

	user->restart			   = 1; 	 // Are we restarting a simulation?
	user->save_breakpoints	   = 10;		 // After how many steps do we make a breakpoint?
	if (user->AnalyticalBenchmark){
		user->restart 		   = 0;		// don't restart and don't create breakpoint files if we perform an analytical benchmark
		user->save_breakpoints = 0;
	}
	user->break_point_number   = 0;		 // The number of the breakpoint file
	user->incr_breakpoints 	   = 100;		 // After how many steps should an incremental breakpointfile be created?
	user->fileno   			   = 0;

	user->LowerViscosityCutoff  = 1.0;		// lowermost viscosity cutoff in model
	user->UpperViscosityCutoff  = 1e30;		// uppermost viscosity cutoff in model

	user->num_phases 		   = 2;  // the default # of phases. In case we set props from the command line, and we use a folding setup combined with FDSTAG we might have to do something smarter (check the setup and increase this)

	user->GravityAngle 		   	= 	0.0;		// angle of gravity with z-axis (can be changed in the x-z plane)
	user->FSSA					=	0.0;		// Free Surface Stabilization Algorithm parameter used in FDSTAG, to stabilize an internal free surface [should be between 0-1]

	user->ArtTemp               = PETSC_FALSE;

	// Default erosion parameters
	user->ApplyErosion 		  	=  0; 	 	// 0-no; 1-cascade
	user->SurfaceAngle 			=  10.0;      // Angle that the surface makes with horizontal    [degree]
	user->SurfaceNoiseAmplitude =  100;  	// Amplitude of noise on free surface 				[m     ]
	user->fluvial_erosion 		=  1.6e-2; 	// Fluvial erosion rate
	user->diffusion_erosion 	=  2.0e-2; 	// Diffusion erosion rate
	user->baselevelx0 			=  1;      	//	Outlet at front [if 1; no outlet=0]
	user->baselevelx1 			=  0;      	//
	user->baselevely0 			=  0;     	//
	user->baselevely1 			=  0;      	//
	user->ErosionParameters.UseInternalFreeSurface 	= 0;	// don't use by default
	user->ErosionParameters.StickyAirPhase 			= 1;	// sticky air phase
	user->ErosionParameters.ErosionModel 			= 0;	// none by default
	user->ErosionParameters.SedimentationModel 		= 0;	// no sedimentation by default; 1-constant rate sedimentation

	//Default FD erosion code parameters
	user->ErosionParameters.FE_ErosionCode.ResolutionFactorX 		=	2;				// how much larger is the FD erosion code resolution compared to the mechanical code?
	user->ErosionParameters.FE_ErosionCode.ResolutionFactorY 		=	2;				// how much larger is the FD erosion code resolution compared to the mechanical code?
	user->ErosionParameters.FE_ErosionCode.InitialRandomNoise_m 	=	1;				// what is the amplitude of the random noise on the erosion surface [in meters]?
	user->ErosionParameters.FE_ErosionCode.InitialUpliftedSide_m 	=	100;			// By how much meters is the right side of the eroded surface uplifed [in meters]?
	user->ErosionParameters.FE_ErosionCode.dt 						=	100.0*SecYear;	// ideal erosion timestep
	user->ErosionParameters.FE_ErosionCode.rain_m_year 				=	0.3;			// rain in m/year
	user->ErosionParameters.FE_ErosionCode.k0 						=	3.2e-12;		// base erodability
	user->ErosionParameters.FE_ErosionCode.c 						=	1.0;			// prefactor for fluvial discharge
	user->ErosionParameters.FE_ErosionCode.n 						=	2.0;			// streampower exponent
    user->ErosionParameters.FE_ErosionCode.BC                       =   1; 
    user->ErosionParameters.FE_ErosionCode.fill_lake                =   0; 
    user->ErosionParameters.FE_ErosionCode.mode_river               =   0;
    user->ErosionParameters.FE_ErosionCode.nbre_river               =   1; 
    user->ErosionParameters.FE_ErosionCode.rain_river_year 			=	0.0;			// rain in m/year
    
    
	// Default Output parameters
	user->Output.velocity				=   1;
	user->Output.temperature			=   1;
	user->Output.surface_topography		=   1;
	user->Output.bottom_topography		=   1;
	user->Output.quadrature				=   1;

	// Default Pushing BC parameters
	user->AddPushing					=	0;
	user->Pushing.reset_pushing_coord	= 	0;
	user->Pushing.theta					= 	0.0;

	// Default gravity field parameters
	user->GravityField.GetIt		=	0;
	user->GravityField.SaveDebug	=	0;
	user->GravityField.SaveRef		=	0;
	user->GravityField.SaveVTK		=	0;
	user->GravityField.UseNumerics	= 	PETSC_FALSE;
	user->GravityField.UseAnalytics	= 	PETSC_FALSE;
	user->GravityField.num_intp 	=	2;
	user->GravityField.survey_nx 	=	11;
	user->GravityField.survey_ny 	=	11;
	user->GravityField.survey_xs 	=	0.0;
	user->GravityField.survey_xm 	=	1.0;
	user->GravityField.survey_ys 	=	0.0;
	user->GravityField.survey_ym 	=	1.0;
	user->GravityField.survey_z		=	0.0;
	user->GravityField.ReferenceDensity = 2700.0;
	user->GravityField.StdDev	=	1.0;
	sprintf(user->GravityField.RefDatFile2load, "ReferenceData/GravityField_REF.bin");
	user->GravityField.LithColNum      =   1;
	user->GravityField.LithColDens[0]  =   2670.0;
	user->GravityField.LithColDepth[0] =   30.0e3;

	// Default isostasy parameters
	user->Isostasy.GetIt            = 0;
	user->Isostasy.SaveRef          = 0;
	user->Isostasy.corr_topo        = 0.0;
	user->Isostasy.ref_rho          = 3215.0; // kg/m3
	user->Isostasy.ref_xi          = 0; // []
	user->Isostasy.ref_yi          = 0; // []

	// Default surface velocity field parameters
	user->SurfVelField.GetIt		=	0;
	user->SurfVelField.SaveRef		=	0;
	user->SurfVelField.VxStdDev		=	1.0;
	user->SurfVelField.VyStdDev		=	1.0;
	user->SurfVelField.VzStdDev		=	1.0;	
	sprintf(user->GravityField.RefDatFile2load, "ReferenceData/SurfVelField_REF.bin");


	// Default optimisation parameters
	user->Optimisation.GetIt		=	0;
	user->Optimisation.MisfitGravity=	0.0;
	user->Optimisation.MisfitSurfVel=	0.0;
	user->Optimisation.MisfitTiso   =	0.0;
	user->Optimisation.mpi_group_id	=	0;


	user->DA_Materials = PETSC_NULL;
	user->Materials    = PETSC_NULL;

	/* Initialize material properties */
	LaMEMInitializeMaterialProperties(user);

	/* Read an input file if required */
	if (user->InputParamFile){
		//	ReadInputFile(user);
		LaMEMReadInputFile(user);		// use the parser, to read the input file
	}

	/* Change values @ command prompt */

	PetscOptionsGetReal(PETSC_NULL,"-W"      ,	&user->W 			, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-L"      ,	&user->L 			, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-H"      ,	&user->H 			, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-He"      ,	&user->H 			, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-y_front" ,	&user->y_front 			, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-ampl2D" ,	&user->ampl2D 		, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-ampl3D" ,	&user->ampl3D 		, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-amplNoise",&user->amplNoise	, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-mumax"  ,	&user->mumax 		, PETSC_NULL);

	PetscOptionsGetInt(PETSC_NULL ,"-nel_x",	&user->nel_x 		, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL ,"-nel_y",	&user->nel_y 		, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL ,"-nel_z",	&user->nel_z 		, PETSC_NULL);

	nel_input_max=3;
	ierr = PetscOptionsGetIntArray(PETSC_NULL,"-nel", nel_array, &nel_input_max, &found); CHKERRQ(ierr);
	// this allows us to also specify the # of elements as -nel 8,16,32  which gives nel_x=8, nel_y=16, nel_z=32
	if (found==PETSC_TRUE) {
		user->nel_x = nel_array[0];
		user->nel_y = nel_array[1];
		user->nel_z = nel_array[2];
	}

	PetscOptionsGetInt(PETSC_NULL ,"-refinex",		&user->refinex		, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL ,"-refiney",		&user->refiney		, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL ,"-refinez",		&user->refinez		, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL ,"-time_end",	&user->time_end 	, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL ,"-time_end_temp",	&user->time_end_temp 	, PETSC_NULL);

	PetscOptionsGetInt(PETSC_NULL ,"-save_timesteps",	&user->save_timesteps 	, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-CFL"    ,		&user->CFL 			, PETSC_NULL);

	PetscOptionsGetReal(PETSC_NULL,"-dt"    ,		&user->dt 			, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-dt_max"    ,	&user->dt_max 			, PETSC_NULL);

	// FDSTAG Canonical Model Setup
	PetscOptionsGetString(PETSC_NULL,"-msetup", setup_name, PETSC_MAX_PATH_LEN, &found);
	if(found == PETSC_TRUE)
	{	if     (!strcmp(setup_name, "parallel"))   user->msetup = PARALLEL;
		else if(!strcmp(setup_name, "redundant"))  user->msetup = REDUNDANT;
		else if(!strcmp(setup_name, "diapir"))     user->msetup = DIAPIR;
		else if(!strcmp(setup_name, "block"))      user->msetup = BLOCK;
		else if(!strcmp(setup_name, "subduction")) user->msetup = SUBDUCTION;
		else if(!strcmp(setup_name, "folding"))    user->msetup = FOLDING;
		else if(!strcmp(setup_name, "detachment")) user->msetup = DETACHMENT;
		else if(!strcmp(setup_name, "slab"))       user->msetup = SLAB;
		else if(!strcmp(setup_name, "spheres"))    user->msetup = SPHERES;
		else if(!strcmp(setup_name, "bands"))      user->msetup = BANDS;
		else SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER,"#ERROR! Incorrect model setup: %s \n", setup_name);
	}

	PetscOptionsGetInt(PETSC_NULL ,"-Setup.Model",	&user->Setup.Model 	, PETSC_NULL);	// 		0-diapir, 1-single layer folding
	PetscOptionsGetInt(PETSC_NULL ,"-GridAdvectionMethod",	&user->GridAdvectionMethod 	, PETSC_NULL);	// 		0 - eulerian, 2-ALE, 1-Lagrangian
	PetscOptionsGetBool( PETSC_NULL,"-SkipStokesSolver",&user->SkipStokesSolver,PETSC_NULL );

	PetscOptionsGetInt(PETSC_NULL ,"-internalBC_coord",&user->internalBC_coord	, PETSC_NULL);	//
	PetscOptionsGetInt(PETSC_NULL ,"-internalBC_frontel",&user->internalBC_frontel	, PETSC_NULL);	//
	PetscOptionsGetInt(PETSC_NULL ,"-internalBC_backel",&user->internalBC_backel	, PETSC_NULL);	//
	PetscOptionsGetInt(PETSC_NULL ,"-internalBC_node",&user->internalBC_node	, PETSC_NULL);	//
	PetscOptionsGetInt(PETSC_NULL ,"-zdepth_BC_el",&user->zdepth_BC_el	, PETSC_NULL);	//
	PetscOptionsGetInt(PETSC_NULL ,"-zdepth_BC_node",&user->zdepth_BC_node	, PETSC_NULL);	//
	PetscOptionsGetInt(PETSC_NULL ,"-BC.InternalBound",&user->BC.InternalBound	, PETSC_NULL);	//
	PetscOptionsGetInt(PETSC_NULL ,"-BC.UpperBound",&user->BC.UpperBound	, PETSC_NULL);	//
	PetscOptionsGetInt(PETSC_NULL ,"-BC.LowerBound",&user->BC.LowerBound	, PETSC_NULL);	//
	PetscOptionsGetInt(PETSC_NULL ,"-BC.LeftBound",	&user->BC.LeftBound	, PETSC_NULL);	//
	PetscOptionsGetInt(PETSC_NULL ,"-BC.RightBound",&user->BC.RightBound	, PETSC_NULL);	//
	PetscOptionsGetInt(PETSC_NULL ,"-BC.FrontBound",&user->BC.FrontBound	, PETSC_NULL);	//
	PetscOptionsGetInt(PETSC_NULL ,"-BC.BackBound" ,&user->BC.BackBound		, PETSC_NULL);	//
	PetscOptionsGetReal(PETSC_NULL ,"-BC.Vy_front",	&user->BC.Vy_front 	, PETSC_NULL);	//  	y-velocity @ front boundary
	PetscOptionsGetReal(PETSC_NULL ,"-BC.Vy_back",	&user->BC.Vy_back 	, PETSC_NULL);	//  	y-velocity @ front boundary
	PetscOptionsGetReal(PETSC_NULL ,"-BC.Vz_top",	&user->BC.Vz_top 	, PETSC_NULL);	//  	y-velocity @ front boundary
	PetscOptionsGetReal(PETSC_NULL ,"-BC.Vz_bot",	&user->BC.Vz_bot 	, PETSC_NULL);	//  	y-velocity @ front boundary
	PetscOptionsGetReal(PETSC_NULL ,"-BC.Vx_left",	&user->BC.Vx_left 	, PETSC_NULL);	//  	y-velocity @ front boundary
	PetscOptionsGetReal(PETSC_NULL ,"-BC.Vx_right",	&user->BC.Vx_right 	, PETSC_NULL);	//  	y-velocity @ front boundary

	PetscOptionsGetReal(PETSC_NULL ,"-BC.Exx",			&user->BC.Exx 		, PETSC_NULL);	//  	y-velocity @ front boundary
	PetscOptionsGetReal(PETSC_NULL ,"-BC.Eyy",			&user->BC.Eyy 		, PETSC_NULL);	//  	y-velocity @ front boundary

	PetscOptionsGetReal(PETSC_NULL ,"-Vx_Front",	&user->Vx_Front		, PETSC_NULL);	// Vx_Front in cm/year
	PetscOptionsGetReal(PETSC_NULL ,"-Vx_Back",		&user->Vx_Back		, PETSC_NULL);	// Vx_Back in cm/year
	PetscOptionsGetReal(PETSC_NULL ,"-Vy_Front",	&user->Vx_Front		, PETSC_NULL);	// Vy_Front in cm/year
	PetscOptionsGetReal(PETSC_NULL ,"-Vy_Back",		&user->Vx_Back		, PETSC_NULL);	// Vy_Back in cm/year
	PetscOptionsGetReal(PETSC_NULL ,"-Vy_Partx",	&user->Vy_Partx		, PETSC_NULL);	// Vy_Partx for BC=4
	PetscOptionsGetReal(PETSC_NULL ,"-Vy_Partz",	&user->Vy_Partz		, PETSC_NULL);	// Vy_Partz for BC=4
	PetscOptionsGetInt(PETSC_NULL ,"-MuMeanMethod",	&user->MuMeanMethod	, PETSC_NULL);  // Specify element-averaging of effective material properties.
	PetscOptionsGetInt(PETSC_NULL ,"-NumPartX",	&user->NumPartX	, PETSC_NULL);  		//	# of tracers per cell in x-direction
	PetscOptionsGetInt(PETSC_NULL ,"-NumPartY",	&user->NumPartY	, PETSC_NULL);  		//	# of tracers per cell in y-direction
	PetscOptionsGetInt(PETSC_NULL ,"-NumPartZ",	&user->NumPartZ	, PETSC_NULL);  		//	# of tracers per cell in z-direction
	PetscOptionsGetInt(PETSC_NULL ,"-restart",	&user->restart	, PETSC_NULL);  		//	# restart a simulation if possible?
	PetscOptionsGetInt(PETSC_NULL ,"-fileno",	&user->fileno	, PETSC_NULL);  		//	# restart a simulation with file no ...
	PetscOptionsGetInt(PETSC_NULL ,"-save_breakpoints",	&user->save_breakpoints	, PETSC_NULL);  	//	After how many steps do we create a breakpoint file?
	PetscOptionsGetInt(PETSC_NULL ,"-ApplyErosion",	      &user->ApplyErosion		, PETSC_NULL);  //	Apply erosion or not?
	PetscOptionsGetReal(PETSC_NULL ,"-fluvial_erosion",	  &user->fluvial_erosion		, PETSC_NULL);	//  fluvial erosion rate
	PetscOptionsGetReal(PETSC_NULL ,"-diffusion_erosion", &user->diffusion_erosion			, PETSC_NULL);	//  diffusion erosion rate
	PetscOptionsGetReal(PETSC_NULL ,"-SurfaceNoiseAmplitude", &user->SurfaceNoiseAmplitude	, PETSC_NULL);	// amplitude of noise @ surface [m]
	PetscOptionsGetReal(PETSC_NULL ,"-Hinterface", &user->Hinterface						, PETSC_NULL);	// Horizontal interface for diapir setup
	PetscOptionsGetInt(PETSC_NULL ,"-InitialErosionSurfaceFromFile", &user->InitialErosionSurfaceFromFile, PETSC_NULL);	// load initial erosion surface from file

	PetscOptionsGetReal(PETSC_NULL ,"-NonlinearIterationsAccuracy", &user->NonlinearIterationsAccuracy	, PETSC_NULL);		// accuracy of nonlinear iterations
	PetscOptionsGetInt(PETSC_NULL ,"-MaxNonlinearIterations",		&user->MaxNonlinearIterations		, PETSC_NULL);  	//	maximum number of nonlinear iterations

	PetscOptionsGetReal(PETSC_NULL ,"-LowerViscosityCutoff", &user->LowerViscosityCutoff	, PETSC_NULL);		// lower viscosity cutoff
	PetscOptionsGetReal(PETSC_NULL ,"-UpperViscosityCutoff", &user->UpperViscosityCutoff	, PETSC_NULL);		// upper viscosity cutoff

	PetscOptionsGetReal(PETSC_NULL ,"-GravityAngle", &user->GravityAngle	, PETSC_NULL);		// Gravity angle in x-z plane


	PetscOptionsGetInt(PETSC_NULL ,"-incr_breakpoints",		&user->incr_breakpoints		, PETSC_NULL);  	//	save incr_breakpoints breakpoints
	PetscOptionsGetReal(PETSC_NULL ,"-MaximumSurfaceAngle", &user->MaximumSurfaceAngle	, PETSC_NULL);		// MaximumSurfaceAngle

	PetscOptionsGetInt(PETSC_NULL ,"-ParticleInjectionPhase",		&user->ParticleInjectionPhase	, PETSC_NULL);  		//	which phase do we inject?
	PetscOptionsGetInt(PETSC_NULL ,"-MatlabOutputFiles",			&user->MatlabOutputFiles		, PETSC_NULL);  		//	write matlab output or not?
	PetscOptionsGetInt(PETSC_NULL ,"-VTKOutputFiles",				&user->VTKOutputFiles			, PETSC_NULL);  		//	write VTK output or not?
	PetscOptionsGetInt(PETSC_NULL ,"-AVDPhaseViewer",				&user->AVDPhaseViewer			, PETSC_NULL);  		//	write AVDPhase output or not?
    PetscOptionsGetBool(PETSC_NULL ,"-SavePartitioning",			&user->SavePartitioning        	, PETSC_NULL);

    
	PetscOptionsGetInt(PETSC_NULL ,"-LoadInitialParticlesFromDisc",	&user->LoadInitialParticlesFromDisc			, PETSC_NULL);  		//	Load initial particles from file

	PetscOptionsGetReal(PETSC_NULL ,"-CriticalDiagonalRatio", 		&user->CriticalDiagonalRatio 	, PETSC_NULL);		// CriticalDiagonalRatio to induce remeshing if LaMEM is run in a lagrangian mode
	PetscOptionsGetInt(PETSC_NULL ,"-SaveParticles",				&user->SaveParticles			, PETSC_NULL);  		//	save particles to disk?
	PetscOptionsGetInt(PETSC_NULL ,"-EulerianAfterTimestep",				&user->EulerianAfterTimestep			, PETSC_NULL);  		//	switch to Eulerian mode after timestep ?? - decativated if <0 (default)

	PetscOptionsGetInt(PETSC_NULL  ,"-UseInternalFreeSurface",		&user->ErosionParameters.UseInternalFreeSurface	, PETSC_NULL);  	//	use the internal free surface
	PetscOptionsGetInt(PETSC_NULL  ,"-ErosionModel",				&user->ErosionParameters.ErosionModel	, PETSC_NULL);  			//	which erosion model do we employ? [0-none; 1-fast]


	/*set the FD erosion parameters from the command-line */
	PetscOptionsGetInt(PETSC_NULL  ,"-FE_ErosionCode.ResolutionFactorX",			&user->ErosionParameters.FE_ErosionCode.ResolutionFactorX		, PETSC_NULL);  			// how much larger is the FD erosion code resolution compared to the mechanical code?
	PetscOptionsGetInt(PETSC_NULL  ,"-FE_ErosionCode.ResolutionFactorY",			&user->ErosionParameters.FE_ErosionCode.ResolutionFactorY		, PETSC_NULL);  			// how much larger is the FD erosion code resolution compared to the mechanical code?
	PetscOptionsGetReal(PETSC_NULL ,"-FE_ErosionCode.InitialRandomNoise_m", 	&user->ErosionParameters.FE_ErosionCode.InitialRandomNoise_m 	, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL ,"-FE_ErosionCode.InitialUpliftedSide_m", 	&user->ErosionParameters.FE_ErosionCode.InitialUpliftedSide_m 	, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL ,"-FE_ErosionCode.dt", 						&user->ErosionParameters.FE_ErosionCode.dt 						, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL ,"-FE_ErosionCode.rain_m_year", 				&user->ErosionParameters.FE_ErosionCode.rain_m_year 			, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL ,"-FE_ErosionCode.k0", 						&user->ErosionParameters.FE_ErosionCode.k0 						, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL ,"-FE_ErosionCode.c", 						&user->ErosionParameters.FE_ErosionCode.c 						, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL ,"-FE_ErosionCode.n", 						&user->ErosionParameters.FE_ErosionCode.n 						, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL  ,"-FE_ErosionCode.BC",                       &user->ErosionParameters.FE_ErosionCode.BC                      , PETSC_NULL);
    n_int = 100;
    found = PETSC_FALSE;
    PetscOptionsGetInt(PETSC_NULL  ,"-FE_ErosionCode.mode_river",               &user->ErosionParameters.FE_ErosionCode.mode_river              , &found);
    PetscOptionsGetInt(PETSC_NULL  ,"-FE_ErosionCode.fill_lake",                &user->ErosionParameters.FE_ErosionCode.fill_lake               , PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL  ,"-FE_ErosionCode.nbre_river",               &user->ErosionParameters.FE_ErosionCode.nbre_river              , PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL ,"-FE_ErosionCode.rain_river_year", 			&user->ErosionParameters.FE_ErosionCode.rain_river_year 		, PETSC_NULL);

    PetscOptionsGetRealArray(PETSC_NULL,"-FE_ErosionCode.location_river",user->ErosionParameters.FE_ErosionCode.location_river,&n_int, &flg);         
    if (found){
    	/* print info if we select river parameters */

    	printf("==========================================================================================\n");
    	printf("in file utils.c \n");
    	printf("mode_river:%d, nbre_river:%d \n", user->ErosionParameters.FE_ErosionCode.mode_river, user->ErosionParameters.FE_ErosionCode.nbre_river);


    	for( i=0; i<user->ErosionParameters.FE_ErosionCode.nbre_river; i++ ) {
    		printf("location_river[%d]:%f\n",i,user->ErosionParameters.FE_ErosionCode.location_river[i]);
    	}

    	printf("==========================================================================================\n");
    }

    
    
	PetscOptionsGetInt(PETSC_NULL  ,"-StickyAirPhase",				&user->ErosionParameters.StickyAirPhase			, PETSC_NULL);  														 //	which phase is sticky air?
	PetscOptionsGetReal(PETSC_NULL ,"-FSSA", 						&user->FSSA 									, PETSC_NULL);		// FSSA parameter [should be between 0-1]
	PetscOptionsGetBool(PETSC_NULL, "-ArtificialTemperature", 	    &user->ArtTemp,    PETSC_NULL );

	PetscOptionsGetInt(PETSC_NULL  ,"-SedimentationModel",			&user->ErosionParameters.SedimentationModel	, PETSC_NULL);  			//	which sedimentation model do we employ? [0-none; 1-constant rate]
	PetscOptionsGetReal(PETSC_NULL ,"-InitialFreeSurfaceHeight", 	&user->ErosionParameters.InitialFreeSurfaceHeight 		, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL ,"-SedimentationRate_cmYr", 		&user->ErosionParameters.SedimentationRate_cmYr 		, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL ,"-SedimentLayerThicknessYears", &user->ErosionParameters.SedimentLayerThicknessYears 	, PETSC_NULL);  // we change the sediment phase every ?? years (mainly a visualization issue)
	PetscOptionsGetInt(PETSC_NULL  ,"-PhaseFirstSedimentedLayer",	&user->ErosionParameters.PhaseFirstSedimentedLayer		, PETSC_NULL);  // which is the first sediment layer?
	PetscOptionsGetInt(PETSC_NULL  ,"-PhaseLastSedimentedLayer",	&user->ErosionParameters.PhaseLastSedimentedLayer		, PETSC_NULL);  // which is the last sediment layer? After this, we sediment the first one again

	// --- SurfaceVelocity related input parameters ---
	PetscOptionsGetInt( PETSC_NULL ,"-get_SurfVelField", 			&user->SurfVelField.GetIt, 		PETSC_NULL );
	PetscOptionsGetInt( PETSC_NULL ,"-SurfVelField_SaveRef", 		&user->SurfVelField.SaveRef,	PETSC_NULL );
	PetscOptionsGetString(PETSC_NULL,"-SurfVelField_RefDatFile",  *(&user->SurfVelField.RefDatFile2load),PETSC_MAX_PATH_LEN,PETSC_NULL);
	PetscOptionsGetReal( PETSC_NULL,"-SurfVelField_VxStdDev",         &user->SurfVelField.VxStdDev,PETSC_NULL );	
	PetscOptionsGetReal( PETSC_NULL,"-SurfVelField_VyStdDev",         &user->SurfVelField.VyStdDev,PETSC_NULL );	
	PetscOptionsGetReal( PETSC_NULL,"-SurfVelField_VzStdDev",         &user->SurfVelField.VzStdDev,PETSC_NULL );		
	 // --- Optimization related input parameters ---
	PetscOptionsGetInt( PETSC_NULL ,"-get_Misfit", 					&user->Optimisation.GetIt,		PETSC_NULL );

	// --- Gravity related input parameters ---
	PetscOptionsGetInt( PETSC_NULL ,"-get_GravityField",			&user->GravityField.GetIt, 	   	PETSC_NULL );
	PetscOptionsGetInt( PETSC_NULL ,"-GravityField_SaveRef", 		&user->GravityField.SaveRef,	PETSC_NULL );
	PetscOptionsGetInt( PETSC_NULL ,"-GravityField_SaveDebug",		&user->GravityField.SaveDebug,	PETSC_NULL );
	PetscOptionsGetInt( PETSC_NULL ,"-GravityField_SaveVTK",		&user->GravityField.SaveVTK, 	PETSC_NULL );
	PetscOptionsGetInt( PETSC_NULL ,"-GravityField_survey_nx",		&user->GravityField.survey_nx, 	PETSC_NULL );
	PetscOptionsGetInt( PETSC_NULL ,"-GravityField_survey_ny", 		&user->GravityField.survey_ny, 	PETSC_NULL );
	PetscOptionsGetReal( PETSC_NULL,"-GravityField_survey_xs", 		&user->GravityField.survey_xs, 	PETSC_NULL );
	PetscOptionsGetReal( PETSC_NULL,"-GravityField_survey_xm", 		&user->GravityField.survey_xm, 	PETSC_NULL );
	PetscOptionsGetReal( PETSC_NULL,"-GravityField_survey_ys",	 	&user->GravityField.survey_ys, 	PETSC_NULL );
	PetscOptionsGetReal( PETSC_NULL,"-GravityField_survey_ym", 		&user->GravityField.survey_ym, 	PETSC_NULL );
	PetscOptionsGetReal( PETSC_NULL,"-GravityField_survey_z", 		&user->GravityField.survey_z,  	PETSC_NULL );
	PetscOptionsGetReal( PETSC_NULL,"-GravityField_ReferenceDensity",&user->GravityField.ReferenceDensity,PETSC_NULL );	
	PetscOptionsGetInt( PETSC_NULL ,"-GravityField_num_intp", 		&user->GravityField.num_intp,  	PETSC_NULL );
	PetscOptionsGetBool( PETSC_NULL,"-GravityField_UseNumerics",	&user->GravityField.UseNumerics,PETSC_NULL );
	PetscOptionsGetBool( PETSC_NULL,"-GravityField_UseAnalytics",	&user->GravityField.UseAnalytics,PETSC_NULL );
	PetscOptionsGetString(PETSC_NULL,"-GravityField_RefDatFile",  *(&user->GravityField.RefDatFile2load),PETSC_MAX_PATH_LEN,PETSC_NULL);
	PetscOptionsGetReal( PETSC_NULL,"-GravityField_StdDev",         &user->GravityField.StdDev,PETSC_NULL );	
	ierr = GetLithColumnFromCommandLine(user);CHKERRQ(ierr);

	// --- Isostasy related parameters ---
	PetscOptionsGetInt( PETSC_NULL ,"-ComputeAiryIsostasy", 		&user->Isostasy.GetIt,   	PETSC_NULL );

	// --- Output related parameters ---
	PetscOptionsGetInt(PETSC_NULL,"-Output.velocity" 		    ,	&user->Output.velocity				, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL,"-Output.temperature" 		,	&user->Output.temperature			, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL,"-Output.surface_topography" 	,	&user->Output.surface_topography	, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL,"-Output.bottom_topography" 	,	&user->Output.surface_topography	, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL,"-Output.quadrature" 			,	&user->Output.quadrature			, PETSC_NULL);

	// --- Pushing BC related parameters ---
	PetscOptionsGetInt(PETSC_NULL,"-AddPushing" 			 ,	&user->AddPushing				, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL,"-Pushing.num_changes"	 ,	&user->Pushing.num_changes		, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL,"-Pushing.reset_pushing_coord",	&user->Pushing.reset_pushing_coord		, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-Pushing.L_block"      	 ,	&user->Pushing.L_block 			, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-Pushing.W_block"      	 ,	&user->Pushing.W_block 			, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-Pushing.H_block"      	 ,	&user->Pushing.H_block 			, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-Pushing.x_center_block" ,	&user->Pushing.x_center_block 	, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-Pushing.y_center_block" ,	&user->Pushing.y_center_block 	, PETSC_NULL);
	PetscOptionsGetReal(PETSC_NULL,"-Pushing.z_center_block" ,	&user->Pushing.z_center_block 	, PETSC_NULL);

	
	char 			matprop_opt[PETSC_MAX_PATH_LEN];
	flg = PETSC_FALSE;

	for(i=0;i<user->Pushing.num_changes;i++){
		// V_push in cm/year
		sprintf(matprop_opt,"-Pushing.V_push_%lld",(LLD)i);
		ierr = PetscOptionsGetReal(PETSC_NULL ,matprop_opt,&user->Pushing.V_push[i]	, &flg); 				CHKERRQ(ierr);
		if(flg == PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD,"#    V_push[%lld]	= %g \n",(LLD)i,user->Pushing.V_push[i]);

		// Rate of rotation deg/yr
		sprintf(matprop_opt,"-Pushing.omega_%lld",(LLD)i);
		ierr = PetscOptionsGetReal(PETSC_NULL ,matprop_opt,&user->Pushing.omega[i]	, &flg); 				CHKERRQ(ierr);
		if(flg == PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD,"#    Omega[%lld]	= %g \n",(LLD)i,user->Pushing.omega[i]);

		// Options whether to advect block for each time segment
		sprintf(matprop_opt,"-Pushing.coord_advect_%lld",(LLD)i);
		ierr = PetscOptionsGetInt(PETSC_NULL ,matprop_opt,&user->Pushing.coord_advect[i]	, &flg); 				CHKERRQ(ierr);
		if(flg == PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD,"#    Coord_advect[%lld]	= %g \n",(LLD)i,user->Pushing.coord_advect[i]);

		// Options to define the direction of pushing
				sprintf(matprop_opt,"-Pushing.dir_%lld",(LLD)i);
				ierr = PetscOptionsGetInt(PETSC_NULL ,matprop_opt,&user->Pushing.dir[i]	, &flg); 				CHKERRQ(ierr);
				if(flg == PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD,"#    Direction[%lld]	= %g \n",(LLD)i,user->Pushing.dir[i]);
	}

	for(i=0;i<user->Pushing.num_changes+1;i++){
		//Time is a num_changes+1 array
		sprintf(matprop_opt,"-Pushing.time_%lld",(LLD)i);
		ierr = PetscOptionsGetReal(PETSC_NULL ,matprop_opt,&user->Pushing.time[i]	, &flg); 				CHKERRQ(ierr);
		if(flg == PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD,"#    time[%lld]	= %g \n",(LLD)i,user->Pushing.time[i]);
	}

	// linear solver options
	PetscOptionsGetInt( PETSC_NULL, "-StokesSolver", 				&user->StokesSolver, 			PETSC_NULL );
	PetscOptionsGetInt( PETSC_NULL, "-VelocitySolver", 				&user->VelocitySolver, 			PETSC_NULL );
	PetscOptionsGetBool( PETSC_NULL, "-VelocityTest", 				&user->VelocityTest, 			PETSC_NULL );
	PetscOptionsGetBool( PETSC_NULL, "-ScaleSystem", 				&user->ScaleSystem, 			PETSC_NULL );
	PetscOptionsGetBool( PETSC_NULL, "-use_fdstag_canonical", 	    &user->use_fdstag_canonical,    PETSC_NULL );

	// --- Get material properties from command line
	ierr = GetMaterialPropertiesFromCommandLine(user);

	/* In case we use FDSTAG, we need at least 2 elements (not 1)! */
	// If FDSTAG is used, dt must be different to 0 to avoid zero division at "ComputeStiffnessMatrixRHSTemperature_FDSTAG"
	if( __ELEMENT_TYPE__ == ELEMENT_FDSTAG ){
		if (user->nel_x==1){
			PetscPrintf(PETSC_COMM_WORLD," With FDSTAG we need at least 2 elements in the x-direction! I have increased this. \n");
			user->nel_x=2;
		}
		if (user->nel_y==1){
			PetscPrintf(PETSC_COMM_WORLD," With FDSTAG we need at least 2 elements in the y-direction! I have increased this. \n");
			user->nel_y=2;
		}
		if (user->dt==0.0){
			PetscPrintf(PETSC_COMM_WORLD," With FDSTAG dt cannot be zero. Value has been changed to dt = 1e-3. \n");
			user->dt=1e-3;
		}
	}

	// Ensure that sticky air phase is not larger than the max. number of phases
	if (user->ErosionParameters.StickyAirPhase > (user->num_phases-1)){
		PetscPrintf(PETSC_COMM_WORLD," The sticky air phase is %i but the maximum phase in the model setup is %i. I changed sticky air phase to %i. \n",user->ErosionParameters.StickyAirPhase, (user->num_phases-1),(user->num_phases-1));
		user->ErosionParameters.StickyAirPhase = (user->num_phases-1);
	}


	/* In case nel_x, nel_y or nel_z is specified, recompute the number of nodes here; this ALWAYS overrides a specification of nnode_x,nnode_y and nnode_z */
	if (user->nel_x > 0){
		if( (__ELEMENT_TYPE__ == ELEMENT_Q2P1) ) {
			user->nnode_x = user->nel_x*2+1;
		}
		else {
			user->nnode_x = user->nel_x+1;
		}
	}
	if (user->nel_y > 0){
		if( (__ELEMENT_TYPE__ == ELEMENT_Q2P1) ) {
			user->nnode_y = user->nel_y*2+1;
		}
		else {
			user->nnode_y = user->nel_y+1;
		}
	}
	if (user->nel_z > 0){
		if( (__ELEMENT_TYPE__ == ELEMENT_Q2P1) ) {
			user->nnode_z = user->nel_z*2+1;
		}
		else {
			user->nnode_z = user->nel_z+1;
		}
	}

	/* If we really insist, we can specify nnode_x from the command line */
	PetscOptionsGetInt(PETSC_NULL ,"-nnode_x",	&user->nnode_x 		, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL ,"-nnode_y",	&user->nnode_y 		, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL ,"-nnode_z",	&user->nnode_z 		, PETSC_NULL);

	user->Setup.Diapir_Hi 	 = (user->H-user->z_bot)*user->Hinterface;	// for 0-diapir setup
	user->Setup.SingleFold_H = 1.0;						    			// for  1-Single-layer folding setup

	// only mess with legacy linear solver settings
	if(user->use_fdstag_canonical != PETSC_TRUE)
	{
		ierr = InitializeLinearSolverLegacy( user ); CHKERRQ(ierr);
	}

	// set this option to monitor actual option usage
//	PetscOptionsInsertString("-options_left");

	if (!(user->InputParamFile)){

		ComputeCharacteristicValues(user);

		/* Define phase Material properties in case NO input file is specified (=falling-block test)*/
		/* viscosity                                     								 density 									*/
		user->PhaseProperties.mu[0]  = 1.0*user->Characteristic.Viscosity; 		user->PhaseProperties.rho[0] = 1.0*user->Characteristic.Density;
		user->PhaseProperties.mu[1] = user->mumax; 								user->PhaseProperties.rho[1] = 2.0*user->Characteristic.Density;

		PetscOptionsInsertString("-AddRandomNoiseParticles 0");		// will always give the same results
	}



	{
		// Define material parameters for folding benchmarks, which overrule values from the parameters file
		PetscScalar	R, nL, nM;
		R 	= 	100;
		nL 	=	1;
		nM 	=	1;
		PetscOptionsGetReal(PETSC_NULL ,"-FoldingBenchmark_R",  &R	, &found);	// Viscosity contrast
		if (found== PETSC_TRUE){
			user->PhaseProperties.mu[1] 			=	R;
			user->PhaseProperties.mu[0] 			=	1;
		}
		PetscOptionsGetReal(PETSC_NULL ,"-FoldingBenchmark_nL", &nL	, &found);	// powerlaw exponent of layer
		if (found== PETSC_TRUE){
			user->PhaseProperties.n_exponent[1] 	=	nL;
		}
		PetscOptionsGetReal(PETSC_NULL ,"-FoldingBenchmark_nM", &nM	, &found);	// powerlaw exponent of matrix
		if (found== PETSC_TRUE){
			user->PhaseProperties.n_exponent[0] 	=	nM;
		}

	}

	/* print information about simulation */
	nx = user->nnode_x;
	ny = user->nnode_y;
	nz = user->nnode_z;
	dx = user->W/((double)(nx-1));
	dy = user->L/((double)(ny-1));
	dz = user->H/((double)(nz-1));

	MPI_Comm_size(PETSC_COMM_WORLD,&size);

	PetscPrintf(PETSC_COMM_WORLD," Element-type used         : ");

	if( __ELEMENT_TYPE__ == ELEMENT_Q1P0 ){
		PetscPrintf(PETSC_COMM_WORLD,"Q1P0 ");
	}
	if( __ELEMENT_TYPE__ == ELEMENT_Q1Q1 ){
		PetscPrintf(PETSC_COMM_WORLD,"Q1Q1 ");
	}
	if( (__ELEMENT_TYPE__ == ELEMENT_Q2P1) ) {
		if ( __Q2_TYPE__ == ELEMENT_Q2P1_LOCAL){
			PetscPrintf(PETSC_COMM_WORLD,"Q2Pm1_local");
		}
		if ( __Q2_TYPE__ == ELEMENT_Q2P1_GLOBAL){
			PetscPrintf(PETSC_COMM_WORLD,"Q2Pm1_global");
		}
	}
	if( __ELEMENT_TYPE__ == ELEMENT_FDSTAG ){
		PetscPrintf(PETSC_COMM_WORLD,"FDSTAG");
		user->GridAdvectionMethod  = 0; 			// only works with Eulerian meshes
	}
	PetscPrintf(PETSC_COMM_WORLD,";  [Q1P0 Q2Pm1_global Q2Pm1_local Q1Q1 FDSTAG] (change with -vpt_element)\n");

	PetscPrintf(PETSC_COMM_WORLD," Total # of cpu's          : %lld \n",(LLD)size);
	PetscPrintf(PETSC_COMM_WORLD," Resolution [nx,ny,nz]     : %lld x %lld x %lld \n",(LLD)(nx), (LLD)(ny), (LLD)(nz));
	PetscPrintf(PETSC_COMM_WORLD," Total # of velocity dof's : %lld \n",(LLD)(nx*ny*nz*3));

	/* Info about material averaging method employed */
	if (user->MuMeanMethod==0){PetscPrintf(PETSC_COMM_WORLD," Computing material properties at integration points \n"); }
	if (user->MuMeanMethod==1){PetscPrintf(PETSC_COMM_WORLD," Computing mean material per element using arithmetic averaging from integration points \n"); }
	if (user->MuMeanMethod==2){PetscPrintf(PETSC_COMM_WORLD," Computing mean material per element using geometric averaging from integration points \n"); }
	if (user->MuMeanMethod==3){PetscPrintf(PETSC_COMM_WORLD," Computing mean material per element using harmonic  averaging from integration points \n"); }

	/* Info about grid movement method */
	if (user->GridAdvectionMethod==0){PetscPrintf(PETSC_COMM_WORLD," Employing code in an Eulerian manner \n"); }
	if (user->GridAdvectionMethod==1){PetscPrintf(PETSC_COMM_WORLD," Employing code in a Lagrangian manner \n"); }
	if (user->GridAdvectionMethod==2){PetscPrintf(PETSC_COMM_WORLD," Employing code in an ALE manner, with remeshing near surface \n"); }

	if (user->GravityAngle!=90){
		PetscPrintf(PETSC_COMM_WORLD," Gravity angle with z-axis =  %g \n",user->GravityAngle);
	}

	/* Info about approximate aspect ratio [useful for multigrid, who would like to have aspect ratios ~1] */
	PetscPrintf(PETSC_COMM_WORLD," Approximate Aspect ratio of cells [dx/dz, dy/dz] = [%g %g] \n", dx/dz, dy/dz);

	/* Info about particles if used */
#ifdef PARTICLES
	if (user->ParticleInput==1){
		/* Check whether a sufficient amount of particles are specified; if not increase it.
		 * We need at least as many particles as integration points; thus 2x2x2 for Q1 elements and 3x3x3 for Q2 elements
		 */

		if( (__ELEMENT_TYPE__ == ELEMENT_Q2P1) ){
			if (user->NumPartX<3){ user->NumPartX = 3; PetscPrintf(PETSC_COMM_WORLD," WARNING: At least 3 particles required in x-direction for Q2 elements; I have increased this \n");	}
			if (user->NumPartY<3){ user->NumPartY = 3; PetscPrintf(PETSC_COMM_WORLD," WARNING: At least 3 particles required in y-direction for Q2 elements; I have increased this \n");	}
			if (user->NumPartZ<3){ user->NumPartZ = 3; PetscPrintf(PETSC_COMM_WORLD," WARNING: At least 3 particles required in z-direction for Q2 elements; I have increased this \n");	}
		}
		else if ((__ELEMENT_TYPE__ == ELEMENT_Q1Q1) || (__ELEMENT_TYPE__ == ELEMENT_Q1P0)){
			if (user->NumPartX<2){ user->NumPartX = 2; PetscPrintf(PETSC_COMM_WORLD," WARNING: At least 2 particles required in x-direction for Q1 elements; I have increased this \n");	}
			if (user->NumPartY<2){ user->NumPartY = 2; PetscPrintf(PETSC_COMM_WORLD," WARNING: At least 2 particles required in y-direction for Q1 elements; I have increased this \n");	}
			if (user->NumPartZ<2){ user->NumPartZ = 2; PetscPrintf(PETSC_COMM_WORLD," WARNING: At least 2 particles required in z-direction for Q1 elements; I have increased this \n");	}
		}
	}

	PetscPrintf(PETSC_COMM_WORLD," Number of tracers/cell in [x,y,z]-direction = [%lld,%lld,%lld] \n",(LLD)(user->NumPartX),(LLD)(user->NumPartY),(LLD)(user->NumPartZ));
#endif


	/* Compute some useful stuff */
	user->NumSurfaceNodes = ((PetscInt) (user->FactorSurfaceLayer* ((double) nz)) );
	LaMEMMod(user->NumSurfaceNodes, 2, &mod); if (mod==0){user->NumSurfaceNodes = user->NumSurfaceNodes-1;};
	PetscPrintf(PETSC_COMM_WORLD," NumSurfaceNodes=%lld  nz=%lld \n", (LLD)(user->NumSurfaceNodes), (LLD)nz);

	if (user->GridAdvectionMethod != 2  || user->BC.UpperBound != 0){
		user->num_subdt			  = 1;
	}
	PetscOptionsGetInt(PETSC_NULL ,"-num_subdt",	&user->num_subdt	, PETSC_NULL);  		//	# modify number of sub-dt

	if (user->GridAdvectionMethod==2 && user->FactorSurfaceLayer>0 && user->NumSurfaceNodes<2 ){
		PetscPrintf(PETSC_COMM_WORLD," ***** You selected ALE grid advection, with a surface layer   *****\n");
		PetscPrintf(PETSC_COMM_WORLD," ***** But the surface layer has 1 node only, which is insufficient  *****\n");
		PetscPrintf(PETSC_COMM_WORLD," ***** Increasing value to 3  *****\n");
		user->NumSurfaceNodes=3;
		PetscPrintf(PETSC_COMM_WORLD," \n");
	}

	if  (user->ApplyErosion==0){
		// non zero surface angle is only relevant if erosion is applied to the model
		user->SurfaceAngle 		   = 0;
		user->SurfaceNoiseAmplitude= 0;
	}

	/* Check whether nonlinear iterations are required */
	user->NonlinearIterations = 0;
	for ( i=0; i<user->num_phases; i++ ){
		if ( ((user->PhaseProperties.ViscosityLaw[i] == 4) || (user->PhaseProperties.ViscosityLaw[i] == 2)) && (user->PhaseProperties.n_exponent[i]>1.0) ){
			/* One of the phases has powerlaw rheology with n>1 */
			user->NonlinearIterations = 1;
		}
		if ( ((user->PhaseProperties.ViscosityLaw[i] == 4) || (user->PhaseProperties.ViscosityLaw[i] == 2)) && (user->PhaseProperties.FrictionAngle[i]>0.0) ){
			/* One of the phases has powerlaw rheology and we likely have plasticity acting */
			user->NonlinearIterations = 1;
		}
		if ( ((user->PhaseProperties.ViscosityLaw[i] == 4) || (user->PhaseProperties.ViscosityLaw[i] == 2)) && (user->PhaseProperties.Cohesion[i]<1e100) ){
			/* One of the phases has powerlaw rheology with n>1 */
			user->NonlinearIterations = 1;
		}
	}


	if (	user->NonlinearIterations == 1){
		PetscPrintf(PETSC_COMM_WORLD," # Employing nonlinear iterations \n");
	}
	else{
		PetscPrintf(PETSC_COMM_WORLD," # Not employing nonlinear iterations \n");
	}


	/* Check for periodic boundary conditions */
	if ( (user->BC.LeftBound==3) || (user->BC.RightBound==3) ){
		user->BC.LeftBound=3;  user->BC.RightBound=3;			// if one is periodic, they both should be
		user->BC.BCType_x = DM_BOUNDARY_PERIODIC;
	}
	else{
		user->BC.BCType_x = DM_BOUNDARY_NONE;
	}

	if ( (user->BC.FrontBound==3) || (user->BC.BackBound==3) ){
		user->BC.FrontBound=3;  user->BC.BackBound=3;			// if one is periodic, they both should be
		user->BC.BCType_y = DM_BOUNDARY_PERIODIC;
	}
	else{
		user->BC.BCType_y = DM_BOUNDARY_NONE;
	}

	if ( (user->BC.UpperBound==3) || (user->BC.LowerBound==3) ){
		PetscPrintf(PETSC_COMM_WORLD," Periodic BCs are not yet implemented in these directions! \n");
		MPI_Abort(PETSC_COMM_WORLD,1);
	}
	else{
		user->BC.BCType_z = DM_BOUNDARY_NONE;
	}

	PetscPrintf(PETSC_COMM_WORLD," BC employed: BC.[Left Right; Front Back; Lower Upper]=[%lld %lld; %lld %lld; %lld %lld] \n",
			(LLD)(user->BC.LeftBound), (LLD)(user->BC.RightBound), (LLD)(user->BC.FrontBound), (LLD)(user->BC.BackBound), (LLD)(user->BC.LowerBound), (LLD)(user->BC.UpperBound) );


	/* Set the density of 'sticky-air' to zero if we eliminate the 'sticky-air' from the system */
	if (user->ErosionParameters.UseInternalFreeSurface==1){
		PetscInt 	AirPhase;
		PetscBool	eliminate_stickyair_from_system;

		eliminate_stickyair_from_system = PETSC_FALSE;
		ierr 		= 	PetscOptionsGetBool( PETSC_NULL, "-eliminate_stickyair_from_system", &eliminate_stickyair_from_system,PETSC_NULL ); CHKERRQ(ierr);
		AirPhase 	= 	user->ErosionParameters.StickyAirPhase; 		// sticky air phase

		if ((eliminate_stickyair_from_system) && (user->PhaseProperties.rho[AirPhase]>0)){

			PetscPrintf(PETSC_COMM_WORLD," You eliminate the sticky air phase from the system, but the sticky air density is %f kg/m3 which is inconsistent. I have set air density to 0 kg/m3 !! \n", user->PhaseProperties.rho[AirPhase]);
			user->PhaseProperties.rho[AirPhase]	=	0.0;		// set to zero


		}

	}

	// show initial timestep
	PetscPrintf(PETSC_COMM_WORLD,"# Initial time step: %g years \n",user->dt);

	/* Compute characteristic values */
	ComputeCharacteristicValues(user);

	/* Perform non-dimensionalisation of all input parameters */
	PerformNonDimensionalization(user);

	/* Allocate arrays required for time-dependent data */
	ierr = PetscMalloc( (size_t)user->time_end*sizeof(GlobalTimeDependentData), 	&user->TimeDependentData); CHKERRQ(ierr);
	ierr = PetscMemzero(user->TimeDependentData, (size_t)user->time_end*sizeof(GlobalTimeDependentData)); CHKERRQ(ierr);

	user->ErosionParameters.HorizontalFreeSurfaceHeight = user->ErosionParameters.InitialFreeSurfaceHeight;

	ierr = PetscOptionsInsertString( all_options ); CHKERRQ(ierr); /* force command line args in to override defaults */

	ierr = PetscFree(all_options); CHKERRQ(ierr);

	PetscFunctionReturn(0);

}
/*==========================================================================================================*/
#undef __FUNCT__
#define __FUNCT__ "InitializeLinearSolverLegacy"
PetscErrorCode InitializeLinearSolverLegacy( UserContext *user )
{
	PetscBool SolverDiagnostics;
	PetscMPIInt    nproc;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = MPI_Comm_size( PETSC_COMM_WORLD, &nproc ); CHKERRQ(ierr);

	// Do we want to print info about the solvers?
	ierr = PetscOptionsHasName(PETSC_NULL,"-SolverDiagnostics",&SolverDiagnostics); CHKERRQ(ierr);

	// check solver types selected
	if(user->StokesSolver   < 1 || user->StokesSolver   > 4) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Invalid Stokes solver type");
	if(user->VelocitySolver < 0 || user->VelocitySolver > 5) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Invalid velocity solver type");

	//==========================================
	// Set default options for the Stokes solver
	//==========================================
	if(user->StokesSolver == 1)
	{
		// 1 - Powell-Hesteness iterations
		PetscPrintf(PETSC_COMM_WORLD," Stokes solver             : Powell-Hestenes\n");
		// check velocity solver restrictions
		if(user->VelocitySolver == 2 || user->VelocitySolver == 3)
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Multigrid methods can not be used with Powell-Hesteness solver");
		PetscOptionsInsertString("-ph_kappa 1e4");
		PetscOptionsInsertString("-ph_max_it 50");
		PetscOptionsInsertString("-ph_tol 1e-8");
	}
	else if(user->StokesSolver == 2)
	{
		// 2 - Schur Complement Reduction
		PetscPrintf(PETSC_COMM_WORLD," Stokes solver             : Schur Complement Reduction\n");
		PetscOptionsInsertString("-scr_ksp_type fgmres");
		PetscOptionsInsertString("-scr_ksp_rtol 1.0e-6");
		PetscOptionsInsertString("-scr_ksp_max_it 1000");
		PetscOptionsInsertString("-scr_ksp_monitor");
		PetscOptionsInsertString("-scr_ksp_converged_reason");
		PetscOptionsInsertString("-scr_pc_type sor");
	}
	else if(user->StokesSolver == 3)
	{
		// 3 - Fully Coupled Solver
		PetscPrintf(PETSC_COMM_WORLD," Stokes solver             : Fully Coupled\n");
		PetscOptionsInsertString("-fc_ksp_type gcr");
		PetscOptionsInsertString("-fc_ksp_max_it 1000");
		PetscOptionsInsertString("-fc_ksp_min_it 4");
		PetscOptionsInsertString("-fc_ksp_converged_reason");

//==============================================================================================================
// *** this shit is necessary, because setting same options via function interface does not fucking work !!! ***
		PetscOptionsInsertString("-fc_pc_type fieldsplit");
		PetscOptionsInsertString("-fc_pc_fieldsplit_type SCHUR");
		PetscOptionsInsertString("-fc_pc_fieldsplit_schur_factorization_type UPPER");
//==============================================================================================================

		PetscOptionsInsertString("-use_stokes_residual");
		PetscOptionsInsertString("-use_stokes_monitor");
		PetscOptionsInsertString("-fc_ksp_v_atol 1.0e-16");
		PetscOptionsInsertString("-fc_ksp_p_atol 1.0e-16");
		PetscOptionsInsertString("-fc_ksp_v_rtol 1.0e-6");
		PetscOptionsInsertString("-fc_ksp_p_rtol 1.0e-6");
		PetscOptionsInsertString("-fc_ksp_v_min_atol 1.0e-5");
		PetscOptionsInsertString("-fc_ksp_p_min_atol 1.0e-5");
		PetscOptionsInsertString("-use_stokes_norm Linf");
		// pressure solver parameters
		PetscOptionsInsertString("-ps_ksp_type preonly");
		PetscOptionsInsertString("-ps_pc_type jacobi");

		// Add debugging options if requested
		if (SolverDiagnostics){
			// add some additional info about solvers to output
			PetscOptionsInsertString("-fc_ksp_view");
		}

	}
	else if (user->StokesSolver == 4)
	{
		// 4 - Only performing Matrix-Vector products (for scaling tests)
		PetscPrintf(PETSC_COMM_WORLD," Stokes solver             : NONE - Performing Matrix-Vector products instead \n");
	}

	//============================================
	// Set default options for the velocity solver
	//============================================
	if(user->VelocitySolver == 0)
	{
		// 0 - User defined
		PetscPrintf(PETSC_COMM_WORLD," Velocity solver           : User-defined (set options with prefix vs_)\n");
		PetscPrintf(PETSC_COMM_WORLD," Note: Galerkin geometric multigrid is deactivated for user-defined solver \n");
	}
	else if(user->VelocitySolver == 1)
	{
		// 1 - Direct (MUMPS in case we run on >1 procs)
		PetscOptionsInsertString("-vs_ksp_type preonly");
		PetscOptionsInsertString("-vs_pc_type lu");
		if(nproc != 1)
		{	PetscPrintf(PETSC_COMM_WORLD," Velocity solver           : Direct (MUMPS)\n");
			// MUMPS in case we run in parallel.
			// If you want to use a different parallel direct solver (e.g. superlu_dist), you should specify that on the command line.
			PetscOptionsInsertString("-vs_pc_factor_mat_solver_package mumps");
		}
		else PetscPrintf(PETSC_COMM_WORLD," Velocity solver           : Direct (Build-in direct solver)\n");
	}
	else if(user->VelocitySolver == 2)
	{
		// 2 - Galerkin geometric multigrid
		PetscPrintf(PETSC_COMM_WORLD," Velocity solver           : Galerkin geometric multigrid\n");
		// check discretization restrictions
		if(__ELEMENT_TYPE__ == ELEMENT_FDSTAG)
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Galerkin geometric multigrid in not implemented for FDSTAG discretization");
		// Krylov options
		PetscOptionsInsertString("-vs_ksp_type fgmres");
		// choose proper accuracies for different Stokes solver types
		if(user->StokesSolver == 2)
		{	PetscOptionsInsertString("-vs_ksp_rtol 1.0e-6");
			PetscOptionsInsertString("-vs_ksp_max_it 1000");
		}
		else if(user->StokesSolver == 3)
		{	PetscOptionsInsertString("-vs_ksp_rtol 1.0e-3");
			PetscOptionsInsertString("-vs_ksp_max_it 1000");
		}
		// multigrid options
		PetscOptionsInsertString("-vs_pc_type mg");
		PetscOptionsInsertString("-vs_pc_mg_levels 2");
		PetscOptionsInsertString("-vs_pc_mg_galerkin");
		PetscOptionsInsertString("-vs_pc_mg_type multiplicative");
		PetscOptionsInsertString("-vs_pc_mg_cycle_type v");
		// level smoothers
		// note: default solver for fieldspit is preonly + bjacoby + ilu
		PetscOptionsInsertString("-vs_mg_levels_ksp_type gmres");
		PetscOptionsInsertString("-vs_mg_levels_ksp_max_it 4");
		PetscOptionsInsertString("-vs_mg_levels_pc_type fieldsplit");
		PetscOptionsInsertString("-vs_mg_levels_pc_fieldsplit_type ADDITIVE");
		PetscOptionsInsertString("-vs_mg_levels_pc_fieldsplit_block_size 3");
		// coarse solver
		PetscOptionsInsertString("-vs_mg_coarse_ksp_type preonly");
		if(nproc != 1)
		{	PetscOptionsInsertString("-vs_mg_coarse_pc_type redundant");
			PetscOptionsInsertString("-vs_mg_coarse_pc_redundant_number 2");
			PetscOptionsInsertString("-vs_mg_coarse_redundant_pc_factor_mat_solver_package mumps");
		}
		else PetscOptionsInsertString("-vs_mg_coarse_pc_type lu");
	}
	else if(user->VelocitySolver == 3)
	{

		// 3 - Fieldsplit + Algebraic Multigrid (GAMG)
		PetscPrintf(PETSC_COMM_WORLD," Velocity solver           : Fieldsplit + Algebraic Multigrid (ML)\n");
		// Krylov options
		PetscOptionsInsertString("-vs_ksp_type gcr");
		// choose proper accuracies for different Stokes solver types
		if(user->StokesSolver == 2)			// schur complement
		{	PetscOptionsInsertString("-vs_ksp_rtol 1.0e-6");
			PetscOptionsInsertString("-vs_ksp_max_it 500");
		}
		else if(user->StokesSolver == 3)	// fully coupled
		{	PetscOptionsInsertString("-vs_ksp_rtol 1.0e-3");
			PetscOptionsInsertString("-vs_ksp_max_it 1000");
		}
		// fieldsplit preconditioner
		PetscOptionsInsertString("-vs_pc_type fieldsplit");
		PetscOptionsInsertString("-vs_pc_fieldsplit_block_size 3");
		PetscOptionsInsertString("-vs_pc_fieldsplit_type ADDITIVE");

		// algebraic multigrid solver with default options for each velocity block
		PetscOptionsInsertString("-vs_fieldsplit_ksp_type preonly");
		PetscOptionsInsertString("-vs_fieldsplit_pc_type ml");

		// Add debugging options if requested
		if (SolverDiagnostics){
			// add some additional info about solvers to output
			PetscOptionsInsertString("-vs_ksp_converged_reason");
			PetscOptionsInsertString("-vs_fieldsplit_ksp_converged_reason");
		}

	}
	else if(user->VelocitySolver == 4)
		{

			// 4 - GCR + Algebraic Multigrid on full velocity block (HYPRE)
			PetscPrintf(PETSC_COMM_WORLD," Velocity solver           : GCR + Algebraic Multigrid Preconditioner on full velocity block (HYPRE) \n");
			// Krylov options
			PetscOptionsInsertString("-vs_ksp_type gcr");
			// choose proper accuracies for different Stokes solver types
			if(user->StokesSolver == 2)			// schur complement
			{	PetscOptionsInsertString("-vs_ksp_rtol 1.0e-3");
				PetscOptionsInsertString("-vs_ksp_max_it 1000");
			}
			else if(user->StokesSolver == 3)	// fully coupled
			{	PetscOptionsInsertString("-vs_ksp_rtol 1.0e-3");
				PetscOptionsInsertString("-vs_ksp_max_it 1000");
			}

			// AMG preconditioner on full velocity block
			PetscOptionsInsertString("-vs_pc_type hypre");


			// Add debugging options if requested
			if (SolverDiagnostics){
				// add some additional info about solvers to output
				PetscOptionsInsertString("-vs_ksp_converged_reason");
			}

		}
	else if(user->VelocitySolver == 5)
		{

			// 5 - Fieldsplit + Algebraic Multigrid (GAMG)
			PetscPrintf(PETSC_COMM_WORLD," Velocity solver           : Fieldsplit + Algebraic Multigrid (GAMG)\n");
			// Krylov options
			PetscOptionsInsertString("-vs_ksp_type gcr");
			// choose proper accuracies for different Stokes solver types
			if(user->StokesSolver == 2)			// schur complement
			{	PetscOptionsInsertString("-vs_ksp_rtol 1.0e-6");
				PetscOptionsInsertString("-vs_ksp_max_it 500");
			}
			else if(user->StokesSolver == 3)	// fully coupled
			{	PetscOptionsInsertString("-vs_ksp_rtol 1.0e-3");
				PetscOptionsInsertString("-vs_ksp_max_it 1000");
			}
			// fieldsplit preconditioner
			PetscOptionsInsertString("-vs_pc_type fieldsplit");
			PetscOptionsInsertString("-vs_pc_fieldsplit_block_size 3");
			PetscOptionsInsertString("-vs_pc_fieldsplit_type ADDITIVE");

			// algebraic multigrid solver with default options for each velocity block
			PetscOptionsInsertString("-vs_fieldsplit_ksp_type preonly");
			//PetscOptionsInsertString("-vs_fieldsplit_pc_type ml");


			// Use identical parameters as ML, but with GAMG [options are from an entry in the usergroup by Mark Adam,
			//   plus from comparing the solvers that were set by running LaMEM with ML and GAMG]:
			PetscOptionsInsertString("-vs_fieldsplit_pc_type gamg");
			PetscOptionsInsertString("-vs_fieldsplit_pc_gamg_type agg");
			PetscOptionsInsertString("-vs_fieldsplit_pc_gamg_agg_nsmooths 1");
			//PetscOptionsInsertString("	");
			PetscOptionsInsertString("-vs_fieldsplit_pc_gamg_threshold .05");
			PetscOptionsInsertString("-vs_fieldsplit_pc_gamg_coarse_eq_limit 50");
			PetscOptionsInsertString("-vs_fieldsplit_mg_coarse_pc_type redundant");
			PetscOptionsInsertString("-vs_fieldsplit_mg_levels_ksp_type richardson");
			PetscOptionsInsertString(" -vs_fieldsplit_mg_levels_pc_type sor");


			// Add debugging options if requested
			if (SolverDiagnostics){
				// add some additional info about solvers to output
				PetscOptionsInsertString("-vs_ksp_converged_reason");
				PetscOptionsInsertString("-vs_fieldsplit_ksp_converged_reason");
			}

	}

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/
/* Compute global properties that are interesting to record */
#undef __FUNCT__
#define __FUNCT__ "ComputeGlobalProperties"
PetscErrorCode ComputeGlobalProperties( DM da, UserContext *user, PetscInt itime, Vec Velocity, LaMEMVelPressureDA C)
{
	PetscMPIInt                 size;
	PetscErrorCode 				ierr;
	GlobalTimeDependentData		TimeDepData;
	PetscBool					flg=PETSC_FALSE;
	PetscInt					ix, iy, iz, xs, ys, zs, xm, ym, zm, iphase, ipart, phase;
	PetscInt		 			xmp,ymp,zmp,xsp,ysp,zsp, iel_x, iel_y, iel_z, intp, istress;
	PetscScalar					Vx_min_loc, Vx_max_loc, Vy_min_loc, Vy_max_loc, Vz_min_loc, Vz_max_loc;
	PetscScalar					Vx_min,     Vx_max,     Vy_min,     Vy_max,     Vz_min,     Vz_max;
	PetscScalar					Vrms_loc, Vrms;
	PetscScalar					MinXCoordPhase_Local[max_num_phases], 	MaxXCoordPhase_Local[max_num_phases];
	PetscScalar					MinYCoordPhase_Local[max_num_phases], 	MaxYCoordPhase_Local[max_num_phases];
	PetscScalar					MinZCoordPhase_Local[max_num_phases], 	MaxZCoordPhase_Local[max_num_phases];
	PetscScalar					MinXCoordPhase[max_num_phases], 		MaxXCoordPhase[max_num_phases];
	PetscScalar					MinYCoordPhase[max_num_phases], 		MaxYCoordPhase[max_num_phases];
	PetscScalar					MinZCoordPhase[max_num_phases], 		MaxZCoordPhase[max_num_phases];
	PetscScalar					DevStress_Local[6], 					DevStress[6];
	PetscScalar					DevStrainrate_Local[6], 				DevStrainrate[6];
	PetscScalar					Volume_Local,							Volume;
	Field 						***velocity;
	Vec							local_Vel;
	Particles					ParticleLocal;
	//MaterialsElement			***materials;
	PetscScalar 				***materials_array;
	MaterialsElementDynamic 	material_data;


	/* Perform some initializations etc, to retrieve velocity etc @ every cpu */
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);



	/* CPU has a portion of velocity */
	ierr = DMGetLocalVector(da,&local_Vel); CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da, Velocity, INSERT_VALUES, local_Vel); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(da,   Velocity, INSERT_VALUES, local_Vel); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, local_Vel,	&velocity);		CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);	CHKERRQ(ierr);

	Vx_min_loc = 0; Vx_max_loc = 0;
	Vy_min_loc = 0; Vy_max_loc = 0;
	Vz_min_loc = 0; Vz_max_loc = 0;
	Vrms_loc=0;
	for (iz=zs; iz<zs+zm; iz++){
		for (iy=ys; iy<ys+ym; iy++){
			for(ix=xs; ix<xs+xm; ix++){

				// Track max and min velocities on local CPU
				if (velocity[iz][iy][ix].Vx<Vx_min_loc){	Vx_min_loc = velocity[iz][iy][ix].Vx;	}
				if (velocity[iz][iy][ix].Vy<Vy_min_loc){	Vy_min_loc = velocity[iz][iy][ix].Vy;	}
				if (velocity[iz][iy][ix].Vz<Vz_min_loc){	Vz_min_loc = velocity[iz][iy][ix].Vz;	}
				if (velocity[iz][iy][ix].Vx>Vx_max_loc){	Vx_max_loc = velocity[iz][iy][ix].Vx;	}
				if (velocity[iz][iy][ix].Vy>Vy_max_loc){	Vy_max_loc = velocity[iz][iy][ix].Vy;	}
				if (velocity[iz][iy][ix].Vz>Vz_max_loc){	Vz_max_loc = velocity[iz][iy][ix].Vz;	}

				Vrms_loc = Vrms_loc + PetscSqr(  velocity[iz][iy][ix].Vx*velocity[iz][iy][ix].Vx +
						velocity[iz][iy][ix].Vy*velocity[iz][iy][ix].Vy +
						velocity[iz][iy][ix].Vz*velocity[iz][iy][ix].Vz  );

			}
		}
	}
	ierr = DMDAVecRestoreArray(da,local_Vel,&velocity);	CHKERRQ(ierr);
	DMRestoreLocalVector(da,&local_Vel);

	/* Copy the repective maximum and minimum values to CPU zero */
	ierr = MPI_Allreduce(&Vx_min_loc,&Vx_min,1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr = MPI_Allreduce(&Vy_min_loc,&Vy_min,1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr = MPI_Allreduce(&Vz_min_loc,&Vz_min,1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr = MPI_Allreduce(&Vx_max_loc,&Vx_max,1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr = MPI_Allreduce(&Vy_max_loc,&Vy_max,1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr = MPI_Allreduce(&Vz_max_loc,&Vz_max,1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr = MPI_Allreduce(&Vrms_loc,  &Vrms  ,1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);

	Vrms               = 	Vrms/( (double) size);
	TimeDepData.Vrms   = 	Vrms*user->Characteristic.Velocity;

	/* Time */
	TimeDepData.Time = user->time*user->Characteristic.Time;	// time

	TimeDepData.Vx_max = 	Vx_max*user->Characteristic.Velocity;
	TimeDepData.Vx_min = 	Vx_min*user->Characteristic.Velocity;
	TimeDepData.Vy_max = 	Vy_max*user->Characteristic.Velocity;
	TimeDepData.Vy_min = 	Vy_min*user->Characteristic.Velocity;
	TimeDepData.Vz_max = 	Vz_max*user->Characteristic.Velocity;
	TimeDepData.Vz_min = 	Vz_min*user->Characteristic.Velocity;


	/* Find the maximum and minimum coordinate of each particle within a specific phase */
	for (iphase=0; iphase<max_num_phases; iphase++){
		MinXCoordPhase_Local[iphase] 	= 	user->x_left  + 2*user->W;
		MaxXCoordPhase_Local[iphase] 	= 	user->x_left;
		MinYCoordPhase_Local[iphase] 	= 	user->y_front + 2*user->L;
		MaxYCoordPhase_Local[iphase] 	= 	user->y_front;
		MinZCoordPhase_Local[iphase] 	= 	user->z_bot   + 2*user->H;
		MaxZCoordPhase_Local[iphase] 	= 	user->z_bot;
	}
#ifdef PARTICLES
	/* Loop over all particles */
	for (ipart=0; ipart<user->num_particle_local; ipart++){
		ParticleLocal		=	user->ParticlesLocal[ipart];
		phase 				= 	((PetscInt) ParticleLocal.phase);
		if (ParticleLocal.z<MinZCoordPhase_Local[phase]){
			MinZCoordPhase_Local[phase] = ParticleLocal.z;
		}
		if (ParticleLocal.z>MaxZCoordPhase_Local[phase]){
			MaxZCoordPhase_Local[phase] = ParticleLocal.z;
		}
		if (ParticleLocal.x<MinXCoordPhase_Local[phase]){
			MinXCoordPhase_Local[phase] = ParticleLocal.x;
		}
		if (ParticleLocal.x>MaxXCoordPhase_Local[phase]){
			MaxXCoordPhase_Local[phase] = ParticleLocal.x;
		}
		if (ParticleLocal.y<MinYCoordPhase_Local[phase]){
			MinYCoordPhase_Local[phase] = ParticleLocal.y;
		}
		if (ParticleLocal.y>MaxYCoordPhase_Local[phase]){
			MaxYCoordPhase_Local[phase] = ParticleLocal.y;
		}
	}
	/* Collect global min. and maximum */
	ierr = MPI_Allreduce(&MinXCoordPhase_Local,&MinXCoordPhase,max_num_phases, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr = MPI_Allreduce(&MaxXCoordPhase_Local,&MaxXCoordPhase,max_num_phases, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr = MPI_Allreduce(&MinYCoordPhase_Local,&MinYCoordPhase,max_num_phases, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr = MPI_Allreduce(&MaxYCoordPhase_Local,&MaxYCoordPhase,max_num_phases, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr = MPI_Allreduce(&MinZCoordPhase_Local,&MinZCoordPhase,max_num_phases, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr = MPI_Allreduce(&MaxZCoordPhase_Local,&MaxZCoordPhase,max_num_phases, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);

#endif
	for (iphase=0; iphase<max_num_phases; iphase++){
		TimeDepData.MinXCoordPhase[iphase] = MinXCoordPhase[iphase]*user->Characteristic.Length;
		TimeDepData.MaxXCoordPhase[iphase] = MaxXCoordPhase[iphase]*user->Characteristic.Length;
		TimeDepData.MinYCoordPhase[iphase] = MinYCoordPhase[iphase]*user->Characteristic.Length;
		TimeDepData.MaxYCoordPhase[iphase] = MaxYCoordPhase[iphase]*user->Characteristic.Length;
		TimeDepData.MinZCoordPhase[iphase] = MinZCoordPhase[iphase]*user->Characteristic.Length;
		TimeDepData.MaxZCoordPhase[iphase] = MaxZCoordPhase[iphase]*user->Characteristic.Length;
	}

	/* Compute maximum & minimum topography */
	PetscScalar meanVec;
	PetscInt	length;
	VecMin(user->SurfaceTopography,PETSC_NULL, &TimeDepData.MinTopography);
	VecMax(user->SurfaceTopography,PETSC_NULL, &TimeDepData.MaxTopography);
	VecNorm(user->SurfaceTopography,NORM_1, &meanVec);
	VecGetSize(user->SurfaceTopography, &length);
	meanVec = meanVec/((PetscScalar) length);

	TimeDepData.MinTopography 	= TimeDepData.MinTopography *user->Characteristic.Length;
	TimeDepData.MaxTopography 	= TimeDepData.MaxTopography *user->Characteristic.Length;
	TimeDepData.MeanTopography 	= meanVec					*user->Characteristic.Length;
	//PetscPrintf(PETSC_COMM_WORLD,"Topography: min=%f, max=%f \n",TimeDepData.MinTopography,TimeDepData.MaxTopography);


	/* volumetrically averaged stresses & strain rates */
	ierr = DMDAVecGetArray(user->DA_Materials, user->Materials, &materials_array); CHKERRQ(ierr);
	ierr = DMDAGetCorners(user->DA_Processors,&xsp,&ysp,&zsp,&xmp,&ymp,&zmp);CHKERRQ(ierr);

	// Loop over all elements
	for (istress=0; istress<6; istress++){
		DevStress_Local[istress] = 0.0;
		DevStrainrate_Local[istress] = 0.0;
	}
	Volume_Local = 0.0;

	for (iel_z=zsp; iel_z<zsp+zmp; iel_z++){
		for (iel_y=ysp; iel_y<ysp+ymp; iel_y++){
			for(iel_x=xsp; iel_x<xsp+xmp; iel_x++){
				LaMEMSetMaterialDataMemoryFromArray( &material_data, iel_x-xsp,iel_y-ysp,iel_z-zsp, C->ngp_vel, materials_array );

				/* Loop over integration points */
				for (intp=0; intp<C->ngp_vel; intp++){
					for (istress=0; istress<6; istress++){

						// Sum stresses
						DevStress_Local[istress] 	= DevStress_Local[istress] + material_data.DevStress[istress][intp]*material_data.ElementVolume[intp]/((double) C->ngp_vel);

						// Sum strainrates
						DevStrainrate_Local[istress] = DevStrainrate_Local[istress] + material_data.DevStrainrate[istress][intp]*material_data.ElementVolume[intp]/((double) C->ngp_vel);

					}

					// Volume
					Volume_Local 			= 	Volume_Local + material_data.ElementVolume[intp];
				}

			}
		}
	}
	ierr = DMDAVecRestoreArray(user->DA_Materials, user->Materials, &materials_array); CHKERRQ(ierr);

	// collect stress data from all CPU's
	ierr = MPI_Allreduce(&DevStress_Local,		&DevStress,		6, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr = MPI_Allreduce(&DevStrainrate_Local,	&DevStrainrate,	6, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr = MPI_Allreduce(&Volume_Local,   &Volume,   1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);

	// Create volumetric average and store data
	for (istress=0; istress<6; istress++){
		DevStress[istress] 		= DevStress[istress]/Volume;
		DevStrainrate[istress] 	= DevStrainrate[istress]/Volume;


		TimeDepData.DevStress[istress] 		= DevStress[istress]*user->Characteristic.Stress;
		TimeDepData.DevStrainrate[istress] 	= DevStrainrate[istress]*user->Characteristic.Strainrate;
	}

	/* Print stuff to screen if required */
	flg = PETSC_FALSE;	// initialize
	PetscOptionsGetBool( PETSC_NULL, "-DisplayAverageQuantities", &flg, PETSC_NULL );
	if (flg==PETSC_TRUE){
		PetscPrintf(PETSC_COMM_WORLD,"# Volume-averaged quantities: \n");
		istress = 0; PetscPrintf(PETSC_COMM_WORLD,"# Txx = %g Exx = %g Volume = %g Mu_xx=%g \n",TimeDepData.DevStress[istress], TimeDepData.DevStrainrate[istress], Volume, TimeDepData.DevStress[istress]/TimeDepData.DevStrainrate[istress]/2.0);
		istress = 1; PetscPrintf(PETSC_COMM_WORLD,"# Tyy = %g Eyy = %g Volume = %g Mu_yy=%g \n",TimeDepData.DevStress[istress], TimeDepData.DevStrainrate[istress], Volume, TimeDepData.DevStress[istress]/TimeDepData.DevStrainrate[istress]/2.0);
		istress = 2; PetscPrintf(PETSC_COMM_WORLD,"# Tzz = %g Ezz = %g Volume = %g Mu_zz=%g \n",TimeDepData.DevStress[istress], TimeDepData.DevStrainrate[istress], Volume, TimeDepData.DevStress[istress]/TimeDepData.DevStrainrate[istress]/2.0);
		istress = 3; PetscPrintf(PETSC_COMM_WORLD,"# Txz = %g Exz = %g Volume = %g Mu_xz=%g \n",TimeDepData.DevStress[istress], TimeDepData.DevStrainrate[istress], Volume, TimeDepData.DevStress[istress]/TimeDepData.DevStrainrate[istress]/2.0);
		istress = 4; PetscPrintf(PETSC_COMM_WORLD,"# Txy = %g Exy = %g Volume = %g Mu_xy=%g \n",TimeDepData.DevStress[istress], TimeDepData.DevStrainrate[istress], Volume, TimeDepData.DevStress[istress]/TimeDepData.DevStrainrate[istress]/2.0);
		istress = 5; PetscPrintf(PETSC_COMM_WORLD,"# Tyz = %g Eyz = %g Volume = %g Mu_yz=%g \n",TimeDepData.DevStress[istress], TimeDepData.DevStrainrate[istress], Volume, TimeDepData.DevStress[istress]/TimeDepData.DevStrainrate[istress]/2.0);
	}

	/* Store array */
	user->TimeDependentData[itime] = TimeDepData;


	/* Print some info to the screen (useful for monitoring some simulations) */
	PetscPrintf(PETSC_COMM_WORLD,"Minimum z-level phase 1 %g \n",TimeDepData.MinZCoordPhase[1]);


	PetscFunctionReturn(0);
}
/*==========================================================================================================*/

/*==========================================================================================================*/

#undef __FUNCT__
#define __FUNCT__ "DASetGhostedCoordinates"
/*@
   DASetGhostedCoordinates - Sets into the DMDA a vector that indicates the coordinates of the local nodes (including ghost nodes).

   Not collective

   Input Parameter:
.  da - the distributed array
.  c - coordinate vector (including the ghost nodes)

   Note:
   Each process has the coordinates for its local AND ghost nodes

   For two and three dimensions coordinates are interlaced (x_0,y_0,x_1,y_1,...)
   and (x_0,y_0,z_0,x_1,y_1,z_1...)

   The vector, c should be obtained from either
     DMGetCoordinatesLocal(da,&gc);
     VecDuplicate(gc,&c);
     VecDestroy(gc);
   or
     DMCreateLocalVector(da,&c);

   The user is responsible for destroying the vector c.

   The ghosted coordinates are inserted into the global vector using INSERT_VALUES.
   Hence, determining which ghost values will appear in the global vector is determined
   by the order in which the messages are recieved. It's the case of last in, best dressed...

  Level: intermediate

.keywords: distributed array, get, corners, nodes, local indices, coordinates

.seealso: DMGetCoordinatesLocal(), DASetCoordinates(), DASetUniformCoordinates(), DMGetCoordinates(), DMGetCoordinateDM()
@*/
PetscErrorCode DASetGhostedCoordinates(DM da,Vec c)
{
	PetscErrorCode ierr;
	Vec global,gc;
	DM cda;


	ierr = DMGetCoordinateDM(da,&cda); CHKERRQ(ierr);

	ierr = DMGetCoordinates(da,&global); CHKERRQ(ierr);

	ierr = DMGetCoordinatesLocal(da,&gc); CHKERRQ(ierr);
	ierr = VecCopy(c,gc); CHKERRQ(ierr);
	//	VecDestroy(gc);  // THIS STEP IS CRUCIAL (not sure why) //
	ierr = DMLocalToGlobalBegin(cda,gc,INSERT_VALUES,global); CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd  (cda,gc,INSERT_VALUES,global); CHKERRQ(ierr);


	ierr = DMGlobalToLocalBegin(cda,global,INSERT_VALUES,gc); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(cda,global,INSERT_VALUES,gc); CHKERRQ(ierr);

	//ierr = VecDestroy(gc); CHKERRQ(ierr);
	//ierr = VecDestroy(global); CHKERRQ(ierr);
	//ierr = DMDestroy(cda); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}




/*==========================================================================================================*/
/*




 * Writes the stiffness matrix to disk
 */
#undef __FUNCT__
#define __FUNCT__ "WriteStiffnessMatrixToDisk"
PetscErrorCode WriteStiffnessMatrixToDisk( DM da, Mat A11, Mat A12, Mat A21, Mat A22, Mat approx_S, Vec f, Vec h, Vec sol_vel, Vec Pressure, UserContext *user)
{
	PetscErrorCode	ierr;
	PetscViewer		viewer;
	PetscBool		DumpStiffnessMatrixes;
	PetscInt 		nnode_x, nnode_y, nnode_z;
	Vec 			InfoVec;

	// Check whether a command-line option is present to dump the stiffness matrixes
	DumpStiffnessMatrixes = PETSC_FALSE;
	ierr = PetscOptionsHasName(PETSC_NULL,"-DumpStiffnessMatrixes",&DumpStiffnessMatrixes); CHKERRQ(ierr);

	if(DumpStiffnessMatrixes == PETSC_TRUE)
	{
		// Create vector with info
		ierr = DMDAGetInfo(da, 0, &nnode_x, &nnode_y, &nnode_z, 0,0,0,0,0,0,0,0,0); CHKERRQ(ierr); // # of nodes in all directions

		ierr = VecCreate(PETSC_COMM_WORLD, &InfoVec); CHKERRQ(ierr);
		ierr = VecSetSizes(InfoVec,PETSC_DECIDE,9); CHKERRQ(ierr);
		ierr = VecSetFromOptions(InfoVec); CHKERRQ(ierr);

		ierr = VecSetValue(InfoVec,0,((double) user->nel_x  ), INSERT_VALUES); CHKERRQ(ierr);
		ierr = VecSetValue(InfoVec,1,((double) user->nel_y  ), INSERT_VALUES); CHKERRQ(ierr);
		ierr = VecSetValue(InfoVec,2,((double) user->nel_z  ), INSERT_VALUES); CHKERRQ(ierr);
		ierr = VecSetValue(InfoVec,3,((double) nnode_x), INSERT_VALUES); CHKERRQ(ierr);
		ierr = VecSetValue(InfoVec,4,((double) nnode_y), INSERT_VALUES); CHKERRQ(ierr);
		ierr = VecSetValue(InfoVec,5,((double) nnode_z), INSERT_VALUES); CHKERRQ(ierr);

		ierr = VecSetValue(InfoVec,6,user->W, INSERT_VALUES); CHKERRQ(ierr);
		ierr = VecSetValue(InfoVec,7,user->L, INSERT_VALUES); CHKERRQ(ierr);
		ierr = VecSetValue(InfoVec,8,user->H, INSERT_VALUES); CHKERRQ(ierr);

		ierr = VecAssemblyBegin(InfoVec); 	CHKERRQ(ierr);
		ierr = VecAssemblyEnd(InfoVec); 	CHKERRQ(ierr);

		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"stiffness.dat",FILE_MODE_WRITE,&viewer); CHKERRQ(ierr);
		ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_BINARY_MATLAB); CHKERRQ(ierr);

		ierr = MatView(A11,viewer); CHKERRQ(ierr);
		ierr = MatView(A12,viewer); CHKERRQ(ierr);
		ierr = MatView(A21,viewer); CHKERRQ(ierr);
		ierr = MatView(A22,viewer); CHKERRQ(ierr);
		ierr = MatView(approx_S,viewer); CHKERRQ(ierr);
		ierr = VecView(InfoVec,viewer); CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

		ierr = VecDestroy(&InfoVec); CHKERRQ(ierr);

		PetscPrintf(PETSC_COMM_WORLD,"*************  Wrote A11, A12, A21, A22 & approx_Schur stiffness matrixes in stiffness.dat *************\n");

		if(f != PETSC_NULL || h != PETSC_NULL || sol_vel != PETSC_NULL || Pressure != PETSC_NULL)
		{
			ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"rhs_vector.dat",FILE_MODE_WRITE,&viewer); CHKERRQ(ierr);
			ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_BINARY_MATLAB); CHKERRQ(ierr);

			if(f != PETSC_NULL)        { ierr = VecView(f, viewer); CHKERRQ(ierr); }

			if(h != PETSC_NULL)        { ierr = VecView(h, viewer); CHKERRQ(ierr); }

			if(sol_vel != PETSC_NULL)  { ierr = VecView(sol_vel, viewer); CHKERRQ(ierr); }

			if(Pressure != PETSC_NULL) { ierr = VecView(Pressure, viewer); CHKERRQ(ierr); }

			ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
			PetscPrintf(PETSC_COMM_WORLD,"*************  Wrote f,h & sol_vel rhs-vectors & V/P solution vectors in rhs_vector.dat ********************************************\n");
			PetscPrintf(PETSC_COMM_WORLD,"*************  Read the files into MATLAB with Read_StiffnessMatrix_ForDebugging.m *****************\n");
		}

		// write A11,A12, A21, A22, f & h into a single binary file that can be read with

		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"StiffnessMatrixesRhs_Binary.dat",FILE_MODE_WRITE,&viewer); CHKERRQ(ierr);

		ierr = MatView(A11,viewer); CHKERRQ(ierr);
		ierr = MatView(A12,viewer); CHKERRQ(ierr);
		ierr = MatView(A21,viewer); CHKERRQ(ierr);
		ierr = MatView(A22,viewer); CHKERRQ(ierr);
		if(f != PETSC_NULL) { ierr = VecView(f,viewer); CHKERRQ(ierr); }
		if(h != PETSC_NULL) { ierr = VecView(h,viewer); CHKERRQ(ierr); }

		ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"*************  Wrote A11,A12, A21, A22, b1 & b2 into StiffnessMatrixesRhs_Binary.dat using PETSC binary format ************************************\n");

		MPI_Barrier(PETSC_COMM_WORLD);
	}

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/



/*==========================================================================================================*/
/* WriteStiffnessMatrixToDisk
 * Writes the stiffness matrix to disk
 */
#undef __FUNCT__
#define __FUNCT__ "WriteTemperatureStiffnessMatrixToDisk"
PetscErrorCode WriteTemperatureStiffnessMatrixToDisk( DM da_temp, Mat T_MAT, Vec rhs_temp, Vec Temp, UserContext *user)
{
	PetscErrorCode	ierr;
	PetscViewer		viewer;
	PetscBool		DumpStiffnessMatrixes;
	PetscInt 		nnode_x, nnode_y, nnode_z;
	Vec 			InfoVec;

	/* Check whether a command-line option is present to dump the stiffness matrixes */
	DumpStiffnessMatrixes = PETSC_FALSE;
	ierr = PetscOptionsHasName(PETSC_NULL,"-DumpStiffnessMatrixes",&DumpStiffnessMatrixes); CHKERRQ(ierr);

	if (DumpStiffnessMatrixes==PETSC_TRUE){
		// Create vector with info
		ierr = DMDAGetInfo(da_temp,		0,&nnode_x,  &nnode_y,	&nnode_z,	0,0,0,0,0,0,0,0,0);	CHKERRQ(ierr);  // # of nodes 	 in all directions

		//ierr = VecCreateSeq(PETSC_COMM_WORLD, 9, &InfoVec); CHKERRQ(ierr);
		ierr = VecCreate(PETSC_COMM_WORLD, &InfoVec); CHKERRQ(ierr);
		ierr = VecSetSizes(InfoVec,PETSC_DECIDE,9); CHKERRQ(ierr);
		ierr = VecSetFromOptions(InfoVec); CHKERRQ(ierr);

		ierr = VecSetValue(InfoVec,0,((double) user->nel_x  ), INSERT_VALUES); CHKERRQ(ierr);
		ierr = VecSetValue(InfoVec,1,((double) user->nel_y  ), INSERT_VALUES); CHKERRQ(ierr);
		ierr = VecSetValue(InfoVec,2,((double) user->nel_z  ), INSERT_VALUES); CHKERRQ(ierr);
		ierr = VecSetValue(InfoVec,3,((double) nnode_x), INSERT_VALUES); CHKERRQ(ierr);
		ierr = VecSetValue(InfoVec,4,((double) nnode_y), INSERT_VALUES); CHKERRQ(ierr);
		ierr = VecSetValue(InfoVec,5,((double) nnode_z), INSERT_VALUES); CHKERRQ(ierr);

		ierr = VecSetValue(InfoVec,6,user->W, INSERT_VALUES); CHKERRQ(ierr);
		ierr = VecSetValue(InfoVec,7,user->L, INSERT_VALUES); CHKERRQ(ierr);
		ierr = VecSetValue(InfoVec,8,user->H, INSERT_VALUES); CHKERRQ(ierr);


		ierr = VecAssemblyBegin(InfoVec); 	CHKERRQ(ierr);
		ierr = VecAssemblyEnd(InfoVec); 	CHKERRQ(ierr);


//		ierr = PetscViewerBinaryMatlabOpen(PETSC_COMM_WORLD,"stiffness_T.dat",&viewer); CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"stiffness_T.dat",FILE_MODE_WRITE,&viewer); CHKERRQ(ierr);
		ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_BINARY_MATLAB); CHKERRQ(ierr);

		ierr = MatView(T_MAT,viewer); CHKERRQ(ierr);
		ierr = VecView(InfoVec,viewer); CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

		ierr = VecDestroy(&InfoVec); CHKERRQ(ierr);

		PetscPrintf(PETSC_COMM_WORLD,"*************  Wrote T_MAT stiffness matrix in stiffness_T.dat *************\n");

		if(rhs_temp!=PETSC_NULL){
//			ierr = PetscViewerBinaryMatlabOpen(PETSC_COMM_WORLD,"rhs_vector_T.dat",&viewer); CHKERRQ(ierr);
			ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"rhs_vector_T.dat",FILE_MODE_WRITE,&viewer); CHKERRQ(ierr);
			ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_BINARY_MATLAB); CHKERRQ(ierr);

			if(rhs_temp != PETSC_NULL){
				ierr = VecView(rhs_temp,viewer); CHKERRQ(ierr);
			}
			ierr = VecView(Temp ,viewer); CHKERRQ(ierr);

			ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
			PetscPrintf(PETSC_COMM_WORLD,"*************  Wrote rhs_temp & T rhs-vectors & temperature solution vectors in rhs_vector_T.dat ********************************************\n");
			PetscPrintf(PETSC_COMM_WORLD,"*************  Read the files into MATLAB with Read_StiffnessMatrix_ForDebugging.m *****************\n");


		}

		MPI_Barrier(PETSC_COMM_WORLD);
	}

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/*
Creates an M x N sized array.
The entries in the array can be accessed via two mechanisms;
 a) A2[i][j]  , 0 <= i < M,  i <= j < N
 b) A[k]      , 0 <= k < M * N
Either A2 or A can be PETSC_NULL;
 */
#undef __FUNCT__
#define __FUNCT__ "LaMEMCreate2dArray"
PetscErrorCode LaMEMCreate2dArray( const PetscInt M, const PetscInt N, PetscScalar ***_A2, PetscScalar **_A )
{
	PetscErrorCode ierr;
	PetscInt i;
	PetscScalar *A, **A2;

	if( (_A==PETSC_NULL) && (_A2==PETSC_NULL) ) {
		PetscPrintf( PETSC_COMM_WORLD, "WARNING(%s): Both args cannot be PETSC_NULL\n", __func__ );
		PetscFunctionReturn(0);
	}

	/* create flat array */
	ierr = PetscMalloc( sizeof(PetscScalar)*(size_t)(M*N), &A ); CHKERRQ(ierr);
	ierr = PetscMemzero( A, sizeof(PetscScalar)*(size_t)(M*N) ); CHKERRQ(ierr);

	/* create multidimensional representation if needed */
	if( _A2 != PETSC_NULL ) {
		ierr = PetscMalloc( sizeof(PetscScalar**)*(size_t)M, &A2 ); CHKERRQ(ierr);
		for( i=0; i<M; i++ ) {
			A2[i] = &A[i*N];
		}
		/* set the pointer */
		*_A2 = A2;
	}

	if( _A != PETSC_NULL ) {
		/* set the pointer */
		*_A  = A;
	}

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/

/*==========================================================================================================*/
#undef __FUNCT__
#define __FUNCT__ "LaMEMDestroy2dArray"
PetscErrorCode LaMEMDestroy2dArray(PetscScalar ***_A2, PetscScalar **_A )
{
	PetscErrorCode ierr;
	PetscScalar *A, **A2;
	PetscScalar *flat;

	/* get the flat memory space and free it */
	if( _A2 != PETSC_NULL ) {
		A2 = *_A2;

		flat = &A2[0][0];

		/* free flat array */
		ierr = PetscFree( flat ); CHKERRQ(ierr);
	}
	else if( _A != PETSC_NULL ) {
		A = *_A;

		flat = A;

		/* free flat array */
		ierr = PetscFree( flat ); CHKERRQ(ierr);
	}
	else {
		PetscPrintf( PETSC_COMM_WORLD, "WARNING(%s): Both args cannot be PETSC_NULL\n", __func__ );
		PetscFunctionReturn(0);
	}


	/* free multidimensional representation if it exists */
	if( _A2 != PETSC_NULL ) {
		ierr = PetscFree( A2 ); CHKERRQ(ierr);
		*_A2 = PETSC_NULL;
	}
	if( _A != PETSC_NULL ) {
		*_A  = PETSC_NULL;
	}

	PetscFunctionReturn(0);
}

/*==========================================================================================================*/
#undef __FUNCT__
#define __FUNCT__ "GetMaterialPropertiesFromCommandLine"
PetscErrorCode GetMaterialPropertiesFromCommandLine(UserContext *user)
{
	PetscErrorCode 	ierr;
	PetscBool		flg,get_options;
	PetscInt 		i;
	char 			matprop_opt[PETSC_MAX_PATH_LEN];


	PetscFunctionBegin;

	flg = PETSC_FALSE;
	get_options = PETSC_FALSE;

	ierr = PetscOptionsGetBool( PETSC_NULL, "-SetMaterialProperties", &get_options, PETSC_NULL ); 					CHKERRQ(ierr);

	if(get_options) {
		PetscPrintf(PETSC_COMM_WORLD,"# --------------------------------------------------------------------------\n");
		PetscPrintf(PETSC_COMM_WORLD,"# Material properties set from command line: \n");

		for(i=0;i<user->num_phases;i++){

			// viscosity
			sprintf(matprop_opt,"-mu_%lld",(LLD)i);
			ierr = PetscOptionsGetReal(PETSC_NULL ,matprop_opt,&user->PhaseProperties.mu[i]	, &flg); 				CHKERRQ(ierr);
			if(flg == PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD,"#    mu[%lld]	= %g \n",(LLD)i,user->PhaseProperties.mu[i]);

			// density
			sprintf(matprop_opt,"-rho_%lld",(LLD)i);
			ierr = PetscOptionsGetReal(PETSC_NULL ,matprop_opt,&user->PhaseProperties.rho[i]	, &flg);			CHKERRQ(ierr);
			if(flg == PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD,"#    rho[%lld]	= %g \n",(LLD)i,user->PhaseProperties.rho[i]);

			// feel free to implement others ...
			//user->PhaseProperties.mu[i]
			//user->PhaseProperties.FrankKamenetskii[i]
			//user->PhaseProperties.Powerlaw_e0[i]
			//user->PhaseProperties.A[i]
			//user->PhaseProperties.E[i]
			//
			//user->PhaseProperties.ElasticShearModule[i]
			//user->PhaseProperties.Cohesion[i]
			//user->PhaseProperties.CohesionAfterWeakening[i]
			//
			//user->PhaseProperties.T_Conductivity[i]
			//user->PhaseProperties.HeatCapacity[i]
			//user->PhaseProperties.RadioactiveHeat[i]
			//
			//user->PhaseProperties.rho[i]
			//user->PhaseProperties.Density_T0[i]
			//user->PhaseProperties.ThermalExpansivity[i]

		}
		PetscPrintf(PETSC_COMM_WORLD,"# --------------------------------------------------------------------------\n");
	}


	PetscFunctionReturn(0);

}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "GetLithColumnFromCommandLine"
PetscErrorCode GetLithColumnFromCommandLine(UserContext *user)
{
	PetscErrorCode 	ierr;
	PetscBool		flg,get_options;
	PetscInt 		i;
	char 			matprop_opt[PETSC_MAX_PATH_LEN];


	PetscFunctionBegin;

	flg = PETSC_FALSE;
	get_options = PETSC_FALSE;

	ierr = PetscOptionsGetBool( PETSC_NULL, "-SetLithColumn", &get_options, PETSC_NULL ); 					CHKERRQ(ierr);

	if(get_options) {



		PetscPrintf(PETSC_COMM_WORLD,"# --------------------------------------------------------------------------\n");
		PetscPrintf(PETSC_COMM_WORLD,"# Background density set from command line: \n");

		for(i=0;i<=user->GravityField.LithColNum;i++){

			// depth
			if(i<user->GravityField.LithColNum)
			{
				sprintf(matprop_opt,"-lithdepth_%lld",(LLD)i);
				ierr = PetscOptionsGetReal(PETSC_NULL ,matprop_opt,&user->GravityField.LithColDepth[i]	, &flg); 				CHKERRQ(ierr);
				if(flg == PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD,"#    depth[%lld]  = %g \n",(LLD)i,user->GravityField.LithColDepth[i]);
			}
			// density
			sprintf(matprop_opt,"-lithrho_%lld",(LLD)i);
			ierr = PetscOptionsGetReal(PETSC_NULL ,matprop_opt,&user->GravityField.LithColDens[i]	, &flg);			CHKERRQ(ierr);
			if(flg == PETSC_TRUE) PetscPrintf(PETSC_COMM_WORLD,"#    rho[%lld]  = %g \n",(LLD)i,user->GravityField.LithColDens[i]);


		}
		PetscPrintf(PETSC_COMM_WORLD,"# --------------------------------------------------------------------------\n");
	}


	PetscFunctionReturn(0);

}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DMDAGetProcessorRank"
PetscErrorCode DMDAGetProcessorRank(DM da, PetscInt *rank_x, PetscInt *rank_y, PetscInt *rank_z, PetscInt *rank_col)
{
	PetscMPIInt	rank;
	PetscInt	m, n, p, i, j, k, colind;
	PetscFunctionBegin;
	// get MPI processor rank
	MPI_Comm_rank(((PetscObject)da)->comm, &rank);
	// get number processors in each coordinate direction
	DMDAGetInfo(da, 0, 0, 0, 0, &m, &n, &p, 0, 0, 0, 0, 0, 0);
	// determine i-j-k coordinates of processor
	// x-index runs first (i), then y (j), followed by z (k)
	getLocalRank(&i, &j, &k, rank, m, n);
	// compute index of x-y column of processors
	// (same rule as above for x and y coordinates)
	colind = i + j*m;
	// assign output
	if(rank_x)   (*rank_x)   = i;
	if(rank_y)   (*rank_y)   = j;
	if(rank_z)   (*rank_z)   = k;
	if(rank_col) (*rank_col) = colind;
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DMDAView"
PetscErrorCode DMDAView(const char * name, DM da)
{
	PetscMPIInt     rank, size;
	PetscErrorCode  ierr;
	PetscInt        x, y, z, m, n, p;
	PetscInt        gx, gy, gz, gm, gn, gp;
	PetscFunctionBegin;
	// print basic DMDA local information on every processor
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank); CHKERRQ(ierr);
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);
	DMDAGetCorners     (da, &x,  &y,  &z,  &m,  &n,  &p);
	DMDAGetGhostCorners(da, &gx, &gy, &gz, &gm, &gn, &gp);
	ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,
		"name=%s; rank=%2lld; x=%2lld; y=%2lld; z=%2lld; m=%2lld; n=%2lld; p=%2lld; gx=%2lld; gy=%2lld; gz=%2lld; gm=%2lld; gn=%2lld; gp=%2lld\n",
		name, (LLD)rank, (LLD)x, (LLD)y, (LLD)z, (LLD)m, (LLD)n, (LLD)p, (LLD)gx, (LLD)gy, (LLD)gz, (LLD)gm, (LLD)gn, (LLD)gp); CHKERRQ(ierr);
	ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DMDAViewVTK"
PetscErrorCode DMDAViewVTK(const char * filename, DM da)
{
	PetscViewer	    viewer;
	PetscErrorCode  ierr;
	PetscFunctionBegin;
	ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF, filename, &viewer); CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_VTK);     CHKERRQ(ierr);
	ierr = DMView(da, viewer);                                       CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&viewer );                             CHKERRQ(ierr);
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "makeIntArray"
PetscErrorCode makeIntArray(PetscInt **arr, const PetscInt *init, const PetscInt n)
{
	PetscInt       *tmp;
	size_t          sz;
	PetscErrorCode 	ierr;
	PetscFunctionBegin;
	// compute size in bytes
	sz = (size_t)n*sizeof(PetscInt);
	// allocate space
	ierr = PetscMalloc(sz, &tmp); CHKERRQ(ierr);
	// initialize memory from vector (if required)
	if(init) { ierr = PetscMemcpy(tmp, init, sz); CHKERRQ(ierr); }
	// or just clear memory
	else { ierr = PetscMemzero(tmp, sz); CHKERRQ(ierr); }
	// return pointer
	*arr = tmp;
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "makeScalArray"
PetscErrorCode makeScalArray(PetscScalar **arr, const PetscScalar *init, const PetscInt n)
{
	PetscScalar    *tmp;
	size_t          sz;
	PetscErrorCode 	ierr;
	PetscFunctionBegin;
	// compute size in bytes
	sz = (size_t)n*sizeof(PetscScalar);
	// allocate space
	ierr = PetscMalloc(sz, &tmp); CHKERRQ(ierr);
	// initialize memory from vector (if required)
	if(init) { ierr = PetscMemcpy(tmp, init, sz); CHKERRQ(ierr); }
	// or just clear memory
	else { ierr = PetscMemzero(tmp, sz); CHKERRQ(ierr); }
	// return pointer
	*arr = tmp;
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DebugSave2Bin_GlobalVec"
PetscErrorCode DebugSave2Bin_GlobalVec(Vec vec_data,const char *filename,PetscInt itime, const char *folder)
{
	PetscViewer 	view_out;
	PetscErrorCode 	ierr;
	char 			filename_in[PETSC_MAX_PATH_LEN];

	PetscFunctionBegin;

	ierr = LaMEMCreateOutputDirectory(folder); CHKERRQ(ierr);

	sprintf(filename_in,"./%s/Debug_%s_%lld.bin",folder,filename,(LLD)itime);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"--- DEBUG: Write %s \n",filename_in);				CHKERRQ(ierr);
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename_in, FILE_MODE_WRITE, &view_out);	CHKERRQ(ierr);
    ierr = VecView(vec_data,view_out);														CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&view_out);													CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "DebugSave2Bin_GlobalMat"
PetscErrorCode DebugSave2Bin_GlobalMat(Mat mat_data,const char *filename,PetscInt itime, const char *folder)
{
	PetscViewer 	view_out;
	PetscErrorCode 	ierr;
	char 			filename_in[PETSC_MAX_PATH_LEN];

	PetscFunctionBegin;

	ierr = LaMEMCreateOutputDirectory(folder); CHKERRQ(ierr);

	sprintf(filename_in,"./%s/Debug_%s_%lld.bin",folder, filename,(LLD)itime);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"--- DEBUG: Write %s \n",filename_in);				CHKERRQ(ierr);
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename_in, FILE_MODE_WRITE, &view_out);	CHKERRQ(ierr);
    ierr = MatView(mat_data,view_out);														CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&view_out);													CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

//=====================================================================================================================
#undef __FUNCT__
#define __FUNCT__ "DMDACreate2dFrom3d"
PetscErrorCode DMDACreate2dFrom3d(DM da,PetscInt gp,DMDADirection dir,DM *da_plane,MPI_Comm *PLANE_COMM)
{
/*
 * 	Extracts a DMDA2d from a DMDA3d depending on a valid grid index, which is part of the distributed plane
 *
 * 	INPUT:
 * 	- da:			DMDA3d
 *	- gp:			global grid point number in Z direction
 *	- dir:			coordinate direction
 *
 * 	OUTPUT:
 * 	- da_plane		DADM of the plane
 * 	- PLANE_COMM 	... its MPI communicator (is needed to work with da_plane and destroy petsc objects later on...)
 *
 * 	Contributed by Tobias Baumann (Mainz, Oct 2012)
 */
//	--- Declarations --------------------------------------------------------------------------------------------------
	PetscErrorCode		ierr;
	DMDALocalInfo 		linfo;
	PetscInt			cpu[3],cpu_x,cpu_y,mx,my;
	const PetscInt 		*lx,*ly,*lz,*lxp,*lyp;

//	--- Code ----------------------------------------------------------------------------------------------------------
	PetscFunctionBegin;

	// Get information of 3dDMDA
	ierr = DMDAGetLocalInfo(da,&linfo);CHKERRQ(ierr);
	ierr = DMDAGetInfo(da, PETSC_NULL, PETSC_NULL,PETSC_NULL, PETSC_NULL, &cpu[0],&cpu[1],&cpu[2],PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
	ierr = DMDAGetOwnershipRanges(da,&lx,&ly,&lz);CHKERRQ(ierr);

	// Get communicator only for this plane
	ierr = DMDAGetProcessorSubset(da,dir,gp,PLANE_COMM);CHKERRQ(ierr);
	// PLANE_COMM == MPI_COMM_NULL for all ranks that do not own gp!

	// Create 2d DMDA + global vector
	if(dir==DMDA_X){
		mx=linfo.my;	my=linfo.mz;
		cpu_x=cpu[1];	cpu_y=cpu[2];
		lxp=ly;			lyp=lz;
	}else if(dir==DMDA_Y){
		mx=linfo.mx;	my=linfo.mz;
		cpu_x=cpu[0];	cpu_y=cpu[2];
		lxp=lx;			lyp=lz;
	}else if(dir==DMDA_Z){
		mx=linfo.mx;	my=linfo.my;
		cpu_x=cpu[0];	cpu_y=cpu[1];
		lxp=lx;			lyp=ly;
	}else SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Invalid direction");

	if(*PLANE_COMM != MPI_COMM_NULL)
		ierr = DMDACreate2d(*PLANE_COMM,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,mx,my,cpu_x,cpu_y,1,1,lxp,lyp,da_plane);CHKERRQ(ierr);
	// arguments 2,3,4,9,10 should be specified too

	PetscFunctionReturn(0);
}
//=====================================================================================================================
#undef __FUNCT__
#define __FUNCT__ "DMExtractGlobalVec2dFromGlobalVec3d"
PetscErrorCode DMExtractGlobalVec2dFromGlobalVec3d(DM da,Vec gvec_da,PetscInt gp, DM da_plane,MPI_Comm PLANE_COMM,Vec *gvec_plane)
{
/*
 * 	INPUT:
 * 	- da:			DMDA3d
 *	- gvec_da:		... its associated global vector
 *	- dir:			coordinate direction
 *	- gp:			global grid point number in Z direction
 *	- da_plane:		DMDA of the plane
 *
 * 	OUTPUT:
 * 	- gvec_plane	The associated global vector of da_plane that has the data of gvec_da => needs to be destroyed
 *
 *	Contributed by Tobias Baumann (Mainz, Oct 2012)
 */
//	--- Declarations --------------------------------------------------------------------------------------------------
	PetscErrorCode		ierr;
	DMDALocalInfo 		linfo;
	PetscInt			ix,iy;
	PetscScalar			***array_da;
	PetscScalar			**array_plane;

//	--- Code ----------------------------------------------------------------------------------------------------------
	PetscFunctionBegin;

	// Get DMDA3d array
	ierr = DMDAVecGetArray(da,gvec_da,&array_da); CHKERRQ(ierr);

	if(PLANE_COMM != MPI_COMM_NULL){
		// Get DMDA2d linfo + vector
		ierr = DMDAGetLocalInfo(da_plane,&linfo);CHKERRQ(ierr);
		ierr = DMCreateGlobalVector(da_plane,gvec_plane);CHKERRQ(ierr); // -> gvec needs to be destroyed outside !!

		// Get DMDA2d array
		ierr = DMDAVecGetArray(da_plane,*gvec_plane,&array_plane);CHKERRQ(ierr);
		for(iy=linfo.ys; iy<linfo.ys+linfo.ym; iy++){
			for(ix=linfo.xs; ix<linfo.xs+linfo.xm; ix++){
					array_plane[iy][ix] = array_da[gp][iy][ix];
			}
		}
		// Restore DMDA2d array
		ierr = DMDAVecRestoreArray(da_plane,*gvec_plane,&array_plane);CHKERRQ(ierr);
	}

	// Restore DMDA3d array
	ierr = DMDAVecRestoreArray(da,gvec_da,&array_da); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//=====================================================================================================================
#undef __FUNCT__
#define __FUNCT__ "DMDADestroy2dFrom3d"
PetscErrorCode DMDADestroy2dFrom3d(DM da_plane,MPI_Comm PLANE_COMM)
{
/*
 * 	Destroys DMDA2d + its communicator that was created from DMDA3d
 * 	INPUT:
 * 	- da_plane		DADM2d
 * 	- PLANE_COMM	... its MPI communicator
 *
 * 	Contributed by Tobias Baumann (Mainz, Oct 2012)
 */
//	--- Declarations --------------------------------------------------------------------------------------------------
	PetscErrorCode		ierr;
//	--- Code ----------------------------------------------------------------------------------------------------------
	PetscFunctionBegin;

	if(PLANE_COMM!=MPI_COMM_NULL){
		ierr = DMDestroy(&da_plane);CHKERRQ(ierr);
		MPI_Comm_free(&PLANE_COMM);
	}
	PetscFunctionReturn(0);
}
//=====================================================================================================================
#undef __FUNCT__
#define __FUNCT__ "DMDestroyGlobalVec2dFromGlobalVec3d"
PetscErrorCode DMDestroyGlobalVec2dFromGlobalVec3d(Vec gvec_plane,MPI_Comm PLANE_COMM)
{
/*
 * 	Destroys DMDA2d-gvec that was created from DMDA3d
 * 	INPUT:
 * 	- gvec_plane	The associated global vector of da_plane
 * 	- PLANE_COMM	The associated MPI communicator
 *
 * 	Contributed by Tobias Baumann (Mainz, Oct 2012)
 */
//	--- Declarations --------------------------------------------------------------------------------------------------
	PetscErrorCode		ierr;
//	--- Code ----------------------------------------------------------------------------------------------------------
	PetscFunctionBegin;

	if(PLANE_COMM!=MPI_COMM_NULL){
		ierr = VecDestroy(&gvec_plane);CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}
//=====================================================================================================================
#undef __FUNCT__
#define __FUNCT__ "SaveProcessorPartitioning"
PetscErrorCode SaveProcessorPartitioning(UserContext *user)
{
/*
 * 	Saves the processor partitioning of the domain to file
 * 	INPUT:
 * 	- user struct
 *
 * 	Contributed by Tobias Baumann (Mainz, Jan 2013)
 */
//	--- Declarations --------------------------------------------------------------------------------------------------
	PetscErrorCode		ierr;
	DM                  cda;
	Vec                 gc;
	PetscScalar         *xc,*yc,*zc;
	DMDACoor3d       ***coors;
	PetscInt            rx,ry,rz,xs,xm,ys,ym,zs,zm,M,N,P,Nmsg,Nprocx,Nprocy,Nprocz,onaxis=0,flag;
	PetscScalar         indSE[3];
    MPI_Request         request;
    MPI_Status          status;
	int                fileID;
	char 				*fname;

	//	--- Code ----------------------------------------------------------------------------------------------------------
	PetscFunctionBegin;

    PetscPrintf(PETSC_COMM_WORLD,"# --- Save processor partitioning ---\n");

	// Get info
	ierr = DMDAGetInfo(user->DA_Vel,0,&M,&N,&P,&Nprocx,&Nprocy,&Nprocz,0,0,0,0,0,0);
	ierr = DMDAGetCorners(user->DA_Vel,&xs,&ys,&zs,&xm,&ym,&zm);  CHKERRQ(ierr);

	// Get coordinates of velocity DMDA
	ierr = DMGetCoordinates(user->DA_Vel , &gc);                CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(user->DA_Vel, &cda);               CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda, gc, &coors);                      CHKERRQ(ierr);

	if(!(xs==0 && ys==0 && zs==0)){

		// X
		if(ys==0 && zs==0){
			indSE[0]  = coors[0][0][xs   ].x;
			indSE[1]  = coors[0][0][xs+xm-1].x;
			onaxis    = 1;
		}
		// Y
		if(xs==0 && zs==0){
			indSE[0]  = coors[0][ys   ][0].y;
			indSE[1]  = coors[0][ys+ym-1][0].y;
			onaxis    = 1;
		}
		// Z
		if(xs==0 && ys==0){
			indSE[0]  = coors[zs   ][0][0].z;
			indSE[1]  = coors[zs+zm-1][0][0].z;
			onaxis    = 1;
		}

		if(onaxis){
			//send indSE to rank zero
			ierr = MPI_Isend(indSE,2,MPI_DOUBLE,0,21,PETSC_COMM_WORLD,&request); CHKERRQ(ierr);
			// wait until sending process has finished
			ierr = MPI_Wait(&request, &status); CHKERRQ(ierr);
		}

	}
	else{ // This is rank 0

		// Allocate memory for three vectors x,y,z
		ierr = PetscMalloc((size_t)(Nprocx+1)*sizeof(PetscScalar),&xc);   CHKERRQ(ierr);
		ierr = PetscMalloc((size_t)(Nprocy+1)*sizeof(PetscScalar),&yc);   CHKERRQ(ierr);
		ierr = PetscMalloc((size_t)(Nprocz+1)*sizeof(PetscScalar),&zc);   CHKERRQ(ierr);

		// coordinate of origin
		xc[0] = coors[0][0][0].x;
		yc[0] = coors[0][0][0].y;
		zc[0] = coors[0][0][0].z;

		// number of messages
		Nmsg = Nprocx-1 + Nprocy-1 + Nprocz-1;

		// receive coordinates
		while(Nmsg!=0){
			MPI_Iprobe(MPI_ANY_SOURCE, 21, PETSC_COMM_WORLD,&flag, &status);
			if(flag){
				ierr = MPI_Recv(indSE,2,MPI_DOUBLE,MPI_ANY_SOURCE,21,PETSC_COMM_WORLD,&status); CHKERRQ(ierr);
                // get local ranks of sending-processor
				getLocalRank(&rx,&ry,&rz,status.MPI_SOURCE,Nprocx,Nprocy);

                if(ry==0 && rz==0){ // X-axis
					xc[rx]=indSE[0];
					if(rx==Nprocx-1)
						xc[Nprocx]=indSE[1];
				}
				if(rx==0 && rz==0){ // Y-axis
					yc[ry]=indSE[0];
					if(ry==Nprocy-1)
						yc[Nprocy]=indSE[1];
				}
				if(rx==0 && ry==0){ // Z-axis
					zc[rz]=indSE[0];
					if(rz==Nprocz-1)
						zc[Nprocz]=indSE[1];
				}
				// reduce counter of number of msgs to receive
				Nmsg=Nmsg-1;
			}
		}

		// Save to disk
		asprintf(&fname,"ProcessorPartitioning_%lldcpu_%lld.%lld.%lld.bin",(LLD)(Nprocx*Nprocy*Nprocz),(LLD)Nprocx,(LLD)Nprocy,(LLD)Nprocz);
		PetscBinaryOpen(fname,FILE_MODE_WRITE,&fileID);
		PetscBinaryWrite(fileID,&Nprocx,1,PETSC_INT,PETSC_FALSE);
		PetscBinaryWrite(fileID,&Nprocy,1,PETSC_INT,PETSC_FALSE);
		PetscBinaryWrite(fileID,&Nprocz,1,PETSC_INT,PETSC_FALSE);
		PetscBinaryWrite(fileID,xc,Nprocx+1,PETSC_SCALAR,PETSC_FALSE);
		PetscBinaryWrite(fileID,yc,Nprocy+1,PETSC_SCALAR,PETSC_FALSE);
		PetscBinaryWrite(fileID,zc,Nprocz+1,PETSC_SCALAR,PETSC_FALSE);
		PetscBinaryWrite(fileID,&user->Characteristic.Length,1,PETSC_SCALAR,PETSC_FALSE);
		PetscBinaryClose(fileID);

		free(fname);

		// Free allocated memory
		ierr = PetscFree(xc);                                         CHKERRQ(ierr);
		ierr = PetscFree(yc);                                         CHKERRQ(ierr);
		ierr = PetscFree(zc);                                         CHKERRQ(ierr);
	}

	// close access
	ierr = DMDAVecRestoreArray(cda, gc, &coors);                  CHKERRQ(ierr);


	PetscFunctionReturn(0);
}
//==========================================================================================================
PetscInt ISRankZero(MPI_Comm comm)
{
	PetscMPIInt rank;

	MPI_Comm_rank(comm, &rank);

	return (rank == 0);
}
//==========================================================================================================
PetscInt ISParallel(MPI_Comm comm)
{
	PetscMPIInt size;

	MPI_Comm_size(comm, &size);

	return (size > 1);
}
//==========================================================================================================
// Creates an output directory
#undef __FUNCT__
#define __FUNCT__ "LaMEMCreateOutputDirectory"
PetscErrorCode LaMEMCreateOutputDirectory(const char *DirectoryName)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// generate a new directory on rank zero
	if(ISRankZero(PETSC_COMM_WORLD))
	{
		if(mkdir(DirectoryName, S_IRWXU))
		{
			PetscPrintf(PETSC_COMM_WORLD," Writing output to existing directory %s \n", DirectoryName);
		}
		else
		{
			PetscPrintf(PETSC_COMM_WORLD," Created output directory %s \n", DirectoryName);
		}
	}

	// all other ranks should wait
	ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//==========================================================================================================

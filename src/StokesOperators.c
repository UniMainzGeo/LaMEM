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

StokesOperators.c, contains the following subroutines:

ComputeKinmtx					-	Compute kinematic matrix on element level
ComputeMaterialMatrix			-	Compute material matrix (currently: isotropic incompressible viscous)
Multiply_Kin_Times_Mat_Times_Kin-	Performs multiplication
ComputeInversePP				-	Computes inverse of pressure mass matrix
VP_invPP_PV						-	Performs matrix multiplications
CreateStencilInGlobalStiffness	-	Creates local2global numbering
ComputeDivergence				-	Compute divergence in all elements
UpdateRHS_Div					-	Update RHS to take into account divergence (required by Uzawa)

$Id$
$Date$

This version is compatible with PETSc version 3.3-p2
The code cannot be distributed, used for commercial purposes or handed on,
without the explicit agreement of Boris Kaus.
*/


#include "LaMEM.h"
#include "StokesOperators.h"
#include "Utils.h"
#include "Elements.h"
#include "Quadrature.h"
#include "Assembly.h"

/*==========================================================================================================*/
/* Computes the kinematic matrix */
void ComputeKinmtx(const PetscInt nnel, const PetscInt edof, double KinMtx[6][MAX_edof], double **dhdPhys, double **dhdsVel, double **InvJacob )
{
	/* It is assumed that the velocity, solution vector has the following form:
	*   [Vx(0) Vy(0) Vz(0)  Vx(1) Vy(1) Vz(1)   ...]  where the order is exactly as expected
	*   for the shape function */
	PetscInt	i, j, k;			//loop control

	// Compute dhdPhys - derivative versus natural coordinates
	for( i = 0; i < 3; i++){ for( j = 0; j < nnel; j++){ dhdPhys[i][j] = 0; }}
	for( i = 0; i < 3; i++){ for( j = 0; j < nnel; j++){ for( k = 0; k < 3; k++){
				dhdPhys[i][j] +=  InvJacob[i][k]*dhdsVel[k][j];
	}}}
	// Use dhdPhys to compute factors in front of strainrate
	for( i = 0; i < 6; i++){ for( j = 0; j < edof; j++){ KinMtx[i][j] = 0; }} // initialize
	for( i = 0; i < nnel; i++){ KinMtx[0][3*i] 	= 	dhdPhys[0][i]; 	}// Exx component
	for( i = 0; i < nnel; i++){ KinMtx[1][3*i+1] = 	dhdPhys[1][i];  }// Eyy component
	for( i = 0; i < nnel; i++){ KinMtx[2][3*i+2] = 	dhdPhys[2][i];  }// Ezz component
	for( i = 0; i < nnel; i++){ KinMtx[3][3*i  ] = 	dhdPhys[2][i]; 	KinMtx[3][3*i+2] = 	dhdPhys[0][i]; }// Exz component
	for( i = 0; i < nnel; i++){ KinMtx[4][3*i  ] = 	dhdPhys[1][i];	KinMtx[4][3*i+1] = 	dhdPhys[0][i]; } // Exy component
	for( i = 0; i < nnel; i++){ KinMtx[5][3*i+1] = 	dhdPhys[2][i];	KinMtx[5][3*i+2] = 	dhdPhys[1][i]; } // Eyz component
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Material matrix */
void ComputeMaterialMatrix( const double mu, double MatMtx[6][6] )
{
	PetscInt i,j;

	for (i=0;i<6;i++){ for (j=0;j<6;j++){ MatMtx[i][j]=0.0;	}}
	MatMtx[0][0] = 2*mu;	MatMtx[1][1] = 2*mu;	MatMtx[2][2] = 2*mu;
	MatMtx[3][3] =   mu;	MatMtx[4][4] =   mu;	MatMtx[5][5] =   mu;
}

/*==========================================================================================================*/
/* performs VV_add = Kinmtx*Matmtx*Kinmtx */
void Multiply_Kin_Times_Mat_Times_Kin( const PetscInt edof, double VV_add[MAX_edof][MAX_edof], double Kinmtx[6][MAX_edof], double Matmtx[6][6] )
{
	PetscInt 	i, j, k;			//loop control
	double	MatKin[6][MAX_edof];	//Intermediate matrix

	/* Compute Intermediate matrix */
	for( i = 0; i < 6; i++){ for( j = 0; j < edof; j++){ MatKin[i][j] = 0; }}
	for( i = 0; i < 6; i++){ for( j = 0; j < edof; j++){for( k = 0; k < 6; k++){
				MatKin[i][j] +=  Matmtx[i][k]*Kinmtx[k][j];
	}}}
	/* Compute VV_add */
	for( i = 0; i < edof; i++){ for( j = 0; j < edof; j++){ VV_add[i][j] = 0; }}
	for( i = 0; i < edof; i++){
		for( j = 0; j < edof; j++){
			for( k = 0; k < 6; k++){
				VV_add[i][j] +=  Kinmtx[k][i]*MatKin[k][j];
			}
		}
	}
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Compute inverse of PP matrix */
PetscErrorCode ComputeInversePP( double PP[MAX_npres][MAX_npres], double InversePP[MAX_npres][MAX_npres] )
{
	if( __ELEMENT_TYPE__ == ELEMENT_Q1P0 ) {
		InversePP[0][0] = 1.0/PP[0][0];
		PetscFunctionReturn(0);
	}

	if( __ELEMENT_TYPE__ == ELEMENT_Q1Q1 ) {
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER,"DO NOT COMPUTE invPP IF USING Q1Q1 elements!!! \n");
		PetscFunctionReturn(0);
	}

	if( __ELEMENT_TYPE__ == ELEMENT_Q2P1 ){
		double a, b, c, d,e,f,g,h,i,j,k,l,m,n,o,p;

			/* extract some data */
		a=PP[0][0];	b=PP[0][1]; c=PP[0][2]; d=PP[0][3];	e=PP[1][0];	f=PP[1][1]; g=PP[1][2]; h=PP[1][3];
		i=PP[2][0];	j=PP[2][1]; k=PP[2][2]; l=PP[2][3];	m=PP[3][0];	n=PP[3][1]; o=PP[3][2]; p=PP[3][3];

		/* inverse of matrix from MAPLE - yes I know it looks uggly*/
		InversePP[0][0] = (-f*k*p+f*l*o+j*g*p-j*h*o-n*g*l+n*h*k)/(-a*f*k*p+a*f*l*o+a*j*g*p-a*j*h*o-a*n*g*l+a*n*h*k+e*b*k*p-e*b*l*o-e*j*c*p+e*j*d*o+e*n*c*l-e*n*d*k-i*b*g*p+i*b*h*o+i*f*c*p-i*f*d*o-i*n*c*h+i*n*d*g+m*b*g*l-m*b*h*k-m*f*c*l+m*f*d*k+m*j*c*h-m*j*d*g);
		InversePP[0][1] = (b*k*p-b*l*o-j*c*p+j*d*o+n*c*l-n*d*k)/(-a*f*k*p+a*f*l*o+a*j*g*p-a*j*h*o-a*n*g*l+a*n*h*k+e*b*k*p-e*b*l*o-e*j*c*p+e*j*d*o+e*n*c*l-e*n*d*k-i*b*g*p+i*b*h*o+i*f*c*p-i*f*d*o-i*n*c*h+i*n*d*g+m*b*g*l-m*b*h*k-m*f*c*l+m*f*d*k+m*j*c*h-m*j*d*g);
		InversePP[0][2] = (-b*g*p+b*h*o+f*c*p-f*d*o-n*c*h+n*d*g)/(-a*f*k*p+a*f*l*o+a*j*g*p-a*j*h*o-a*n*g*l+a*n*h*k+e*b*k*p-e*b*l*o-e*j*c*p+e*j*d*o+e*n*c*l-e*n*d*k-i*b*g*p+i*b*h*o+i*f*c*p-i*f*d*o-i*n*c*h+i*n*d*g+m*b*g*l-m*b*h*k-m*f*c*l+m*f*d*k+m*j*c*h-m*j*d*g);
		InversePP[0][3] = -(-b*g*l+b*h*k+f*c*l-f*d*k-j*c*h+j*d*g)/(-a*f*k*p+a*f*l*o+a*j*g*p-a*j*h*o-a*n*g*l+a*n*h*k+e*b*k*p-e*b*l*o-e*j*c*p+e*j*d*o+e*n*c*l-e*n*d*k-i*b*g*p+i*b*h*o+i*f*c*p-i*f*d*o-i*n*c*h+i*n*d*g+m*b*g*l-m*b*h*k-m*f*c*l+m*f*d*k+m*j*c*h-m*j*d*g);

		InversePP[1][0] = -(-e*k*p+e*l*o+i*g*p-i*h*o-m*g*l+m*h*k)/(-a*f*k*p+a*f*l*o+a*j*g*p-a*j*h*o-a*n*g*l+a*n*h*k+e*b*k*p-e*b*l*o-e*j*c*p+e*j*d*o+e*n*c*l-e*n*d*k-i*b*g*p+i*b*h*o+i*f*c*p-i*f*d*o-i*n*c*h+i*n*d*g+m*b*g*l-m*b*h*k-m*f*c*l+m*f*d*k+m*j*c*h-m*j*d*g);
		InversePP[1][1] = -(a*k*p-a*l*o-i*c*p+i*d*o+m*c*l-m*d*k)/(-a*f*k*p+a*f*l*o+a*j*g*p-a*j*h*o-a*n*g*l+a*n*h*k+e*b*k*p-e*b*l*o-e*j*c*p+e*j*d*o+e*n*c*l-e*n*d*k-i*b*g*p+i*b*h*o+i*f*c*p-i*f*d*o-i*n*c*h+i*n*d*g+m*b*g*l-m*b*h*k-m*f*c*l+m*f*d*k+m*j*c*h-m*j*d*g);
		InversePP[1][2] = -(-a*g*p+a*h*o+e*c*p-e*d*o-m*c*h+m*d*g)/(-a*f*k*p+a*f*l*o+a*j*g*p-a*j*h*o-a*n*g*l+a*n*h*k+e*b*k*p-e*b*l*o-e*j*c*p+e*j*d*o+e*n*c*l-e*n*d*k-i*b*g*p+i*b*h*o+i*f*c*p-i*f*d*o-i*n*c*h+i*n*d*g+m*b*g*l-m*b*h*k-m*f*c*l+m*f*d*k+m*j*c*h-m*j*d*g);
		InversePP[1][3] = (-a*g*l+a*h*k+e*c*l-e*d*k-i*c*h+i*d*g)/(-a*f*k*p+a*f*l*o+a*j*g*p-a*j*h*o-a*n*g*l+a*n*h*k+e*b*k*p-e*b*l*o-e*j*c*p+e*j*d*o+e*n*c*l-e*n*d*k-i*b*g*p+i*b*h*o+i*f*c*p-i*f*d*o-i*n*c*h+i*n*d*g+m*b*g*l-m*b*h*k-m*f*c*l+m*f*d*k+m*j*c*h-m*j*d*g);

		InversePP[2][0] = -(e*j*p-e*l*n-i*f*p+i*h*n+m*f*l-m*h*j)/(-a*f*k*p+a*f*l*o+a*j*g*p-a*j*h*o-a*n*g*l+a*n*h*k+e*b*k*p-e*b*l*o-e*j*c*p+e*j*d*o+e*n*c*l-e*n*d*k-i*b*g*p+i*b*h*o+i*f*c*p-i*f*d*o-i*n*c*h+i*n*d*g+m*b*g*l-m*b*h*k-m*f*c*l+m*f*d*k+m*j*c*h-m*j*d*g);
		InversePP[2][1] = (a*j*p-a*l*n-i*b*p+i*d*n+m*b*l-m*d*j)/(-a*f*k*p+a*f*l*o+a*j*g*p-a*j*h*o-a*n*g*l+a*n*h*k+e*b*k*p-e*b*l*o-e*j*c*p+e*j*d*o+e*n*c*l-e*n*d*k-i*b*g*p+i*b*h*o+i*f*c*p-i*f*d*o-i*n*c*h+i*n*d*g+m*b*g*l-m*b*h*k-m*f*c*l+m*f*d*k+m*j*c*h-m*j*d*g);
		InversePP[2][2] = -(a*f*p-a*h*n-e*b*p+e*d*n+m*b*h-m*d*f)/(-a*f*k*p+a*f*l*o+a*j*g*p-a*j*h*o-a*n*g*l+a*n*h*k+e*b*k*p-e*b*l*o-e*j*c*p+e*j*d*o+e*n*c*l-e*n*d*k-i*b*g*p+i*b*h*o+i*f*c*p-i*f*d*o-i*n*c*h+i*n*d*g+m*b*g*l-m*b*h*k-m*f*c*l+m*f*d*k+m*j*c*h-m*j*d*g);
		InversePP[2][3] = (a*f*l-a*h*j-e*b*l+e*d*j+i*b*h-i*d*f)/(-a*f*k*p+a*f*l*o+a*j*g*p-a*j*h*o-a*n*g*l+a*n*h*k+e*b*k*p-e*b*l*o-e*j*c*p+e*j*d*o+e*n*c*l-e*n*d*k-i*b*g*p+i*b*h*o+i*f*c*p-i*f*d*o-i*n*c*h+i*n*d*g+m*b*g*l-m*b*h*k-m*f*c*l+m*f*d*k+m*j*c*h-m*j*d*g);

		InversePP[3][0] = (e*j*o-e*k*n-i*f*o+i*g*n+m*f*k-m*g*j)/(-a*f*k*p+a*f*l*o+a*j*g*p-a*j*h*o-a*n*g*l+a*n*h*k+e*b*k*p-e*b*l*o-e*j*c*p+e*j*d*o+e*n*c*l-e*n*d*k-i*b*g*p+i*b*h*o+i*f*c*p-i*f*d*o-i*n*c*h+i*n*d*g+m*b*g*l-m*b*h*k-m*f*c*l+m*f*d*k+m*j*c*h-m*j*d*g);
		InversePP[3][1] = -(a*j*o-a*k*n-i*b*o+i*c*n+m*b*k-m*c*j)/(-a*f*k*p+a*f*l*o+a*j*g*p-a*j*h*o-a*n*g*l+a*n*h*k+e*b*k*p-e*b*l*o-e*j*c*p+e*j*d*o+e*n*c*l-e*n*d*k-i*b*g*p+i*b*h*o+i*f*c*p-i*f*d*o-i*n*c*h+i*n*d*g+m*b*g*l-m*b*h*k-m*f*c*l+m*f*d*k+m*j*c*h-m*j*d*g);
		InversePP[3][2] = (a*f*o-a*g*n-e*b*o+e*c*n+m*b*g-m*c*f)/(-a*f*k*p+a*f*l*o+a*j*g*p-a*j*h*o-a*n*g*l+a*n*h*k+e*b*k*p-e*b*l*o-e*j*c*p+e*j*d*o+e*n*c*l-e*n*d*k-i*b*g*p+i*b*h*o+i*f*c*p-i*f*d*o-i*n*c*h+i*n*d*g+m*b*g*l-m*b*h*k-m*f*c*l+m*f*d*k+m*j*c*h-m*j*d*g);
		InversePP[3][3] = -(a*f*k-a*g*j-e*b*k+e*c*j+i*b*g-i*c*f)/(-a*f*k*p+a*f*l*o+a*j*g*p-a*j*h*o-a*n*g*l+a*n*h*k+e*b*k*p-e*b*l*o-e*j*c*p+e*j*d*o+e*n*c*l-e*n*d*k-i*b*g*p+i*b*h*o+i*f*c*p-i*f*d*o-i*n*c*h+i*n*d*g+m*b*g*l-m*b*h*k-m*f*c*l+m*f*d*k+m*j*c*h-m*j*d*g);
	}

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Compute VP*inv(PP)*PV */
void VP_invPP_PV( const PetscInt edof, const PetscInt npres, double VV_correction[MAX_edof][MAX_edof], double VP[MAX_edof][MAX_npres], double PP[MAX_npres][MAX_npres], double PV[MAX_npres][MAX_edof] )
{
	PetscInt 	i, j, k;				//loop control
	double	InversePP[MAX_npres][MAX_npres], invPP_PV[MAX_npres][MAX_edof];			//Intermediate matrix

	/* Compute Intermediate matrix */
	ComputeInversePP(PP, InversePP);
	for( i = 0; i < npres; i++){ for( j = 0; j < edof; j++){ invPP_PV[i][j] = 0; }}
	for( i = 0; i < npres; i++){ for( j = 0; j < edof; j++){for( k = 0; k < npres; k++){
				invPP_PV[i][j] +=  InversePP[i][k]*PV[k][j];
	}}}

	/* Compute VV_add */
	for( i = 0; i < edof; i++){ for( j = 0; j < edof; j++){ VV_correction[i][j] = 0; }}
	for( i = 0; i < edof; i++){ for( j = 0; j < edof; j++){for( k = 0; k < npres; k++){
				VV_correction[i][j] +=  VP[i][k]*invPP_PV[k][j];
	}}}
}
/*==========================================================================================================*/
/* Computes divergence for every element */
#undef __FUNCT__
#define __FUNCT__ "ComputeDivergence"
PetscErrorCode ComputeDivergence( LaMEMVelPressureDA C, DM da, Vec Velocity, DM da_pres, Vec Divergence,
		PetscScalar *MaximumDivergence)
{
	PetscErrorCode	ierr;
	PetscInt		i,j,k, iel_x, iel_y, iel_z;
	PetscInt		 xmp,ymp,zmp,xsp,ysp,zsp;
	PetscInt		 ii,jj,kk, intp, num_elem;
	DMDACoor3d		 ***coords, coord_elem[MAX_nnel], CoordIntp;
	DM				cda;
	Vec				gc;
	Field			 ***velocity;
	PetscScalar		IntPoint[3][MAX_ngp_vel], IntWeight[MAX_ngp_vel], ShapeVel[MAX_nnel], **dhdsVel, Point[3], ShapeP[MAX_npres];
	PetscScalar		**Jacob, **InvJacob, DetJacob, KinMtx[6][MAX_edof], **dhdPhys;
	PetscScalar		PV[MAX_npres][MAX_edof];
	PetscScalar		V_element[MAX_edof], Div[MAX_npres];
	PetscScalar		maxDiv, MaximumDiv, maxVel;
	PressureElem	 ***div_local;
	Vec				local_Vel;
	DAVPElementType element_type;
	PetscInt		ngp_vel, nnel, nintp_1D, npres, edof;
	PetscScalar		***div_local_array;
	PressureElemDynamic div_local_data;


	PetscFunctionBegin;


	element_type = C->type;
	ngp_vel  = C->ngp_vel;
	nnel     = C->nnel;
	edof     = C->edof;
	npres    = C->npres;
	nintp_1D = C->nintp_1D;

	LaMEMCreate2dArray( 3, nnel, &dhdsVel, PETSC_NULL );
	LaMEMCreate2dArray( 3, nnel, &dhdPhys, PETSC_NULL );
	LaMEMCreate2dArray( 3, 3, &Jacob, PETSC_NULL );
	LaMEMCreate2dArray( 3, 3, &InvJacob, PETSC_NULL );



	VecSet( Divergence, 0.0 );

	//	DMDASetUniformCoordinates(da,user->x_left,user->x_left+user->W,user->y_front,user->y_front + user->L,user->z_bot,user->z_bot + user->H);
	DMGetCoordinateDM(da,&cda);	//coordinates
	DMGetCoordinatesLocal(da,&gc);
	DMDAVecGetArray(cda,gc,&coords);

	/* The divergence vector has the same size as the pressure vector */
	ierr = DMDAVecGetArray(da_pres, Divergence,	&div_local);		CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da_pres, Divergence,	&div_local_array);		CHKERRQ(ierr);


	// Copy global velocity solution to local processor, including ghostpoints
	DMGetLocalVector(da,&local_Vel);
	DMGlobalToLocalBegin(da, Velocity, INSERT_VALUES, local_Vel);
	DMGlobalToLocalEnd(da,   Velocity, INSERT_VALUES, local_Vel);
	ierr = DMDAVecGetArray(da, local_Vel,	&velocity);		CHKERRQ(ierr);


	num_elem = 0;
	maxDiv   = 0.0;
	ierr = DMDAGetCorners(da,&xsp,&ysp,&zsp,&xmp,&ymp,&zmp);CHKERRQ(ierr);
	i=0; j=0; k=0;
	for (iel_z=zsp; iel_z<zsp+zmp; iel_z++){
		for (iel_y=ysp; iel_y<ysp+ymp; iel_y++){
			for(iel_x=xsp; iel_x<xsp+xmp; iel_x++){
				LaMEMSetPressueElemDataMemoryFromArray( &div_local_data, iel_x,iel_y,iel_z, npres,div_local_array );

				if( element_type==DAVP_Q1P0 ) {
					i = iel_x;		j = iel_y;		k = iel_z;
				}
				else if( (element_type==DAVP_Q2PM1L) || (element_type==DAVP_Q2PM1G) ){
					i = 2*iel_x;  	j = 2*iel_y; 	k = 2*iel_z;
				}
				else {
					SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Unknown element type" );
				}

				//------------------------------------------------------
				// Form the local stiffness matrix
				// -----------------------------------------------------
				// Extract coordinates of the local element in correct order
				GetElementCoords(coord_elem, coords, i,j,k,1);

				// Stuf that can be precomputed
				IntegrationPoints( nnel, nintp_1D, IntPoint, IntWeight);

				// Initialization
				for(ii=0; ii<npres;ii++){ for(jj=0; jj<edof;  jj++){	PV[ii][jj]=0.0;	}}
				for(ii=0; ii<npres; ii++){	Div[ii]=0;		}

				// Loop over integration points
				for (intp=0; intp<ngp_vel; intp++){

					Point[0]= IntPoint[0][intp];Point[1]= IntPoint[1][intp];Point[2]= IntPoint[2][intp];
					ComputeShapeFunctionVelocity(ShapeVel, dhdsVel, Point);				 // Velocity shape function
					ComputeCoordIntp(nnel,ShapeVel, coord_elem,  &CoordIntp);				 // Real coordinate of integration point
					ComputeShapeFunctionPressure(ShapeP, 	CoordIntp, Point);					 // Pressure shape function
					ComputeJacobianElement(nnel,dhdsVel,coord_elem,Jacob,InvJacob,&DetJacob);  // Jacobian etc. of element
					ComputeKinmtx( nnel, edof, KinMtx, dhdPhys, dhdsVel, InvJacob);					 // Kinematic matrix

					/* Catch errors (for debugging only) */
					if (DetJacob<0){
						PetscPrintf(PETSC_COMM_WORLD,"In routine UpdateRHSDiv: Negative Jacobian on DetJacob=%g  \n",DetJacob);
						PetscPrintf(PETSC_COMM_WORLD,"  element: [%lld,%lld,%lld]  \n",(LLD)iel_x,(LLD)iel_y,(LLD)iel_z);
						PetscPrintf(PETSC_COMM_WORLD," x-coords of local element: [%g,%g,%g,%g,%g,%g,%g,%g] \n",coord_elem[0].x,coord_elem[1].x,coord_elem[2].x,coord_elem[3].x,coord_elem[4].x,coord_elem[5].x,coord_elem[6].x,coord_elem[7].x);
						PetscPrintf(PETSC_COMM_WORLD," y-coords of local element: [%g,%g,%g,%g,%g,%g,%g,%g] \n",coord_elem[0].y,coord_elem[1].y,coord_elem[2].y,coord_elem[3].y,coord_elem[4].y,coord_elem[5].y,coord_elem[6].y,coord_elem[7].y);
						PetscPrintf(PETSC_COMM_WORLD," z-coords of local element: [%g,%g,%g,%g,%g,%g,%g,%g] \n",coord_elem[0].z,coord_elem[1].z,coord_elem[2].z,coord_elem[3].z,coord_elem[4].z,coord_elem[5].z,coord_elem[6].z,coord_elem[7].z);

						MPI_Abort(PETSC_COMM_WORLD,1);
					}


					// Compute PV matrix
					for (ii=0; ii<nnel; ii++){for (jj=0; jj<npres; jj++){for (kk=0; kk<3; kk++){
								PV[jj][3*ii+kk]	=	PV[jj][3*ii+kk]	+	dhdPhys[kk][ii]*ShapeP[jj]*DetJacob*IntWeight[intp];
					}}}

				}
				//------------------------------------------------------
				// End of forming the local stiffness matrix
				// -----------------------------------------------------

				// Extract Velocity of the current element
				GetVelocityElement(velocity, V_element, i,j,k);

				for (ii=0; ii<edof; ii++){	if (PetscAbs(V_element[ii])>maxVel){
						maxVel = PetscAbs(V_element[ii]);
				}}


				// Compute divergence of the current element
				for (ii=0; ii<npres; ii++){	for (jj=0; jj<edof; jj++){
						Div[ii] = Div[ii] - PV[ii][jj]*V_element[jj];
				}}


				// Update and store Divergence (which is defined element-wise)
				for (ii=0; ii<npres; ii++){
					div_local_data.P[ii] 	= 	Div[ii];
				}

				// Track the maximum divergence
				for (ii=0; ii<npres; ii++){	if (PetscAbs(Div[ii])>maxDiv){
						maxDiv = PetscAbs(Div[ii]);
				}}

				num_elem = num_elem+1;

			}
		}
	}
	DMDAVecRestoreArray(da_pres,Divergence,	&div_local);
	ierr = DMDAVecRestoreArray(da_pres, Divergence,	&div_local_array);		CHKERRQ(ierr);
	//VecAssemblyBegin(Divergence);			VecAssemblyEnd(Divergence);


	ierr = DMDAVecRestoreArray(da,	local_Vel,&velocity);	CHKERRQ(ierr);
	DMRestoreLocalVector(da,		&local_Vel);
	ierr = DMDAVecRestoreArray(cda,gc,&coords);	CHKERRQ(ierr);


	if( MaximumDivergence != NULL ) {
		// Collect the maximum divergence of all processors
		ierr = MPI_Allreduce(&maxDiv,&MaximumDiv,1, MPIU_SCALAR,MPI_MAX,PETSC_COMM_WORLD); CHKERRQ(ierr);

		// Normalize the maximum divergence over velocity
		VecNorm(Velocity, NORM_INFINITY, &maxVel);

		*MaximumDivergence =   MaximumDiv/maxVel;
	}

	// Clean up
	//DMDestroy( cda );
	//VecDestroy(gc);

	LaMEMDestroy2dArray(&dhdsVel, PETSC_NULL );
	LaMEMDestroy2dArray(&dhdPhys, PETSC_NULL );
	LaMEMDestroy2dArray(&Jacob, PETSC_NULL );
	LaMEMDestroy2dArray(&InvJacob, PETSC_NULL );

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


void LaMEMSetPressueElemDataMemoryFromArray(
		PressureElemDynamic *data,
		const PetscInt i, const PetscInt j, const PetscInt k,
		const PetscInt block_size,
		PetscScalar ***list )
{
	PetscScalar *block_ijk;

	block_ijk = &list[k][j][i*block_size];

	data->P  = &block_ijk[0];
}

/*==========================================================================================================*/
/*
Computes:
  vel = G pressure
in a matrix-free manner.

The following routine can be computed more efficiently by pre-computing the matrixes during formation
of the stiffness matrix and by retrieving them when required
*/
#undef __FUNCT__
#define __FUNCT__ "ComputeGradient"
PetscErrorCode ComputeGradient( LaMEMVelPressureDA C, DM da_vel, DM da_pres, Vec pressure, Vec vel, Vec vel_local, UserContext *user)
{
	PetscErrorCode	ierr;
	PetscInt		i,j,k,ipres, iel_x, iel_y, iel_z;
	PetscInt		xmp,ymp,zmp,xsp,ysp,zsp;
	PetscInt		ii,jj,kk, intp, num_elem;
	DMDACoor3d		***coords, coord_elem[MAX_nnel], CoordIntp;
	DM				cda_vel;
	Vec				gc;
	Field			***rhs;
	PetscScalar		IntPoint[3][MAX_ngp_vel], IntWeight[MAX_ngp_vel], ShapeVel[MAX_nnel], **dhdsVel, Point[3], ShapeP[MAX_npres];
	PetscScalar		**Jacob, **InvJacob, DetJacob, KinMtx[6][MAX_edof], **dhdPhys;
	PetscScalar		VP[MAX_edof][MAX_npres];
	PetscScalar		element_pressure[MAX_npres];
	PetscScalar		dVRHS[MAX_edof];
	PressureElem	***pressure_local;
	DAVPElementType element_type;
	PetscInt ngp_vel, nnel, nintp_1D, npres, edof;
	PetscScalar ***pressure_local_array;
	PressureElemDynamic pressure_local_data;


	PetscFunctionBegin;


	element_type = C->type;
	ngp_vel  = C->ngp_vel;
	nnel     = C->nnel;
	edof     = C->edof;
	npres    = C->npres;
	nintp_1D = C->nintp_1D;


	LaMEMCreate2dArray( 3, nnel, &dhdsVel, PETSC_NULL );
	LaMEMCreate2dArray( 3, nnel, &dhdPhys, PETSC_NULL );
	LaMEMCreate2dArray( 3, 3, &Jacob, PETSC_NULL );
	LaMEMCreate2dArray( 3, 3, &InvJacob, PETSC_NULL );


	ierr = VecSet(vel,0.0); CHKERRQ(ierr);
	ierr = VecSet(vel_local,0.0); CHKERRQ(ierr);


	ierr = DMGetCoordinateDM(da_vel,&cda_vel); CHKERRQ(ierr);//coordinates
	ierr = DMGetCoordinatesLocal(da_vel,&gc); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda_vel,gc,&coords); CHKERRQ(ierr);

	ierr = DMDAVecGetArray(da_pres, pressure, &pressure_local ); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da_pres, pressure, &pressure_local_array ); CHKERRQ(ierr);


	/* Make a loop over all local elements, construct the element stiffness matrix,
	* correct for bounda_velry conditions and  */
	ierr = DMDAVecGetArray( da_vel, vel_local, &rhs ); CHKERRQ(ierr);

	num_elem = 0; i=0; j=0; k=0;
	ierr = DMDAGetCorners(user->DA_Processors,&xsp,&ysp,&zsp,&xmp,&ymp,&zmp);CHKERRQ(ierr);
	for (iel_z=zsp; iel_z<zsp+zmp; iel_z++){
		for (iel_y=ysp; iel_y<ysp+ymp; iel_y++){
			for(iel_x=xsp; iel_x<xsp+xmp; iel_x++){
				LaMEMSetPressueElemDataMemoryFromArray( &pressure_local_data, iel_x,iel_y,iel_z, npres, pressure_local_array );


				if( element_type==DAVP_Q1P0 ) {
					i = iel_x;		j = iel_y;		k = iel_z;
				}
				else if( (element_type==DAVP_Q2PM1L) || (element_type==DAVP_Q2PM1G) ){
					i = 2*iel_x;  	j = 2*iel_y; 	k = 2*iel_z;
				}
				else {
					SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Unknown element type" );
				}

				/*------------------------------------------------------
				* Form the local stiffness matrix
				* -----------------------------------------------------*/
				/* Extract coordinates of the local element in correct order */
				GetElementCoords(coord_elem, coords, i,j,k,1);

				/* Stuff that can be precomputed */
				IntegrationPoints( nnel, nintp_1D, IntPoint, IntWeight);

				/* Initialization */
				for(ii=0; ii<edof; ii++){ for(jj=0; jj<npres; jj++){	VP[ii][jj]=0.0;	}}
				for(ii=0; ii<npres; ii++){	element_pressure[ii]=0;		}
				for(ii=0; ii<edof; ii++) {	dVRHS[ii]=0;	}
//				for(ii=0; ii<npres; ii++){	dP[ii]=0;		}

				/* Loop over integration points */
				for (intp=0; intp<ngp_vel; intp++){

					Point[0]= IntPoint[0][intp];Point[1]= IntPoint[1][intp];Point[2]= IntPoint[2][intp];
					ComputeShapeFunctionVelocity(ShapeVel, dhdsVel, Point);				 // Velocity shape function
					ComputeCoordIntp(nnel,ShapeVel, coord_elem,  &CoordIntp);				 // Real coordinate of integration point
					ComputeShapeFunctionPressure(ShapeP, 	CoordIntp, Point);					 // Pressure shape function
					ComputeJacobianElement(nnel,dhdsVel,coord_elem,Jacob,InvJacob,&DetJacob);  // Jacobian etc. of element
					ComputeKinmtx( nnel, edof, KinMtx, dhdPhys, dhdsVel, InvJacob);					 // Kinematic matrix

					/* Catch errors (for debugging only) */
					if (DetJacob<0){
						PetscPrintf(PETSC_COMM_WORLD,"In routine UpdateRHSDiv: Negative Jacobian on DetJacob=%g  \n",DetJacob);
						PetscPrintf(PETSC_COMM_WORLD,"  element: [%lld,%lld,%lld]  \n",(LLD)iel_x,(LLD)iel_y,(LLD)iel_z);
						PetscPrintf(PETSC_COMM_WORLD," x-coords of local element: [%g,%g,%g,%g,%g,%g,%g,%g] \n",coord_elem[0].x,coord_elem[1].x,coord_elem[2].x,coord_elem[3].x,coord_elem[4].x,coord_elem[5].x,coord_elem[6].x,coord_elem[7].x);
						PetscPrintf(PETSC_COMM_WORLD," y-coords of local element: [%g,%g,%g,%g,%g,%g,%g,%g] \n",coord_elem[0].y,coord_elem[1].y,coord_elem[2].y,coord_elem[3].y,coord_elem[4].y,coord_elem[5].y,coord_elem[6].y,coord_elem[7].y);
						PetscPrintf(PETSC_COMM_WORLD," z-coords of local element: [%g,%g,%g,%g,%g,%g,%g,%g] \n",coord_elem[0].z,coord_elem[1].z,coord_elem[2].z,coord_elem[3].z,coord_elem[4].z,coord_elem[5].z,coord_elem[6].z,coord_elem[7].z);

						MPI_Abort(PETSC_COMM_WORLD,1);
					}

					// Compute VP matrix
					for (ii=0; ii<nnel;ii++){
						for (jj=0; jj<npres; jj++){
							for (kk=0; kk<3; kk++){
								VP[3*ii+kk][jj]	=	VP[3*ii+kk][jj]	-	dhdPhys[kk][ii]*ShapeP[jj]*DetJacob*IntWeight[intp];
					}}}

				}
				/*------------------------------------------------------
				* End of forming the VP matrix
				* -----------------------------------------------------*/


				/* Extract pressure of the current element */
				for (ipres=0; ipres<npres; ipres++){
					element_pressure[ipres]  = pressure_local_data.P[ipres];
				}

				/* Compute new rhs */
				for (ii=0; ii<edof; ii++){
					for (jj=0; jj<npres; jj++){
						dVRHS[ii] = dVRHS[ii] - VP[ii][jj]*element_pressure[jj];
					}
				}

				/* Add dVRHS to the rhs vector (not taking care about BCs yet) */
				SetValuesRHS(rhs, dVRHS, i, j, k);

				num_elem = num_elem+1;

			}
		}
	}


	ierr = DMDAVecRestoreArray(cda_vel,gc,&coords);					CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da_pres, pressure, &pressure_local);	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da_pres, pressure, &pressure_local_array ); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray( da_vel, vel_local, &rhs );			CHKERRQ(ierr);

	ierr = DMLocalToGlobalBegin( da_vel, vel_local, ADD_VALUES, vel ); CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd( da_vel, vel_local, ADD_VALUES, vel );   CHKERRQ(ierr);

	/* Clean up */
	//DMDestroy( cda_vel );
	//VecDestroy(gc);

	LaMEMDestroy2dArray(&dhdsVel, PETSC_NULL );
	LaMEMDestroy2dArray(&dhdPhys, PETSC_NULL );
	LaMEMDestroy2dArray(&Jacob, PETSC_NULL );
	LaMEMDestroy2dArray(&InvJacob, PETSC_NULL );

	/* Set Boundary conditions to RHS - set them to zero !!*/
	SetBoundaryConditionsRHS(da_vel, user, vel, 0);


	PetscFunctionReturn(0);
}
/*==========================================================================================================*/

struct _p_MatrixFreeData_Divergence {
	Vec x_local;
	DM da_vel, da_pres;
	PetscInt m,n, M,N;
	LaMEMVelPressureDA C;
};

// Citation from PETSc documentation for MatShellSetOperation:
// All user-provided functions (execept for MATOP_DESTROY) should have the same
// calling sequence, as the usual matrix interface routines.
// So we can safely remove, unused parameters in the following function:

// PetscErrorCode MatDestroy_MFOperatorDivergence( Mat A, Vec x, Vec y )
PetscErrorCode MatDestroy_MFOperatorDivergence( Mat A )
{
	PetscErrorCode ierr;
	struct _p_MatrixFreeData_Divergence *mf_data;


	ierr = MatShellGetContext( A, (void**)&mf_data ); CHKERRQ(ierr);

	VecDestroy( &mf_data->x_local );
	PetscFree( mf_data );

	mf_data = PETSC_NULL;

	PetscFunctionReturn(0);
}

// PetscErrorCode MatDestroy_MFOperatorDivergence( Mat A, Vec x, Vec y )
PetscErrorCode MatMult_MFOperatorDivergence( Mat A, Vec x, Vec y )
{
	PetscErrorCode ierr;
	struct _p_MatrixFreeData_Divergence *mf_data;


	ierr = MatShellGetContext( A, (void**)&mf_data ); CHKERRQ(ierr);

	ierr = VecSet(y,0.0); CHKERRQ(ierr);
	ComputeDivergence( mf_data->C, mf_data->da_vel, x, mf_data->da_pres, y, PETSC_NULL );

	PetscFunctionReturn(0);
}

PetscErrorCode MatCreate_MFOperatorDivergence( LaMEMVelPressureDA C, DM da_vel,Vec vel, DM da_pres,Vec pres, Mat *A )
{
	PetscErrorCode ierr;
	struct _p_MatrixFreeData_Divergence *mf_data;
	MPI_Comm comm;
	PetscInt m,n, M,N;


	ierr = PetscMalloc( sizeof(struct _p_MatrixFreeData_Divergence), &mf_data ); CHKERRQ(ierr);

	/* da */
	mf_data->da_vel = da_vel;
	mf_data->da_pres = da_pres;
	mf_data->C = C;

	/* sizes */
	VecGetLocalSize( pres, &m );
	VecGetSize( pres, &M );

	VecGetLocalSize( vel, &n );
	VecGetSize( vel, &N );

	mf_data->m = m;
	mf_data->n = n;

	mf_data->M = M;
	mf_data->N = N;

	/* stash */
	DMCreateLocalVector( da_vel, &mf_data->x_local );

	PetscObjectGetComm( (PetscObject)da_vel, &comm );

	/* create matrix */
	ierr = MatCreateShell( comm, m,n, M,N, (void*)mf_data, A ); CHKERRQ(ierr);

	/* set operations */
	MatShellSetOperation( *A, MATOP_MULT,    (void(*)(void))MatMult_MFOperatorDivergence );
	MatShellSetOperation( *A, MATOP_DESTROY, (void(*)(void))MatDestroy_MFOperatorDivergence );


	PetscFunctionReturn(0);
}

void test_MFOperatorDivergence( LaMEMVelPressureDA C, DM da, Vec Velocity, DM da_pres, Vec Divergence, PetscScalar *MaximumDivergence )
{
	PetscScalar norm;
	Mat A;


	ComputeDivergence( C, da, Velocity, da_pres, Divergence, MaximumDivergence );
	VecNorm( Divergence, NORM_2, &norm );
	PetscPrintf( PETSC_COMM_WORLD, "orig: ||p||_2 = %1.10e \n", norm );


	MatCreate_MFOperatorDivergence( C, da,Velocity,da_pres,Divergence, &A );
	MatMult( A, Velocity, Divergence );
	VecNorm( Divergence, NORM_2, &norm );
	PetscPrintf( PETSC_COMM_WORLD, "  MF: ||p||_2 = %1.10e \n", norm );

	MatDestroy(&A);
}

/*==========================================================================================================*/

struct _p_MatrixFreeData_Gradient {
	Vec y_local;
	DM da_vel, da_pres;
	PetscInt m,n, M,N;
	UserContext *user;
	LaMEMVelPressureDA C;
};

// Citation from PETSc documentation for MatShellSetOperation:
// All user-provided functions (execept for MATOP_DESTROY) should have the same
// calling sequence, as the usual matrix interface routines.
// So we can safely remove, unused parameters in the following function:

PetscErrorCode MatDestroy_MFOperatorGradient( Mat A )
{
	PetscErrorCode ierr;
	struct _p_MatrixFreeData_Gradient *mf_data;


	ierr = MatShellGetContext( A, (void**)&mf_data ); CHKERRQ(ierr);

	VecDestroy( &mf_data->y_local );
	PetscFree( mf_data );

	mf_data = PETSC_NULL;

	PetscFunctionReturn(0);
}

PetscErrorCode MatMult_MFOperatorGradient( Mat A, Vec x, Vec y )
{
	PetscErrorCode ierr;
	struct _p_MatrixFreeData_Gradient *mf_data;


	ierr = MatShellGetContext( A, (void**)&mf_data ); CHKERRQ(ierr);

	ierr = VecSet(y,0.0); CHKERRQ(ierr);
	ComputeGradient( mf_data->C, mf_data->da_vel, mf_data->da_pres, x, y, mf_data->y_local, mf_data->user );

	PetscFunctionReturn(0);
}

PetscErrorCode MatCreate_MFOperatorGradient( LaMEMVelPressureDA C, DM da_vel,Vec vel, DM da_pres,Vec pres, UserContext *data, Mat *A )
{
	PetscErrorCode ierr;
	struct _p_MatrixFreeData_Gradient *mf_data;
	MPI_Comm comm;
	PetscInt m,n, M,N;


	ierr = PetscMalloc( sizeof(struct _p_MatrixFreeData_Gradient), &mf_data ); CHKERRQ(ierr);
		   
	/* da */
	mf_data->da_vel = da_vel;
	mf_data->da_pres = da_pres;
	mf_data->C = C;

	/* sizes */
	VecGetLocalSize( vel, &m );
	VecGetSize( vel, &M );

	VecGetLocalSize( pres, &n );
	VecGetSize( pres, &N );


	mf_data->m = m;
	mf_data->n = n;

	mf_data->M = M;
	mf_data->N = N;

	/* stash */
	DMCreateLocalVector( da_vel, &mf_data->y_local );

	mf_data->user = data;


	PetscObjectGetComm( (PetscObject)da_vel, &comm );

	/* create matrix */
	ierr = MatCreateShell( comm, m,n, M,N, (void*)mf_data, A ); CHKERRQ(ierr);

	/* set operations */
	MatShellSetOperation( *A, MATOP_MULT,    (void(*)(void))MatMult_MFOperatorGradient );
	MatShellSetOperation( *A, MATOP_DESTROY, (void(*)(void))MatDestroy_MFOperatorGradient );


	PetscFunctionReturn(0);
}

void test_MFOperatorGradient( LaMEMVelPressureDA C, DM da_vel, DM da_pres, Vec pressure, Vec vel, Vec vel_local, UserContext *user )
{
	PetscScalar norm;
	Mat A;


	ComputeGradient( C, da_vel, da_pres, pressure, vel, vel_local, user );
	VecNorm( vel, NORM_2, &norm );
	PetscPrintf( PETSC_COMM_WORLD, "orig: ||u||_2 = %1.10e \n", norm );


	MatCreate_MFOperatorGradient( C, da_vel,vel,da_pres,pressure, user, &A );
	MatMult( A, pressure, vel );
	VecNorm( vel, NORM_2, &norm );
	PetscPrintf( PETSC_COMM_WORLD, "  MF: ||u||_2 = %1.10e \n", norm );

	MatDestroy(&A);

}





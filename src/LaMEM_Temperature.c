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

LaMEM_Temperature.c

Contains temperature-related routines that are used in StokesMG

$Id$
$Date$

This version is compatible with PETSc version 3.3-p2
The code cannot be distributed, used for commercial purposes or handed on,
without the explicit agreement of Boris Kaus.
*/

#include "LaMEM.h"
#include "LaMEM_Temperature.h"
#include "Quadrature.h"
#include "Elements.h"
#include "Assembly_FDSTAG.h"
#include "Utils.h"
#include "Material.h"

#ifdef TEMPERATURE

/*==========================================================================================================*/
/* Computes the kinematic matrix for temperature */
void ComputeKinmtxTemp( const PetscInt nnel, const PetscInt edof_temp, double **KinMtx, double **dhdPhys, double **dhdsVel,  double **InvJacob )
{
	/* It is assumed that the velocity, solution vector has the following form:
	*   [Vx(0) Vy(0) Vz(0)  Vx(1) Vy(1) Vz(1)   ...]  where the order is exactly as expected
	*   for the shape function */
	PetscInt	i, j, k;			//loop control
	// Compute dhdPhys - derivative versus natural coordinates
	for( i = 0; i < 3; i++)
	{	for( j = 0; j < nnel; j++)
		 	dhdPhys[i][j] = 0;
	}
	for( i = 0; i < 3; i++)
	{	for( j = 0; j < nnel; j++)
		{ 	for( k = 0; k < 3; k++)
				dhdPhys[i][j] += InvJacob[i][k]*dhdsVel[k][j];
		}
	}
	// Use dhdPhys to compute factors in front of strainrate
	for( i = 0; i < 3; i++){ for( j = 0; j < edof_temp; j++){ KinMtx[i][j] = 0; }} // initialize
	for( i = 0; i < nnel; i++){ KinMtx[0][i] = 	dhdPhys[0][i]; 	}// d/dx component
	for( i = 0; i < nnel; i++){ KinMtx[1][i] = 	dhdPhys[1][i];  }// d/dy component
	for( i = 0; i < nnel; i++){ KinMtx[2][i] = 	dhdPhys[2][i];  }// d/dz component
}

/*==========================================================================================================*/
/* Material matrix for temperature */
void ComputeMaterialMatrixTemp( const double k, double MatMtx[3][3] )
{
	PetscInt i,j;

	for (i=0;i<3;i++){ for (j=0;j<3;j++){ MatMtx[i][j]=0.0;	}}
	MatMtx[0][0] = k;	MatMtx[1][1] = k;	MatMtx[2][2] = k;
}
/*==========================================================================================================*/
/* performs TT_add = Kinmtx*Matmtx*Kinmtx */
void Multiply_Kin_Times_Mat_Times_Kin_Temp( const PetscInt edof_temp, double TT_add[MAX_edof_temp][MAX_edof_temp], double **Kinmtx, double Matmtx[3][3] )
{
	PetscInt 	i, j, k;					//loop control
	double	MatKin[3][MAX_edof_temp];	//Intermediate matrix

	/* Compute Intermediate matrix */
	for( i = 0; i < 3; i++){ for( j = 0; j < edof_temp; j++){ MatKin[i][j] = 0; }}
	for( i = 0; i < 3; i++){ for( j = 0; j < edof_temp; j++){for( k = 0; k < 3; k++){
				MatKin[i][j] +=  Matmtx[i][k]*Kinmtx[k][j];
	}}}
	/* Compute VV_add */
	for( i = 0; i < edof_temp; i++){ for( j = 0; j < edof_temp; j++){ TT_add[i][j] = 0; }}
	for( i = 0; i < edof_temp; i++){
		for( j = 0; j < edof_temp; j++){
			for( k = 0; k < 3; k++){
				TT_add[i][j] +=  Kinmtx[k][i]*MatKin[k][j];
			}
		}
	}
}
/*==========================================================================================================*/

/*==========================================================================================================*/
/* performs TT_mass = ShapeTemp*transpose(ShapeTemp) */
void Multiply_ShapeT_Times_ShapeT( const PetscInt edof_temp, double ShapeTemp[], double TT_mass[MAX_edof_temp][MAX_edof_temp] )
{
	PetscInt 	i, j;			//loop control

	/* Compute VV_add */
	for( i = 0; i < edof_temp; i++){ for( j = 0; j < edof_temp; j++){ TT_mass[i][j] = 0; }}
	for( i = 0; i < edof_temp; i++){
		for( j = 0; j < edof_temp; j++){
			TT_mass[i][j] +=  ShapeTemp[i]*ShapeTemp[j];
		}
	}

}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Create local2global for temperature equations */
#undef __FUNCT__
#define __FUNCT__ "CreateStencilInGlobalStiffnessTemp"
PetscErrorCode CreateStencilInGlobalStiffnessTemp( const PetscInt edof_temp, MatStencil *row,MatStencil *col,const PetscInt i,const PetscInt j, const PetscInt k )
{
	PetscInt ii;

	if(  (__ELEMENT_TYPE__ == ELEMENT_Q1P0) ||  (__ELEMENT_TYPE__ == ELEMENT_Q1Q1) ) {
		row[0 ].i = i;   row[0 ].j = j;   row[0 ].k = k; 	row[0 ].c = 0;  // T
		row[1 ].i = i+1; row[1 ].j = j;   row[1 ].k = k; 	row[1 ].c = 0;  // T
		row[2 ].i = i+1; row[2 ].j = j+1; row[2 ].k = k; 	row[2 ].c = 0;  // T
		row[3 ].i = i;   row[3 ].j = j+1; row[3 ].k = k; 	row[3 ].c = 0;  // T
		row[4 ].i = i;   row[4 ].j = j;   row[4 ].k = k+1; 	row[4 ].c = 0;  // T
		row[5 ].i = i+1; row[5 ].j = j;   row[5 ].k = k+1; 	row[5 ].c = 0;  // T
		row[6 ].i = i+1; row[6 ].j = j+1; row[6 ].k = k+1; 	row[6 ].c = 0;  // T
		row[7 ].i = i;   row[7 ].j = j+1; row[7 ].k = k+1; 	row[7 ].c = 0;  // T
		for (ii=0; ii<edof_temp; ii++){
			col[ii].i = row[ii].i;
			col[ii].j = row[ii].j;
			col[ii].k = row[ii].k;
			col[ii].c = row[ii].c;
		}
	}
	else if( __ELEMENT_TYPE__ == ELEMENT_Q2P1 ) {
		// Q2P1 element!
		row[0 ].i = i;     row[0 ].j = j;     row[0 ].k = k;   row[0 ].c = 0;  	// T
		row[1 ].i = i+2;   row[1 ].j = j;     row[1 ].k = k;   row[1 ].c = 0;  	// T
		row[2 ].i = i+2;   row[2 ].j = j+2;   row[2 ].k = k;   row[2 ].c = 0;  	// T
		row[3 ].i = i;     row[3 ].j = j+2;   row[3 ].k = k;   row[3 ].c = 0;  	// T
		row[4 ].i = i;     row[4 ].j = j;     row[4 ].k = k+2; row[4 ].c = 0;  	// T
		row[5 ].i = i+2;   row[5 ].j = j;     row[5 ].k = k+2; row[5 ].c = 0;  	// T
		row[6 ].i = i+2;   row[6 ].j = j+2;   row[6 ].k = k+2; row[6 ].c = 0;  	// T
		row[7 ].i = i;     row[7 ].j = j+2;   row[7 ].k = k+2; row[7 ].c = 0;  	// T
		row[8 ].i = i+1;   row[8 ].j = j;     row[8 ].k = k;   row[8 ].c = 0;  	// T
		row[9 ].i = i+2;   row[9 ].j = j+1;   row[9 ].k = k;   row[9 ].c = 0;  	// T
		row[10].i = i+1;   row[10].j = j+2;   row[10].k = k;   row[10].c = 0;  	// T
		row[11].i = i;     row[11].j = j+1;   row[11].k = k;   row[11].c = 0;  	// T
		row[12].i = i+1;   row[12].j = j;     row[12].k = k+2; row[12].c = 0;  	// T
		row[13].i = i+2;   row[13].j = j+1;   row[13].k = k+2; row[13].c = 0;  	// T
		row[14].i = i+1;   row[14].j = j+2;   row[14].k = k+2; row[14].c = 0;  	// T
		row[15].i = i;     row[15].j = j+1;   row[15].k = k+2; row[15].c = 0;  	// T
		row[16].i = i;     row[16].j = j;     row[16].k = k+1; row[16].c = 0;  	// T
		row[17].i = i+2;   row[17].j = j;     row[17].k = k+1; row[17].c = 0;  	// T
		row[18].i = i+2;   row[18].j = j+2;   row[18].k = k+1; row[18].c = 0;  	// T
		row[19].i = i;     row[19].j = j+2;   row[19].k = k+1; row[19].c = 0;  	// T
		row[20].i = i+1;   row[20].j = j+1;   row[20].k = k;   row[20].c = 0;  	// T
		row[21].i = i+1;   row[21].j = j+1;   row[21].k = k+2; row[21].c = 0;  	// T
		row[22].i = i;     row[22].j = j+1;   row[22].k = k+1; row[22].c = 0;  	// T
		row[23].i = i+2;   row[23].j = j+1;   row[23].k = k+1; row[23].c = 0;  	// T
		row[24].i = i+1;   row[24].j = j;     row[24].k = k+1; row[24].c = 0;  	// T
		row[25].i = i+1;   row[25].j = j+2;   row[25].k = k+1; row[25].c = 0;  	// T
		row[26].i = i+1;   row[26].j = j+1;   row[26].k = k+1; row[26].c = 0;  	// T


		for (ii=0; ii<edof_temp; ii++){
			col[ii].i = row[ii].i;
			col[ii].j = row[ii].j;
			col[ii].k = row[ii].k;
			col[ii].c = row[ii].c;
		}
	}
	else{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Unknown element type!" );
	}

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/
/* Get the Temperature for an element*/
#undef __FUNCT__
#define __FUNCT__ "GetTemperatureElement"
PetscErrorCode  GetTemperatureElement( PetscScalar ***temperature, PetscScalar Temp_element[], PetscInt i, PetscInt j, PetscInt k )
{

	if ( (__ELEMENT_TYPE__ == ELEMENT_Q1P0) || (__ELEMENT_TYPE__ == ELEMENT_Q1Q1) ||  (__ELEMENT_TYPE__ == ELEMENT_FDSTAG)) {
		Temp_element[0 ]	=	temperature[k  ][j  ][i  ];
		Temp_element[1 ]	=	temperature[k  ][j  ][i+1];
		Temp_element[2 ]	=	temperature[k  ][j+1][i+1];
		Temp_element[3 ]	=	temperature[k  ][j+1][i  ];
		Temp_element[4 ]	=	temperature[k+1][j  ][i  ];
		Temp_element[5 ]	=	temperature[k+1][j  ][i+1];
		Temp_element[6 ]	=	temperature[k+1][j+1][i+1];
		Temp_element[7 ]	=	temperature[k+1][j+1][i  ];
	}
	else if( __ELEMENT_TYPE__ == ELEMENT_Q2P1 ) {
		Temp_element[0 ]	=	temperature[k  ][j  ][i  ];
		Temp_element[1 ]	=	temperature[k  ][j  ][i+2];
		Temp_element[2 ]	=	temperature[k  ][j+2][i+2];
		Temp_element[3 ]	=	temperature[k  ][j+2][i  ];
		Temp_element[4 ]	=	temperature[k+2][j  ][i  ];
		Temp_element[5 ]	=	temperature[k+2][j  ][i+2];
		Temp_element[6 ]	=	temperature[k+2][j+2][i+2];
		Temp_element[7 ]	=	temperature[k+2][j+2][i  ];
		Temp_element[8 ]	=	temperature[k  ][j  ][i+1];
		Temp_element[9 ]	=	temperature[k  ][j+1][i+2];
		Temp_element[10]	=	temperature[k  ][j+2][i+1];
		Temp_element[11]	=	temperature[k  ][j+1][i  ];
		Temp_element[12]	=	temperature[k+2][j  ][i+1];
		Temp_element[13]	=	temperature[k+2][j+1][i+2];
		Temp_element[14]	=	temperature[k+2][j+2][i+1];
		Temp_element[15]	=	temperature[k+2][j+1][i  ];
		Temp_element[16]	=	temperature[k+1][j  ][i  ];
		Temp_element[17]	=	temperature[k+1][j  ][i+2];
		Temp_element[18]	=	temperature[k+1][j+2][i+2];
		Temp_element[19]	=	temperature[k+1][j+2][i  ];
		Temp_element[20]	=	temperature[k  ][j+1][i+1];
		Temp_element[21]	=	temperature[k+2][j+1][i+1];
		Temp_element[22]	=	temperature[k+1][j+1][i  ];
		Temp_element[23]	=	temperature[k+1][j+1][i+2];
		Temp_element[24]	=	temperature[k+1][j  ][i+1];
		Temp_element[25]	=	temperature[k+1][j+2][i+1];
		Temp_element[26]	=	temperature[k+1][j+1][i+1];
	}
	else{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Unknown element type!" );
	}

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/

/*==========================================================================================================*/
/* Set values to the right-hand-side */
#undef __FUNCT__
#define __FUNCT__ "SetValuesRHS_Temp"
PetscErrorCode   SetValuesRHS_Temp( PetscScalar ***rhs, PetscScalar T_RHS[], PetscInt i, PetscInt j, PetscInt k )
{

	if( (__ELEMENT_TYPE__ == ELEMENT_Q1P0) ||  (__ELEMENT_TYPE__ == ELEMENT_Q1Q1) ) {
		rhs[k  ][j  ][i  ]  = rhs[k  ][j  ][i  ] + T_RHS[0 ];
		rhs[k  ][j  ][i+1]  = rhs[k  ][j  ][i+1] + T_RHS[1 ];
		rhs[k  ][j+1][i+1]  = rhs[k  ][j+1][i+1] + T_RHS[2 ];
		rhs[k  ][j+1][i  ]  = rhs[k  ][j+1][i  ] + T_RHS[3 ];
		rhs[k+1][j  ][i  ]  = rhs[k+1][j  ][i  ] + T_RHS[4 ];
		rhs[k+1][j  ][i+1]  = rhs[k+1][j  ][i+1] + T_RHS[5 ];
		rhs[k+1][j+1][i+1]  = rhs[k+1][j+1][i+1] + T_RHS[6 ];
		rhs[k+1][j+1][i  ]  = rhs[k+1][j+1][i  ] + T_RHS[7 ];
	}
	else if( __ELEMENT_TYPE__ == ELEMENT_Q2P1 ) {
		rhs[k  ][j  ][i  ]  = rhs[k  ][j  ][i  ]  + T_RHS[0 ];
		rhs[k  ][j  ][i+2]  = rhs[k  ][j  ][i+2]  + T_RHS[1 ];
		rhs[k  ][j+2][i+2]  = rhs[k  ][j+2][i+2]  + T_RHS[2 ];
		rhs[k  ][j+2][i  ]  = rhs[k  ][j+2][i  ]  + T_RHS[3 ];
		rhs[k+2][j  ][i  ]  = rhs[k+2][j  ][i  ]  + T_RHS[4 ];
		rhs[k+2][j  ][i+2]  = rhs[k+2][j  ][i+2]  + T_RHS[5 ];
		rhs[k+2][j+2][i+2]  = rhs[k+2][j+2][i+2]  + T_RHS[6 ];
		rhs[k+2][j+2][i  ]  = rhs[k+2][j+2][i  ]  + T_RHS[7 ];
		rhs[k  ][j  ][i+1]  = rhs[k  ][j  ][i+1]  + T_RHS[8 ];
		rhs[k  ][j+1][i+2]  = rhs[k  ][j+1][i+2]  + T_RHS[9 ];
		rhs[k  ][j+2][i+1]  = rhs[k  ][j+2][i+1]  + T_RHS[10];
		rhs[k  ][j+1][i  ]  = rhs[k  ][j+1][i  ]  + T_RHS[11];
		rhs[k+2][j  ][i+1]  = rhs[k+2][j  ][i+1]  + T_RHS[12];
		rhs[k+2][j+1][i+2]  = rhs[k+2][j+1][i+2]  + T_RHS[13];
		rhs[k+2][j+2][i+1]  = rhs[k+2][j+2][i+1]  + T_RHS[14];
		rhs[k+2][j+1][i  ]  = rhs[k+2][j+1][i  ]  + T_RHS[15];
		rhs[k+1][j  ][i  ]  = rhs[k+1][j  ][i  ]  + T_RHS[16];
		rhs[k+1][j  ][i+2]  = rhs[k+1][j  ][i+2]  + T_RHS[17];
		rhs[k+1][j+2][i+2]  = rhs[k+1][j+2][i+2]  + T_RHS[18];
		rhs[k+1][j+2][i  ]  = rhs[k+1][j+2][i  ]  + T_RHS[19];
		rhs[k  ][j+1][i+1]  = rhs[k  ][j+1][i+1]  + T_RHS[20];
		rhs[k+2][j+1][i+1]  = rhs[k+2][j+1][i+1]  + T_RHS[21];
		rhs[k+1][j+1][i  ]  = rhs[k+1][j+1][i  ]  + T_RHS[22];
		rhs[k+1][j+1][i+2]  = rhs[k+1][j+1][i+2]  + T_RHS[23];
		rhs[k+1][j  ][i+1]  = rhs[k+1][j  ][i+1]  + T_RHS[24];
		rhs[k+1][j+2][i+1]  = rhs[k+1][j+2][i+1]  + T_RHS[25];
		rhs[k+1][j+1][i+1]  = rhs[k+1][j+1][i+1]  + T_RHS[26];
	}
	else{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Unknown element type!" );
	}

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/




/*==========================================================================================================*/
/* Compute the stiffness matrix for the temperature equation
 *
 * This function determines whether it is a FEM or a FDSTAG problem and calls the appropriate subroutines
 */
#undef __FUNCT__
#define __FUNCT__ "ComputeStiffnessMatrixRHSTemperature"
PetscErrorCode ComputeStiffnessMatrixRHSTemperature( LaMEMVelPressureDA C, DM da_temp, DM da, Mat J, UserContext user,
		Vec Temp_local, Vec Temp, Vec rhs_Temp_local, Vec rhs_Temp, PetscScalar dt )
{
	PetscErrorCode 		ierr;
	DAVPElementType 	element_type;

	PetscFunctionBegin;

	element_type = C->type;

	if (element_type==DAVP_FDSTAG) {
		/* FDSTAG formulation */
		ierr = ComputeStiffnessMatrixRHSTemperature_FDSTAG(da_temp, da, J, &user,
			Temp, rhs_Temp_local, rhs_Temp, dt ); CHKERRQ(ierr); // 	Compute stiffness matrix
	}
	else {
		/* FEM formulation */
		ierr = ComputeStiffnessMatrixRHSTemperature_FEM( C, da_temp, da, J, &user,
			Temp_local, Temp, rhs_Temp_local, rhs_Temp, dt ); CHKERRQ(ierr); // 	Compute stiffness matrix

	}

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Compute the temperature stiffness matrix using the FEM	*/
#undef __FUNCT__
#define __FUNCT__ "ComputeStiffnessMatrixRHSTemperature_FEM"
PetscErrorCode ComputeStiffnessMatrixRHSTemperature_FEM( LaMEMVelPressureDA C, DM da_temp, DM da, Mat J, UserContext *user,
	Vec Temp_local, Vec Temp, Vec rhs_Temp_local, Vec rhs_Temp, PetscScalar dt )
{
	PetscMPIInt     rank;
	PetscErrorCode  ierr;
	//  DMDA          da = (DM) dmmg->dm;
	PetscInt		i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs;
	PetscInt		xmp,ymp,zmp,xsp,ysp,zsp, xsg, ysg, zsg, xmg, ymg, zmg;
	PetscInt		ii,jj, intp,iiii, row_BC[1], *BC_row_vec, num_bc, num_elem;
	PetscInt		iel_x, iel_y, iel_z;
	PetscInt		*BC_row_vec_glob, phase;
	MatStencil		row[MAX_edof_temp], col[MAX_edof_temp];
	DMDACoor3d		***coords, coord_elem[MAX_nnel];
	DM			 cda;
	Vec			 gc;
	PetscScalar	 MatMtx[3][3], rho, cp,  cond, Q, ShearHeat;
	PetscScalar	 IntPoint[3][MAX_ngp_vel], IntWeight[MAX_ngp_vel], ShapeTemp[MAX_nnel], **dhdsTemp, Point[3];
	PetscScalar	 **Jacob, **InvJacob, DetJacob, **KinMtx, **dhdPhys;
	PetscScalar	 TT[MAX_edof_temp][MAX_edof_temp], TT_add[MAX_edof_temp][MAX_edof_temp], TT_mass[MAX_edof_temp][MAX_edof_temp], T_RHS[MAX_edof_temp], T_old[MAX_edof_temp];
	PetscScalar	 v[MAX_edof_temp*MAX_edof_temp], T_intp;
	ISLocalToGlobalMapping ltog;
	//MaterialsElement	***materials;
	PetscScalar		***temperature, ***rhs;
	PetscInt				edof_temp, nnel, nintp_1D, ngp_vel;
	DAVPElementType element_type;
	PetscScalar ***materials_array;
	MaterialsElementDynamic material_data;


	PetscFunctionBegin;

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	nnel        = C->nnel;
	nintp_1D    = C->nintp_1D;
	edof_temp   = C->edof_temp;
	ngp_vel     = C->ngp_vel;
	element_type = C->type;


	/* setup some memory */
	LaMEMCreate2dArray( 3, edof_temp, &KinMtx, PETSC_NULL );
	LaMEMCreate2dArray( 3, nnel, &dhdPhys, PETSC_NULL );
	LaMEMCreate2dArray( 3, 3, &Jacob, PETSC_NULL );
	LaMEMCreate2dArray( 3, 3, &InvJacob, PETSC_NULL );

	LaMEMCreate2dArray( 3, nnel, &dhdsTemp, PETSC_NULL );


	//	DASetUniformCoordinates(da,user->x_left,user->x_left+user->W,user->y_front,user->y_front + user->L,user->z_bot,user->z_bot + user->H);
	DMGetCoordinateDM(da,&cda);	//coordinates
	DMGetCoordinatesLocal(da,&gc);
	DMDAVecGetArray(cda,gc,&coords);

	// Copy global temperature solution ("old T") to local processor, including ghostpoints
	DMGlobalToLocalBegin(da_temp, Temp, INSERT_VALUES, Temp_local);
	DMGlobalToLocalEnd(da_temp,   Temp, INSERT_VALUES, Temp_local);
	ierr = DMDAVecGetArray(da_temp, Temp_local,	&temperature);		CHKERRQ(ierr);


	//	VecNorm(Temp, NORM_INFINITY, &max_temp);
	//	PetscPrintf(PETSC_COMM_WORLD,"max_temp = %g \n", max_temp);
	//	VecNorm(Temp_local, NORM_INFINITY, &max_temp);
	//	PetscPrintf(PETSC_COMM_WORLD,"max_temp = %g \n", max_temp);


	// Initialize rhs vector
	VecSet(rhs_Temp,0.0);  	VecSet(rhs_Temp_local,0.0);
	DMGlobalToLocalBegin(da_temp, rhs_Temp, INSERT_VALUES, rhs_Temp_local);
	DMGlobalToLocalEnd(da_temp,   rhs_Temp, INSERT_VALUES, rhs_Temp_local);
	ierr = DMDAVecGetArray(da_temp, rhs_Temp_local,		 &rhs);		CHKERRQ(ierr);

	/* print some info */
	ierr = DMDAGetInfo(da_temp,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0);	CHKERRQ(ierr);  	// # of nodes in all directions
	PetscPrintf(PETSC_COMM_WORLD,"#  Forming temperature stiffness matrix mx,my,mz=(%lld, %lld, %lld) ... ",(LLD)mx,(LLD)my,(LLD)mz);

	ierr = DMDAGetCorners(user->DA_Processors,&xsp,&ysp,&zsp,&xmp,&ymp,&zmp);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da_temp, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);

	DMDAVecGetArray(user->DA_Materials, user->Materials, &materials_array);

	/* Make a loop over all local elements, construct the element stiffness matrix */
	i=0; j=0; k=0;
	num_elem 	= 0;
	for (iel_z=zsp; iel_z<zsp+zmp; iel_z++){
		for (iel_y=ysp; iel_y<ysp+ymp; iel_y++){
			for(iel_x=xsp; iel_x<xsp+xmp; iel_x++){
				LaMEMSetMaterialDataMemoryFromArray( &material_data, iel_x-xsp,iel_y-ysp,iel_z-zsp, C->ngp_vel, materials_array );

				if( (element_type==DAVP_Q1P0) || (element_type==DAVP_Q1Q1) ){ // q1p0
					i = iel_x;
					j = iel_y;
					k = iel_z;
				}
				else if( (element_type==DAVP_Q2PM1L) || (element_type==DAVP_Q2PM1G) ){ // q2p1
					i = 2*iel_x;
					j = 2*iel_y;
					k = 2*iel_z;
				}
				else if (element_type==DAVP_FDSTAG){
					SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Routine does not  work for FDSTAG!" );
					ierr = MPI_Abort(PETSC_COMM_WORLD,1);  CHKERRQ(ierr);
				}

				/*
				i = element_counter_scale*iel_x;
				j = element_counter_scale*iel_y;
				k = element_counter_scale*iel_z;
				*/
				/*------------------------------------------------------
				* Form the local stiffness matrix
				* -----------------------------------------------------*/
				// 	if (num_elem==0){
				/* Extract coordinates of the local element in correct order */
				GetElementCoords(coord_elem, coords, i,j,k, 1);
				CorrectElementCoordsForPeriodicity( coord_elem, i, j, user,1 );

				/* Stuff that can be precomputed */
				IntegrationPoints( nnel, nintp_1D, IntPoint, IntWeight);

				/* Extract old temperature */
				GetTemperatureElement(temperature, T_old, i,j,k);


				/* Initialization */
				for(ii=0; ii<edof_temp; ii++){ for(jj=0; jj<edof_temp;  jj++){	TT[ii][jj]		=0.0;	}}

				for(ii=0; ii<edof_temp; ii++){ 									T_RHS[ii]		=0.0;	}

				//	    	/if (num_elem==0){
				/* Loop over integration points */
				for (intp=0; intp<ngp_vel; intp++){


					phase 		= 	( (PetscInt) material_data.Phases[intp]);				//  Dominant phase @ current integration point

					cond	 	=	user->PhaseProperties.T_Conductivity[phase];		//  Thermal conductivity
					Q 	 		=	user->PhaseProperties.RadioactiveHeat[phase];		//  Radioactive heat
					cp 	 		=	user->PhaseProperties.HeatCapacity[phase];			//  Heat capacity
					rho 	 	=	material_data.Density[intp];						// 	Density (might be T-dependent etc.)
					ShearHeat 	=	material_data.ShearHeat[intp];						//	Shear heating


					if (dt==0){
						dt=1e-10;		// hack, in case dt=0
					}

					Point[0]= IntPoint[0][intp];Point[1]= IntPoint[1][intp];Point[2]= IntPoint[2][intp];
					ComputeShapeFunctionVelocity(ShapeTemp, dhdsTemp, Point);				  // Temperature shape function is same as velocity shape function
					ComputeJacobianElement(nnel,dhdsTemp,coord_elem,Jacob,InvJacob,&DetJacob);  	  // Jacobian etc. of element
					ComputeKinmtxTemp( nnel, edof_temp, KinMtx, dhdPhys, dhdsTemp, InvJacob);					  // Kinematic matrix
					ComputeMaterialMatrixTemp(cond, MatMtx);								  // Material matrix for temperature


					/* Catch errors (for debugging only) */
					if (DetJacob<0){
						PetscPrintf(PETSC_COMM_WORLD,"Element too deformed in ComputeJacobianRHSTemperature\n");
						PetscPrintf(PETSC_COMM_WORLD,"In routine ComputeJacobian: Negative Jacobian on rank=%lld , DetJacob=%g  \n",(LLD)rank,DetJacob);
						PetscPrintf(PETSC_COMM_WORLD,"  element: [%lld,%lld,%lld]  \n",(LLD)iel_x,(LLD)iel_y,(LLD)iel_z);
						PetscPrintf(PETSC_COMM_WORLD," x-coords of local element: [%g,%g,%g,%g,%g,%g,%g,%g] \n",coord_elem[0].x,coord_elem[1].x,coord_elem[2].x,coord_elem[3].x,coord_elem[4].x,coord_elem[5].x,coord_elem[6].x,coord_elem[7].x);
						PetscPrintf(PETSC_COMM_WORLD," y-coords of local element: [%g,%g,%g,%g,%g,%g,%g,%g] \n",coord_elem[0].y,coord_elem[1].y,coord_elem[2].y,coord_elem[3].y,coord_elem[4].y,coord_elem[5].y,coord_elem[6].y,coord_elem[7].y);
						PetscPrintf(PETSC_COMM_WORLD," z-coords of local element: [%g,%g,%g,%g,%g,%g,%g,%g] \n",coord_elem[0].z,coord_elem[1].z,coord_elem[2].z,coord_elem[3].z,coord_elem[4].z,coord_elem[5].z,coord_elem[6].z,coord_elem[7].z);
						MPI_Abort(PETSC_COMM_WORLD,1);
					}

					// Compute TT matrix
					for(ii=0; ii<edof_temp; ii++){ for(jj=0; jj<edof_temp;  jj++){	TT_add[ii][jj]	=0.0;	}}
					for(ii=0; ii<edof_temp; ii++){ for(jj=0; jj<edof_temp;  jj++){	TT_mass[ii][jj]	=0.0;	}}
					Multiply_Kin_Times_Mat_Times_Kin_Temp( edof_temp, TT_add, KinMtx, MatMtx);			  //
					Multiply_ShapeT_Times_ShapeT( edof_temp, ShapeTemp, TT_mass);						  // ShapeT*ShapeT'
					for (ii=0; ii<edof_temp; ii++){
						for (jj=0; jj<edof_temp; jj++){
							TT[ii][jj]	=	TT[ii][jj] + ( TT_add[ii][jj] +	rho*cp/dt*TT_mass[ii][jj])*DetJacob*IntWeight[intp];
						}
					}

					// Compute T_RHS vector
					T_intp = 0.0;
					for (jj=0; jj<edof_temp; jj++){
						T_intp = T_intp + T_old[jj]*ShapeTemp[jj];
					}
					for (ii=0; ii<edof_temp; ii++){
						T_RHS[ii]	=	T_RHS[ii] + (rho*cp/dt*T_intp + Q + ShearHeat)*ShapeTemp[ii]*DetJacob*IntWeight[intp];
					}

					if (1==0){	//  debugging part
						PetscPrintf(PETSC_COMM_WORLD,"Integration point = [%g,%g,%g]\n", Point[0] , Point[1], Point[2]);


						PetscPrintf(PETSC_COMM_WORLD,"Jacob[0][:]= %g %g %g\n", Jacob[0][0] , Jacob[0][1], Jacob[0][2]);
						PetscPrintf(PETSC_COMM_WORLD,"Jacob[1][:]= %g %g %g\n", Jacob[1][0] , Jacob[1][1], Jacob[1][2]);
						PetscPrintf(PETSC_COMM_WORLD,"Jacob[2][:]= %g %g %g\n", Jacob[2][0] , Jacob[2][1], Jacob[2][2]);

						PetscPrintf(PETSC_COMM_WORLD,"InvJacob[0][:]= %g %g %g\n", InvJacob[0][0] , InvJacob[0][1], InvJacob[0][2]);
						PetscPrintf(PETSC_COMM_WORLD,"InvJacob[1][:]= %g %g %g\n", InvJacob[1][0] , InvJacob[1][1], InvJacob[1][2]);
						PetscPrintf(PETSC_COMM_WORLD,"InvJacob[2][:]= %g %g %g\n", InvJacob[2][0] , InvJacob[2][1], InvJacob[2][2]);

						PetscPrintf(PETSC_COMM_WORLD,"DetJacob= %g \n", DetJacob);

						PetscPrintf(PETSC_COMM_WORLD," x-coords of local element: [%g,%g,%g,%g,%g,%g,%g,%g] \n",coord_elem[0].x,coord_elem[1].x,coord_elem[2].x,coord_elem[3].x,coord_elem[4].x,coord_elem[5].x,coord_elem[6].x,coord_elem[7].x);
						PetscPrintf(PETSC_COMM_WORLD," y-coords of local element: [%g,%g,%g,%g,%g,%g,%g,%g] \n",coord_elem[0].y,coord_elem[1].y,coord_elem[2].y,coord_elem[3].y,coord_elem[4].y,coord_elem[5].y,coord_elem[6].y,coord_elem[7].y);
						PetscPrintf(PETSC_COMM_WORLD," z-coords of local element: [%g,%g,%g,%g,%g,%g,%g,%g] \n",coord_elem[0].z,coord_elem[1].z,coord_elem[2].z,coord_elem[3].z,coord_elem[4].z,coord_elem[5].z,coord_elem[6].z,coord_elem[7].z);

						for (iiii = 0; iiii<3; iiii++){
							PetscPrintf(PETSC_COMM_WORLD,"Kinmtx(%lld,:)= [%g %g %g %g %g %g %g %g ]\n",(LLD)(iiii+1),KinMtx[iiii][0], KinMtx[iiii][1],KinMtx[iiii][2], KinMtx[iiii][3],KinMtx[iiii][4], KinMtx[iiii][5],KinMtx[iiii][6], KinMtx[iiii][7]);
						}

						for (iiii = 0; iiii<8; iiii++){
							PetscPrintf(PETSC_COMM_WORLD,"TT(%lld,:)= [%g %g %g %g %g %g %g %g  ]\n",(LLD)iiii,TT[iiii][0], TT[iiii][1],TT[iiii][2], TT[iiii][3],TT[iiii][4], TT[iiii][5], TT[iiii][6], TT[iiii][7]);
						}


						MPI_Abort(PETSC_COMM_WORLD,1);
					}	// end of debugging part

				}

				/*------------------------------------------------------
				* End of forming the local stiffness matrix
				* -----------------------------------------------------*/
				// 	} // if  num_elem==0

				for (ii=0; ii<edof_temp; ii++){for (jj=0; jj<edof_temp; jj++){
						v[ii*edof_temp + jj] = TT[ii][jj];
				}}

				SetValuesRHS_Temp(rhs, T_RHS, i, j, k);


				// Add local stiffness to global stiffness matrix
				CreateStencilInGlobalStiffnessTemp(edof_temp,row,col,i,j,k);										// get local2global
				ierr = MatSetValuesStencil(J,edof_temp,row,edof_temp,col,v,ADD_VALUES);CHKERRQ(ierr);	// add to global stiffness

				num_elem = num_elem+1;
			}
		}
	}
	ierr = MatAssemblyBegin(J,	MAT_FINAL_ASSEMBLY);	CHKERRQ(ierr);	// Assemble global stiffness matrix
	ierr = MatAssemblyEnd(J,	MAT_FINAL_ASSEMBLY);	CHKERRQ(ierr);  // Assemble global stiffness matrix
	ierr = DMDAVecRestoreArray(da_temp, Temp_local, &temperature);	CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(da_temp,rhs_Temp_local,&rhs);	CHKERRQ(ierr);
	DMLocalToGlobalBegin(da_temp,rhs_Temp_local,ADD_VALUES,rhs_Temp);			DMLocalToGlobalEnd(da_temp,rhs_Temp_local,ADD_VALUES,rhs_Temp);

	DMDAVecRestoreArray(cda,gc,&coords);
	DMDAVecRestoreArray(user->DA_Materials, user->Materials, &materials_array);


	/* Set the boundary conditions in the matrix ---------------------------------------------------- */
	ierr = DMDAGetCorners(da_temp,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(da_temp,&xsg,&ysg,&zsg,&xmg,&ymg,&zmg);CHKERRQ(ierr);

	PetscMalloc( 3*(size_t)(2*xmg*ymg + 2*xmg*zmg + 2*ymg*zmg)*sizeof(PetscInt), 	&BC_row_vec);
	num_bc = 0;
	for (k=zs; k<zs+zm; k++){
		for (j=ys; j<ys+ym; j++){
			for(i=xs; i<xs+xm; i++){

				if (k==0    ){	// Bottom boundary
					GetGlobalIndexTemp(xmg, ymg, i-xsg  ,j-ysg  ,k-zsg  ,0  , &row_BC[0]);
					BC_row_vec[num_bc] 	 = row_BC[0];
					num_bc = num_bc+1;
				}
				if (k==mz-1    ){	// Top boundary
					GetGlobalIndexTemp(xmg, ymg, i-xsg  ,j-ysg  ,k-zsg  ,0  , &row_BC[0]);
					BC_row_vec[num_bc] 	 = row_BC[0];
					num_bc = num_bc+1;

				}

			}
		}
	}
	// set BCs in matrix
	PetscMalloc( (size_t)(num_bc)*sizeof(PetscInt), 	&BC_row_vec_glob);
	DMGetLocalToGlobalMapping(da_temp, &ltog);
	ISLocalToGlobalMappingApply(ltog,num_bc,BC_row_vec,BC_row_vec_glob);
	MatZeroRows(J,num_bc,BC_row_vec_glob, 1.0,PETSC_NULL,PETSC_NULL); 								// Delete row and put 1 at diagonal

	PetscFree(BC_row_vec);
	PetscFree(BC_row_vec_glob);
	/* End of setting the boundary conditions in the matrix ---------------------------------------------------- */

	/* Set BCs to the RHS vector ------------------------------------------------------------------------------ */
	// Initialize rhs vector
	ierr = DMDAVecGetArray(da_temp, rhs_Temp,		 &rhs);		CHKERRQ(ierr);
	num_bc = 0;
	for (k=zs; k<zs+zm; k++){
		for (j=ys; j<ys+ym; j++){
			for(i=xs; i<xs+xm; i++){
				// Only constant temperature conditions have to be set - zero flux follows automatically
				//  if nothing is set
				if (k==0    ){	// Bottom boundary - fixed T
					rhs[k][j][i]  = user->Temp_bottom;
				}
				if (k==mz-1    ){	// Top boundary - fixed T
					rhs[k][j][i]  = user->Temp_top;
				}

			}
		}
	}
	ierr = DMDAVecRestoreArray(da_temp,rhs_Temp,&rhs);	CHKERRQ(ierr);
	/* End of setting BCs to the RHS vector ------------------------------------------------------------------ */
	PetscPrintf(PETSC_COMM_WORLD," finished. \n");


//	DMDestroy(cda);
//	VecDestroy(gc);

	LaMEMDestroy2dArray(&KinMtx, PETSC_NULL );
	LaMEMDestroy2dArray(&dhdPhys, PETSC_NULL );
	LaMEMDestroy2dArray(&Jacob, PETSC_NULL );
	LaMEMDestroy2dArray(&InvJacob, PETSC_NULL );

	LaMEMDestroy2dArray(&dhdsTemp, PETSC_NULL );

	//	VecNorm(rhs_Temp,NORM_INFINITY,&test);
	//	PetscPrintf(PETSC_COMM_WORLD," test=%g \n",test);


	/* Hack for the manner in which BC's are set here */
	MatSetOption(J,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE);		//no new nonzeros

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/

/*==========================================================================================================*/
/* Find the nearest nodal point for a given particle with known natural coordinates [eta,zetha,phi]  */
#undef __FUNCT__
#define __FUNCT__ "FindNearestNode"
PetscErrorCode    FindNearestNode( PetscScalar eta, PetscScalar zetha, PetscScalar phi, PetscInt *ix_add_out, PetscInt *iy_add_out, PetscInt *iz_add_out )
{
	PetscInt	ix_add=0, iy_add=0, iz_add=0;

	if(  (__ELEMENT_TYPE__ == ELEMENT_Q1P0) ||  (__ELEMENT_TYPE__ == ELEMENT_Q1Q1) ||  (__ELEMENT_TYPE__ == ELEMENT_FDSTAG) ) {
		if (eta  <=0){ ix_add = 0; }
		if (eta  > 0){ ix_add = 1; }
		if (zetha<=0){ iy_add = 0; }
		if (zetha> 0){ iy_add = 1; }
		if (phi  <=0){ iz_add = 0; }
		if (phi  > 0){ iz_add = 1; }
	}
	else if( __ELEMENT_TYPE__ == ELEMENT_Q2P1 ){
		if ( (eta  <=-0.33)       		   	 ){ ix_add = 0; }
		if ( (eta  > -0.33) && (eta  < 0.33) ){ ix_add = 1; }
		if ( (eta  >  0.33) 				 ){ ix_add = 2; }
		if ( (zetha<=-0.33)       		   	 ){ iy_add = 0; }
		if ( (zetha> -0.33) && (zetha< 0.33) ){ iy_add = 1; }
		if ( (zetha>  0.33) 				 ){ iy_add = 2; }
		if ( (phi  <=-0.33)       		   	 ){ iz_add = 0; }
		if ( (phi  > -0.33) && (phi  < 0.33) ){ iz_add = 1; }
		if ( (phi  >  0.33) 				 ){ iz_add = 2; }
	}
	else{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Unknown element type!" );
	}
	*ix_add_out = ix_add;
	*iy_add_out = iy_add;
	*iz_add_out = iz_add;

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/* ====================================================================================================
* Interpolate particle temperatures to the nodes */
#undef __FUNCT__
#define __FUNCT__ "ParticleTemperatureToNodes"
PetscErrorCode ParticleTemperatureToNodes( LaMEMVelPressureDA C, UserContext *user, DM da_temp, Vec Temp_local, Vec Temp )
{
	PetscMPIInt     rank, size;
	PetscErrorCode	ierr;
	PetscInt		ipart, ielx,iely,ielz, ix_add, iy_add, iz_add;
	PetscInt		i,j,k,xs,ys,zs,xm,ym,zm,mz, ii,jj,kk, PrintErrorMessage;
	Particles		ParticleLocal;
	PetscScalar		***temperature, ***number, eta, zetha, phi;
	PetscScalar		MeanTemp, Number;
	Vec				Number_Global, Number_Local;
	PetscInt 		mod;
	DAVPElementType vpt_element_type;

	vpt_element_type = C->type;

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&size);


	// Create arrays to count the # of particles that have been added
	VecSet(Temp      ,0.0);						VecSet(Temp_local      ,0.0);
	VecDuplicate(Temp, 	     &Number_Global);	VecSet(Number_Global,0.0);
	VecDuplicate(Temp_local, &Number_Local);	VecSet(Number_Local ,0.0);

	// Copy numbers to local processor
	DMGlobalToLocalBegin(da_temp, Number_Global, INSERT_VALUES, Number_Local);
	DMGlobalToLocalEnd(da_temp,   Number_Global, INSERT_VALUES, Number_Local);
	DMDAVecGetArray(da_temp, Number_Local,	&number);

	// Copy global temperature solution ("old T") to local processor, including ghostpoints
	DMGlobalToLocalBegin(da_temp, Temp, INSERT_VALUES, Temp_local);
	DMGlobalToLocalEnd(da_temp,   Temp, INSERT_VALUES, Temp_local);
	DMDAVecGetArray(da_temp, Temp_local,	&temperature);
	DMDAGetGhostCorners(da_temp,&xs,&ys,&zs,&xm,&ym,&zm);

	/* Loop over all local particles */
	for (ipart=0; ipart<user->num_particle_local; ipart++){
		ParticleLocal = user->ParticlesLocal[ipart];

		/* Compute the element it belongs to  */
		ParticleLocal	=	user->ParticlesLocal[ipart];
		eta             =   ParticleLocal.eta;
		zetha           =   ParticleLocal.zetha;
		phi             =   ParticleLocal.phi;

		if( (vpt_element_type==DAVP_Q2PM1L) || (vpt_element_type==DAVP_Q2PM1G) ) {
			ielx = ((PetscInt) ParticleLocal.ix);
			LaMEMMod(ielx, 2, &mod);
			if (mod>0){  ielx = ielx-1;  };

			iely = ((PetscInt) ParticleLocal.iy);
			LaMEMMod(iely, 2, &mod);
			if (mod>0){  iely = iely-1;  };

			ielz = ((PetscInt) ParticleLocal.iz);
			LaMEMMod(ielz, 2, &mod);
			if (mod>0){  ielz = ielz-1;  };
		}
		else if( (vpt_element_type == DAVP_Q1P0) || (vpt_element_type == DAVP_Q1Q1) || (vpt_element_type == DAVP_FDSTAG)) {
			ielx = ((PetscInt) ParticleLocal.ix);
			iely = ((PetscInt) ParticleLocal.iy);
			ielz = ((PetscInt) ParticleLocal.iz);
		}
		else{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Unknown element type!" );
		}

		/* Find closest node */
		FindNearestNode(eta, zetha, phi, &ix_add, &iy_add, &iz_add);

		/* Add temperature to node */
		if ( ((ielx+ix_add)>=xs) && ((ielx+ix_add)<(xs+xm)) &&
				((iely+iy_add)>=ys) && ((iely+iy_add)<(ys+ym)) &&
				((ielz+iz_add)>=zs) && ((ielz+iz_add)<(zs+zm)) ){

			temperature[ielz+iz_add][iely+iy_add][ielx+ix_add] = temperature[ielz+iz_add][iely+iy_add][ielx+ix_add] + ParticleLocal.T;

			/* Keep track of how many temperatures have been added to a particular node */
			number[ielz+iz_add][iely+iy_add][ielx+ix_add] = number[ielz+iz_add][iely+iy_add][ielx+ix_add] + 1.0;
		}

	}


	DMDAVecRestoreArray(da_temp, Temp_local,	&temperature);
	DMLocalToGlobalBegin(da_temp,Temp_local,ADD_VALUES,Temp);
	DMLocalToGlobalEnd(da_temp,Temp_local,ADD_VALUES,Temp);

	DMDAVecRestoreArray(da_temp, Number_Local,	&number);
	DMLocalToGlobalBegin(da_temp,Number_Local, ADD_VALUES,Number_Global);
	DMLocalToGlobalEnd(da_temp,  Number_Local, ADD_VALUES,Number_Global);


	/* Average temperatures @ the nodes */
	// Copy numbers to local processor
	VecSet(Number_Local,0.0);
	DMGlobalToLocalBegin(da_temp, Number_Global, INSERT_VALUES, Number_Local);
	DMGlobalToLocalEnd(da_temp,   Number_Global, INSERT_VALUES, Number_Local);
	DMDAVecGetArray(da_temp, Number_Local,	&number);

	// Copy global temperature solution ("old T") to local processor, including ghostpoints
	DMGlobalToLocalBegin(da_temp, Temp, INSERT_VALUES, Temp_local);
	DMGlobalToLocalEnd(da_temp,   Temp, INSERT_VALUES, Temp_local);
	DMDAVecGetArray(da_temp, Temp_local,	&temperature);
	DMDAGetInfo(da_temp,0,0,0,&mz,0,0,0,0,0,0,0,0,0);				// # of nodes in all directions
	DMDAGetCorners(da_temp,&xs,&ys,&zs,&xm,&ym,&zm);
	for (k=zs; k<zs+zm; k++){
		for (j=ys; j<ys+ym; j++){
			for(i=xs; i<xs+xm; i++){

				if (number[k][j][i]>0){
					temperature[k][j][i] = temperature[k][j][i]/number[k][j][i];
				}

				// Set boundary conditions
				if (k==0){
					temperature[k][j][i] = user->Temp_bottom;
				}
				if (k==mz-1){
					temperature[k][j][i] = user->Temp_top;
				}

			}
		}
	}

	/* Make a loop over all nodes to find nodes w/o particles. Set temperature at these nodes
	*  to the mean temperature of the surrounding cells. Empty cells should not occur, since
	*  we inject particles in those cells.*/
	PrintErrorMessage = 1;
	for (k=zs; k<zs+zm-1; k++){
		for (j=ys; j<ys+ym-1; j++){
			for(i=xs; i<xs+xm-1; i++){
				if (number[k][j][i]==0){
					/* Empty cell detected */

					/* Compute mean T */
					MeanTemp=0; Number = 0;
					for (ii=i; ii<i+2; ii++){
						for (jj=j; jj<j+2; jj++){
							for (kk=k; kk<k+2; kk++){
								if (number[kk][jj][ii]>0){
									MeanTemp =  MeanTemp +	temperature[kk][jj][ii];
									Number   =	Number   +  1.0;
								}
							}
						}
					}

					if (Number>0){
						MeanTemp=MeanTemp/Number;
					} else {
						if (PrintErrorMessage==1){
							PetscPrintf(PETSC_COMM_WORLD," Serious problem: surrounding T nodes have also no particles !! \n");
							PrintErrorMessage = 0;  // print this message only once/timestep
						}
					}


					/* Set mean T */
					temperature[k][j][i] = MeanTemp;
				}
			}
		}
	}


	ierr = DMDAVecRestoreArray(da_temp, Temp_local,	&temperature);		CHKERRQ(ierr);
	VecSet(Temp,0);
	DMLocalToGlobalBegin(da_temp,Temp_local,INSERT_VALUES,Temp);
	DMLocalToGlobalEnd  (da_temp,Temp_local,INSERT_VALUES,Temp);
	ierr = DMDAVecRestoreArray(da_temp, Number_Local,	&number);		CHKERRQ(ierr);

	/* Clean up stuff */
	VecDestroy(&Number_Global);
	VecDestroy(&Number_Local);

	PetscFunctionReturn(0);

}
/* ==================================================================================================== */

/* ====================================================================================================
* Interpolate temperature from nodal points to particles.
* 	Two approaches are followed:
* 		(1) Interpolate complete temperature back to particles
* 		(2) Add temperature difference back to particle  												*/
#undef __FUNCT__
#define __FUNCT__ "TemperatureDiffNodesToParticles"
PetscErrorCode TemperatureDiffNodesToParticles( LaMEMVelPressureDA C, UserContext *user, DM da_temp, Vec Temp, Vec Temp_old )
{
	PetscInt		ipart, ix, iy, iz, i;
	PetscInt		xs,ys,zs,xm,ym,zm;
	PetscScalar	 	eta, zetha, phi;
	Particles		ParticleLocal;
	PetscScalar		***temperature, ***temperature_old, T_old[MAX_edof_temp], T_new[MAX_edof_temp];
	PetscScalar		ShapeTemp[MAX_edof_temp], **dhdsTemp, Point[3];
	PetscScalar		Tnew_point, Tdiff_point;
	Vec				Temp_old_local, Temp_local;
	PetscInt			mod;
	PetscInt					nnel, edof_temp;
	DAVPElementType 	element_type;

	nnel        = C->nnel;
	edof_temp   = C->edof_temp;
	element_type = C->type;


	LaMEMCreate2dArray( 3, edof_temp, &dhdsTemp, PETSC_NULL );


	/* Create vectors to store local temperature */
	DMCreateLocalVector(da_temp,&Temp_local);
	DMCreateLocalVector(da_temp,&Temp_old_local);

	DMGlobalToLocalBegin(da_temp, Temp, INSERT_VALUES, Temp_local);
	DMGlobalToLocalEnd(da_temp,   Temp, INSERT_VALUES, Temp_local);
	DMDAVecGetArray(da_temp, Temp_local,	&temperature);

	DMGlobalToLocalBegin(da_temp, Temp_old, INSERT_VALUES, Temp_old_local);
	DMGlobalToLocalEnd(da_temp,   Temp_old, INSERT_VALUES, Temp_old_local);
	DMDAVecGetArray(da_temp, Temp_old_local,	&temperature_old);

	DMDAGetCorners(da_temp,&xs,&ys,&zs,&xm,&ym,&zm);

	/* Loop over all local particles */
	for (ipart=0; ipart<user->num_particle_local; ipart++){
		ParticleLocal = user->ParticlesLocal[ipart];

		/* Compute the element it belongs to  */
		ParticleLocal	=	user->ParticlesLocal[ipart];

		if (ParticleLocal.phase>-4){
			if( (element_type==DAVP_Q1P0) || (element_type==DAVP_Q1Q1) || (element_type==DAVP_FDSTAG) ){
				ix				=	((PetscInt) ParticleLocal.ix);
				iy				=	((PetscInt) ParticleLocal.iy);
				iz				=	((PetscInt) ParticleLocal.iz);
			}
			else if( (element_type==DAVP_Q2PM1L) || (element_type==DAVP_Q2PM1G) ) {
				ix = ((PetscInt) ParticleLocal.ix);
				LaMEMMod(ix, 2, &mod);
				if (mod>0){  ix = ix-1;  };

				iy = ((PetscInt) ParticleLocal.iy);
				LaMEMMod(iy, 2, &mod);
				if (mod>0){  iy = iy-1;  };

				iz = ((PetscInt) ParticleLocal.iz);
				LaMEMMod(iz, 2, &mod);
				if (mod>0){  iz = iz-1;  };
			}
			else{
				SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Unknown element type!" );
			}

			if ( (ix>=xs) && (ix<(xs+xm)) &&  (iy>=ys) && (iy<(ys+ym)) &&  (iz>=zs) && (iz<(zs+zm))){
				eta             =   ParticleLocal.eta;
				zetha           =   ParticleLocal.zetha;
				phi             =   ParticleLocal.phi;

				// Extract old and new temperatures for given element
				GetTemperatureElement(temperature_old, T_old, ix,iy,iz);
				GetTemperatureElement(temperature    , T_new, ix,iy,iz);

				// Compute shape function at given natural coordinates
				Point[0] = eta; Point[1]= zetha; Point[2]=phi;
				ComputeShapeFunctionVelocity(ShapeTemp, dhdsTemp, Point);				 // Velocity shape function

				Tnew_point  = 0;
				Tdiff_point = 0;
				for (i=0; i<nnel; i++){
					Tnew_point 	= Tnew_point 	+ 	ShapeTemp[i]*T_new[i];
					Tdiff_point = Tdiff_point 	+ 	ShapeTemp[i]*(T_new[i]-T_old[i]);
				}

				// Possibility 1:
				ParticleLocal.T = Tnew_point;

				// Possibility 2:
				//ParticleLocal.T = ParticleLocal.T + Tdiff_point;

				user->ParticlesLocal[ipart] = ParticleLocal;
			}

		}
	}
	DMDAVecRestoreArray(da_temp, Temp_local,		&temperature	);
	DMDAVecRestoreArray(da_temp, Temp_old_local,	&temperature_old);

	LaMEMDestroy2dArray(&dhdsTemp, PETSC_NULL );


	VecDestroy(&Temp_local);
	VecDestroy(&Temp_old_local);

	PetscFunctionReturn(0);
}
/* ==================================================================================================== */


#endif

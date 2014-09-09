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

Assembly.c, contains the following subroutines:

ComputePressureLocal2Global 			-	Creates local2global numbering for pressure
ComputeVelocityLocal2Global				-	Creates local2global numbering for velocity
SetValuesRHS							-	Adds values to the RHS vectorGetVelocityElement
ComputeStiffnessMatrixes				-	Compute stiffness matrix -- steering routine (calls the FEM or FDSTAG subroutines to do computations)
ComputeStiffnessMatrixes_FEM 			-	Compute stiffness matrix using the FEM method
ComputeRHS								-	Compute RHS vector
SetBC_ElementStiffnessMatrix_Symmetric	-
Add_BC_ToLocalStiffness					-
SetBoundaryConditionsRHS				-	Set BC values to RHS vector

$Id$
$Date$

This version is compatible with PETSc version 3.3-p2
The code cannot be distributed, used for commercial purposes or handed on,
without the explicit agreement of Boris Kaus.
*/


#include "LaMEM.h"
#include "Assembly.h"
#include "Utils.h"
#include "Assembly_FDSTAG.h"
#include "Elements.h"
#include "Material.h"
#include "Quadrature.h"
#include "StokesOperators.h"

/*==========================================================================================================*/
/* Create Pressure local2global  */
#undef __FUNCT__
#define __FUNCT__ "ComputePressureLocal2Global"
PetscErrorCode ComputePressureLocal2Global(const PetscInt npres, PetscInt i, PetscInt j, PetscInt  k, PetscInt PRES_local2global[], Mat MATRIX)
{
	MatStencil 	row[MAX_npres];
	PetscInt 	kk;

	if  ( (__ELEMENT_TYPE__ == ELEMENT_Q1P0) | (__ELEMENT_TYPE__ == ELEMENT_Q2P1) ) {
		/* Discontinuous pressure shape function */
		for (kk=0; kk<npres; kk++){
			row[kk].i = i;
			row[kk].j = j;
			row[kk].k = k;
			row[kk].c = kk;
		}

	}
	else if (__ELEMENT_TYPE__ == ELEMENT_Q1Q1) {
		/* Continuous pressure shape function */

		row[0].i = i;    row[0].j = j;    row[0].k = k;   row[0].c = 0;
		row[1].i = i+1;  row[1].j = j;    row[1].k = k;   row[1].c = 0;
		row[2].i = i+1;  row[2].j = j+1;  row[2].k = k;   row[2].c = 0;
		row[3].i = i;    row[3].j = j+1;  row[3].k = k;   row[3].c = 0;
		row[4].i = i;    row[4].j = j;    row[4].k = k+1; row[4].c = 0;
		row[5].i = i+1;  row[5].j = j;    row[5].k = k+1; row[5].c = 0;
		row[6].i = i+1;  row[6].j = j+1;  row[6].k = k+1; row[6].c = 0;
		row[7].i = i;    row[7].j = j+1;  row[7].k = k+1; row[7].c = 0;

	}

	// Compute a local2global vector from the MatStencil values given above
	StencilToLocalNumbering(MATRIX,npres, row, PRES_local2global);

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Create local2global  */
#undef __FUNCT__
#define __FUNCT__ "ComputeVelocityLocal2Global"
PetscErrorCode ComputeVelocityLocal2Global( const PetscInt edof, MatStencil *row,MatStencil *col,PetscInt local2global[],
		const PetscInt i ,const PetscInt j , const PetscInt k, Mat MATRIX)
{
	PetscInt ii;

	if(  (__ELEMENT_TYPE__ == ELEMENT_Q1P0) |  (__ELEMENT_TYPE__ == ELEMENT_Q1Q1) ) {
	row[0 ].i = i;   row[0 ].j = j;   row[0 ].k = k; row[0 ].c = 0;  // Vx
	row[1 ].i = i;   row[1 ].j = j;   row[1 ].k = k; row[1 ].c = 1;  // Vy
	row[2 ].i = i;   row[2 ].j = j;   row[2 ].k = k; row[2 ].c = 2;  // Vz
	row[3 ].i = i+1; row[3 ].j = j;   row[3 ].k = k; row[3 ].c = 0;  // Vx
	row[4 ].i = i+1; row[4 ].j = j;   row[4 ].k = k; row[4 ].c = 1;  // Vy
	row[5 ].i = i+1; row[5 ].j = j;   row[5 ].k = k; row[5 ].c = 2;  // Vz
	row[6 ].i = i+1; row[6 ].j = j+1; row[6 ].k = k; row[6 ].c = 0;  // Vx
	row[7 ].i = i+1; row[7 ].j = j+1; row[7 ].k = k; row[7 ].c = 1;  // Vy
	row[8 ].i = i+1; row[8 ].j = j+1; row[8 ].k = k; row[8 ].c = 2;  // Vz
	row[9 ].i = i;   row[9 ].j = j+1; row[9 ].k = k; row[9 ].c = 0;  // Vx
	row[10].i = i;   row[10].j = j+1; row[10].k = k; row[10].c = 1;  // Vy
	row[11].i = i;   row[11].j = j+1; row[11].k = k; row[11].c = 2;  // Vz

	row[12].i = i;   row[12].j = j;   row[12].k = k+1; row[12].c = 0;  //Vx
	row[13].i = i;   row[13].j = j;   row[13].k = k+1; row[13].c = 1;  //Vy
	row[14].i = i;   row[14].j = j;   row[14].k = k+1; row[14].c = 2;  //Vz
	row[15].i = i+1; row[15].j = j;   row[15].k = k+1; row[15].c = 0;  //Vx
	row[16].i = i+1; row[16].j = j;   row[16].k = k+1; row[16].c = 1;  //Vy
	row[17].i = i+1; row[17].j = j;   row[17].k = k+1; row[17].c = 2;  //Vz
	row[18].i = i+1; row[18].j = j+1; row[18].k = k+1; row[18].c = 0;  //Vx
	row[19].i = i+1; row[19].j = j+1; row[19].k = k+1; row[19].c = 1;  //Vy
	row[20].i = i+1; row[20].j = j+1; row[20].k = k+1; row[20].c = 2;  //Vz
	row[21].i = i;   row[21].j = j+1; row[21].k = k+1; row[21].c = 0;  //Vx
	row[22].i = i;   row[22].j = j+1; row[22].k = k+1; row[22].c = 1;  //Vy
	row[23].i = i;   row[23].j = j+1; row[23].k = k+1; row[23].c = 2;  //Vz
	for (ii=0; ii<edof; ii++){
		col[ii].i = row[ii].i;
		col[ii].j = row[ii].j;
		col[ii].k = row[ii].k;
		col[ii].c = row[ii].c;
	}
	}

	if( __ELEMENT_TYPE__ == ELEMENT_Q2P1 ) {
	// Q2P1 element!
	row[0 ].i = i;     row[0 ].j = j;     row[0 ].k = k; row[0 ].c = 0;  // Vx
	row[1 ].i = i;     row[1 ].j = j;     row[1 ].k = k; row[1 ].c = 1;  // Vy
	row[2 ].i = i;     row[2 ].j = j;     row[2 ].k = k; row[2 ].c = 2;  // Vz
	row[3 ].i = i+2;   row[3 ].j = j;     row[3 ].k = k; row[3 ].c = 0;  // Vx
	row[4 ].i = i+2;   row[4 ].j = j;     row[4 ].k = k; row[4 ].c = 1;  // Vy
	row[5 ].i = i+2;   row[5 ].j = j;     row[5 ].k = k; row[5 ].c = 2;  // Vz
	row[6 ].i = i+2;   row[6 ].j = j+2;   row[6 ].k = k; row[6 ].c = 0;  // Vx
	row[7 ].i = i+2;   row[7 ].j = j+2;   row[7 ].k = k; row[7 ].c = 1;  // Vy
	row[8 ].i = i+2;   row[8 ].j = j+2;   row[8 ].k = k; row[8 ].c = 2;  // Vz
	row[9 ].i = i;     row[9 ].j = j+2;   row[9 ].k = k; row[9 ].c = 0;  // Vx
	row[10].i = i;     row[10].j = j+2;   row[10].k = k; row[10].c = 1;  // Vy
	row[11].i = i;     row[11].j = j+2;   row[11].k = k; row[11].c = 2;  // Vz
	row[12].i = i;     row[12].j = j;     row[12].k = k+2; row[12].c = 0;  // Vx
	row[13].i = i;     row[13].j = j;     row[13].k = k+2; row[13].c = 1;  // Vy
	row[14].i = i;     row[14].j = j;     row[14].k = k+2; row[14].c = 2;  // Vz
	row[15].i = i+2;   row[15].j = j;     row[15].k = k+2; row[15].c = 0;  // Vx
	row[16].i = i+2;   row[16].j = j;     row[16].k = k+2; row[16].c = 1;  // Vy
	row[17].i = i+2;   row[17].j = j;     row[17].k = k+2; row[17].c = 2;  // Vz
	row[18].i = i+2;   row[18].j = j+2;   row[18].k = k+2; row[18].c = 0;  // Vx
	row[19].i = i+2;   row[19].j = j+2;   row[19].k = k+2; row[19].c = 1;  // Vy
	row[20].i = i+2;   row[20].j = j+2;   row[20].k = k+2; row[20].c = 2;  // Vz
	row[21].i = i;     row[21].j = j+2;   row[21].k = k+2; row[21].c = 0;  // Vx
	row[22].i = i;     row[22].j = j+2;   row[22].k = k+2; row[22].c = 1;  // Vy
	row[23].i = i;     row[23].j = j+2;   row[23].k = k+2; row[23].c = 2;  // Vz
	row[24].i = i+1;   row[24].j = j;     row[24].k = k; row[24].c = 0;  // Vx
	row[25].i = i+1;   row[25].j = j;     row[25].k = k; row[25].c = 1;  // Vy
	row[26].i = i+1;   row[26].j = j;     row[26].k = k; row[26].c = 2;  // Vz
	row[27].i = i+2;   row[27].j = j+1;   row[27].k = k; row[27].c = 0;  // Vx
	row[28].i = i+2;   row[28].j = j+1;   row[28].k = k; row[28].c = 1;  // Vy
	row[29].i = i+2;   row[29].j = j+1;   row[29].k = k; row[29].c = 2;  // Vz
	row[30].i = i+1;   row[30].j = j+2;   row[30].k = k; row[30].c = 0;  // Vx
	row[31].i = i+1;   row[31].j = j+2;   row[31].k = k; row[31].c = 1;  // Vy
	row[32].i = i+1;   row[32].j = j+2;   row[32].k = k; row[32].c = 2;  // Vz
	row[33].i = i;     row[33].j = j+1;   row[33].k = k; row[33].c = 0;  // Vx
	row[34].i = i;     row[34].j = j+1;   row[34].k = k; row[34].c = 1;  // Vy
	row[35].i = i;     row[35].j = j+1;   row[35].k = k; row[35].c = 2;  // Vz
	row[36].i = i+1;   row[36].j = j;     row[36].k = k+2; row[36].c = 0;  // Vx
	row[37].i = i+1;   row[37].j = j;     row[37].k = k+2; row[37].c = 1;  // Vy
	row[38].i = i+1;   row[38].j = j;     row[38].k = k+2; row[38].c = 2;  // Vz
	row[39].i = i+2;   row[39].j = j+1;   row[39].k = k+2; row[39].c = 0;  // Vx
	row[40].i = i+2;   row[40].j = j+1;   row[40].k = k+2; row[40].c = 1;  // Vy
	row[41].i = i+2;   row[41].j = j+1;   row[41].k = k+2; row[41].c = 2;  // Vz
	row[42].i = i+1;   row[42].j = j+2;   row[42].k = k+2; row[42].c = 0;  // Vx
	row[43].i = i+1;   row[43].j = j+2;   row[43].k = k+2; row[43].c = 1;  // Vy
	row[44].i = i+1;   row[44].j = j+2;   row[44].k = k+2; row[44].c = 2;  // Vz
	row[45].i = i;     row[45].j = j+1;   row[45].k = k+2; row[45].c = 0;  // Vx
	row[46].i = i;     row[46].j = j+1;   row[46].k = k+2; row[46].c = 1;  // Vy
	row[47].i = i;     row[47].j = j+1;   row[47].k = k+2; row[47].c = 2;  // Vz
	row[48].i = i;     row[48].j = j;   row[48].k = k+1; row[48].c = 0;  // Vx
	row[49].i = i;     row[49].j = j;   row[49].k = k+1; row[49].c = 1;  // Vy
	row[50].i = i;     row[50].j = j;   row[50].k = k+1; row[50].c = 2;  // Vz
	row[51].i = i+2;     row[51].j = j;   row[51].k = k+1; row[51].c = 0;  // Vx
	row[52].i = i+2;     row[52].j = j;   row[52].k = k+1; row[52].c = 1;  // Vy
	row[53].i = i+2;     row[53].j = j;   row[53].k = k+1; row[53].c = 2;  // Vz
	row[54].i = i+2;     row[54].j = j+2; row[54].k = k+1; row[54].c = 0;  // Vx
	row[55].i = i+2;     row[55].j = j+2; row[55].k = k+1; row[55].c = 1;  // Vy
	row[56].i = i+2;     row[56].j = j+2; row[56].k = k+1; row[56].c = 2;  // Vz
	row[57].i = i;       row[57].j = j+2; row[57].k = k+1; row[57].c = 0;  // Vx
	row[58].i = i;       row[58].j = j+2; row[58].k = k+1; row[58].c = 1;  // Vy
	row[59].i = i;       row[59].j = j+2; row[59].k = k+1; row[59].c = 2;  // Vz
	row[60].i = i+1;     row[60].j = j+1; row[60].k = k; row[60].c = 0;  // Vx
	row[61].i = i+1;     row[61].j = j+1; row[61].k = k; row[61].c = 1;  // Vy
	row[62].i = i+1;     row[62].j = j+1; row[62].k = k; row[62].c = 2;  // Vz
	row[63].i = i+1;     row[63].j = j+1; row[63].k = k+2; row[63].c = 0;  // Vx
	row[64].i = i+1;     row[64].j = j+1; row[64].k = k+2; row[64].c = 1;  // Vy
	row[65].i = i+1;     row[65].j = j+1; row[65].k = k+2; row[65].c = 2;  // Vz
	row[66].i = i;       row[66].j = j+1; row[66].k = k+1; row[66].c = 0;  // Vx
	row[67].i = i;       row[67].j = j+1; row[67].k = k+1; row[67].c = 1;  // Vy
	row[68].i = i;       row[68].j = j+1; row[68].k = k+1; row[68].c = 2;  // Vz
	row[69].i = i+2;     row[69].j = j+1; row[69].k = k+1; row[69].c = 0;  // Vx
	row[70].i = i+2;     row[70].j = j+1; row[70].k = k+1; row[70].c = 1;  // Vy
	row[71].i = i+2;     row[71].j = j+1; row[71].k = k+1; row[71].c = 2;  // Vz
	row[72].i = i+1;     row[72].j = j;   row[72].k = k+1; row[72].c = 0;  // Vx
	row[73].i = i+1;     row[73].j = j;   row[73].k = k+1; row[73].c = 1;  // Vy
	row[74].i = i+1;     row[74].j = j;   row[74].k = k+1; row[74].c = 2;  // Vz
	row[75].i = i+1;     row[75].j = j+2; row[75].k = k+1; row[75].c = 0;  // Vx
	row[76].i = i+1;     row[76].j = j+2; row[76].k = k+1; row[76].c = 1;  // Vy
	row[77].i = i+1;     row[77].j = j+2; row[77].k = k+1; row[77].c = 2;  // Vz
	row[78].i = i+1;     row[78].j = j+1; row[78].k = k+1; row[78].c = 0;  // Vx
	row[79].i = i+1;     row[79].j = j+1; row[79].k = k+1; row[79].c = 1;  // Vy
	row[80].i = i+1;     row[80].j = j+1; row[80].k = k+1; row[80].c = 2;  // Vz

	for (ii=0; ii<edof; ii++){
		col[ii].i = row[ii].i;
		col[ii].j = row[ii].j;
		col[ii].k = row[ii].k;
		col[ii].c = row[ii].c;
	}
	}

	// Compute a local2global vector from the MatStencil values given above
	StencilToLocalNumbering(MATRIX,edof,row, local2global);

	PetscFunctionReturn(0);
}


/*==========================================================================================================*/
/* Set values to the right-hand-side */
#undef __FUNCT__
#define __FUNCT__ "SetValuesRHS"
PetscErrorCode SetValuesRHS( Field ***rhs, PetscScalar V_RHS[], PetscInt i, PetscInt j, PetscInt k )
{
	if(  (__ELEMENT_TYPE__ == ELEMENT_Q1P0) |  (__ELEMENT_TYPE__ == ELEMENT_Q1Q1) ) {
	rhs[k  ][j  ][i  ].Vx = rhs[k  ][j  ][i  ].Vx + V_RHS[0 ];
	rhs[k  ][j  ][i  ].Vy = rhs[k  ][j  ][i  ].Vy + V_RHS[1 ];
	rhs[k  ][j  ][i  ].Vz = rhs[k  ][j  ][i  ].Vz + V_RHS[2 ];
	rhs[k  ][j  ][i+1].Vx = rhs[k  ][j  ][i+1].Vx + V_RHS[3 ];
	rhs[k  ][j  ][i+1].Vy = rhs[k  ][j  ][i+1].Vy + V_RHS[4 ];
	rhs[k  ][j  ][i+1].Vz = rhs[k  ][j  ][i+1].Vz + V_RHS[5 ];
	rhs[k  ][j+1][i+1].Vx = rhs[k  ][j+1][i+1].Vx + V_RHS[6 ];
	rhs[k  ][j+1][i+1].Vy = rhs[k  ][j+1][i+1].Vy + V_RHS[7 ];
	rhs[k  ][j+1][i+1].Vz = rhs[k  ][j+1][i+1].Vz + V_RHS[8 ];
	rhs[k  ][j+1][i  ].Vx = rhs[k  ][j+1][i  ].Vx + V_RHS[9 ];
	rhs[k  ][j+1][i  ].Vy = rhs[k  ][j+1][i  ].Vy + V_RHS[10];
	rhs[k  ][j+1][i  ].Vz = rhs[k  ][j+1][i  ].Vz + V_RHS[11];

	rhs[k+1][j  ][i  ].Vx = rhs[k+1][j  ][i  ].Vx + V_RHS[12];
	rhs[k+1][j  ][i  ].Vy = rhs[k+1][j  ][i  ].Vy + V_RHS[13];
	rhs[k+1][j  ][i  ].Vz = rhs[k+1][j  ][i  ].Vz + V_RHS[14];
	rhs[k+1][j  ][i+1].Vx = rhs[k+1][j  ][i+1].Vx + V_RHS[15];
	rhs[k+1][j  ][i+1].Vy = rhs[k+1][j  ][i+1].Vy + V_RHS[16];
	rhs[k+1][j  ][i+1].Vz = rhs[k+1][j  ][i+1].Vz + V_RHS[17];
	rhs[k+1][j+1][i+1].Vx = rhs[k+1][j+1][i+1].Vx + V_RHS[18];
	rhs[k+1][j+1][i+1].Vy = rhs[k+1][j+1][i+1].Vy + V_RHS[19];
	rhs[k+1][j+1][i+1].Vz = rhs[k+1][j+1][i+1].Vz + V_RHS[20];
	rhs[k+1][j+1][i  ].Vx = rhs[k+1][j+1][i  ].Vx + V_RHS[21];
	rhs[k+1][j+1][i  ].Vy = rhs[k+1][j+1][i  ].Vy + V_RHS[22];
	rhs[k+1][j+1][i  ].Vz = rhs[k+1][j+1][i  ].Vz + V_RHS[23];
	}

	if( __ELEMENT_TYPE__ == ELEMENT_Q2P1 ) {
	rhs[k  ][j  ][i  ].Vx = rhs[k  ][j  ][i  ].Vx + V_RHS[0 ];
	rhs[k  ][j  ][i  ].Vy = rhs[k  ][j  ][i  ].Vy + V_RHS[1 ];
	rhs[k  ][j  ][i  ].Vz = rhs[k  ][j  ][i  ].Vz + V_RHS[2 ];
	rhs[k  ][j  ][i+2].Vx = rhs[k  ][j  ][i+2].Vx + V_RHS[3 ];
	rhs[k  ][j  ][i+2].Vy = rhs[k  ][j  ][i+2].Vy + V_RHS[4 ];
	rhs[k  ][j  ][i+2].Vz = rhs[k  ][j  ][i+2].Vz + V_RHS[5 ];
	rhs[k  ][j+2][i+2].Vx = rhs[k  ][j+2][i+2].Vx + V_RHS[6 ];
	rhs[k  ][j+2][i+2].Vy = rhs[k  ][j+2][i+2].Vy + V_RHS[7 ];
	rhs[k  ][j+2][i+2].Vz = rhs[k  ][j+2][i+2].Vz + V_RHS[8 ];
	rhs[k  ][j+2][i  ].Vx = rhs[k  ][j+2][i  ].Vx + V_RHS[9 ];
	rhs[k  ][j+2][i  ].Vy = rhs[k  ][j+2][i  ].Vy + V_RHS[10];
	rhs[k  ][j+2][i  ].Vz = rhs[k  ][j+2][i  ].Vz + V_RHS[11];
	rhs[k+2][j  ][i  ].Vx = rhs[k+2][j  ][i  ].Vx + V_RHS[12];
	rhs[k+2][j  ][i  ].Vy = rhs[k+2][j  ][i  ].Vy + V_RHS[13];
	rhs[k+2][j  ][i  ].Vz = rhs[k+2][j  ][i  ].Vz + V_RHS[14];
	rhs[k+2][j  ][i+2].Vx = rhs[k+2][j  ][i+2].Vx + V_RHS[15];
	rhs[k+2][j  ][i+2].Vy = rhs[k+2][j  ][i+2].Vy + V_RHS[16];
	rhs[k+2][j  ][i+2].Vz = rhs[k+2][j  ][i+2].Vz + V_RHS[17];
	rhs[k+2][j+2][i+2].Vx = rhs[k+2][j+2][i+2].Vx + V_RHS[18];
	rhs[k+2][j+2][i+2].Vy = rhs[k+2][j+2][i+2].Vy + V_RHS[19];
	rhs[k+2][j+2][i+2].Vz = rhs[k+2][j+2][i+2].Vz + V_RHS[20];
	rhs[k+2][j+2][i  ].Vx = rhs[k+2][j+2][i  ].Vx + V_RHS[21];
	rhs[k+2][j+2][i  ].Vy = rhs[k+2][j+2][i  ].Vy + V_RHS[22];
	rhs[k+2][j+2][i  ].Vz = rhs[k+2][j+2][i  ].Vz + V_RHS[23];
	rhs[k  ][j  ][i+1].Vx = rhs[k  ][j  ][i+1].Vx + V_RHS[24];
	rhs[k  ][j  ][i+1].Vy = rhs[k  ][j  ][i+1].Vy + V_RHS[25];
	rhs[k  ][j  ][i+1].Vz = rhs[k  ][j  ][i+1].Vz + V_RHS[26];
	rhs[k  ][j+1][i+2].Vx = rhs[k  ][j+1][i+2].Vx + V_RHS[27];
	rhs[k  ][j+1][i+2].Vy = rhs[k  ][j+1][i+2].Vy + V_RHS[28];
	rhs[k  ][j+1][i+2].Vz = rhs[k  ][j+1][i+2].Vz + V_RHS[29];
	rhs[k  ][j+2][i+1].Vx = rhs[k  ][j+2][i+1].Vx + V_RHS[30];
	rhs[k  ][j+2][i+1].Vy = rhs[k  ][j+2][i+1].Vy + V_RHS[31];
	rhs[k  ][j+2][i+1].Vz = rhs[k  ][j+2][i+1].Vz + V_RHS[32];
	rhs[k  ][j+1][i  ].Vx = rhs[k  ][j+1][i  ].Vx + V_RHS[33];
	rhs[k  ][j+1][i  ].Vy = rhs[k  ][j+1][i  ].Vy + V_RHS[34];
	rhs[k  ][j+1][i  ].Vz = rhs[k  ][j+1][i  ].Vz + V_RHS[35];
	rhs[k+2][j  ][i+1].Vx = rhs[k+2][j  ][i+1].Vx + V_RHS[36];
	rhs[k+2][j  ][i+1].Vy = rhs[k+2][j  ][i+1].Vy + V_RHS[37];
	rhs[k+2][j  ][i+1].Vz = rhs[k+2][j  ][i+1].Vz + V_RHS[38];
	rhs[k+2][j+1][i+2].Vx = rhs[k+2][j+1][i+2].Vx + V_RHS[39];
	rhs[k+2][j+1][i+2].Vy = rhs[k+2][j+1][i+2].Vy + V_RHS[40];
	rhs[k+2][j+1][i+2].Vz = rhs[k+2][j+1][i+2].Vz + V_RHS[41];
	rhs[k+2][j+2][i+1].Vx = rhs[k+2][j+2][i+1].Vx + V_RHS[42];
	rhs[k+2][j+2][i+1].Vy = rhs[k+2][j+2][i+1].Vy + V_RHS[43];
	rhs[k+2][j+2][i+1].Vz = rhs[k+2][j+2][i+1].Vz + V_RHS[44];
	rhs[k+2][j+1][i  ].Vx = rhs[k+2][j+1][i  ].Vx + V_RHS[45];
	rhs[k+2][j+1][i  ].Vy = rhs[k+2][j+1][i  ].Vy + V_RHS[46];
	rhs[k+2][j+1][i  ].Vz = rhs[k+2][j+1][i  ].Vz + V_RHS[47];
	rhs[k+1][j  ][i  ].Vx = rhs[k+1][j  ][i  ].Vx + V_RHS[48];
	rhs[k+1][j  ][i  ].Vy = rhs[k+1][j  ][i  ].Vy + V_RHS[49];
	rhs[k+1][j  ][i  ].Vz = rhs[k+1][j  ][i  ].Vz + V_RHS[50];
	rhs[k+1][j  ][i+2].Vx = rhs[k+1][j  ][i+2].Vx + V_RHS[51];
	rhs[k+1][j  ][i+2].Vy = rhs[k+1][j  ][i+2].Vy + V_RHS[52];
	rhs[k+1][j  ][i+2].Vz = rhs[k+1][j  ][i+2].Vz + V_RHS[53];
	rhs[k+1][j+2][i+2].Vx = rhs[k+1][j+2][i+2].Vx + V_RHS[54];
	rhs[k+1][j+2][i+2].Vy = rhs[k+1][j+2][i+2].Vy + V_RHS[55];
	rhs[k+1][j+2][i+2].Vz = rhs[k+1][j+2][i+2].Vz + V_RHS[56];
	rhs[k+1][j+2][i  ].Vx = rhs[k+1][j+2][i  ].Vx + V_RHS[57];
	rhs[k+1][j+2][i  ].Vy = rhs[k+1][j+2][i  ].Vy + V_RHS[58];
	rhs[k+1][j+2][i  ].Vz = rhs[k+1][j+2][i  ].Vz + V_RHS[59];
	rhs[k  ][j+1][i+1].Vx = rhs[k  ][j+1][i+1].Vx + V_RHS[60];
	rhs[k  ][j+1][i+1].Vy = rhs[k  ][j+1][i+1].Vy + V_RHS[61];
	rhs[k  ][j+1][i+1].Vz = rhs[k  ][j+1][i+1].Vz + V_RHS[62];
	rhs[k+2][j+1][i+1].Vx = rhs[k+2][j+1][i+1].Vx + V_RHS[63];
	rhs[k+2][j+1][i+1].Vy = rhs[k+2][j+1][i+1].Vy + V_RHS[64];
	rhs[k+2][j+1][i+1].Vz = rhs[k+2][j+1][i+1].Vz + V_RHS[65];
	rhs[k+1][j+1][i  ].Vx = rhs[k+1][j+1][i  ].Vx + V_RHS[66];
	rhs[k+1][j+1][i  ].Vy = rhs[k+1][j+1][i  ].Vy + V_RHS[67];
	rhs[k+1][j+1][i  ].Vz = rhs[k+1][j+1][i  ].Vz + V_RHS[68];
	rhs[k+1][j+1][i+2].Vx = rhs[k+1][j+1][i+2].Vx + V_RHS[69];
	rhs[k+1][j+1][i+2].Vy = rhs[k+1][j+1][i+2].Vy + V_RHS[70];
	rhs[k+1][j+1][i+2].Vz = rhs[k+1][j+1][i+2].Vz + V_RHS[71];
	rhs[k+1][j  ][i+1].Vx = rhs[k+1][j  ][i+1].Vx + V_RHS[72];
	rhs[k+1][j  ][i+1].Vy = rhs[k+1][j  ][i+1].Vy + V_RHS[73];
	rhs[k+1][j  ][i+1].Vz = rhs[k+1][j  ][i+1].Vz + V_RHS[74];
	rhs[k+1][j+2][i+1].Vx = rhs[k+1][j+2][i+1].Vx + V_RHS[75];
	rhs[k+1][j+2][i+1].Vy = rhs[k+1][j+2][i+1].Vy + V_RHS[76];
	rhs[k+1][j+2][i+1].Vz = rhs[k+1][j+2][i+1].Vz + V_RHS[77];
	rhs[k+1][j+1][i+1].Vx = rhs[k+1][j+1][i+1].Vx + V_RHS[78];
	rhs[k+1][j+1][i+1].Vy = rhs[k+1][j+1][i+1].Vy + V_RHS[79];
	rhs[k+1][j+1][i+1].Vz = rhs[k+1][j+1][i+1].Vz + V_RHS[80];
	}

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Compute the stiffness matrix
 *
 * This function determines whether it is a FEM or a FDSTAG problem and calls the appropriate subroutines
 */
#undef __FUNCT__
#define __FUNCT__ "ComputeStiffnessMatrixes"
PetscErrorCode ComputeStiffnessMatrixes(
		LaMEMVelPressureDA C, DM da, DM da_pres,
		Mat VV_MAT, 		Mat PP_MAT, 	Mat PV_MAT, 	Mat VP_MAT,
		UserContext user, PetscScalar dt, Vec rhs_add_BC,
		Vec rhs_add_BC_loc, Vec rhs_add_pBC, Vec pv_rhs_push, Vec ViscosityScaling, Mat approx_S, PetscInt *remesh )
{
	PetscErrorCode 		ierr;
	DAVPElementType 	element_type;
	PetscInt remesh1;

	PetscFunctionBegin;

	element_type = C->type;

	if (element_type==DAVP_FDSTAG) {
		/* FDSTAG formulation */
		ierr 	= ComputeStiffnessMatrixes_FDSTAG(da, da_pres, VV_MAT, PP_MAT, PV_MAT, VP_MAT, &user, dt,
				rhs_add_BC, rhs_add_BC_loc, pv_rhs_push, ViscosityScaling, approx_S ); CHKERRQ(ierr); // 	Compute stiffness matrix
		*remesh = 0;
	}
	else {
		/* FEM formulation */
		ierr 	= ComputeStiffnessMatrixes_FEM( C, 	da, da_pres, VV_MAT, PP_MAT, PV_MAT, VP_MAT, &user,
				rhs_add_BC, rhs_add_BC_loc, rhs_add_pBC, ViscosityScaling, approx_S, &remesh1 ); CHKERRQ(ierr); // 	Compute stiffness matrix
		*remesh = remesh1;
	}




	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Compute the stiffness matrix using FEM formulation */
#undef __FUNCT__
#define __FUNCT__ "ComputeStiffnessMatrixes_FEM"
PetscErrorCode ComputeStiffnessMatrixes_FEM(
		LaMEMVelPressureDA C, DM da, DM da_pres,
		Mat VV_MAT, 		Mat PP_MAT, 	Mat PV_MAT, 	Mat VP_MAT,
		UserContext *user, 	Vec rhs_add_BC,
		Vec rhs_add_BC_loc, Vec rhs_add_pBC, Vec ViscosityScaling, Mat approx_S, PetscInt *remesh)
{
	PetscMPIInt         rank;
	PetscErrorCode      ierr;
	//DMDA             	da = (DA) dmmg->dm;
	PetscInt       	 	i,j,k,xm,ym,zm,xs,ys,zs,phase;
	PetscInt		 	xmp,ymp,zmp,xsp,ysp,zsp;
	PetscInt			nel_x, nel_y, nel_z, nnode_x, nnode_y, nnode_z;
	PetscInt		 	ii,jj,intp, num_elem;
	PetscInt	     	iel_x, iel_y, iel_z, VEL_local2global[MAX_edof], PRES_local2global[MAX_npres];
	PetscInt 			VEL_local2global_localNumbering[MAX_edof], PRES_local2global_localNumbering[MAX_npres];
	const PetscInt		*P_globalindices, *Vel_globalindices;
	//PetscInt			xsg,ysg,zsg,xmg,ymg,zmg;
	MatStencil     		row[MAX_edof], col[MAX_edof];
	Field			 	***rhs_for_BoundaryConditions;
	DMDACoor3d		 	***coords, coord_elem[MAX_nnel];
	sBC 			 	BoundaryConditions;
	DM			 		cda;
	Vec			 		gc;
	PetscScalar  		mu_average, mu_viscous, mu_plastic, mu_eff, rho_eff, IntPoint[3][MAX_ngp_vel], IntWeight[MAX_ngp_vel], F_BC[MAX_edof], P_BC[MAX_npres];;
	PetscScalar			mu[C->ngp_vel], rho[C->ngp_vel], G[C->ngp_vel], T, mu_p[C->ngp_vel], mu_v[C->ngp_vel];
	//PetscScalar	 IntPoint[3][ngp_vel], IntWeight[ngp_vel], ShapeVel[nnel], dhdsVel[3][nnel], Point[3], ShapeP[npres];
	//PetscScalar	 Jacob[3][3], InvJacob[3][3], DetJacob[1], KinMtx[6][edof], dhdPhys[3][nnel], VV_add[edof][edof];
	PetscScalar	 		VV[MAX_edof][MAX_edof], VP[MAX_edof][MAX_npres], PV[MAX_npres][MAX_edof], PP[MAX_npres][MAX_npres], PP_prec[MAX_npres][MAX_npres];
	PetscScalar	 		vv[MAX_edof*MAX_edof], pp[MAX_npres*MAX_npres], vp[MAX_npres*MAX_edof], pv[MAX_npres*MAX_edof], o_on_mu[MAX_npres];
	PetscScalar 		P_lithos, minDiagonalRatio, minDiagonalRatio_local;
	PetscLogDouble		cputime_start, 	cputime_end;
	//MaterialsElement		***materials;
	PetscInt npres, edof, ngp_vel, nnel, nintp_1D;
	DAVPElementType element_type;
#ifdef WRITE_STIFFNESSMATRIX_FOR_DEBUGGING
	Vec 				InfoVec;
#endif
	PetscScalar ***materials_array;
	MaterialsElementDynamic material_data;
	PetscBool assemble_inverse_mass_matrix,flg;

	ISLocalToGlobalMapping P_ltogm, Vel_ltogm;

	PetscFunctionBegin;


	npres        = C->npres;
	edof         = C->edof;
	ngp_vel      = C->ngp_vel;
	nnel         = C->nnel;
	nintp_1D     = C->nintp_1D;
	element_type = C->type;

	assemble_inverse_mass_matrix = PETSC_FALSE;
	ierr = PetscOptionsGetBool( PETSC_NULL, "-assemble_inverse_mass_matrix", &assemble_inverse_mass_matrix, &flg ); CHKERRQ(ierr);

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	//	DASetUniformCoordinates(da,user->x_left,user->x_left+user->W,user->y_front,user->y_front + user->L,user->z_bot,user->z_bot + user->H);
	ierr = DMGetCoordinateDM(da,&cda);       CHKERRQ(ierr);	//coordinates
	ierr = DMGetCoordinatesLocal(da,&gc);    CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,gc,&coords);  CHKERRQ(ierr);

	BoundaryConditions = user->BC;

	
	ierr = DMGetLocalVector(da,&rhs_add_BC_loc);            CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da, rhs_add_BC, INSERT_VALUES, rhs_add_BC_loc); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(da,rhs_add_BC, INSERT_VALUES, rhs_add_BC_loc);    CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da,rhs_add_BC_loc,&rhs_for_BoundaryConditions);	    CHKERRQ(ierr);
	

	/* print some info */
	ierr = DMDAGetInfo(user->DA_Processors,	0,&nel_x	,&nel_y,	&nel_z,		0,0,0,0,0,0,0,0,0);	CHKERRQ(ierr); 	// # of elements in all directions
	ierr = DMDAGetInfo(da,				    0,&nnode_x,  &nnode_y,	&nnode_z,	0,0,0,0,0,0,0,0,0);	CHKERRQ(ierr);  // # of nodes 	 in all directions
	PetscTime(&cputime_start);
	PetscPrintf(PETSC_COMM_WORLD,"#  Forming FEM stiffness matrix mx,my,mz=(%lld, %lld, %lld) ... ",(LLD)nnode_x,(LLD)nnode_y,(LLD)nnode_z);


	ierr = DMDAGetCorners(user->DA_Processors,&xsp,&ysp,&zsp,&xmp,&ymp,&zmp);CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,      			   &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);

	/* Obtain mapping from current proc to global numbering for both velocity and pressure */
//	ierr = DMDAGetGlobalIndices(da_pres,&np  ,&P_globalindices); CHKERRQ(ierr);
//	ierr = DMDAGetGlobalIndices(da,     &nvel,&Vel_globalindices);  CHKERRQ(ierr);

	ierr = DMGetLocalToGlobalMapping(da_pres, &P_ltogm);   CHKERRQ(ierr);
	ierr = DMGetLocalToGlobalMapping(da,      &Vel_ltogm); CHKERRQ(ierr);

	ierr = ISLocalToGlobalMappingGetIndices(P_ltogm,   &P_globalindices);   CHKERRQ(ierr);
	ierr = ISLocalToGlobalMappingGetIndices(Vel_ltogm, &Vel_globalindices); CHKERRQ(ierr);


	ierr = DMDAVecGetArray(user->DA_Materials, user->Materials, &materials_array); CHKERRQ(ierr);

	// Make a loop over all local elements, construct the element stiffness matrix
	minDiagonalRatio_local = 	1e6;
	num_elem 	= 	0;
	for (iel_z=zsp; iel_z<zsp+zmp; iel_z++){
		for (iel_y=ysp; iel_y<ysp+ymp; iel_y++){
			for(iel_x=xsp; iel_x<xsp+xmp; iel_x++){
				/* load data */
				LaMEMSetMaterialDataMemoryFromArray( &material_data, iel_x-xsp,iel_y-ysp,iel_z-zsp, C->ngp_vel, materials_array );

				i = iel_x;	j = iel_y; k = iel_z;		// Initialize variables
				if( (element_type==DAVP_Q1P0) | (element_type==DAVP_Q1Q1) ) {
					i = iel_x;
					j = iel_y;
					k = iel_z;
				}
				else if( (element_type==DAVP_Q2PM1L) || (element_type==DAVP_Q2PM1G) ){
					i = 2*iel_x;
					j = 2*iel_y;
					k = 2*iel_z;
				}

				/*
				i = element_counter_scale * iel_x;
				j = element_counter_scale * iel_y;
				k = element_counter_scale * iel_z;
				*/

				/*------------------------------------------------------
				* Form the local stiffness matrix
				* -----------------------------------------------------*/
				//  	if (num_elem==0){
				/* Extract coordinates of the local element in correct order */
				GetElementCoords(coord_elem, coords, i,j,k, 1);
				CorrectElementCoordsForPeriodicity( coord_elem, i, j, user,1 );


				/* Stuff that can be precomputed */
				IntegrationPoints(nnel, nintp_1D, IntPoint, IntWeight);

				/* Initialization */
				for(ii=0; ii<edof;  ii++){ for(jj=0; jj<edof;  jj++){	     VV[ii][jj]=0.0;	}}
				for(ii=0; ii<edof;  ii++){ for(jj=0; jj<npres; jj++){	     VP[ii][jj]=0.0;	}}
				for(ii=0; ii<npres; ii++){ for(jj=0; jj<edof;  jj++){	     PV[ii][jj]=0.0;	}}
				for(ii=0; ii<npres; ii++){ for(jj=0; jj<npres; jj++){	     PP[ii][jj]=0.0;	}}
				for(ii=0; ii<npres; ii++){ for(jj=0; jj<npres; jj++){	PP_prec[ii][jj]=0.0;	}}




				/* Compute effective viscosity,elastic shear module and density at each integration point
				 * This is only done if we do a particle-based averaging.
				 * If we compute an analytical solution, these values are already explicitly hardcoded @ integration points
				 * */
				if (user->AnalyticalBenchmark==PETSC_FALSE){

					for (intp=0; intp<ngp_vel; intp++){
						phase 		= 	( (PetscInt) material_data.Phases[intp]);				//  Dominant phase @ current integration point
						G[intp] 	=	user->PhaseProperties.ElasticShearModule[phase];	//	Elastic shear module from lookup table

						T 			=	material_data.Temperature[intp];					// temperature

						// Compute effective density @ integration point
						ComputeEffectiveDensity(phase, user, T, &rho_eff);
						rho[intp] 	=	rho_eff;					//	Density

						// Estimate lithostatic pressure based on the depth to the top of the model (needed for plasticity cutoff's)
						P_lithos = (user->z_bot + user->H - coord_elem[0].z)*rho[0]*user->Gravity;

						mu_plastic 	=	material_data.PlasticViscosity[intp];
						mu_eff 		=	material_data.Viscosity[intp];

						// Effective viscosity might depend on various parameters
						ComputeEffectiveViscosity(phase, user, T, material_data.SecondInvariantDevStrainrate[intp],
								material_data.SecondInvariantDevStress[intp],  material_data.Pressure[intp], P_lithos,
								material_data.PlasticStrain[intp], user->PlasticityCutoff,
								&mu_viscous, &mu_plastic, &mu_eff);
						mu[intp] 	= 	mu_eff;
						mu_p[intp]	=	mu_plastic;
						mu_v[intp]	=	mu_viscous;


					}


					/* Compute an average element viscosity [required for stabilizing the element], and an average density (in a similar manner) */
					AverageViscosityDensityPerElement(mu, mu_v, mu_p, &mu_average, rho, user, ngp_vel);

					/* Store the true effective viscosity, the elasticity and the density that are used in the computations
					 * [since they might change @ every iteration step ] */
					for (intp=0; intp<ngp_vel; intp++){
						material_data.Viscosity[intp]			=	mu[intp];
						material_data.PlasticViscosity[intp]	=	mu_p[intp];
						material_data.TrueViscosity[intp] 		=	mu_v[intp];
						material_data.ElasticShearModule[intp]	=	G[intp];
						material_data.Density[intp]				=	rho[intp];
					}

				}
				else{
					mu_average=0;
					for (intp=0; intp<ngp_vel; intp++){
						mu[intp] 	= 	material_data.Viscosity[intp];
						mu_average 	=  mu_average + mu[intp]/((double) ngp_vel);
					}
				}

				/* Compute the element stiffness matrix */
				ComputeElementStiffnessMatrixes(C, VV, VP, PV, PP, PP_prec, IntPoint, IntWeight, mu, coord_elem, &minDiagonalRatio_local );

				//WriteElementStiffnessMatrixToDisk(C, VV, VP, PV, V_RHS, PP, mu_average, rho_g, coord_elem, num_elem, iel_x, iel_y, iel_z);
				/*------------------------------------------------------
				* End of forming the local stiffness matrix
				* -----------------------------------------------------*/


				/* Sofar, boundary conditions have not been taken into account
				*
				* If, however, the default PETSC method is employed, the resulting matrix is not symmetric.
				* That is a problem for the CG based Uzawa method, used in this code, which assumes the inner system
				* to be symmetric. It also unnecessarily slows down the code.
				*
				* The remedy is to make the element stiffness matrix symmetric, and add known contributions
				* to the RHS vector.
				*
				*/

#ifdef SYMMETRIZE_MATRIX_AT_ELEMENTLEVEL
				SetBC_ElementStiffnessMatrix_Symmetric( C, VV, VP, PV, F_BC, P_BC, coord_elem, iel_x, iel_y, iel_z, nel_x, nel_y, nel_z,  BoundaryConditions, user);

				// Add RHS additions to global vector
				SetValuesRHS(rhs_for_BoundaryConditions, F_BC, i, j, k);
#endif

				/* Obtain velocity and pressure local2global indices in LOCAL numbering */
				ComputeVelocityLocal2Global( edof, row,col,VEL_local2global_localNumbering,i,j,k,VV_MAT);
				ComputePressureLocal2Global( npres, iel_x,iel_y,iel_z,PRES_local2global_localNumbering, PP_MAT);

				for (ii=0; ii<edof; ii++){
					VEL_local2global[ii] = Vel_globalindices[VEL_local2global_localNumbering[ii]];
				}
				for (ii=0; ii<npres; ii++){
					PRES_local2global[ii] = P_globalindices[PRES_local2global_localNumbering[ii]];
				}

				// Put coefficients into the VV matrix
				for (ii=0; ii<edof; ii++){for (jj=0; jj<edof; jj++){
						vv[ii*edof + jj] = VV[ii][jj];
				}}
				ierr = MatSetValues(VV_MAT,edof,VEL_local2global,edof,VEL_local2global,vv,ADD_VALUES);CHKERRQ(ierr);


				// Put coefficients into PP matrix
				for (ii=0; ii<npres; ii++){for (jj=0; jj<npres; jj++){
					pp[ii*npres + jj] = PP[ii][jj];

					if (element_type==DAVP_Q1Q1) {
						/* scale with 1/mu */
						pp[ii*npres + jj] = -pp[ii*npres + jj]*(1.0/mu_average);
					}

				}}

				//
				if (assemble_inverse_mass_matrix) {
					double invPP[MAX_npres][MAX_npres];

					ComputeInversePP(PP,invPP);
					for (ii=0; ii<npres; ii++){
						for (jj=0; jj<npres; jj++){
							pp[ii*npres + jj] = invPP[ii][jj];
						}
					}
				}
				//


				ierr 		= MatSetValues(PP_MAT,npres,PRES_local2global,npres,PRES_local2global,pp,ADD_VALUES);CHKERRQ(ierr);

				// Put coefficients into VP matrix
				for (ii=0; ii<edof; ii++){
					for (jj=0; jj<npres; jj++){
						vp[ii*npres + jj] = VP[ii][jj];
					}
				}
				ierr 		= MatSetValues(VP_MAT,edof,VEL_local2global,npres,PRES_local2global,vp,ADD_VALUES);CHKERRQ(ierr);

				// Put coefficients into PV matrix
				for (ii=0; ii<npres; ii++){for (jj=0; jj<edof; jj++){
					pv[ii*edof + jj] = PV[ii][jj];
				}}
				ierr 		= MatSetValues(PV_MAT,npres,PRES_local2global,edof,VEL_local2global,pv,ADD_VALUES);CHKERRQ(ierr);

				/* Store the inverse of the element viscosity (for building a preconditioner @ a later stage) */
				/* Assemble viscosity scaling all at once */
				for (ii=0; ii<npres; ii++){
					o_on_mu[ii] = 1.0/mu_average;
				}
				ierr = VecSetValues( ViscosityScaling, npres, PRES_local2global, o_on_mu, INSERT_VALUES ); CHKERRQ(ierr);


				// compose 1/eta*M to be used as a preconditioning matrix
				for (ii=0; ii<npres; ii++)
				{	for (jj=0; jj<npres; jj++)
					{	pp[ii*npres + jj] = PP_prec[ii][jj];
						// scale with 1/mu
						pp[ii*npres + jj] = pp[ii*npres + jj]*(1.0/mu_average);
					}
				}
				ierr = MatSetValues(approx_S,npres,PRES_local2global,npres,PRES_local2global,pp,ADD_VALUES); CHKERRQ(ierr);


				#ifdef SYMMETRIZE_MATRIX_AT_ELEMENTLEVEL
				ierr = VecSetValues(rhs_add_pBC,npres,PRES_local2global,P_BC,ADD_VALUES); CHKERRQ(ierr);
				#endif

				num_elem 	= num_elem+1;

			}
		}
	}

	ierr = MatAssemblyBegin(VV_MAT,	MAT_FINAL_ASSEMBLY);	CHKERRQ(ierr);	// Assemble global velocity stiffness matrix
	ierr = MatAssemblyEnd(VV_MAT,	MAT_FINAL_ASSEMBLY);	CHKERRQ(ierr);  // Assemble global velocity stiffness matrix

	ierr = MatAssemblyBegin(PP_MAT,	MAT_FINAL_ASSEMBLY);	CHKERRQ(ierr);	// Assemble global pressure mass matrix
	ierr = MatAssemblyEnd(PP_MAT,	MAT_FINAL_ASSEMBLY);	CHKERRQ(ierr);  // Assemble global pressure mass matrix

	ierr = MatAssemblyBegin(VP_MAT,	MAT_FINAL_ASSEMBLY);	CHKERRQ(ierr);	// Assemble global VP matrix
	ierr = MatAssemblyEnd(VP_MAT,	MAT_FINAL_ASSEMBLY);	CHKERRQ(ierr);  // Assemble global VP matrix

	ierr = MatAssemblyBegin(PV_MAT,	MAT_FINAL_ASSEMBLY);	CHKERRQ(ierr);	// Assemble global PV matrix
	ierr = MatAssemblyEnd(PV_MAT,	MAT_FINAL_ASSEMBLY);	CHKERRQ(ierr);  // Assemble global PV matrix

	ierr = MatAssemblyBegin(approx_S,	MAT_FINAL_ASSEMBLY);	CHKERRQ(ierr);	// Assemble schur preconditioning matrix
	ierr = MatAssemblyEnd(approx_S,		MAT_FINAL_ASSEMBLY);	CHKERRQ(ierr);  // Assemble schur preconditioning matrix

	ierr = VecAssemblyBegin(ViscosityScaling);	CHKERRQ(ierr);	//Assemble ViscosityScaling
	ierr = VecAssemblyEnd(ViscosityScaling);	CHKERRQ(ierr);  //Assemble ViscosityScaling

	#ifdef SYMMETRIZE_MATRIX_AT_ELEMENTLEVEL
	ierr = VecAssemblyBegin(rhs_add_pBC);	CHKERRQ(ierr);		//Assemble rhs_add_pBC
	ierr = VecAssemblyEnd(rhs_add_pBC);		CHKERRQ(ierr);  	//Assemble rhs_add_pBC
	#endif

	ierr = DMDAVecRestoreArray(da,rhs_add_BC_loc,&rhs_for_BoundaryConditions);	CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da,rhs_add_BC_loc,ADD_VALUES,rhs_add_BC);	DMLocalToGlobalEnd(da,  rhs_add_BC_loc,ADD_VALUES, rhs_add_BC); CHKERRQ(ierr);	// update including ghost points
	ierr = DMRestoreLocalVector(da,&rhs_add_BC_loc); CHKERRQ(ierr);
	

	ierr = DMDAVecRestoreArray(cda,gc,&coords); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(user->DA_Materials, user->Materials, &materials_array); CHKERRQ(ierr);

	//ierr = DMDestroy(cda); CHKERRQ(ierr);	//coordinates
	//ierr = VecDestroy(gc); CHKERRQ(ierr);


	/* Set the boundary conditions - Note that the BC values are set in a separate routine! */
#if !defined(SYMMETRIZE_MATRIX_AT_ELEMENTLEVEL)

	PetscInt 				xsg, ysg, zsg, xmg, ymg, zmg, *BC_row_vec_glob, num_bc, *BC_row_vec,  row_BC[1];
	ISLocalToGlobalMapping 	ltog;


	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(da,&xsg,&ysg,&zsg,&xmg,&ymg,&zmg);CHKERRQ(ierr);

	ierr = PetscMalloc( 3*(2*xmg*ymg + 2*xmg*zmg + 2*ymg*zmg)*sizeof(PetscInt), 	&BC_row_vec); CHKERRQ(ierr);
	num_bc = 0;
	for (k=zs; k<zs+zm; k++){
		for (j=ys; j<ys+ym; j++){
			for(i=xs; i<xs+xm; i++){


				if (i==0 ){	// Left boundary
					if 		((user->BC.LeftBound==1) || (user->BC.LeftBound==6) || (user->BC.LeftBound==7)){ // free slip w. specified BG strainrate, thin-sheet, 3D growth rate
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,0  , &row_BC[0]);	// Vx (dof==0)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
					}
					else if 	(user->BC.LeftBound==2){ // no slip
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,0  , &row_BC[0]);	// Vx (dof==0)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,1  , &row_BC[0]);	// Vy (dof==1)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,2  , &row_BC[0]);	// Vz (dof==2)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
					}
					else if 	(user->BC.FrontBound==3){ // periodic

					}
					else if 	(user->BC.LeftBound==5){ // growthrate BC
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,0  , &row_BC[0]);	// Vx (dof==0)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,1  , &row_BC[0]);	// Vy (dof==1)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
					}
				}
				if (i==nnode_x-1 ){	// Right boundary
					if 		((user->BC.RightBound==1) || (user->BC.RightBound==7)){ // free slip w. specified BG strainrate, 3D growth rate
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,0  , &row_BC[0]);	// Vx (dof==0)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
					}
					else if 	(user->BC.RightBound==2){ // no slip
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,0  , &row_BC[0]);	// Vx (dof==0)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,1  , &row_BC[0]);	// Vy (dof==1)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,2  , &row_BC[0]);	// Vz (dof==2)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
					}
					else if 	(user->BC.RightBound==3){ // periodic

					}
					else if 	((user->BC.RightBound==5) || (user->BC.RightBound==6)){ // growthrate BC or thinsheet BC
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,0  , &row_BC[0]);	// Vx (dof==0)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,1  , &row_BC[0]);	// Vy (dof==1)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
					}
				}

				if (j==0 	){	// Front boundary
					if 		((user->BC.FrontBound==1) || (user->BC.FrontBound==5) || (user->BC.FrontBound==6) || (user->BC.FrontBound==7)){ // free slip w. specified BG strainrate, growth rate, thin-sheet
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,1  , &row_BC[0]);	// Vy (dof==1)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
					}
					else if 	(user->BC.FrontBound==2){ // no slip
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,0  , &row_BC[0]);	// Vx (dof==0)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,1  , &row_BC[0]);	// Vy (dof==1)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,2  , &row_BC[0]);	// Vz (dof==2)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
					}
					else if 	(user->BC.FrontBound==3){ // periodic

					}
					else if 	(user->BC.FrontBound==4){ // free slip with specified velocity over certain area
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,1  , &row_BC[0]);	// Vy (dof==1)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;

						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,2  , &row_BC[0]);	// Vz (dof==2)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
					}
					else if 	(user->BC.FrontBound==8){ // free slip with specified velocity over certain area
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,0  , &row_BC[0]);	// Vx (dof==0)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;

						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,1  , &row_BC[0]);	// Vy (dof==1)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
					}
				}
				if (j==nnode_y-1 ){	// Back boundary
					if 		((user->BC.BackBound==1) || (user->BC.BackBound==5) || (user->BC.BackBound==6) || (user->BC.BackBound==7)){ // free slip w. specified BG strainrate, growth rate, thin-sheet
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,1  , &row_BC[0]);	// Vy (dof==1)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
					}
					else if 	(user->BC.BackBound==2){ // no slip
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,0  , &row_BC[0]);	// Vx (dof==0)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,1  , &row_BC[0]);	// Vy (dof==1)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,2  , &row_BC[0]);	// Vz (dof==2)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
					}
					else if 		(user->BC.BackBound==3){ // free slip w. specified shear velocity

					}
				}
				if (k==0    ){	// Bottom boundary
					if 		((user->BC.LowerBound==1) || (user->BC.LowerBound==6) || (user->BC.LowerBound==7)){ // free slip w. specified BG strainrate
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,2  , &row_BC[0]);	// Vz (dof==2)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
					}
					else if 	(user->BC.LowerBound==2){ // no slip
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,0  , &row_BC[0]);	// Vx (dof==0)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,1  , &row_BC[0]);	// Vy (dof==0)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,2  , &row_BC[0]);	// Vz (dof==0)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
					}
					else if 	(user->BC.LowerBound==5){ // growthrate BC
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,1  , &row_BC[0]);	// Vy (dof==1)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,2  , &row_BC[0]);	// Vz (dof==2)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
					}
					else if 	(user->BC.LowerBound==8){ // no slip with prescribed horizontal BG strainrates
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,0  , &row_BC[0]);	// Vx (dof==0)  - prescribed
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,1  , &row_BC[0]);	// Vy (dof==1)- prescribed
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,2  , &row_BC[0]);	// Vz (dof==2)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
					}

				}
				if (k==nnode_z-1    ){	// Top boundary
					if 		(user->BC.UpperBound==1){ // free slip w. specified BG strainrate
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,2  , &row_BC[0]);	// Vz (dof==2)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
					}
					else if 	(user->BC.UpperBound==2){ // no slip
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,0  , &row_BC[0]);	// Vx (dof==0)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,1  , &row_BC[0]);	// Vy (dof==1)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,2  , &row_BC[0]);	// Vz (dof==2)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
					}
					else if 	(user->BC.UpperBound==5){ // growthrate BC
						GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,1  , &row_BC[0]);	// Vy (dof==1)
						BC_row_vec[num_bc] 	 = row_BC[0];
						num_bc = num_bc+1;
					}
				}
				if (user->BC.InternalBound>0){
					if (j==user->internalBC_node 	){	// Internal boundary
						if 	(user->BC.InternalBound==9){ // free slip with specified velocity over certain area
							GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,0  , &row_BC[0]);	// Vx (dof==0)
							BC_row_vec[num_bc] 	 = row_BC[0];
							num_bc = num_bc+1;

							GetGlobalIndex(xmg, ymg, zmg, i-xsg  ,j-ysg  ,k-zsg  ,1  , &row_BC[0]);	// Vy (dof==1)
							BC_row_vec[num_bc] 	 = row_BC[0];
							num_bc = num_bc+1;
						}
					}
				}

			}
		}
	}

	// set matrix
	ierr = PetscMalloc( (num_bc)*sizeof(PetscInt), 	&BC_row_vec_glob); CHKERRQ(ierr);
	ierr = DAGetISLocalToGlobalMapping(da, &ltog); CHKERRQ(ierr);
	ierr = ISLocalToGlobalMappingApply(ltog,num_bc,BC_row_vec,BC_row_vec_glob); CHKERRQ(ierr);
	ierr = MatZeroRows(VV_MAT,num_bc,BC_row_vec_glob, 1.0); CHKERRQ(ierr); 	// Delete row and put 1 at diagonal
	ierr = MatZeroRows(VP_MAT,num_bc,BC_row_vec_glob, 0.0); CHKERRQ(ierr); 	// Delete row altogether
	ierr = MatTranspose(VP_MAT,MAT_INITIAL_MATRIX,&PV_MAT); CHKERRQ(ierr);	// Ensure that PV is set accordingly



	ierr = PetscFree(BC_row_vec);
	ierr = PetscFree(BC_row_vec_glob);
#endif

	// indicate that no new nonzero's will be formed
	ierr = MatSetOption(VV_MAT,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE); CHKERRQ(ierr);	//
	ierr = MatSetOption(VP_MAT,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE); CHKERRQ(ierr);		//
	ierr = MatSetOption(PV_MAT,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE); CHKERRQ(ierr);		//
	ierr = MatSetOption(PP_MAT,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE); CHKERRQ(ierr);		//

	//MatSetOption(VV_MAT,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);		// generate error if new nonzero is created



	PetscTime(&cputime_end);
	PetscPrintf(PETSC_COMM_WORLD,"  finished. [%g s]\n",cputime_end-cputime_start);

	/* Compute the minimum diagonal ratio in the whole computational domain */
	ierr = MPI_Allreduce(&minDiagonalRatio_local,&minDiagonalRatio,1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD); CHKERRQ(ierr);
	
	PetscPrintf(PETSC_COMM_WORLD," minDiagonalRatio=%g [Critical = %g; change with -CriticalDiagonalRatio]  \n",minDiagonalRatio, user->CriticalDiagonalRatio);
	if ( (minDiagonalRatio<user->CriticalDiagonalRatio) && (user->GridAdvectionMethod==1) ){
		*remesh = 1;

		/* We should really remesh as the lagrangian mesh is too deformed. */
		PetscPrintf(PETSC_COMM_WORLD," Inducing a remesh step as minDiagonalRatio=%g [Critical=%g; change with -CriticalDiagonalRatio]  \n",minDiagonalRatio, user->CriticalDiagonalRatio);

	}
	else{
		*remesh = 0;
	}

	// FOR DEBUGGING: PRINT THE MATRIX TO DISK IF REQUIRED -------------------------
	//WriteStiffnessMatrixToDisk(da, VV_MAT, VP_MAT, PV_MAT, PP_MAT, PETSC_NULL, PETSC_NULL, user);

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/

/*==========================================================================================================*/
/* Compute the rhs vector
 *
 * This function determines whether it is a FEM or a FDSTAG problem and calls the appropriate subroutines
 */
#undef __FUNCT__
#define __FUNCT__ "ComputeRHS"
PetscErrorCode ComputeRHS( LaMEMVelPressureDA C, DM da, DM da_p, Vec b, Vec Pressure, UserContext user)
{
	PetscErrorCode 		ierr;
	DAVPElementType 	element_type;

	PetscFunctionBegin;

	element_type = C->type;

	if (element_type==DAVP_FDSTAG) {
		/* FDSTAG formulation */
		ierr = ComputeRHS_FDSTAG(da, b, &user);	CHKERRQ(ierr);
	}
	else {
		/* FEM formulation */
		ierr = ComputeRHS_FEM(C, da, da_p, b, Pressure, &user);	CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Compute the right-hand-side vector */
#undef __FUNCT__
#define __FUNCT__ "ComputeRHS_FEM"
PetscErrorCode ComputeRHS_FEM( LaMEMVelPressureDA C, DM da, DM da_p, Vec b, Vec Pressure, UserContext *user)
{
	DMDACoor3d		 		***coords, coord_elem[MAX_nnel], CoordIntp;
	DM			 			cda;
	Vec			 			gc, bloc;
	Field					***rhs;
	PetscErrorCode 			ierr;
	PetscInt       			i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs,ii, jj, kk, intp, nel_x, nel_y, nel_z;
	PetscInt		 		iel_x, iel_y, iel_z;
	PetscScalar	 			IntPoint[3][MAX_ngp_vel], IntWeight[MAX_ngp_vel],ShapeVel[MAX_nnel], **dhdsVel, Point[3];
	PetscScalar	 			**Jacob, **InvJacob, DetJacob,V_RHS[MAX_edof], rho_g;
	PetscScalar 	 		PP[MAX_npres][MAX_npres], ShapeP[MAX_npres], VP[MAX_edof][MAX_npres];
	PetscScalar	 			KinMtx[6][MAX_edof], **dhdPhys, Gravity_x, Gravity_z;
	PressureElem	 		***PressureNew;
	//MaterialsElement		***materials;
	PetscInt 					npres, edof, ngp_vel, nnel, nintp_1D;
	DAVPElementType 		element_type;
	PetscScalar 			***materials_array;
	MaterialsElementDynamic material_data;
	PetscScalar 			***PressureNew_array;
	PressureElemDynamic 	PressureNew_data;


	PetscFunctionBegin;

	nnel     		= C->nnel;
	nintp_1D 		= C->nintp_1D;
	npres   		= C->npres;
	edof    		= C->edof;
	ngp_vel 		= C->ngp_vel;
	element_type 	= C->type;

	Gravity_z 		=	 cos(user->GravityAngle/180*M_PI);
	Gravity_x 		=	-sin(user->GravityAngle/180*M_PI);

	LaMEMCreate2dArray( 3, nnel, &dhdsVel, PETSC_NULL );
	LaMEMCreate2dArray( 3, nnel, &dhdPhys, PETSC_NULL );
	LaMEMCreate2dArray( 3, 3, &Jacob, PETSC_NULL );
	LaMEMCreate2dArray( 3, 3, &InvJacob, PETSC_NULL );

	/* Get Coordinates */
	// 	DASetUniformCoordinates(da,user->x_left,user->x_left+user->W,user->y_front,user->y_front + user->L,user->z_bot,user->z_bot + user->H);
	ierr = DMGetCoordinateDM(da,&cda); 												CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(da,&gc); 										CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,gc,&coords); 											CHKERRQ(ierr);

	/* print some info  */
	ierr = DMDAGetInfo(user->DA_Processors,	0,&nel_x	,&nel_y,	&nel_z,		0,0,0,0,0,0,0,0,0);	CHKERRQ(ierr); 	// # of elements in all directions
	ierr = DMDAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0);									CHKERRQ(ierr);  	// # of nodes in all directions
	ierr = DMDAGetInfo(user->DA_Processors, 0, &mx, &my, &mz, 0,0,0,0,0,0,0,0,0);		CHKERRQ(ierr);
	ierr = DMDAGetCorners(user->DA_Processors,&xs,&ys,&zs,&xm,&ym,&zm);			CHKERRQ(ierr);

	ierr = DMGetLocalVector(da,&bloc); 													CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da, b, INSERT_VALUES, bloc); 							CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(da, b, INSERT_VALUES, bloc); 								CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da,bloc,&rhs);	CHKERRQ(ierr);

	ierr = DMDAVecGetArray(da_p, 			Pressure, 				&PressureNew); 			CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da_p, 			Pressure, 				&PressureNew_array); 	CHKERRQ(ierr);
	ierr = DMDAVecGetArray(user->DA_Materials, user->Materials, &materials_array); CHKERRQ(ierr);

	/* Compute RHS for each element; ignoring BC's */
	for (iel_z=zs; iel_z<zs+zm; iel_z++){
		for (iel_y=ys; iel_y<ys+ym; iel_y++){
			for(iel_x=xs; iel_x<xs+xm; iel_x++){
				LaMEMSetMaterialDataMemoryFromArray( &material_data, iel_x-xs,iel_y-ys,iel_z-zs, C->ngp_vel, materials_array );
				LaMEMSetPressueElemDataMemoryFromArray( &PressureNew_data, iel_x,iel_y,iel_z, npres, PressureNew_array );

				i = iel_x;	j = iel_y; k = iel_z;		// Initialize variables
				if( (element_type==DAVP_Q1P0) |  (element_type==DAVP_Q1Q1)  ){
					i = iel_x;
					j = iel_y;
					k = iel_z;
				}
				else if( (element_type==DAVP_Q2PM1L) || (element_type==DAVP_Q2PM1G) ){
					i = 2*iel_x;
					j = 2*iel_y;
					k = 2*iel_z;
				}

				/*
				i = element_counter_scale*iel_x;
				j = element_counter_scale*iel_y;
				k = element_counter_scale*iel_z;
				*/

				/* Extract coordinates of the local element in correct order */
				GetElementCoords(coord_elem, coords, i,j,k, 1);
				CorrectElementCoordsForPeriodicity( coord_elem, i, j, user, 1);


				/* Stuff that can be precomputed */
				IntegrationPoints( nnel, nintp_1D, IntPoint, IntWeight);

				/* Initialization */
				for(ii=0; ii<edof; ii++){ 								V_RHS[ii] =0.0;	 }
				for(ii=0; ii<npres;ii++){ for(jj=0; jj<npres; jj++){	PP[ii][jj]=0.0;	}}
				for(ii=0; ii<edof; ii++){ for(jj=0; jj<npres; jj++){	VP[ii][jj]=0.0;	}}



				/*---------------------------------------------------------------------*/
				for (intp=0; intp<ngp_vel; intp++){
					rho_g = material_data.Density[intp]*user->Gravity;	//  density*g


					if ((PetscAbsScalar(rho_g)<1e-11)  && (PetscAbsScalar(user->Gravity)>0) ){
						PetscPrintf(PETSC_COMM_SELF,"*** Emergency: abs(rho_g)<1e-11 & abs(g)>0 in forming rhs! material_data.Density[intp] = %g *** \n",material_data.Density[intp]);
					}
					if (isnan(rho_g)){
						PetscPrintf(PETSC_COMM_SELF,"*** Emergency: isnan(rho_g) in forming rhs! [ielx,iely,ielz,intp]=[%lld,%lld,%lld,%lld]*** \n", (LLD)iel_x, (LLD)iel_y, (LLD)iel_z, (LLD)intp);
					}

					Point[0]= IntPoint[0][intp];Point[1]= IntPoint[1][intp];Point[2]= IntPoint[2][intp];
					ComputeShapeFunctionVelocity(ShapeVel, dhdsVel, Point);				 // Velocity shape function
					ComputeCoordIntp(nnel,ShapeVel, coord_elem,  &CoordIntp);				 // Real coordinate of integration point
					ComputeJacobianElement(nnel,dhdsVel,coord_elem,Jacob,InvJacob,&DetJacob);  // Jacobian etc. of element
					ComputeShapeFunctionPressure(ShapeP, 	CoordIntp, Point);					 // Pressure shape function
					ComputeKinmtx(nnel,edof,KinMtx, dhdPhys, dhdsVel, InvJacob);					 // Kinematic matrix

					/* Catch errors (for debugging only) */
					if (DetJacob<0){

						PetscPrintf(PETSC_COMM_WORLD,"In routine ComputeRHS: Negative Jacobian on DetJacob=%g  \n",DetJacob);
						PetscPrintf(PETSC_COMM_WORLD,"  element: [%lld,%lld,%lld]  \n",(LLD)iel_x,(LLD)iel_y,(LLD)iel_z);
						if( element_type==DAVP_Q1P0 ){
							PetscPrintf(PETSC_COMM_WORLD," x-coords of local element: [%g,%g,%g,%g,%g,%g,%g,%g] \n",coord_elem[0].x,coord_elem[1].x,coord_elem[2].x,coord_elem[3].x,coord_elem[4].x,coord_elem[5].x,coord_elem[6].x,coord_elem[7].x);
							PetscPrintf(PETSC_COMM_WORLD," y-coords of local element: [%g,%g,%g,%g,%g,%g,%g,%g] \n",coord_elem[0].y,coord_elem[1].y,coord_elem[2].y,coord_elem[3].y,coord_elem[4].y,coord_elem[5].y,coord_elem[6].y,coord_elem[7].y);
							PetscPrintf(PETSC_COMM_WORLD," z-coords of local element: [%g,%g,%g,%g,%g,%g,%g,%g] \n",coord_elem[0].z,coord_elem[1].z,coord_elem[2].z,coord_elem[3].z,coord_elem[4].z,coord_elem[5].z,coord_elem[6].z,coord_elem[7].z);
						}
						else if( (element_type==DAVP_Q2PM1L) || (element_type==DAVP_Q2PM1G) ){
							PetscPrintf(PETSC_COMM_WORLD," x-coords of local element: [%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g] \n",coord_elem[0].x,coord_elem[1].x,coord_elem[2].x,coord_elem[3].x,coord_elem[4].x,coord_elem[5].x,coord_elem[6].x,coord_elem[7].x,
									coord_elem[8].x,coord_elem[9].x,coord_elem[10].x,coord_elem[11].x,coord_elem[12].x,coord_elem[13].x,coord_elem[14].x,coord_elem[15].x);

						}

						MPI_Abort(PETSC_COMM_WORLD,1);
					}


					// Compute PP matrix
					for (ii=0; ii<npres; ii++){for (jj=0; jj<npres; jj++){
							PP[ii][jj]	=	PP[ii][jj]	+	ShapeP[ii]*ShapeP[jj]*DetJacob*IntWeight[intp];
					}}

					// Compute VP matrix
					for (ii=0; ii<nnel;ii++){for (jj=0; jj<npres; jj++){for (kk=0; kk<3; kk++){
								VP[3*ii+kk][jj]	=	VP[3*ii+kk][jj]	+	dhdPhys[kk][ii]*ShapeP[jj]*DetJacob*IntWeight[intp];
					}}}


					// Compute contribution of gravity to the RHS
					for (ii=0; ii<nnel; ii++){
						V_RHS[3*ii+0] = V_RHS[3*ii+0]	-	Gravity_x*rho_g*ShapeVel[ii]*DetJacob*IntWeight[intp];  // x-component

						V_RHS[3*ii+2] = V_RHS[3*ii+2]	-	Gravity_z*rho_g*ShapeVel[ii]*DetJacob*IntWeight[intp];  // z-component
					}

				}
				/*---------------------------------------------------------------------*/

				/* Extract Pressure of the current element */
/*				for (ipres=0; ipres<npres; ipres++){
					P_element[ipres]  = PressureNew_data.P[ipres];
				}
*/

				/* Add pressure from last timestep to RHS */
				/*
				for (ii=0; ii<edof; ii++){
					for (jj=0; jj<npres; jj++){
						V_RHS[ii] = V_RHS[ii] - VP[ii][jj]*P_element[jj];
					}
				}
				*/

				// Add V_RHS to global vector
				SetValuesRHS(rhs, V_RHS, i, j, k);

			}
		}
	}
	ierr = DMDAVecRestoreArray(da,bloc,&rhs);	CHKERRQ(ierr); CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da,bloc,ADD_VALUES,b);	DMLocalToGlobalEnd(da,  bloc,ADD_VALUES,b); CHKERRQ(ierr);	// update including ghost points
	ierr = DMRestoreLocalVector(da,&bloc); CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(user->DA_Materials,   user->Materials, &materials_array); 	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(cda,							gc,						&coords); 			CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(da_p,				Pressure,				&PressureNew); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da_p, 			Pressure, 				&PressureNew_array); CHKERRQ(ierr);

	//ierr = DMDestroy( cda ); CHKERRQ(ierr);
	//ierr = VecDestroy( gc ); CHKERRQ(ierr);

	/* Set boundary conditions to RHS */
	SetBoundaryConditionsRHS(da, user, b, 1);

	LaMEMDestroy2dArray(&dhdsVel, PETSC_NULL );
	LaMEMDestroy2dArray(&dhdPhys, PETSC_NULL );
	LaMEMDestroy2dArray(&Jacob, PETSC_NULL );
	LaMEMDestroy2dArray(&InvJacob, PETSC_NULL );

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Set BC's on the element stiffness matrix and make this matrix symmetric */
#undef __FUNCT__
#define __FUNCT__ "SetBC_ElementStiffnessMatrix_Symmetric"
PetscErrorCode SetBC_ElementStiffnessMatrix_Symmetric( LaMEMVelPressureDA C, double VV[MAX_edof][MAX_edof], double VP[MAX_edof][MAX_npres],
		double PV[MAX_npres][MAX_edof], double F_BC[], double P_BC[],
		DMDACoor3d coord_elem[],
		PetscInt iel_x, PetscInt iel_y, PetscInt iel_z,
		PetscInt nel_x, PetscInt nel_y, PetscInt nel_z, sBC BoundaryConditions, UserContext *user)
{
	PetscInt 			i, ndof, nodes[MAX_nnel_1D*MAX_nnel_1D], k, SetDOF;
	PetscInt 		nn,mm,oo,pp, ApplyKinematicCondition_FromElement_Z;
	PetscScalar 	BoundaryVelocity[MAX_nnel_1D*MAX_nnel_1D][3], fac;
	PetscInt 		npres, edof, nnel_1D;
	DAVPElementType element_type;
	PetscBool 		flg;


	npres   = C->npres;
	edof    = C->edof;
	nnel_1D = C->nnel_1D;
	element_type = C->type;


	/* Initialize rhs vector, in which additions due to BC's will be added */
	for(i=0; i<edof;  i++){	F_BC[i] = 0;	}
	for(i=0; i<npres; i++){	P_BC[i] = 0;	}



	/* Compute a factor with which we multiply the strainrate boundary conditions
	 * This allows us to only apply the strainrate condition at the upper part of the domain
	 */
	fac = 1.0;
	PetscOptionsGetInt(PETSC_NULL ,"-ApplyKinematicCondition_FromElement_Z",&ApplyKinematicCondition_FromElement_Z	, &flg);	//
	if ((flg) && (ApplyKinematicCondition_FromElement_Z>0)){
		/* Compute the factor based on the current z-element */
		if (iel_z<=ApplyKinematicCondition_FromElement_Z){
			fac = ((PetscScalar) iel_z)/( (PetscScalar) ApplyKinematicCondition_FromElement_Z);
		}
	}


	/* Note: the way I number the nodes in an element is far from optimal.
	 * I'm sure one can come up with a better (and more general) technique;
	 * feel free to implement that and let me know.  BK
	 *
	 */

	// --------------------- Left BC -----------------------------
	if (iel_x==0){

		for(i=0; i<nnel_1D*nnel_1D; i++){ for(k=0; k<3; k++){
				BoundaryVelocity[i][k] = 0;		// Initialize all velocities to zero
		}}

		// Define left boundary nodes
		if( (element_type==DAVP_Q1P0) |  (element_type==DAVP_Q1Q1) ){
		nodes[0] = 0; nodes[1] = 3; nodes[2] = 4; nodes[3] = 7;
		}
		else if( (element_type==DAVP_Q2PM1L) || (element_type==DAVP_Q2PM1G) ){
		nodes[0] = 0;  nodes[1] = 11; nodes[2] = 3; nodes[3] = 16;
		nodes[4] = 22; nodes[5] = 19; nodes[6] = 4; nodes[7] = 15;
		nodes[8] = 7;
		}

		ndof 	= 0;
		SetDOF  = -1;
		if	(BoundaryConditions.LeftBound==1){ 			// free slip w. specified BG strainrate
			ndof 	= 1;
			SetDOF  = 0;								// Vx

			// define velocity @ all x-nodes on left side of model, based on BG strainrate
			for(i=0; i<nnel_1D*nnel_1D; i++){


				BoundaryVelocity[i][0] = -BoundaryConditions.Exx*coord_elem[nodes[i]].x*fac;
			}
		}
		else if 	(BoundaryConditions.LeftBound==2){ 	// no slip
			ndof = 3; 									// Vx,Vy,Vz
		}
		else if 	(BoundaryConditions.LeftBound==3){ 	// periodic

		}
		else if 	(BoundaryConditions.LeftBound==5){ 	// growth rate
			ndof 	= 2; 								// Vx,Vy
		}
		else if 	((BoundaryConditions.LeftBound==6) || (BoundaryConditions.LeftBound==7)){ // thin-sheet, growth rate
			ndof 	= 1;
			SetDOF  = 0; 								// Vx
		}
		else{
			SETERRQ1(PETSC_COMM_WORLD, 1,"Unknown left boundary condition BC.LeftBound=%lld",(LLD)(BoundaryConditions.LeftBound));
		}

		if (ndof>0){
			// Zero rows and collumns of local stiffness matrix
			Add_BC_ToLocalStiffness(C,VV, VP, PV, F_BC, P_BC, nodes, ndof, SetDOF, BoundaryVelocity);
		}

	}
	/* ----------------------------------------------------------- */


	// --------------------- Right BC -----------------------------
	if (iel_x==nel_x-1){

		for(i=0; i<nnel_1D*nnel_1D; i++){ for(k=0; k<3; k++){
				BoundaryVelocity[i][k] = 0;		// Initialize all velocities to zero
		}}

		// Define right boundary nodes
		if( (element_type==DAVP_Q1P0) |  (element_type==DAVP_Q1Q1)  ){
		nodes[0] = 1; nodes[1] = 2; nodes[2] = 5; nodes[3] = 6;
		}
		else if( (element_type==DAVP_Q2PM1L) || (element_type==DAVP_Q2PM1G) ){
		nodes[0] = 1;  nodes[1] = 9;  nodes[2] = 2; nodes[3] = 17;
		nodes[4] = 23; nodes[5] = 18; nodes[6] = 5; nodes[7] = 13;
		nodes[8] = 6;
		}

		ndof 	= 0;
		SetDOF  = -1;
		if	((BoundaryConditions.RightBound==1) || (BoundaryConditions.RightBound==7)){ 			// free slip w. specified BG strainrate, growth rate
			ndof 	= 1;
			SetDOF  = 0; 									// Vx

			// define velocity @ all x-nodes on right side of model, based on BG strainrate
			for(i=0; i<nnel_1D*nnel_1D; i++){
				BoundaryVelocity[i][0] = -BoundaryConditions.Exx*coord_elem[nodes[i]].x;
			}
		}
		else if 	(BoundaryConditions.RightBound==2){ 	// no slip
			ndof = 3;										// Vx,Vy,Vz
		}
		else if 	(BoundaryConditions.RightBound==3){ 	// periodic

		}
		else if 	(BoundaryConditions.RightBound==5){ 	// growth rate
			ndof 	= 2; 									// Vx,Vy
			for(i=0; i<nnel_1D*nnel_1D; i++){
				BoundaryVelocity[i][0] = -BoundaryConditions.Exx*coord_elem[nodes[i]].x*fac;
				if(user->Setup.Model==7){ BoundaryVelocity[i][0] = -BoundaryConditions.Exx*fac;}
			}
			for(i=0; i<nnel_1D*nnel_1D; i++){
				BoundaryVelocity[i][1] = 0.0;
			}
		}
		else if 	(BoundaryConditions.RightBound==6){ 	// thin-sheet
			ndof 	= 2; 									// Vx,Vy
			nn = (PetscInt)round(((PetscScalar)nel_y-1.0)/4.0);
			mm = (PetscInt)round(((PetscScalar)nel_y-1.0)/2.0);

			for(i=0; i<nnel_1D*nnel_1D; i++){
				if(iel_y<=nn){
					BoundaryVelocity[i][0] = -BoundaryConditions.Exx*fac*user->W;
				}
				else if(iel_y>nn && iel_y<mm){
					BoundaryVelocity[i][0] = -fac*BoundaryConditions.Exx*user->W*(cos(2.0*M_PI*((coord_elem[nodes[i]].y)/user->L-(1.0/4.0))))*(cos(2.0*M_PI*((coord_elem[nodes[i]].y)/user->L-(1.0/4.0))));
				}
				else if(iel_y>=mm){
					BoundaryVelocity[i][0] = 0.0;
				}
			}
			for(i=0; i<nnel_1D*nnel_1D; i++){
				BoundaryVelocity[i][1] = 0.0;
			}
		}
		else{
			SETERRQ1(PETSC_COMM_WORLD, 1,"Unknown right boundary condition BC.RightBound=%lld",(LLD)(BoundaryConditions.RightBound));
		}
		if (ndof>0){
			// Zero rows and collumns of local stiffness matrix
			Add_BC_ToLocalStiffness(C, VV, VP, PV, F_BC, P_BC, nodes, ndof, SetDOF, BoundaryVelocity);
		}

	}
	// -----------------------------------------------------------


	// --------------------- Front BC -----------------------------
	if (iel_y==0){

		for(i=0; i<nnel_1D*nnel_1D; i++){ for(k=0; k<3; k++){
				BoundaryVelocity[i][k] = 0;		// Initialize all velocities to zero
		}}

		// Define front boundary nodes
		if( (element_type==DAVP_Q1P0) |  (element_type==DAVP_Q1Q1)  ){
		nodes[0] = 0; nodes[1] = 1; nodes[2] = 4; nodes[3] = 5;
		}
		else if( (element_type==DAVP_Q2PM1L) || (element_type==DAVP_Q2PM1G) ){
		nodes[0] = 0;  nodes[1] = 8;  nodes[2] = 1; nodes[3] = 16;
		nodes[4] = 24; nodes[5] = 17; nodes[6] = 4; nodes[7] = 12;
		nodes[8] = 5;
		}

		ndof 	=	0;
		SetDOF  =  -1;
		if	(BoundaryConditions.FrontBound==1){ 			// free slip w. specified BG strainrate
			ndof 	= 	1;
			SetDOF  =   1; 									// Vy
			// define velocity @ all y-nodes on front side of model, based on BG strainrate
			for(i=0; i<nnel_1D*nnel_1D; i++){
				BoundaryVelocity[i][0] = -fac*BoundaryConditions.Eyy*coord_elem[nodes[i]].y;
			}
		}
		else if 	(BoundaryConditions.FrontBound==2){ 	// no slip
			ndof = 3; 										// Vx,Vy,Vz
		}
		else if 	(BoundaryConditions.FrontBound==3){ 	// periodic

		}
		else if 	((BoundaryConditions.FrontBound==5) || (BoundaryConditions.FrontBound==6) || (BoundaryConditions.FrontBound==7)){ // growth rate, thin-sheet
			ndof 	= 	1;
			SetDOF  =   1; 									// Vy
		}
		else if 	(BoundaryConditions.FrontBound==8){ 	// thin-sheet
			if(user->zdepth_BC_el == 0){
				user->zdepth_BC_el = nel_z;
			}
			if(iel_z>(nel_z-1-(user->zdepth_BC_el))){
				ndof 	= 2; 									// Vx,Vy
				//nn = round((nel_y-1.0)/4); mm = round((nel_y-1)/2);
				nn = 0;
				mm = (PetscInt)round(((PetscScalar)nel_x-1.0)/8.0);
				oo = (PetscInt)round(3.0*((PetscScalar)nel_x-1.0)/8.0);
				pp = (PetscInt)round(((PetscScalar)nel_x-1.0)/2.0);

				for(i=0; i<nnel_1D*nnel_1D; i++){
					if(iel_x>=nn && iel_x<mm){
						//BoundaryVelocity[i][1] = -BoundaryConditions.Eyy*user->W*(cos(2.0*M_PI*((coord_elem[nodes[i]].x)/user->W-(1.0/8.0))))*(cos(2.0*M_PI*((coord_elem[nodes[i]].x)/user->W-(1.0/8.0))));
						BoundaryVelocity[i][1] = -fac*BoundaryConditions.Eyy*user->W;
					}
					else if(iel_x<=oo && iel_x>=mm){
						BoundaryVelocity[i][1] = -fac*BoundaryConditions.Eyy*user->W;
					}
					else if(iel_x>oo && iel_x<pp){
						BoundaryVelocity[i][1] = -fac*BoundaryConditions.Eyy*user->W*(sin(2.0*M_PI*((coord_elem[nodes[i]].x)/user->W-(1.0/8.0))))*(sin(2.0*M_PI*((coord_elem[nodes[i]].x)/user->W-(1.0/8.0))));
					}
					else if(iel_x>=pp){
						BoundaryVelocity[i][1] = 0.0;
					}
				}
				for(i=0; i<nnel_1D*nnel_1D; i++){
					BoundaryVelocity[i][0] = 0.0;
				}
			}
			else{
				ndof 	= 2;
				for(i=0; i<nnel_1D*nnel_1D; i++){
					BoundaryVelocity[i][0] = 0.0;
				}
				for(i=0; i<nnel_1D*nnel_1D; i++){
					BoundaryVelocity[i][1] = 0.0;
				}
			}
		}
		else{
			SETERRQ1(PETSC_COMM_WORLD, 1,"Unknown front boundary condition BC.FrontBound=%lld",(LLD)(BoundaryConditions.FrontBound));
		}

		if (ndof>0){
			// Zero rows and collumns of local stiffness matrix
			Add_BC_ToLocalStiffness(C, VV, VP, PV, F_BC, P_BC, nodes, ndof, SetDOF, BoundaryVelocity);
		}

	}
	// -----------------------------------------------------------


	// --------------------- Back BC -----------------------------
	if (iel_y==nel_y-1){

		for(i=0; i<nnel_1D*nnel_1D; i++){ for(k=0; k<3; k++){
				BoundaryVelocity[i][k] = 0;		// Initialize all velocities to zero
		}}

		// Define back boundary nodes
		if( (element_type==DAVP_Q1P0) |  (element_type==DAVP_Q1Q1)  ){
		nodes[0] = 2; nodes[1] = 3; nodes[2] = 6; nodes[3] = 7;
		}
		else if( (element_type==DAVP_Q2PM1L) || (element_type==DAVP_Q2PM1G) ){
		nodes[0] = 3;  nodes[1] = 10; nodes[2] = 2; nodes[3] = 19;
		nodes[4] = 25; nodes[5] = 18; nodes[6] = 7; nodes[7] = 14;
		nodes[8] = 6;
		}

		ndof 	=	0;
		SetDOF  =   -1;
		if	((BoundaryConditions.BackBound==1) || (BoundaryConditions.BackBound==7)){ 			// free slip w. specified BG strainrate, growth rate
			ndof 	= 	1;
			SetDOF  =   1; 								// Vy

			// define velocity @ all y-nodes on back side of model, based on BG strainrate
			for(i=0; i<nnel_1D*nnel_1D; i++){
				BoundaryVelocity[i][0] = -fac*BoundaryConditions.Eyy*coord_elem[nodes[i]].y;
			}
		}
		else if 	(BoundaryConditions.BackBound==2){ 	// no slip
			ndof = 3; 									// Vx,Vy,Vz
		}
		else if 	(BoundaryConditions.BackBound==3){ 	// periodic

		}
		else if 	((BoundaryConditions.BackBound==5) || (BoundaryConditions.BackBound==6) || (BoundaryConditions.BackBound==8)){ // growth rate, thin-sheet
			ndof 	= 	1;
			SetDOF  =   1; 							 	// Vy
		}
		else{
			SETERRQ1(PETSC_COMM_WORLD, 1,"Unknown back boundary condition BC.BackBound=%lld",(LLD)(BoundaryConditions.BackBound));
		}

		if (ndof>0){
			// Zero rows and collumns of local stiffness matrix
			Add_BC_ToLocalStiffness(C, VV, VP, PV, F_BC, P_BC, nodes, ndof, SetDOF, BoundaryVelocity);
		}

	}
	// -----------------------------------------------------------


	// --------------------- Lower BC -----------------------------
	if (iel_z==0){

		for(i=0; i<nnel_1D*nnel_1D; i++){ for(k=0; k<3; k++){
				BoundaryVelocity[i][k] = 0;		// Initialize all velocities to zero
		}}

		// Define lower boundary nodes
		if( (element_type==DAVP_Q1P0) |  (element_type==DAVP_Q1Q1)  ){
		nodes[0] = 0; nodes[1] = 1; nodes[2] = 2; nodes[3] = 3;
		}
		else if( (element_type==DAVP_Q2PM1L) || (element_type==DAVP_Q2PM1G) ){
		nodes[0] = 0;  nodes[1] = 8;  nodes[2] = 1; nodes[3] = 11;
		nodes[4] = 20; nodes[5] = 9; nodes[6] = 3; nodes[7] = 10;
		nodes[8] = 2;
		}

		ndof 	=	0;
		SetDOF  =   -1;
		if	((BoundaryConditions.LowerBound==1) || (BoundaryConditions.LowerBound==6) || (BoundaryConditions.LowerBound==7)){ 			// free slip w. specified BG strainrate, thin-sheet, growth rate
			ndof 	=	 1;
			SetDOF  =    2;  								// Vz
			// define velocity @ all z-nodes on bottomx` side of model, based on BG strainrate
			for(i=0; i<nnel_1D*nnel_1D; i++){
				BoundaryVelocity[i][0] = fac*(BoundaryConditions.Eyy+BoundaryConditions.Exx)*coord_elem[nodes[i]].z;
			}

		}
		else if 	(BoundaryConditions.LowerBound==2){ 	// no slip
			ndof = 3; 										// Vx,Vy,Vz
		}
		else if 	(BoundaryConditions.LowerBound==8){ 	// no slip but with prescribed x and y velocity
			ndof = 3; 										// Vx,Vy,Vz
			for(i=0; i<nnel_1D*nnel_1D; i++){
				BoundaryVelocity[i][0] = -BoundaryConditions.Exx*coord_elem[nodes[i]].x;		// Vx velocity
				BoundaryVelocity[i][1] = -BoundaryConditions.Eyy*coord_elem[nodes[i]].y ;		// Vy velocity
			}


		}
		else if 	(BoundaryConditions.LowerBound==5){
			ndof 	=	 1;
			//SetDOF  =    1; 								// Vy
			//Add_BC_ToLocalStiffness(C, VV, VP, PV, F_BC, P_BC, nodes, ndof, SetDOF, BoundaryVelocity);
			SetDOF  =    2; 								// Vz
			Add_BC_ToLocalStiffness(C, VV, VP, PV, F_BC, P_BC, nodes, ndof, SetDOF, BoundaryVelocity);
		}
		else{
			SETERRQ1(PETSC_COMM_WORLD, 1,"Unknown lower boundary condition BC.LowerBound=%lld",(LLD)(BoundaryConditions.LowerBound));
		}

		if (ndof>0){
			// Zero rows and collumns of local stiffness matrix
			Add_BC_ToLocalStiffness(C, VV, VP, PV, F_BC, P_BC, nodes, ndof, SetDOF, BoundaryVelocity);
		}


	}
	// -----------------------------------------------------------


	// --------------------- Upper BC -----------------------------
	if (iel_z==nel_z-1){

		for(i=0; i<nnel_1D*nnel_1D; i++){ for(k=0; k<3; k++){
				BoundaryVelocity[i][k] = 0;		// Initialize all velocities to zero
		}}

		// Define upper boundary nodes
		if( (element_type==DAVP_Q1P0) |  (element_type==DAVP_Q1Q1)  ){
		nodes[0] = 4; nodes[1] = 5; nodes[2] = 6; nodes[3] = 7;
		}
		else if( (element_type==DAVP_Q2PM1L) || (element_type==DAVP_Q2PM1G) ){
		nodes[0] = 4;  nodes[1] = 12; nodes[2] = 5; nodes[3] = 15;
		nodes[4] = 21; nodes[5] = 13; nodes[6] = 7; nodes[7] = 14;
		nodes[8] = 6;
		}
		ndof 	=	0;
		SetDOF 	= 	-1;
		if	(BoundaryConditions.UpperBound==1){ 			// free slip w. specified BG strainrate
			ndof 	= 	1;
			SetDOF  =   2;  								// Vz

			// define velocity @ all z-nodes on top side of model, based on BG strainrate
			for(i=0; i<nnel_1D*nnel_1D; i++){
				BoundaryVelocity[i][0] = fac*(BoundaryConditions.Eyy+BoundaryConditions.Exx)*coord_elem[nodes[i]].z;
			}

		}
		else if 	(BoundaryConditions.UpperBound==2){ 	// no slip
			ndof = 3; 										// Vx,Vy,Vz
		}
		else if 	(BoundaryConditions.UpperBound==5){
			ndof 	= 	1;
			SetDOF  =    1; 								// Vy
		}
		else if 	(BoundaryConditions.UpperBound==0){
			// Free surface; nothing to be done
		}
		else{
			SETERRQ1(PETSC_COMM_WORLD, 1,"Unknown upper boundary condition BC.UpperBound=%lld",(LLD)(BoundaryConditions.UpperBound));
		}

		if (ndof>0){
			// Zero rows and collumns of local stiffness matrix
			Add_BC_ToLocalStiffness(C, VV, VP, PV, F_BC, P_BC, nodes, ndof, SetDOF, BoundaryVelocity);
		}

	}
	// -----------------------------------------------------------


	// --------------------- Internal BC -------------------------
	if (user->BC.InternalBound>0){
		if (iel_y==user->internalBC_frontel){

			//for(i=0; i<nnel_1D*nnel_1D; i++){ for(k=0; k<3; k++){
			//		BoundaryVelocity[i][k] = 0;		// Initialize all velocities to zero
			//}}

			// Define front boundary nodes
			if( (element_type==DAVP_Q1P0) |  (element_type==DAVP_Q1Q1)  ){
			nodes[0] = 0; nodes[1] = 1; nodes[2] = 4; nodes[3] = 5;
			}
			else if( (element_type==DAVP_Q2PM1L) || (element_type==DAVP_Q2PM1G) ){
			nodes[0] = 0;  nodes[1] = 8;  nodes[2] = 1; nodes[3] = 16;
			nodes[4] = 24; nodes[5] = 17; nodes[6] = 4; nodes[7] = 12;
			nodes[8] = 5;
			}

			ndof 	=	0;
			SetDOF  =  -1;
			if 	(BoundaryConditions.InternalBound==9){ 	// thin-sheet
				if(user->zdepth_BC_el == 0){
					user->zdepth_BC_el = nel_z;
				}
				if(iel_z>(nel_z-1-(user->zdepth_BC_el))){
					ndof 	= 2; 					    // Vx,Vy
					nn = 0;
					mm = (PetscInt)round((PetscScalar)nel_x/8.0);
					oo = (PetscInt)round((PetscScalar)nel_x/2.0);
					pp = (PetscInt)round((PetscScalar)nel_x/4.0) + (PetscInt)round((PetscScalar)nel_x/2.0);

					for(i=0; i<nnel_1D*nnel_1D; i++){
						if(iel_x>=nn && iel_x<oo){
							BoundaryVelocity[i][1] = -fac*BoundaryConditions.Eyy*user->W;
							//BoundaryVelocity[i][1] = -BoundaryConditions.Eyy*user->W*(cos(2.0*M_PI*((coord_elem[nodes[i]].x)/user->W-(1.0/8.0))))*(cos(2.0*M_PI*((coord_elem[nodes[i]].x)/user->W-(1.0/8.0))));
						}
						else if(iel_x>=oo && iel_x<pp){
							BoundaryVelocity[i][1] = -fac*BoundaryConditions.Eyy*user->W*(sin(2.0*M_PI*((coord_elem[nodes[i]].x)/user->W-(1.0/4.0))))*(sin(2.0*M_PI*((coord_elem[nodes[i]].x)/user->W-(1.0/4.0))));
						}
						else if(iel_x>=pp){
							BoundaryVelocity[i][1] = 0.0;
						}
					}
					for(i=0; i<nnel_1D*nnel_1D; i++){
						BoundaryVelocity[i][0] = 0.0;
					}
				}
				/*else{
					ndof 	= 2;
					for(i=0; i<nnel_1D*nnel_1D; i++){
						BoundaryVelocity[i][0] = 0.0;
					}
					for(i=0; i<nnel_1D*nnel_1D; i++){
						BoundaryVelocity[i][1] = 0.0;
					}
				}*/
			}
			else{
				SETERRQ1(PETSC_COMM_WORLD, 1,"Unknown front boundary condition BC.InternalBound=%lld",(LLD)(BoundaryConditions.InternalBound));
			}

			if (ndof>0){
				// Zero rows and collumns of local stiffness matrix
				Add_BC_ToLocalStiffness(C, VV, VP, PV, F_BC, P_BC, nodes, ndof, SetDOF, BoundaryVelocity);
			}

		}

		if (iel_y==user->internalBC_backel){

			//for(i=0; i<nnel_1D*nnel_1D; i++){ for(k=0; k<3; k++){
			//		BoundaryVelocity[i][k] = 0;		// Initialize all velocities to zero
			//}}

			// Define back boundary nodes
			if( (element_type==DAVP_Q1P0) |  (element_type==DAVP_Q1Q1)  ){
			nodes[0] = 2; nodes[1] = 3; nodes[2] = 6; nodes[3] = 7;
			}
			else if( (element_type==DAVP_Q2PM1L) || (element_type==DAVP_Q2PM1G) ){
			nodes[0] = 3;  nodes[1] = 10; nodes[2] = 2; nodes[3] = 19;
			nodes[4] = 25; nodes[5] = 18; nodes[6] = 7; nodes[7] = 14;
			nodes[8] = 6;
			}

			ndof 	=	0;
			SetDOF  =   -1;
			if 	(BoundaryConditions.InternalBound==9){ 	// thin-sheet
				if(user->zdepth_BC_el == 0){
					user->zdepth_BC_el = nel_z;
				}
				if(iel_z>(nel_z-1-(user->zdepth_BC_el))){
					ndof 	= 2; 					    // Vx,Vy
					nn = 0;
					mm = (PetscInt)round((PetscScalar)nel_x/8.0);
					oo = (PetscInt)round((PetscScalar)nel_x/2.0);
					pp = (PetscInt)round((PetscScalar)nel_x/4.0) + (PetscInt)round((PetscScalar)nel_x/2.0);

					for(i=0; i<nnel_1D*nnel_1D; i++){
						if(iel_x>=nn && iel_x<oo){
							BoundaryVelocity[i][1] = -fac*BoundaryConditions.Eyy*user->W;
							//BoundaryVelocity[i][1] = -BoundaryConditions.Eyy*user->W*(cos(2.0*M_PI*((coord_elem[nodes[i]].x)/user->W-(1.0/8.0))))*(cos(2.0*M_PI*((coord_elem[nodes[i]].x)/user->W-(1.0/8.0))));
						}
						else if(iel_x>=oo && iel_x<pp){
							BoundaryVelocity[i][1] = -fac*BoundaryConditions.Eyy*user->W*(sin(2.0*M_PI*((coord_elem[nodes[i]].x)/user->W-(1.0/4.0))))*(sin(2.0*M_PI*((coord_elem[nodes[i]].x)/user->W-(1.0/4.0))));
						}
						else if(iel_x>=pp){
							BoundaryVelocity[i][1] = 0.0;
						}
					}
					for(i=0; i<nnel_1D*nnel_1D; i++){
						BoundaryVelocity[i][0] = 0.0;
					}
				}
				/*else{
					ndof 	= 2;
					for(i=0; i<nnel_1D*nnel_1D; i++){
						BoundaryVelocity[i][0] = 0.0;
					}
					for(i=0; i<nnel_1D*nnel_1D; i++){
						BoundaryVelocity[i][1] = 0.0;
					}
				}*/
			}
			else{
				SETERRQ1(PETSC_COMM_WORLD, 1,"Unknown back boundary condition BC.InternalBound=%lld",(LLD)(BoundaryConditions.InternalBound));
			}

			if (ndof>0){
				// Zero rows and collumns of local stiffness matrix
				Add_BC_ToLocalStiffness(C, VV, VP, PV, F_BC, P_BC, nodes, ndof, SetDOF, BoundaryVelocity);
			}
		}
	}
	// -----------------------------------------------------------


	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Help function to zero rows and collumns in element stiffness matrix [for setting boundary conditions]*/
#undef __FUNCT__
#define __FUNCT__ "Add_BC_ToLocalStiffness"
PetscErrorCode Add_BC_ToLocalStiffness(LaMEMVelPressureDA C, double VV[MAX_edof][MAX_edof],  double VP[MAX_edof][MAX_npres],
		double PV[MAX_npres][MAX_edof], double F_BC[], double P_BC[],
		PetscInt nodes[], PetscInt ndof, PetscInt SetDOF,
		double BoundaryVelocity[MAX_nnel_1D*MAX_nnel_1D][3] )
{
	PetscInt 	idof, dof, kk;
	PetscInt 		i;
	PetscInt npres, edof, nnel_1D;


	npres   = C->npres;
	edof    = C->edof;
	nnel_1D = C->nnel_1D;

	/* At an element level, we have the following ordering
	 *
	 * | VV VP | |V|  |F_BC|
	 * | PV PP | |P|  |P_BC|
	 *
	 */



	for (idof=0; idof<ndof; idof++){
		for (i=0; i<nnel_1D*nnel_1D; i++){
			if (ndof>1){
				dof = 3*nodes[i] + idof;
			}
			else{
				dof = 3*nodes[i] + SetDOF;
			}

			// zero the row
			for (kk=0; kk<edof; kk++){
				VV[dof][kk] = 0.0;
			}

			// Zero collumns of VV matrix
			for (kk=0; kk<edof; kk++){

				// Add contribution to RHS. Also the real boundary value is added here, to be
				// consistent with the manner in which diagonal terms are added on the diagonal of the
				// global stiffness matrix (it might end up as an 8 on the diagonal instead of a 1, as
				// you might expect).
				if (kk != dof){
					F_BC[kk] =  F_BC[kk] - BoundaryVelocity[i][idof]*VV[kk][dof] ;
				}
				else{
					F_BC[kk] =  BoundaryVelocity[i][idof];
				}

				// Zero collumn
				VV[kk][dof] = 0.0;
			}

			// Zero row of VP matrix
			for (kk=0; kk<npres; kk++){
				VP[dof][kk] = 0.0;
			}

			// Zero collumn of PV matrix & push values to RHS
			for (kk=0; kk<npres; kk++){
				P_BC[kk] 	=   P_BC[kk] - BoundaryVelocity[i][idof]*PV[kk][dof];
				PV[kk][dof] = 	0.0;
			}



			// Add 1 in diagonal
			VV[dof][dof] = 1.0;

		}
	}

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Set's the boundary conditions in the RHS vector */
#undef __FUNCT__
#define __FUNCT__ "SetBoundaryConditionsRHS"
PetscErrorCode SetBoundaryConditionsRHS( DM da, UserContext *user, Vec b, PetscInt SetNonZeroValuesFlag)
{
	PetscErrorCode 		ierr;
	PetscInt			i,j,k,xs,xm,ys,ym,zs,zm, mx,my,mz,maxi,maxk,nn,mm, oo, pp;
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
#if !defined(SYMMETRIZE_MATRIX_AT_ELEMENTLEVEL)
	factor = 1.0;
#endif
#if defined(SYMMETRIZE_MATRIX_AT_ELEMENTLEVEL)
	factor = 0.0;
#endif



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
					else if (user->BC.LeftBound==5){ // growth rate
						if(SetNonZeroValuesFlag==1){ rhs[k][j][i].Vx = 0.0; rhs[k][j][i].Vy = 0.0;}
						else{ rhs[k][j][i].Vx = 0.0; rhs[k][j][i].Vy = 0.0;}
					}
					else if ((user->BC.LeftBound==6) || (user->BC.LeftBound==7)){ // growth rate, thin-sheet
						if(SetNonZeroValuesFlag==1){ rhs[k][j][i].Vx = 0.0; }
						else{ rhs[k][j][i].Vx = 0.0; }
					}
					else if (user->BC.LeftBound==3){
						// periodic - taken care of already
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
					else if (user->BC.RightBound==5){ // growth rate
						if(SetNonZeroValuesFlag==1){ rhs[k][j][i].Vx = (-user->BC.Exx)*coords[k  ][j  ][i  ].x*factor; rhs[k][j][i].Vy = 0.0;}
						else{ rhs[k][j][i].Vx = 0.0; rhs[k][j][i].Vy = 0.0;}
					}
					else if (user->BC.RightBound==6){ // thin-sheet
						if(SetNonZeroValuesFlag==1){ rhs[k][j][i].Vy = 0.0;}
						else{ rhs[k][j][i].Vy = 0.0;}
						nn = (PetscInt)round((PetscScalar)my/4.0); //nn = 3;
						mm = (PetscInt)round((PetscScalar)my/2.0); //mm = 5;
						if(j<nn){
							if(SetNonZeroValuesFlag==1){ rhs[k][j][i].Vx = -user->BC.Exx*factor;}
							else{ rhs[k][j][i].Vx = 0.0;}
						}
						else if(j>=nn && j<=mm){
							if(SetNonZeroValuesFlag==1){ rhs[k][j][i].Vx = -user->BC.Exx*(cos(2.0*M_PI*((coords[k][j][i].y)/user->L-(1.0/4.0))))*(cos(2.0*M_PI*((coords[k][j][i].y)/user->L-(1.0/4.0))))*factor;}
							else{ rhs[k][j][i].Vx = 0.0;}
						}
						else if(j>mm){
							if(SetNonZeroValuesFlag==1){ rhs[k][j][i].Vx = 0.0;}
							else{ rhs[k][j][i].Vx = 0.0;}
						}
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
					else if (user->BC.FrontBound==3){ // free slip w. specified BG strainrate and shear velocity
						if (SetNonZeroValuesFlag==1){
							rhs[k][j][i].Vy  = -user->BC.Eyy*coords[k  ][j  ][i  ].y*factor;
							rhs[k][j][i].Vx  = user->Vx_Front*factor;;
						}
						else{
							rhs[k][j][i].Vy  = 0.0;
							rhs[k][j][i].Vx  = 0.0;
						}

					}
					else if (user->BC.FrontBound==4){ // free slip w. specified BG strainrate and shear velocity over part.area
						if (SetNonZeroValuesFlag==1){
							maxi = (PetscInt)(user->Vy_Partx*xm);
							maxk = (PetscInt)(user->Vy_Partz*zm);
							if ((i<maxi+1) && (k>(maxk-1))){
								rhs[k][j][i].Vy  = user->Vy_Front*factor;
							}
							else{
								rhs[k][j][i].Vy  = 0.0;
							}
						}
						else{
							rhs[k][j][i].Vy  = 0.0;
						}
						rhs[k][j][i].Vz  = 0.0;		// vertical motion contraint

					}
					else if ((user->BC.FrontBound==5) || (user->BC.FrontBound==6) || (user->BC.FrontBound==7)){ // growth rate, thin-sheet
						if(SetNonZeroValuesFlag==1){ rhs[k][j][i].Vy = 0.0;}
						else{ rhs[k][j][i].Vy = 0.0;}
					}
					else if (user->BC.FrontBound==8){
						if (SetNonZeroValuesFlag==1){

							if(user->zdepth_BC_node == 0){
								user->zdepth_BC_node = mz-1;
							}
							if(k>=(mz-1-(user->zdepth_BC_node))){
								//nn = round((nel_y-1.0)/4); mm = round((nel_y-1)/2);
								nn = 0;
								mm = (PetscInt)round((PetscScalar)mx/8.0);
								oo = (PetscInt)round(3.0*(PetscScalar)mx/8.0);
								pp = (PetscInt)round((PetscScalar)mx/2.0);

								if(mx>=nn && mx<mm){
									//BoundaryVelocity[i][1] = -BoundaryConditions.Eyy*user->W*(cos(2.0*M_PI*((coord_elem[nodes[i]].x)/user->W-(1.0/8.0))))*(cos(2.0*M_PI*((coord_elem[nodes[i]].x)/user->W-(1.0/8.0))));
									rhs[k][j][i].Vy = -user->BC.Eyy*user->W*factor;
								}
								else if(mx<=oo && mx>=mm){
									rhs[k][j][i].Vy = -user->BC.Eyy*user->W*factor;
								}
								else if(mx>oo && mx<pp){
									rhs[k][j][i].Vy = -user->BC.Eyy*user->W*(sin(2.0*M_PI*((coords[k][j][i].x)/user->W-(1.0/8.0))))*(sin(2.0*M_PI*((coords[k][j][i].x)/user->W-(1.0/8.0))))*factor;
								}
								else if(mx>=pp){
									rhs[k][j][i].Vy = 0.0;
								}
								rhs[k][j][i].Vx = 0.0;
							}
							else{
								rhs[k][j][i].Vx = 0.0;
								rhs[k][j][i].Vy = 0.0;
							}

						}
						else{
							rhs[k][j][i].Vx = 0.0;
							rhs[k][j][i].Vy = 0.0;
						}

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
					else if (user->BC.BackBound==3){ // free slip w. specified BG strainrate & shear velocity
						if (SetNonZeroValuesFlag==1){
							rhs[k][j][i].Vy  = -user->BC.Eyy*coords[k  ][j  ][i  ].y*factor;
							rhs[k][j][i].Vx  = user->Vx_Back*factor;
						}
						else {
							rhs[k][j][i].Vy  = 0.0;
							rhs[k][j][i].Vx  = 0.0;
						}

					}
					else if ((user->BC.BackBound==5) || (user->BC.BackBound==6)){ // growth rate, thin-sheet
						if(SetNonZeroValuesFlag==1){ rhs[k][j][i].Vy = 0.0;}
						else{ rhs[k][j][i].Vy = 0.0;}
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
					else if (user->BC.LowerBound==8){ // no slip with prescribed velocity
						rhs[k][j][i].Vx  = -(user->BC.Exx)*coords[k  ][j  ][i  ].x*factor;
						rhs[k][j][i].Vy  = -(user->BC.Eyy)*coords[k  ][j  ][i  ].y*factor;
						rhs[k][j][i].Vz  = 0.0;
					}
					else if (user->BC.LowerBound==5){ // growth rate
						if(SetNonZeroValuesFlag==1){
							rhs[k][j][i].Vz = 0.0;
							rhs[k][j][i].Vy = 0.0;
						}
						else{
							rhs[k][j][i].Vz = 0.0;
							rhs[k][j][i].Vy = 0.0;
						}
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
					else if (user->BC.UpperBound==5){ // growth rate
						if(SetNonZeroValuesFlag==1){ rhs[k][j][i].Vy = 0.0;}
						else{ rhs[k][j][i].Vy = 0.0;}
					}
					else if (user->BC.UpperBound==0){	// free surface

					}
					else{
						SETERRQ1(PETSC_COMM_WORLD, 1,"Unknown upper boundary condition BC.UpperBound=%lld",(LLD)(user->BC.UpperBound));
					}

				}

				if (user->BC.InternalBound>0){		// If we selected internal boundary conditions as an option, set them here.
					if (j==user->internalBC_node){					/* Internal boundary */
						if (user->BC.InternalBound==9){
							if (SetNonZeroValuesFlag==1){
								if(user->zdepth_BC_node == 0){
									user->zdepth_BC_node = mz-1;
								}
								if(k>=(mz-1-(user->zdepth_BC_node))){
									nn = 0;
									mm = (PetscInt)round(((PetscScalar)mx-1.0)/8.0);
									oo = (PetscInt)round(((PetscScalar)mx-1.0)/2.0);
									pp = (PetscInt)round(((PetscScalar)mx-1.0)/4.0) + (PetscInt)round(((PetscScalar)mx-1.0)/2.0);

									if(mx>=nn && mx<oo){
										//rhs[k][j][i].Vy = -user->BC.Eyy*user->W*(sin(2.0*M_PI*((coords[k][j][i].x)/user->W-(1.0/8.0))))*(sin(2.0*M_PI*((coords[k][j][i].x)/user->W-(1.0/8.0))))*factor;
										rhs[k][j][i].Vy = -user->BC.Eyy*user->W*factor;
									}
									else if(mx>=oo && mx<pp){
										rhs[k][j][i].Vy = -user->BC.Eyy*user->W*(sin(2.0*M_PI*((coords[k][j][i].x)/user->W-(1.0/4.0))))*(sin(2.0*M_PI*((coords[k][j][i].x)/user->W-(1.0/4.0))))*factor;
									}
									else if(mx>=pp){
										rhs[k][j][i].Vy = 0.0;
									}
									rhs[k][j][i].Vx = 0.0;
								}
								/*else{
									rhs[k][j][i].Vx = 0.0;
									rhs[k][j][i].Vy = 0.0;
								}*/

							}
							else{
								rhs[k][j][i].Vx = 0.0;
								rhs[k][j][i].Vy = 0.0;
							}

						}
						else{
							SETERRQ1(PETSC_COMM_WORLD, 1,"Unknown internal boundary condition BC.InternalBound=%lld",(LLD)(user->BC.InternalBound));
						}
					}

				}

			}
		}
	}
	ierr = DMDAVecRestoreArray(da ,b,	   &rhs);	CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(cda,gc,&coords);	CHKERRQ(ierr);

	//DMDestroy( cda );
	//VecDestroy( gc );


	PetscFunctionReturn(0);
}
/*==========================================================================================================*/



/*==========================================================================================================*/
/* Computes the element stiffness matrixes VV, PP, VP and PV, assuming
 * 	(1) a  constant viscosity for the element.
 * 	(2) an isotropic Newtonian viscous rheology.
 *
 * The routine performs numerical integration.
 */
#undef __FUNCT__
#define __FUNCT__ "ComputeElementStiffnessMatrixes"

PetscErrorCode ComputeElementStiffnessMatrixes( LaMEMVelPressureDA C, double VV[MAX_edof][MAX_edof], double VP[MAX_edof][MAX_npres],
		double PV[MAX_npres][MAX_edof], 		double PP[MAX_npres][MAX_npres], double PP_prec[MAX_npres][MAX_npres],
		double IntPoint[MAX_ndim][MAX_ngp_vel],	double IntWeight[],	double *mu,
		DMDACoor3d coord_elem[], PetscScalar *minDiagonalRatio)
{
	//PetscErrorCode 		ierr;
	PetscInt			ii, jj, kk, intp, iiii;
	PetscScalar 		Point[3], ShapeVel[MAX_nnel], **dhdsVel,ShapeP[MAX_npres], DiagonalLength_1, DiagonalLength_2, DiagonalLength_3, maxDiagonalLength, minDiagonalLength;
	PetscScalar	 		**Jacob, **InvJacob, DetJacob;
	PetscScalar 		KinMtx[6][MAX_edof], **dhdPhys, VV_add[MAX_edof][MAX_edof];
	PetscScalar			MatMtx[6][6], DiagonalRatio;
	DMDACoor3d 			CoordIntp;
	PetscInt 				npres, edof, nnel, ngp_vel;
	DAVPElementType 	element_type;


	nnel    = C->nnel;
	npres   = C->npres;
	edof    = C->edof;
	ngp_vel = C->ngp_vel;
	element_type = C->type;



	LaMEMCreate2dArray( 3, nnel, &dhdsVel, PETSC_NULL );
	LaMEMCreate2dArray( 3, nnel, &dhdPhys, PETSC_NULL );
	LaMEMCreate2dArray( 3, 3, &Jacob, PETSC_NULL );
	LaMEMCreate2dArray( 3, 3, &InvJacob, PETSC_NULL );


	/* Initialization */
	for(ii=0; ii<edof; ii++){ for(jj=0; jj<edof;  jj++){	 VV[ii][jj]=0.0;	}}
	for(ii=0; ii<edof; ii++){ for(jj=0; jj<npres; jj++){	 VP[ii][jj]=0.0;	}}
	for(ii=0; ii<npres;ii++){ for(jj=0; jj<edof;  jj++){	 PV[ii][jj]=0.0;	}}
	for(ii=0; ii<npres;ii++){ for(jj=0; jj<npres; jj++){	 PP[ii][jj]=0.0;	}}
	for(ii=0; ii<npres;ii++){ for(jj=0; jj<npres; jj++){PP_prec[ii][jj]=0.0;	}}

	//for(ii=0; ii<edof; ii++){ 								V_RHS[ii]=0.0;	}

	/* Loop over integration points */
	for (intp=0; intp<ngp_vel; intp++){


		if (PetscAbsScalar(mu[0])<1e-8){
			PetscPrintf(PETSC_COMM_SELF,"*** Emergency: abs(mu)<1e-8 in forming stiffness matrix! *** \n");
			//MPI_Abort(PETSC_COMM_WORLD,1);
		}
		if (isnan(mu[0])){
			PetscPrintf(PETSC_COMM_SELF,"*** Emergency: isnan(mu) in forming stiffness matrix! *** \n");
			//MPI_Abort(PETSC_COMM_WORLD,1);
		}

		Point[0]= IntPoint[0][intp];Point[1]= IntPoint[1][intp];Point[2]= IntPoint[2][intp];


		ComputeShapeFunctionVelocity(ShapeVel, dhdsVel, Point);				 		// Velocity shape function
		ComputeCoordIntp(nnel,ShapeVel, coord_elem,  &CoordIntp);				 	// Real coordinate of integration point

		ComputeShapeFunctionPressure(ShapeP, 	CoordIntp, Point);		     		// Pressure shape function

		ComputeJacobianElement(nnel,dhdsVel,coord_elem,Jacob,InvJacob,&DetJacob);  	// Jacobian etc. of element
		ComputeKinmtx(nnel,edof,KinMtx, dhdPhys, dhdsVel, InvJacob);				// Kinematic matrix

		ComputeMaterialMatrix(mu[intp], MatMtx);									 // Material matrix




		// Catch errors (for debugging only)
		if (DetJacob<0){
			PetscPrintf(PETSC_COMM_WORLD,"In routine ComputeElementStiffnessMatrixes: Negative Jacobian, DetJacob=%g  \n",DetJacob);
			if( (element_type==DAVP_Q1P0) |  (element_type==DAVP_Q1Q1)  ){
				PetscPrintf(PETSC_COMM_WORLD," x-coords of local element: [%g,%g,%g,%g,%g,%g,%g,%g] \n",coord_elem[0].x,coord_elem[1].x,coord_elem[2].x,coord_elem[3].x,coord_elem[4].x,coord_elem[5].x,coord_elem[6].x,coord_elem[7].x);
				PetscPrintf(PETSC_COMM_WORLD," y-coords of local element: [%g,%g,%g,%g,%g,%g,%g,%g] \n",coord_elem[0].y,coord_elem[1].y,coord_elem[2].y,coord_elem[3].y,coord_elem[4].y,coord_elem[5].y,coord_elem[6].y,coord_elem[7].y);
				PetscPrintf(PETSC_COMM_WORLD," z-coords of local element: [%g,%g,%g,%g,%g,%g,%g,%g] \n",coord_elem[0].z,coord_elem[1].z,coord_elem[2].z,coord_elem[3].z,coord_elem[4].z,coord_elem[5].z,coord_elem[6].z,coord_elem[7].z);
			}
			else if( (element_type==DAVP_Q2PM1L) || (element_type==DAVP_Q2PM1G) ){
				PetscPrintf(PETSC_COMM_WORLD," x-coords of local element: [%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g] \n",coord_elem[0].x,coord_elem[1].x,coord_elem[2].x,coord_elem[3].x,coord_elem[4].x,coord_elem[5].x,coord_elem[6].x,coord_elem[7].x,
						coord_elem[8].x,coord_elem[9].x,coord_elem[10].x,coord_elem[11].x,coord_elem[12].x,coord_elem[13].x,coord_elem[14].x,coord_elem[15].x);

			}

			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP,"Jacobian is negative! Mesh is too deformed!");

		}

		// Compute VV matrix
		Multiply_Kin_Times_Mat_Times_Kin( edof, VV_add, KinMtx, MatMtx);			 // VV_add = Kinmtx*MatMtx*Kinmtx
		for(ii=0; ii<edof; ii++){for (jj=0; jj<edof; jj++) {
			VV[ii][jj]=VV[ii][jj] + VV_add[ii][jj]*DetJacob*IntWeight[intp];
		}}

		// Compute VP matrix
		for (ii=0; ii<nnel;ii++){for (jj=0; jj<npres; jj++){for (kk=0; kk<3; kk++){
			VP[3*ii+kk][jj]	=	VP[3*ii+kk][jj]	-	dhdPhys[kk][ii]*ShapeP[jj]*DetJacob*IntWeight[intp];
		}}}

		// Compute PV matrix
		for (ii=0; ii<nnel; ii++){for (jj=0; jj<npres; jj++){for (kk=0; kk<3; kk++){
			PV[jj][3*ii+kk]	=	PV[jj][3*ii+kk]	-	dhdPhys[kk][ii]*ShapeP[jj]*DetJacob*IntWeight[intp];
		}}}

		// Compute PP and PP_prec matrix
		{
			PetscScalar shift;

			if (element_type==DAVP_Q1Q1){	shift = 0.0156;	}  // 1/16 in 2D, 1/64 in 3D
			else{							shift = 0.0;	}

			for (ii=0; ii<npres; ii++){for (jj=0; jj<npres; jj++){
				PP[ii][jj]		=	PP[ii][jj]		+	DetJacob*IntWeight[intp]*(ShapeP[ii]*ShapeP[jj] - shift);
			}}

			for (ii=0; ii<npres; ii++){for (jj=0; jj<npres; jj++){
				PP_prec[ii][jj]	=	PP_prec[ii][jj]	+	DetJacob*IntWeight[intp]*ShapeP[ii]*ShapeP[jj];
			}}

		}
		// Compute contribution of gravity to the RHS
		//for (ii=0; ii<nnel; ii++){
		//	V_RHS[3*ii+2] = V_RHS[3*ii+2]	-	rho_g*ShapeVel[ii]*DetJacob*IntWeight[intp];
		//}

		if (1==0){	//  debugging part
			PetscPrintf(PETSC_COMM_WORLD,"Integration point = [%g,%g,%g]\n", Point[0] , Point[1], Point[2]);
			PetscPrintf(PETSC_COMM_WORLD,"	CoordIntp = [%g,%g,%g]\n",	CoordIntp.x , 	CoordIntp.y, 	CoordIntp.z);


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


			iiii = 0;
			PetscPrintf(PETSC_COMM_WORLD,"Kinmtx[:][%lld]= %g %g %g %g %g %g \n",(LLD)iiii,KinMtx[0][iiii], KinMtx[1][iiii],KinMtx[2][iiii], KinMtx[3][iiii],KinMtx[4][iiii], KinMtx[5][iiii]);

			iiii=0;
			PetscPrintf(PETSC_COMM_WORLD,"PV[%lld][:]= %g %g %g %g %g %g  \n", (LLD)iiii,PV[0][1],PV[0][2], PV[0][3], PV[0][4], PV[0][5], PV[0][6]);

			if (npres>1){
				PetscPrintf(PETSC_COMM_WORLD,"ShapeP= %g %g %g %g   \n", ShapeP[0], ShapeP[1], ShapeP[2], ShapeP[3], ShapeP[4]);
			}

		//	MPI_Abort(PETSC_COMM_WORLD,1);
		}	// end of debugging part

	}

	LaMEMDestroy2dArray(&dhdsVel, PETSC_NULL );
	LaMEMDestroy2dArray(&dhdPhys, PETSC_NULL );
	LaMEMDestroy2dArray(&Jacob, PETSC_NULL );
	LaMEMDestroy2dArray(&InvJacob, PETSC_NULL );

	/* Perform some element metrics, to see if the element becomes too distorted */
	DiagonalLength_1 = (coord_elem[6].x-coord_elem[0].x)*(coord_elem[6].x-coord_elem[0].x) + (coord_elem[6].y-coord_elem[0].y)*(coord_elem[6].y-coord_elem[0].y) + (coord_elem[6].z-coord_elem[0].z)*(coord_elem[6].z-coord_elem[0].z);
	DiagonalLength_2 = (coord_elem[7].x-coord_elem[1].x)*(coord_elem[7].x-coord_elem[1].x) + (coord_elem[7].y-coord_elem[1].y)*(coord_elem[7].y-coord_elem[1].y) + (coord_elem[7].z-coord_elem[1].z)*(coord_elem[7].z-coord_elem[1].z);
	DiagonalLength_3 = (coord_elem[4].x-coord_elem[2].x)*(coord_elem[4].x-coord_elem[2].x) + (coord_elem[4].y-coord_elem[2].y)*(coord_elem[4].y-coord_elem[2].y) + (coord_elem[4].z-coord_elem[2].z)*(coord_elem[4].z-coord_elem[2].z);

	maxDiagonalLength = DiagonalLength_1;
	if (DiagonalLength_2>maxDiagonalLength){ maxDiagonalLength = DiagonalLength_2;}
	if (DiagonalLength_3>maxDiagonalLength){ maxDiagonalLength = DiagonalLength_3;}

	minDiagonalLength = DiagonalLength_1;
	if (DiagonalLength_2<minDiagonalLength){ minDiagonalLength = DiagonalLength_2;}
	if (DiagonalLength_3<minDiagonalLength){ minDiagonalLength = DiagonalLength_3;}
	DiagonalRatio = minDiagonalLength/maxDiagonalLength;



	if (DiagonalRatio<*minDiagonalRatio){
		*minDiagonalRatio = DiagonalRatio; 			 // Store the minimum jacobian
	}

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/



/*==========================================================================================================*/
/* Writes an element stiffness matrix to disk (for debugging purposes)
 */
#undef __FUNCT__
#define __FUNCT__ "WriteElementStiffnessMatrixToDisk"
PetscErrorCode WriteElementStiffnessMatrixToDisk( LaMEMVelPressureDA C, double VV[MAX_edof][MAX_edof], 			double VP[MAX_edof][MAX_npres],
		double PV[MAX_npres][MAX_edof], 					double V_RHS[], 			double PP[MAX_npres][MAX_npres],
		double mu_average,
		double rho_g, 								DMDACoor3d coord_elem[], 		PetscInt number,
		PetscInt ielx, 								PetscInt iely, 					PetscInt ielz)
{

	PetscErrorCode ierr;
	Mat 			PV_MAT, VP_MAT, VV_MAT, PP_MAT, COORD_MAT, INFO_MAT, VRHS_MAT;
	PetscInt 		i,j;
	char			SaveFileName[PETSC_MAX_PATH_LEN];
	PetscViewer		view_out;
	PetscInt npres, edof, nnel;


	nnel    = C->nnel;
	npres   = C->npres;
	edof    = C->edof;

	/* ----------------------------------------------------
	 * Transform all vectors & matrixes into PETSC format
	 * ----------------------------------------------------*/

	// PV
	ierr = MatCreate(PETSC_COMM_WORLD,&PV_MAT); CHKERRQ(ierr);
	ierr = MatSetSizes(PV_MAT,PETSC_DECIDE,PETSC_DECIDE,npres,edof); CHKERRQ(ierr);
	ierr = MatSetType(PV_MAT,MATSEQDENSE); CHKERRQ(ierr);
	for (i=0; i<npres; i++){
		for (j=0; j<edof; j++){
			 ierr = MatSetValue(PV_MAT,i,j,PV[i][j],INSERT_VALUES); CHKERRQ(ierr);
		}
	}
	ierr = MatAssemblyBegin(PV_MAT,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(PV_MAT,	MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	// VP
	ierr = MatCreate(PETSC_COMM_WORLD,&VP_MAT); CHKERRQ(ierr);
	ierr = MatSetSizes(VP_MAT,PETSC_DECIDE,PETSC_DECIDE,edof,npres); CHKERRQ(ierr);
	ierr = MatSetType(VP_MAT,MATSEQDENSE); CHKERRQ(ierr);
	for (i=0; i<edof; i++){
		for (j=0; j<npres; j++){
			ierr =  MatSetValue(VP_MAT,i,j,VP[i][j],INSERT_VALUES); CHKERRQ(ierr);
		}
	}
	ierr = MatAssemblyBegin(VP_MAT,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(VP_MAT,	MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	// VV
	ierr = MatCreate(PETSC_COMM_WORLD,&VV_MAT); CHKERRQ(ierr);
	ierr = MatSetSizes(VV_MAT,PETSC_DECIDE,PETSC_DECIDE,edof,edof); CHKERRQ(ierr);
	ierr = MatSetType(VV_MAT,MATSEQDENSE); CHKERRQ(ierr);
	for (i=0; i<edof; i++){
		for (j=0; j<edof; j++){
			 ierr = MatSetValue(VV_MAT,i,j,VV[i][j],INSERT_VALUES); CHKERRQ(ierr);
		}
	}
	ierr = MatAssemblyBegin(VV_MAT,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(VV_MAT,	MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	// V_RHS
	ierr = MatCreate(PETSC_COMM_WORLD,&VRHS_MAT); CHKERRQ(ierr);
	ierr = MatSetSizes(VRHS_MAT,PETSC_DECIDE,PETSC_DECIDE,edof,1); CHKERRQ(ierr);
	ierr = MatSetType(VRHS_MAT,MATSEQDENSE); CHKERRQ(ierr);
	for (i=0; i<edof; i++){
		ierr = MatSetValue(VRHS_MAT,i,0,V_RHS[i],INSERT_VALUES); CHKERRQ(ierr);
	}
	ierr = MatAssemblyBegin(VRHS_MAT,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(VRHS_MAT,	MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);


	// PP
	ierr = MatCreate(PETSC_COMM_WORLD,&PP_MAT); CHKERRQ(ierr);
	ierr = MatSetSizes(PP_MAT,PETSC_DECIDE,PETSC_DECIDE,npres,npres); CHKERRQ(ierr);
	ierr = MatSetType(PP_MAT,MATSEQDENSE); CHKERRQ(ierr);
	for (i=0; i<npres; i++){
		for (j=0; j<npres; j++){
			ierr =  MatSetValue(PP_MAT,i,j,PP[i][j],INSERT_VALUES); CHKERRQ(ierr);
		}
	}
	ierr = MatAssemblyBegin(PP_MAT,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(PP_MAT,	MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	// COORD_MAT
	ierr = MatCreate(PETSC_COMM_WORLD,&COORD_MAT); CHKERRQ(ierr);
	ierr = MatSetSizes(COORD_MAT,PETSC_DECIDE,PETSC_DECIDE,3,nnel); CHKERRQ(ierr);
	ierr = MatSetType(COORD_MAT,MATSEQDENSE); CHKERRQ(ierr);
	for (j=0; j<nnel; j++){
		ierr =  MatSetValue(COORD_MAT,0,j,coord_elem[j].x,INSERT_VALUES); CHKERRQ(ierr);
		 ierr = MatSetValue(COORD_MAT,1,j,coord_elem[j].y,INSERT_VALUES); CHKERRQ(ierr);
		 ierr = MatSetValue(COORD_MAT,2,j,coord_elem[j].z,INSERT_VALUES); CHKERRQ(ierr);
	}
	ierr = MatAssemblyBegin(COORD_MAT,	MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(COORD_MAT,	MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	// INFO_MAT
	ierr = MatCreate(PETSC_COMM_WORLD,&INFO_MAT); CHKERRQ(ierr);
	ierr = MatSetSizes(INFO_MAT,PETSC_DECIDE,PETSC_DECIDE,5,1); CHKERRQ(ierr);
	ierr = MatSetType(INFO_MAT,MATSEQDENSE); CHKERRQ(ierr);
	ierr = MatSetValue(INFO_MAT,0,0, mu_average,INSERT_VALUES); CHKERRQ(ierr);
	ierr = MatSetValue(INFO_MAT,1,0, rho_g,INSERT_VALUES); CHKERRQ(ierr);
	ierr = MatSetValue(INFO_MAT,2,0,((double) ielx),INSERT_VALUES); CHKERRQ(ierr);
	ierr = MatSetValue(INFO_MAT,3,0,((double) iely),INSERT_VALUES); CHKERRQ(ierr);
	ierr = MatSetValue(INFO_MAT,4,0,((double) ielz),INSERT_VALUES); CHKERRQ(ierr);
	ierr = MatAssemblyBegin(INFO_MAT,	MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(INFO_MAT,	MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);


	/* Write data to disk */
	sprintf(SaveFileName,"ElementStiffnessMatrix%lld.out",(LLD)number);  // construct the filename

//	ierr = PetscViewerBinaryMatlabOpen(PETSC_COMM_SELF,SaveFileName,&view_out); CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,SaveFileName,FILE_MODE_WRITE,&view_out); CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(view_out, PETSC_VIEWER_BINARY_MATLAB); CHKERRQ(ierr);

	ierr = MatView(COORD_MAT, 				view_out); CHKERRQ(ierr);
	ierr = MatView(VV_MAT, 				view_out); CHKERRQ(ierr);
	ierr = MatView(VP_MAT, 				view_out); CHKERRQ(ierr);
	ierr = MatView(PV_MAT, 				view_out); CHKERRQ(ierr);
	ierr = MatView(PP_MAT, 				view_out); CHKERRQ(ierr);
	ierr = MatView(INFO_MAT, 				view_out); CHKERRQ(ierr);
	ierr = MatView(VRHS_MAT, 				view_out); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&view_out); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_SELF," Saved output in file %s \n",SaveFileName); CHKERRQ(ierr);


	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*
Use the DA's to generate matrices and perform efficient memory allocation for A11 and A22.
For the mixed terms, A12, we use a generic routine.
For A21, we could simply call MatTranspose, but its faster to re-use the same generic routine
*/

// this function is not used
/*
#undef __FUNCT__
#define __FUNCT__ "GenerateStokesOperators"
PetscErrorCode GenerateStokesOperators( DM u_da, DM p_da, Mat *_A11, Mat *_A12, Mat *_A21, Mat *_A22 )
{
	PetscErrorCode ierr;
	//Vec u,p;  need these to determine the parallel layout 
	Mat A11,A12,A21,A22;
	PetscInt rStart, rEnd, cStart, cEnd;
	PetscInt nrows_local, ncols_local;
	PetscMPIInt size;
	PetscInt *dnnz, *onnz;


	MPI_Comm_size( PETSC_COMM_WORLD, &size );

	// A11 (V x V)
	ierr = DMGetMatrix( u_da, MATAIJ, &A11 ); CHKERRQ(ierr);

	// A22 (P x P) 
	ierr = DMGetMatrix( p_da, MATAIJ, &A22 ); CHKERRQ(ierr);

	// A12 (V x P)
	ierr = MatGetOwnershipRange(       A11, &rStart, &rEnd ); CHKERRQ(ierr);
	ierr = MatGetOwnershipRangeColumn( A22, &cStart, &cEnd ); CHKERRQ(ierr);

	nrows_local = rEnd - rStart;
	ncols_local = cEnd - cStart;

	ierr = MatCreate( PETSC_COMM_WORLD, &A12 ); CHKERRQ(ierr);
	ierr = MatSetSizes( A12, nrows_local, ncols_local, PETSC_DECIDE, PETSC_DECIDE); CHKERRQ(ierr);
	ierr = MatSetType( A12, MATAIJ ); CHKERRQ(ierr);

	ierr = PetscMalloc( sizeof(PetscInt)*nrows_local, &dnnz ); CHKERRQ(ierr);
	ierr = PetscMalloc( sizeof(PetscInt)*nrows_local, &onnz ); CHKERRQ(ierr);

	if( size==1 ){
//		MatSeqAIJSetPreallocation( A12, npres*8, PETSC_NULL );
		ierr = MatSeqAIJSetPreallocation( A12, PETSC_NULL, dnnz ); CHKERRQ(ierr);
	}
	else {
//		MatMPIAIJSetPreallocation( A12, npres*8, PETSC_NULL, npres*8, PETSC_NULL );
		ierr = MatMPIAIJSetPreallocation( A12, PETSC_NULL, dnnz, PETSC_NULL, onnz ); CHKERRQ(ierr);
	}

	ierr = PetscFree( dnnz ); CHKERRQ(ierr);
	ierr = PetscFree( onnz ); CHKERRQ(ierr);

	// A21 (P x V)
	ierr = MatGetOwnershipRange(       A22, &rStart, &rEnd ); CHKERRQ(ierr);
	ierr = MatGetOwnershipRangeColumn( A11, &cStart, &cEnd ); CHKERRQ(ierr);

	nrows_local = rEnd - rStart;
	ncols_local = cEnd - cStart;

	ierr = MatCreate( PETSC_COMM_WORLD, &A21 ); CHKERRQ(ierr);
	ierr = MatSetSizes( A21, nrows_local, ncols_local, PETSC_DECIDE, PETSC_DECIDE); CHKERRQ(ierr);
	ierr = MatSetType( A21, MATAIJ ); CHKERRQ(ierr);

	ierr = PetscMalloc( sizeof(PetscInt)*nrows_local, &dnnz ); CHKERRQ(ierr);
	ierr = PetscMalloc( sizeof(PetscInt)*nrows_local, &onnz ); CHKERRQ(ierr);

	if( size==1 ){
//		MatSeqAIJSetPreallocation( A21, edof, PETSC_NULL );
		ierr = MatSeqAIJSetPreallocation( A21, PETSC_NULL, dnnz ); CHKERRQ(ierr);
	}
	else {
//		MatMPIAIJSetPreallocation( A21, edof, PETSC_NULL, edof, PETSC_NULL );
		ierr = MatMPIAIJSetPreallocation( A21, PETSC_NULL, dnnz, PETSC_NULL, onnz ); CHKERRQ(ierr);
	}

	ierr = PetscFree( dnnz ); CHKERRQ(ierr);
	ierr = PetscFree( onnz ); CHKERRQ(ierr);


	PetscFunctionReturn(0);

}
*/


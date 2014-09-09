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

$Id:$
$Date:$

This version is compatible with PETSc version 3.3-p2
The code cannot be distributed, used for commercial purposes or handed on,
without the explicit agreement of Boris Kaus.
*/

#include "LaMEM.h"
#include "Utils.h"

/*

< LOCAL NODE NUMBERING FOR Q1 >

zeta = -1.0;

3-------2
|       |
|       |
|       |
0-------1


zeta =  1.0;

7-------6
|       |
|       |
|       |
4-------5

*/

void ComputeElementIndex_Q1Q1( LaMEMVelPressureDA C, Particles ParticleLocal, PetscInt *_ix, PetscInt *_iy, PetscInt *_iz )
{
	PetscInt ix,iy,iz;

	if(C) C = NULL; // lazy ad-hoc

	ix = ((PetscInt) ParticleLocal.ix);
	iy = ((PetscInt) ParticleLocal.iy);
	iz = ((PetscInt) ParticleLocal.iz);

	*_ix = ix;
	*_iy = iy;
	*_iz = iz;
}

void CreateStencilInGlobalStiffnessTemp_Q1Q1( LaMEMVelPressureDA C, MatStencil *row,MatStencil *col,const PetscInt i,const PetscInt j, const PetscInt k )
{
	PetscInt ii;
	PetscInt edof_temp = C->edof_temp;


	// Q1Q1 element!
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

void GetTemperatureElement_Q1Q1( LaMEMVelPressureDA C, PetscScalar ***temperature, PetscScalar Temp_element[], PetscInt i, PetscInt j, PetscInt k )
{
	if(C) C = NULL; // lazy ad-hoc

	Temp_element[0 ]	=	temperature[k  ][j  ][i  ];
	Temp_element[1 ]	=	temperature[k  ][j  ][i+1];
	Temp_element[2 ]	=	temperature[k  ][j+1][i+1];
	Temp_element[3 ]	=	temperature[k  ][j+1][i  ];
	Temp_element[4 ]	=	temperature[k+1][j  ][i  ];
	Temp_element[5 ]	=	temperature[k+1][j  ][i+1];
	Temp_element[6 ]	=	temperature[k+1][j+1][i+1];
	Temp_element[7 ]	=	temperature[k+1][j+1][i  ];
}

void SetValuesRHS_Temp_Q1Q1( LaMEMVelPressureDA C, PetscScalar ***rhs, PetscScalar T_RHS[], PetscInt i, PetscInt j, PetscInt k )
{
	if(C) C = NULL; // lazy ad-hoc

	rhs[k  ][j  ][i  ]  = rhs[k  ][j  ][i  ] + T_RHS[0 ];
	rhs[k  ][j  ][i+1]  = rhs[k  ][j  ][i+1] + T_RHS[1 ];
	rhs[k  ][j+1][i+1]  = rhs[k  ][j+1][i+1] + T_RHS[2 ];
	rhs[k  ][j+1][i  ]  = rhs[k  ][j+1][i  ] + T_RHS[3 ];
	rhs[k+1][j  ][i  ]  = rhs[k+1][j  ][i  ] + T_RHS[4 ];
	rhs[k+1][j  ][i+1]  = rhs[k+1][j  ][i+1] + T_RHS[5 ];
	rhs[k+1][j+1][i+1]  = rhs[k+1][j+1][i+1] + T_RHS[6 ];
	rhs[k+1][j+1][i  ]  = rhs[k+1][j+1][i  ] + T_RHS[7 ];
}

void FindNearestNode_Q1Q1( LaMEMVelPressureDA C, PetscScalar eta, PetscScalar zetha, PetscScalar phi, PetscInt *ix_add_out, PetscInt *iy_add_out, PetscInt *iz_add_out )
{
	PetscInt	ix_add=0, iy_add=0, iz_add=0;

	if(C) C = NULL; // lazy ad-hoc

	if (eta  <=0){ ix_add = 0; }
	if (eta  > 0){ ix_add = 1; }
	if (zetha<=0){ iy_add = 0; }
	if (zetha> 0){ iy_add = 1; }
	if (phi  <=0){ iz_add = 0; }
	if (phi  > 0){ iz_add = 1; }

	*ix_add_out = ix_add;
	*iy_add_out = iy_add;
	*iz_add_out = iz_add;

}

void ComputeVelocityLocal2Global_Q1Q1(
		LaMEMVelPressureDA C,
		MatStencil *row, MatStencil *col,
		PetscInt local2global[],
		const PetscInt i ,const PetscInt j , const PetscInt k,
		Mat MATRIX )
{
	PetscInt ii;
	PetscInt edof = C->edof;

	// Q2P1 element!
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
	// Compute a local2global vector from the MatStencil values given above
	StencilToLocalNumbering(MATRIX,edof,row, local2global);

}

void CreateStencilInGlobalRhs_Q1Q1( LaMEMVelPressureDA C, PetscInt row_rhs[],const PetscInt i,const PetscInt j, const PetscInt k, const PetscInt mx, const PetscInt my, const PetscInt mz )
{
	if(C)  C = NULL; // lazy ad-hoc
	if(mz) C = NULL; // lazy ad-hoc

	// Q1Q1 element
	GetGlobalIndex(mx, my, i  ,j  ,k  ,0  , &row_rhs[0]);
	GetGlobalIndex(mx, my, i  ,j  ,k  ,1  , &row_rhs[1 ]);
	GetGlobalIndex(mx, my, i  ,j  ,k  ,2  , &row_rhs[2 ]);

	GetGlobalIndex(mx, my, i+1,j  ,k  ,0  , &row_rhs[3 ]);
	GetGlobalIndex(mx, my, i+1,j  ,k  ,1  , &row_rhs[4 ]);
	GetGlobalIndex(mx, my, i+1,j  ,k  ,2  , &row_rhs[5 ]);

	GetGlobalIndex(mx, my, i+1,j+1,k  ,0  , &row_rhs[6 ]);
	GetGlobalIndex(mx, my, i+1,j+1,k  ,1  , &row_rhs[7 ]);
	GetGlobalIndex(mx, my, i+1,j+1,k  ,2  , &row_rhs[8 ]);
	GetGlobalIndex(mx, my, i  ,j+1,k  ,0  , &row_rhs[9 ]);
	GetGlobalIndex(mx, my, i  ,j+1,k  ,1  , &row_rhs[10]);
	GetGlobalIndex(mx, my, i  ,j+1,k  ,2  , &row_rhs[11]);

	GetGlobalIndex(mx, my, i  ,j  ,k+1,0  , &row_rhs[12]);
	GetGlobalIndex(mx, my, i  ,j  ,k+1,1  , &row_rhs[13]);
	GetGlobalIndex(mx, my, i  ,j  ,k+1,2  , &row_rhs[14]);
	GetGlobalIndex(mx, my, i+1,j  ,k+1,0  , &row_rhs[15]);
	GetGlobalIndex(mx, my, i+1,j  ,k+1,1  , &row_rhs[16]);
	GetGlobalIndex(mx, my, i+1,j  ,k+1,2  , &row_rhs[17]);
	GetGlobalIndex(mx, my, i+1,j+1,k+1,0  , &row_rhs[18]);
	GetGlobalIndex(mx, my, i+1,j+1,k+1,1  , &row_rhs[19]);
	GetGlobalIndex(mx, my, i+1,j+1,k+1,2  , &row_rhs[20]);
	GetGlobalIndex(mx, my, i  ,j+1,k+1,0  , &row_rhs[21]);
	GetGlobalIndex(mx, my, i  ,j+1,k+1,1  , &row_rhs[22]);
	GetGlobalIndex(mx, my, i  ,j+1,k+1,2  , &row_rhs[23]);
}

void SetValuesRHS_Q1Q1( LaMEMVelPressureDA C, Field ***rhs, PetscScalar V_RHS[], PetscInt i, PetscInt j, PetscInt k )
{
	if(C) C = NULL; // lazy ad-hoc

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


PetscErrorCode ComputeShapeFunctionVelocity_Q1Q1( LaMEMVelPressureDA C, double ShapeVel[], double **dhdsVel, const double Point[] )
{
	double x,y,z;

	if(C) C = NULL; // lazy ad-hoc

	x = Point[0];
	y = Point[1];
	z = Point[2];

	/* Velocity shape function for Q1 element */
	ShapeVel[0] = 0.125*(1.0-x)*(1.0-y)*(1.0-z);
	ShapeVel[1] = 0.125*(1.0+x)*(1.0-y)*(1.0-z);
	ShapeVel[2] = 0.125*(1.0+x)*(1.0+y)*(1.0-z);
	ShapeVel[3] = 0.125*(1.0-x)*(1.0+y)*(1.0-z);
	ShapeVel[4] = 0.125*(1.0-x)*(1.0-y)*(1.0+z);
	ShapeVel[5] = 0.125*(1.0+x)*(1.0-y)*(1.0+z);
	ShapeVel[6] = 0.125*(1.0+x)*(1.0+y)*(1.0+z);
	ShapeVel[7] = 0.125*(1.0-x)*(1.0+y)*(1.0+z);

	dhdsVel[0][0] = -0.125*(1-y)*(1-z);
	dhdsVel[0][1] =  0.125*(1-y)*(1-z);
	dhdsVel[0][2] =  0.125*(1+y)*(1-z);
	dhdsVel[0][3] = -0.125*(1+y)*(1-z);
	dhdsVel[0][4] = -0.125*(1-y)*(1+z);
	dhdsVel[0][5] =  0.125*(1-y)*(1+z);
	dhdsVel[0][6] =  0.125*(1+y)*(1+z);
	dhdsVel[0][7] = -0.125*(1+y)*(1+z);

	dhdsVel[1][0] = -0.125*(1-x)*(1-z);
	dhdsVel[1][1] = -0.125*(1+x)*(1-z);
	dhdsVel[1][2] =  0.125*(1+x)*(1-z);
	dhdsVel[1][3] =  0.125*(1-x)*(1-z);
	dhdsVel[1][4] = -0.125*(1-x)*(1+z);
	dhdsVel[1][5] = -0.125*(1+x)*(1+z);
	dhdsVel[1][6] =  0.125*(1+x)*(1+z);
	dhdsVel[1][7] =  0.125*(1-x)*(1+z);

	dhdsVel[2][0] = -0.125*(1-x)*(1-y);
	dhdsVel[2][1] = -0.125*(1+x)*(1-y);
	dhdsVel[2][2] = -0.125*(1+x)*(1+y);
	dhdsVel[2][3] = -0.125*(1-x)*(1+y);
	dhdsVel[2][4] =  0.125*(1-x)*(1-y);
	dhdsVel[2][5] =  0.125*(1+x)*(1-y);
	dhdsVel[2][6] =  0.125*(1+x)*(1+y);
	dhdsVel[2][7] =  0.125*(1-x)*(1+y);

	PetscFunctionReturn(0);
}


void ComputeShapeFunctionPressure_Q1Q1( LaMEMVelPressureDA C, double ShapeP[], const DMDACoor3d CoordIntp, double Point[] )
{
	if(C)           C = NULL; // lazy ad-hoc
	if(CoordIntp.x) C = NULL; // lazy ad-hoc
	if(Point[0])    C = NULL; // lazy ad-hoc

	// Global pressure shape function:
	ShapeP[0] = 1.0;
}

PetscErrorCode GetVelocityElement_Q1Q1( LaMEMVelPressureDA C, Field ***velocity, PetscScalar Vel_element[], PetscInt i, PetscInt j, PetscInt k )
{
	if(C) C = NULL; // lazy ad-hoc

	Vel_element[0 ]	=	velocity[k  ][j  ][i  ].Vx;
	Vel_element[1 ]	=	velocity[k  ][j  ][i  ].Vy;
	Vel_element[2 ]	=	velocity[k  ][j  ][i  ].Vz;
	Vel_element[3 ]	=	velocity[k  ][j  ][i+1].Vx;
	Vel_element[4 ]	=	velocity[k  ][j  ][i+1].Vy;
	Vel_element[5 ]	=	velocity[k  ][j  ][i+1].Vz;
	Vel_element[6 ]	=	velocity[k  ][j+1][i+1].Vx;
	Vel_element[7 ]	=	velocity[k  ][j+1][i+1].Vy;
	Vel_element[8 ]	=	velocity[k  ][j+1][i+1].Vz;
	Vel_element[9 ]	=	velocity[k  ][j+1][i  ].Vx;
	Vel_element[10]	=	velocity[k  ][j+1][i  ].Vy;
	Vel_element[11]	=	velocity[k  ][j+1][i  ].Vz;

	Vel_element[12]	=	velocity[k+1][j  ][i  ].Vx;
	Vel_element[13]	=	velocity[k+1][j  ][i  ].Vy;
	Vel_element[14]	=	velocity[k+1][j  ][i  ].Vz;
	Vel_element[15]	=	velocity[k+1][j  ][i+1].Vx;
	Vel_element[16]	=	velocity[k+1][j  ][i+1].Vy;
	Vel_element[17]	=	velocity[k+1][j  ][i+1].Vz;
	Vel_element[18]	=	velocity[k+1][j+1][i+1].Vx;
	Vel_element[19]	=	velocity[k+1][j+1][i+1].Vy;
	Vel_element[20]	=	velocity[k+1][j+1][i+1].Vz;
	Vel_element[21]	=	velocity[k+1][j+1][i  ].Vx;
	Vel_element[22]	=	velocity[k+1][j+1][i  ].Vy;
	Vel_element[23]	=	velocity[k+1][j+1][i  ].Vz;

	PetscFunctionReturn(0);
}


PetscErrorCode GetElementCoords_Q1Q1( LaMEMVelPressureDA C, DMDACoor3d *coord_elem, DMDACoor3d ***coords, PetscInt i,PetscInt j, PetscInt k, PetscInt Flag )
{
	if(C)    C = NULL; // lazy ad-hoc
	if(Flag) C = NULL; // lazy ad-hoc

	coord_elem[0] = coords[k  ][j  ][i  ];
	coord_elem[1] = coords[k  ][j  ][i+1];
	coord_elem[2] = coords[k  ][j+1][i+1];
	coord_elem[3] = coords[k  ][j+1][i  ];
	coord_elem[4] = coords[k+1][j  ][i  ];
	coord_elem[5] = coords[k+1][j  ][i+1];
	coord_elem[6] = coords[k+1][j+1][i+1];
	coord_elem[7] = coords[k+1][j+1][i  ];

	PetscFunctionReturn(0);

}


PetscErrorCode LaMEMVelPressureDACreate_Q1Q1( LaMEMVelPressureDA C )
{

	/* set the element parameters */
	C->nnel        = 8;
	C->ngp_vel     = 8;
	C->nintp_1D    = 2;
	C->ElementType = 3;
	C->npres       = 8;
	C->edof        = 24;
	C->edof_temp   = 8;
	C->nnode_el_1D = 1;
	C->nnel_1D     = 2;
	C->ndim        = 3;

	/* face indices (taken from SetBC_ElementStiffnessMatrix_Symmetric()) */
	/* left (xi = -1.0) */
	C->boundary_nodes_left[0] = 0;
	C->boundary_nodes_left[1] = 3;
	C->boundary_nodes_left[2] = 4;
	C->boundary_nodes_left[3] = 7;
	/* right (xi = 1.0) */
	C->boundary_nodes_right[0] = 1;
	C->boundary_nodes_right[1] = 2;
	C->boundary_nodes_right[2] = 5;
	C->boundary_nodes_right[3] = 6;

	/* front (eta = -1.0) */
	C->boundary_nodes_front[0] = 0;
	C->boundary_nodes_front[1] = 1;
	C->boundary_nodes_front[2] = 4;
	C->boundary_nodes_front[3] = 5;
	/* back (eta = 1.0) */
	C->boundary_nodes_back[0] = 2;
	C->boundary_nodes_back[1] = 3;
	C->boundary_nodes_back[2] = 6;
	C->boundary_nodes_back[3] = 7;

	/* lower (zeta = -1.0) */
	C->boundary_nodes_lower[0] = 0;
	C->boundary_nodes_lower[1] = 1;
	C->boundary_nodes_lower[2] = 2;
	C->boundary_nodes_lower[3] = 3;
	/* upper (zeta = 1.0) */
	C->boundary_nodes_upper[0] = 4;
	C->boundary_nodes_upper[1] = 5;
	C->boundary_nodes_upper[2] = 6;
	C->boundary_nodes_upper[3] = 7;

	/* Set the operations */
	C->fp_ComputeElementIndex = &ComputeElementIndex_Q1Q1;

	C->fp_CreateStencilInGlobalStiffnessTemp = &CreateStencilInGlobalStiffnessTemp_Q1Q1;
	C->fp_GetTemperatureElement              = &GetTemperatureElement_Q1Q1;
	C->fp_SetValuesRHS_Temp                  = &SetValuesRHS_Temp_Q1Q1;
	C->fp_FindNearestNode                    = &FindNearestNode_Q1Q1;

	C->fp_ComputeVelocityLocal2Global = &ComputeVelocityLocal2Global_Q1Q1;
	C->fp_CreateStencilInGlobalRhs    = &CreateStencilInGlobalRhs_Q1Q1;
	C->fp_SetValuesRHS                = &SetValuesRHS_Q1Q1;

	C->fp_ComputeShapeFunctionVelocity = &ComputeShapeFunctionVelocity_Q1Q1;
	C->fp_ComputeShapeFunctionPressure = &ComputeShapeFunctionPressure_Q1Q1;


	C->fp_GetVelocityElement = &GetVelocityElement_Q1Q1;
	C->fp_GetElementCoords   = &GetElementCoords_Q1Q1;


	PetscFunctionReturn(0);
}


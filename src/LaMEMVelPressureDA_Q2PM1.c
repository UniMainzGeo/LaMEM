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

< LOCAL NODE NUMBERING FOR Q2 >

zeta = -1.0;

3---10----2
|         |
|         |
11   20   9
|         |
|         |
0----8----1



zeta = 0.0;

19---25--18
|         |
|         |
22  26   23
|         |
|         |
16---24--17



zeta =  1.0;

7---14----6
|         |
|         |
15  21   13
|         |
|         |
4---12----5

*/
void ComputeElementIndex_Q2PM1( LaMEMVelPressureDA C, Particles ParticleLocal, PetscInt *_ix, PetscInt *_iy, PetscInt *_iz )
{
	PetscInt mod, ix,iy,iz;
	
	if(C) C = NULL; // lazy ad-hoc

	ix = ((PetscInt) ParticleLocal.ix); 
	LaMEMMod(ix, 2, &mod); 
	if (mod>0){  ix = ix-1;  };
	
	iy = ((PetscInt) ParticleLocal.iy); 
	LaMEMMod(iy, 2, &mod); 
	if (mod>0){  iy = iy-1;  };
	
	iz = ((PetscInt) ParticleLocal.iz); 
	LaMEMMod(iz, 2, &mod); 
	if (mod>0){  iz = iz-1;  };
	
	*_ix = ix;
	*_iy = iy;
	*_iz = iz;
}

void CreateStencilInGlobalStiffnessTemp_Q2PM1( LaMEMVelPressureDA C, MatStencil *row,MatStencil *col,const PetscInt i,const PetscInt j, const PetscInt k )
{
	PetscInt ii;
	PetscInt edof_temp = C->edof_temp;

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

void GetTemperatureElement_Q2PM1( LaMEMVelPressureDA C, PetscScalar ***temperature, PetscScalar Temp_element[], PetscInt i, PetscInt j, PetscInt k )
{
	if(C) C = NULL; // lazy ad-hoc

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

void SetValuesRHS_Temp_Q2PM1( LaMEMVelPressureDA C, PetscScalar ***rhs, PetscScalar T_RHS[], PetscInt i, PetscInt j, PetscInt k )
{
	if(C) C = NULL; // lazy ad-hoc

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

void FindNearestNode_Q2PM1( LaMEMVelPressureDA C, PetscScalar eta, PetscScalar zetha, PetscScalar phi, PetscInt *ix_add_out, PetscInt *iy_add_out, PetscInt *iz_add_out )
{
	PetscInt	ix_add=0, iy_add=0, iz_add=0;
	
	if(C) C = NULL; // lazy ad-hoc

	if ( (eta  <=-0.33)       		   	 ){ ix_add = 0; }
	if ( (eta  > -0.33) && (eta  < 0.33) ){ ix_add = 1; }
	if ( (eta  >  0.33) 				 ){ ix_add = 2; }
	if ( (zetha<=-0.33)       		   	 ){ iy_add = 0; }
	if ( (zetha> -0.33) && (zetha< 0.33) ){ iy_add = 1; }
	if ( (zetha>  0.33) 				 ){ iy_add = 2; }
	if ( (phi  <=-0.33)       		   	 ){ iz_add = 0; }
	if ( (phi  > -0.33) && (phi  < 0.33) ){ iz_add = 1; }
	if ( (phi  >  0.33) 				 ){ iz_add = 2; }
	
	*ix_add_out = ix_add;
	*iy_add_out = iy_add;
	*iz_add_out = iz_add;
	
}

void ComputeVelocityLocal2Global_Q2PM1( 
		LaMEMVelPressureDA C,
		MatStencil *row, MatStencil *col, 
		PetscInt local2global[],
		const PetscInt i ,const PetscInt j , const PetscInt k, 
		Mat MATRIX )
{
	PetscInt ii;
	PetscInt edof = C->edof;
	
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
	// Compute a local2global vector from the MatStencil values given above
	StencilToLocalNumbering(MATRIX,edof,row, local2global);
	
}

void CreateStencilInGlobalRhs_Q2PM1( LaMEMVelPressureDA C, PetscInt row_rhs[],const PetscInt i,const PetscInt j, const PetscInt k, const PetscInt mx, const PetscInt my, const PetscInt mz )
{
	if(C)  C = NULL; // lazy ad-hoc
	if(mz) C = NULL; // lazy ad-hoc

	// Q2P1 element	
	GetGlobalIndex(mx, my, i  ,j  ,k  ,0  , &row_rhs[0 ]);
	GetGlobalIndex(mx, my, i  ,j  ,k  ,1  , &row_rhs[1 ]);
	GetGlobalIndex(mx, my, i  ,j  ,k  ,2  , &row_rhs[2 ]);
	GetGlobalIndex(mx, my, i+2,j  ,k  ,0  , &row_rhs[3 ]);
	GetGlobalIndex(mx, my, i+2,j  ,k  ,1  , &row_rhs[4 ]);
	GetGlobalIndex(mx, my, i+2,j  ,k  ,2  , &row_rhs[5 ]);
	GetGlobalIndex(mx, my, i+2,j+2,k  ,0  , &row_rhs[6 ]);
	GetGlobalIndex(mx, my, i+2,j+2,k  ,1  , &row_rhs[7 ]);
	GetGlobalIndex(mx, my, i+2,j+2,k  ,2  , &row_rhs[8 ]);
	GetGlobalIndex(mx, my, i  ,j+2,k  ,0  , &row_rhs[9 ]);
	GetGlobalIndex(mx, my, i  ,j+2,k  ,1  , &row_rhs[10]);
	GetGlobalIndex(mx, my, i  ,j+2,k  ,2  , &row_rhs[11]);
	GetGlobalIndex(mx, my, i  ,j  ,k+2,0  , &row_rhs[12]);
	GetGlobalIndex(mx, my, i  ,j  ,k+2,1  , &row_rhs[13]);
	GetGlobalIndex(mx, my, i  ,j  ,k+2,2  , &row_rhs[14]);
	GetGlobalIndex(mx, my, i+2,j  ,k+2,0  , &row_rhs[15]);
	GetGlobalIndex(mx, my, i+2,j  ,k+2,1  , &row_rhs[16]);
	GetGlobalIndex(mx, my, i+2,j  ,k+2,2  , &row_rhs[17]);
	GetGlobalIndex(mx, my, i+2,j+2,k+2,0  , &row_rhs[18]);
	GetGlobalIndex(mx, my, i+2,j+2,k+2,1  , &row_rhs[19]);
	GetGlobalIndex(mx, my, i+2,j+2,k+2,2  , &row_rhs[20]);
	GetGlobalIndex(mx, my, i  ,j+2,k+2,0  , &row_rhs[21]);
	GetGlobalIndex(mx, my, i  ,j+2,k+2,1  , &row_rhs[22]);
	GetGlobalIndex(mx, my, i  ,j+2,k+2,2  , &row_rhs[23]);
	GetGlobalIndex(mx, my, i+1,j  ,k  ,0  , &row_rhs[24]);
	GetGlobalIndex(mx, my, i+1,j  ,k  ,1  , &row_rhs[25]);
	GetGlobalIndex(mx, my, i+1,j  ,k  ,2  , &row_rhs[26]);
	GetGlobalIndex(mx, my, i+2,j+1,k  ,0  , &row_rhs[27]);
	GetGlobalIndex(mx, my, i+2,j+1,k  ,1  , &row_rhs[28]);
	GetGlobalIndex(mx, my, i+2,j+1,k  ,2  , &row_rhs[29]);
	GetGlobalIndex(mx, my, i+1,j+2,k  ,0  , &row_rhs[30]);
	GetGlobalIndex(mx, my, i+1,j+2,k  ,1  , &row_rhs[31]);
	GetGlobalIndex(mx, my, i+1,j+2,k  ,2  , &row_rhs[32]);
	GetGlobalIndex(mx, my, i  ,j+1,k  ,0  , &row_rhs[33]);
	GetGlobalIndex(mx, my, i  ,j+1,k  ,1  , &row_rhs[34]);
	GetGlobalIndex(mx, my, i  ,j+1,k  ,2  , &row_rhs[35]);
	GetGlobalIndex(mx, my, i+1,j  ,k+2,0  , &row_rhs[36]);
	GetGlobalIndex(mx, my, i+1,j  ,k+2,1  , &row_rhs[37]);
	GetGlobalIndex(mx, my, i+1,j  ,k+2,2  , &row_rhs[38]);
	GetGlobalIndex(mx, my, i+2,j+1,k+2,0  , &row_rhs[39]);
	GetGlobalIndex(mx, my, i+2,j+1,k+2,1  , &row_rhs[40]);
	GetGlobalIndex(mx, my, i+2,j+1,k+2,2  , &row_rhs[41]);
	GetGlobalIndex(mx, my, i+1,j+2,k+2,0  , &row_rhs[42]);
	GetGlobalIndex(mx, my, i+1,j+2,k+2,1  , &row_rhs[43]);
	GetGlobalIndex(mx, my, i+1,j+2,k+2,2  , &row_rhs[44]);
	GetGlobalIndex(mx, my, i  ,j+1,k+2,0  , &row_rhs[45]);
	GetGlobalIndex(mx, my, i  ,j+1,k+2,1  , &row_rhs[46]);
	GetGlobalIndex(mx, my, i  ,j+1,k+2,2  , &row_rhs[47]);
	GetGlobalIndex(mx, my, i  ,j  ,k+1,0  , &row_rhs[48]);
	GetGlobalIndex(mx, my, i  ,j  ,k+1,1  , &row_rhs[49]);
	GetGlobalIndex(mx, my, i  ,j  ,k+1,2  , &row_rhs[50]);
	GetGlobalIndex(mx, my, i+2,j  ,k+1,0  , &row_rhs[51]);
	GetGlobalIndex(mx, my, i+2,j  ,k+1,1  , &row_rhs[52]);
	GetGlobalIndex(mx, my, i+2,j  ,k+1,2  , &row_rhs[53]);
	GetGlobalIndex(mx, my, i+2,j+2,k+1,0  , &row_rhs[54]);
	GetGlobalIndex(mx, my, i+2,j+2,k+1,1  , &row_rhs[55]);
	GetGlobalIndex(mx, my, i+2,j+2,k+1,2  , &row_rhs[56]);
	GetGlobalIndex(mx, my, i  ,j+2,k+1,0  , &row_rhs[57]);
	GetGlobalIndex(mx, my, i  ,j+2,k+1,1  , &row_rhs[58]);
	GetGlobalIndex(mx, my, i  ,j+2,k+1,2  , &row_rhs[59]);
	GetGlobalIndex(mx, my, i+1,j+1,k  ,0  , &row_rhs[60]);
	GetGlobalIndex(mx, my, i+1,j+1,k  ,1  , &row_rhs[61]);
	GetGlobalIndex(mx, my, i+1,j+1,k  ,2  , &row_rhs[62]);
	GetGlobalIndex(mx, my, i+1,j+1,k+2,0  , &row_rhs[63]);
	GetGlobalIndex(mx, my, i+1,j+1,k+2,1  , &row_rhs[64]);
	GetGlobalIndex(mx, my, i+1,j+1,k+2,2  , &row_rhs[65]);
	GetGlobalIndex(mx, my, i  ,j+1,k+1,0  , &row_rhs[66]);
	GetGlobalIndex(mx, my, i  ,j+1,k+1,1  , &row_rhs[67]);
	GetGlobalIndex(mx, my, i  ,j+1,k+1,2  , &row_rhs[68]);
	GetGlobalIndex(mx, my, i+2,j+1,k+1,0  , &row_rhs[69]);
	GetGlobalIndex(mx, my, i+2,j+1,k+1,1  , &row_rhs[70]);
	GetGlobalIndex(mx, my, i+2,j+1,k+1,2  , &row_rhs[71]);
	GetGlobalIndex(mx, my, i+1,j  ,k+1,0  , &row_rhs[72]);
	GetGlobalIndex(mx, my, i+1,j  ,k+1,1  , &row_rhs[73]);
	GetGlobalIndex(mx, my, i+1,j  ,k+1,2  , &row_rhs[74]);
	GetGlobalIndex(mx, my, i+1,j+2,k+1,0  , &row_rhs[75]);
	GetGlobalIndex(mx, my, i+1,j+2,k+1,1  , &row_rhs[76]);
	GetGlobalIndex(mx, my, i+1,j+2,k+1,2  , &row_rhs[77]);
	GetGlobalIndex(mx, my, i+1,j+1,k+1,0  , &row_rhs[78]);
	GetGlobalIndex(mx, my, i+1,j+1,k+1,1  , &row_rhs[79]);
	GetGlobalIndex(mx, my, i+1,j+1,k+1,2  , &row_rhs[80]);
	
}

void SetValuesRHS_Q2PM1( LaMEMVelPressureDA C, Field ***rhs, PetscScalar V_RHS[], PetscInt i, PetscInt j, PetscInt k )
{
	if(C) C = NULL; // lazy ad-hoc

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


PetscErrorCode ComputeShapeFunctionVelocity_Q2PM1( LaMEMVelPressureDA C, double ShapeVel[], double **dhdsVel, const double Point[] )
{
	double x,y,z;
	
	if(C) C = NULL; // lazy ad-hoc

	x = Point[0];
	y = Point[1];
	z = Point[2];
	
	/* Velocity shape function for Q2 element */
	/* shape function */
	
	/* vertices of the cube */
	ShapeVel[0]		=	0.125*x*y*z*(x-1.0)*(y-1.0)*(z-1.0);
	ShapeVel[1]		=	0.125*x*y*z*(x+1.0)*(y-1.0)*(z-1.0);
	ShapeVel[2]		=	0.125*x*y*z*(x+1.0)*(y+1.0)*(z-1.0);
	ShapeVel[3]		=	0.125*x*y*z*(x-1.0)*(y+1.0)*(z-1.0);
	ShapeVel[4]		=	0.125*x*y*z*(x-1.0)*(y-1.0)*(z+1.0);
	ShapeVel[5]		=	0.125*x*y*z*(x+1.0)*(y-1.0)*(z+1.0);
	ShapeVel[6]		=	0.125*x*y*z*(x+1.0)*(y+1.0)*(z+1.0);
	ShapeVel[7]		=	0.125*x*y*z*(x-1.0)*(y+1.0)*(z+1.0);
	
	/* back face (zeta = -1.0) midside nodes (lower, right, upper, left) */
	ShapeVel[8]		=	0.25*y*z*(1.0-x*x)*(y-1.0)*(z-1.0);
	ShapeVel[9]		=	0.25*x*z*(1.0-y*y)*(x+1.0)*(z-1.0);
	ShapeVel[10]	=	0.25*y*z*(1.0-x*x)*(y+1.0)*(z-1.0);
	ShapeVel[11]	=	0.25*x*z*(1.0-y*y)*(x-1.0)*(z-1.0);
	
	/* front face (zeta = -1.0) midside nodes (lower, right, upper, left) */
	ShapeVel[12]	=	0.25*y*z*(1.0-x*x)*(y-1.0)*(z+1.0);
	ShapeVel[13]	=	0.25*x*z*(1.0-y*y)*(x+1.0)*(z+1.0);
	ShapeVel[14]	=	0.25*y*z*(1.0-x*x)*(y+1.0)*(z+1.0);
	ShapeVel[15]	=	0.25*x*z*(1.0-y*y)*(x-1.0)*(z+1.0);
	
	/* vertices of the quad located at zeta = 0.0 */
	ShapeVel[16]	=	0.25*x*y*(x-1.0)*(y-1.0)*(1.0-z*z);
	ShapeVel[17]	=	0.25*x*y*(x+1.0)*(y-1.0)*(1.0-z*z);
	ShapeVel[18]	=	0.25*x*y*(x+1.0)*(y+1.0)*(1.0-z*z);
	ShapeVel[19]	=	0.25*x*y*(x-1.0)*(y+1.0)*(1.0-z*z);
	
	ShapeVel[20]	=	0.5*z*(1.0-x*x)*(1.0-y*y)*(z-1.0); /* center of back face  (zeta = -1.0) */
	ShapeVel[21]	=	0.5*z*(1.0-x*x)*(1.0-y*y)*(z+1.0); /* center of front face (zeta =  1.0) */
	
	/* midside nodes (left, right, lower, upper) for zeta = 0.0 */
	ShapeVel[22]	=	0.5*x*(x-1.0  )*(1.0-y*y)*(1.0-z*z);
	ShapeVel[23]	=	0.5*x*(x+1.0  )*(1.0-y*y)*(1.0-z*z);
	ShapeVel[24]	=	0.5*y*(1.0-x*x)*(y-1.0  )*(1.0-z*z);
	ShapeVel[25]	=	0.5*y*(1.0-x*x)*(y+1.0  )*(1.0-z*z);
	
	ShapeVel[26]	=	(1.0-x*x)*(1.0-y*y)*(1.0-z*z); /* center of cube */
	
	/* Compute derivative of shape function */
	dhdsVel[0][0]		=   0.125*y*z*(x-1.0)*(y-1.0)*(z-1.0)+0.125*x*y*z*(y-1.0)*(z-1.0);
	dhdsVel[0][1]		=   0.125*y*z*(x+1.0)*(y-1.0)*(z-1.0)+0.125*x*y*z*(y-1.0)*(z-1.0);
	dhdsVel[0][2]		=   0.125*y*z*(x+1.0)*(y+1.0)*(z-1.0)+0.125*x*y*z*(y+1.0)*(z-1.0);
	dhdsVel[0][3]		=   0.125*y*z*(x-1.0)*(y+1.0)*(z-1.0)+0.125*x*y*z*(y+1.0)*(z-1.0);
	dhdsVel[0][4]		=   0.125*y*z*(x-1.0)*(y-1.0)*(z+1.0)+0.125*x*y*z*(y-1.0)*(z+1.0);
	dhdsVel[0][5]		=   0.125*y*z*(x+1.0)*(y-1.0)*(z+1.0)+0.125*x*y*z*(y-1.0)*(z+1.0);
	dhdsVel[0][6]		=   0.125*y*z*(x+1.0)*(y+1.0)*(z+1.0)+0.125*x*y*z*(y+1.0)*(z+1.0);
	dhdsVel[0][7]		=   0.125*y*z*(x-1.0)*(y+1.0)*(z+1.0)+0.125*x*y*z*(y+1.0)*(z+1.0);
	
	dhdsVel[0][8]		=   -0.5*x*y*z*(y-1.0)*(z-1.0);
	dhdsVel[0][9]		=   0.25*z*(1.0-y*y)*(x+1.0)*(z-1.0)+0.25*x*z*(1.0-y*y)*(z-1.0);
	dhdsVel[0][10]		=   -0.5*x*y*z*(y+1.0)*(z-1.0);
	dhdsVel[0][11]		=   0.25*z*(1.0-y*y)*(x-1.0)*(z-1.0)+0.25*x*z*(1.0-y*y)*(z-1.0);
	dhdsVel[0][12]		=   -0.5*x*y*z*(y-1.0)*(z+1.0);
	dhdsVel[0][13]		=   0.25*z*(1.0-y*y)*(x+1.0)*(z+1.0)+0.25*x*z*(1.0-y*y)*(z+1.0);
	dhdsVel[0][14]		=   -0.5*x*y*z*(y+1.0)*(z+1.0);
	dhdsVel[0][15]		=   0.25*z*(1.0-y*y)*(x-1.0)*(z+1.0)+0.25*x*z*(1.0-y*y)*(z+1.0);
	dhdsVel[0][16]		=   0.25*y*(x-1.0)*(y-1.0)*(1.0-z*z)+0.25*x*y*(y-1.0)*(1.0-z*z);
	dhdsVel[0][17]		=   0.25*y*(x+1.0)*(y-1.0)*(1.0-z*z)+0.25*x*y*(y-1.0)*(1.0-z*z);
	dhdsVel[0][18]		=   0.25*y*(x+1.0)*(y+1.0)*(1.0-z*z)+0.25*x*y*(y+1.0)*(1.0-z*z);
	dhdsVel[0][19]		=   0.25*y*(x-1.0)*(y+1.0)*(1.0-z*z)+0.25*x*y*(y+1.0)*(1.0-z*z);
	dhdsVel[0][20]		=   -z*x*(1.0-y*y)*(z-1.0);
	dhdsVel[0][21]		=   -z*x*(1.0-y*y)*(z+1.0);
	dhdsVel[0][22]		=   0.5*(x-1.0)*(1.0-y*y)*(1.0-z*z)	+	0.5*x*(1.0-y*y)*(1.0-z*z);
	dhdsVel[0][23]		=   0.5*(x+1.0)*(1.0-y*y)*(1.0-z*z)	+	0.5*x*(1.0-y*y)*(1.0-z*z);
	dhdsVel[0][24]		=   -y*x*(y-1.0)*(1.0-z*z);
	dhdsVel[0][25]		=   -x*y*(y+1.0)*(1.0-z*z);
	dhdsVel[0][26]		=   -2.0*x*(1.0-y*y)*(1.0-z*z);
	
	dhdsVel[1][0]		=	0.125*x*z*(x-1.0)*(y-1.0)*(z-1.0)+0.125*x*y*z*(x-1.0)*(z-1.0);
	dhdsVel[1][1]		=	0.125*x*z*(x+1.0)*(y-1.0)*(z-1.0)+0.125*x*y*z*(x+1.0)*(z-1.0);
	dhdsVel[1][2]		=	0.125*x*z*(x+1.0)*(y+1.0)*(z-1.0)+0.125*x*y*z*(x+1.0)*(z-1.0);
	dhdsVel[1][3]		=	0.125*x*z*(x-1.0)*(y+1.0)*(z-1.0)+0.125*x*y*z*(x-1.0)*(z-1.0);
	dhdsVel[1][4]		=	0.125*x*z*(x-1.0)*(y-1.0)*(z+1.0)+0.125*x*y*z*(x-1.0)*(z+1.0);
	dhdsVel[1][5]		=	0.125*x*z*(x+1.0)*(y-1.0)*(z+1.0)+0.125*x*y*z*(x+1.0)*(z+1.0);
	dhdsVel[1][6]		=	0.125*x*z*(x+1.0)*(y+1.0)*(z+1.0)+0.125*x*y*z*(x+1.0)*(z+1.0);
	dhdsVel[1][7]		=	0.125*x*z*(x-1.0)*(y+1.0)*(z+1.0)+0.125*x*y*z*(x-1.0)*(z+1.0);
	dhdsVel[1][8]    	=	0.25*z*(1.0-x*x)*(y-1.0)*(z-1.0)+0.25*y*z*(1.0-x*x)*(z-1.0);
	dhdsVel[1][9]		=	-0.5*x*y*z*(x+1.0)*(z-1.0);
	dhdsVel[1][10]		=	0.25*z*(1.0-x*x)*(y+1.0)*(z-1.0)+0.25*y*z*(1.0-x*x)*(z-1.0);
	dhdsVel[1][11]		=	-0.5*x*y*z*(x-1.0)*(z-1.0);
	dhdsVel[1][12]		=	0.25*z*(1.0-x*x)*(y-1.0)*(z+1.0)+0.25*y*z*(1.0-x*x)*(z+1.0);
	dhdsVel[1][13]		=	-0.5*x*y*z*(x+1.0)*(z+1.0);
	dhdsVel[1][14]		=	0.25*z*(1.0-x*x)*(y+1.0)*(z+1.0)+0.25*y*z*(1.0-x*x)*(z+1.0);
	dhdsVel[1][15]		=	-0.5*x*y*z*(x-1.0)*(z+1.0);
	dhdsVel[1][16]		=	0.25*x*(x-1.0)*(y-1.0)*(1.0-z*z)+0.25*x*y*(x-1.0)*(1.0-z*z);
	dhdsVel[1][17]		=	0.25*x*(x+1.0)*(y-1.0)*(1.0-z*z)+0.25*x*y*(x+1.0)*(1.0-z*z);
	dhdsVel[1][18]		=	0.25*x*(x+1.0)*(y+1.0)*(1.0-z*z)+0.25*x*y*(x+1.0)*(1.0-z*z);
	dhdsVel[1][19]		=	0.25*x*(x-1.0)*(y+1.0)*(1.0-z*z)+0.25*x*y*(x-1.0)*(1.0-z*z);
	dhdsVel[1][20]		=	-z*(1.0-x*x)*y*(z-1.0);
	dhdsVel[1][21]		=	-z*(1.0-x*x)*y*(z+1.0);
	dhdsVel[1][22]		=	-x*(x-1.0)*y*(1.0-z*z);
	dhdsVel[1][23]		=	-x*(x+1.0)*y*(1.0-z*z);
	dhdsVel[1][24]		=	0.5*(1.0-x*x)*(y-1.0)*(1.0-z*z)+0.5*y*(1.0-x*x)*(1.0-z*z);
	dhdsVel[1][25]		=	0.5*(1.0-x*x)*(y+1.0)*(1.0-z*z)+0.5*y*(1.0-x*x)*(1.0-z*z);
	dhdsVel[1][26]		=	-2.0*(1.0-x*x)*y*(1.0-z*z);
	
	dhdsVel[2][0] 		=    0.125*x*y*(x-1.0)*(y-1.0)*(z-1.0)+0.125*x*y*z*(x-1.0)*(y-1.0);
	dhdsVel[2][1]		=    0.125*x*y*(x+1.0)*(y-1.0)*(z-1.0)+0.125*x*y*z*(x+1.0)*(y-1.0);
	dhdsVel[2][2]		=	 0.125*x*y*(x+1.0)*(y+1.0)*(z-1.0)+0.125*x*y*z*(x+1.0)*(y+1.0);
	dhdsVel[2][3]		=    0.125*x*y*(x-1.0)*(y+1.0)*(z-1.0)+0.125*x*y*z*(x-1.0)*(y+1.0);
	dhdsVel[2][4]		=    0.125*x*y*(x-1.0)*(y-1.0)*(z+1.0)+0.125*x*y*z*(x-1.0)*(y-1.0);
	dhdsVel[2][5]		=    0.125*x*y*(x+1.0)*(y-1.0)*(z+1.0)+0.125*x*y*z*(x+1.0)*(y-1.0);
	dhdsVel[2][6]		=	 0.125*x*y*(x+1.0)*(y+1.0)*(z+1.0)+0.125*x*y*z*(x+1.0)*(y+1.0);
	dhdsVel[2][7]		=    0.125*x*y*(x-1.0)*(y+1.0)*(z+1.0)+0.125*x*y*z*(x-1.0)*(y+1.0);
	dhdsVel[2][8]		=    0.25*y*(1.0-x*x)*(y-1.0)*(z-1.0)+0.25*y*z*(1.0-x*x)*(y-1.0);
	dhdsVel[2][9]		=	 0.25*x*(1.0-y*y)*(x+1.0)*(z-1.0)+0.25*x*z*(1.0-y*y)*(x+1.0);
	dhdsVel[2][10]		=   0.25*y*(1.0-x*x)*(y+1.0)*(z-1.0)+0.25*y*z*(1.0-x*x)*(y+1.0);
	dhdsVel[2][11]		=   0.25*x*(1.0-y*y)*(x-1.0)*(z-1.0)+0.25*x*z*(1.0-y*y)*(x-1.0);
	dhdsVel[2][12]		=   0.25*y*(1.0-x*x)*(y-1.0)*(z+1.0)+0.25*y*z*(1.0-x*x)*(y-1.0);
	dhdsVel[2][13]		=   0.25*x*(1.0-y*y)*(x+1.0)*(z+1.0)+0.25*x*z*(1.0-y*y)*(x+1.0);
	dhdsVel[2][14]		=   0.25*y*(1.0-x*x)*(y+1.0)*(z+1.0)+0.25*y*z*(1.0-x*x)*(y+1.0);
	dhdsVel[2][15]		=   0.25*x*(1.0-y*y)*(x-1.0)*(z+1.0)+0.25*x*z*(1.0-y*y)*(x-1.0);
	dhdsVel[2][16]		=   -0.5*x*y*(x-1.0)*(y-1.0)*z;
	dhdsVel[2][17]		=   -0.5*x*y*z*(x+1.0)*(y-1.0);
	dhdsVel[2][18]		=   -0.5*x*y*(x+1.0)*(y+1.0)*z;
	dhdsVel[2][19]		=   -0.5*x*y*(x-1.0)*(y+1.0)*z;
	dhdsVel[2][20]		=	0.5*(1.0-x*x)*(1.0-y*y)*(z-1.0)+0.5*z*(1.0-x*x)*(1.0-y*y);
	dhdsVel[2][21]		=	0.5*(1.0-x*x)*(1.0-y*y)*(z+1.0)+0.5*z*(1.0-x*x)*(1.0-y*y);
	dhdsVel[2][22]		=   -x*(x-1.0)*(1.0-y*y)*z;
	dhdsVel[2][23]		=   -x*(x+1.0)*(1.0-y*y)*z;
	dhdsVel[2][24]		=	-y*(1.0-x*x)*(y-1.0)*z;
	dhdsVel[2][25]		=   -y*(1.0-x*x)*(y+1.0)*z;
	dhdsVel[2][26]		=	-2*(1.0-x*x)*(1.0-y*y)*z;	
	
	PetscFunctionReturn(0);
} 


void ComputeShapeFunctionPressure_Q2PM1_GLOBAL( LaMEMVelPressureDA C, double ShapeP[], const DMDACoor3d CoordIntp, double Point[] )
{
	if(C)        C = NULL; // lazy ad-hoc
	if(Point[0]) C = NULL; // lazy ad-hoc

	// Global pressure shape function: (WARNING!!! WHERE IS SCLAING???)
	ShapeP[0] = 1.0;
	ShapeP[1] = CoordIntp.x;
	ShapeP[2] = CoordIntp.y;
	ShapeP[3] = CoordIntp.z; 
} 

void ComputeShapeFunctionPressure_Q2PM1_LOCAL( LaMEMVelPressureDA C, double ShapeP[], const DMDACoor3d CoordIntp, double Point[] )
{
	if(C)           C = NULL; // lazy ad-hoc
	if(CoordIntp.x) C = NULL; // lazy ad-hoc

	// Local pressure shape function:
	ShapeP[0] = 1.0; 
	ShapeP[1] = Point[0]; 
	ShapeP[2] = Point[1]; 
	ShapeP[3] = Point[2];
} 

PetscErrorCode GetVelocityElement_Q2PM1( LaMEMVelPressureDA C, Field ***velocity, PetscScalar Vel_element[], PetscInt i, PetscInt j, PetscInt k )
{
	if(C) C = NULL; // lazy ad-hoc

	Vel_element[0 ]	=	velocity[k  ][j  ][i  ].Vx;
	Vel_element[1 ]	=	velocity[k  ][j  ][i  ].Vy;
	Vel_element[2 ]	=	velocity[k  ][j  ][i  ].Vz;
	Vel_element[3 ]	=	velocity[k  ][j  ][i+2].Vx;
	Vel_element[4 ]	=	velocity[k  ][j  ][i+2].Vy;
	Vel_element[5 ]	=	velocity[k  ][j  ][i+2].Vz;
	Vel_element[6 ]	=	velocity[k  ][j+2][i+2].Vx;
	Vel_element[7 ]	=	velocity[k  ][j+2][i+2].Vy;
	Vel_element[8 ]	=	velocity[k  ][j+2][i+2].Vz;
	Vel_element[9 ]	=	velocity[k  ][j+2][i  ].Vx;
	Vel_element[10]	=	velocity[k  ][j+2][i  ].Vy;
	Vel_element[11]	=	velocity[k  ][j+2][i  ].Vz;
	Vel_element[12]	=	velocity[k+2][j  ][i  ].Vx;
	Vel_element[13]	=	velocity[k+2][j  ][i  ].Vy;
	Vel_element[14]	=	velocity[k+2][j  ][i  ].Vz;
	Vel_element[15]	=	velocity[k+2][j  ][i+2].Vx;
	Vel_element[16]	=	velocity[k+2][j  ][i+2].Vy;
	Vel_element[17]	=	velocity[k+2][j  ][i+2].Vz;
	
	Vel_element[18]	=	velocity[k+2][j+2][i+2].Vx;
	Vel_element[19]	=	velocity[k+2][j+2][i+2].Vy;
	Vel_element[20]	=	velocity[k+2][j+2][i+2].Vz;
	
	Vel_element[21]	=	velocity[k+2][j+2][i  ].Vx;
	Vel_element[22]	=	velocity[k+2][j+2][i  ].Vy;
	Vel_element[23]	=	velocity[k+2][j+2][i  ].Vz;
	
	Vel_element[24]	=	velocity[k  ][j  ][i+1].Vx;
	Vel_element[25]	=	velocity[k  ][j  ][i+1].Vy;
	Vel_element[26]	=	velocity[k  ][j  ][i+1].Vz;
	
	Vel_element[27]	=	velocity[k  ][j+1][i+2].Vx;
	Vel_element[28]	=	velocity[k  ][j+1][i+2].Vy;
	Vel_element[29]	=	velocity[k  ][j+1][i+2].Vz;
	
	Vel_element[30]	=	velocity[k  ][j+2][i+1].Vx;
	Vel_element[31]	=	velocity[k  ][j+2][i+1].Vy;
	Vel_element[32]	=	velocity[k  ][j+2][i+1].Vz;
	
	Vel_element[33]	=	velocity[k  ][j+1][i  ].Vx;
	Vel_element[34]	=	velocity[k  ][j+1][i  ].Vy;
	Vel_element[35]	=	velocity[k  ][j+1][i  ].Vz;
	
	Vel_element[36]	=	velocity[k+2][j  ][i+1].Vx;
	Vel_element[37]	=	velocity[k+2][j  ][i+1].Vy;
	Vel_element[38]	=	velocity[k+2][j  ][i+1].Vz;
	
	Vel_element[39]	=	velocity[k+2][j+1][i+2].Vx;
	Vel_element[40]	=	velocity[k+2][j+1][i+2].Vy;
	Vel_element[41]	=	velocity[k+2][j+1][i+2].Vz;
	
	Vel_element[42]	=	velocity[k+2][j+2][i+1].Vx;
	Vel_element[43]	=	velocity[k+2][j+2][i+1].Vy;
	Vel_element[44]	=	velocity[k+2][j+2][i+1].Vz;
	
	Vel_element[45]	=	velocity[k+2][j+1][i  ].Vx;
	Vel_element[46]	=	velocity[k+2][j+1][i  ].Vy;
	Vel_element[47]	=	velocity[k+2][j+1][i  ].Vz;
	
	Vel_element[48]	=	velocity[k+1][j  ][i  ].Vx;
	Vel_element[49]	=	velocity[k+1][j  ][i  ].Vy;
	Vel_element[50]	=	velocity[k+1][j  ][i  ].Vz;
	
	Vel_element[51]	=	velocity[k+1][j  ][i+2].Vx;
	Vel_element[52]	=	velocity[k+1][j  ][i+2].Vy;
	Vel_element[53]	=	velocity[k+1][j  ][i+2].Vz;
	
	Vel_element[54]	=	velocity[k+1][j+2][i+2].Vx;
	Vel_element[55]	=	velocity[k+1][j+2][i+2].Vy;
	Vel_element[56]	=	velocity[k+1][j+2][i+2].Vz;
	
	Vel_element[57]	=	velocity[k+1][j+2][i  ].Vx;
	Vel_element[58]	=	velocity[k+1][j+2][i  ].Vy;
	Vel_element[59]	=	velocity[k+1][j+2][i  ].Vz;
	
	Vel_element[60]	=	velocity[k  ][j+1][i+1].Vx;
	Vel_element[61]	=	velocity[k  ][j+1][i+1].Vy;
	Vel_element[62]	=	velocity[k  ][j+1][i+1].Vz;
	
	Vel_element[63]	=	velocity[k+2][j+1][i+1].Vx;
	Vel_element[64]	=	velocity[k+2][j+1][i+1].Vy;
	Vel_element[65]	=	velocity[k+2][j+1][i+1].Vz;
	
	Vel_element[66]	=	velocity[k+1][j+1][i  ].Vx;
	Vel_element[67]	=	velocity[k+1][j+1][i  ].Vy;
	Vel_element[68]	=	velocity[k+1][j+1][i  ].Vz;
	
	Vel_element[69]	=	velocity[k+1][j+1][i+2].Vx;
	Vel_element[70]	=	velocity[k+1][j+1][i+2].Vy;
	Vel_element[71]	=	velocity[k+1][j+1][i+2].Vz;
	
	Vel_element[72]	=	velocity[k+1][j  ][i+1].Vx;
	Vel_element[73]	=	velocity[k+1][j  ][i+1].Vy;
	Vel_element[74]	=	velocity[k+1][j  ][i+1].Vz;
	
	Vel_element[75]	=	velocity[k+1][j+2][i+1].Vx;
	Vel_element[76]	=	velocity[k+1][j+2][i+1].Vy;
	Vel_element[77]	=	velocity[k+1][j+2][i+1].Vz;
	
	Vel_element[78]	=	velocity[k+1][j+1][i+1].Vx;
	Vel_element[79]	=	velocity[k+1][j+1][i+1].Vy;
	Vel_element[80]	=	velocity[k+1][j+1][i+1].Vz;
	
	PetscFunctionReturn(0);
}


PetscErrorCode GetElementCoords_Q2PM1( LaMEMVelPressureDA C, DMDACoor3d *coord_elem, DMDACoor3d ***coords, PetscInt i,PetscInt j, PetscInt k, PetscInt Flag )
{
	if(C) C = NULL; // lazy ad-hoc
	
	if (Flag==1){
		coord_elem[0 ]=coords[k  ][j  ][i  ];	coord_elem[1 ]=coords[k  ][j  ][i+2];
		coord_elem[2 ]=coords[k  ][j+2][i+2];	coord_elem[3 ]=coords[k  ][j+2][i  ];
		coord_elem[4 ]=coords[k+2][j  ][i  ];	coord_elem[5 ]=coords[k+2][j  ][i+2];
		coord_elem[6 ]=coords[k+2][j+2][i+2];	coord_elem[7 ]=coords[k+2][j+2][i  ];
		coord_elem[8 ]=coords[k  ][j  ][i+1];	coord_elem[9 ]=coords[k  ][j+1][i+2];
		coord_elem[10]=coords[k  ][j+2][i+1];	coord_elem[11]=coords[k  ][j+1][i  ];
		coord_elem[12]=coords[k+2][j  ][i+1];	coord_elem[13]=coords[k+2][j+1][i+2];
		coord_elem[14]=coords[k+2][j+2][i+1];	coord_elem[15]=coords[k+2][j+1][i  ];
		coord_elem[16]=coords[k+1][j  ][i  ];	coord_elem[17]=coords[k+1][j  ][i+2];
		coord_elem[18]=coords[k+1][j+2][i+2];	coord_elem[19]=coords[k+1][j+2][i  ];
		coord_elem[20]=coords[k  ][j+1][i+1];	coord_elem[21]=coords[k+2][j+1][i+1];
		coord_elem[22]=coords[k+1][j+1][i  ];	coord_elem[23]=coords[k+1][j+1][i+2];
		coord_elem[24]=coords[k+1][j  ][i+1];	coord_elem[25]=coords[k+1][j+2][i+1];
		coord_elem[26]=coords[k+1][j+1][i+1];
	}
	else if (Flag==0) { // useful for tracing particles
		coord_elem[0].x=coords[k  ][j  ][i  ].x;	coord_elem[1].x=coords[k  ][j  ][i+1].x; 
		coord_elem[2].x=coords[k  ][j+1][i+1].x;	coord_elem[3].x=coords[k  ][j+1][i  ].x;
		coord_elem[4].x=coords[k+1][j  ][i  ].x;	coord_elem[5].x=coords[k+1][j  ][i+1].x; 
		coord_elem[6].x=coords[k+1][j+1][i+1].x;	coord_elem[7].x=coords[k+1][j+1][i  ].x;
		
		coord_elem[0].y=coords[k  ][j  ][i  ].y;	coord_elem[1].y=coords[k  ][j  ][i+1].y; 
		coord_elem[2].y=coords[k  ][j+1][i+1].y;	coord_elem[3].y=coords[k  ][j+1][i  ].y;
		coord_elem[4].y=coords[k+1][j  ][i  ].y;	coord_elem[5].y=coords[k+1][j  ][i+1].y; 
		coord_elem[6].y=coords[k+1][j+1][i+1].y;	coord_elem[7].y=coords[k+1][j+1][i  ].y;
		
		coord_elem[0].z=coords[k  ][j  ][i  ].z;	coord_elem[1].z=coords[k  ][j  ][i+1].z; 
		coord_elem[2].z=coords[k  ][j+1][i+1].z;	coord_elem[3].z=coords[k  ][j+1][i  ].z;
		coord_elem[4].z=coords[k+1][j  ][i  ].z;	coord_elem[5].z=coords[k+1][j  ][i+1].z; 
		coord_elem[6].z=coords[k+1][j+1][i+1].z;	coord_elem[7].z=coords[k+1][j+1][i  ].z;
	}
	
	PetscFunctionReturn(0);
	
}


PetscErrorCode LaMEMVelPressureDACreate_Q2PM1( LaMEMVelPressureDA C )
{
	
	/* set the element parameters */
	C->nnel        = 27;
	C->ngp_vel     = 27;
	C->nintp_1D    = 3;
	C->ElementType = 1;
	C->npres       = 4;
	C->edof        = 81;
	C->edof_temp   = 27;
	C->nnode_el_1D = 2;
	C->nnel_1D     = 3;
	C->ndim        = 3;
	
	
	
	/* face indices (taken from SetBC_ElementStiffnessMatrix_Symmetric()) */
	/* left (xi = -1.0) */
	C->boundary_nodes_left[0] = 0;  C->boundary_nodes_left[1] = 11; C->boundary_nodes_left[2] = 3;
	C->boundary_nodes_left[3] = 16; C->boundary_nodes_left[4] = 22; C->boundary_nodes_left[5] = 19;
	C->boundary_nodes_left[6] = 4; C->boundary_nodes_left[7] = 15; C->boundary_nodes_left[8] = 7;
	/* right (xi = 1.0) */
	C->boundary_nodes_right[0] = 1;  C->boundary_nodes_right[1] = 9;  C->boundary_nodes_right[2] = 2;
	C->boundary_nodes_right[3] = 17; C->boundary_nodes_right[4] = 23; C->boundary_nodes_right[5] = 18;
	C->boundary_nodes_right[6] = 5; C->boundary_nodes_right[7] = 13; C->boundary_nodes_right[8] = 6;
	
	/* front (eta = -1.0) */
	C->boundary_nodes_front[0] = 0;  C->boundary_nodes_front[1] = 8;  C->boundary_nodes_front[2] = 1; 
	C->boundary_nodes_front[3] = 16; C->boundary_nodes_front[4] = 24; C->boundary_nodes_front[5] = 17; 
	C->boundary_nodes_front[6] = 4; C->boundary_nodes_front[7] = 12; C->boundary_nodes_front[8] = 5;
	/* back (eta = 1.0) */
	C->boundary_nodes_back[0] = 3;  C->boundary_nodes_back[1] = 10; C->boundary_nodes_back[2] = 2; 
	C->boundary_nodes_back[3] = 19; C->boundary_nodes_back[4] = 25; C->boundary_nodes_back[5] = 18; 
	C->boundary_nodes_back[6] = 7; C->boundary_nodes_back[7] = 14; C->boundary_nodes_back[8] = 6;
	
	/* lower (zeta = -1.0) */
	C->boundary_nodes_lower[0] = 0;  C->boundary_nodes_lower[1] = 8;  C->boundary_nodes_lower[2] = 1; 
	C->boundary_nodes_lower[3] = 11; C->boundary_nodes_lower[4] = 20; C->boundary_nodes_lower[5] = 9; 
	C->boundary_nodes_lower[6] = 3; C->boundary_nodes_lower[7] = 10; C->boundary_nodes_lower[8] = 2;
	/* upper (zeta = 1.0) */
	C->boundary_nodes_upper[0] = 4;  C->boundary_nodes_upper[1] = 12; C->boundary_nodes_upper[2] = 5; 
	C->boundary_nodes_upper[3] = 15; C->boundary_nodes_upper[4] = 21; C->boundary_nodes_upper[5] = 13; 
	C->boundary_nodes_upper[6] = 7; C->boundary_nodes_upper[7] = 14; C->boundary_nodes_upper[8] = 6;
	
	/* Set the operations */
	C->fp_ComputeElementIndex = &ComputeElementIndex_Q2PM1;
	
	C->fp_CreateStencilInGlobalStiffnessTemp = &CreateStencilInGlobalStiffnessTemp_Q2PM1;
	C->fp_GetTemperatureElement              = &GetTemperatureElement_Q2PM1;
	C->fp_SetValuesRHS_Temp                  = &SetValuesRHS_Temp_Q2PM1;
	C->fp_FindNearestNode                    = &FindNearestNode_Q2PM1;
	
	C->fp_ComputeVelocityLocal2Global = &ComputeVelocityLocal2Global_Q2PM1;
	C->fp_CreateStencilInGlobalRhs    = &CreateStencilInGlobalRhs_Q2PM1;
	C->fp_SetValuesRHS                = &SetValuesRHS_Q2PM1;
	
	C->fp_ComputeShapeFunctionVelocity = &ComputeShapeFunctionVelocity_Q2PM1;
	/* select local or global basis function for pressure */
	if( C->type == DAVP_Q2PM1G ) { 
		C->fp_ComputeShapeFunctionPressure = &ComputeShapeFunctionPressure_Q2PM1_GLOBAL;
	}
	else if( C->type == DAVP_Q2PM1L ) {
		C->fp_ComputeShapeFunctionPressure = &ComputeShapeFunctionPressure_Q2PM1_LOCAL;
	}
	else {
		PetscPrintf( PETSC_COMM_WORLD, "Unrecognised type for Q2PM1. Choose one of { DAVP_Q2PM1L, DAVP_Q2PM1G }\n" );
		SETERRQ( PETSC_COMM_WORLD, PETSC_ERR_SUP, "Unrecognised type for Q2PM1");
	}
	
	C->fp_GetVelocityElement = &GetVelocityElement_Q2PM1;
	C->fp_GetElementCoords   = &GetElementCoords_Q2PM1;
	
	
	PetscFunctionReturn(0);
}


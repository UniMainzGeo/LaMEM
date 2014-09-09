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

Elements.c, contains the following subroutines:

ComputeShapeFunctionVelocity	-	Compute linear or quadratic velocity shape functions and derivatives
ComputeShapeFunctionPressure	-	Computes pressure shape function
ComputeJacobianElement			-	Computes Jacobian, its inverse and determinant
ComputeCoordIntp				-	Compute real coordinate of a point, given its natural coordinate
GetVelocityElement				-	Extract the velocity from an element into a vector
GetElementCoords				-	Gets coordinates for an element

$Id$
$Date$

This version is compatible with PETSc version 3.3-p2
The code cannot be distributed, used for commercial purposes or handed on,
without the explicit agreement of Boris Kaus.
*/

#include "LaMEM.h"
#include "Elements.h"

/*==========================================================================================================*/
/* Compute shape function for Velocity @ integration point */
#undef __FUNCT__
#define __FUNCT__ "ComputeShapeFunctionVelocity"
PetscErrorCode ComputeShapeFunctionVelocity( double ShapeVel[], double **dhdsVel, const double Point[] )
{
	double x,y,z;

	x		=	Point[0];	y=Point[1];	z=Point[2];

	if( __ELEMENT_TYPE__ == ELEMENT_Q2P1 ) {
		/* Velocity shape function for Q2 element */
		/* shape function */
		ShapeVel[0]		=	0.125*x*y*z*(x-1.0)*(y-1.0)*(z-1.0);	ShapeVel[1]		=	0.125*x*y*z*(x+1.0)*(y-1.0)*(z-1.0);
		ShapeVel[2]		=	0.125*x*y*z*(x+1.0)*(y+1.0)*(z-1.0);	ShapeVel[3]		=	0.125*x*y*z*(x-1.0)*(y+1.0)*(z-1.0);
		ShapeVel[4]		=	0.125*x*y*z*(x-1.0)*(y-1.0)*(z+1.0);    ShapeVel[5]		=	0.125*x*y*z*(x+1.0)*(y-1.0)*(z+1.0);
		ShapeVel[6]		=	0.125*x*y*z*(x+1.0)*(y+1.0)*(z+1.0); 	ShapeVel[7]		=	0.125*x*y*z*(x-1.0)*(y+1.0)*(z+1.0);
		ShapeVel[8]		=	0.25*y*z*(1.0-x*x)*(y-1.0)*(z-1.0);	    ShapeVel[9]		=	0.25*x*z*(1.0-y*y)*(x+1.0)*(z-1.0);
		ShapeVel[10]	=	0.25*y*z*(1.0-x*x)*(y+1.0)*(z-1.0);	    ShapeVel[11]	=	0.25*x*z*(1.0-y*y)*(x-1.0)*(z-1.0);
		ShapeVel[12]	=	0.25*y*z*(1.0-x*x)*(y-1.0)*(z+1.0);		ShapeVel[13]	=	0.25*x*z*(1.0-y*y)*(x+1.0)*(z+1.0);
		ShapeVel[14]	=	0.25*y*z*(1.0-x*x)*(y+1.0)*(z+1.0);		ShapeVel[15]	=	0.25*x*z*(1.0-y*y)*(x-1.0)*(z+1.0);
		ShapeVel[16]	=	0.25*x*y*(x-1.0)*(y-1.0)*(1.0-z*z);		ShapeVel[17]	=	0.25*x*y*(x+1.0)*(y-1.0)*(1.0-z*z);
		ShapeVel[18]	=	0.25*x*y*(x+1.0)*(y+1.0)*(1.0-z*z);		ShapeVel[19]	=	0.25*x*y*(x-1.0)*(y+1.0)*(1.0-z*z);
		ShapeVel[20]	=	0.5*z*(1.0-x*x)*(1.0-y*y)*(z-1.0);		ShapeVel[21]	=	0.5*z*(1.0-x*x)*(1.0-y*y)*(z+1.0);
		ShapeVel[22]	=	0.5*x*(x-1.0  )*(1.0-y*y)*(1.0-z*z);	ShapeVel[23]	=	0.5*x*(x+1.0  )*(1.0-y*y)*(1.0-z*z);
		ShapeVel[24]	=	0.5*y*(1.0-x*x)*(y-1.0  )*(1.0-z*z);	ShapeVel[25]	=	0.5*y*(1.0-x*x)*(y+1.0  )*(1.0-z*z);
		ShapeVel[26]	=	(1.0-x*x)*(1.0-y*y)*(1.0-z*z);

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
	}
	else if( (__ELEMENT_TYPE__ == ELEMENT_Q1P0) ||  (__ELEMENT_TYPE__ == ELEMENT_Q1Q1)  ||  (__ELEMENT_TYPE__ == ELEMENT_FDSTAG)  ){
		/* Velocity shape function for Q1 element
		 *
		 * Note: the FDSTAG grid does not have this shape function for velocity (but a different one instead). Yet, for locating tracers,
		 * we do assume elements (or control volumes) to be similar in the FDSTAG & Q1P0/Q1Q1 cases, which is why we use linear shape functions.
		 * */
		ShapeVel[0]		=	0.125*(1.0-x)*(1.0-y)*(1.0-z);  ShapeVel[1]=0.125*(1.0+x)*(1.0-y)*(1.0-z);
		ShapeVel[2]		=	0.125*(1.0+x)*(1.0+y)*(1.0-z);  ShapeVel[3]=0.125*(1.0-x)*(1.0+y)*(1.0-z);
		ShapeVel[4]		=	0.125*(1.0-x)*(1.0-y)*(1.0+z);  ShapeVel[5]=0.125*(1.0+x)*(1.0-y)*(1.0+z);
		ShapeVel[6]		=	0.125*(1.0+x)*(1.0+y)*(1.0+z);  ShapeVel[7]=0.125*(1.0-x)*(1.0+y)*(1.0+z);

		dhdsVel[0][0]	=  -0.125*(1-y)*(1-z);	dhdsVel[0][1]=0.125*(1-y)*(1-z);dhdsVel[0][2]	=0.125*(1+y)*(1-z);
		dhdsVel[0][3]	=  -0.125*(1+y)*(1-z);	dhdsVel[0][4]=-0.125*(1-y)*(1+z);dhdsVel[0][5]=0.125*(1-y)*(1+z);
		dhdsVel[0][6]	=	0.125*(1+y)*(1+z);	dhdsVel[0][7]=-0.125*(1+y)*(1+z);

		dhdsVel[1][0]	=  -0.125*(1-x)*(1-z);	dhdsVel[1][1]=-0.125*(1+x)*(1-z);dhdsVel[1][2]=0.125*(1+x)*(1-z);
		dhdsVel[1][3]	=   0.125*(1-x)*(1-z);	dhdsVel[1][4]=-0.125*(1-x)*(1+z);dhdsVel[1][5]=-0.125*(1+x)*(1+z);
		dhdsVel[1][6]	=   0.125*(1+x)*(1+z);	dhdsVel[1][7]=0.125*(1-x)*(1+z);

		dhdsVel[2][0]	=  -0.125*(1-x)*(1-y);	dhdsVel[2][1]=  -0.125*(1+x)*(1-y);	dhdsVel[2][2]	=  -0.125*(1+x)*(1+y);
		dhdsVel[2][3]	=  -0.125*(1-x)*(1+y); 	dhdsVel[2][4]=	0.125*(1-x)*(1-y); 	dhdsVel[2][5]	=	0.125*(1+x)*(1-y);
		dhdsVel[2][6]	=	0.125*(1+x)*(1+y);	dhdsVel[2][7]=	0.125*(1-x)*(1+y);
	}
	else {
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Element type not implemented " );
	}

	PetscFunctionReturn(0);
}

/*==========================================================================================================*/
/* Compute the pressure shape function */
#undef __FUNCT__
#define __FUNCT__ "ComputeShapeFunctionPressure"
PetscErrorCode ComputeShapeFunctionPressure( double ShapeP[], const DMDACoor3d CoordIntp, double Point[] )
{
	if( __ELEMENT_TYPE__ == ELEMENT_Q1P0 ) {
		ShapeP[0] = 1.0;
	}
	else if( __ELEMENT_TYPE__ == ELEMENT_Q1Q1 ) {
		{
			double x,y,z;
			x		=	Point[0];	y=Point[1];	z=Point[2];

			/* Linear, continuous P shape function (same as velocity shape fct.) */
			ShapeP[0]		=	0.125*(1.0-x)*(1.0-y)*(1.0-z);  ShapeP[1]=0.125*(1.0+x)*(1.0-y)*(1.0-z);
			ShapeP[2]		=	0.125*(1.0+x)*(1.0+y)*(1.0-z);  ShapeP[3]=0.125*(1.0-x)*(1.0+y)*(1.0-z);
			ShapeP[4]		=	0.125*(1.0-x)*(1.0-y)*(1.0+z);  ShapeP[5]=0.125*(1.0+x)*(1.0-y)*(1.0+z);
			ShapeP[6]		=	0.125*(1.0+x)*(1.0+y)*(1.0+z);  ShapeP[7]=0.125*(1.0-x)*(1.0+y)*(1.0+z);

		}

	}
	else if( __ELEMENT_TYPE__ == ELEMENT_Q2P1 ) {
		if( __Q2_TYPE__ == ELEMENT_Q2P1_GLOBAL ) {
			// Global pressure shape function:
			ShapeP[0]=1.0; ShapeP[1] = CoordIntp.x; ShapeP[2]=CoordIntp.y; ShapeP[3]=CoordIntp.z;
		}
		else if( __Q2_TYPE__ == ELEMENT_Q2P1_LOCAL ) {
			ShapeP[0]=1.0; ShapeP[1] = Point[0]; ShapeP[2]=Point[1]; ShapeP[3]=Point[2];
		}
		//PetscPrintf(PETSC_COMM_WORLD,"ShapeP= %g %g %g %g   \n", ShapeP[0], ShapeP[1], ShapeP[2], ShapeP[3], ShapeP[4]);
	}
	else {
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Element type not implemented for pressure shape function" );
	}

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Compute the Jacobian, its inverse and the determinant */
#undef __FUNCT__
#define __FUNCT__ "ComputeJacobianElement"
PetscErrorCode ComputeJacobianElement( const PetscInt nnel, double **dhdsVel, const DMDACoor3d coord_elem[], double **Jacob,double **InvJacob, double *DetJacob )
{
	PetscInt 	ii, j, k;	//loop control
	double 	a, b, c, d, e,f,g,h,i, coord[MAX_nnel][3];

	for(ii=0;ii<3;ii++){ for(j=0;j<3;j++){
			Jacob[ii][j] = 0;
	}}
	for(ii=0;ii<nnel;ii++){
		coord[ii][0] = coord_elem[ii].x;
		coord[ii][1] = coord_elem[ii].y;
		coord[ii][2] = coord_elem[ii].z;
	}
	for( ii = 0; ii < 3; ii++){
		for( j = 0; j < 3; j++){
			for( k = 0; k < nnel; k++){
				Jacob[ii][j] +=  dhdsVel[ii][k]*coord[k][j];
			}
		}
	}
	/* Compute inverse and determinant (derived with MAPLE) */
	a = Jacob[0][0]; b = Jacob[0][1]; c = Jacob[0][2];
	d = Jacob[1][0]; e = Jacob[1][1]; f = Jacob[1][2];
	g = Jacob[2][0]; h = Jacob[2][1]; i = Jacob[2][2];
	InvJacob[0][0]	=	-(e*i-f*h)/(-a*e*i+a*f*h+d*b*i-d*c*h-g*b*f+g*c*e);
	InvJacob[0][1]	=	 (b*i-c*h)/(-a*e*i+a*f*h+d*b*i-d*c*h-g*b*f+g*c*e);
	InvJacob[0][2]	=	-(b*f-c*e)/(-a*e*i+a*f*h+d*b*i-d*c*h-g*b*f+g*c*e);
	InvJacob[1][0]	=	 (d*i-f*g)/(-a*e*i+a*f*h+d*b*i-d*c*h-g*b*f+g*c*e);
	InvJacob[1][1]	=	-(a*i-c*g)/(-a*e*i+a*f*h+d*b*i-d*c*h-g*b*f+g*c*e);
	InvJacob[1][2]	=	 (a*f-c*d)/(-a*e*i+a*f*h+d*b*i-d*c*h-g*b*f+g*c*e);
	InvJacob[2][0]	=	-(d*h-e*g)/(-a*e*i+a*f*h+d*b*i-d*c*h-g*b*f+g*c*e);
	InvJacob[2][1]	=	-(-a*h+b*g)/(-a*e*i+a*f*h+d*b*i-d*c*h-g*b*f+g*c*e);
	InvJacob[2][2]	=	 (-a*e+b*d)/(-a*e*i+a*f*h+d*b*i-d*c*h-g*b*f+g*c*e);
	*DetJacob    	=	a*e*i-a*f*h-d*b*i+d*c*h+g*b*f-g*c*e;
	if (*DetJacob<0){
		PetscPrintf(PETSC_COMM_WORLD,"Watch out: element too deformed!!! \n");
	//	MPI_Abort(PETSC_COMM_WORLD,1);
	}

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Compute the real coordinates of an integration point */
#undef __FUNCT__
#define __FUNCT__ "ComputeCoordIntp"
PetscErrorCode ComputeCoordIntp( const PetscInt nnel, const double ShapeVel[], const  DMDACoor3d coord_elem[], DMDACoor3d *CoordIntp_out )
{
	PetscInt 	 i;
	DMDACoor3d CoordIntp;

	CoordIntp.x = 0.0; 	CoordIntp.y = 0.0; 	CoordIntp.z = 0.0;
	for (i=0; i<nnel; i++){
		CoordIntp.x = CoordIntp.x + ShapeVel[i]*coord_elem[i].x;
		CoordIntp.y = CoordIntp.y + ShapeVel[i]*coord_elem[i].y;
		CoordIntp.z = CoordIntp.z + ShapeVel[i]*coord_elem[i].z;
	}
	*CoordIntp_out = CoordIntp;

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Get the velocity for an element */
#undef __FUNCT__
#define __FUNCT__ "GetVelocityElement"
PetscErrorCode GetVelocityElement( Field ***velocity, PetscScalar Vel_element[], PetscInt i, PetscInt j, PetscInt k )
{
	if( (__ELEMENT_TYPE__ == ELEMENT_Q1P0) ||   (__ELEMENT_TYPE__ == ELEMENT_Q1Q1) ||   (__ELEMENT_TYPE__ == ELEMENT_FDSTAG) ) {
		/*Q1P0 or Q1Q1 FEM
		 * FDSTAG is 'tricked' here in that we first interpolate the staggered velocities to corner nodes, in order to be able to
		 * use alle existing FE routines. It could be improved
		 * */
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
	}
	else if( __ELEMENT_TYPE__ == ELEMENT_Q2P1 ) {
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
	}
	else {
		SETERRQ( PETSC_COMM_WORLD,PETSC_ERR_SUP, "Element type not implemented" );
	}

	PetscFunctionReturn(0);
}
/*==========================================================================================================*/

/*==========================================================================================================*/
/* Get the velocity for an element */
#undef __FUNCT__
#define __FUNCT__ "GetPressureElement"
PetscErrorCode GetPressureElement( P_array ***pressure, PetscScalar P_element[], PetscInt i, PetscInt j, PetscInt k )
{
	if(  (__ELEMENT_TYPE__ == ELEMENT_Q1Q1) ) {
		P_element[0 ]	=	pressure[k  ][j  ][i  ].p;
		P_element[1 ]	=	pressure[k  ][j  ][i+1].p;
		P_element[2 ]	=	pressure[k  ][j+1][i+1].p;
		P_element[3 ]	=	pressure[k  ][j+1][i  ].p;
		P_element[4 ]	=	pressure[k+1][j  ][i  ].p;
		P_element[5 ]	=	pressure[k+1][j  ][i+1].p;
		P_element[6 ]	=	pressure[k+1][j+1][i  ].p;
		P_element[7]   =    pressure[k+1][j+1][i+1].p;
	}


	PetscFunctionReturn(0);
}
/*==========================================================================================================*/

/*==========================================================================================================*/
/* Get element coordinates */
#undef __FUNCT__
#define __FUNCT__ "GetElementCoords"
PetscErrorCode GetElementCoords( DMDACoor3d *coord_elem, DMDACoor3d ***coords, PetscInt i,PetscInt j, PetscInt k, PetscInt Flag)
{
	if( __ELEMENT_TYPE__ == ELEMENT_Q2P1 ) {	/* Q2P1 element */
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
	}
	else if( (__ELEMENT_TYPE__ == ELEMENT_Q1P0) || (__ELEMENT_TYPE__ == ELEMENT_Q1Q1) || (__ELEMENT_TYPE__ == ELEMENT_FDSTAG)  ) {   /* Q1P0/Q1Q1 element or FDSTAG */
		coord_elem[0]=coords[k  ][j  ][i  ];	coord_elem[1]=coords[k  ][j  ][i+1];
		coord_elem[2]=coords[k  ][j+1][i+1];	coord_elem[3]=coords[k  ][j+1][i  ];
		coord_elem[4]=coords[k+1][j  ][i  ];	coord_elem[5]=coords[k+1][j  ][i+1];
		coord_elem[6]=coords[k+1][j+1][i+1];	coord_elem[7]=coords[k+1][j+1][i  ];
	}
	else {
		SETERRQ( PETSC_COMM_WORLD, PETSC_ERR_SUP, "Element type not implemented" );
	}

	PetscFunctionReturn(0);

}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Get dx, dy & dz for the current element */
#undef __FUNCT__
#define __FUNCT__ "ElementSpacing"
PetscErrorCode ElementSpacing( DMDACoor3d *coord_elem, PetscScalar *dx, PetscScalar *dy, PetscScalar *dz)
{

	*dx = coord_elem[6].x - coord_elem[0].x;
	*dy = coord_elem[6].y - coord_elem[0].y;
	*dz = coord_elem[6].z - coord_elem[0].z;


	PetscFunctionReturn(0);

}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Compute all dx,dy,dz spacing that is required for the FDSTAG discretization (and that includes the neighboring cells). */
#undef __FUNCT__
#define __FUNCT__ "FDSTAG_ComputeSpacing"
PetscErrorCode FDSTAG_ComputeSpacing(DMDACoor3d ***coords, PetscInt ielx,PetscInt iely, PetscInt ielz, UserContext *user,
		PetscScalar dx_NormalVelocity[7], PetscScalar dy_NormalVelocity[7], PetscScalar dz_NormalVelocity[7],
		PetscScalar dx_P[2], PetscScalar dy_P[2], PetscScalar dz_P[2], PetscScalar *z_center)
{
	// NOTE: THIS ROUTINE HAS BUGS AND WILL NOT WORK FOR VARIABLE SPACING IN X,Y,Z DIRECTION !
	// BECAUSE OF TIME ISSUES THIS HAS NOT BEEN FIXED YET, BUT A (PRIMITIVE) ERROR-CHECK HAS BEEN ADDED

	PetscErrorCode 	ierr;
	PetscInt 		num, i, j, k;
	DMDACoor3d		coord_elem[MAX_nnel];
	PetscScalar 	dx,dy,dz, eps=1e-12;
	//PetscScalar 	dx_NormalVelocity[7];
	//PetscScalar 	dy_NormalVelocity[7], dz_NormalVelocity[7];		// spacing required for normal velocities/normal stresses
	//PetscScalar 	dx_P[2], dy_P[2], dz_P[2];

	/* The numbering scheme is as follows:
	 *
	 */


	/* lower-left corner */
	num 					= 	0;
	i 						=	ielx-1; j = iely-1; k = ielz-1;
	if (k< 0){	k = k+1;}							// bottom boundary (ghost node)
	if (i<0){	i = i+1;}	// left boundary (ghost node)
	if (j< 0){	j = j+1;}							// bottom boundary (ghost node)
	ierr 					= 	GetElementCoords(coord_elem, coords, i,j,k, 1);					CHKERRQ(ierr);
	ierr 					= 	CorrectElementCoordsForPeriodicity( coord_elem, i, j, user,1 );	CHKERRQ(ierr);
	ierr 					= 	ElementSpacing(coord_elem, &dx, &dy, &dz);						CHKERRQ(ierr);
	dx_NormalVelocity[num] 	= 	dx;
	dy_NormalVelocity[num] 	= 	dy;
	dz_NormalVelocity[num] 	= 	dz;

	/* Center */
	num 					= 	1;
	i 						=	ielx; j = iely; k = ielz;
	ierr 					= 	GetElementCoords(coord_elem, coords, i,j,k, 1);					CHKERRQ(ierr);
	ierr 					= 	CorrectElementCoordsForPeriodicity( coord_elem, i, j, user,1 );	CHKERRQ(ierr);
	ierr 					= 	ElementSpacing(coord_elem, &dx, &dy, &dz);						CHKERRQ(ierr);
	dx_NormalVelocity[num] 	= 	dx;
	dy_NormalVelocity[num] 	= 	dy;
    dz_NormalVelocity[num] 	= 	dz;

	/* z-coordinate of center point */
	*z_center				=	(coords[k  ][j  ][i  ].z + coords[k  ][j+1][i  ].z + coords[k  ][j  ][i+1].z + coords[k  ][j+1][i+1].z +
								 coords[k+1][j  ][i  ].z + coords[k+1][j+1][i  ].z + coords[k+1][j  ][i+1].z + coords[k+1][j+1][i+1].z)/8.0;


	/* error detection */
	if ( (PetscAbsScalar(dx_NormalVelocity[1] - dx_NormalVelocity[0])>eps) || (PetscAbsScalar(dy_NormalVelocity[1] - dy_NormalVelocity[0])>eps) || (PetscAbsScalar(dz_NormalVelocity[1] - dz_NormalVelocity[0])>eps) ){

		PetscPrintf(PETSC_COMM_WORLD," dx=[%g,%g], dy=[%g,%g], dz=[%g,%g]  \n",dx_NormalVelocity[0],dx_NormalVelocity[1],dy_NormalVelocity[0],dy_NormalVelocity[1],dz_NormalVelocity[0],dz_NormalVelocity[1] );
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Code with FDSTAG currently only works with constant dx,dy,dz! [%g, %g] ");

	}

	/* Set constant dx/dy/dz everywhere */
	for (num=0; num<7; num++){
		dx_NormalVelocity[num] 	= 	dx;
		dy_NormalVelocity[num] 	= 	dy;
		dz_NormalVelocity[num] 	= 	dz;
	}

	
#if 0
	/* West */
	num 					= 	0;
	i 						=	ielx-1; j = iely; k = ielz;
	if (i<0){	i = i+1;}	// left boundary (ghost node)
	ierr 					= 	GetElementCoords(coord_elem, coords, i,j,k, 1);					CHKERRQ(ierr);
	ierr 					= 	CorrectElementCoordsForPeriodicity( coord_elem, i, j, user,1 );	CHKERRQ(ierr);
	ierr 					= 	ElementSpacing(coord_elem, &dx, &dy, &dz);						CHKERRQ(ierr);
	dx_NormalVelocity[num] 	= 	dx;
	dy_NormalVelocity[num] 	= 	dy;
	dz_NormalVelocity[num] 	= 	dz;

	/* Center */
	num 					= 	1;
	i 						=	ielx; j = iely; k = ielz;
	ierr 					= 	GetElementCoords(coord_elem, coords, i,j,k, 1);					CHKERRQ(ierr);
	ierr 					= 	CorrectElementCoordsForPeriodicity( coord_elem, i, j, user,1 );	CHKERRQ(ierr);
	ierr 					= 	ElementSpacing(coord_elem, &dx, &dy, &dz);						CHKERRQ(ierr);
	dx_NormalVelocity[num] 	= 	dx;
	dy_NormalVelocity[num] 	= 	dy;
	dz_NormalVelocity[num] 	= 	dz;


	/* East */
	num 					= 	2;
	i 						=	ielx+1; j = iely; k = ielz;
	if (i== user->finest_nnode_x-1){	i = i-1;}	// right boundary (ghost node)
	ierr 					= 	GetElementCoords(coord_elem, coords, i,j,k, 1);					CHKERRQ(ierr);
	ierr 					= 	CorrectElementCoordsForPeriodicity( coord_elem, i, j, user,1 );	CHKERRQ(ierr);
	ierr 					= 	ElementSpacing(coord_elem, &dx, &dy, &dz);						CHKERRQ(ierr);
	dx_NormalVelocity[num] 	= 	dx;
	dy_NormalVelocity[num] 	= 	dy;
	dz_NormalVelocity[num] 	= 	dz;

	/* South */
	num 					= 	3;
	i 						=	ielx; j = iely-1; k = ielz;
	if (j< 0){	j = j+1;}							// south boundary (ghost node)
	ierr 					= 	GetElementCoords(coord_elem, coords, i,j,k, 1);					CHKERRQ(ierr);
	ierr 					= 	CorrectElementCoordsForPeriodicity( coord_elem, i, j, user,1 );	CHKERRQ(ierr);
	ierr 					= 	ElementSpacing(coord_elem, &dx, &dy, &dz);						CHKERRQ(ierr);
	dx_NormalVelocity[num] 	= 	dx;
	dy_NormalVelocity[num] 	= 	dy;
	dz_NormalVelocity[num] 	= 	dz;

	/* North */
	num 					= 	4;
	i 						=	ielx; j = iely+1; k = ielz;
	if (j== user->finest_nnode_y-1){	j = j-1;}	// north boundary (ghost node)
	ierr 					= 	GetElementCoords(coord_elem, coords, i,j,k, 1);					CHKERRQ(ierr);
	ierr 					= 	CorrectElementCoordsForPeriodicity( coord_elem, i, j, user,1 );	CHKERRQ(ierr);
	ierr 					= 	ElementSpacing(coord_elem, &dx, &dy, &dz);						CHKERRQ(ierr);
	dx_NormalVelocity[num] 	= 	dx;
	dy_NormalVelocity[num] 	= 	dy;
	dz_NormalVelocity[num] 	= 	dz;


	/* Bottom */
	num 					= 	5;
	i 						=	ielx; j = iely; k = ielz-1;
	if (k< 0){	k = k+1;}							// bottom boundary (ghost node)
	ierr 					= 	GetElementCoords(coord_elem, coords, i,j,k, 1);					CHKERRQ(ierr);
	ierr 					= 	CorrectElementCoordsForPeriodicity( coord_elem, i, j, user,1 );	CHKERRQ(ierr);
	ierr 					= 	ElementSpacing(coord_elem, &dx, &dy, &dz);						CHKERRQ(ierr);
	dx_NormalVelocity[num] 	= 	dx;
	dy_NormalVelocity[num] 	= 	dy;
	dz_NormalVelocity[num] 	= 	dz;

	/* Top */
	num 					= 	6;
	i 						=	ielx; j = iely; k = ielz+1;
	if (k== user->finest_nnode_z-1){	k = k-1;}	// top boundary (ghost node)
	ierr 					= 	GetElementCoords(coord_elem, coords, i,j,k, 1);					CHKERRQ(ierr);
	ierr 					= 	CorrectElementCoordsForPeriodicity( coord_elem, i, j, user,1 );	CHKERRQ(ierr);
	ierr 					= 	ElementSpacing(coord_elem, &dx, &dy, &dz);						CHKERRQ(ierr);
	dx_NormalVelocity[num] 	= 	dx;
	dy_NormalVelocity[num] 	= 	dy;
	dz_NormalVelocity[num] 	= 	dz;

#endif


	/* Compute spacing for pressure points */
	dx_P[0] 				= 	(dx_NormalVelocity[1]+dx_NormalVelocity[0])/2.0; 	// west
	dx_P[1] 				= 	(dx_NormalVelocity[1]+dx_NormalVelocity[2])/2.0; 	// east
	dy_P[0] 				= 	(dy_NormalVelocity[1]+dy_NormalVelocity[3])/2.0; 	// south
	dy_P[1] 				= 	(dy_NormalVelocity[1]+dy_NormalVelocity[4])/2.0; 	// north
	dz_P[0] 				= 	(dz_NormalVelocity[1]+dz_NormalVelocity[5])/2.0; 	// bottom
	dz_P[1] 				= 	(dz_NormalVelocity[1]+dz_NormalVelocity[6])/2.0; 	// top



	PetscFunctionReturn(0);

}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Extracts viscosity at center and corner points (which we need to form the FD stencil) */
#undef __FUNCT__
#define __FUNCT__ "FDSTAG_ExtractMaterialParameters"
PetscErrorCode FDSTAG_ExtractMaterialParameters(PetscScalar ***viscosity_center, PetscScalar ***viscosity_XY,
		PetscScalar ***viscosity_YZ, PetscScalar ***viscosity_XZ, PetscScalar ***density_center, PetscScalar ***PhaseProportionsAir_Center, PetscScalar ***LocalSurfaceTopography,
		PetscInt ielx,PetscInt iely, PetscInt ielz, PetscInt zs_FreeSurface, PetscScalar dx_P[2], PetscScalar dy_P[2], PetscScalar dz_P[2], UserContext *user,
		PetscScalar ViscosityCenter[7], PetscScalar Viscosity_XY[4],PetscScalar Viscosity_YZ[4],PetscScalar Viscosity_XZ[4], PetscScalar dRho_dxdydz[3],
		PetscInt FreeSurfaceCells[7], PetscScalar *z_FreeSurface)
{
	PetscInt 	i,j,k, num;

	/* Extract viscosities and label cells as free surface @ centers ----------- */
	for (num=0; num<7; num++){FreeSurfaceCells[num] 	=	0; ViscosityCenter[num] = 0.0;}		// initialize

	// West
	num 					= 	0;
	i 						=	ielx-1; j = iely; k = ielz;
	if (i<0){	i = i+1;}							// left boundary (ghost node)
	ViscosityCenter[num] 	=	viscosity_center[k][j][i];
	if ((PhaseProportionsAir_Center[k][j][i]>0.99) && (user->ErosionParameters.UseInternalFreeSurface==1)){
		FreeSurfaceCells[num] 	=	1;
	}

	// Center
	num 					= 	1;
	i 						=	ielx; j = iely; k = ielz;
	ViscosityCenter[num] 	=	viscosity_center[k][j][i];
	if ((PhaseProportionsAir_Center[k][j][i]>0.99) && (user->ErosionParameters.UseInternalFreeSurface==1)){
		FreeSurfaceCells[num] 	=	1;
	}

	// East
	num 					= 	2;
	i 						=	ielx; j = iely; k = ielz;
	if (i== user->finest_nnode_x-1){	i = i-1;}	// right boundary (ghost node)
	ViscosityCenter[num] 	=	viscosity_center[k][j][i];
	if ((PhaseProportionsAir_Center[k][j][i]>0.99) && (user->ErosionParameters.UseInternalFreeSurface==1)){
		FreeSurfaceCells[num] 	=	1;
	}

	/* South */
	num 					= 	3;
	i 						=	ielx; j = iely-1; k = ielz;
	if (j< 0){	j = j+1;}							// south boundary (ghost node)
	ViscosityCenter[num] 	=	viscosity_center[k][j][i];
	if ((PhaseProportionsAir_Center[k][j][i]>0.99) && (user->ErosionParameters.UseInternalFreeSurface==1)){
		FreeSurfaceCells[num] 	=	1;
	}

	/* North */
	num 					= 	4;
	i 						=	ielx; j = iely+1; k = ielz;
	if (j== user->finest_nnode_y-1){	j = j-1;}	// north boundary (ghost node)
	ViscosityCenter[num] 	=	viscosity_center[k][j][i];
	if ((PhaseProportionsAir_Center[k][j][i]>0.99) && (user->ErosionParameters.UseInternalFreeSurface==1)){
		FreeSurfaceCells[num] 	=	1;
	}

	/* Bottom */
	num 					= 	5;
	i 						=	ielx; j = iely; k = ielz-1;
	if (k< 0){	k = k+1;}							// bottom boundary (ghost node)
	ViscosityCenter[num] 	=	viscosity_center[k][j][i];
	if ((PhaseProportionsAir_Center[k][j][i]>0.99) && (user->ErosionParameters.UseInternalFreeSurface==1)){
		FreeSurfaceCells[num] 	=	1;
	}
	/* Top */
	num 					= 	6;
	i 						=	ielx; j = iely; k = ielz+1;
	if (k== user->finest_nnode_z-1){	k = k-1;}	// top boundary (ghost node)
	ViscosityCenter[num] 	=	viscosity_center[k][j][i];
	if ((PhaseProportionsAir_Center[k][j][i]>0.99) && (user->ErosionParameters.UseInternalFreeSurface==1)){
		FreeSurfaceCells[num] 	=	1;
	}
	/* ------------------------------------------------------------------------------- */

	/* Extract viscosities @ Sxy points ---------------------------------------------- */
	i 					=	ielx; j = iely; k = ielz;
	Viscosity_XY[0]		=	viscosity_XY[k  ][j  ][i  ];	// SouthWest
	Viscosity_XY[1]		=	viscosity_XY[k  ][j  ][i+1];	// SouthEast
	Viscosity_XY[2]		=	viscosity_XY[k  ][j+1][i  ];	// NorthWest
	Viscosity_XY[3]		=	viscosity_XY[k  ][j+1][i+1];	// NorthEast
	/* ---------------------------------------------------------------------------- 	*/

	/* Extract viscosities @ Syz points ---------------------------------------------- */
	i 					=	ielx; j = iely; k = ielz;
	Viscosity_YZ[0]		=	viscosity_YZ[k  ][j  ][i  ];	// BottomSouth
	Viscosity_YZ[1]		=	viscosity_YZ[k+1][j  ][i  ];	// TopSouth
	Viscosity_YZ[2]		=	viscosity_YZ[k  ][j+1][i  ];	// BottomNorth
	Viscosity_YZ[3]		=	viscosity_YZ[k+1][j+1][i  ];	// TopNorth
	/* ---------------------------------------------------------------------------- 	*/

	/* Extract viscosities @ Sxz points ---------------------------------------------- */
	i 					=	ielx; j = iely; k = ielz;
	Viscosity_XZ[0]		=	viscosity_XZ[k  ][j  ][i  ];	// BottomWest
	Viscosity_XZ[1]		=	viscosity_XZ[k+1][j  ][i  ];	// TopWest
	Viscosity_XZ[2]		=	viscosity_XZ[k  ][j  ][i+1];	// BottomEast
	Viscosity_XZ[3]		=	viscosity_XZ[k+1][j  ][i+1];	// TopEast
	/* ---------------------------------------------------------------------------- 	*/

	/* Compute derivatives of density vs x,y & z direction -------------------------	*/
	/* This is required for the FSSA algorithm */
	i 						=	ielx; j = iely; k = ielz;
	if ((i>0) && (i<user->nel_x-1)){
		dRho_dxdydz[0]		=	(density_center[k  ][j  ][i+1]-density_center[k  ][j  ][i-1])/(dx_P[1]+dx_P[0]);	// drho/dx, computed using central derivatives
	}
	else if (i==0){
		dRho_dxdydz[0]		=	(density_center[k  ][j  ][i+1]-density_center[k  ][j  ][i  ])/(dx_P[1]);			// drho/dx, computed using forward derivatives
	}
	else if (i==(user->nel_x-1)){
		dRho_dxdydz[0]		=	(density_center[k  ][j  ][i  ]-density_center[k  ][j  ][i-1])/(dx_P[0]);			// drho/dx, computed using backwards derivatives
	}

	if ((j>0) && (j<user->nel_y-1)){
		dRho_dxdydz[1]		=	(density_center[k  ][j+1][i  ]-density_center[k  ][j-1][i  ])/(dy_P[1]+dy_P[0]);	// drho/dy, computed using central derivatives
	}
	else if (j==0){
		dRho_dxdydz[1]		=	(density_center[k  ][j+1][i  ]-density_center[k  ][j  ][i  ])/(dy_P[1]);			// drho/dy, computed using forward derivatives
	}
	else if (j==(user->nel_y-1)){
		dRho_dxdydz[1]		=	(density_center[k  ][j-1][i  ]-density_center[k  ][j  ][i  ])/(dy_P[0]);			// drho/dy, computed using backwards derivatives
	}

	if ((k>0) && (k<user->nel_z-1)){
		dRho_dxdydz[2]		=	(density_center[k+1][j  ][i  ]-density_center[k-1][j  ][i  ])/(dz_P[1]+dz_P[0]);	// drho/dz, computed using central derivatives
	}
	else if (k==0){
		dRho_dxdydz[2]		=	(density_center[k+1][j  ][i  ]-density_center[k  ][j  ][i  ])/(dz_P[1]);			// drho/dz, computed using forward derivatives
	}
	else if (k==(user->nel_z-1)){
		dRho_dxdydz[2]		=	(density_center[k-1][j  ][i  ]-density_center[k  ][j  ][i  ])/(dz_P[0]);			// drho/dz, computed using backwards derivatives
	}
	/* ---------------------------------------------------------------------------- 	*/


	/* Height of the free surface above the center of the current cell*/
	*z_FreeSurface 			=	(	LocalSurfaceTopography[zs_FreeSurface][iely  ][ielx  ] + LocalSurfaceTopography[zs_FreeSurface][iely+1][ielx  ]	+
									LocalSurfaceTopography[zs_FreeSurface][iely+1][ielx+1] + LocalSurfaceTopography[zs_FreeSurface][iely+1][ielx  ]	)/4.0;





	PetscFunctionReturn(0);

}
/*==========================================================================================================*/


/*==========================================================================================================*/
/* Get element coordinates */
#undef __FUNCT__
#define __FUNCT__ "CorrectElementCoordsForPeriodicity"
PetscErrorCode CorrectElementCoordsForPeriodicity( DMDACoor3d *coord_elem, PetscInt ix, PetscInt iy, UserContext *user, PetscInt Flag)
{
	if( __ELEMENT_TYPE__ == ELEMENT_Q2P1 ) {	/* Q2P1 element */
		if (Flag==1){
			if ((ix==user->finest_nnode_x-3) && (user->BC.LeftBound==3)){
				coord_elem[1].x  = coord_elem[1].x	+	user->W ;
				coord_elem[2].x  = coord_elem[2].x	+	user->W ;
				coord_elem[5].x  = coord_elem[5].x	+	user->W ;
				coord_elem[6].x  = coord_elem[6].x	+	user->W ;
				coord_elem[9].x  = coord_elem[9].x	+	user->W ;
				coord_elem[13].x = coord_elem[13].x	+	user->W ;
				coord_elem[17].x = coord_elem[17].x	+	user->W ;
				coord_elem[18].x = coord_elem[18].x	+	user->W ;
				coord_elem[23].x = coord_elem[23].x	+	user->W ;
			}
			if ((iy==user->finest_nnode_y-3) && (user->BC.FrontBound==3)){

				coord_elem[2].y  = coord_elem[2].y	+	user->L ;
				coord_elem[3].y  = coord_elem[3].y	+	user->L ;
				coord_elem[6].y  = coord_elem[6].y	+	user->L ;
				coord_elem[7].y  = coord_elem[7].y	+	user->L ;
				coord_elem[10].y = coord_elem[10].y	+	user->L ;
				coord_elem[14].y = coord_elem[14].y	+	user->L ;
				coord_elem[18].y = coord_elem[18].y	+	user->L ;
				coord_elem[19].y = coord_elem[19].y +	user->L ;
				coord_elem[25].y = coord_elem[25].y	+	user->L ;

			}

		}
		else if (Flag==0){
			if ((ix==user->finest_nnode_x-2) && (user->BC.LeftBound==3)){
				coord_elem[1].x = coord_elem[1].x	+	user->W ;
				coord_elem[2].x = coord_elem[2].x	+	user->W ;
				coord_elem[5].x = coord_elem[5].x	+	user->W ;
				coord_elem[6].x = coord_elem[6].x	+	user->W ;
			}
			if ((iy==user->finest_nnode_y-2) && (user->BC.FrontBound==3)){
				coord_elem[2].y = coord_elem[2].y	+	user->L ;
				coord_elem[3].y = coord_elem[3].y	+	user->L ;
				coord_elem[6].y = coord_elem[6].y	+	user->L ;
				coord_elem[7].y = coord_elem[7].y	+	user->L ;
			}
		}
	}
	else if( (__ELEMENT_TYPE__ == ELEMENT_Q1P0) ||  (__ELEMENT_TYPE__ == ELEMENT_Q1Q1)  ||  (__ELEMENT_TYPE__ == ELEMENT_FDSTAG)) {   /* Q1 element or FDSTAG */
		if ((ix==user->finest_nnode_x-2) && (user->BC.LeftBound==3)){
			coord_elem[1].x = coord_elem[1].x	+	user->W ;
			coord_elem[2].x = coord_elem[2].x	+	user->W ;
			coord_elem[5].x = coord_elem[5].x	+	user->W ;
			coord_elem[6].x = coord_elem[6].x	+	user->W ;
		}
		if ((iy==user->finest_nnode_y-2) && (user->BC.FrontBound==3)){
			coord_elem[2].y = coord_elem[2].y	+	user->L ;
			coord_elem[3].y = coord_elem[3].y	+	user->L ;
			coord_elem[6].y = coord_elem[6].y	+	user->L ;
			coord_elem[7].y = coord_elem[7].y	+	user->L ;
		}


	}
	else {

	}

	PetscFunctionReturn(0);

}
/*==========================================================================================================*/

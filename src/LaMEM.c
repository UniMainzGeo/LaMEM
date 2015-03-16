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

LaMEM.c, contains the following subroutine:

main -	Main routine

$Id$
$Date$

This version is compatible with PETSc version 3.3-p2
The code cannot be distributed, used for commercial purposes or handed on,
without the explicit agreement of Boris Kaus.
*/
static char help[] = "Solves 3D Stokes equations using multigrid .\n\n";

#include "LaMEM.h"

extern PetscErrorCode PCCreate_SemiRedundant(PC);

/*==========================================================================================================*/
/* Main routine */
#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv )
{
	PetscErrorCode 	ierr;
	PetscScalar		LaMEM_OutputParameters=0.0;
	PetscInt		mpi_group_id=0;

	// Initialize PETSC
	ierr = PetscInitialize(&argc,&argv,(char *)0,help); CHKERRQ(ierr);
	
    ierr = PCRegister("pc_semiredundant",PCCreate_SemiRedundant);
    
    
	// call LaMEM main library function
	ierr = LaMEMLib(&LaMEM_OutputParameters,&mpi_group_id); CHKERRQ(ierr);

	// cleanup PETSC
	ierr = PetscFinalize(); CHKERRQ(ierr);
	return 0;
}
/* End of code */
/*==========================================================================================================*/


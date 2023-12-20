/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : LaMEM
 **   License      : MIT, see LICENSE file for details
 **   Contributors : Anton Popov, Boris Kaus, see AUTHORS file for complete list
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : kaus@uni-mainz.de, popov@uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
#include "LaMEM.h"
#include "scaling.h"
#include "objFunct.h"
#include "parsing.h"
#include "adjoint.h"
#include "phase.h"
//---------------------------------------------------------------------------
static char help[] = "Solves 3D Stokes equations using multigrid .\n\n";
//---------------------------------------------------------------------------
int main(int argc, char **argv)
{
	PetscErrorCode 	ierr;

	// Initialize PETSC
	ierr = PetscInitialize(&argc,&argv,(char *)0, help); CHKERRQ(ierr);
	ModParam IOparam;
	char      str[_str_len_];

	// set default to be a forward run and overwrite it with input file options
	ierr = PetscMalloc(sizeof(ModParam), &IOparam);  CHKERRQ(ierr);
	ierr = PetscMemzero(&IOparam, sizeof(ModParam)); CHKERRQ(ierr);

	IOparam.use = _none_;
	ierr = FBLoad(&IOparam.fb, PETSC_FALSE); CHKERRQ(ierr);
	ierr = getStringParam(IOparam.fb, _OPTIONAL_, "Adjoint_mode", str, "None"); CHKERRQ(ierr);
	if     (!strcmp(str, "None"))                   IOparam.use = _none_;
	else if(!strcmp(str, "GenericInversion"))       IOparam.use = _inversion_;
	else if(!strcmp(str, "AdjointGradients"))       IOparam.use = _adjointgradients_;
	else if(!strcmp(str, "GradientDescent"))        IOparam.use = _gradientdescent_;
	else if(!strcmp(str, "SyntheticForwardRun"))    IOparam.use = _syntheticforwardrun_;
	else{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Unknown parameter for 'Adjoint_mode'. Possibilities are [None; GenericInversion; AdjointGradients; GradientDescent or SyntheticForwardRun]");
	} 
	
	/* Name stages */
	PetscCall(PetscLogStageRegister("Initial guess",  &IOparam.stages[0])); 
	PetscCall(PetscLogStageRegister("SNES solve",     &IOparam.stages[1]));
	PetscCall(PetscLogStageRegister("Advect markers", &IOparam.stages[2])); 
	PetscCall(PetscLogStageRegister("I/O",            &IOparam.stages[3]));

	if(IOparam.use == 0)
	{
		// Forward simulation	
		ierr = LaMEMLibMain(NULL,IOparam.stages); CHKERRQ(ierr);
	}
	else
	{
		// Inversion or adjoint gradient computation
		ierr = LaMEMAdjointMain(&IOparam); CHKERRQ(ierr);
	}

	// destroy file buffer
	ierr = FBDestroy(&IOparam.fb); CHKERRQ(ierr);

	// cleanup PETSC
	ierr = PetscFinalize(); CHKERRQ(ierr);

	return 0;
}
//--------------------------------------------------------------------------

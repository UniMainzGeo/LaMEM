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
#include "phase.h"
#include "objFunct.h"
#include "parsing.h"
#include "options.h"
#include "adjoint.h"

//---------------------------------------------------------------------------
static char help[] = "Solves 3D Stokes equations using multigrid .\n\n";
//---------------------------------------------------------------------------
int main(int argc, char **argv)
{
	FB       *fb;
	ModParam IOparam;
	char     str[_str_len_];

	

	// Initialize PETSC
	PetscCall(PetscInitialize(&argc,&argv,(char *)0, help));

	// set default to be a forward run and overwrite it with input file options
	IOparam.use = _none_;

	// load and parse input file
	PetscCall(FBLoad(&fb));

	// set solver options
	PetscCall(setSolverOptions(fb));

	IOparam.fb = fb;

	PetscCall(getStringParam(IOparam.fb, _OPTIONAL_, "Adjoint_mode", str, "None"));

	if     (!strcmp(str, "None"))                   IOparam.use = _none_;
	else if(!strcmp(str, "GenericInversion"))       IOparam.use = _inversion_;
	else if(!strcmp(str, "AdjointGradients"))       IOparam.use = _adjointgradients_;
	else if(!strcmp(str, "GradientDescent"))        IOparam.use = _gradientdescent_;
	else if(!strcmp(str, "SyntheticForwardRun"))    IOparam.use = _syntheticforwardrun_;
	else
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Unknown parameter for 'Adjoint_mode'. Possibilities are [None; GenericInversion; AdjointGradients; GradientDescent or SyntheticForwardRun]");
	} 
	
	if(IOparam.use == _none_)
	{
		// forward simulation
		PetscCall(LaMEMLibMain(NULL, fb));
	}
	else
	{
		// inversion or adjoint gradient computation
		PetscCall(LaMEMAdjointMain(&IOparam));
	}

	// destroy file buffer
	PetscCall(FBDestroy(&fb));

	// cleanup PETSC
	PetscCall(PetscFinalize());

	return 0;
}
//--------------------------------------------------------------------------

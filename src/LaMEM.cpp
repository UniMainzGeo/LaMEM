/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2011-2018, JGU Mainz, Anton Popov, Boris Kaus
 **    All rights reserved.
 **
 **    This software was developed at:
 **
 **         Institute of Geosciences
 **         Johannes-Gutenberg University, Mainz
 **         Johann-Joachim-Becherweg 21
 **         55128 Mainz, Germany
 **
 **    project:    LaMEM
 **    filename:   LaMEM.cpp
 **
 **    LaMEM is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, version 3 of the License.
 **
 **    LaMEM is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with LaMEM. If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    Contact:
 **        Boris Kaus       [kaus@uni-mainz.de]
 **        Anton Popov      [popov@uni-mainz.de]
 **
 **
 **    Main development team:
 **         Anton Popov      [popov@uni-mainz.de]
 **         Boris Kaus       [kaus@uni-mainz.de]
 **         Tobias Baumann
 **			Georg Reuber
 **         Adina Pusok
 **         Arthur Bauville
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

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
	PetscErrorCode 	ierr;

	// Initialize PETSC
	ierr = PetscInitialize(&argc,&argv,(char *)0, help); CHKERRQ(ierr);

	ModParam IOparam;
	FB       *fb;
	char      str[_str_len_];

	// set default to be a forward run and overwrite it with input file options

	IOparam.use = _none_;
	ierr = FBLoad(&fb, PETSC_FALSE); CHKERRQ(ierr);
	ierr = getStringParam(fb, _OPTIONAL_, "use", str, NULL); CHKERRQ(ierr);
	if     (!strcmp(str, "none"))                   IOparam.use = _none_;
	else if(!strcmp(str, "inversion"))              IOparam.use = _inversion_;
	else if(!strcmp(str, "adjointgradients"))       IOparam.use = _adjointgradients_;
	else if(!strcmp(str, "gradientdescent"))        IOparam.use = _gradientdescent_;
	else if(!strcmp(str, "syntheticforwardrun"))    IOparam.use = _syntheticforwardrun_;
	
	if(IOparam.use == 0)
	{
		// Forward simulation
		ierr = LaMEMLibMain(NULL); CHKERRQ(ierr);
	}
	else
	{
		// Inversion or gradient computation
		ierr = LaMEMAdjointMain(&IOparam, fb); CHKERRQ(ierr);
	}

	// destroy file buffer
	ierr = FBDestroy(&fb); CHKERRQ(ierr);

	ierr = PetscMemzero(&IOparam, sizeof(ModParam)); CHKERRQ(ierr);

	// cleanup PETSC
	ierr = PetscFinalize(); CHKERRQ(ierr);

	return 0;
}

//--------------------------------------------------------------------------

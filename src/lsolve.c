/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2011-2015, JGU Mainz, Anton Popov, Boris Kaus
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
 **    filename:   lsolve.c
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
 **         Adina Pusok
 **         Arthur Bauville
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
//---------------------------------------------------------------------------
//...................   LINEAR PRECONDITIONER ROUTINES   ....................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "parsing.h"
#include "scaling.h"
#include "tssolve.h"
#include "fdstag.h"
#include "solVar.h"
#include "bc.h"
#include "JacRes.h"
#include "multigrid.h"
#include "matrix.h"
#include "lsolve.h"
//---------------------------------------------------------------------------
// * implement preconditioners in PETSc
// * add default solver options
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesSetFromOptions"
PetscErrorCode PCStokesSetFromOptions(PCStokes pc)
{
	PetscBool found;
	char      pname[MAX_NAME_LEN];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = PetscOptionsGetString(NULL, NULL,"-jp_type", pname, MAX_NAME_LEN, &found); CHKERRQ(ierr);

	if(found == PETSC_TRUE)
	{
		if(!strcmp(pname, "bf"))
		{
			PetscPrintf(PETSC_COMM_WORLD, " Preconditioner type            : block factorization\n");
            
			pc->type = _STOKES_BF_;
		}
		else if(!strcmp(pname, "mg"))
		{
			PetscPrintf(PETSC_COMM_WORLD, " Preconditioner type            : coupled Galerkin geometric multigrid\n");
			pc->type = _STOKES_MG_;
		}
		else if(!strcmp(pname, "user"))
		{
			PetscPrintf(PETSC_COMM_WORLD, " Preconditioner type            : user-defined\n");
			pc->type = _STOKES_USER_;
		}
		else SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER,"#Incorrect Jacobian preconditioner type: %s", pname);
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, " Preconditioner type            : user-defined\n");
		pc->type = _STOKES_USER_;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesCreate"
PetscErrorCode PCStokesCreate(PCStokes *p_pc, PMat pm)
{
	//========================================================================
	// create Stokes preconditioner context
	//========================================================================

	PCStokes pc;
	PMatType pm_type;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// allocate space
	ierr = PetscMalloc(sizeof(p_PCStokes), &pc); CHKERRQ(ierr);

	// clear object
	ierr = PetscMemzero(pc, sizeof(p_PCStokes)); CHKERRQ(ierr);

	// read options
	ierr = PCStokesSetFromOptions(pc); CHKERRQ(ierr);

	if(pc->type == _STOKES_BF_)
	{
		// Block Factorization
		pc->Create  = PCStokesBFCreate;
		pc->Setup   = PCStokesBFSetup;
		pc->Destroy = PCStokesBFDestroy;
		pc->Apply   = PCStokesBFApply;
		pm_type     = _BLOCK_;
	}
	else if(pc->type == _STOKES_MG_)
	{
		// Galerkin multigrid
		pc->Create  = PCStokesMGCreate;
		pc->Setup   = PCStokesMGSetup;
		pc->Destroy = PCStokesMGDestroy;
		pc->Apply   = PCStokesMGApply;
		pm_type     = _MONOLITHIC_;
	}
	else if(pc->type == _STOKES_USER_)
	{
		// user-defined
		pc->Create  = PCStokesUserCreate;
		pc->Setup   = PCStokesUserSetup;
		pc->Destroy = PCStokesUserDestroy;
		pc->Apply   = PCStokesUserApply;
		pm_type     = _MONOLITHIC_;
	}

	// check matrix type
	if(pm->type != pm_type) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Incorrect Stokes preconditioner matrix type used");

	// set matrix
	pc->pm = pm;

	// create preconditioner
	ierr = pc->Create(pc); CHKERRQ(ierr);

	// return preconditioner
	(*p_pc) = pc;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesSetup"
PetscErrorCode PCStokesSetup(PCStokes pc)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = pc->Setup(pc); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesDestroy"
PetscErrorCode PCStokesDestroy(PCStokes pc)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = pc->Destroy(pc); CHKERRQ(ierr);
	ierr = PetscFree(pc);   CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//........................... BLOCK FACTORIZATION ...........................
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesBFCreate"
PetscErrorCode PCStokesBFCreate(PCStokes pc)
{
	PC          vpc;
	PCStokesBF *bf;
	JacRes     *jr;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// allocate space
	ierr = PetscMalloc(sizeof(PCStokesBF), (void**)&bf); CHKERRQ(ierr);

	// clear object
	ierr = PetscMemzero(bf, sizeof(PCStokesBF)); CHKERRQ(ierr);

	// store context
	pc->data = (void*)bf;

	// read options
	ierr = PCStokesBFSetFromOptions(pc); CHKERRQ(ierr);

	// access context
	jr = pc->pm->jr;

	// create velocity solver
	ierr = KSPCreate(PETSC_COMM_WORLD, &bf->vksp); CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(bf->vksp,"vs_");    CHKERRQ(ierr);
	ierr = KSPSetFromOptions(bf->vksp);            CHKERRQ(ierr);

	// create & set velocity multigrid preconditioner
	if(bf->vtype == _VEL_MG_)
	{
		ierr = MGCreate(&bf->vmg, jr);           CHKERRQ(ierr);
		ierr = KSPGetPC(bf->vksp, &vpc);         CHKERRQ(ierr);
		ierr = PCSetType(vpc, PCSHELL);          CHKERRQ(ierr);
		ierr = PCShellSetContext(vpc, &bf->vmg); CHKERRQ(ierr);
		ierr = PCShellSetApply(vpc, MGApply);    CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesBFSetFromOptions"
PetscErrorCode PCStokesBFSetFromOptions(PCStokes pc)
{
	PCStokesBF *bf;

	PetscBool   flg;
	char        pname[MAX_NAME_LEN];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	bf = (PCStokesBF*)pc->data;

	// set factorization type
	ierr = PetscOptionsGetString(NULL, NULL,"-bf_type", pname, MAX_NAME_LEN, &flg); CHKERRQ(ierr);

	if(flg == PETSC_TRUE)
	{
		if(!strcmp(pname, "upper"))
		{
			PetscPrintf(PETSC_COMM_WORLD, " Block factorization type       : upper \n");

			bf->type = _UPPER_;
		}
		else if(!strcmp(pname, "lower"))
		{
			PetscPrintf(PETSC_COMM_WORLD, " Block factorization type       : lower \n");

			bf->type = _LOWER_;
		}
		else SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER,"#Incorrect block factorization type: %s", pname);
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, " Block factorization type       : upper \n");

		bf->type = _UPPER_;
	}

	// set velocity solver type
	ierr = PetscOptionsGetString(NULL, NULL,"-bf_vs_type", pname, MAX_NAME_LEN, &flg); CHKERRQ(ierr);

	if(flg == PETSC_TRUE)
	{
		if(!strcmp(pname, "mg"))
		{
			PetscPrintf(PETSC_COMM_WORLD, " Velocity preconditioner        : Galerkin geometric multigrid\n");

			bf->vtype = _VEL_MG_;
		}
		else if(!strcmp(pname, "user"))
		{
			PetscPrintf(PETSC_COMM_WORLD, " Velocity preconditioner        : user-defined\n");

			bf->vtype = _VEL_USER_;
		}
		else SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_USER,"#Incorrect velocity solver type: %s", pname);
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, " Velocity preconditioner        : user-defined\n");

		bf->vtype = _VEL_USER_;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesBFDestroy"
PetscErrorCode PCStokesBFDestroy(PCStokes pc)
{
	PCStokesBF *bf;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	bf = (PCStokesBF*)pc->data;

	ierr = KSPDestroy(&bf->vksp);  CHKERRQ(ierr);

	if(bf->vtype == _VEL_MG_)
	{
		ierr = MGDestroy(&bf->vmg); CHKERRQ(ierr);
	}

	ierr = PetscFree(bf); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesBFSetup"
PetscErrorCode PCStokesBFSetup(PCStokes pc)
{
	PCStokesBF *bf;
	PMatBlock  *P;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	bf = (PCStokesBF*)pc->data;
	P  = (PMatBlock*) pc->pm->data;

	ierr = KSPSetOperators(bf->vksp, P->Avv, P->Avv); CHKERRQ(ierr);

	if(bf->vtype == _VEL_MG_)
	{
		ierr = MGSetup(&bf->vmg, P->Avv); CHKERRQ(ierr);
	}

	ierr = KSPSetUp(bf->vksp); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesBFApply"
PetscErrorCode PCStokesBFApply(Mat JP, Vec r, Vec x)
{
	//======================================================================
	// r - residual vector      (input)
	// x - approximate solution (output)
	//======================================================================

	PCStokes    pc;
	PCStokesBF *bf;
	PMatBlock  *P;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	ierr = MatShellGetContext(JP, (void**)&pc); CHKERRQ(ierr);

	bf = (PCStokesBF*)pc->data;
	P  = (PMatBlock*) pc->pm->data;

	// extract residual blocks
	ierr = VecScatterBlockToMonolithic(P->rv, P->rp, r, SCATTER_REVERSE); CHKERRQ(ierr);

	if(bf->type == _UPPER_)
	{
		//=======================
		// BLOCK UPPER TRIANGULAR
		//=======================

		// Schur complement already contains negative sign (no negative sign here)
		ierr = MatMult(P->iS, P->rp, P->xp);     CHKERRQ(ierr); // xp = (S^-1)*rp

		ierr = MatMult(P->Avp, P->xp, P->wv);    CHKERRQ(ierr); // wv = Avp*xp

		ierr = VecAXPY(P->rv, -1.0, P->wv);      CHKERRQ(ierr); // rv = rv - wv

		ierr = KSPSolve(bf->vksp, P->rv, P->xv); CHKERRQ(ierr); // xv = (Avv^-1)*rv
	}
	else if(bf->type == _LOWER_)
	{
		//=======================
		// BLOCK LOWER TRIANGULAR
		//=======================

		ierr = KSPSolve(bf->vksp, P->rv, P->xv); CHKERRQ(ierr); // xv = (Avv^-1)*rv

		ierr = MatMult(P->Apv, P->xv, P->wp);    CHKERRQ(ierr); // wp = Apv*xv

		ierr = VecAXPY(P->rp, -1.0, P->wp);      CHKERRQ(ierr); // rp = rp - wp

		// Schur complement already contains negative sign (no negative sign here)
		ierr = MatMult(P->iS, P->rp, P->xp);     CHKERRQ(ierr); // xp = (S^-1)*rp
	}

	// compose approximate solution
	ierr = VecScatterBlockToMonolithic(P->xv, P->xp, x, SCATTER_FORWARD); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//....................... COUPLED GALERKIN MULTIGRID ........................
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesMGCreate"
PetscErrorCode PCStokesMGCreate(PCStokes pc)
{
	PCStokesMG *mg;
	JacRes     *jr;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// allocate space
	ierr = PetscMalloc(sizeof(PCStokesMG), (void**)&mg); CHKERRQ(ierr);

	// store context
	pc->data = (void*)mg;

	// access context
	jr = pc->pm->jr;

	// create context
	ierr = MGCreate(&mg->mg, jr); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesMGDestroy"
PetscErrorCode PCStokesMGDestroy(PCStokes pc)
{
	PCStokesMG *mg;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	mg = (PCStokesMG*)pc->data;

	ierr = MGDestroy(&mg->mg); CHKERRQ(ierr);
	ierr = PetscFree(mg);      CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesMGSetup"
PetscErrorCode PCStokesMGSetup(PCStokes pc)
{
	PCStokesMG *mg;
	PMatMono   *P;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// acces context
	mg = (PCStokesMG*)pc->data;
	P  = (PMatMono*)  pc->pm->data;

	ierr = MGSetup(&mg->mg, P->A); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesMGApply"
PetscErrorCode PCStokesMGApply(Mat JP, Vec x, Vec y)
{
	PCStokes    pc;
	PCStokesMG *mg;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = MatShellGetContext(JP, (void**)&pc); CHKERRQ(ierr);

	mg = (PCStokesMG*)pc->data;

	// apply multigrid preconditioner
	ierr = PCApply(mg->mg.pc, x, y); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//............................. USER-DEFINED ................................
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesUserCreate"
PetscErrorCode PCStokesUserCreate(PCStokes pc)
{
	PCStokesUser *user;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// allocate space
	ierr = PetscMalloc(sizeof(PCStokesUser), (void**)&user); CHKERRQ(ierr);

	// store context
	pc->data = (void*)user;

	// create user-defined preconditioner
	// attach index sets in case fieldsplit preconditioner needs to be used
	// set additional options (whatever with -jp_ prefix)
	ierr = PCCreate(PETSC_COMM_WORLD, &user->pc); CHKERRQ(ierr);
	ierr = PCSetOptionsPrefix(user->pc, "jp_");   CHKERRQ(ierr);
	ierr = PCStokesUserAttachIS(pc);              CHKERRQ(ierr);
	ierr = PCSetFromOptions(user->pc);            CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesUserAttachIS"
PetscErrorCode PCStokesUserAttachIS(PCStokes pc)
{
	PCStokesUser *user;
	JacRes       *jr;
	DOFIndex     *dof;
	PetscInt      st, lnv, lnp;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	user = (PCStokesUser*)pc->data;
	jr   =  pc->pm->jr;
	dof  = &jr->fs->dof;
	st   =  dof->st;
	lnv  =  dof->lnv;
	lnp  =  dof->lnp;

	// create index sets
	ierr = ISCreateStride(PETSC_COMM_WORLD, lnv, st,     1, &user->isv); CHKERRQ(ierr);
	ierr = ISCreateStride(PETSC_COMM_WORLD, lnp, st+lnv, 1, &user->isp); CHKERRQ(ierr);

	// this needs to be defined before index sets can be attached
	ierr = PCSetType(user->pc, PCFIELDSPLIT); CHKERRQ(ierr);

	// attach index sets
	ierr = PCFieldSplitSetIS(user->pc, "v", user->isv); CHKERRQ(ierr);
	ierr = PCFieldSplitSetIS(user->pc, "p", user->isp); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesUserDestroy"
PetscErrorCode PCStokesUserDestroy(PCStokes pc)
{
	PCStokesUser *user;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get context
	user = (PCStokesUser*)pc->data;

	// cleanup
	ierr = PCDestroy(&user->pc);  CHKERRQ(ierr);
	ierr = ISDestroy(&user->isv); CHKERRQ(ierr);
	ierr = ISDestroy(&user->isp); CHKERRQ(ierr);
	ierr = PetscFree(user);       CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesUserSetup"
PetscErrorCode PCStokesUserSetup(PCStokes pc)
{
	PCStokesUser *user;
	PMatMono     *P;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	user = (PCStokesUser*)pc->data;
	P    = (PMatMono*)    pc->pm->data;

	// compute preconditioner
	ierr = PCSetOperators(user->pc, P->A, P->A);  CHKERRQ(ierr);
	ierr = PCSetUp(user->pc);                     CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PCStokesUserApply"
PetscErrorCode PCStokesUserApply(Mat JP, Vec x, Vec y)
{
	PCStokes      pc;
	PCStokesUser *user;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = MatShellGetContext(JP, (void**)&pc); CHKERRQ(ierr);

	user = (PCStokesUser*)pc->data;

	// apply user-defined preconditioner
	ierr = PCApply(user->pc, x, y); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------


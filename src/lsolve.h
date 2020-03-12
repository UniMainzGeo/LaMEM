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
 **    filename:   lsolve.h
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
#ifndef __lsolve_h__
#define __lsolve_h__
//---------------------------------------------------------------------------
// Stokes preconditioner type
enum PCStokesType
{
	_STOKES_BF_,  // block factorization
	_STOKES_MG_,  // Galerkin multigrid
	_STOKES_USER_ // user-defined

};

//---------------------------------------------------------------------------

// Stokes preconditioner type
enum PCBFType
{
	_UPPER_,  // upper triangular factorization
	_LOWER_,  // lower triangular factorization
	_BFBT_ 	  // BFBT

};

//---------------------------------------------------------------------------

// velocity block preconditioner type (bf only)
enum PCVelType
{
	_VEL_MG_,  // Galerkin multigrid
	_VEL_USER_ // user-defined

};
//---------------------------------------------------------------------------

typedef struct _p_PCStokes *PCStokes;

typedef struct _p_PCStokes
{
	PCStokesType  type;
	PMat          pm;   // preconditioner matrix
	void         *data; // type-specific context

	// operations
	PetscErrorCode (*Create)  (PCStokes pc);
	PetscErrorCode (*Setup)   (PCStokes pc);
	PetscErrorCode (*Destroy) (PCStokes pc);
	PetscErrorCode (*Apply)   (Mat P, Vec x, Vec y);

} p_PCStokes;

// PCStokes - pointer to an opaque structure (to be used in declarations)
// sizeof(p_PCStokes) - size of the opaque structure

//---------------------------------------------------------------------------

PetscErrorCode PCStokesCreate(PCStokes *p_pc, PMat pm);

PetscErrorCode PCStokesSetFromOptions(PCStokes pc);

PetscErrorCode PCStokesSetup(PCStokes pc);

PetscErrorCode PCStokesDestroy(PCStokes pc);

//---------------------------------------------------------------------------

// Block Factorization preconditioner context
struct PCStokesBF
{
	PCVelType vtype; // velocity solver type
	KSP       vksp;  // velocity solver
	KSP 	  pksp;  // pressure solver
	MG        vmg;   // velocity multigrid context
	MG		  pmg; 	 // pressure multigrid context
	PCBFType  type;  // factorization type

};

//---------------------------------------------------------------------------

PetscErrorCode PCStokesBFCreate(PCStokes pc);

PetscErrorCode PCStokesBFSetFromOptions(PCStokes pc);

PetscErrorCode PCStokesBFDestroy(PCStokes pc);

PetscErrorCode PCStokesBFSetup(PCStokes pc);

PetscErrorCode PCStokesBFApply(Mat JP, Vec x, Vec y);

//---------------------------------------------------------------------------

// Galerkin multigrid preconditioner context
struct PCStokesMG
{
	MG mg; // coupled multigrid context

};

//---------------------------------------------------------------------------

PetscErrorCode PCStokesMGCreate(PCStokes pc);

PetscErrorCode PCStokesMGDestroy(PCStokes pc);

PetscErrorCode PCStokesMGSetup(PCStokes pc);

PetscErrorCode PCStokesMGApply(Mat JP, Vec x, Vec y);

//---------------------------------------------------------------------------

// User-defined
struct PCStokesUser
{
	PC pc;       // general preconditioner object
	IS isv, isp; // velocity and pressure index sets

} ;

//---------------------------------------------------------------------------

PetscErrorCode PCStokesUserCreate(PCStokes pc);

PetscErrorCode PCStokesUserAttachIS(PCStokes pc);

PetscErrorCode PCStokesUserDestroy(PCStokes pc);

PetscErrorCode PCStokesUserSetup(PCStokes pc);

PetscErrorCode PCStokesUserApply(Mat JP, Vec x, Vec y);

//---------------------------------------------------------------------------

#endif

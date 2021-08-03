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
 **    filename:   dike.h
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
//.................. MATERIAL PARAMETERS READING ROUTINES....................                                                                                                      
//---------------------------------------------------------------------------                                                                                                      
#ifndef __dike_h__
#define __dike_h__
//---------------------------------------------------------------------------   

struct FB; 
//struct JacRes;   // necessary????? try to remove
struct ConstEqCtx;
struct DBMat;
struct TSSol;
struct Scaling;

//---------------------------------------------------------------------------                                                                                                      
//.......................   Dike Parameters  .......................                                                                                                      
//---------------------------------------------------------------------------                                                                                                      
struct Dike
{
public:
  PetscInt    ID;        // dike ID
  PetscScalar Mf;        // amount of magma-accomodated extension in front of box 
  PetscScalar Mb;        // amount of magma-accommodated extension in back of box
  PetscInt PhaseID;      // associated material phase ID
  PetscScalar t0_dike;    // starting time for moving the dike
  PetscScalar t1_dike;   // end time for moving the dike
  PetscScalar v_dike;     // velocity with which the dike move

  PetscScalar dikeRHS;   // output, added divergence to RHS of continuity equation, should it be private? s
};

      
struct DBPropDike
{
  Scaling  *scal;                     // scaling parameters 
  //  DBMat    *dbm;  
  PetscInt numDike;                   // number of dikes
  Dike     matDike[_max_num_dike_];   // dike properties per dike ID
};

// create the dike strutures for read-in 
PetscErrorCode DBDikeCreate(DBPropDike *dbdike, DBMat *dbm, FB *fb, PetscBool PrintOutput);

// read in dike parameters
PetscErrorCode DBReadDike(DBPropDike *dbdike, DBMat *dbm, FB *fb, PetscBool PrintOutput);

// compute the added RHS of the dike for the continuity equation
PetscErrorCode GetDikeContr(ConstEqCtx *ctx, PetscScalar *phRat, PetscScalar &dikeRHS);

// compute the new locations of the dikes in case they move with a specified velocity
PetscErrorCode MovingDike(ConstEqCtx *ctx, TSSol *ts, PetscScalar &left_new, PetscScalar &right_new);

//---------------------------------------------------------------------------
#endif

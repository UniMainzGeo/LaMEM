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

//struct Scaling;
struct FB;
struct JacRes;
struct ModParam;

//---------------------------------------------------------------------------                                                                                                      
//.......................   Dike Parameters  .......................                                                                                                      
//---------------------------------------------------------------------------                                                                                                      
struct Dike
{
public:

        PetscInt    ID;   // dike ID
        PetscScalar Mf;   // amount of magma-accomodated extension in front of box 
        PetscScalar Mb;   // amount of magma-accommodated extension in back of box
  PetscInt Phase;         // associated material phase ID

  //private:  

    PetscScalar dikeRHS;
};


struct DBPropDike
{
  PetscInt numDike;                   // number of dikes
  Dike     matDike[_max_num_dike_];   // dike properties per dike ID
};

// read dike properties
PetscErrorCode DBDikeCreate(DBPropDike *dbdike, FB *fb, PetscBool PrintOutput);

// read-indike parameters
PetscErrorCode DBReadDike(DBPropDike *dbdike, FB *fb, PetscBool PrintOutput);

//---------------------------------------------------------------------------
#endif

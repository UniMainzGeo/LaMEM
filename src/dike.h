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
struct ConstEqCtx;
struct DBMat;
struct FDSTAG;
struct TSSol;
struct JacRes;
struct Controls;

//---------------------------------------------------------------------------       
//.......................   Dike Parameters  .......................                                                                                                      
//---------------------------------------------------------------------------

struct Dike
{
public:
  PetscInt    ID;        // dike ID
  PetscInt    dyndike;
  PetscScalar Mf;        // amount of magma-accomodated extension in front of box 
  PetscScalar Mb;        // amount of magma-accommodated extension in back of box
  PetscScalar Mc;        // amount of magma-accommodated extension in center of box
  PetscScalar y_Mc;      // location in y direction of Mc, if in x-direction x_Mc needs to be given or in z-direction z_Mc
  PetscScalar x_Mc;
  PetscScalar z_Mc;
  PetscInt PhaseID;      // associated material phase ID
  PetscInt PhaseTransID; // associated phase transition ID (necessary for moving dike)
  Vec devxx_mean;
  Vec dPm;
  Vec lthickness;
};
      
struct DBPropDike
{
  PetscInt numDike;                   // number of dikes
  Dike     matDike[_max_num_dike_];   // dike properties per dike ID
};

// create the dike strutures for read-in 
PetscErrorCode DBDikeCreate(DBPropDike *dbdike, DBMat *dbm, FB *fb, JacRes *jr, PetscBool PrintOutput);

// read in dike parameters
PetscErrorCode DBReadDike(DBPropDike *dbdike, DBMat *dbm, FB *fb, JacRes *jr, PetscBool PrintOutput);

// compute the added RHS of the dike for the continuity equation
PetscErrorCode GetDikeContr(ConstEqCtx *ctx, PetscScalar *phRat, PetscInt &Airphase, PetscScalar &dikeRHS, PetscScalar &y_c, PetscInt J);

// compute dike heat after Behn & Ito, 2008
PetscErrorCode Dike_k_heatsource(JacRes *jr,
                                Material_t *phases,
                                PetscScalar &Tc,
                                PetscScalar *phRat,          // phase ratios in the control volume
                                PetscScalar &k,
                                PetscScalar &rho_A,
                                PetscScalar &y_c,
                                PetscInt J); 

PetscErrorCode Locate_Dike_Zones(JacRes *jr);

//---------------------------------------------------------------------------
#endif

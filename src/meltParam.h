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
 **    filename:   meltParam.h
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
//...................   PARAMETERIZED MELT FRACTION   .......................
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Implementation of melt parameterization described in:
// A new parameterization of hydrous mantle melting
// Katz, Richard F.; Spiegelman, Marc; Langmuir, Charles H.
// Geochem. Geophys. Geosyst.  Vol. 4, No. 9, 1073  DOI 10.1029/2002GC000433
// 09 September 2003

// Available at http://www.ldeo.columbia.edu/~katz/meltParam/
//---------------------------------------------------------------------------

#ifndef __meltParam_h__
#define __meltParam_h__

#define UNUSED    1e20
#define MAXITS    60
#define X_ACC     0.00001
#define SIG(x,y) ( (y)>=0.0 ? fabs(x) : -fabs(x) )

//---------------------------------------------------------------------------
// melting parameter
//---------------------------------------------------------------------------
typedef struct melt_parameters_s {
  PetscScalar A1, A2, A3, B1, B2, B3, C1, C2, C3;
  PetscScalar r1, r2, beta1, beta2, K,gamma;
  PetscScalar D_water,chi1, chi2, lambda;
  PetscScalar Cp, DS;
} meltPar_Katz;

//---------------------------------------------------------------------------
//   private function prototypes inline function declare
//---------------------------------------------------------------------------
PetscScalar calcDT(PetscScalar P, PetscScalar X, PetscScalar F,meltPar_Katz *mp);
PetscScalar calcF (PetscScalar T, PetscScalar dT,PetscScalar P,PetscScalar Fcpx,meltPar_Katz *mp);
PetscScalar FX_bal(PetscScalar x1,PetscScalar x2,PetscScalar T,PetscScalar P,PetscScalar X,PetscScalar Fcpx,meltPar_Katz *mp);
PetscScalar FT_bal(PetscScalar x1,PetscScalar x2,PetscScalar T,PetscScalar P,PetscScalar X,PetscScalar M,meltPar_Katz *mp);
PetscScalar FZero (PetscScalar F, PetscScalar T, PetscScalar P,PetscScalar X,PetscScalar Fcpx, meltPar_Katz *mp);
PetscScalar HZero (PetscScalar F, PetscScalar T, PetscScalar P,PetscScalar X,PetscScalar M,meltPar_Katz *mp);

//---------------------------------------------------------------------------
// public function prototypes
//---------------------------------------------------------------------------
void        setMeltParamsToDefault_Katz(meltPar_Katz *mp);
PetscScalar MPgetFEquilib (PetscScalar P,PetscScalar T, PetscScalar X, PetscScalar M,meltPar_Katz *mp);
PetscScalar MPgetFReactive(PetscScalar P,PetscScalar T, PetscScalar Cf,PetscScalar M,meltPar_Katz *mp);
PetscScalar MPgetTEquilib (PetscScalar P,PetscScalar F, PetscScalar X, PetscScalar M,meltPar_Katz *mp);
PetscScalar MPgetFconsH   (PetscScalar P,PetscScalar Ti,PetscScalar X, PetscScalar M,PetscScalar *Tf,meltPar_Katz *mp);
PetscScalar MPgetTSolidus (PetscScalar P,PetscScalar X, meltPar_Katz *mp);

//---------------------------------------------------------------------------
//  default values of parameters from the paper, table 2
//  these are used when setMeltParamsToDefault() is called
//---------------------------------------------------------------------------
#define DEFAULT_A1        1085.7
#define DEFAULT_A2         132.9
#define DEFAULT_A3          -5.1
#define DEFAULT_B1        1475.0
#define DEFAULT_B2          80.0
#define DEFAULT_B3          -3.2
#define DEFAULT_C1        1780.0 
#define DEFAULT_C2          45.0
#define DEFAULT_C3          -2.0
#define DEFAULT_R1           0.5
#define DEFAULT_R2           0.08
#define DEFAULT_BETA1        1.5
#define DEFAULT_BETA2        1.5
#define DEFAULT_K           43.0
#define DEFAULT_GAMMA        0.75
#define DEFAULT_D_WATER      0.01
#define DEFAULT_CHI1         0.12
#define DEFAULT_CHI2         0.01
#define DEFAULT_LAMBDA       0.6
#define DEFAULT_CP        1000.0
#define DEFAULT_DS         300.0

//---------------------------------------------------------------------------

#endif

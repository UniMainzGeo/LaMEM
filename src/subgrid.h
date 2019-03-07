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
 **    filename:   subgrid.h
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
//.................   SUBGRID MARKER RESAMPLING ROUTINES   ..................
//---------------------------------------------------------------------------
#ifndef __subgrid_h__
#define __subgrid_h__
//---------------------------------------------------------------------------






struct AdvCtx;
/*
struct FDSTAG;
struct JacRes;
struct FreeSurf;
struct DBMat;
*/

//---------------------------------------------------------------------------

// resample markers
PetscErrorCode ADVMarkSubGrid(AdvCtx *actx);

/*

// exchange markers between the processors resulting from the position change
PetscErrorCode ADVExchange(AdvCtx *actx);

// project history INCREMENTS from grid to markers
PetscErrorCode ADVProjHistGridToMark(AdvCtx *actx);

// interpolate field history increments to markers
PetscErrorCode ADVInterpFieldToMark(AdvCtx *actx, InterpCase icase);

// update marker positions from current velocities & time step
PetscErrorCode ADVAdvectMark(AdvCtx *actx);

// count number of markers to be sent to each neighbor domain
PetscErrorCode ADVMapMarkToDomains(AdvCtx *actx);

// communicate number of markers with neighbor processes
PetscErrorCode ADVExchangeNumMark(AdvCtx *actx);

// create send and receive buffers for asynchronous MPI communication
PetscErrorCode ADVCreateMPIBuff(AdvCtx *actx);

// communicate markers with neighbor processes
PetscErrorCode ADVExchangeMark(AdvCtx *actx);

// store received markers, collect garbage
PetscErrorCode ADVCollectGarbage(AdvCtx *actx);

// free communication buffer
PetscErrorCode ADVDestroyMPIBuff(AdvCtx *actx);

// find host cells for local markers
PetscErrorCode ADVMapMarkToCells(AdvCtx *actx);

// creates arrays to optimize marker-cell interaction
PetscErrorCode ADVUpdateMarkCell(AdvCtx *actx);

// delete marker outflow
PetscErrorCode ADVMarkDeleteOutflow(AdvCtx *actx);

// change marker phase when crossing free surface
PetscErrorCode ADVMarkCrossFreeSurf(AdvCtx *actx);

// check marker phases
PetscErrorCode ADVCheckMarkPhases(AdvCtx *actx);

// update history variables without advection
PetscErrorCode ADVUpdateHistADVNone(AdvCtx *actx);

// get maximum inverse time step (CFL)
PetscErrorCode ADVSelectTimeStep(AdvCtx *actx, PetscInt *restart);

*/

//---------------------------------------------------------------------------
#endif

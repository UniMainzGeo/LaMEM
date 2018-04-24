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
 **    filename:   adjoint.c
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

/*
 * TODO:
 * - Pay attention local and global might still be inconsistent in Exchange Volume
 *
 * - So far only tested for exactly 1 extraction step
 * - What happens if the melt quantity actually decreases instead of increase? This might be not treated here yet
 * - What pressure and T should the new marker have?
 *   #The temperature is fixed and given as material parameter. The pressure must be interpolated from the surrounding.
 * - Two things are hard-coded (ctrl+F "hard-coded" to find them)
 * - How to handle the actual volume change in injection?:
 * 		1) Delete all markers from cell and inject equally distributed (and volume equal)
 *
 *      2) Reduce volume of existing markers in the cell such that volume is again correct?
 * 		3) Is it phyiscally  correct if it is just corrected for the inflated volumeas is done right now? <?>
 * 		4) Artificially massively increase the injected volume until it is correct? No
 *
 * 	POTENTIAL HELP:
 * 	- Implement a timestep criterium for the melt extraction. Only incrementally move melt up such that velocity doesn't run away.
 */

#include "LaMEM.h"
#include "cvi.h"
#include "phase.h"
#include "fdstag.h"
#include "constEq.h"
#include "advect.h"
#include "tools.h"
#include "JacRes.h"
#include "meltextraction.h"
#include "scaling.h"
#include "parsing.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MeltExtractionCreate"
PetscErrorCode MeltExtractionCreate(JacRes *jr, FB *fb)
{ // First functions that is called
	Scaling *scal;
	Material_t *mat;
	PetscInt maxPhaseID;


	PetscErrorCode ierr;
	PetscFunctionBegin;

	scal       = jr->scal;
    maxPhaseID = jr->dbm->numPhases-1;
    mat        = jr->dbm->phases;


	ierr = DMCreateGlobalVector(jr->fs->DA_CEN, &jr->gdMV)       ; CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(jr->fs->DA_CEN, &jr->gdMVmerge)  ; CHKERRQ(ierr);
	ierr = DMCreateLocalVector (jr->fs->DA_CEN, &jr->ldMV)       ; CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(jr->fs->DA_CEN, &jr->gdc)        ; CHKERRQ(ierr);
	ierr = DMCreateLocalVector (jr->fs->DA_CEN, &jr->ldc)        ; CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(jr->fs->DA_CEN, &jr->gdMoho)     ; CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(jr->fs->DA_CEN, &jr->gdMohomerge); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(jr->fs->DA_CEN, &jr->Miphase)    ; CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(jr->fs->DA_CEN, &jr->Vol); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "PhInt", &mat->PhInt,   1,  maxPhaseID); CHKERRQ(ierr);
	ierr = getIntParam   (fb,_OPTIONAL_,"PhExt",&mat->PhExt,1, maxPhaseID); CHKERRQ(ierr);





	PetscFunctionReturn(0);

}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MeltExtractionDestroy"
PetscErrorCode MeltExtractionDestroy(JacRes *jr)
{ // Last

	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = VecDestroy(&jr->gdMV);           CHKERRQ(ierr);
	ierr = VecDestroy(&jr->gdMVmerge);      CHKERRQ(ierr);
	ierr = VecDestroy(&jr->ldMV);           CHKERRQ(ierr);
	ierr = VecDestroy(&jr->gdc);            CHKERRQ(ierr);
	ierr = VecDestroy(&jr->ldc);            CHKERRQ(ierr);
	ierr = VecDestroy(&jr->ldMoho);         CHKERRQ(ierr);
	ierr = VecDestroy(&jr->Miphase);		CHKERRQ(ierr);
	ierr = VecDestroy(&jr->Vol);		CHKERRQ(ierr);


	PetscFunctionReturn(0);

}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MeltExtractionSave"
PetscErrorCode MeltExtractionSave(AdvCtx *actx, JacRes *jr)
{    /* 2:Create the structure data that must be used to inject new particles, and compute the associated sink and source term.Copy dMF (see Consteq.cpp)
        into a grid whose coordinate are based on the center of cell, and which save the variable from the bulk variables */

	PetscErrorCode ierr;
	PetscFunctionBegin;

	PetscScalar  ***ldMV,***Mipbuff,***Volume                  ;   // 3D structure storing the melt that has been extracted
	PetscInt      i, j, k, nx, ny, nz, sx, sy, sz, iter, iphase;   // Iteration
	SolVarBulk   *svBulk                                       ;   // pointer 2 the solution variable defined as "bulk cell properties"
	SolVarCell   *svCell                                       ;   // pointer 2 the Solution variable defined in the cell
	FDSTAG       *fs                                           ;   // pointer 2 the grid structure and variables
	PData        *pd                                           ;   // pointer 2 the phase diagram structure
	PetscInt     numPhases                                     ;   // number of phases
    PetscScalar  *phRat                                        ;   // Phase Ratio
    PetscScalar  ***p,pc                                       ;   // pressure
    PetscScalar  ***T,Tc                                       ;   // temperature
    PetscScalar  mfeff                                         ;   // Effective melt extracted
    Material_t   *mat                                          ;   // Material properties structure
    Material_t   *phases                                       ;   // Phases
    //AdvCtx       *actx                                         ;

    // Access to the context (?)
	fs        = jr->fs                                         ;   // take the structured grid data from jr. The out put is a pointer structure
	numPhases = jr->dbm->numPhases                             ;   // take the number of phases from dbm structures
    pd        = jr->Pd                                         ;   // take the structure associated to the phase diagram
    phases    = jr->dbm->phases                                ;   // take the phases

    //
    // Calling Moho_Tracking (working on)
    //
    // Initialize & get array

	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp,      &p)        ;      CHKERRQ(ierr);
    ierr = DMDAVecGetArray(fs->DA_CEN, jr->lT,      &T)        ;      CHKERRQ(ierr);

	iter = 0                                                   ;
	GET_CELL_RANGE(nx, sx, fs->dsx)                            ;
	GET_CELL_RANGE(ny, sy, fs->dsy)                            ;
	GET_CELL_RANGE(nz, sz, fs->dsz)                            ;
	ierr = VecZeroEntries(jr->Vol)                       ; CHKERRQ(ierr);


    for(iphase=0;iphase<numPhases;iphase++)
    {
      	  mat=&phases[iphase];
    	  if(mat->Pd_rho == 1)
    	  {
    	//Create the buffer
    	ierr = VecZeroEntries(jr->Miphase)                       ; CHKERRQ(ierr);
    	ierr = DMDAVecGetArray(fs->DA_CEN, jr->Miphase, &Mipbuff); CHKERRQ(ierr);
        START_STD_LOOP
          {
    		  svCell    = &jr->svCell[iter]                 ;     // take the central node based properties
     		  svBulk    = &svCell->svBulk                     ;     // take the bulk solution variables
    		  phRat     = jr->svCell[iter++].phRat            ;     // take phase ratio on the central node

        	if(phRat[iphase])
        	{
    	           // Temperature&Pressure
    	           // access current pressure
    		    	pc    = p[k][j][i]                                                  ;
    			   // current temperature
    			    Tc    = T[k][j][i]                                                  ;
    		        ierr  = SetDataPhaseDiagram(pd, pc, Tc, 0, mat->pdn); CHKERRQ(ierr) ;
    			    mfeff = pd->mf-svBulk->mfextot                                      ;// historical variables
    		        if (mfeff<0) mfeff=0                                                ;//Correction
    			    if( mfeff>phases[iphase].Mtrs)
    				  {
    			    	Mipbuff[k][j][i]   = - phRat[iphase] * (mfeff-phases[iphase].Mleft);
    				   }
    			    else
    			    {
    			    	Mipbuff[k][j][i]   = 0                                           ;
    			    }
        	   }

           }END_STD_LOOP
           ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->Miphase, &Mipbuff)                 ;        CHKERRQ(ierr);
           // Send the data to Melt Extraction Exchange volume & compute the injection

           // Update Miphase
           ierr =  MeltExtractionExchangeVolume(jr,iphase)	;	CHKERRQ(ierr);
           // Update the marker properties (Interpolate the properties back to the marker, then eventually inject)
           ierr =  MeltExtractionInterpMarker(actx,iphase)	;	CHKERRQ(ierr);		// Issue1: How to handle the extrusion&free surface?

           ierr = DMDAVecGetArray(fs->DA_CEN,jr->Vol,&Volume); CHKERRQ(ierr);
           ierr = DMDAVecGetArray(fs->DA_CEN,jr->Miphase,&Mipbuff); CHKERRQ(ierr);
            iter=0;
           START_STD_LOOP{
               Volume[k][j][i] +=Mipbuff[k][j][i] ;

           }END_STD_LOOP

           ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->Miphase, &Mipbuff)                 ;        CHKERRQ(ierr);
           ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->Vol, &Volume)                      ;        CHKERRQ(ierr);




           }
     }

    ierr = DMDAVecGetArray(fs->DA_CEN,jr->Vol,&Volume); CHKERRQ(ierr);

	PetscFunctionReturn(0);

}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MeltExtractionInterpMarker"
PetscErrorCode MeltExtractionInterpMarker(AdvCtx *actx, PetscInt iphase)
{
	//=======================================================================
	// interpolate increments of the history field of melt extraction to markers
	//=======================================================================

		FDSTAG      *fs ;
		JacRes      *jr ;
		Marker      *P  ;
		AdvVelCtx   vi  ;
		PetscScalar UP  ;
		PetscInt    nx, ny, sx, sy, sz;
		PetscInt    jj, ID, I, J, K   ;
		PetscScalar *vgdc, ***vldc, *Mipbuff ;

		PetscErrorCode ierr;
		PetscFunctionBegin;

		ierr = ADVelCreate(actx, &vi);  CHKERRQ(ierr);

		fs = actx->fs;
		jr = actx->jr;

		// starting indices & number of cells
		sx = fs->dsx.pstart; nx = fs->dsx.ncels;
		sy = fs->dsy.pstart; ny = fs->dsy.ncels;
		sz = fs->dsz.pstart;

        ierr = VecGetArray(jr->gdc, &vgdc);  CHKERRQ(ierr);
		// access 1D layouts of global vectors
    	ierr = VecGetArray(jr->Miphase, &Mipbuff); CHKERRQ(ierr);

        // Loop over the cell, to take dMF 
		for(jj = 0; jj < fs->nCells; jj++) vgdc[jj] = Mipbuff[jj]; // Save (again?) the dMF in a vector without takint into account k,j,i but the global indexing

		// restore access
		ierr = VecRestoreArray(jr->gdc, &vgdc);  CHKERRQ(ierr);
        ierr = VecRestoreArray(jr->Miphase, &Mipbuff)                 ;        CHKERRQ(ierr);
		// communicate boundary values
		GLOBAL_TO_LOCAL(fs->DA_CEN, jr->gdc, jr->ldc);   // Local/Global Grid and MPI communication routines

		// access 3D layouts of local vectors
		ierr = DMDAVecGetArray(fs->DA_CEN, jr->ldc, &vldc); CHKERRQ(ierr); // Create the local vector containing all the information related to the dM

		// scan ALL markers
		for(jj = 0; jj < actx->nummark; jj++)
		{
		    // access next marker
		    	P = &actx->markers[jj];

			// get consecutive index of the host cell
		    	ID = actx->cellnum[jj];

			// expand I, J, K cell indices
		    	GET_CELL_IJK(ID, I, J, K, nx, ny)

			// access buffer
		    	UP = vldc[sz+K][sy+J][sx+I];

		    	if(UP > 0)
		    	{
				ierr = MeltExtractionInject(jr,actx,&vi, ID, I, J, K, UP,iphase);  CHKERRQ(ierr);
				vldc[sz+K][sy+J][sx+I] = 0; // It avoid to repeat the injection.
		    	}
		    	else if(UP<0)
		    	{
		    		if(P->phase==iphase)
		    		{
		    		P->Mtot += -UP;  // has to increase by the amount of melt change
		    		P->Mvol +=  UP;  // has to decrease
		    		}
		    	}
		}

		// restore access
		ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->ldc, &vldc); CHKERRQ(ierr);

		ierr = ADVelDestroy(&vi);      CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//-----------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MeltExtractionInterpMarkerBackToGrid"
PetscErrorCode MeltExtractionInterpMarkerBackToGrid(AdvCtx *actx)
{
	FDSTAG      *fs;
	JacRes      *jr;
	Marker      *P;
	SolVarCell  *svCell;
	PetscInt     ID, I, J, K, II, JJ, KK;
	PetscInt     ii, jj, numPhases;
	PetscInt     nx, ny, sx, sy, sz, nCells;
	PetscScalar  xp, yp, zp, xc, yc, zc, wxn, wyn, wzn, wxc, wyc, wzc, w = 0.0;
	PetscScalar  UPXY, UPXZ, UPYZ;
	PetscScalar *gxy, *gxz, *gyz, ***lxy, ***lxz, ***lyz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs        = actx->fs;
	jr        = actx->jr;
	numPhases = actx->dbm->numPhases;

	// check marker phases
	ierr = ADVCheckMarkPhases(actx); CHKERRQ(ierr);

	//======
	// CELLS
	//======

	// marker-to-grid projection (cell nodes)

	// number of cells
	nx     = fs->dsx.ncels;
	ny     = fs->dsy.ncels;
	nCells = fs->nCells;

	// clear history variables
	for(jj = 0; jj < nCells; jj++)
	{
		// access solution variable
		svCell = &jr->svCell[jj];

		// clear phase ratios
		for(ii = 0; ii < numPhases; ii++) svCell->phRat[ii] = 0.0;

		// clear history variables
		svCell->svBulk.mfextot = 0.0;
		svCell->svBulk.mfVol  = 0.0;
	}

	// scan ALL markers
	for(jj = 0; jj < actx->nummark; jj++)
	{
		// access next marker
		P = &actx->markers[jj];

		// get consecutive index of the host cell
		ID = actx->cellnum[jj];

		// expand I, J, K cell indices
		GET_CELL_IJK(ID, I, J, K, nx, ny)

		// get marker coordinates
		xp = P->X[0];
		yp = P->X[1];
		zp = P->X[2];

		// get interpolation weights in cell control volumes
		wxc = WEIGHT_POINT_CELL(I, xp, fs->dsx);
		wyc = WEIGHT_POINT_CELL(J, yp, fs->dsy);
		wzc = WEIGHT_POINT_CELL(K, zp, fs->dsz);

		// get total interpolation weight
		w = wxc*wyc*wzc;

		// access solution variable of the host cell
		svCell = &jr->svCell[ID];

		// update phase ratios
		svCell->phRat[P->phase] += w;

		// update history variables
		svCell->svBulk.mfextot += w*P->Mtot;
		svCell->svBulk.mfVol   += w*P->Mvol;

	}

	// normalize interpolated values
	for(jj = 0; jj < nCells; jj++)
	{
		// access solution variable
		svCell = &jr->svCell[jj];

		// normalize phase ratios
		ierr = getPhaseRatio(numPhases, svCell->phRat, &w); CHKERRQ(ierr);

		// normalize history variables
		svCell->svBulk.mfextot /= w;
		svCell->svBulk.mfVol   /= w;
	}

	//======
	// EDGES
	//======

	// starting indices & number of cells
	sx = fs->dsx.pstart; nx = fs->dsx.ncels;
	sy = fs->dsy.pstart; ny = fs->dsy.ncels;
	sz = fs->dsz.pstart;

	// clear local vectors
	ierr = VecZeroEntries(jr->ldxy); CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->ldxz); CHKERRQ(ierr);
	ierr = VecZeroEntries(jr->ldyz); CHKERRQ(ierr);

	// access 3D layouts of local vectors
	ierr = DMDAVecGetArray(fs->DA_XY, jr->ldxy, &lxy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_XZ, jr->ldxz, &lxz); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_YZ, jr->ldyz, &lyz); CHKERRQ(ierr);

	// set interpolated fields to defaults
	UPXY = 1.0; UPXZ = 1.0; UPYZ = 1.0;

	// scan ALL markers
	for(jj = 0; jj < actx->nummark; jj++)
	{
		// access next marker
		P = &actx->markers[jj];

		// get consecutive index of the host cell
		ID = actx->cellnum[jj];

		// expand I, J, K cell indices
		GET_CELL_IJK(ID, I, J, K, nx, ny)

		// get marker coordinates
		xp = P->X[0];
		yp = P->X[1];
		zp = P->X[2];

		// get coordinates of cell center
		xc = fs->dsx.ccoor[I];
		yc = fs->dsy.ccoor[J];
		zc = fs->dsz.ccoor[K];

		// map marker on the control volumes of edge nodes
		if(xp > xc) { II = I+1; } else { II = I; }
		if(yp > yc) { JJ = J+1; } else { JJ = J; }
		if(zp > zc) { KK = K+1; } else { KK = K; }

		// get interpolation weights in cell control volumes
		wxc = WEIGHT_POINT_CELL(I, xp, fs->dsx);
		wyc = WEIGHT_POINT_CELL(J, yp, fs->dsy);
		wzc = WEIGHT_POINT_CELL(K, zp, fs->dsz);

		// get interpolation weights in node control volumes
		wxn = WEIGHT_POINT_NODE(II, xp, fs->dsx);
		wyn = WEIGHT_POINT_NODE(JJ, yp, fs->dsy);
		wzn = WEIGHT_POINT_NODE(KK, zp, fs->dsz);

		UPXY = P->Mtot;  UPXZ = P->Mtot;  UPYZ = P->Mtot;

		// update required fields from marker to edge nodes
		lxy[sz+K ][sy+JJ][sx+II] += wxn*wyn*wzc*UPXY;
		lxz[sz+KK][sy+J ][sx+II] += wxn*wyc*wzn*UPXZ;
		lyz[sz+KK][sy+JJ][sx+I ] += wxc*wyn*wzn*UPYZ;
	}

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_XY, jr->ldxy, &lxy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_XZ, jr->ldxz, &lxz); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_YZ, jr->ldyz, &lyz); CHKERRQ(ierr);

	// assemble global vectors
	LOCAL_TO_GLOBAL(fs->DA_XY, jr->ldxy, jr->gdxy)
	LOCAL_TO_GLOBAL(fs->DA_XZ, jr->ldxz, jr->gdxz)
	LOCAL_TO_GLOBAL(fs->DA_YZ, jr->ldyz, jr->gdyz)

	// access 1D layouts of global vectors
	ierr = VecGetArray(jr->gdxy, &gxy);  CHKERRQ(ierr);
	ierr = VecGetArray(jr->gdxz, &gxz);  CHKERRQ(ierr);
	ierr = VecGetArray(jr->gdyz, &gyz);  CHKERRQ(ierr);

	// copy (normalized) data to the residual context
	for(jj = 0; jj < fs->nXYEdg; jj++) jr->svXYEdge[jj].svDev.mfextot = gxy[jj]/jr->svXYEdge[jj].ws;
	for(jj = 0; jj < fs->nXZEdg; jj++) jr->svXZEdge[jj].svDev.mfextot = gxz[jj]/jr->svXZEdge[jj].ws;
	for(jj = 0; jj < fs->nYZEdg; jj++) jr->svYZEdge[jj].svDev.mfextot = gyz[jj]/jr->svYZEdge[jj].ws;

	// restore access
	ierr = VecRestoreArray(jr->gdxy, &gxy); CHKERRQ(ierr);
	ierr = VecRestoreArray(jr->gdxz, &gxz); CHKERRQ(ierr);
	ierr = VecRestoreArray(jr->gdyz, &gyz); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MeltExtractionExchangeVolume"
PetscErrorCode MeltExtractionExchangeVolume(JacRes *jr, PetscInt iphase)
{
	    FDSTAG      *fs                                                          ;
		Discret1D   *dsz                                                         ;
		SolVarCell  *svCell                                                      ;
		Controls    *ctrl                                                        ;
		PetscInt     i, j, k, K, sx, sy, sz, nx, ny, nz, iter,L                  ;
		Vec          dgmvvec, dgmvvecmerge                                       ;
		PetscScalar  bz, ez                                                      ;
		PetscScalar  level, IR                                                   ;
		PetscScalar  *vdgmvvec, *vdgmvvecmerge, ***vdgmvvecmerge2, ***vdgmvvec2, ***vdgmv, ***Mipbuff ;
		Material_t   *phases                                                     ;         // Phases
		PetscErrorCode ierr;
		PetscFunctionBegin;

		// access context
		fs     = jr->fs                       ;
		dsz    = &fs->dsz                     ;
		ctrl   = &jr->ctrl                    ;
		L      = (PetscInt)fs->dsz.rank       ; // rank of the processor
		phases = jr->dbm->phases              ;
		level  = phases[iphase].DInt/jr->scal->length; // It has to be modified
		IR     = phases[iphase].RelInt        ; // Amount of intrusion that has to be injected within the crust, at level

		// get local coordinate bounds
		ierr = FDSTAGGetLocalBox(fs, NULL, NULL, &bz, NULL, NULL, &ez); CHKERRQ(ierr);

		// create column communicator
		ierr = Discret1DGetColumnComm(dsz); CHKERRQ(ierr);

		ierr = DMDAVecGetArray  (fs->DA_CEN, jr->Miphase, &Mipbuff)   ; CHKERRQ(ierr);
		ierr = DMGetGlobalVector(jr->DA_CELL_2D, &dgmvvec)            ; CHKERRQ(ierr);
		ierr = DMGetGlobalVector(jr->DA_CELL_2D, &dgmvvecmerge)       ; CHKERRQ(ierr);
		ierr = VecZeroEntries   (dgmvvec)                             ; CHKERRQ(ierr);
		ierr = VecZeroEntries   (dgmvvecmerge)                        ; CHKERRQ(ierr);
		// scan all local cells
		GET_CELL_RANGE(nx, sx, fs->dsx)
		GET_CELL_RANGE(ny, sy, fs->dsy)
		GET_CELL_RANGE(nz, sz, fs->dsz)
		ierr = FDSTAGGetLocalBox(fs, NULL, NULL, &bz, NULL, NULL, &ez) ; CHKERRQ(ierr);
		ierr = DMDAVecGetArray(jr->DA_CELL_2D, dgmvvec, &vdgmvvec2)    ; CHKERRQ(ierr);

		iter = 0 ;

		START_STD_LOOP
		{
			vdgmvvec2[L][j][i] += Mipbuff[k][j][i];
		}
		END_STD_LOOP
		ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->Miphase  , &Mipbuff)     ; CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dgmvvec, &vdgmvvec2)     ; CHKERRQ(ierr);

		if(dsz->nproc != 1 )
		{
			ierr = VecGetArray(dgmvvec, &vdgmvvec)           ; CHKERRQ(ierr);
			ierr = VecGetArray(dgmvvecmerge, &vdgmvvecmerge) ; CHKERRQ(ierr);

			ierr = MPI_Allreduce(vdgmvvec, vdgmvvecmerge, (PetscMPIInt)(nx*ny), MPIU_SCALAR, MPI_SUM, dsz->comm); CHKERRQ(ierr);

			ierr = VecRestoreArray(dgmvvec, &vdgmvvec); CHKERRQ(ierr);
			ierr = VecRestoreArray(dgmvvecmerge, &vdgmvvecmerge); CHKERRQ(ierr);
		}
		else
		{
			ierr = VecCopy(dgmvvec,dgmvvecmerge);  CHKERRQ(ierr);
		}
		ierr = DMDAVecGetArray  (fs->DA_CEN, jr->Miphase, &Mipbuff)   ; CHKERRQ(ierr);
		ierr = DMDAVecGetArray  (jr->DA_CELL_2D, dgmvvecmerge, &vdgmvvecmerge2)   ; CHKERRQ(ierr);
		// scan all local cells
		iter=0;
		GET_CELL_RANGE(nx, sx, fs->dsx)
		GET_CELL_RANGE(ny, sy, fs->dsy)
		GET_CELL_RANGE(nz, sz, fs->dsz)
		START_PLANE_LOOP
		{
			if(vdgmvvecmerge2[L][j][i] < 0)
			{
				// check whether point belongs to domain
				if(level >= bz && level < ez)
				{
					// find containing cell
					K = FindPointInCell(dsz->ncoor, 0, dsz->ncels, level);

					// interpolate velocity
					Mipbuff[sz+K][j][i] = -IR*vdgmvvecmerge2[L][j][i];
					PetscPrintf(PETSC_COMM_WORLD, "Mipbuff=%.20f ",Mipbuff[sz+K][j][i]);
				}
			}
		}
		END_PLANE_LOOP
		ierr = DMDAVecRestoreArray(jr->DA_CELL_2D, dgmvvecmerge, &vdgmvvecmerge2); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->Miphase ,    &Mipbuff)        ; CHKERRQ(ierr);
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MeltExtractionInject"
PetscErrorCode MeltExtractionInject(JacRes *jr,AdvCtx *actx, AdvVelCtx *vi, PetscInt ID, PetscInt I, PetscInt J, PetscInt K, PetscScalar UP, PetscInt iphase)
{ /* 3_b
   */

	PetscInt    jj, ipn, n, ninj, pind, found, PhInject,ii,sind=0;
	PetscScalar xs[3], xe[3], xp[3],*X,dx[3];
	PetscScalar x,sumind=0.0;
	PetscScalar cf_rand;
	PetscRandom    rctx;
	Marker      *P;
	FDSTAG      *fs;
	Material_t  *phases;
	PetscErrorCode ierr;
	Scaling    *scal;
	PetscFunctionBegin;

	found = 0;
	fs        = actx->fs;
    phases    = jr->dbm->phases                                ;   // take the phases

	// ATTENTION hard-coded
	PhInject = phases[iphase].PhInt;
	// ATTENTION hard-coded


	// get markers in cell
	n = vi->markstart[ID+1] - vi->markstart[ID];

	// scan cell markers
	for(jj = 0; jj < n; jj++)
	{
		// get marker index
		pind = vi->markind[vi->markstart[ID] + jj];

		P = &actx->markers[pind];

		if(P->phase == PhInject && UP > 0)
		{
			UP = UP - (1-P->Mvol);
			if (UP < 0)
			{
				P->Mvol += UP + (1-P->Mvol);   // has to increase
				UP = 0;
			}
			else
			{
				P->Mvol += (1-P->Mvol);   // has to increase
			}
			found = 1;
		}
	}

	// We have not found a marker of the correct phase or there is still melt to be injected
	if(found == 0 || UP > 0)
	{
		ninj = (PetscInt)ceil(UP);  // Amount of markers we have to inject

		// allocate memory for new markers
		actx->nrecv = ninj;
		ierr = PetscMalloc((size_t)actx->nrecv*sizeof(Marker), &actx->recvbuf); CHKERRQ(ierr);
		ierr = PetscMemzero(actx->recvbuf, (size_t)actx->nrecv*sizeof(Marker)); CHKERRQ(ierr);

		// initialize the random number generator
		ierr = PetscRandomCreate(PETSC_COMM_SELF, &rctx); CHKERRQ(ierr);
		ierr = PetscRandomSetFromOptions(rctx);            CHKERRQ(ierr);

		// get cell coordinates
		xs[0] = fs->dsx.ncoor[I]; xe[0] = fs->dsx.ncoor[I+1];
		xs[1] = fs->dsy.ncoor[J]; xe[1] = fs->dsy.ncoor[J+1];
		xs[2] = fs->dsz.ncoor[K]; xe[2] = fs->dsz.ncoor[K+1];

		for(ipn = 0; ipn<ninj; ipn++)
		{
			// create random coordinate within this cell
			ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
			xp[0] = (xe[0] - xs[0]) * cf_rand + xs[0];
			ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
			xp[1] = (xe[1] - xs[1]) * cf_rand + xs[1];
			ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
			xp[2] = (xe[2] - xs[2]) * cf_rand + xs[2];

			// calculate the closest (parent marker)
			for (ii = 0; ii < n; ii++)
			{
				X  = actx->markers[ii].X;
				dx[0] = X[0] - xp[0];
				dx[1] = X[1] - xp[1];
				dx[2] = X[2] - xp[2];
				x  = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];

				if (ii == 0)    { sumind = x; sind = ii; }
				if (x < sumind) { sumind = x; sind = ii; }
			}

			// create new marker
			// actx->recvbuf[ipn]      = markers[sind];

			// hard-coded new marker properties for debugging
			actx->recvbuf[ipn].phase = PhInject;
			actx->recvbuf[ipn].p = actx->markers[sind].p;
			actx->recvbuf[ipn].T = phases[iphase].TInt;
			actx->recvbuf[ipn].APS = 5;
			actx->recvbuf[ipn].Mtot = 0;

			actx->recvbuf[ipn].S.xx = 0;
			actx->recvbuf[ipn].S.xy = 0;
			actx->recvbuf[ipn].S.xz = 0;
			actx->recvbuf[ipn].S.yy = 0;
			actx->recvbuf[ipn].S.yz = 0;
			actx->recvbuf[ipn].S.zz = 0;

			if (UP > 1)
			{
				actx->recvbuf[ipn].Mvol = 2;
				UP -= 1;
			}
			else
			{
				actx->recvbuf[ipn].Mvol = 1+UP;
				UP = 0;
			}

			actx->recvbuf[ipn].U[0] = 0;
			actx->recvbuf[ipn].U[1] = 0;
			actx->recvbuf[ipn].U[2] = 0;
			actx->recvbuf[ipn].X[0] = xp[0];
			actx->recvbuf[ipn].X[1] = xp[1];
			actx->recvbuf[ipn].X[2] = xp[2];
		}

		// destroy random context
		ierr = PetscRandomDestroy(&rctx); CHKERRQ(ierr);

		// store new markers
		actx->ndel = 0;
		ierr = ADVCollectGarbage(actx); CHKERRQ(ierr);

		// compute host cells for all the markers
		ierr = ADVMapMarkToCells(actx); CHKERRQ(ierr);

		// update arrays for marker-cell interaction
		ierr = ADVUpdateMarkCell(actx); CHKERRQ(ierr);

		// print info
		// PetscPrintf(PETSC_COMM_WORLD,"Melt extraction injected %i markers.\n", ninj);

		// clear
		ierr = PetscFree(actx->recvbuf); CHKERRQ(ierr);
	}
	PetscFunctionReturn(0);
}
//----------------------------------------------------------------------------------------------------------------------------------------//
/*#undef __FUNCT__
#define __FUNCT__ "Moho_Tracking"
PetscErrorCode Moho_Tracking(JacRes *jr)
{
	// based on Melt exchange
	FDSTAG      *fs                                                          ;
	Discret1D   *dsz                                                         ;
	SolVarCell  *svCell                                                      ;
	Controls    *ctrl                                                        ;
	PetscInt     i, j, k, K, sx, sy, sz, nx, ny, nz, iter,L                  ;
	Vec          dgmvvec, dgmvvecmerge                                       ;
	PetscScalar  bz, ez                                                      ;
	PetscScalar  level, IR                                                   ;
	PetscScalar  *vdgmvvec, *vdgmvvecmerge, **vdgmvvecmerge2, ***vdgmvvec2, ***vdgmv, ***Mipbuff ;
	Material_t   *phases  ;
	PetscScalar  *phRat	;

	// access context
	fs     = jr->fs	;
	dsz    = &fs->dsz	;
	ctrl   = &jr->ctrl	;
	L      = (PetscInt)fs->dsz.rank ; // rank of the processor
	phases = jr->dbm->phases;

	// Loop over the nodes to find where there is mantle phases

	// Loop to find threshold mantle phase = 0

	// Merge the information into a 2D grid








	PetscFunctionReturn(0);
}

*/

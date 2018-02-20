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

#include "LaMEM.h"
#include "JacRes.h"
#include "fdstag.h"
#include "constEq.h"
#include "phase.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "MeltExtraction"
PetscErrorCode MeltExtraction(JacRes *jr)
{
	FDSTAG     *fs;
	SolVarCell *svCell;
	SolVarEdge *svEdge;
	SolVarDev  *svDev;
	SolVarBulk *svBulk;
	Material_t *phases;

	PetscErrorCode ierr;
	PetscInt    i, j, k, ii, nx, ny, nz, sx, sy, sz, mx, my, mz, mcx, mcy, mcz;
	PetscInt    I1, I2, J1, J2, K1, K2;
	PetscScalar pShift;
	PetscInt    iter, numPhases;
	PetscScalar ***p, ***T;

	PetscFunctionBegin;

	fs = jr->fs;

	phases    =  jr->dbm->phases;    // phase parameters
	numPhases =  jr->dbm->numPhases; // number phases

	// Get some global solution
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lp,      &p);      CHKERRQ(ierr);
	ierr = DMDAVecGetArray(fs->DA_CEN, jr->lT,      &T);      CHKERRQ(ierr);

	//-------------------------------
	// central points
	//-------------------------------
	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// access solution variables
		svCell = &jr->svCell[iter++];
		svDev  = &svCell->svDev;
		svBulk = &svCell->svBulk;

		// scan all phases
		for(ii = 0; ii < numPhases; ii++)
		{
			// update present phases only
			if(svCell->phRat[ii])
			{
				// Get PD data
				if(phases[ii].Pd_rho == 1)
				{
					// Get the data from phase diagram
					SetDataPhaseDiagram(jr->Pd, p[k][j][i]-pShift, T[k][j][i], 0, phases[ii].pdn);  CHKERRQ(ierr);
					if (jr->Pd->mf > phases[ii].Mtrs)
					{
						jr->Pd->mfext    = jr->Pd->mf - phases[ii].Mleft;
						jr->Pd->mf       = phases[ii].Mleft;
						jr->Pd->mfextot += jr->Pd->mfext;
					}
					svBulk->mf       += svCell->phRat[ii]*jr->Pd->mf;
					svBulk->mfext    += svCell->phRat[ii]*jr->Pd->mfext;
					svBulk->mfextot += svCell->phRat[ii]*jr->Pd->mfextot;
					svDev->mf         = jr->Pd->mf;
					svDev->mfext      = jr->Pd->mfext;
					svDev->mfextot    = jr->Pd->mfextot;
				}
			}
		}
	}
	END_STD_LOOP

	//-------------------------------
	// xy edge points
	//-------------------------------
	iter = 0;
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_CELL_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// access solution variables
		svEdge = &jr->svXYEdge[iter++];
		svDev  = &svEdge->svDev;

		//=================
		// SECOND INVARIANT
		//=================

		// check index bounds
		I1 = i;   if(I1 == mx) I1--;
		I2 = i-1; if(I2 == -1) I2++;
		J1 = j;   if(J1 == my) J1--;
		J2 = j-1; if(J2 == -1) J2++;

		// scan all phases
		for(ii = 0; ii < numPhases; ii++)
		{
			// update present phases only
			if(svEdge->phRat[ii])
			{

				// Get PD data
				if(phases[ii].Pd_rho == 1)
				{
					// Get the data from phase diagram
					SetDataPhaseDiagram(jr->Pd, 0.25*(p[k][j][i] + p[k][j][i-1] + p[k][j-1][i] + p[k][j-1][i-1])-pShift, 0.25*(T[k][j][i] + T[k][j][i-1] + T[k][j-1][i] + T[k][j-1][i-1]), 0, phases[ii].pdn);  CHKERRQ(ierr);
					if (jr->Pd->mf > phases[ii].Mtrs)
					{
						jr->Pd->mfext    = jr->Pd->mf - phases[ii].Mleft;
						jr->Pd->mf       = phases[ii].Mleft;
						jr->Pd->mfextot += jr->Pd->mfext;
					}
					svDev->mf        = jr->Pd->mf;
				    svDev->mfext     = jr->Pd->mfext;
				    svDev->mfextot   = jr->Pd->mfextot;
				}
			}
		}
	}
	END_STD_LOOP

	//-------------------------------
	// xz edge points
	//-------------------------------
	iter = 0;
	GET_NODE_RANGE(nx, sx, fs->dsx)
	GET_CELL_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// access solution variables
		svEdge = &jr->svXZEdge[iter++];
		svDev  = &svEdge->svDev;

		//=================
		// SECOND INVARIANT
		//=================

		// check index bounds
		I1 = i;   if(I1 == mx) I1--;
		I2 = i-1; if(I2 == -1) I2++;
		K1 = k;   if(K1 == mz) K1--;
		K2 = k-1; if(K2 == -1) K2++;

		// scan all phases
		for(ii = 0; ii < numPhases; ii++)
		{
			// update present phases only
			if(svEdge->phRat[ii])
			{

				// Get PD data
				if(phases[ii].Pd_rho == 1)
				{
					// Get the data from phase diagram
					SetDataPhaseDiagram(jr->Pd, 0.25*(p[k][j][i] + p[k][j][i-1] + p[k-1][j][i] + p[k-1][j][i-1])-pShift, 0.25*(T[k][j][i] + T[k][j][i-1] + T[k-1][j][i] + T[k-1][j][i-1]), 0, phases[ii].pdn);  CHKERRQ(ierr);
					if (jr->Pd->mf > phases[ii].Mtrs)
					{
						jr->Pd->mfext    = jr->Pd->mf - phases[ii].Mleft;
						jr->Pd->mf       = phases[ii].Mleft;
						jr->Pd->mfextot += jr->Pd->mfext;
					}
					svDev->mf       = jr->Pd->mf;
				    svDev->mfext    = jr->Pd->mfext;
					svDev->mfextot  = jr->Pd->mfextot;
				}
			}
		}
	}
	END_STD_LOOP

	//-------------------------------
	// yz edge points
	//-------------------------------
	iter = 0;
	GET_CELL_RANGE(nx, sx, fs->dsx)
	GET_NODE_RANGE(ny, sy, fs->dsy)
	GET_NODE_RANGE(nz, sz, fs->dsz)

	START_STD_LOOP
	{
		// access solution variables
		svEdge = &jr->svYZEdge[iter++];
		svDev  = &svEdge->svDev;

		//=================
		// SECOND INVARIANT
		//=================

		// check index bounds
		J1 = j;   if(J1 == my) J1--;
		J2 = j-1; if(J2 == -1) J2++;
		K1 = k;   if(K1 == mz) K1--;
		K2 = k-1; if(K2 == -1) K2++;

		// scan all phases
		for(ii = 0; ii < numPhases; ii++)
		{
			// update present phases only
			if(svEdge->phRat[ii])
			{

				// Get PD data
				if(phases[ii].Pd_rho == 1)
				{
					// Get the data from phase diagram
					SetDataPhaseDiagram(jr->Pd, 0.25*(p[k][j][i] + p[k][j-1][i] + p[k-1][j][i] + p[k-1][j-1][i])-pShift, 0.25*(T[k][j][i] + T[k][j-1][i] + T[k-1][j][i] + T[k-1][j-1][i]), 0, phases[ii].pdn);  CHKERRQ(ierr);
					if (jr->Pd->mf > phases[ii].Mtrs)
					{
						jr->Pd->mfext    = jr->Pd->mf - phases[ii].Mleft;
						jr->Pd->mf       = phases[ii].Mleft;
						jr->Pd->mfextot += jr->Pd->mfext;
					}
					svDev->mf        = jr->Pd->mf;
					svDev->mfext     = jr->Pd->mfext;
				    svDev->mfextot   = jr->Pd->mfextot;
				}
			}
		}
	}
	END_STD_LOOP

	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lp,      &p);      CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(fs->DA_CEN, jr->lT,      &T);      CHKERRQ(ierr);


	PetscFunctionReturn(0);
}

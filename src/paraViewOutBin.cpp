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
 **    filename:   paraViewOutBin.c
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
//.................   FDSTAG PARAVIEW XML OUTPUT ROUTINES   .................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "paraViewOutBin.h"
#include "scaling.h"
#include "parsing.h"
#include "fdstag.h"
#include "JacRes.h"
#include "phase.h"
#include "outFunct.h"
#include "tools.h"
//---------------------------------------------------------------------------
// * phase-ratio output
// * integrate AVD phase viewer
//---------------------------------------------------------------------------
//............................. Output buffer ...............................
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "OutBufCreate"
PetscErrorCode OutBufCreate(OutBuf *outbuf, JacRes *jr)
{
	FDSTAG   *fs;
	PetscInt rx, ry, rz, sx, sy, sz, nx, ny, nz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	fs = jr->fs;

	// initialize parameters
	outbuf->fs = fs;
	outbuf->fp = NULL;
	outbuf->cn = 0;

	// get local output grid sizes
	GET_OUTPUT_RANGE(rx, nx, sx, fs->dsx)
	GET_OUTPUT_RANGE(ry, ny, sy, fs->dsy)
	GET_OUTPUT_RANGE(rz, nz, sz, fs->dsz)

	// allocate output buffer
	ierr = PetscMalloc((size_t)(_max_num_comp_*nx*ny*nz)*sizeof(float), &outbuf->buff); CHKERRQ(ierr);

	// set pointers to center, corner & edge buffers (reuse from JacRes object)
	outbuf->lbcen = jr->ldxx;
	outbuf->lbcor = jr->lbcor;
	outbuf->lbxy  = jr->ldxy;
	outbuf->lbxz  = jr->ldxz;
	outbuf->lbyz  = jr->ldyz;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "OutBufDestroy"
PetscErrorCode OutBufDestroy(OutBuf *outbuf)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// free output buffer
	ierr = PetscFree(outbuf->buff); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
void OutBufConnectToFile(OutBuf *outbuf, FILE *fp)
{
	// set file pointer
	outbuf->fp = fp;

	// clear buffer
	outbuf->cn = 0;
}
//---------------------------------------------------------------------------
void OutBufDump(OutBuf *outbuf)
{
	// dump output buffer contents to disk

	int nbytes;

	// compute number of bytes
	nbytes = (int)outbuf->cn*(int)sizeof(float);

	// dump number of bytes
	fwrite(&nbytes, sizeof(int), 1, outbuf->fp);

	// dump buffer contents
	fwrite(outbuf->buff, sizeof(float), (size_t)outbuf->cn, outbuf->fp);

	// clear buffer
	outbuf->cn = 0;
}
//---------------------------------------------------------------------------
void OutBufPutCoordVec(
	OutBuf      *outbuf,
	Discret1D   *ds,
	PetscScalar  cf)  // scaling coefficient
{
	// put FDSTAG coordinate vector to output buffer

	PetscInt    i, r, n, s;
	float       *buff;
	PetscScalar *ncoor;

	// get number of node points for output
	GET_OUTPUT_RANGE(r, n, s, (*ds))

	// access output buffer and coordinate array
	buff  = outbuf->buff;
	ncoor = ds->ncoor;

	// scale & copy to buffer
	for(i = 0; i < n; i++) buff[i] = (float) (cf*ncoor[i]);

	// update number of elements in the buffer
	outbuf->cn += n;

}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "OutBufPut3DVecComp"
PetscErrorCode OutBufPut3DVecComp(
	OutBuf      *outbuf,
	PetscInt     ncomp,  // number of components
	PetscInt     dir,    // component identifier
	PetscScalar  cf,     // scaling coefficient
	PetscScalar  shift)  // shift parameter (subtracted from scaled values)
{
	// put component of 3D vector to output buffer
	// component data is taken from obuf->gbcor vector

	FDSTAG      *fs;
	float       *buff;
	PetscScalar ***arr;
	PetscInt    i, j, k, rx, ry, rz, sx, sy, sz, nx, ny, nz, cnt;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access grid layout & buffer
	fs   = outbuf->fs;
	buff = outbuf->buff;

	// scatter ghost points to local buffer vector from global source vector
	LOCAL_TO_LOCAL(fs->DA_COR, outbuf->lbcor)

	// access local buffer vector
	ierr = DMDAVecGetArray(fs->DA_COR, outbuf->lbcor, &arr); CHKERRQ(ierr);

	// get sub-domain ranks, starting node IDs, and number of nodes
	GET_OUTPUT_RANGE(rx, nx, sx, fs->dsx)
	GET_OUTPUT_RANGE(ry, ny, sy, fs->dsy)
	GET_OUTPUT_RANGE(rz, nz, sz, fs->dsz)

	// set counter
	cnt = dir;

	// copy vector component to buffer
	if(cf < 0.0)
	{
		// negative scaling -> logarithmic output
		cf = -cf;

		START_STD_LOOP
		{
			// write
			buff[cnt] = (float) PetscLog10Real(cf*arr[k][j][i] - shift);

			// update counter
			cnt += ncomp;
		}
		END_STD_LOOP
	}
	else
	{
		// positive scaling -> standard output

		START_STD_LOOP
		{
			// write
			buff[cnt] = (float) (cf*arr[k][j][i] - shift);

			// update counter
			cnt += ncomp;
		}
		END_STD_LOOP

	}

	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_COR, outbuf->lbcor, &arr); CHKERRQ(ierr);

	// update number of elements in the buffer
	outbuf->cn += nx*ny*nz;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "OutBufZero3DVecComp"
PetscErrorCode OutBufZero3DVecComp(
	OutBuf      *outbuf,
	PetscInt     ncomp,  // number of components
	PetscInt     dir)    // component identifier
{
	// put zero component to output buffer

	FDSTAG      *fs;
	float       *buff;
	PetscInt    ii, nn, rx, ry, rz, sx, sy, sz, nx, ny, nz, cnt;

	PetscFunctionBegin;

	// access grid layout & buffer
	fs   = outbuf->fs;
	buff = outbuf->buff;

	// get sub-domain ranks, starting node IDs, and number of nodes
	GET_OUTPUT_RANGE(rx, nx, sx, fs->dsx)
	GET_OUTPUT_RANGE(ry, ny, sy, fs->dsy)
	GET_OUTPUT_RANGE(rz, nz, sz, fs->dsz)

	// set counter
	cnt = dir;

	// get total number of output nodes
	nn = nx*ny*nz;

	for(ii = 0; ii < nn; ii++)
	{
		buff[cnt] = 0.0;

		cnt += ncomp;

	}

	// update number of elements in the buffer
	outbuf->cn += nn;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//.......................... Vector output mask .............................
//---------------------------------------------------------------------------
void OutMaskSetDefault(OutMask *omask)
{
	omask->phase      = 1;
	omask->visc_total = 1;
	omask->visc_creep = 1;
	omask->visc_plast = 1;
	omask->velocity   = 1;
	omask->pressure   = 1;
}
//---------------------------------------------------------------------------
PetscInt OutMaskCountActive(OutMask *omask)
{
	PetscInt cnt = 0;

	if(omask->phase)          cnt++; // phase
	if(omask->density)        cnt++; // density
	if(omask->visc_total)     cnt++; // total effective viscosity
	if(omask->visc_creep)     cnt++; // creep effective viscosity
	if(omask->visc_plast)     cnt++; // viscoplastic viscosity
	if(omask->velocity)       cnt++; // velocity
	if(omask->pressure)       cnt++; // pressure
	if(omask->eff_press)      cnt++; // effective pressure
	if(omask->over_press)     cnt++; // overpressure
	if(omask->litho_press)    cnt++; // lithostatic pressure
	if(omask->pore_press)     cnt++; // pore pressure
	if(omask->temperature)    cnt++; // temperature
	if(omask->dev_stress)     cnt++; // deviatoric stress tensor
	if(omask->j2_dev_stress)  cnt++; // deviatoric stress second invariant
	if(omask->strain_rate)    cnt++; // deviatoric strain rate tensor
	if(omask->j2_strain_rate) cnt++; // deviatoric strain rate second invariant
	if(omask->melt_fraction)  cnt++; // melt fraction
	if(omask->fluid_density)  cnt++; // fluid density
	if(omask->vol_rate)       cnt++; // volumetric strain rate
	if(omask->vorticity)      cnt++; // vorticity vector
	if(omask->ang_vel_mag)    cnt++; // average angular velocity magnitude
	if(omask->tot_strain)     cnt++; // total strain
	if(omask->plast_strain)   cnt++; // accumulated plastic strain
	if(omask->plast_dissip)   cnt++; // plastic dissipation
	if(omask->tot_displ)      cnt++; // total displacements
	if(omask->SHmax)          cnt++; // maximum horizontal stress
	if(omask->EHmax)          cnt++; // maximum horizontal stress
	if(omask->yield)          cnt++; // yield stress
	if(omask->DIIdif)         cnt++; // diffusion creep relative strain rate
	if(omask->DIIdis)         cnt++; // dislocation creep relative strain rate
	if(omask->DIIprl)         cnt++; // Peierls creep relative strain rate

	// === debugging vectors ===============================================
	if(omask->moment_res)     cnt++; // momentum residual
	if(omask->cont_res)       cnt++; // continuity residual
	if(omask->energ_res)      cnt++; // energy residual

	// phase aggregates
	cnt += omask->num_agg;

	return cnt;
}
//---------------------------------------------------------------------------
//...................... ParaView output driver object ......................
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutCreate"
PetscErrorCode PVOutCreate(PVOut *pvout, FB *fb)
{
	OutMask *omask;
	PetscInt i, j, np, numPhases, maxPhaseID;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	omask      = &pvout->omask;
	numPhases  = pvout->jr->dbm->numPhases;
	maxPhaseID = numPhases-1;

	// initialize
	pvout->outpvd = 1;

	OutMaskSetDefault(omask);

	// read
	ierr = getStringParam(fb, _OPTIONAL_, "out_file_name",       pvout->outfile, "output");       CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_pvd",            &pvout->outpvd,            1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_phase",          &omask->phase,             1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_density",        &omask->density,           1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_visc_total",     &omask->visc_total,        1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_visc_creep",     &omask->visc_creep,        1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_visc_plast",     &omask->visc_plast,        1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_velocity",       &omask->velocity,          1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_pressure",       &omask->pressure,          1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_eff_press",      &omask->eff_press,         1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_over_press",     &omask->over_press,        1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_litho_press",    &omask->litho_press,       1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_pore_press",     &omask->pore_press,        1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_temperature",    &omask->temperature,       1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_dev_stress",     &omask->dev_stress,        1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_j2_dev_stress",  &omask->j2_dev_stress,     1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_strain_rate",    &omask->strain_rate,       1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_j2_strain_rate", &omask->j2_strain_rate,    1, 1); CHKERRQ(ierr);
//	ierr = getIntParam   (fb, _OPTIONAL_, "out_vol_rate",       &omask->vol_rate,          1, 1); CHKERRQ(ierr);
//	ierr = getIntParam   (fb, _OPTIONAL_, "out_vorticity",      &omask->vorticity,         1, 1); CHKERRQ(ierr);
//	ierr = getIntParam   (fb, _OPTIONAL_, "out_ang_vel_mag",    &omask->ang_vel_mag,       1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_shmax",          &omask->SHmax,             1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_ehmax",          &omask->EHmax,             1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_yield",          &omask->yield,             1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_rel_dif_rate",   &omask->DIIdif,            1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_rel_dis_rate",   &omask->DIIdis,            1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_rel_prl_rate",   &omask->DIIprl,            1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_tot_strain",     &omask->tot_strain,        1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_plast_strain",   &omask->plast_strain,      1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_plast_dissip",   &omask->plast_dissip,      1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_tot_displ",      &omask->tot_displ,         1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_moment_res",     &omask->moment_res,        1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_cont_res",       &omask->cont_res,          1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_energ_res",      &omask->energ_res,         1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_melt_fraction",  &omask->melt_fraction,     1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_fluid_density",  &omask->fluid_density,     1, 1); CHKERRQ(ierr);

	// read phase aggregates
	ierr = FBFindBlocks(fb, _OPTIONAL_, "<PhaseAggStart>", "<PhaseAggEnd>"); CHKERRQ(ierr);

	if(fb->nblocks > _max_num_phase_agg_)
	{
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many phase aggregates specified! Max allowed: %lld", (LLD)_max_num_phase_agg_);
	}

	omask->num_agg = fb->nblocks;

	for(i = 0; i < fb->nblocks; i++)
	{
		ierr = getStringParam(fb, _REQUIRED_, "name",      omask->agg_name[i],      NULL);            CHKERRQ(ierr);
		ierr = getIntParam   (fb, _REQUIRED_, "numPhase", &np,                      1,  numPhases);   CHKERRQ(ierr);
		ierr = getIntParam   (fb, _REQUIRED_, "phaseID",   omask->agg_phase_ID[i],  np, maxPhaseID);  CHKERRQ(ierr);

		// store number of phases
		omask->agg_num_phase[i] = np;

		fb->blockID++;

	}

	ierr = FBFreeBlocks(fb); CHKERRQ(ierr);

	// check
	if(!pvout->jr->ctrl.actTemp)             omask->energ_res = 0; // heat diffusion is deactivated
	if( pvout->jr->ctrl.gwType == _GW_NONE_) omask->eff_press = 0; // pore pressure is deactivated

	// print summary
	PetscPrintf(PETSC_COMM_WORLD, "Output parameters:\n");
	PetscPrintf(PETSC_COMM_WORLD, "   Output file name                        : %s \n", pvout->outfile);
	PetscPrintf(PETSC_COMM_WORLD, "   Write .pvd file                         : %s \n", pvout->outpvd ? "yes" : "no");

	if(omask->phase)          PetscPrintf(PETSC_COMM_WORLD, "   Phase                                   @ \n");
	if(omask->density)        PetscPrintf(PETSC_COMM_WORLD, "   Density                                 @ \n");
	if(omask->visc_total)     PetscPrintf(PETSC_COMM_WORLD, "   Total effective viscosity               @ \n");
	if(omask->visc_creep)     PetscPrintf(PETSC_COMM_WORLD, "   Creep effective viscosity               @ \n");
	if(omask->visc_plast)     PetscPrintf(PETSC_COMM_WORLD, "   Viscoplastic viscosity                  @ \n");
	if(omask->velocity)       PetscPrintf(PETSC_COMM_WORLD, "   Velocity                                @ \n");
	if(omask->pressure)       PetscPrintf(PETSC_COMM_WORLD, "   Pressure                                @ \n");
	if(omask->eff_press)      PetscPrintf(PETSC_COMM_WORLD, "   Effective pressure                      @ \n");
	if(omask->over_press)     PetscPrintf(PETSC_COMM_WORLD, "   Overpressure                            @ \n");
	if(omask->litho_press)    PetscPrintf(PETSC_COMM_WORLD, "   Lithostatic pressure                    @ \n");
	if(omask->pore_press)     PetscPrintf(PETSC_COMM_WORLD, "   Pore pressure                           @ \n");
	if(omask->temperature)    PetscPrintf(PETSC_COMM_WORLD, "   Temperature                             @ \n");
	if(omask->dev_stress)     PetscPrintf(PETSC_COMM_WORLD, "   Deviatoric stress tensor                @ \n");
	if(omask->j2_dev_stress)  PetscPrintf(PETSC_COMM_WORLD, "   Deviatoric stress second invariant      @ \n");
	if(omask->strain_rate)    PetscPrintf(PETSC_COMM_WORLD, "   Deviatoric strain rate tensor           @ \n");
	if(omask->j2_strain_rate) PetscPrintf(PETSC_COMM_WORLD, "   Deviatoric strain rate second invariant @ \n");
	if(omask->SHmax)          PetscPrintf(PETSC_COMM_WORLD, "   Maximum horizontal stress               @ \n");
	if(omask->EHmax)          PetscPrintf(PETSC_COMM_WORLD, "   Maximum horizontal extension            @ \n");
	if(omask->yield)          PetscPrintf(PETSC_COMM_WORLD, "   Yield stress                            @ \n");
	if(omask->DIIdif)         PetscPrintf(PETSC_COMM_WORLD, "   Diffusion creep relative strain rate    @ \n");
	if(omask->DIIdis)         PetscPrintf(PETSC_COMM_WORLD, "   Dislocation creep relative strain rate  @ \n");
	if(omask->DIIprl)         PetscPrintf(PETSC_COMM_WORLD, "   Peierls creep relative strain rate      @ \n");
	if(omask->tot_strain)     PetscPrintf(PETSC_COMM_WORLD, "   Accumulated Total Strain (ATS)          @ \n");
	if(omask->plast_strain)   PetscPrintf(PETSC_COMM_WORLD, "   Accumulated Plastic Strain (APS)        @ \n");
	if(omask->plast_dissip)   PetscPrintf(PETSC_COMM_WORLD, "   Plastic dissipation                     @ \n");
	if(omask->tot_displ)      PetscPrintf(PETSC_COMM_WORLD, "   Total displacements                     @ \n");
	if(omask->moment_res)     PetscPrintf(PETSC_COMM_WORLD, "   Momentum residual                       @ \n");
	if(omask->cont_res)       PetscPrintf(PETSC_COMM_WORLD, "   Continuity residual                     @ \n");
	if(omask->energ_res)      PetscPrintf(PETSC_COMM_WORLD, "   energy residual                         @ \n");
	if(omask->melt_fraction)  PetscPrintf(PETSC_COMM_WORLD, "   Melt fraction                           @ \n");
	if(omask->fluid_density)  PetscPrintf(PETSC_COMM_WORLD, "   Fluid density                           @ \n");

	for(i = 0; i < omask->num_agg; i++)
	{
		PetscPrintf(PETSC_COMM_WORLD, "   Aggregate: < %s >   Phases: < ", omask->agg_name[i]);

		for(j = 0; j < omask->agg_num_phase[i]; j++)
		{
			PetscPrintf(PETSC_COMM_WORLD, "%lld ", (LLD)omask->agg_phase_ID[i][j]);
		}

		PetscPrintf(PETSC_COMM_WORLD, ">\n");
	}

	PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");

	// count active output vectors
	pvout->nvec = OutMaskCountActive(omask);

	// create output buffer and vectors
	ierr = PVOutCreateData(pvout); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutCreateData"
PetscErrorCode PVOutCreateData(PVOut *pvout)
{
	JacRes   *jr;
	OutBuf   *outbuf;
	Scaling  *scal;
	OutMask  *omask;
	PetscInt  i, iter;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	jr     =  pvout->jr;
	outbuf = &pvout->outbuf;
	omask  = &pvout->omask;
	scal   =  jr->scal;
	iter   =  0;

	// create output buffer
	ierr = OutBufCreate(&pvout->outbuf, jr); CHKERRQ(ierr);

	// create vectors
	ierr = PetscMalloc(sizeof(OutVec)*(size_t)pvout->nvec, &pvout->outvecs); CHKERRQ(ierr);
	ierr = PetscMemzero(pvout->outvecs, sizeof(OutVec)*(size_t)pvout->nvec); CHKERRQ(ierr);

	if(omask->phase)          OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "phase",          scal->lbl_unit,             &PVOutWritePhase,        1, NULL);
	if(omask->density)        OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "density",        scal->lbl_density,          &PVOutWriteDensity,      1, NULL);
	if(omask->visc_total)     OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "visc_total",     scal->lbl_viscosity,        &PVOutWriteViscTotal,    1, NULL);
	if(omask->visc_creep)     OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "visc_creep",     scal->lbl_viscosity,        &PVOutWriteViscCreep,    1, NULL);
	if(omask->visc_plast)     OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "visc_plast",     scal->lbl_viscosity,        &PVOutWriteViscoPlastic, 1, NULL);
	if(omask->velocity)       OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "velocity",       scal->lbl_velocity,         &PVOutWriteVelocity,     3, NULL);
	if(omask->pressure)       OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "pressure",       scal->lbl_stress,           &PVOutWritePressure,     1, NULL);
	if(omask->eff_press)      OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "eff_press",      scal->lbl_stress,           &PVOutWriteEffPress,     1, NULL);
	if(omask->over_press)     OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "over_press",     scal->lbl_stress,           &PVOutWriteOverPress,    1, NULL);
	if(omask->litho_press)    OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "litho_press",    scal->lbl_stress,           &PVOutWriteLithoPress,   1, NULL);
	if(omask->pore_press)     OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "pore_press",     scal->lbl_stress,           &PVOutWritePorePress,    1, NULL);
	if(omask->temperature)    OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "temperature",    scal->lbl_temperature,      &PVOutWriteTemperature,  1, NULL);
	if(omask->dev_stress)     OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "dev_stress",     scal->lbl_stress,           &PVOutWriteDevStress,    9, NULL);
	if(omask->strain_rate)    OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "strain_rate",    scal->lbl_strain_rate,      &PVOutWriteStrainRate,   9, NULL);
	if(omask->j2_dev_stress)  OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "j2_dev_stress",  scal->lbl_stress,           &PVOutWriteJ2DevStress,  1, NULL);
	if(omask->j2_strain_rate) OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "j2_strain_rate", scal->lbl_strain_rate,      &PVOutWriteJ2StrainRate, 1, NULL);
	if(omask->vol_rate)       OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "vol_rate",       scal->lbl_strain_rate,      &PVOutWriteVolRate,      1, NULL);
	if(omask->vorticity)      OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "vorticity",      scal->lbl_strain_rate,      &PVOutWriteVorticity,    3, NULL);
	if(omask->ang_vel_mag)    OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "ang_vel_mag",    scal->lbl_angular_velocity, &PVOutWriteAngVelMag,    1, NULL);
	if(omask->tot_strain)     OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "tot_strain",     scal->lbl_unit,             &PVOutWriteTotStrain,    1, NULL);
	if(omask->plast_strain)   OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "plast_strain",   scal->lbl_unit,             &PVOutWritePlastStrain,  1, NULL);
	if(omask->plast_dissip)   OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "plast_dissip",   scal->lbl_dissipation_rate, &PVOutWritePlastDissip,  1, NULL);
	if(omask->tot_displ)      OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "tot_displ",      scal->lbl_length,           &PVOutWriteTotDispl,     3, NULL);
	if(omask->SHmax)          OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "SHmax",          scal->lbl_unit,             &PVOutWriteSHmax,        3, NULL);
	if(omask->EHmax)          OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "EHmax",          scal->lbl_unit,             &PVOutWriteEHmax,        3, NULL);
	if(omask->yield)          OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "yield",          scal->lbl_stress,           &PVOutWriteYield,        1, NULL);
	if(omask->DIIdif)         OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "rel_dif_rate",   scal->lbl_unit,             &PVOutWriteRelDIIdif,    1, NULL);
	if(omask->DIIdis)         OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "rel_dis_rate",   scal->lbl_unit,             &PVOutWriteRelDIIdis,    1, NULL);
	if(omask->DIIprl)         OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "rel_prl_rate",   scal->lbl_unit,             &PVOutWriteRelDIIprl,    1, NULL);
	// === debugging vectors ===============================================
	if(omask->melt_fraction)  OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "melt_fraction",  scal->lbl_unit,             &PVOutWriteMeltFraction, 1, NULL);
	if(omask->fluid_density)  OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "fluid_density",  scal->lbl_density,	      &PVOutWriteFluidDensity, 1, NULL);
	if(omask->moment_res)     OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "moment_res",     scal->lbl_volumetric_force, &PVOutWriteMomentRes,    3, NULL);
	if(omask->cont_res)       OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "cont_res",       scal->lbl_strain_rate,      &PVOutWriteContRes,      1, NULL);
	if(omask->energ_res)      OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, "energ_res",      scal->lbl_dissipation_rate, &PVOutWritEnergRes,      1, NULL);

	// setup phase aggregate output vectors
	for(i = 0; i < omask->num_agg; i++)
	{
		OutVecCreate(&pvout->outvecs[iter++], jr, outbuf, omask->agg_name[i], scal->lbl_unit, &PVOutWritePhaseAgg, omask->agg_num_phase[i], omask->agg_phase_ID[i]);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutDestroy"
PetscErrorCode PVOutDestroy(PVOut *pvout)
{

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// output vectors
	PetscFree(pvout->outvecs);

	// output buffer
	ierr = OutBufDestroy(&pvout->outbuf); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteTimeStep"
PetscErrorCode PVOutWriteTimeStep(PVOut *pvout, const char *dirName, PetscScalar ttime)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// update .pvd file if necessary
	ierr = UpdatePVDFile(dirName, pvout->outfile, "pvtr", &pvout->offset, ttime, pvout->outpvd); CHKERRQ(ierr);

	// write parallel data .pvtr file
	ierr = PVOutWritePVTR(pvout, dirName); CHKERRQ(ierr);

	// write sub-domain data .vtr files
	ierr = PVOutWriteVTR(pvout, dirName); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWritePVTR"
PetscErrorCode PVOutWritePVTR(PVOut *pvout, const char *dirName)
{
	FILE        *fp;
	FDSTAG      *fs;
	char        *fname;
	OutVec      *outvecs;
	PetscInt     i, rx, ry, rz;
	PetscMPIInt  nproc, iproc;

	PetscFunctionBegin;

	// only first process generates this file (WARNING! Bottleneck!)
	if(!ISRankZero(PETSC_COMM_WORLD)) PetscFunctionReturn(0);

	// access staggered grid layout
	fs = pvout->outbuf.fs;

	// open outfile.pvtr file in the output directory (write mode)
	asprintf(&fname, "%s/%s.pvtr", dirName, pvout->outfile);
	fp = fopen(fname,"w");
	if(fp == NULL) SETERRQ1(PETSC_COMM_SELF, 1,"cannot open file %s", fname);
	free(fname);

	// write header
	WriteXMLHeader(fp, "PRectilinearGrid");

	// open rectilinear grid data block (write total grid size)
	fprintf(fp, "\t<PRectilinearGrid GhostLevel=\"0\" WholeExtent=\"%lld %lld %lld %lld %lld %lld\">\n",
		1LL, (LLD)fs->dsx.tnods,
		1LL, (LLD)fs->dsy.tnods,
		1LL, (LLD)fs->dsz.tnods);

	// write cell data block (empty)
	fprintf(fp, "\t\t<PCellData>\n");
	fprintf(fp, "\t\t</PCellData>\n");

	// write coordinate block
	fprintf(fp, "\t\t<PCoordinates>\n");
	fprintf(fp, "\t\t\t<PDataArray type=\"Float32\" Name=\"Coordinates_X\" NumberOfComponents=\"1\" format=\"appended\"/>\n");
	fprintf(fp, "\t\t\t<PDataArray type=\"Float32\" Name=\"Coordinates_Y\" NumberOfComponents=\"1\" format=\"appended\"/>\n");
	fprintf(fp, "\t\t\t<PDataArray type=\"Float32\" Name=\"Coordinates_Z\" NumberOfComponents=\"1\" format=\"appended\"/>\n");
	fprintf(fp, "\t\t</PCoordinates>\n");

	// write description of output vectors (parameterized)
	outvecs = pvout->outvecs;
	fprintf(fp, "\t\t<PPointData>\n");
	for(i = 0; i < pvout->nvec; i++)
	{	fprintf(fp,"\t\t\t<PDataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"%lld\" format=\"appended\"/>\n",
			outvecs[i].name, (LLD)outvecs[i].ncomp);
	}
	fprintf(fp, "\t\t</PPointData>\n");

	// get total number of sub-domains
	MPI_Comm_size(PETSC_COMM_WORLD, &nproc);

	// write local grid sizes (extents) and data file names for all sub-domains
	for(iproc = 0; iproc < nproc; iproc++)
	{
		// get sub-domain ranks in all coordinate directions
		getLocalRank(&rx, &ry, &rz, iproc, fs->dsx.nproc, fs->dsy.nproc);

		// write data
		fprintf(fp, "\t\t<Piece Extent=\"%lld %lld %lld %lld %lld %lld\" Source=\"%s_p%1.8lld.vtr\"/>\n",
			(LLD)(fs->dsx.starts[rx] + 1), (LLD)(fs->dsx.starts[rx+1] + 1),
			(LLD)(fs->dsy.starts[ry] + 1), (LLD)(fs->dsy.starts[ry+1] + 1),
			(LLD)(fs->dsz.starts[rz] + 1), (LLD)(fs->dsz.starts[rz+1] + 1), pvout->outfile, (LLD)iproc);
	}

	// close rectilinear grid data block
	fprintf(fp, "\t</PRectilinearGrid>\n");
	fprintf(fp, "</VTKFile>\n");

	// close file
	fclose(fp);
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteVTR"
PetscErrorCode PVOutWriteVTR(PVOut *pvout, const char *dirName)
{
	FILE          *fp;
	FDSTAG        *fs;
	JacRes        *jr;
	char          *fname;
	OutBuf        *outbuf;
	OutVec        *outvecs;
	PetscInt       i, rx, ry, rz, sx, sy, sz, nx, ny, nz;
	PetscMPIInt    rank;
	size_t         offset = 0;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get global sub-domain rank
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

	// access output buffer object & staggered grid layout
	outbuf = &pvout->outbuf;
	fs     =  outbuf->fs;
	jr     =  pvout->jr;

	// get sizes of output grid
	GET_OUTPUT_RANGE(rx, nx, sx, fs->dsx)
	GET_OUTPUT_RANGE(ry, ny, sy, fs->dsy)
	GET_OUTPUT_RANGE(rz, nz, sz, fs->dsz)

	// open outfile_p_XXXXXX.vtr file in the output directory (write mode)
	asprintf(&fname, "%s/%s_p%1.8lld.vtr", dirName, pvout->outfile, (LLD)rank);
	fp = fopen(fname,"w");
	if(fp == NULL) SETERRQ1(PETSC_COMM_SELF, 1,"cannot open file %s", fname);
	free(fname);

	// link output buffer to file
	OutBufConnectToFile(outbuf, fp);

	// write header
	WriteXMLHeader(fp, "RectilinearGrid");

	// open rectilinear grid data block (write total grid size)
	fprintf(fp, "\t<RectilinearGrid WholeExtent=\"%lld %lld %lld %lld %lld %lld\">\n",
		(LLD)(fs->dsx.starts[rx] + 1), (LLD)(fs->dsx.starts[rx+1] + 1),
		(LLD)(fs->dsy.starts[ry] + 1), (LLD)(fs->dsy.starts[ry+1] + 1),
		(LLD)(fs->dsz.starts[rz] + 1), (LLD)(fs->dsz.starts[rz+1] + 1));

	// open sub-domain (piece) description block
	fprintf(fp, "\t\t<Piece Extent=\"%lld %lld %lld %lld %lld %lld\">\n",
		(LLD)(fs->dsx.starts[rx] + 1), (LLD)(fs->dsx.starts[rx+1] + 1),
		(LLD)(fs->dsy.starts[ry] + 1), (LLD)(fs->dsy.starts[ry+1] + 1),
		(LLD)(fs->dsz.starts[rz] + 1), (LLD)(fs->dsz.starts[rz+1] + 1));

	// write cell data block (empty)
	fprintf(fp, "\t\t\t<CellData>\n");
	fprintf(fp, "\t\t\t</CellData>\n");

	// write coordinate block
	fprintf(fp, "\t\t\t<Coordinates>\n");

	fprintf(fp, "\t\t\t\t<DataArray type=\"Float32\" Name=\"Coordinates_X\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)offset);
	offset += sizeof(int) + sizeof(float)*(size_t)nx;

	fprintf(fp, "\t\t\t\t<DataArray type=\"Float32\" Name=\"Coordinates_Y\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)offset);
	offset += sizeof(int) + sizeof(float)*(size_t)ny;

	fprintf(fp, "\t\t\t\t<DataArray type=\"Float32\" Name=\"Coordinates_Z\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)offset);
	offset += sizeof(int) + sizeof(float)*(size_t)nz;

	fprintf(fp, "\t\t\t</Coordinates>\n");

	// write description of output vectors (parameterized)
	outvecs = pvout->outvecs;
	fprintf(fp, "\t\t\t<PointData>\n");
	for(i = 0; i < pvout->nvec; i++)
	{	fprintf(fp, "\t\t\t\t<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"%lld\" format=\"appended\" offset=\"%lld\"/>\n",
			outvecs[i].name, (LLD)outvecs[i].ncomp, (LLD)offset);
		// update offset
		offset += sizeof(int) + sizeof(float)*(size_t)(nx*ny*nz*outvecs[i].ncomp);
	}
	fprintf(fp, "\t\t\t</PointData>\n");

	// close sub-domain and grid blocks
	fprintf(fp, "\t\t</Piece>\n");
	fprintf(fp, "\t</RectilinearGrid>\n");

	// write appended data section
	fprintf(fp, "\t<AppendedData encoding=\"raw\">\n");
	fprintf(fp,"_");

	// coordinate vectors
	OutBufPutCoordVec(outbuf, &fs->dsx, jr->scal->length); OutBufDump(outbuf);
	OutBufPutCoordVec(outbuf, &fs->dsy, jr->scal->length); OutBufDump(outbuf);
	OutBufPutCoordVec(outbuf, &fs->dsz, jr->scal->length); OutBufDump(outbuf);

	for(i = 0; i < pvout->nvec; i++)
	{
		// compute each output vector using its own setup function
		ierr = outvecs[i].OutVecWrite(&outvecs[i]); CHKERRQ(ierr);
		// write vector to output file
		OutBufDump(outbuf);
	}

	// close appended data section and file
	fprintf(fp, "\n\t</AppendedData>\n");
	fprintf(fp, "</VTKFile>\n");

	// close file
	fclose(fp);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//........................... Service Functions .............................
//---------------------------------------------------------------------------
void WriteXMLHeader(FILE *fp, const char *file_type)
{
	// write standard header to ParaView XML file
	fprintf(fp,"<?xml version=\"1.0\"?>\n");
#ifdef PETSC_WORDS_BIGENDIAN
	fprintf(fp,"<VTKFile type=\"%s\" version=\"0.1\" byte_order=\"BigEndian\">\n", file_type);
#else
	fprintf(fp,"<VTKFile type=\"%s\" version=\"0.1\" byte_order=\"LittleEndian\">\n", file_type);
#endif
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "UpdatePVDFile"
PetscErrorCode UpdatePVDFile(
		const char *dirName, const char *outfile, const char *ext,
		long int *offset, PetscScalar ttime, PetscInt outpvd)
{
	FILE        *fp;
	char        *fname;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// check whether pvd is requested
	if(!outpvd) PetscFunctionReturn(0);

	// only first process generates this file (WARNING! Bottleneck!)
	if(!ISRankZero(PETSC_COMM_WORLD)) PetscFunctionReturn(0);

	// open outfile.pvd file (write or update mode)
	asprintf(&fname, "%s.pvd", outfile);
	if(!ttime) fp = fopen(fname,"w");
	else       fp = fopen(fname,"r+");
	free(fname);

	if(fp == NULL) SETERRQ1(PETSC_COMM_SELF, 1,"cannot open file %s", fname);

	if(!ttime)
	{
		// write header
		WriteXMLHeader(fp, "Collection");

		// open time step collection
		fprintf(fp,"<Collection>\n");
	}
	else
	{
		// put the file pointer on the next entry
		ierr = fseek(fp, (*offset), SEEK_SET); CHKERRQ(ierr);
	}

	// add entry to .pvd file
	fprintf(fp,"\t<DataSet timestep=\"%1.6e\" file=\"%s/%s.%s\"/>\n",
		ttime, dirName, outfile, ext);

	// store current position in the file
	(*offset) = ftell(fp);

	// close time step collection
	fprintf(fp,"</Collection>\n");
	fprintf(fp,"</VTKFile>\n");

	// close file
	fclose(fp);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

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
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "paraViewOutBin.h"
#include "outFunct.h"
#include "tools.h"
#include "nlsolveExplicit.h"
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

	// clear object
	ierr = PetscMemzero(outbuf, sizeof(OutBuf)); CHKERRQ(ierr);

	// initialize parameters
	outbuf->fs    = fs;
	outbuf->fp    = NULL;
	outbuf->cn    = 0;

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
//...........  Multi-component output vector data structure .................
//---------------------------------------------------------------------------
void OutVecCreate(
	OutVec         *outvec,
	const char     *name,
	const char     *label,
	OutVecFunctPtr  OutVecFunct,
	PetscInt        ncomp)
{
	// store name
	asprintf(&outvec->name, "%s %s", name, label);

	// output function
	outvec->OutVecFunct = OutVecFunct;

	// number of components
	outvec->ncomp = ncomp;

}
//---------------------------------------------------------------------------
void OutVecDestroy(OutVec *outvec)
{
	LAMEM_FREE(outvec->name);
}
//---------------------------------------------------------------------------
//.......................... Vector output mask .............................
//---------------------------------------------------------------------------
void OutMaskSetDefault(OutMask *omask)
{
	// clear
	PetscMemzero(omask, sizeof(OutMask));

	omask->phase          		= 1;
	omask->visc_total     		= 1;
	omask->visc_creep     		= 1;
	omask->visc_viscoplastic 	= 1;
	
	omask->velocity       		= 1;
	omask->pressure       		= 1;
}
//---------------------------------------------------------------------------
PetscInt OutMaskCountActive(OutMask *omask)
{
	PetscInt cnt = 0;

	if(omask->phase)          		cnt++; // phase
	if(omask->density)        		cnt++; // density
	if(omask->visc_total)     		cnt++; // total effective viscosity
	if(omask->visc_creep)     		cnt++; // creep effective viscosity
	if(omask->visc_viscoplastic) 	cnt++; // viscoplastic viscosity
	if(omask->velocity)       		cnt++; // velocity
	if(omask->pressure)       		cnt++; // pressure
	if(omask->overpressure)   		cnt++; // overpressure
	if(omask->lithospressure) 	    cnt++; // lithostatic pressure
	if(omask->temperature)    		cnt++; // temperature
	if(omask->dev_stress)     		cnt++; // deviatoric stress tensor
	if(omask->j2_dev_stress)  		cnt++; // deviatoric stress second invariant
	if(omask->strain_rate)    		cnt++; // deviatoric strain rate tensor
	if(omask->j2_strain_rate) 		cnt++; // deviatoric strain rate second invariant
	if(omask->vol_rate)       		cnt++; // volumetric strain rate
	if(omask->vorticity)      		cnt++; // vorticity vector
	if(omask->ang_vel_mag)    		cnt++; // average angular velocity magnitude
	if(omask->tot_strain)     		cnt++; // total strain
	if(omask->plast_strain)   		cnt++; // accumulated plastic strain
	if(omask->plast_dissip)   		cnt++; // plastic dissipation
	if(omask->tot_displ)     		cnt++; // total displacements
	if(omask->SHmax)          		cnt++; // maximum horizontal stress
	if(omask->EHmax)          		cnt++; // maximum horizontal stress
	if(omask->ISA)            		cnt++; // Infinite Strain Axis
	if(omask->GOL)            		cnt++; // Grain Orientation Lag
	// === debugging vectors ===============================================
	if(omask->moment_res)     		cnt++; // momentum residual
	if(omask->cont_res)      	 	cnt++; // continuity residual
	if(omask->energ_res)      		cnt++; // energy residual
	if(omask->jac_test)       		cnt++; // matrix-vector Jacobian test

	return cnt;
}
//---------------------------------------------------------------------------
//...................... ParaView output driver object ......................
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutClear"
PetscErrorCode PVOutClear(PVOut *pvout)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// clear object
	ierr = PetscMemzero(pvout, sizeof(PVOut)); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutCreate"
PetscErrorCode PVOutCreate(PVOut *pvout, JacRes *jr, const char *filename)
{
	Scaling  *scal;
	OutMask  *omask;
	OutVec   *outvecs;
	PetscInt  cnt;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	scal = &jr->scal;

	// set file name
	asprintf(&pvout->outfile, "%s", filename);

	// get output mask
	omask = &pvout->omask;

	// activate default vectors
	OutMaskSetDefault(omask);

	// create output buffer object
	ierr = OutBufCreate(&pvout->outbuf, jr); CHKERRQ(ierr);

	// set pvd file flag & offset
	pvout->offset = 0;
	pvout->outpvd = 0;

	// read options
	ierr = PVOutReadFromOptions(pvout); CHKERRQ(ierr);

	//===============
	// OUTPUT VECTORS
	//===============

	// count active output vectors
	pvout->nvec = OutMaskCountActive(omask);

	if(jr->actTemp != PETSC_TRUE) omask->energ_res = 0;

	// allocate space
	ierr = PetscMalloc(sizeof(OutVec)*(size_t)pvout->nvec, &pvout->outvecs); CHKERRQ(ierr);

	// access output vectors
	outvecs = pvout->outvecs;

	// set all output functions & collect information to allocate buffers
	cnt = 0;

	if(omask->phase)          		OutVecCreate(&outvecs[cnt++], "phase",         		scal->lbl_unit,             &PVOutWritePhase,        1);
	if(omask->density)        		OutVecCreate(&outvecs[cnt++], "density",       		scal->lbl_density,          &PVOutWriteDensity,      1);
	if(omask->visc_total)     		OutVecCreate(&outvecs[cnt++], "visc_total",    		scal->lbl_viscosity,        &PVOutWriteViscTotal,    1);
	if(omask->visc_creep)     		OutVecCreate(&outvecs[cnt++], "visc_creep",    		scal->lbl_viscosity,        &PVOutWriteViscCreep,    1);
	if(omask->visc_viscoplastic) 	OutVecCreate(&outvecs[cnt++], "visc_viscoplastic",	scal->lbl_viscosity,        &PVOutWriteViscoPlastic, 1);
	if(omask->velocity)       		OutVecCreate(&outvecs[cnt++], "velocity",       	scal->lbl_velocity,         &PVOutWriteVelocity,     3);
	if(omask->pressure)       		OutVecCreate(&outvecs[cnt++], "pressure",       	scal->lbl_stress,           &PVOutWritePressure,     1);
	if(omask->overpressure)   		OutVecCreate(&outvecs[cnt++], "overpressure",   	scal->lbl_stress,           &PVOutWriteOverPressure, 1);
	if(omask->lithospressure)       OutVecCreate(&outvecs[cnt++], "lithospressure",     scal->lbl_stress,           &PVOutWriteLithosPressure,1);
	if(omask->temperature)    		OutVecCreate(&outvecs[cnt++], "temperature",    	scal->lbl_temperature,      &PVOutWriteTemperature,  1);
	if(omask->dev_stress)     		OutVecCreate(&outvecs[cnt++], "dev_stress",     	scal->lbl_stress,           &PVOutWriteDevStress,    9);
	if(omask->strain_rate)    		OutVecCreate(&outvecs[cnt++], "strain_rate",    	scal->lbl_strain_rate,      &PVOutWriteStrainRate,   9);
	if(omask->j2_dev_stress) 		OutVecCreate(&outvecs[cnt++], "j2_dev_stress",  	scal->lbl_stress,           &PVOutWriteJ2DevStress,  1);
	if(omask->j2_strain_rate) 		OutVecCreate(&outvecs[cnt++], "j2_strain_rate", 	scal->lbl_strain_rate,      &PVOutWriteJ2StrainRate, 1);
	if(omask->vol_rate)       		OutVecCreate(&outvecs[cnt++], "vol_rate",       	scal->lbl_strain_rate,      &PVOutWriteVolRate,      1);
	if(omask->vorticity)      		OutVecCreate(&outvecs[cnt++], "vorticity",      	scal->lbl_strain_rate,      &PVOutWriteVorticity,    3);
	if(omask->ang_vel_mag)    		OutVecCreate(&outvecs[cnt++], "ang_vel_mag",    	scal->lbl_angular_velocity, &PVOutWriteAngVelMag,    1);
	if(omask->tot_strain)     		OutVecCreate(&outvecs[cnt++], "tot_strain",     	scal->lbl_unit,             &PVOutWriteTotStrain,    1);
	if(omask->plast_strain)   		OutVecCreate(&outvecs[cnt++], "plast_strain",   	scal->lbl_unit,             &PVOutWritePlastStrain,  1);
	if(omask->plast_dissip)   		OutVecCreate(&outvecs[cnt++], "plast_dissip",   	scal->lbl_dissipation_rate, &PVOutWritePlastDissip,  1);
	if(omask->tot_displ)      		OutVecCreate(&outvecs[cnt++], "tot_displ",      	scal->lbl_length,           &PVOutWriteTotDispl,     3);
	if(omask->SHmax)          		OutVecCreate(&outvecs[cnt++], "SHmax",          	scal->lbl_unit,             &PVOutWriteSHmax,        3);
	if(omask->EHmax)          		OutVecCreate(&outvecs[cnt++], "EHmax",          	scal->lbl_unit,             &PVOutWriteEHmax,        3);
	if(omask->ISA)            		OutVecCreate(&outvecs[cnt++], "ISA",            	scal->lbl_unit,             &PVOutWriteISA,          3);
	if(omask->GOL)            		OutVecCreate(&outvecs[cnt++], "GOL",            	scal->lbl_unit,             &PVOutWriteGOL,          1);
	// === debugging vectors ===============================================
	if(omask->moment_res)     		OutVecCreate(&outvecs[cnt++], "moment_res",     	scal->lbl_volumetric_force, &PVOutWriteMomentRes,    3);
	if(omask->cont_res)       		OutVecCreate(&outvecs[cnt++], "cont_res",       	scal->lbl_strain_rate,      &PVOutWriteContRes,      1);
	if(omask->energ_res)      		OutVecCreate(&outvecs[cnt++], "energ_res",      	scal->lbl_dissipation_rate, &PVOutWritEnergRes,      1);
	if(omask->jac_test)       		OutVecCreate(&outvecs[cnt++], "jac_test",       	scal->lbl_unit,             &PVOutWriteJacTest,      3);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutReadFromOptions"
PetscErrorCode PVOutReadFromOptions(PVOut *pvout)
{
	OutMask  *omask;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get output mask
	omask = &pvout->omask;

	ierr = PetscOptionsGetInt(NULL, NULL, "-out_pvd",            		&pvout->outpvd,         	NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, NULL, "-out_phase",          		&omask->phase,          	NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, NULL, "-out_density",        		&omask->density,        	NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, NULL, "-out_visc_total",     		&omask->visc_total,     	NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, NULL, "-out_visc_creep",     		&omask->visc_creep,     	NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, NULL, "-out_visc_viscoplastic",     &omask->visc_viscoplastic,  NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, NULL, "-out_velocity",       		&omask->velocity,      	 	NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, NULL, "-out_pressure",       		&omask->pressure,       	NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, NULL, "-out_overpressure",   		&omask->overpressure,   	NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, NULL, "-out_lithospressure",        &omask->lithospressure,     NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, NULL, "-out_temperature",    		&omask->temperature,    	NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, NULL, "-out_dev_stress",     		&omask->dev_stress,     	NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, NULL, "-out_j2_dev_stress",  		&omask->j2_dev_stress,  	NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, NULL, "-out_strain_rate",    		&omask->strain_rate,    	NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, NULL, "-out_j2_strain_rate", 		&omask->j2_strain_rate, 	NULL); CHKERRQ(ierr);
//	ierr = PetscOptionsGetInt(NULL, NULL, "-out_vol_rate",       		&omask->vol_rate,       	NULL); CHKERRQ(ierr);
//	ierr = PetscOptionsGetInt(NULL, NULL, "-out_vorticity",      		&omask->vorticity,      	NULL); CHKERRQ(ierr);
//	ierr = PetscOptionsGetInt(NULL, NULL, "-out_ang_vel_mag",    		&omask->ang_vel_mag,    	NULL); CHKERRQ(ierr);
//	ierr = PetscOptionsGetInt(NULL, NULL, "-out_tot_strain",     		&omask->tot_strain,     	NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, NULL, "-out_shmax",          		&omask->SHmax,          	NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, NULL, "-out_ehmax",          		&omask->EHmax,          	NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, NULL, "-out_isa",            		&omask->ISA,            	NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, NULL, "-out_gol",            		&omask->GOL,            	NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, NULL, "-out_plast_strain",   		&omask->plast_strain,   	NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, NULL, "-out_plast_dissip",   		&omask->plast_dissip,   	NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, NULL, "-out_tot_displ",      		&omask->tot_displ,      	NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, NULL, "-out_moment_res",     		&omask->moment_res,     	NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, NULL, "-out_cont_res",       		&omask->cont_res,       	NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-out_energ_res",      		&omask->energ_res,      	NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-out_jac_test",       		&omask->jac_test,       	NULL); CHKERRQ(ierr);

	if(pvout->outpvd)
	{
		PetscPrintf(PETSC_COMM_WORLD, " Writing grid .pvd file to disk\n");
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutDestroy"
PetscErrorCode PVOutDestroy(PVOut *pvout)
{
	PetscInt i;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// file name
	LAMEM_FREE(pvout->outfile);

	// output vectors
	for(i = 0; i < pvout->nvec; i++)
		OutVecDestroy(&pvout->outvecs[i]);

	PetscFree(pvout->outvecs);

	// output buffer
	ierr = OutBufDestroy(&pvout->outbuf); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWriteTimeStep"
PetscErrorCode PVOutWriteTimeStep(PVOut *pvout, JacRes *jr, const char *dirName, PetscScalar ttime, PetscInt tindx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// update .pvd file if necessary
	if(pvout->outpvd)
	{
		ierr = UpdatePVDFile(dirName, pvout->outfile, "pvtr", &pvout->offset, ttime, tindx); CHKERRQ(ierr);
	}

	// write parallel data .pvtr file
	ierr = PVOutWritePVTR(pvout, dirName); CHKERRQ(ierr);

//ierr = ShowValues(jr,7); CHKERRQ(ierr);

	// write sub-domain data .vtr files
	ierr = PVOutWriteVTR(pvout, jr, dirName); CHKERRQ(ierr);

//ierr = ShowValues(jr,8); CHKERRQ(ierr);

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
		fprintf(fp, "\t\t<Piece Extent=\"%lld %lld %lld %lld %lld %lld\" Source=\"%s_p%1.6lld.vtr\"/>\n",
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
PetscErrorCode PVOutWriteVTR(PVOut *pvout, JacRes *jr, const char *dirName)
{
	FILE          *fp;
	FDSTAG        *fs;
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
	fs     = outbuf->fs;

	// get sizes of output grid
	GET_OUTPUT_RANGE(rx, nx, sx, fs->dsx)
	GET_OUTPUT_RANGE(ry, ny, sy, fs->dsy)
	GET_OUTPUT_RANGE(rz, nz, sz, fs->dsz)

	// open outfile_p_XXXXXX.vtr file in the output directory (write mode)
	asprintf(&fname, "%s/%s_p%1.6lld.vtr", dirName, pvout->outfile, (LLD)rank);
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
	OutBufPutCoordVec(outbuf, &fs->dsx, jr->scal.length); OutBufDump(outbuf);
	OutBufPutCoordVec(outbuf, &fs->dsy, jr->scal.length); OutBufDump(outbuf);
	OutBufPutCoordVec(outbuf, &fs->dsz, jr->scal.length); OutBufDump(outbuf);

//PetscPrintf(PETSC_COMM_WORLD, "    7aa ------------------------------------------\n");
//ierr = ShowValues(jr); CHKERRQ(ierr);

	for(i = 0; i < pvout->nvec; i++)
	{
		// compute each output vector using its own setup function
		ierr = outvecs[i].OutVecFunct(jr, outbuf); CHKERRQ(ierr);

		// write vector to output file
		OutBufDump(outbuf);
	}

//PetscPrintf(PETSC_COMM_WORLD, "    7ab ------------------------------------------\n");
//ierr = ShowValues(jr); CHKERRQ(ierr);

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
		long int *offset, PetscScalar ttime, PetscInt tindx)
{
	FILE        *fp;
	char        *fname;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// only first process generates this file (WARNING! Bottleneck!)
	if(!ISRankZero(PETSC_COMM_WORLD)) PetscFunctionReturn(0);

	// open outfile.pvd file (write or update mode)
	asprintf(&fname, "%s.pvd", outfile);
	if(!tindx) fp = fopen(fname,"w");
	else       fp = fopen(fname,"r+");
	if(fp == NULL) SETERRQ1(PETSC_COMM_SELF, 1,"cannot open file %s", fname);
	free(fname);

	if(!tindx)
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

//---------------------------------------------------------------------------
//.................   FDSTAG PARAVIEW XML OUTPUT ROUTINES   .................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "fdstag.h"
#include "solVar.h"
#include "scaling.h"
#include "bc.h"
#include "JacRes.h"
#include "paraViewOutBin.h"
#include "outFunct.h"
#include "Utils.h"

//---------------------------------------------------------------------------
//............................. Output buffer ...............................
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "OutBufCreate"
PetscErrorCode OutBufCreate(
	OutBuf   *outbuf,
	FDSTAG   *fs,
	PetscInt  maxnc,
	PetscBool mkcen,
	PetscBool mkedg)
{
	PetscInt rx, ry, rz, sx, sy, sz, nx, ny, nz;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// initialize parameters
	outbuf->fs    = fs;
	outbuf->fp    = NULL;
	outbuf->maxnc = maxnc;
	outbuf->mkcen = mkcen;
	outbuf->mkedg = mkedg;
	outbuf->cn    = 0;

	// get local output grid sizes
	GET_OUTPUT_RANGE(rx, nx, sx, fs->dsx)
	GET_OUTPUT_RANGE(ry, ny, sy, fs->dsy)
	GET_OUTPUT_RANGE(rz, nz, sz, fs->dsz)

	// allocate output buffer
	ierr = PetscMalloc((size_t)(nx*ny*nz*maxnc)*sizeof(float), &outbuf->buff); CHKERRQ(ierr);

	// allocate corner buffers
	ierr = DMCreateGlobalVector(fs->DA_COR, &outbuf->gbcor); CHKERRQ(ierr);
	ierr = DMCreateLocalVector (fs->DA_COR, &outbuf->lbcor); CHKERRQ(ierr);

	// allocate center buffers
	if(mkcen)
	{
		ierr = DMCreateGlobalVector(fs->DA_CEN, &outbuf->gbcen); CHKERRQ(ierr);
		ierr = DMCreateLocalVector (fs->DA_CEN, &outbuf->lbcen); CHKERRQ(ierr);
	}

	// allocate edge buffers
	if(mkedg)
	{
		ierr = DMCreateGlobalVector(fs->DA_XY,  &outbuf->gbxy); CHKERRQ(ierr);
		ierr = DMCreateLocalVector (fs->DA_XY,  &outbuf->lbxy); CHKERRQ(ierr);

		ierr = DMCreateGlobalVector(fs->DA_XZ,  &outbuf->gbxz); CHKERRQ(ierr);
		ierr = DMCreateLocalVector (fs->DA_XZ,  &outbuf->lbxz); CHKERRQ(ierr);

		ierr = DMCreateGlobalVector(fs->DA_YZ,  &outbuf->gbyz); CHKERRQ(ierr);
		ierr = DMCreateLocalVector (fs->DA_YZ,  &outbuf->lbyz); CHKERRQ(ierr);
	}

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

	// free corner buffers
	ierr = VecDestroy(&outbuf->gbcor); CHKERRQ(ierr);
	ierr = VecDestroy(&outbuf->lbcor); CHKERRQ(ierr);

	// free center buffers
	if(outbuf->mkcen)
	{
		ierr = VecDestroy(&outbuf->gbcen); CHKERRQ(ierr);
		ierr = VecDestroy(&outbuf->lbcen); CHKERRQ(ierr);
	}

	// free edge buffers
	if(outbuf->mkedg)
	{
		ierr = VecDestroy(&outbuf->gbxy); CHKERRQ(ierr);
		ierr = VecDestroy(&outbuf->lbxy); CHKERRQ(ierr);

		ierr = VecDestroy(&outbuf->gbxz); CHKERRQ(ierr);
		ierr = VecDestroy(&outbuf->lbxz); CHKERRQ(ierr);

		ierr = VecDestroy(&outbuf->gbyz); CHKERRQ(ierr);
		ierr = VecDestroy(&outbuf->lbyz); CHKERRQ(ierr);
	}

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
	nbytes = outbuf->cn*(int)sizeof(float);

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
	PetscInt     ncomp, // number of components
	PetscInt     dir,   // component identifier
	PetscScalar  cf)    // scaling coefficient
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
	GLOBAL_TO_LOCAL(fs->DA_COR, outbuf->gbcor, outbuf->lbcor)

	// access local buffer vector
	ierr = DMDAVecGetArray(fs->DA_COR, outbuf->lbcor, &arr); CHKERRQ(ierr);

	// get sub-domain ranks, starting node IDs, and number of nodes
	GET_OUTPUT_RANGE(rx, nx, sx, fs->dsx)
	GET_OUTPUT_RANGE(ry, ny, sy, fs->dsy)
	GET_OUTPUT_RANGE(rz, nz, sz, fs->dsz)

	// set counter
	cnt = dir;
	// copy vector component to buffer
	START_STD_LOOP
	{	// write
		buff[cnt] = (float) (cf*arr[k][j][i]);
		// update counter
		cnt += ncomp;
	}
	END_STD_LOOP
	// restore access
	ierr = DMDAVecRestoreArray(fs->DA_COR, outbuf->lbcor, &arr); CHKERRQ(ierr);

	// update number of elements in the buffer
	outbuf->cn += nx*ny*nz;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
//...........  Multi-component output vector data structure .................
//---------------------------------------------------------------------------
void OutVecCreate(
		OutVec         *outvec,
		const char     *name,
		OutVecFunctPtr  OutVecFunct,
		PetscInt        ncomp,
		PetscInt        cen,
		PetscInt        edg,
		PetscInt       *maxnc,
		PetscInt       *mkcen,
		PetscInt       *mkedg)
{
	// store name
	asprintf(&outvec->name, "%s", name);

	// output function
	outvec->OutVecFunct = OutVecFunct;

	// number of components
	outvec->ncomp = ncomp;

	// track maximum number of components
	if(*maxnc < ncomp) *maxnc = ncomp;

	// track center buffer request
	if(cen) *mkcen = 1;

	// track edge buffer request
	if(edg) *mkedg = 1;
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
	memset(omask, 0, sizeof(OutMask));

// ADHOC
/*
	// activate default output vectors
	omask->phase          = 1;
	omask->density        = 1;
	omask->viscosity      = 1;
	omask->velocity       = 1;
	omask->pressure       = 1;
	omask->dev_stress     = 1;
	omask->j2_dev_stress  = 1;
	omask->strain_rate    = 1;
	omask->j2_strain_rate = 1;
	omask->moment_res     = 1;
	omask->cont_res       = 1;
*/
	omask->phase          = 1;
	omask->viscosity      = 1;
	omask->velocity       = 1;
	omask->pressure       = 1;
}
//---------------------------------------------------------------------------
PetscInt OutMaskCountActive(OutMask *omask)
{
	PetscInt cnt = 0;
	if(omask->phase)          cnt++; // phase
	if(omask->density)        cnt++; // density
	if(omask->viscosity)      cnt++; // effective viscosity
	if(omask->velocity)       cnt++; // velocity
	if(omask->pressure)       cnt++; // pressure
	if(omask->temperature)    cnt++; // temperature
	if(omask->moment_res)     cnt++; // momentum residual   (debugging)
	if(omask->cont_res)       cnt++; // continuity residual (debugging)
	if(omask->energ_res)      cnt++; // energy residual     (debugging)
	if(omask->dev_stress)     cnt++; // deviatoric stress tensor
	if(omask->j2_dev_stress)  cnt++; // deviatoric stress second invariant
	if(omask->strain_rate)    cnt++; // deviatoric strain rate tensor
	if(omask->j2_strain_rate) cnt++; // deviatoric strain rate second invariant
	if(omask->vol_rate)       cnt++; // volumetric strain rate
	if(omask->vorticity)      cnt++; // vorticity vector
	if(omask->ang_vel_mag)    cnt++; // average angular velocity magnitude
	if(omask->tot_strain)     cnt++; // total strain
	if(omask->plast_strain)   cnt++; // accumulated plastic strain
	if(omask->plast_dissip)   cnt++; // plastic dissipation
	if(omask->tot_displ)      cnt++; // total displacements
	if(omask->DII_CEN)        cnt++; // effective strain rate invariant in center  (debugging)
	if(omask->DII_XY)         cnt++; // effective strain rate invariant on xy-edge (debugging)
	if(omask->DII_XZ)         cnt++; // effective strain rate invariant on xz-edge (debugging)
	if(omask->DII_YZ)         cnt++; // effective strain rate invariant on yz-edge (debugging)

	return cnt;
}
//---------------------------------------------------------------------------
//...................... ParaView output driver object ......................
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutCreate"
PetscErrorCode PVOutCreate(PVOut *pvout, FDSTAG *fs, Scaling *scal, const char *filename)
{
	PetscInt  cnt, maxnc, mkcen, mkedg;
	OutMask  *omask;
	OutVec   *outvecs;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// set file name
	asprintf(&pvout->outfile, "%s", filename);

	// set output scaling for coordinates
	pvout->crdScal = scal->out_length;

	// get output mask
	omask = &pvout->omask;

	// activate default vectors
	OutMaskSetDefault(omask);

	// read output mask from file here
	// ... NOT IMPLEMENTED YET

// ADHOC (uncomment output flags)

//	memset(omask, 0, sizeof(OutMask));

//	omask->DII_CEN        = 1;
//	omask->DII_XY         = 1;
//	omask->DII_XZ         = 1;
//	omask->DII_YZ         = 1;
//	omask->velocity       = 1;
//	omask->j2_strain_rate = 1;
//	omask->strain_rate    = 1;

	// count active output vectors
	pvout->nvec = OutMaskCountActive(omask);

	// allocate space
	ierr = PetscMalloc(sizeof(OutVec)*(size_t)pvout->nvec, &pvout->outvecs); CHKERRQ(ierr);

	// access output vectors
	outvecs = pvout->outvecs;

	// set all output functions & collect information to allocate buffers
	cnt   = 0;
	maxnc = 0;
	mkcen = 0;
	mkedg = 0;

	if(omask->phase)          OutVecCreate(&outvecs[cnt++], "phase",          &PVOutWritePhase,        1, 1, 0, &maxnc, &mkcen, &mkedg);
	if(omask->density)        OutVecCreate(&outvecs[cnt++], "density",        &PVOutWriteDensity,      1, 1, 0, &maxnc, &mkcen, &mkedg);
	if(omask->viscosity)      OutVecCreate(&outvecs[cnt++], "viscosity",      &PVOutWriteViscosity,    1, 1, 0, &maxnc, &mkcen, &mkedg);
	if(omask->velocity)       OutVecCreate(&outvecs[cnt++], "velocity",       &PVOutWriteVelocity,     3, 0, 0, &maxnc, &mkcen, &mkedg);
	if(omask->pressure)       OutVecCreate(&outvecs[cnt++], "pressure",       &PVOutWritePressure,     1, 0, 0, &maxnc, &mkcen, &mkedg);
	if(omask->temperature)    OutVecCreate(&outvecs[cnt++], "temperature",    &PVOutWriteTemperature,  1, 0, 0, &maxnc, &mkcen, &mkedg);
	if(omask->moment_res)     OutVecCreate(&outvecs[cnt++], "moment_res",     &PVOutWriteMomentRes,    3, 0, 0, &maxnc, &mkcen, &mkedg);
	if(omask->cont_res)       OutVecCreate(&outvecs[cnt++], "cont_res",       &PVOutWriteContRes,      1, 0, 0, &maxnc, &mkcen, &mkedg);
	if(omask->energ_res)      OutVecCreate(&outvecs[cnt++], "energ_res",      &PVOutWritEnergRes,      1, 0, 0, &maxnc, &mkcen, &mkedg);
	if(omask->dev_stress)     OutVecCreate(&outvecs[cnt++], "dev_stress",     &PVOutWriteDevStress,    6, 1, 1, &maxnc, &mkcen, &mkedg);
	if(omask->j2_dev_stress)  OutVecCreate(&outvecs[cnt++], "j2_dev_stress",  &PVOutWriteJ2DevStress,  1, 1, 0, &maxnc, &mkcen, &mkedg);
	if(omask->strain_rate)    OutVecCreate(&outvecs[cnt++], "strain_rate",    &PVOutWriteStrainRate,   6, 1, 1, &maxnc, &mkcen, &mkedg);
	if(omask->j2_strain_rate) OutVecCreate(&outvecs[cnt++], "j2_strain_rate", &PVOutWriteJ2StrainRate, 1, 1, 1, &maxnc, &mkcen, &mkedg);
	if(omask->vol_rate)       OutVecCreate(&outvecs[cnt++], "vol_rate",       &PVOutWriteVolRate,      1, 1, 0, &maxnc, &mkcen, &mkedg);
	if(omask->vorticity)      OutVecCreate(&outvecs[cnt++], "vorticity",      &PVOutWriteVorticity,    3, 0, 1, &maxnc, &mkcen, &mkedg);
	if(omask->ang_vel_mag)    OutVecCreate(&outvecs[cnt++], "ang_vel_mag",    &PVOutWriteAngVelMag,    1, 0, 1, &maxnc, &mkcen, &mkedg);
	if(omask->tot_strain)     OutVecCreate(&outvecs[cnt++], "tot_strain",     &PVOutWriteTotStrain,    1, 0, 0, &maxnc, &mkcen, &mkedg);
	if(omask->plast_strain)   OutVecCreate(&outvecs[cnt++], "plast_strain",   &PVOutWritePlastStrain,  1, 0, 0, &maxnc, &mkcen, &mkedg);
	if(omask->plast_dissip)   OutVecCreate(&outvecs[cnt++], "plast_dissip",   &PVOutWritePlastDissip,  1, 0, 0, &maxnc, &mkcen, &mkedg);
	if(omask->tot_displ)      OutVecCreate(&outvecs[cnt++], "tot_displ",      &PVOutWriteTotDispl,     3, 0, 0, &maxnc, &mkcen, &mkedg);
	if(omask->DII_CEN)        OutVecCreate(&outvecs[cnt++], "DII_CEN",        &PVOutWriteDII_CEN,      1, 1, 0, &maxnc, &mkcen, &mkedg);
	if(omask->DII_XY)         OutVecCreate(&outvecs[cnt++], "DII_XY",         &PVOutWriteDII_XY,       1, 0, 1, &maxnc, &mkcen, &mkedg);
	if(omask->DII_XZ)         OutVecCreate(&outvecs[cnt++], "DII_XZ",         &PVOutWriteDII_XZ,       1, 0, 1, &maxnc, &mkcen, &mkedg);
	if(omask->DII_YZ)         OutVecCreate(&outvecs[cnt++], "DII_YZ",         &PVOutWriteDII_YZ,       1, 0, 1, &maxnc, &mkcen, &mkedg);

	// create output buffer object
	ierr = OutBufCreate(&pvout->outbuf, fs, maxnc, mkcen, mkedg);

	//================================
	// time step buffer
	// to be replaced by python script
	//================================

	// number of collected time steps
	pvout->nstep = 0;

	// buffer counter
	pvout->bcnt = 0;

	// time step database name
	asprintf(&pvout->dbname, "%s_timestep.db", pvout->outfile);

	// buffer for output time stamps
	ierr = makeScalArray(&pvout->bstamps, NULL, _timestep_buff_size_); CHKERRQ(ierr);

	// buffer for indices of saved time steps
	ierr = makeIntArray(&pvout->bindexs, NULL, _timestep_buff_size_); CHKERRQ(ierr);

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

	// write .pvd output file to disk
	ierr = PVOutWritePVD(pvout);

	// file name
	LAMEM_FREE(pvout->outfile);

	// output vectors
	for(i = 0; i < pvout->nvec; i++)
		OutVecDestroy(&pvout->outvecs[i]);

	PetscFree(pvout->outvecs);

	// output buffer
	ierr = OutBufDestroy(&pvout->outbuf); CHKERRQ(ierr);

	//================================
	// time step buffer
	// to be replaced by python script
	//===============================

	LAMEM_FREE(pvout->dbname);

	ierr = PetscFree(pvout->bstamps); CHKERRQ(ierr);

	ierr = PetscFree(pvout->bindexs); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
void PVOutWriteXMLHeader(FILE *fp, const char *file_type)
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
#define __FUNCT__ "PVOutWriteTimeStep"
PetscErrorCode PVOutWriteTimeStep(PVOut *pvout, JacResCtx *jrctx, PetscScalar ttime, PetscInt tindx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// update time step database (required for .pvd file)
	ierr = PVOutUpdateTimeStepBuffer(pvout, ttime, tindx); CHKERRQ(ierr);

	// write parallel data .pvtr file
	ierr = PVOutWritePVTR(pvout, tindx); CHKERRQ(ierr);

	// write sub-domain data .vtr files
	ierr = PVOutWriteVTR(pvout, jrctx, tindx); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutUpdateTimeStepBuffer"
PetscErrorCode PVOutUpdateTimeStepBuffer(PVOut *pvout, PetscScalar ttime, PetscInt tindx)
{
	PetscMPIInt rank;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// only first process generates this file
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
	if(rank != 0) PetscFunctionReturn(0);

	// dump time step buffer if necessary
	if(pvout->bcnt == _timestep_buff_size_)
	{	ierr = PVOutDumpTimeStepBuffer(pvout); CHKERRQ(ierr); }

	// load time step data to buffer
	pvout->bstamps[pvout->bcnt] = ttime;
	pvout->bindexs[pvout->bcnt] = tindx;

	// increment time step & buffer counter
	pvout->nstep++;
	pvout->bcnt++;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutDumpTimeStepBuffer"
PetscErrorCode PVOutDumpTimeStepBuffer(PVOut *pvout)
{
	FILE *db;

	PetscFunctionBegin;

	// return if buffer is empty
	if(pvout->bcnt == 0) PetscFunctionReturn(0);

	// open time step database file
	if (pvout->nstep <= _timestep_buff_size_) db = fopen(pvout->dbname,"w"); // new file (write mode)
	else                                      db = fopen(pvout->dbname,"a"); // existing file (append mode)
	if(db == NULL) SETERRQ1(PETSC_COMM_WORLD, 1,"cannot open file %s", pvout->dbname);

	// dump contents of the buffer to database
	fwrite(pvout->bstamps, sizeof(PetscScalar), (size_t)pvout->bcnt, db);
	fwrite(pvout->bindexs, sizeof(PetscInt),    (size_t)pvout->bcnt, db);

	// zero out buffer counter
	pvout->bcnt = 0;

	// close file
	fclose(db);
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWritePVD"
PetscErrorCode PVOutWritePVD(PVOut *pvout)
{
	FILE          *fp, *db;
	char          *name;
	PetscScalar    ttime;
	PetscInt       tindx, chsz, j, cnt;
	PetscMPIInt    rank;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// only first process generates this file
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
	if(rank != 0) PetscFunctionReturn(0);

	// only proceed if at least one time step was stored
	if(pvout->nstep == 0) PetscFunctionReturn(0);

	// dump the rest of the time step buffer
	ierr = PVOutDumpTimeStepBuffer(pvout); CHKERRQ(ierr);

	// open outfile.pvd file (write, mode)
	asprintf(&name, "%s.pvd", pvout->outfile);
	fp = fopen(name,"w");
	if(fp == NULL) SETERRQ1(PETSC_COMM_WORLD, 1,"cannot open file %s", name);
	free(name);

	// open time step database file (read mode)
	db = fopen(pvout->dbname,"r");
	if(db == NULL) SETERRQ1(PETSC_COMM_WORLD, 1,"cannot open file %s", pvout->dbname);

	// write header
	PVOutWriteXMLHeader(fp, "Collection");

	// open time step collection
	fprintf(fp,"<Collection>\n");

	// write links to time step output files
	cnt = 0;
	do
	{	// get size of next chunk to read
		if((pvout->nstep - cnt) > _timestep_buff_size_) chsz = _timestep_buff_size_;
		else 	                                        chsz = (size_t) (pvout->nstep - cnt);

		// load data to buffer
		fread(pvout->bstamps, sizeof(PetscScalar), (size_t)chsz, db);
		fread(pvout->bindexs, sizeof(PetscInt),    (size_t)chsz, db);

		// process buffer
		for(j = 0; j < chsz; j++)
		{
			// retrieve next time stamp and index
			ttime = pvout->bstamps[j];
			tindx = pvout->bindexs[j];

			// add record to .pvd file (outfile.pvtr in the Temestep_XXXXXX directory)
			fprintf(fp,"\t<DataSet timestep=\"%1.6e\" file=\"Timestep_%1.6lld/%s.pvtr\"/>\n",
				ttime, (LLD)tindx, pvout->outfile);
		}
		cnt += chsz;
	} while (cnt != pvout->nstep);

	// close time step collection
	fprintf(fp,"</Collection>\n");
	fprintf(fp,"</VTKFile>\n");

	// close files
	fclose(fp);
	fclose(db);

	// delete time step database file
	remove(pvout->dbname);
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVOutWritePVTR"
PetscErrorCode PVOutWritePVTR(PVOut *pvout, PetscInt tindx)
{
	FILE        *fp;
	FDSTAG      *fs;
	char        *fname;
	OutVec      *outvecs;
	PetscInt     i, rx, ry, rz;
	PetscMPIInt  rank, nproc, iproc;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// only first process generates this file (WARNING! Potential Bottleneck!)
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
	if(rank != 0)  PetscFunctionReturn(0);

	// access staggered grid layout
	fs = pvout->outbuf.fs;

	// open outfile.pvtr file in the Temestep_XXXXXX directory (write mode)
	asprintf(&fname, "Timestep_%1.6lld/%s.pvtr", (LLD)tindx, pvout->outfile);
	fp = fopen(fname,"w");
	if(fp == NULL) SETERRQ1(PETSC_COMM_WORLD, 1,"cannot open file %s", fname);
	free(fname);

	// write header
	PVOutWriteXMLHeader(fp, "PRectilinearGrid");

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
PetscErrorCode PVOutWriteVTR(PVOut *pvout, JacResCtx *jrctx, PetscInt tindx)
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

	// open outfile_p_XXXXXX.vtr file in the Temestep_XXXXXX directory (write mode)
	asprintf(&fname, "Timestep_%1.6lld/%s_p%1.6lld.vtr", (LLD)tindx, pvout->outfile, (LLD)rank);
	fp = fopen(fname,"w");
	if(fp == NULL) SETERRQ1(PETSC_COMM_WORLD, 1,"cannot open file %s", fname);
	free(fname);

	// link output buffer to file
	OutBufConnectToFile(outbuf, fp);

	// write header
	PVOutWriteXMLHeader(fp, "RectilinearGrid");

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
	OutBufPutCoordVec(outbuf, &fs->dsx, pvout->crdScal); OutBufDump(outbuf);
	OutBufPutCoordVec(outbuf, &fs->dsy, pvout->crdScal); OutBufDump(outbuf);
	OutBufPutCoordVec(outbuf, &fs->dsz, pvout->crdScal); OutBufDump(outbuf);

	for(i = 0; i < pvout->nvec; i++)
	{
		// compute each output vector using its own setup function
		ierr = outvecs[i].OutVecFunct(jrctx, outbuf); CHKERRQ(ierr);
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

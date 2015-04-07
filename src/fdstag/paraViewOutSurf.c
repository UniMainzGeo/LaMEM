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
#include "interpolate.h"
#include "surf.h"
#include "paraViewOutBin.h"
#include "paraViewOutSurf.h"
#include "outFunct.h"
#include "Utils.h"
/*
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
	PetscInt     ncomp,  // number of components
	PetscInt     dir,    // component identifier
	PetscScalar  cf,     // scaling coefficient
	PetscScalar  shift) // shift parameter (subtracted from scaled values)
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

*/
//---------------------------------------------------------------------------
//................ ParaView free surface output driver object ...............
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVSurfClear"
PetscErrorCode PVSurfClear(PVSurf *pvout)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// clear object
	ierr = PetscMemzero(pvout, sizeof(PVSurf)); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVSurfCreate"
PetscErrorCode PVSurfCreate(PVSurf *pvout, FreeSurf *surf, const char *filename)
{
	FDSTAG   *fs;
	PetscInt  rx, ry, sx, sy, nx, ny;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access staggered grid layout
	fs = surf->jr->fs;

	// get local output grid sizes
	GET_OUTPUT_RANGE(rx, nx, sx, fs->dsx)
	GET_OUTPUT_RANGE(ry, ny, sy, fs->dsy)

	// set context
	pvout->surf = surf;

	// set file name
	ierr = asprintf(&pvout->outfile, "%s_surf", filename); CHKERRQ(ierr);

	// allocate output buffer
	ierr = PetscMalloc((size_t)(_max_num_comp_surf_*nx*ny)*sizeof(float), &pvout->buff); CHKERRQ(ierr);

	// set pvd file offset
	pvout->offset = 0;

	// set output mask
	pvout->outpvd     = 0;
	pvout->velocity   = 0;
	pvout->topography = 1;
	pvout->amplitude  = 0;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVSurfReadFromOptions"
PetscErrorCode PVSurfReadFromOptions(PVSurf *pvout)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// read output flags
	ierr = PetscOptionsGetInt(NULL, "-out_surf_pvd",        &pvout->outpvd,     NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, "-out_surf_velocity",   &pvout->velocity,   NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, "-out_surf_topography", &pvout->topography, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, "-out_surf_amplitude",  &pvout->amplitude,  NULL); CHKERRQ(ierr);

	if(pvout->outpvd)
	{
		PetscPrintf(PETSC_COMM_WORLD, " Writing surface .pvd file to disk\n");
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVSurfDestroy"
PetscErrorCode PVSurfDestroy(PVSurf *pvout)
{
	PetscFunctionBegin;

	LAMEM_FREE(pvout->outfile);
	PetscFree (pvout->buff);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVSurfWriteTimeStep"
PetscErrorCode PVSurfWriteTimeStep(PVSurf *pvout, const char *dirName, PetscScalar ttime, PetscInt tindx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// update .pvd file if necessary
	if(pvout->outpvd)
	{
		ierr = UpdatePVDFile(dirName, pvout->outfile, "pvts", &pvout->offset, ttime, tindx); CHKERRQ(ierr);
	}

	// write parallel data .pvts file
	ierr = PVSurfWritePVTS(pvout, dirName); CHKERRQ(ierr);

	// write sub-domain data .vts files
	ierr = PVSurfWriteVTS(pvout, dirName); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVSurfWritePVTS"
PetscErrorCode PVSurfWritePVTS(PVSurf *pvout, const char *dirName)
{
	FILE        *fp;
	FDSTAG      *fs;
	char        *fname;
	Scaling     *scal;
	PetscInt     rx, ry, rz;
	PetscMPIInt  nproc, iproc;

	PetscFunctionBegin;

	// only first process generates this file (WARNING! Bottleneck!)
	if(!ISRankZero(PETSC_COMM_WORLD)) PetscFunctionReturn(0);

	// access context
	fs   =  pvout->surf->jr->fs;
	scal = &pvout->surf->jr->scal;

	// open outfile.pvts file in the output directory (write mode)
	asprintf(&fname, "%s/%s.pvts", dirName, pvout->outfile);
	fp = fopen(fname,"w");
	if(fp == NULL) SETERRQ1(PETSC_COMM_SELF, 1,"cannot open file %s", fname);
	free(fname);

	// write header
	WriteXMLHeader(fp, "PStructuredGrid");

	// open structured grid data block (write total grid size)
	fprintf(fp, "\t<PStructuredGrid GhostLevel=\"0\" WholeExtent=\"%lld %lld %lld %lld 1 1\">\n",
		1LL, (LLD)fs->dsx.tnods,
		1LL, (LLD)fs->dsy.tnods);

	// write cell data block (empty)
	fprintf(fp, "\t\t<PCellData>\n");
	fprintf(fp, "\t\t</PCellData>\n");

	// write coordinate block
	fprintf(fp, "\t\t<PPoints>\n");

	fprintf(fp,"\t\t\t<PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\"/>\n");

	fprintf(fp, "\t\t</PPoints>\n");

	// write description of output vectors
	fprintf(fp, "\t\t<PPointData>\n");

	if(pvout->velocity)
	{
		fprintf(fp,"\t\t\t<PDataArray type=\"Float32\" Name=\"velocity %s\" NumberOfComponents=\"3\" format=\"appended\"/>\n", scal->lbl_velocity);
	}

	if(pvout->topography)
	{
		fprintf(fp,"\t\t\t<PDataArray type=\"Float32\" Name=\"topography %s\" NumberOfComponents=\"1\" format=\"appended\"/>\n", scal->lbl_length);
	}

	if(pvout->amplitude)
	{
		fprintf(fp,"\t\t\t<PDataArray type=\"Float32\" Name=\"amplitude %s\" NumberOfComponents=\"1\" format=\"appended\"/>\n", scal->lbl_length);
	}

	fprintf(fp, "\t\t</PPointData>\n");

	// get total number of free surface sub-domains
	nproc = fs->dsx.nproc*fs->dsy.nproc;

	// write local grid sizes (extents) and data file names for all sub-domains
	for(iproc = 0; iproc < nproc; iproc++)
	{
		// get sub-domain ranks in all coordinate directions
		getLocalRank(&rx, &ry, &rz, iproc, fs->dsx.nproc, fs->dsy.nproc);

		// write data
		fprintf(fp, "\t\t<Piece Extent=\"%lld %lld %lld %lld 1 1\" Source=\"%s_p%1.6lld.vts\"/>\n",
			(LLD)(fs->dsx.starts[rx] + 1), (LLD)(fs->dsx.starts[rx+1] + 1),
			(LLD)(fs->dsy.starts[ry] + 1), (LLD)(fs->dsy.starts[ry+1] + 1),
			pvout->outfile, (LLD)iproc);
	}

	// close structured grid data block
	fprintf(fp, "\t</PStructuredGrid>\n");
	fprintf(fp, "</VTKFile>\n");

	// close file
	fclose(fp);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVSurfWriteVTS"
PetscErrorCode PVSurfWriteVTS(PVSurf *pvout, const char *dirName)
{

	FILE          *fp;
	FDSTAG        *fs;
	Scaling       *scal;
	char          *fname;
	PetscInt       i, rx, ry, rz, sx, sy, sz, nx, ny, nz;
	size_t         offset = 0;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get staggered grid layout
	fs = pvout->surf->jr->fs;

	// only ranks zero in z direction generate this file
	if(fs->dsz.rank) PetscFunctionReturn(0);

	// open outfile_p_XXXXXX.vts file in the output directory (write mode)
	asprintf(&fname, "%s/%s_p%1.6lld.vts", dirName, pvout->outfile, (LLD)fs->dsz.color);
	fp = fopen(fname,"w");
	if(fp == NULL) SETERRQ1(PETSC_COMM_SELF, 1,"cannot open file %s", fname);
	free(fname);

	// get scaling object
	scal = &pvout->surf->jr->scal;

	// get sizes of output grid
	GET_OUTPUT_RANGE(rx, nx, sx, fs->dsx)
	GET_OUTPUT_RANGE(ry, ny, sy, fs->dsy)

/*
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


	if(omask->velocity)       OutVecCreate(&outvecs[cnt++], "velocity",       scal->lbl_velocity,         &PVOutWriteVelocity,     3);


	// write description of output vectors

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

	for(i = 0; i < pvout->nvec; i++)
	{
		// compute each output vector using its own setup function
		ierr = outvecs[i].OutVecFunct(jr, outbuf); CHKERRQ(ierr);
		// write vector to output file
		OutBufDump(outbuf);
	}
*/
	// close appended data section and file
	fprintf(fp, "\n\t</AppendedData>\n");
	fprintf(fp, "</VTKFile>\n");

	// close file
	fclose(fp);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------





/*

	// Compute average height & subtract it
	VecNorm(FIELD,NORM_1, &AverageHeight);
	VecGetSize(FIELD,&size);
	AverageHeight = AverageHeight/((PetscScalar) size);
	{
		const char *field_name;
		ierr = DMDAGetFieldName( da, 0, &field_name );CHKERRQ(ierr);
		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", "SurfaceAmplitude" );
		for( k=sk+nz-1; k<sk+nz; k++ ) {
			for( j=sj; j<sj+ny; j++ ) {
				for( i=si; i<si+nx; i++ ) {
					double val;
					val = (ZCOORD[k][j][i] - AverageHeight)	* scaling_length;
					fprintf( vtk_fp, "\t\t\t\t\t%1.8f\n", val );
				}
			}
		}

	}
	fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");


*/


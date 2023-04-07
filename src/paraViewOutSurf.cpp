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
 **    filename:   paraViewOutSurf.c
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
//..............   FREE SRUFACE PARAVIEW XML OUTPUT ROUTINES   ..............
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "paraViewOutSurf.h"
#include "paraViewOutBin.h"
#include "parsing.h"
#include "scaling.h"
#include "fdstag.h"
#include "surf.h"
#include "JacRes.h"
#include "tools.h"
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVSurfCreate"
PetscErrorCode PVSurfCreate(PVSurf *pvsurf, FB *fb)
{
	char filename[_str_len_];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// free surface cases only
	if(!pvsurf->surf->UseFreeSurf) PetscFunctionReturn(0);

	// check activation
	ierr = getIntParam(fb, _OPTIONAL_, "out_surf", &pvsurf->outsurf, 1, 1); CHKERRQ(ierr);

	if(!pvsurf->outsurf) PetscFunctionReturn(0);

	// initialize
	pvsurf->outpvd     = 1;
	pvsurf->topography = 1;

	// read
	ierr = getStringParam(fb, _OPTIONAL_, "out_file_name",       filename,        "output"); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_surf_pvd",        &pvsurf->outpvd,     1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_surf_velocity",   &pvsurf->velocity,   1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_surf_topography", &pvsurf->topography, 1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_surf_amplitude",  &pvsurf->amplitude,  1, 1); CHKERRQ(ierr);

	// print summary
	PetscPrintf(PETSC_COMM_WORLD, "Surface output parameters:\n");
	PetscPrintf(PETSC_COMM_WORLD, "   Write .pvd file : %s \n", pvsurf->outpvd ? "yes" : "no");

	if(pvsurf->velocity)   PetscPrintf(PETSC_COMM_WORLD, "   Velocity        @ \n");
	if(pvsurf->topography) PetscPrintf(PETSC_COMM_WORLD, "   Topography      @ \n");
	if(pvsurf->amplitude)  PetscPrintf(PETSC_COMM_WORLD, "   Amplitude       @ \n");

	PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");

	// set file name
	sprintf(pvsurf->outfile, "%s_surf", filename);

	// create output buffer
	ierr = PVSurfCreateData(pvsurf); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVSurfCreateData"
PetscErrorCode PVSurfCreateData(PVSurf *pvsurf)
{
	FDSTAG   *fs;
	PetscInt  rx, ry, sx, sy, nx, ny;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// check activation
	if(!pvsurf->outsurf) PetscFunctionReturn(0);

	// access staggered grid layout
	fs = pvsurf->surf->jr->fs;

	// get local output grid sizes
	GET_OUTPUT_RANGE(rx, nx, sx, fs->dsx)
	GET_OUTPUT_RANGE(ry, ny, sy, fs->dsy)

	// buffer is only necessary on ranks zero in z direction
	if(!fs->dsz.rank)
	{
		// allocate output buffer
		ierr = PetscMalloc((size_t)(_max_num_comp_surf_*nx*ny)*sizeof(float), &pvsurf->buff); CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVSurfDestroy"
PetscErrorCode PVSurfDestroy(PVSurf *pvsurf)
{
	PetscFunctionBegin;

	// check activation
	if(!pvsurf->outsurf) PetscFunctionReturn(0);

	PetscFree(pvsurf->buff);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVSurfWriteTimeStep"
PetscErrorCode PVSurfWriteTimeStep(PVSurf *pvsurf, const char *dirName, PetscScalar ttime)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// check activation
	if(!pvsurf->outsurf) PetscFunctionReturn(0);

	// update .pvd file if necessary
	ierr = UpdatePVDFile(dirName, pvsurf->outfile, "pvts", &pvsurf->offset, ttime, pvsurf->outpvd); CHKERRQ(ierr);

	// write parallel data .pvts file
	ierr = PVSurfWritePVTS(pvsurf, dirName); CHKERRQ(ierr);

	// write sub-domain data .vts files
	ierr = PVSurfWriteVTS(pvsurf, dirName); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVSurfWritePVTS"
PetscErrorCode PVSurfWritePVTS(PVSurf *pvsurf, const char *dirName)
{
	FILE        *fp;
	FDSTAG      *fs;
	char        *fname;
	Scaling     *scal;
	PetscInt     nproc, rx, ry, rz;
	PetscMPIInt  iproc;

	PetscFunctionBegin;

	// only first process generates this file (WARNING! Bottleneck!)
	if(!ISRankZero(PETSC_COMM_WORLD)) PetscFunctionReturn(0);

	// access context
	fs   = pvsurf->surf->jr->fs;
	scal = pvsurf->surf->jr->scal;

	// open outfile.pvts file in the output directory (write mode)
	asprintf(&fname, "%s/%s.pvts", dirName, pvsurf->outfile);
	fp = fopen(fname,"wb");
	if(fp == NULL) SETERRQ1(PETSC_COMM_SELF, 1,"cannot open file %s", fname);
	free(fname);

	// write header
	WriteXMLHeader(fp, "PStructuredGrid");

	// open structured grid data block (write total grid size)
	fprintf(fp, "\t<PStructuredGrid GhostLevel=\"0\" WholeExtent=\"1 %lld 1 %lld 1 1\">\n",
		(LLD)fs->dsx.tnods,
		(LLD)fs->dsy.tnods);

	// write cell data block (empty)
	fprintf(fp, "\t\t<PCellData>\n");
	fprintf(fp, "\t\t</PCellData>\n");

	// write coordinate block
	fprintf(fp, "\t\t<PPoints>\n");

	fprintf(fp,"\t\t\t<PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\"/>\n");

	fprintf(fp, "\t\t</PPoints>\n");

	// write description of output vectors
	fprintf(fp, "\t\t<PPointData>\n");

	if(pvsurf->velocity)
	{
		fprintf(fp,"\t\t\t<PDataArray type=\"Float32\" Name=\"velocity %s\" NumberOfComponents=\"3\" format=\"appended\"/>\n", scal->lbl_velocity);
	}

	if(pvsurf->topography)
	{
		fprintf(fp,"\t\t\t<PDataArray type=\"Float32\" Name=\"topography %s\" NumberOfComponents=\"1\" format=\"appended\"/>\n", scal->lbl_length);
	}

	if(pvsurf->amplitude)
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
		fprintf(fp, "\t\t<Piece Extent=\"%lld %lld %lld %lld 1 1\" Source=\"%s_p%1.8lld.vts\"/>\n",
			(LLD)(fs->dsx.starts[rx] + 1), (LLD)(fs->dsx.starts[rx+1] + 1),
			(LLD)(fs->dsy.starts[ry] + 1), (LLD)(fs->dsy.starts[ry+1] + 1),
			pvsurf->outfile, (LLD)iproc);
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
PetscErrorCode PVSurfWriteVTS(PVSurf *pvsurf, const char *dirName)
{
	FILE      *fp;
	FDSTAG    *fs;
	Scaling   *scal;
	char      *fname;
	PetscInt   rx, ry, sx, sy, nx, ny;
	size_t     offset = 0;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// access context
	fs   = pvsurf->surf->jr->fs;
	scal = pvsurf->surf->jr->scal;

	fp = PETSC_NULL;
	// only ranks zero in z direction generate this file
	if(!fs->dsz.rank)
	{
		// open outfile_p_XXXXXX.vts file in the output directory (write mode)
		asprintf(&fname, "%s/%s_p%1.8lld.vts", dirName, pvsurf->outfile, (LLD)fs->dsz.color);
		fp = fopen(fname,"wb");
		if(fp == NULL) SETERRQ1(PETSC_COMM_SELF, 1,"cannot open file %s", fname);
		free(fname);

		// get sizes of output grid
		GET_OUTPUT_RANGE(rx, nx, sx, fs->dsx)
		GET_OUTPUT_RANGE(ry, ny, sy, fs->dsy)

		// write header
		WriteXMLHeader(fp, "StructuredGrid");

		// open structured grid data block (write total grid size)
		fprintf(fp, "\t<StructuredGrid WholeExtent=\"%lld %lld %lld %lld 1 1\">\n",
			(LLD)(fs->dsx.starts[rx] + 1), (LLD)(fs->dsx.starts[rx+1] + 1),
			(LLD)(fs->dsy.starts[ry] + 1), (LLD)(fs->dsy.starts[ry+1] + 1));

		// open sub-domain (piece) description block
		fprintf(fp, "\t\t<Piece Extent=\"%lld %lld %lld %lld 1 1\">\n",
			(LLD)(fs->dsx.starts[rx] + 1), (LLD)(fs->dsx.starts[rx+1] + 1),
			(LLD)(fs->dsy.starts[ry] + 1), (LLD)(fs->dsy.starts[ry+1] + 1));

		// write cell data block (empty)
		fprintf(fp, "\t\t\t<CellData>\n");
		fprintf(fp, "\t\t\t</CellData>\n");

		// write coordinate block
		fprintf(fp, "\t\t<Points>\n");

		fprintf(fp,"\t\t\t<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%lld\"/>\n",
			(LLD)offset);

		offset += sizeof(uint64_t) + sizeof(float)*(size_t)(nx*ny*3);

		fprintf(fp, "\t\t</Points>\n");

		// write description of output vectors
		fprintf(fp, "\t\t<PointData>\n");

		if(pvsurf->velocity)
		{
			fprintf(fp,"\t\t\t<DataArray type=\"Float32\" Name=\"velocity %s\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%lld\"/>\n",
				scal->lbl_velocity, (LLD)offset);

			offset += sizeof(uint64_t) + sizeof(float)*(size_t)(nx*ny*3);
		}

		if(pvsurf->topography)
		{
			fprintf(fp,"\t\t\t<DataArray type=\"Float32\" Name=\"topography %s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",
				scal->lbl_length, (LLD)offset);

			offset += sizeof(uint64_t) + sizeof(float)*(size_t)(nx*ny);
		}

		if(pvsurf->amplitude)
		{
			fprintf(fp,"\t\t\t<DataArray type=\"Float32\" Name=\"amplitude %s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",
				scal->lbl_length, (LLD)offset);

			offset += sizeof(uint64_t) + sizeof(float)*(size_t)(nx*ny);
		}

		fprintf(fp, "\t\t</PointData>\n");

		// close sub-domain and grid blocks
		fprintf(fp, "\t\t</Piece>\n");
		fprintf(fp, "\t</StructuredGrid>\n");

		// write appended data section
		fprintf(fp, "\t<AppendedData encoding=\"raw\">\n");
		fprintf(fp,"_");
	}

	// write point coordinates
	ierr = PVSurfWriteCoord (pvsurf, fp); CHKERRQ(ierr);

	// write output vectors
	if(pvsurf->velocity)   { ierr = PVSurfWriteVel      (pvsurf, fp); CHKERRQ(ierr); }
	if(pvsurf->topography) { ierr = PVSurfWriteTopo     (pvsurf, fp); CHKERRQ(ierr); }
	if(pvsurf->amplitude)  { ierr = PVSurfWriteAmplitude(pvsurf, fp); CHKERRQ(ierr); }

	if(!fs->dsz.rank)
	{
		// close appended data section and file
		fprintf(fp, "\n\t</AppendedData>\n");
		fprintf(fp, "</VTKFile>\n");

		// close file
		fclose(fp);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
void OutputBufferWrite(
	FILE     *fp,
	float    *buff,
	PetscInt  cn)
{
	if(!cn) return;

	// dump output buffer contents to disk
	uint64_t nbytes;

	// compute number of bytes
	nbytes = (uint64_t)cn*(int)sizeof(float);

	// dump number of bytes
	fwrite(&nbytes, sizeof(uint64_t), 1, fp);

	// dump buffer contents
	fwrite(buff, sizeof(float), (size_t)cn, fp);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVSurfWriteCoord"
PetscErrorCode PVSurfWriteCoord(PVSurf *pvsurf, FILE *fp)
{
	FreeSurf    *surf;
	FDSTAG      *fs;
	float       *buff;
	PetscScalar ***topo, cf;
	PetscInt    i, j, rx, ry, nx, ny, sx, sy, cn, L;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	L    = 0;
	cn   = 0;
	buff = pvsurf->buff;
	surf = pvsurf->surf;
	fs   = surf->jr->fs;
	cf   = surf->jr->scal->length;

	GET_OUTPUT_RANGE(rx, nx, sx, fs->dsx)
	GET_OUTPUT_RANGE(ry, ny, sy, fs->dsy)

	ierr = DMDAVecGetArray(surf->DA_SURF, surf->ltopo, &topo); CHKERRQ(ierr);

	if(!fs->dsz.rank)
	{
		START_PLANE_LOOP
		{
			// store node coordinates
			buff[cn++] = (float)(cf*COORD_NODE(i, sx, fs->dsx));
			buff[cn++] = (float)(cf*COORD_NODE(j, sy, fs->dsy));
			buff[cn++] = (float)(cf*topo[L][j][i]);
		}
		END_PLANE_LOOP
	}

	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->ltopo, &topo); CHKERRQ(ierr);

	OutputBufferWrite(fp, buff, cn);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVSurfWriteVel"
PetscErrorCode PVSurfWriteVel(PVSurf *pvsurf, FILE *fp)
{
	FreeSurf    *surf;
	FDSTAG      *fs;
	float       *buff;
	PetscScalar ***vx, ***vy, ***vz, cf;
	PetscInt    i, j, rx, ry, nx, ny, sx, sy, cn, L;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	L    = 0;
	cn   = 0;
	buff = pvsurf->buff;
	surf = pvsurf->surf;
	fs   = surf->jr->fs;
	cf   = surf->jr->scal->velocity;

	GET_OUTPUT_RANGE(rx, nx, sx, fs->dsx)
	GET_OUTPUT_RANGE(ry, ny, sy, fs->dsy)

	ierr = DMDAVecGetArray(surf->DA_SURF, surf->vx, &vx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(surf->DA_SURF, surf->vy, &vy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(surf->DA_SURF, surf->vz, &vz); CHKERRQ(ierr);

	if(!fs->dsz.rank)
	{
		START_PLANE_LOOP
		{
			// store surface velocities
			buff[cn++] = (float)(cf*vx[L][j][i]);
			buff[cn++] = (float)(cf*vy[L][j][i]);
			buff[cn++] = (float)(cf*vz[L][j][i]);
		}
		END_PLANE_LOOP
	}

	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->vx, &vx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->vy, &vy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->vz, &vz); CHKERRQ(ierr);

	OutputBufferWrite(fp, buff, cn);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVSurfWriteTopo"
PetscErrorCode PVSurfWriteTopo(PVSurf *pvsurf, FILE *fp)
{
	FreeSurf    *surf;
	FDSTAG      *fs;
	float       *buff;
	PetscScalar ***topo, cf;
	PetscInt    i, j, rx, ry, nx, ny, sx, sy, cn, L;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	L    = 0;
	cn   = 0;
	buff = pvsurf->buff;
	surf = pvsurf->surf;
	fs   = surf->jr->fs;
	cf   = surf->jr->scal->length;

	GET_OUTPUT_RANGE(rx, nx, sx, fs->dsx)
	GET_OUTPUT_RANGE(ry, ny, sy, fs->dsy)

	ierr = DMDAVecGetArray(surf->DA_SURF, surf->ltopo, &topo); CHKERRQ(ierr);

	if(!fs->dsz.rank)
	{
		START_PLANE_LOOP
		{
			// store surface topography
			buff[cn++] = (float)(cf*topo[L][j][i]);
		}
		END_PLANE_LOOP
	}

	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->ltopo, &topo); CHKERRQ(ierr);

	OutputBufferWrite(fp, buff, cn);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVSurfWriteAmplitude"
PetscErrorCode PVSurfWriteAmplitude(PVSurf *pvsurf, FILE *fp)
{
	FreeSurf    *surf;
	FDSTAG      *fs;
	float       *buff;
	PetscScalar ***topo, avg_topo, cf;
	PetscInt    i, j, rx, ry, nx, ny, sx, sy, cn, L;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	L    = 0;
	cn   = 0;
	buff = pvsurf->buff;
	surf = pvsurf->surf;
	fs   = surf->jr->fs;
	cf   = surf->jr->scal->length;

	// retrieve average topography
	avg_topo = surf->avg_topo;

	GET_OUTPUT_RANGE(rx, nx, sx, fs->dsx)
	GET_OUTPUT_RANGE(ry, ny, sy, fs->dsy)

	ierr = DMDAVecGetArray(surf->DA_SURF, surf->ltopo, &topo); CHKERRQ(ierr);

	if(!fs->dsz.rank)
	{
		START_PLANE_LOOP
		{
			// store topography amplitude
			buff[cn++] = (float)(cf*(topo[L][j][i] - avg_topo));
		}
		END_PLANE_LOOP
	}

	ierr = DMDAVecRestoreArray(surf->DA_SURF, surf->ltopo, &topo); CHKERRQ(ierr);

	OutputBufferWrite(fp, buff, cn);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

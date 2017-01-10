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
 **    filename:   AVDView.c
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
#include "LaMEM.h"
#include "parsing.h"
#include "scaling.h"
#include "fdstag.h"
#include "solVar.h"
#include "tssolve.h"
#include "bc.h"
#include "JacRes.h"
#include "interpolate.h"
#include "surf.h"
#include "advect.h"
#include "paraViewOutBin.h"
#include "AVDView.h"
#include "tools.h"
//---------------------------------------------------------------------------
// .......................... AVD ParaView Output ...........................
//---------------------------------------------------------------------------

#undef __FUNCT__
#define __FUNCT__ "PVAVDCreate"
PetscErrorCode PVAVDCreate(PVAVD *pvavd, AdvCtx *actx, const char *filename)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	ierr = PVAVDReadFromOptions(pvavd); CHKERRQ(ierr);

	if(!pvavd->outavd) PetscFunctionReturn(0);

	pvavd->actx = actx;

	asprintf(&pvavd->outfile, "%s_phase", filename);

	pvavd->offset = 0;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVAVDDestroy"
PetscErrorCode PVAVDDestroy(PVAVD *pvavd)
{
	PetscFunctionBegin;

	if(!pvavd->outavd) PetscFunctionReturn(0);

	free(pvavd->outfile);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVAVDReadFromOptions"
PetscErrorCode PVAVDReadFromOptions(PVAVD *pvavd)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	pvavd->outavd = 0; // AVD output flag

	ierr = PetscOptionsGetInt(NULL, NULL, "-out_avd", &pvavd->outavd, NULL); CHKERRQ(ierr);

	if(!pvavd->outavd) PetscFunctionReturn(0);

	pvavd->refine = 2; // Voronoi Diagram refinement factor
	pvavd->outpvd = 0; // pvd file output flag

	ierr = PetscOptionsGetInt(NULL, NULL, "-out_avd_ref", &pvavd->refine, NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, NULL, "-out_avd_pvd", &pvavd->outpvd, NULL); CHKERRQ(ierr);

	if(pvavd->outpvd)
	{
		PetscPrintf(PETSC_COMM_WORLD, " Writing AVD3D .pvd file to disk\n");
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVAVDWriteTimeStep"
PetscErrorCode PVAVDWriteTimeStep(PVAVD *pvavd, const char *dirName, PetscScalar ttime, PetscInt tindx)
{
	// Create a 3D Voronoi diagram from particles with phase information
	// write the file to disk and perform scaling/unscaling of the variables

	PetscErrorCode ierr;
	PetscFunctionBegin;

	if(!pvavd->outavd) PetscFunctionReturn(0);

	// update .pvd file if necessary
	if(pvavd->outpvd)
	{
		ierr = UpdatePVDFile(dirName, pvavd->outfile, "pvtr", &pvavd->offset, ttime, tindx); CHKERRQ(ierr);
	}

	ierr = PVAVDWritePVTR(pvavd, dirName); CHKERRQ(ierr);

	ierr = PVAVDWriteVTR(pvavd, dirName); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVAVDWritePVTR"
PetscErrorCode PVAVDWritePVTR(PVAVD *pvavd, const char *dirName)
{
	FILE        *fp;
	char        *fname;
	PetscMPIInt inproc, irank;
	PetscInt    r2d, p, pi, pj, pk, nproc, rank;

	PetscFunctionBegin;

/*
	// only first process generates this file (WARNING! Bottleneck!)
	if(!ISRankZero(PETSC_COMM_WORLD)) PetscFunctionReturn(0);

	MPI_Comm_size(PETSC_COMM_WORLD, &inproc); nproc = (PetscInt)inproc;
	MPI_Comm_rank(PETSC_COMM_WORLD, &irank);  rank  = (PetscInt)irank;

	// open outfile.pvts file in the output directory (write mode)
	asprintf(&fname, "%s/%s.pvtr", dirName, pvavd->outfile);
	fp = fopen(fname,"w");
	if(fp == NULL) SETERRQ1(PETSC_COMM_SELF, 1,"cannot open file %s", fname);
	free(fname);

	pk  = rank/(A->M*A->N);
	r2d = rank - pk*(A->M*A->N);
	pj  = r2d/(A->M);
	pi  = r2d - pj*A->M;

	WriteXMLHeader(fp, "PRectilinearGrid");

	fprintf(fp, "  <PRectilinearGrid WholeExtent=\"%lld %lld %lld %lld %lld %lld\" GhostLevel=\"0\" >\n",
		0LL,(LLD)(A->gmx),
		0LL,(LLD)(A->gmy),
		0LL,(LLD)(A->gmz));

	fprintf(fp, "    <PCoordinates>\n");
	fprintf(fp, "      <PDataArray type=\"Float32\" NumberOfComponents=\"1\" format=\"appended\" />\n");
	fprintf(fp, "      <PDataArray type=\"Float32\" NumberOfComponents=\"1\" format=\"appended\" />\n");
	fprintf(fp, "      <PDataArray type=\"Float32\" NumberOfComponents=\"1\" format=\"appended\" />\n");
	fprintf(fp, "    </PCoordinates>\n");

	fprintf(fp, "    <PCellData>\n");

	fprintf(fp, "      <PDataArray type=\"UInt8\" Name=\"phase\" NumberOfComponents=\"1\" format=\"appended\" />\n");

	fprintf(fp, "    </PCellData>\n");

	fprintf(fp, "    <PPointData>\n");
	fprintf(fp, "    </PPointData>\n");

	for(p=0; p<nproc; p++)
	{
		pk = p/(A->M*A->N);
		r2d = p - pk*(A->M*A->N);
		pj = r2d/(A->M);
		pi = r2d - pj*A->M;

		fprintf(fp, "    <Piece Extent=\"%lld %lld %lld %lld %lld %lld\" Source=\"%s_p%1.6lld.vtr\" />\n",
				(LLD)(A->ownership_ranges_i[pi]),(LLD)(A->ownership_ranges_i[pi+1]),
				(LLD)(A->ownership_ranges_j[pj]),(LLD)(A->ownership_ranges_j[pj+1]),
				(LLD)(A->ownership_ranges_k[pk]),(LLD)(A->ownership_ranges_k[pk+1]),
				pvavd->outfile, (LLD)p );
	}

	fprintf(fp, "  </PRectilinearGrid>\n");

	fprintf(fp, "</VTKFile>\n");

	fclose( fp );
*/

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVAVDWriteVTR"
PetscErrorCode PVAVDWriteVTR(PVAVD *pvavd, const char *dirName)
{
	// WARNING! writing single entry at a time is too slow. Use buffers instead!

	PetscMPIInt   irank;
	FILE          *fp;
	char          *fname;
	PetscInt      i, j, k, ii;
	PetscInt      r2d, pi, pj, pk, rank;
	PetscScalar   chLen;
	float         crd;
	unsigned char phase;
	int           offset, L;

	PetscFunctionBegin;

/*
	// access context
	chLen = pvavd->actx->jr->scal->length;

	MPI_Comm_rank(PETSC_COMM_WORLD, &irank);  rank = (PetscInt)irank;

	// open outfile_p_XXXXXX.vtr file in the output directory (write mode)
	asprintf(&fname, "%s/%s_p%1.6lld.vtr", dirName, pvavd->outfile, (LLD)rank);
	fp = fopen(fname,"w");
	if(fp == NULL) SETERRQ1(PETSC_COMM_SELF, 1,"cannot open file %s", fname);
	free(fname);

	pk  = rank/(A->M*A->N);
	r2d = rank - pk*(A->M*A->N);
	pj  = r2d/(A->M);
	pi  = r2d - pj*A->M;

	// write header
	WriteXMLHeader(fp, "RectilinearGrid");

  fprintf(fp, "  <RectilinearGrid WholeExtent=\"%lld %lld %lld %lld %lld %lld\" >\n",
		  (LLD)(A->ownership_ranges_i[pi]),(LLD)(A->ownership_ranges_i[pi+1]),
		  (LLD)(A->ownership_ranges_j[pj]),(LLD)(A->ownership_ranges_j[pj+1]),
		  (LLD)(A->ownership_ranges_k[pk]),(LLD)(A->ownership_ranges_k[pk+1]));

	fprintf(fp, "    <Piece Extent=\"%lld %lld %lld %lld %lld %lld\" >\n",
			(LLD)(A->ownership_ranges_i[pi]),(LLD)(A->ownership_ranges_i[pi+1]),
			(LLD)(A->ownership_ranges_j[pj]),(LLD)(A->ownership_ranges_j[pj+1]),
			(LLD)(A->ownership_ranges_k[pk]),(LLD)(A->ownership_ranges_k[pk+1]));

	offset = 0;

	fprintf(fp, "    <Coordinates>\n");

	// X
	fprintf(fp, "      <DataArray type=\"Float32\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)offset);
	offset += (int)(sizeof(int) + sizeof(float)*(size_t)(A->mx+1));
	// Y
	fprintf(fp, "      <DataArray type=\"Float32\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)offset);
	offset += (int)(sizeof(int) + sizeof(float)*(size_t)(A->my+1));
	// Z
	fprintf(fp, "      <DataArray type=\"Float32\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)offset);
	offset += (int)(sizeof(int) + sizeof(float)*(size_t)(A->mz+1));

	fprintf(fp, "    </Coordinates>\n");

	fprintf(fp, "    <CellData>\n");

	// phase
	fprintf(fp, "      <DataArray type=\"UInt8\" Name=\"phase\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)offset);

	fprintf(fp, "    </CellData>\n");

	fprintf(fp, "    <PointData>\n");
	fprintf(fp, "    </PointData>\n");

	fprintf(fp, "    </Piece>\n");
  fprintf(fp, "  </RectilinearGrid>\n");

	fprintf(fp,"  <AppendedData encoding=\"raw\">\n");
	fprintf(fp,"_");

	// X
	L = (int)sizeof(float)*(int)(A->mx+1);
	fwrite(&L, sizeof(int), 1, fp);
	for( i=0; i<A->mx+1; i++ ) {
		crd = (float)((A->x0 + (PetscScalar)i*A->dx)*chLen);
		fwrite(&crd,sizeof(float),1,fp);
	}

	// Y
	L = (int)sizeof(float)*(int)(A->my+1);
	fwrite(&L, sizeof(int), 1, fp);
	for( i=0; i<A->my+1; i++ ) {
		crd = (float)((A->y0 + (PetscScalar)i*A->dy)*chLen);
		fwrite(&crd,sizeof(float),1,fp);
	}

	// Z
	L = (int)sizeof(float)*(int)(A->mz+1);
	fwrite(&L, sizeof(int), 1, fp);
	for( i=0; i<A->mz+1; i++ ) {
		crd = (float)((A->z0 + (PetscScalar)i*A->dz)*chLen);
		fwrite(&crd,sizeof(float),1,fp);
	}

	// phase
	L = (int)sizeof(unsigned char)*(int)(A->mz*A->my*A->mx);
	fwrite(&L, sizeof(int), 1, fp);
	for (k=1; k<A->mz+1; k++) {
		for (j=1; j<A->my+1; j++) {
			for (i=1; i<A->mx+1; i++)
			{
				ii    = i + j*A->mx_mesh + k*A->mx_mesh*A->my_mesh;
				phase = (unsigned char)A->points[A->cells[ii].p].phase;
				fwrite(&phase,sizeof(unsigned char),1,fp);
			}
		}
	}
	fprintf(fp,"\n  </AppendedData>\n");

	fprintf(fp, "</VTKFile>\n");

	fclose( fp );
*/

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------



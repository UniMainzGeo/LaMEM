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
 **    filename:   paraViewOutMark.c
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
//..............   MARKER PARAVIEW XML OUTPUT ROUTINES   ....................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "parsing.h"
#include "scaling.h"
#include "tssolve.h"
#include "tools.h"
#include "fdstag.h"
#include "solVar.h"
#include "bc.h"
#include "JacRes.h"
#include "interpolate.h"
#include "surf.h"
#include "advect.h"
#include "paraViewOutBin.h"
#include "paraViewOutMark.h"

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVMarkCreate"
PetscErrorCode PVMarkCreate(PVMark *pvmark, AdvCtx *actx, const char *filename)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// set context
	pvmark->actx = actx;

	// read options
	ierr = PVMarkReadFromOptions(pvmark); CHKERRQ(ierr);

	// set file name
	asprintf(&pvmark->outfile, "%s_mark", filename);

	// set .pvd file offset
	pvmark->offset = 0;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVMarkDestroy"
PetscErrorCode PVMarkDestroy(PVMark *pvmark)
{
	PetscFunctionBegin;

	// file name
	free(pvmark->outfile);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVMarkReadFromOptions"
PetscErrorCode PVMarkReadFromOptions(PVMark *pvmark)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// set output mask
	pvmark->outmark    = 0;
	pvmark->outpvd     = 0;

	// read output flags
	ierr = PetscOptionsGetInt(NULL, NULL, "-out_markers",         &pvmark->outmark,    NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(NULL, NULL, "-out_mark_pvd",        &pvmark->outpvd,     NULL); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVMarkWriteTimeStep"
PetscErrorCode PVMarkWriteTimeStep(PVMark *pvmark, const char *dirName, PetscScalar ttime, PetscInt tindx)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// return if not output
	if(pvmark->outmark==0) PetscFunctionReturn(0);

	// update .pvd file if necessary
	if(pvmark->outpvd)
	{
		ierr = UpdatePVDFile(dirName, pvmark->outfile, "pvtu", &pvmark->offset, ttime, tindx); CHKERRQ(ierr);
	}

	// write parallel data .pvtu file
	ierr = PVMarkWritePVTU(pvmark, dirName); CHKERRQ(ierr);

	// write sub-domain data .vtu files
	ierr = PVMarkWriteVTU(pvmark, dirName); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVMarkWriteVTU"
PetscErrorCode PVMarkWriteVTU(PVMark *pvmark, const char *dirName)
{
	// output markers in .vtu files
	AdvCtx     *actx;
	char       *fname;
	FILE       *fp;
	PetscInt    i, idx, connect, length, phase;
	PetscScalar scal_length;
	float       xp[3];
	size_t      offset = 0;

	PetscFunctionBegin;

	// get context
	actx = pvmark->actx;

	// create file name
	asprintf(&fname, "%s/%s_p%1.8lld.vtu", dirName, pvmark->outfile, (LLD)actx->iproc);

	// open file
	fp = fopen( fname, "w" );
	if(fp == NULL) SETERRQ1(PETSC_COMM_SELF, 1,"cannot open file %s", fname);
	free(fname);

	// write header
	WriteXMLHeader(fp, "UnstructuredGrid");

	// initialize connectivity
	connect = actx->nummark;

	// begin unstructured grid
	fprintf( fp, "\t<UnstructuredGrid>\n" );
	fprintf( fp, "\t\t<Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%lld\">\n",(LLD)actx->nummark,(LLD)connect );

	// cells
	fprintf( fp, "\t\t\t<Cells>\n");

	// connectivity
	fprintf( fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)offset);
	offset += sizeof(int) + sizeof(int)*(size_t)connect;

	// offsets
	fprintf( fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)offset);
	offset += sizeof(int) + sizeof(int)*(size_t)connect;
	// types
	fprintf( fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)offset);
	offset += sizeof(int) + sizeof(int)*(size_t)connect;

	fprintf( fp, "\t\t\t</Cells>\n");

	// write cell data block (empty)
	fprintf( fp, "\t\t\t<CellData>\n");
	fprintf( fp, "\t\t\t</CellData>\n");

	// points
	fprintf( fp, "\t\t\t<Points>\n");

	// point coordinates
	fprintf( fp, "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%lld\" />\n",(LLD)offset);
	offset += sizeof(int) + sizeof(float)*(size_t)(actx->nummark*3);

	fprintf( fp, "\t\t\t</Points>\n");

	// point data - marker phase
	fprintf( fp, "\t\t\t<PointData Scalars=\"\">\n");

	fprintf( fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"Phase\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)offset );
	offset += sizeof(int) + sizeof(int)*(size_t)actx->nummark;

	fprintf( fp, "\t\t\t</PointData>\n");

	fprintf( fp, "\t\t</Piece>\n");
	fprintf( fp, "\t</UnstructuredGrid>\n");

	fprintf( fp,"\t<AppendedData encoding=\"raw\">\n");
	fprintf( fp,"_");

	// -------------------
	// write connectivity
	// -------------------
	length = (int)sizeof(int)*connect;
	fwrite( &length,sizeof(int),1, fp);

	for( i = 0; i < connect; i++)
	{
		idx = i;
		fwrite( &idx, sizeof(int),1, fp );
	}
	// -------------------
	// write offsets
	// -------------------
	length = (int)sizeof(int)*connect;
	fwrite( &length,sizeof(int),1, fp);

	for( i = 0; i < connect; i++)
	{
		idx = i+1;
		fwrite( &idx, sizeof(int),1, fp );
	}
	// -------------------
	// write types
	// -------------------
	length = (int)sizeof(int)*connect;
	fwrite( &length,sizeof(int),1, fp);

	for( i = 0; i < connect; i++)
	{
		idx = 1;
		fwrite( &idx, sizeof(int),1, fp );
	}
	// -------------------
	// write point coordinates
	// -------------------
	length = (int)sizeof(float)*(3*actx->nummark);
	fwrite( &length,sizeof(int),1, fp);

	// scaling length
	scal_length = actx->jr->scal->length;

	for( i = 0; i < actx->nummark; i++)
	{
		xp[0] = (float)(actx->markers[i].X[0]*scal_length);
		xp[1] = (float)(actx->markers[i].X[1]*scal_length);
		xp[2] = (float)(actx->markers[i].X[2]*scal_length);
		fwrite( xp, sizeof(float), (size_t)3, fp );
	}
	// -------------------
	// write field: phases
	// -------------------
	length = (int)sizeof(int)*(actx->nummark);
	fwrite( &length,sizeof(int),1, fp);

	for( i = 0; i < actx->nummark; i++)
	{
		phase = actx->markers[i].phase;
		fwrite( &phase, sizeof(int),1, fp );
	}
	// -------------------

	// end header
	fprintf( fp,"\n\t</AppendedData>\n");
	fprintf( fp, "</VTKFile>\n");

	// close file
	fclose(fp);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVMarkWritePVTU"
PetscErrorCode PVMarkWritePVTU(PVMark *pvmark, const char *dirName)
{
	// create .pvtu file for marker output
	// load the pvtu file in ParaView and apply a Glyph-spheres filter
	AdvCtx   *actx;
	char     *fname;
	FILE     *fp;
	PetscInt i;

	PetscFunctionBegin;

	// only processor 0
	if (!ISRankZero(PETSC_COMM_WORLD)) { PetscFunctionReturn(0); }

	// get context
	actx = pvmark->actx;

	// create file name
	asprintf(&fname, "%s/%s.pvtu", dirName, pvmark->outfile);

	// open file
	fp = fopen( fname, "w" );
	if(fp == NULL) SETERRQ1(PETSC_COMM_SELF, 1,"cannot open file %s", fname);
	free(fname);

	// write header
	WriteXMLHeader(fp, "PUnstructuredGrid");

	// define ghost level
	fprintf( fp, "\t<PUnstructuredGrid GhostLevel=\"0\">\n" );

	// cell data (empty)
	fprintf( fp, "\t\t<PCellData>\n");
	fprintf( fp, "\t\t</PCellData>\n");

	// cells
	fprintf( fp, "\t\t\t<Cells>\n");
	fprintf( fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" />\n");
	fprintf( fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" />\n");
	fprintf( fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"appended\" />\n");
	fprintf( fp, "\t\t\t</Cells>\n");

	// points
	fprintf( fp, "\t\t<PPoints>\n");
	fprintf( fp, "\t\t\t<PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\"/>\n");
	fprintf( fp, "\t\t</PPoints>\n");

	// point data
	fprintf( fp, "\t\t<PPointData>\n");
		fprintf(fp,"\t\t\t<PDataArray type=\"Int32\" Name=\"Phase\" NumberOfComponents=\"1\" format=\"appended\"/>\n");
	fprintf( fp, "\t\t</PPointData>\n");

	for(i = 0; i < actx->nproc; i++){
		fprintf( fp, "\t\t<Piece Source=\"%s_p%1.8lld.vtu\"/>\n",pvmark->outfile,(LLD)i);
	}

	// close the file
	fprintf( fp, "\t</PUnstructuredGrid>\n");
	fprintf( fp, "</VTKFile>\n");

	// close file and free name
	fclose(fp);

	PetscFunctionReturn(0);
}

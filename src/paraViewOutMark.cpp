/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : LaMEM
 **   License      : MIT, see LICENSE file for details
 **   Contributors : Anton Popov, Boris Kaus, see AUTHORS file for complete list
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : kaus@uni-mainz.de, popov@uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
//---------------------------------------------------------------------------
//..............   MARKER PARAVIEW XML OUTPUT ROUTINES   ....................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "paraViewOutMark.h"
#include "paraViewOutBin.h"
#include "parsing.h"
#include "scaling.h"
#include "advect.h"
#include "JacRes.h"
#include "tools.h"
//---------------------------------------------------------------------------
PetscErrorCode PVMarkCreate(PVMark *pvmark, FB *fb)
{
	char filename[_str_len_];

	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// check advection type
	if(pvmark->actx->advect == ADV_NONE) PetscFunctionReturn(0);

	// check activation
	ierr = getIntParam(fb, _OPTIONAL_, "out_mark", &pvmark->outmark, 1, 1); CHKERRQ(ierr);

	if(!pvmark->outmark) PetscFunctionReturn(0);

	// initialize
	pvmark->outpvd = 1;

	// read
	ierr = getStringParam(fb, _OPTIONAL_, "out_file_name", filename,    "output"); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_mark_pvd",  &pvmark->outpvd, 1, 1); CHKERRQ(ierr);

	// print summary
	PetscPrintf(PETSC_COMM_WORLD, "Marker output parameters:\n");
	PetscPrintf(PETSC_COMM_WORLD, "   Write .pvd file : %s \n", pvmark->outpvd ? "yes" : "no");
	PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");

	// set file name
	sprintf(pvmark->outfile, "%s_mark", filename);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PVMarkWriteTimeStep(PVMark *pvmark, const char *dirName, PetscScalar ttime)
{
	PetscErrorCode ierr;
	PetscFunctionBeginUser;

	// check activation
	if(!pvmark->outmark) PetscFunctionReturn(0);

	// update .pvd file if necessary
	ierr = UpdatePVDFile(dirName, pvmark->outfile, "pvtu", &pvmark->offset, ttime, pvmark->outpvd); CHKERRQ(ierr);

	// write parallel data .pvtu file
	ierr = PVMarkWritePVTU(pvmark, dirName); CHKERRQ(ierr);

	// write sub-domain data .vtu files
	ierr = PVMarkWriteVTU(pvmark, dirName); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PVMarkWriteVTU(PVMark *pvmark, const char *dirName)
{
	// output markers in .vtu files
	AdvCtx     *actx;
	char       *fname;
	FILE       *fp;
	PetscInt    i, idx, connect, phase;
	uint64_t	length;
	PetscScalar scal_length;
	float       xp[3];
	size_t      offset = 0;

	PetscFunctionBeginUser;

	// get context
	actx = pvmark->actx;

	// create file name
	asprintf(&fname, "%s/%s_p%1.8lld.vtu", dirName, pvmark->outfile, (LLD)actx->iproc);

	// open file
	fp = fopen( fname, "wb" );
	if(fp == NULL) SETERRQ(PETSC_COMM_SELF, 1,"cannot open file %s", fname);
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
	offset += sizeof(uint64_t) + sizeof(int)*(size_t)connect;

	// offsets
	fprintf( fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)offset);
	offset += sizeof(uint64_t) + sizeof(int)*(size_t)connect;
	// types
	fprintf( fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)offset);
	offset += sizeof(uint64_t) + sizeof(int)*(size_t)connect;

	fprintf( fp, "\t\t\t</Cells>\n");

	// write cell data block (empty)
	fprintf( fp, "\t\t\t<CellData>\n");
	fprintf( fp, "\t\t\t</CellData>\n");

	// points
	fprintf( fp, "\t\t\t<Points>\n");

	// point coordinates
	fprintf( fp, "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%lld\" />\n",(LLD)offset);
	offset += sizeof(uint64_t) + sizeof(float)*(size_t)(actx->nummark*3);

	fprintf( fp, "\t\t\t</Points>\n");

	// point data - marker phase
	fprintf( fp, "\t\t\t<PointData Scalars=\"\">\n");

	fprintf( fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"Phase\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)offset );
	offset += sizeof(uint64_t) + sizeof(int)*(size_t)actx->nummark;

	fprintf( fp, "\t\t\t</PointData>\n");

	fprintf( fp, "\t\t</Piece>\n");
	fprintf( fp, "\t</UnstructuredGrid>\n");

	fprintf( fp,"\t<AppendedData encoding=\"raw\">\n");
	fprintf( fp,"_");

	// -------------------
	// write connectivity
	// -------------------
	length = (uint64_t)sizeof(int)*connect;
	fwrite( &length,sizeof(uint64_t),1, fp);

	for( i = 0; i < connect; i++)
	{
		idx = i;
		fwrite( &idx, sizeof(int),1, fp );
	}
	// -------------------
	// write offsets
	// -------------------
	length = (uint64_t)sizeof(int)*connect;
	fwrite( &length,sizeof(uint64_t),1, fp);

	for( i = 0; i < connect; i++)
	{
		idx = i+1;
		fwrite( &idx, sizeof(int),1, fp );
	}
	// -------------------
	// write types
	// -------------------
	length = (uint64_t)sizeof(int)*connect;
	fwrite( &length,sizeof(uint64_t),1, fp);

	for( i = 0; i < connect; i++)
	{
		idx = 1;
		fwrite( &idx, sizeof(int),1, fp );
	}
	// -------------------
	// write point coordinates
	// -------------------
	length = (uint64_t)sizeof(float)*(3*actx->nummark);
	fwrite( &length,sizeof(uint64_t),1, fp);

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
	length = (uint64_t)sizeof(int)*(actx->nummark);
	fwrite( &length,sizeof(uint64_t),1, fp);

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
PetscErrorCode PVMarkWritePVTU(PVMark *pvmark, const char *dirName)
{
	// create .pvtu file for marker output
	// load the pvtu file in ParaView and apply a Glyph-spheres filter
	AdvCtx   *actx;
	char     *fname;
	FILE     *fp;
	PetscInt i;

	PetscFunctionBeginUser;

	// only processor 0
	if (!ISRankZero(PETSC_COMM_WORLD)) { PetscFunctionReturn(0); }

	// get context
	actx = pvmark->actx;

	// create file name
	asprintf(&fname, "%s/%s.pvtu", dirName, pvmark->outfile);

	// open file
	fp = fopen( fname, "wb" );
	if(fp == NULL) SETERRQ(PETSC_COMM_SELF, 1,"cannot open file %s", fname);
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

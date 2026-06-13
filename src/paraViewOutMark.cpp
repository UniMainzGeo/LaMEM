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
#include "Tensor.h"
#include "advect.h"
#include "JacRes.h"
#include "tools.h"
//---------------------------------------------------------------------------
PetscErrorCode PVMarkCreate(PVMark *pvmark, FB *fb)
{
	char filename[_str_len_];

	
	PetscFunctionBeginUser;

	// check advection type
	if(pvmark->actx->advect == ADV_NONE) PetscFunctionReturn(0);

	// check activation
	PetscCall(getIntParam(fb, _OPTIONAL_, "out_mark", &pvmark->outmark, 1, 1));

	if(!pvmark->outmark) PetscFunctionReturn(0);

	// initialize
	pvmark->outpvd = 1;

	// read
	PetscCall(getStringParam(fb, _OPTIONAL_, "out_file_name", filename,    "output"));
	PetscCall(getIntParam   (fb, _OPTIONAL_, "out_mark_pvd",  &pvmark->outpvd, 1, 1));

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
	
	PetscFunctionBeginUser;

	// check activation
	if(!pvmark->outmark) PetscFunctionReturn(0);

	// update .pvd file if necessary
	PetscCall(UpdatePVDFile(dirName, pvmark->outfile, "pvtu", &pvmark->offset, ttime, pvmark->outpvd));

	// write parallel data .pvtu file
	PetscCall(PVMarkWritePVTU(pvmark, dirName));

	// write sub-domain data .vtu files
	PetscCall(PVMarkWriteVTU(pvmark, dirName));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PVMarkWriteVTU(PVMark *pvmark, const char *dirName)
{
	// output markers in .vtu files
	AdvCtx     *actx;
	char       *fname;
	FILE       *fp;
	PetscInt    i, nummark;
	PetscScalar scal_length;
	float       var_float, Xp[3];
	int         var_int;
	uint64_t    offset = 0;
	uint64_t 	length;

	PetscFunctionBeginUser;

	// get context
	actx = pvmark->actx;

	// create file name
	asprintf(&fname, "%s/%s_p%1.8" PetscInt_FMT ".vtu", dirName, pvmark->outfile, actx->iproc);

	// open file
	fp = fopen( fname, "wb" );
	if(fp == NULL) SETERRQ(PETSC_COMM_SELF, 1,"cannot open file %s", fname);
	free(fname);

	// write header
	WriteXMLHeader(fp, "UnstructuredGrid");

	// initialize connectivity
	nummark = actx->nummark;

	// begin unstructured grid
	fprintf(fp, "\t<UnstructuredGrid>\n" );
	fprintf(fp, "\t\t<Piece NumberOfPoints=\"%" PetscInt_FMT "\" NumberOfCells=\"%" PetscInt_FMT "\">\n", nummark, nummark);

	// cells
	fprintf(fp, "\t\t\t<Cells>\n");

	// connectivity
	fprintf(fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"%" PRIu64 "\"/>\n", offset);
	offset += (uint64_t)(sizeof(uint64_t) + sizeof(int)*(size_t)nummark);

	// offsets
	fprintf(fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"%" PRIu64 "\"/>\n", offset);
	offset += (uint64_t)(sizeof(uint64_t) + sizeof(int)*(size_t)nummark);

	// types
	fprintf(fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"appended\" offset=\"%" PRIu64 "\"/>\n", offset);
	offset += (uint64_t)(sizeof(uint64_t) + sizeof(int)*(size_t)nummark);

	fprintf(fp, "\t\t\t</Cells>\n");

	// write cell data block (empty)
	fprintf(fp, "\t\t\t<CellData>\n");
	fprintf(fp, "\t\t\t</CellData>\n");

	// points
	fprintf(fp, "\t\t\t<Points>\n");

	// point coordinates
	fprintf(fp, "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%" PRIu64 "\" />\n",offset);
	offset += (uint64_t)(sizeof(uint64_t) + sizeof(float)*(size_t)(3*nummark));

	fprintf(fp, "\t\t\t</Points>\n");

	// point data - marker phase
	fprintf(fp, "\t\t\t<PointData Scalars=\"\">\n");

	fprintf(fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"Phase\" format=\"appended\" offset=\"%" PRIu64 "\"/>\n", offset );
	offset += (uint64_t)(sizeof(uint64_t) + sizeof(int)*(size_t)nummark);
	
	fprintf(fp, "\t\t\t\t<DataArray type=\"Float32\" Name=\"APS\" format=\"appended\" offset=\"%" PRIu64 "\"/>\n", offset );
	offset += (uint64_t)(sizeof(uint64_t) + sizeof(float)*(size_t)nummark);

	fprintf(fp, "\t\t\t</PointData>\n");

	fprintf(fp, "\t\t</Piece>\n");
	fprintf(fp, "\t</UnstructuredGrid>\n");

	fprintf(fp,"\t<AppendedData encoding=\"raw\">\n");
	fprintf(fp,"_");

	// -------------------
	// write connectivity
	// -------------------
	length = (uint64_t)(sizeof(int)*(size_t)nummark);
	fwrite(&length, sizeof(uint64_t), 1, fp);

	for(i = 0; i < nummark; i++)
	{
		var_int = (int)i;
		fwrite(&var_int, sizeof(int), 1, fp );
	}

	// -------------------
	// write offsets
	// -------------------
	length = (uint64_t)(sizeof(int)*(size_t)nummark);
	fwrite(&length, sizeof(uint64_t), 1, fp);

	for(i = 0; i < nummark; i++)
	{
		var_int = (int)(i+1);
		fwrite(&var_int, sizeof(int), 1, fp );
	}

	// -------------------
	// write types
	// -------------------
	length = (uint64_t)(sizeof(int)*(size_t)nummark);
	fwrite(&length, sizeof(uint64_t), 1, fp);

	for(i = 0; i < nummark; i++)
	{
		var_int = 1;
		fwrite(&var_int, sizeof(int), 1, fp);
	}

	// -------------------
	// write point coordinates
	// -------------------
	length = (uint64_t)(sizeof(float)*(size_t)(3*nummark));
	fwrite(&length, sizeof(uint64_t), 1, fp);

	// scaling length
	scal_length = actx->jr->scal->length;

	for(i = 0; i < nummark; i++)
	{
		Xp[0] = (float)(actx->markers[i].X[0]*scal_length);
		Xp[1] = (float)(actx->markers[i].X[1]*scal_length);
		Xp[2] = (float)(actx->markers[i].X[2]*scal_length);

		fwrite(Xp, sizeof(float), 3, fp);
	}

	//--------------------
	// write field: phases
	//--------------------
	length = (uint64_t)(sizeof(int)*(size_t)nummark);
	fwrite(&length, sizeof(uint64_t), 1, fp);

	for(i = 0; i < nummark; i++)
	{
		var_int = (int)actx->markers[i].phase;
		fwrite(&var_int, sizeof(int), 1, fp);
	}

	//-----------------
	// write field: APS
	//-----------------
	length = (uint64_t)(sizeof(float)*(size_t)nummark);
	fwrite(&length, sizeof(uint64_t), 1, fp);

	for(i = 0; i < nummark; i++)
	{
		var_float = float(actx->markers[i].APS);
		fwrite(&var_float, sizeof(float),1, fp );
	}

	// end header
	fprintf(fp,"\n\t</AppendedData>\n");
	fprintf(fp, "</VTKFile>\n");

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
	fprintf(fp,"\t\t\t<PDataArray type=\"Float32\" Name=\"APS\" NumberOfComponents=\"1\" format=\"appended\"/>\n");
	fprintf( fp, "\t\t</PPointData>\n");

	for(i = 0; i < actx->nproc; i++)
	{
		fprintf( fp, "\t\t<Piece Source=\"%s_p%1.8" PetscInt_FMT ".vtu\"/>\n",pvmark->outfile,i);
	}

	// close the file
	fprintf( fp, "\t</PUnstructuredGrid>\n");
	fprintf( fp, "</VTKFile>\n");

	// close file and free name
	fclose(fp);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------


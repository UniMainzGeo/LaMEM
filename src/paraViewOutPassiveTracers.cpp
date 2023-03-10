/*
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
#include "paraViewOutBin.h"
#include "parsing.h"
#include "scaling.h"
#include "advect.h"
#include "JacRes.h"
#include "paraViewOutPassiveTracers.h"
#include "tools.h"
#include "passive_tracer.h"

//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVPtrCreate"
PetscErrorCode PVPtrCreate(PVPtr *pvptr, FB *fb)
{
	char filename[_str_len_];

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// check activation
	ierr = getIntParam(fb, _OPTIONAL_, "out_ptr", &pvptr->actx->jr->ctrl.Passive_Tracer, 1, 1); CHKERRQ(ierr);


	// check advection type
	if (pvptr->actx->jr->ctrl.Passive_Tracer == 0) PetscFunctionReturn(0);

	// check activation

    // Default Values ( if passive tracers are used, there is no point to let someone to not
	// visualize them)
	pvptr->ID          = 1;
	pvptr->Pressure    = 1;
	pvptr->Temperature = 1;
	pvptr->outptr      = 1;
	pvptr->outpvd      = 1;
	// read
	ierr = getStringParam(fb, _OPTIONAL_, "out_file_name",           filename,    "output"); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_ptr_ID",              &pvptr->ID,   1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_ptr_Temperature",     &pvptr->Temperature, 1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_ptr_Pressure",        &pvptr->Pressure,  1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_ptr_phase",           &pvptr->Phase,   1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_ptr_MeltFraction",    &pvptr->MeltFraction, 1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_ptr_Active",          &pvptr->Active   , 1, 1); CHKERRQ(ierr);
	ierr = getIntParam   (fb, _OPTIONAL_, "out_ptr_Grid_Mf",          &pvptr->Grid_mf   , 1, 1); CHKERRQ(ierr);

	// print summary
	PetscPrintf(PETSC_COMM_WORLD, "Passive Tracers output parameters:\n");
	if(pvptr->outpvd) PetscPrintf(PETSC_COMM_WORLD, "   Write Passive tracers pvd file  \n");
	PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");

	// set file name
	sprintf(pvptr->outfile, "%s_passive_tracers", filename);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVPtrWriteTimeStep"
PetscErrorCode PVPtrWriteTimeStep(PVPtr *pvptr, const char *dirName, PetscScalar ttime)
{
	PetscErrorCode ierr;
	PetscFunctionBegin;

	// check activation
	if(pvptr->actx->jr->ctrl.Passive_Tracer == 0) PetscFunctionReturn(0);

	// update .pvd file if necessary
	ierr = UpdatePVDFile(dirName, pvptr->outfile, "pvtu", &pvptr->offset, ttime, pvptr->outpvd); CHKERRQ(ierr);

	// write parallel data .pvtu file
	ierr = PVPtrWritePVTU(pvptr, dirName); CHKERRQ(ierr);

	// write sub-domain data .vtu files
	ierr = PVPtrWriteVTU(pvptr, dirName); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVPtrWriteVTU"
PetscErrorCode PVPtrWriteVTU(PVPtr *pvptr, const char *dirName)
{
	// output markers in .vtu files
	P_Tr       *ptr;
	char       *fname;
	FILE       *fp;
	PetscInt    i, idx, connect;
	uint64_t 	length;
	PetscScalar scal_length;
	float       var,Xp[3];
	PetscInt    var_int;
	PetscScalar *xp,*yp,*zp,*buf;
	size_t      offset = 0;

	PetscErrorCode ierr;
	PetscFunctionBegin;

	// get context
	ptr = pvptr->actx->Ptr;

	// create file name
	asprintf(&fname, "%s/%s_p%1.8lld.vtu", dirName, pvptr->outfile, (LLD)pvptr->actx->iproc);

	// open file
	fp = fopen( fname, "wb" );
	if(fp == NULL) SETERRQ1(PETSC_COMM_SELF, 1,"cannot open file %s", fname);
	free(fname);

	// write header
	WriteXMLHeader(fp, "UnstructuredGrid");

	// initialize connectivity
	connect = ptr->nummark;

	// begin unstructured grid
	fprintf( fp, "\t<UnstructuredGrid>\n" );
	fprintf( fp, "\t\t<Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%lld\">\n",(LLD)ptr->nummark,(LLD)connect );

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
	fprintf( fp, "\t\t\t\t<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%lld\" />\n",(LLD)offset);
	offset += sizeof(uint64_t) + sizeof(float)*(size_t)(ptr->nummark*3);

	fprintf( fp, "\t\t\t</Points>\n");

	// point data - marker phase
	fprintf( fp, "\t\t\t<PointData>\n");
	if(pvptr->Phase)
	{
		fprintf( fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"Phase\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)offset );
		offset += sizeof(uint64_t) + sizeof(int)*(size_t)ptr->nummark;
	}

	if(pvptr->Temperature)
	{
		fprintf( fp, "\t\t\t\t<DataArray type=\"Float32\" Name=\"Temperature %s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",pvptr->actx->jr->scal->lbl_temperature, (LLD)offset);
		offset += sizeof(uint64_t) + sizeof(float)*(size_t)ptr->nummark;
	}
	if(pvptr->Pressure)
	{
		fprintf( fp, "\t\t\t\t<DataArray type=\"Float32\" Name=\"Pressure %s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",pvptr->actx->jr->scal->lbl_stress ,(LLD)offset);
		offset += sizeof(uint64_t) + sizeof(float)*(size_t)ptr->nummark;
	}
	if(pvptr->MeltFraction)
	{
		fprintf( fp, "\t\t\t\t<DataArray type=\"Float32\" Name=\"Mf %s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",pvptr->actx->jr->scal->lbl_unit  ,(LLD)offset);

		offset += sizeof(uint64_t) + sizeof(float)*(size_t)ptr->nummark;
	}
	if(pvptr->Grid_mf)
	{
		fprintf( fp, "\t\t\t\t<DataArray type=\"Float32\" Name=\"Mf_Grid %s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n",pvptr->actx->jr->scal->lbl_unit  ,(LLD)offset);
		offset += sizeof(uint64_t) + sizeof(float)*(size_t)ptr->nummark;
		}

	if(pvptr->ID)
	{
		fprintf( fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"ID\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)offset );
		offset += sizeof(uint64_t) + sizeof(int)*(size_t)ptr->nummark;
	}

	if(pvptr->Active)
	{
		fprintf( fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"Active\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)offset );
		offset += sizeof(uint64_t) + sizeof(int)*(size_t)ptr->nummark;
	}

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
	// scaling length
		scal_length = pvptr->actx->jr->scal->length;

		ierr = VecGetArray(ptr->x, &xp)           ; CHKERRQ(ierr);
		ierr = VecGetArray(ptr->y, &yp)           ; CHKERRQ(ierr);
		ierr = VecGetArray(ptr->z, &zp)           ; CHKERRQ(ierr);



	length = (uint64_t)sizeof(float)*(3*ptr->nummark);
	fwrite( &length,sizeof(uint64_t),1, fp);

	for( i = 0; i < ptr->nummark; i++)
	{
		Xp[0] = (float)(xp[i]*scal_length);
		Xp[1] = (float)(yp[i]*scal_length);
		Xp[2] = (float)(zp[i]*scal_length);
		fwrite( Xp, sizeof(float), (size_t)3, fp );
	}
	ierr = VecRestoreArray(ptr->x, &xp)           ; CHKERRQ(ierr);
	ierr = VecRestoreArray(ptr->y, &yp)           ; CHKERRQ(ierr);
	ierr = VecRestoreArray(ptr->z, &zp)           ; CHKERRQ(ierr);




	// -------------------
	// write field: phases
	// -------------------
	if(pvptr->Phase)
	{
		ierr = VecGetArray(ptr->phase, &buf)           ; CHKERRQ(ierr);

		length = (uint64_t)sizeof(int)*(ptr->nummark);
		fwrite( &length,sizeof(uint64_t),1, fp);

		for( i = 0; i < ptr->nummark; i++)
		{
			var_int = PetscInt(buf[i]);
			fwrite( &var_int, sizeof(int),1, fp );
		}
	// -------------------
		ierr = VecRestoreArray(ptr->phase, &buf)           ; CHKERRQ(ierr);


	}

	if(pvptr->Temperature)
		{
			ierr = VecGetArray(ptr->T, &buf)           ; CHKERRQ(ierr);

			length = (uint64_t)sizeof(float)*(ptr->nummark);
			fwrite( &length,sizeof(uint64_t),1, fp);

			for( i = 0; i < ptr->nummark; i++)
			{
				var = float(buf[i]*pvptr->actx->jr->scal->temperature-pvptr->actx->jr->scal->Tshift);
				fwrite( &var, sizeof(float),1, fp );
			}
		// -------------------
			ierr = VecRestoreArray(ptr->T, &buf)           ; CHKERRQ(ierr);

			// end header

		}


	if(pvptr->Pressure)
		{
			ierr = VecGetArray(ptr->p, &buf)           ; CHKERRQ(ierr);

			length = (uint64_t)sizeof(float)*(ptr->nummark);
			fwrite( &length,sizeof(uint64_t),1, fp);

			for( i = 0; i < ptr->nummark; i++)
			{
				var = float(buf[i]*pvptr->actx->jr->scal->stress);
				fwrite( &var, sizeof(float),1, fp );
			}
		// -------------------
			ierr = VecRestoreArray(ptr->p, &buf)           ; CHKERRQ(ierr);


		}

	if(pvptr->MeltFraction)
		{
			ierr = VecGetArray(ptr->Melt_fr, &buf)           ; CHKERRQ(ierr);

			length = (uint64_t)sizeof(float)*(ptr->nummark);
			fwrite( &length,sizeof(uint64_t),1, fp);

			for( i = 0; i < ptr->nummark; i++)
			{
				var = float(buf[i]);
				fwrite( &var, sizeof(float),1, fp );
			}
		// -------------------
			ierr = VecRestoreArray(ptr->Melt_fr, &buf)           ; CHKERRQ(ierr);


		}

	if(pvptr->Grid_mf)
		{
			ierr = VecGetArray(ptr->Melt_Grid, &buf)           ; CHKERRQ(ierr);

			length = (uint64_t)sizeof(float)*(ptr->nummark);
			fwrite( &length,sizeof(uint64_t),1, fp);

			for( i = 0; i < ptr->nummark; i++)
			{
				var = float(buf[i]);
				fwrite( &var, sizeof(float),1, fp );
			}
		// -------------------
			ierr = VecRestoreArray(ptr->Melt_Grid, &buf)           ; CHKERRQ(ierr);


		}

	if(pvptr->ID)
			{
				ierr = VecGetArray(ptr->ID, &buf)           ; CHKERRQ(ierr);

				length = (uint64_t)sizeof(int)*(ptr->nummark);
				fwrite( &length,sizeof(uint64_t),1, fp);

				for( i = 0; i < ptr->nummark; i++)
				{
					var_int = PetscInt(buf[i]);
					fwrite( &var_int, sizeof(int),1, fp );
				}
			// -------------------
				ierr = VecRestoreArray(ptr->ID, &buf)           ; CHKERRQ(ierr);


			}

	if(pvptr->Active)
			{
				ierr = VecGetArray(ptr->C_advection, &buf)           ; CHKERRQ(ierr);

				length = (uint64_t)sizeof(int)*(ptr->nummark);
				fwrite( &length,sizeof(uint64_t),1, fp);

				for( i = 0; i < ptr->nummark; i++)
				{
					var_int = PetscInt(buf[i]);
					fwrite( &var_int, sizeof(int),1, fp );
				}

				ierr = VecRestoreArray(ptr->C_advection, &buf)           ; CHKERRQ(ierr);

			}


	fprintf( fp,"\n\t</AppendedData>\n");
	fprintf( fp, "</VTKFile>\n");
	// close file
	fclose(fp);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "PVPtrWritePVTU"
PetscErrorCode PVPtrWritePVTU(PVPtr *pvptr, const char *dirName)
{
	// create .pvtu file for marker output
	// load the pvtu file in ParaView and apply a Glyph-spheres filter
	char     *fname;
	FILE     *fp;
	PetscInt i;

	PetscFunctionBegin;

	// only processor 0
	if (!ISRankZero(PETSC_COMM_WORLD)) { PetscFunctionReturn(0); }

	// get context

	// create file name
	asprintf(&fname, "%s/%s.pvtu", dirName, pvptr->outfile);

	// open file
	fp = fopen( fname, "wb" );
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

	fprintf( fp, "\t\t<PPointData>\n");

	if(pvptr->Phase)
	{
		// point data
		fprintf(fp,"\t\t\t<PDataArray type=\"Int32\" Name=\"Phase\" NumberOfComponents=\"1\" format=\"appended\"/>\n");
	}
	if(pvptr->Temperature)
		{
			// point data
			fprintf(fp,"\t\t\t<PDataArray type=\"Float32\" Name=\"Temperature %s\" NumberOfComponents=\"1\" format=\"appended\"/>\n",pvptr->actx->jr->scal->lbl_temperature);
		}
	if(pvptr->Pressure)
		{
			// point data
			fprintf(fp,"\t\t\t<PDataArray type=\"Float32\" Name=\"Pressure %s\" NumberOfComponents=\"1\" format=\"appended\"/>\n",pvptr->actx->jr->scal->lbl_stress);
		}
	if(pvptr->MeltFraction)
		{
			// point data
			fprintf(fp,"\t\t\t<PDataArray type=\"Float32\" Name=\"Mf %s\" NumberOfComponents=\"1\" format=\"appended\"/>\n",pvptr->actx->jr->scal->lbl_unit);
		}

	if(pvptr->Grid_mf)
			{
				// point data
				fprintf(fp,"\t\t\t<PDataArray type=\"Float32\" Name=\"Mf_Grid %s\" NumberOfComponents=\"1\" format=\"appended\"/>\n",pvptr->actx->jr->scal->lbl_unit);
			}
	if(pvptr->ID)
		{
			// point data
			fprintf(fp,"\t\t\t<PDataArray type=\"Int32\" Name=\"ID\" NumberOfComponents=\"1\" format=\"appended\"/>\n");
		}
	if(pvptr->Active)
		{
				// point data
			fprintf(fp,"\t\t\t<PDataArray type=\"Int32\" Name=\"Active\" NumberOfComponents=\"1\" format=\"appended\"/>\n");
		}

	fprintf( fp, "\t\t</PPointData>\n");


	for(i = 0; i < 1; i++){
			fprintf( fp, "\t\t<Piece Source=\"%s_p%1.8lld.vtu\"/>\n",pvptr->outfile,(LLD)i);
		}

	// close the file
	fprintf( fp, "\t</PUnstructuredGrid>\n");
	fprintf( fp, "</VTKFile>\n");

	// close file and free name
	fclose(fp);

	PetscFunctionReturn(0);
}





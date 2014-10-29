/*
LaMEM
Lithosphere and Mantle Evolution Model

A 3D numerical code that can be used to model geodynamical processes such as
mantle-lithosphere interaction.

Originally developed by:

Boris J.P. Kaus (Uni-Mainz, Germany). 2007-
kaus@uni-mainz.de

Further-development:

Dave May (ETH Zurich, Switzerland). 2009-
dave.may@erdw.ethz.ch

Tobias Baumann (Uni-Mainz, Germany). 2011-
baumann@uni-mainz.de

Anton Popov (Uni-Mainz, Germany). 2011-
popov@uni-mainz.de

ParaViewOutput.c, contains the following functions:

WriteOutputFileMatlab			-	Write MATLAB output file for each CPU

$Id:$
$Date:$

This version is compatible with PETSc version 3.3-p2
The code cannot be distributed, used for commercial purposes or handed on,
without the explicit agreement of Boris Kaus.
*/

#include "LaMEM.h"
#include "ParaViewOutput.h"
#include "Material.h"
#include "AVDPhaseViewer.h"
#include "Mesh.h"

#undef __FUNCT__
#define __FUNCT__ "LaMEMView_QuadratureFieldsInit"
PetscErrorCode LaMEMView_QuadratureFieldsInit( LaMEMView_QuadratureFields *view )
{
	PetscErrorCode ierr;
	PetscBool flg,view_all_fields,ascii;

	PetscFunctionBegin;
	ierr = PetscMemzero(view,sizeof(LaMEMView_QuadratureFields));CHKERRQ(ierr);

	view->vtk_ascii = PETSC_TRUE;

	view->phase  = PETSC_TRUE;
	view->mu     = PETSC_TRUE;
	view->rho    = PETSC_TRUE;

	ascii = PETSC_FALSE;
	ierr = PetscOptionsGetBool( PETSC_NULL,"-qfield_ascii", &ascii, &flg );CHKERRQ(ierr);
	if (ascii==PETSC_TRUE) { view->vtk_ascii = PETSC_TRUE; }

	ascii = PETSC_FALSE;
	ierr = PetscOptionsGetBool( PETSC_NULL,"-qfield_binary", &ascii, &flg );CHKERRQ(ierr);
	if (ascii==PETSC_TRUE) { view->vtk_ascii = PETSC_FALSE; }

	ierr = PetscOptionsGetBool( PETSC_NULL,"-qfield_mu",  &view->mu, &flg );CHKERRQ(ierr);
	ierr = PetscOptionsGetBool( PETSC_NULL,"-qfield_rho", &view->rho, &flg );CHKERRQ(ierr);
	ierr = PetscOptionsGetBool( PETSC_NULL,"-qfield_G",   &view->G, &flg );CHKERRQ(ierr);
	ierr = PetscOptionsGetBool( PETSC_NULL,"-qfield_numParticles",     &view->numParticles, &flg );CHKERRQ(ierr);
	ierr = PetscOptionsGetBool( PETSC_NULL,"-qfield_pressure",         &view->pressure, &flg );CHKERRQ(ierr);
	ierr = PetscOptionsGetBool( PETSC_NULL,"-qfield_temperature",         &view->temperature, &flg );CHKERRQ(ierr);

	ierr = PetscOptionsGetBool( PETSC_NULL,"-qfield_devStrainrateInv", &view->devStrainrateInv, &flg );CHKERRQ(ierr);
	ierr = PetscOptionsGetBool( PETSC_NULL,"-qfield_devStressInv",     &view->devStressInv, &flg );CHKERRQ(ierr);

	ierr = PetscOptionsGetBool( PETSC_NULL,"-qfield_stress",      &view->stress, &flg );CHKERRQ(ierr);
	ierr = PetscOptionsGetBool( PETSC_NULL,"-qfield_strainrate",  &view->strainrate, &flg );CHKERRQ(ierr);

	ierr = PetscOptionsGetBool( PETSC_NULL,"-qfield_strain",        &view->strain, &flg );CHKERRQ(ierr);
	ierr = PetscOptionsGetBool( PETSC_NULL,"-qfield_plasticStrain", &view->plasticStrain, &flg );CHKERRQ(ierr);

	ierr = PetscOptionsGetBool( PETSC_NULL,"-qfield_phase", &view->phase, &flg );CHKERRQ(ierr);

	view_all_fields = PETSC_FALSE;
	ierr = PetscOptionsGetBool( PETSC_NULL,"-qfield_all", &view_all_fields, &flg );CHKERRQ(ierr);
	if (view_all_fields) {
		view->mu           = PETSC_TRUE;
		view->rho          = PETSC_TRUE;
		view->G            = PETSC_TRUE;
		view->numParticles = PETSC_TRUE;
		view->pressure     = PETSC_TRUE;
		view->temperature  = PETSC_TRUE;
		view->devStrainrateInv = PETSC_TRUE;
		view->devStressInv     = PETSC_TRUE;
		view->stress        = PETSC_TRUE;
		view->strainrate    = PETSC_TRUE;
		view->strain        = PETSC_TRUE;
		view->plasticStrain = PETSC_TRUE;
		view->phase         = PETSC_TRUE;
	}

	PetscFunctionReturn(0);
}

#define LAMEMVTU_LOAD_AND_WRITE(data)\
length = (PetscInt)sizeof(PetscScalar)*NGP;\
fwrite( &length,sizeof(PetscInt),1,vtk_fp);\
for (iz=0; iz<zm; iz++){for (iy=0; iy<ym; iy++){for (ix=0; ix<xm; ix++){\
  LaMEMSetMaterialDataMemoryFromArray( &material_data, ix,iy,iz, ngp_vel, materials_array );\
  for (intp=0; intp<ngp_vel; intp++){\
    scalar_quadrature_field[intp] = data[intp];\
   }\
	 fwrite( scalar_quadrature_field, sizeof(PetscScalar), (size_t)ngp_vel, vtk_fp );\
}}}\

#undef __FUNCT__
#define __FUNCT__ "LaMEMViewQuadraturePoints_3DPVTU_appended"
PetscErrorCode LaMEMViewQuadraturePoints_3DPVTU_appended( LaMEMView_QuadratureFields *view, UserContext *user, PetscInt xm, PetscInt ym, PetscInt zm, DM DA_Materials_fine, Vec Materials_fine, PetscInt ngp_vel, const char file_prefix[] )
{
	char *vtk_filename;
	PetscMPIInt rank;
	FILE *vtk_fp;
	PetscInt k,NGP;
	PetscInt memory_offset;
	PetscInt length;
	PetscScalar ***materials_array;
	MaterialsElementDynamic material_data;
	PetscLogDouble t0,t1;
	PetscErrorCode ierr;

	PetscFunctionBegin;
	ierr = PetscTime(&t0);CHKERRQ(ierr);


	/* create file name */
	MPI_Comm_rank( PETSC_COMM_WORLD, &rank );
	asprintf( &vtk_filename, "subdomain-%s-p%1.6d.vtu", file_prefix, rank );

	/* open file and write header */
	vtk_fp = fopen( vtk_filename, "w" );
	if( vtk_fp == NULL ) {
		printf("ERROR(%s): Cannot open file = %s \n", __FUNCT__,vtk_filename );
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"File opening error");
	}

	fprintf( vtk_fp, "<?xml version=\"1.0\"?>\n");
#ifdef PETSC_WORDS_BIGENDIAN
	fprintf( vtk_fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	fprintf( vtk_fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif

	fprintf( vtk_fp, "\t<UnstructuredGrid>\n" );

	NGP = xm*ym*zm * ngp_vel;
	fprintf( vtk_fp, "\t\t<Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%lld\">\n",(LLD)NGP,(LLD)NGP );


	memory_offset = 0;

	fprintf( vtk_fp, "\n");
	fprintf( vtk_fp, "\t\t\t<Cells>\n");

	// Distinguish 32-bit and 64-bit integers

#if defined(PETSC_USE_64BIT_INDICES)
	// connectivity
	fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Int64\" Name=\"connectivity\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)memory_offset);
	memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscInt)*NGP;
	// offsets
	fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Int64\" Name=\"offsets\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)memory_offset);
	memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscInt)*NGP;
	// types
	fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Int64\" Name=\"types\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)memory_offset);
	memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscInt)*NGP;
#else
	// connectivity
	fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)memory_offset);
	memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscInt)*NGP;
	// offsets
	fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)memory_offset);
	memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscInt)*NGP;
	// types
	fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"appended\" offset=\"%lld\"/>\n",(LLD)memory_offset);
	memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscInt)*NGP;
#endif

	fprintf( vtk_fp, "\t\t\t</Cells>\n");

	fprintf( vtk_fp, "\n");
	fprintf( vtk_fp, "\t\t\t<CellData>\n");
	fprintf( vtk_fp, "\t\t\t</CellData>\n");

	fprintf( vtk_fp, "\n");
	fprintf( vtk_fp, "\t\t\t<Points>\n");

	/* coordinates */
	fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%lld\" />\n",(LLD)memory_offset);
	memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscScalar)*NGP*3;

	fprintf( vtk_fp, "\t\t\t</Points>\n");


	fprintf( vtk_fp, "\n");
	fprintf( vtk_fp, "\t\t\t<PointData Scalars=\" ");
	fprintf( vtk_fp, "\">\n");

	/* fields */
	if (view->mu) {
		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"mu\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)memory_offset );
		memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscScalar)*NGP;
	}
	if (view->rho) {
		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"rho\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)memory_offset );
		memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscScalar)*NGP;
	}
	if (view->phase) {
		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"phase\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)memory_offset );
		memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscScalar)*NGP;
	}
	if (view->G) {
		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"G\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)memory_offset );
		memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscScalar)*NGP;
	}
	if (view->numParticles) {
		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"NumParticles\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)memory_offset );
		memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscScalar)*NGP;
	}
	if (view->pressure) {
		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"Pressure\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)memory_offset );
		memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscScalar)*NGP;
	}
	if (view->temperature) {
		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"Temperature\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)memory_offset );
		memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscScalar)*NGP;
	}
	if (view->devStrainrateInv) {
		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"DevStrainrateInv\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)memory_offset );
		memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscScalar)*NGP;
	}
	if (view->devStressInv) {
		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"DevStressInv\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)memory_offset );
		memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscScalar)*NGP;
	}
	if (view->strain) {
		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"Strain\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)memory_offset );
		memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscScalar)*NGP;
	}
	if (view->plasticStrain) {
		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"PlasticStrain\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)memory_offset );
		memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscScalar)*NGP;
	}

	if (view->n) {
	}
	if (view->C) {
	}
	if (view->phi) {
	}
	if (view->k) {
	}
	if (view->Cp) {
	}
	if (view->Q) {
	}
	if (view->alpha) {
	}
	if (view->FK) {
	}

	if ((view->stress) && (1==0)) {
		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"Txx\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)memory_offset );
		memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscScalar)*NGP;

		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"Tyy\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)memory_offset );
		memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscScalar)*NGP;

		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"Tzz\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)memory_offset );
		memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscScalar)*NGP;

		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"Txy\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)memory_offset );
		memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscScalar)*NGP;

		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"Txz\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)memory_offset );
		memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscScalar)*NGP;

		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"Tyz\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)memory_offset );
		memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscScalar)*NGP;
	}
	if ((view->strainrate) && (1==0)) {
		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"Exx\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)memory_offset );
		memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscScalar)*NGP;

		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"Eyy\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)memory_offset );
		memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscScalar)*NGP;

		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"Ezz\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)memory_offset );
		memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscScalar)*NGP;

		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"Exy\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)memory_offset );
		memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscScalar)*NGP;

		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"Exz\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)memory_offset );
		memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscScalar)*NGP;

		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"Eyz\" format=\"appended\" offset=\"%lld\"/>\n", (LLD)memory_offset );
		memory_offset = memory_offset + (PetscInt)sizeof(PetscInt) + (PetscInt)sizeof(PetscScalar)*NGP;
	}



	/* end fields */

	fprintf( vtk_fp, "\t\t\t</PointData>\n");
	fprintf( vtk_fp, "\n");

	fprintf( vtk_fp, "\t\t</Piece>\n");
	fprintf( vtk_fp, "\t</UnstructuredGrid>\n");

	fprintf( vtk_fp,"\t<AppendedData encoding=\"raw\">\n");
	fprintf( vtk_fp,"_");

	ierr = DMDAVecGetArray(DA_Materials_fine,Materials_fine, &materials_array); CHKERRQ(ierr);
	{
		PetscInt idx;
		PetscScalar xp[3];
		const PetscInt max_quadrature_points = 100;
		PetscScalar scalar_quadrature_field[max_quadrature_points];
		PetscScalar *scalar_quadrature_field_ptr;

		PetscInt ix,iy,iz,intp;
		PetscLogDouble ta,tb;

		if (ngp_vel>max_quadrature_points) {
			SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Using more than 100 quadrature points.. updated code to use dynamic allocation");
		}

		////////////////////////////////////////////////////////
		/* write connectivity */
		length = (PetscInt)sizeof(PetscInt)*NGP;
		fwrite( &length,sizeof(PetscInt),1,vtk_fp);
		for(k=0;k<NGP;k++) {
			idx = k;
			fwrite( &idx, sizeof(PetscInt),1, vtk_fp );
		}
		////////////////////////////////////////////////////////
		/* write offsets */
		length = (PetscInt)sizeof(PetscInt)*NGP;
		fwrite( &length,sizeof(PetscInt),1,vtk_fp);
		for(k=0;k<NGP;k++) {
			idx = k+1;
			fwrite( &idx, sizeof(PetscInt),1, vtk_fp );
		}
		////////////////////////////////////////////////////////
		/* write types */
		length = (PetscInt)sizeof(PetscInt)*NGP;
		fwrite( &length,sizeof(PetscInt),1,vtk_fp);
		for(k=0;k<NGP;k++) {
			idx = 1;
			fwrite( &idx, sizeof(PetscInt),1, vtk_fp );
		}
		////////////////////////////////////////////////////////
		/* write coordinates */
		length = (PetscInt)sizeof(PetscScalar)*NGP*3;
		fwrite( &length,sizeof(PetscInt),1,vtk_fp);
		ierr = PetscTime(&ta);CHKERRQ(ierr);
		for (iz=0; iz<zm; iz++){
			for (iy=0; iy<ym; iy++){
				for (ix=0; ix<xm; ix++){
					LaMEMSetMaterialDataMemoryFromArray( &material_data, ix,iy,iz, ngp_vel, materials_array );

					for (intp=0; intp<ngp_vel; intp++){
						xp[0] = material_data.Coord[0][intp]*user->Characteristic.Length;
						xp[1] = material_data.Coord[1][intp]*user->Characteristic.Length;
						xp[2] = material_data.Coord[2][intp]*user->Characteristic.Length;

						fwrite( xp, sizeof(PetscScalar),3, vtk_fp );
					}
		}}}
		ierr = PetscTime(&tb);CHKERRQ(ierr);

		////////////////////////////////////////////////////////
		/* write field: mu */
		if (view->mu) {
			length = (PetscInt)sizeof(PetscScalar)*NGP;
			fwrite( &length,sizeof(PetscInt),1,vtk_fp);
			for (iz=0; iz<zm; iz++){
				for (iy=0; iy<ym; iy++){
					for (ix=0; ix<xm; ix++){
						LaMEMSetMaterialDataMemoryFromArray( &material_data, ix,iy,iz, ngp_vel, materials_array );
						// load //
						for (intp=0; intp<ngp_vel; intp++){
							scalar_quadrature_field[intp] = material_data.Viscosity[intp]*user->Characteristic.Viscosity;
						}
						// write //
						fwrite( scalar_quadrature_field, sizeof(PetscScalar), (size_t)ngp_vel, vtk_fp );
			}}}
			ierr = PetscTime(&tb);CHKERRQ(ierr);
		}
		//
		////////////////////////////////////////////////////////
		/* write field: rho */
		if (view->rho) {
			length = (PetscInt)sizeof(PetscScalar)*NGP;
			fwrite( &length,sizeof(PetscInt),1,vtk_fp);
			for (iz=0; iz<zm; iz++){
				for (iy=0; iy<ym; iy++){
					for (ix=0; ix<xm; ix++){
						LaMEMSetMaterialDataMemoryFromArray( &material_data, ix,iy,iz, ngp_vel, materials_array );
						// load //
						for (intp=0; intp<ngp_vel; intp++){
							scalar_quadrature_field[intp] = material_data.Density[intp]*user->Characteristic.Density;
						}
						// write //
						fwrite( scalar_quadrature_field, sizeof(PetscScalar), (size_t)ngp_vel, vtk_fp );
					}}}
			ierr = PetscTime(&tb);CHKERRQ(ierr);
		}
		//

		////////////////////////////////////////////////////////
		/* write field: G */
		if (view->G) {
			length = (PetscInt)sizeof(PetscScalar)*NGP;
			fwrite( &length,sizeof(PetscInt),1,vtk_fp);
			for (iz=0; iz<zm; iz++){
				for (iy=0; iy<ym; iy++){
					for (ix=0; ix<xm; ix++){
						LaMEMSetMaterialDataMemoryFromArray( &material_data, ix,iy,iz, ngp_vel, materials_array );
						// load //
						for (intp=0; intp<ngp_vel; intp++){
							scalar_quadrature_field[intp] = material_data.ElasticShearModule[intp]*user->Characteristic.Stress;
						}
						// write //
						fwrite( scalar_quadrature_field, sizeof(PetscScalar), (size_t)ngp_vel, vtk_fp );
					}}}
			ierr = PetscTime(&tb);CHKERRQ(ierr);
		}

		//
		if (view->phase) {
			scalar_quadrature_field_ptr = &material_data.Phases[0];
			LAMEMVTU_LOAD_AND_WRITE(scalar_quadrature_field_ptr);
		}
		if (view->numParticles) {
			scalar_quadrature_field_ptr = &material_data.NumParticles[0];
			LAMEMVTU_LOAD_AND_WRITE(scalar_quadrature_field_ptr);
		}
		if (view->pressure) {
			scalar_quadrature_field_ptr = &material_data.Pressure[0];
			LAMEMVTU_LOAD_AND_WRITE(scalar_quadrature_field_ptr);
		}
		if (view->temperature) {
			scalar_quadrature_field_ptr = &material_data.Temperature[0];
			LAMEMVTU_LOAD_AND_WRITE(scalar_quadrature_field_ptr);
		}
		if (view->devStrainrateInv) {
			scalar_quadrature_field_ptr = &material_data.SecondInvariantDevStrainrate[0];
			LAMEMVTU_LOAD_AND_WRITE(scalar_quadrature_field_ptr);
		}
		if (view->devStressInv) {
			scalar_quadrature_field_ptr = &material_data.SecondInvariantDevStress[0];
			LAMEMVTU_LOAD_AND_WRITE(scalar_quadrature_field_ptr);
		}
		if (view->strain) {
			scalar_quadrature_field_ptr = &material_data.Strain[0];
			LAMEMVTU_LOAD_AND_WRITE(scalar_quadrature_field_ptr);
		}
		if (view->plasticStrain) {
			scalar_quadrature_field_ptr = &material_data.PlasticStrain[0];
			LAMEMVTU_LOAD_AND_WRITE(scalar_quadrature_field_ptr);
		}


		if (view->n) {
		}
		if (view->C) {
		}
		if (view->phi) {
		}
		if (view->k) {
		}
		if (view->Cp) {
		}
		if (view->Q) {
		}
		if (view->alpha) {
		}
		if (view->FK) {
		}

		if (view->stress) {
			scalar_quadrature_field_ptr = &material_data.DevStress[0][0];
			LAMEMVTU_LOAD_AND_WRITE(scalar_quadrature_field_ptr);

			scalar_quadrature_field_ptr = &material_data.DevStress[1][0];
			LAMEMVTU_LOAD_AND_WRITE(scalar_quadrature_field_ptr);

			scalar_quadrature_field_ptr = &material_data.DevStress[2][0];
			LAMEMVTU_LOAD_AND_WRITE(scalar_quadrature_field_ptr);

			scalar_quadrature_field_ptr = &material_data.DevStress[3][0];
			LAMEMVTU_LOAD_AND_WRITE(scalar_quadrature_field_ptr);

			scalar_quadrature_field_ptr = &material_data.DevStress[4][0];
			LAMEMVTU_LOAD_AND_WRITE(scalar_quadrature_field_ptr);

			scalar_quadrature_field_ptr = &material_data.DevStress[5][0];
			LAMEMVTU_LOAD_AND_WRITE(scalar_quadrature_field_ptr);
		}
		if (view->strainrate) {
			scalar_quadrature_field_ptr = &material_data.DevStrainrate[0][0];
			LAMEMVTU_LOAD_AND_WRITE(scalar_quadrature_field_ptr);

			scalar_quadrature_field_ptr = &material_data.DevStrainrate[1][0];
			LAMEMVTU_LOAD_AND_WRITE(scalar_quadrature_field_ptr);

			scalar_quadrature_field_ptr = &material_data.DevStrainrate[2][0];
			LAMEMVTU_LOAD_AND_WRITE(scalar_quadrature_field_ptr);

			scalar_quadrature_field_ptr = &material_data.DevStrainrate[3][0];
			LAMEMVTU_LOAD_AND_WRITE(scalar_quadrature_field_ptr);

			scalar_quadrature_field_ptr = &material_data.DevStrainrate[4][0];
			LAMEMVTU_LOAD_AND_WRITE(scalar_quadrature_field_ptr);

			scalar_quadrature_field_ptr = &material_data.DevStrainrate[5][0];
			LAMEMVTU_LOAD_AND_WRITE(scalar_quadrature_field_ptr);
		}

	}
	ierr = DMDAVecRestoreArray(DA_Materials_fine,Materials_fine, &materials_array); CHKERRQ(ierr);

	fprintf( vtk_fp,"\n\t</AppendedData>\n");

	fprintf( vtk_fp, "</VTKFile>\n");


	if( vtk_fp!= NULL ) {
		fclose( vtk_fp );
		vtk_fp = NULL;
	}

	free(vtk_filename);

	ierr = PetscTime(&t1);CHKERRQ(ierr);
#ifdef PROFILE_TIMING
	PetscPrintf(PETSC_COMM_WORLD,"VTKWriter(%s): Time %1.4e sec\n",__FUNCT__,t1-t0);
#endif
	PetscFunctionReturn(0);
}


/* parallel vts file for quadrature point information */
#undef __FUNCT__
#define __FUNCT__ "LaMEMViewQuadraturePoints_Meta3DPVTS"
PetscErrorCode LaMEMViewQuadraturePoints_Meta3DPVTS( LaMEMView_QuadratureFields *view, UserContext *user, const char file_prefix[], const char local_file_prefix[], DM da, const char DirectoryName[] )
{
	PetscMPIInt nproc,rank;
	char *vtk_filename;
	FILE *vtk_fp;
	PetscInt M, N, P, swidth;

	PetscFunctionBegin;

	/* only master generates this file */
	MPI_Comm_size( PETSC_COMM_WORLD, &nproc );
	MPI_Comm_rank( PETSC_COMM_WORLD, &rank );

	if( rank != 0 ) { PetscFunctionReturn(0); }

	/* create file name */
	asprintf( &vtk_filename, "./%s/%s.pvts", DirectoryName,file_prefix );
	vtk_fp = fopen( vtk_filename, "w" );
	if( vtk_fp == NULL ) {
		printf("ERROR(%s): Cannot open file = %s on rank %lld \n", __FUNCT__, vtk_filename, (LLD)rank );
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"File opening error");
	}
	free( vtk_filename );

	/* (VTK) generate pvts header */
	fprintf( vtk_fp, "<?xml version=\"1.0\"?>\n");

#ifdef PETSC_WORDS_BIGENDIAN
	fprintf( vtk_fp, "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	fprintf( vtk_fp, "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif

	/* define size of the nodal mesh based on the cell DM */
	//fprintf( vtk_fp, "\t<PStructuredGrid GhostLevel=\"0\">\n" ); /* note overlap = 0 */
	DMDAGetInfo( da, 0, &M,&N,&P, 0,0,0, 0,&swidth,0,0,0,0 );
	fprintf( vtk_fp, "\t<PStructuredGrid GhostLevel=\"%lld\" WholeExtent=\"%lld %lld %lld %lld %lld %lld\">\n", (LLD)swidth, 0LL,(LLD)M-1LL, 0LL,(LLD)N-1LL, 0LL,(LLD)P-1LL ); /* note overlap = 1 for Q1 */

	/* DUMP THE CELL REFERENCES */
	fprintf( vtk_fp, "\t\t<PCellData>\n");
	fprintf( vtk_fp, "\t\t</PCellData>\n");


	fprintf( vtk_fp, "\t\t<PPoints>\n");
	fprintf( vtk_fp, "\t\t\t<PDataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\"/>\n");
	fprintf( vtk_fp, "\t\t</PPoints>\n");

	fprintf( vtk_fp, "\t\t<PPointData>\n");

	if (view->mu) {
		if (user->DimensionalUnits==1){
			fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"mu[Pas]\" NumberOfComponents=\"1\"/>\n");
		}
		else{
			fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"mu[]\" NumberOfComponents=\"1\"/>\n");
		}
	}
	if (view->phase) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"phase\" NumberOfComponents=\"1\"/>\n");
	}
	if (view->rho) {
		if (user->DimensionalUnits==1){
			fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"rho[kg/m3]\" NumberOfComponents=\"1\"/>\n");
		}
		else{
			fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"rho[ ]\" NumberOfComponents=\"1\"/>\n");
		}
	}
	if (view->G) {
		if (user->DimensionalUnits==1){
			fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"G [Pa]\" NumberOfComponents=\"1\"/>\n");
		}
		else{
			fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"G [ ]\" NumberOfComponents=\"1\"/>\n");
		}
	}


	if (view->numParticles) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"NumParticles\" NumberOfComponents=\"1\"/>\n");
	}
	if (view->pressure) {
		if (user->DimensionalUnits==1){
			fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Pressure [MPa]\" NumberOfComponents=\"1\"/>\n");
		}
		else{
			fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Pressure [ ]\" NumberOfComponents=\"1\"/>\n");
		}

	}
	if (view->temperature) {
		if (user->DimensionalUnits==1){
			fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Temperature [K]\" NumberOfComponents=\"1\"/>\n");
		}
		else{
			fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Temperature [ ]\" NumberOfComponents=\"1\"/>\n");
		}

	}
	if (view->devStrainrateInv) {
		if (user->DimensionalUnits==1){
			fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"DevStrainrateInv [1/s]\" NumberOfComponents=\"1\"/>\n");
		}
		else{
			fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"DevStrainrateInv [ ]\" NumberOfComponents=\"1\"/>\n");
		}
	}
	if (view->devStressInv) {
		if (user->DimensionalUnits==1){
			fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"DevStressInv [MPa]\" NumberOfComponents=\"1\"/>\n");
		}
		else{
			fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"DevStressInv [ ]\" NumberOfComponents=\"1\"/>\n");
		}
	}

	if (view->strain) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Strain [ ]\" NumberOfComponents=\"1\"/>\n");
	}
	if (view->plasticStrain) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"PlasticStrain [ ]\" NumberOfComponents=\"1\"/>\n");
	}


	if ((view->stress) && (1==0)) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Txx\" NumberOfComponents=\"1\"/>\n");
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Tyy\" NumberOfComponents=\"1\"/>\n");
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Tzz\" NumberOfComponents=\"1\"/>\n");
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Txy\" NumberOfComponents=\"1\"/>\n");
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Txz\" NumberOfComponents=\"1\"/>\n");
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Tyz\" NumberOfComponents=\"1\"/>\n");
	}
	if ((view->strainrate) && (1==0)) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Exx\" NumberOfComponents=\"1\"/>\n");
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Eyy\" NumberOfComponents=\"1\"/>\n");
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Ezz\" NumberOfComponents=\"1\"/>\n");
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Exy\" NumberOfComponents=\"1\"/>\n");
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Exz\" NumberOfComponents=\"1\"/>\n");
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Eyz\" NumberOfComponents=\"1\"/>\n");
	}

	/* These are not supported */
	/*
	if (view->n) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"n\" NumberOfComponents=\"1\"/>\n");
	}
	if (view->C) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"C\" NumberOfComponents=\"1\"/>\n");
	}
	if (view->phi) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"phi\" NumberOfComponents=\"1\"/>\n");
	}
	if (view->k) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"k\" NumberOfComponents=\"1\"/>\n");
	}
	if (view->Cp) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Cp\" NumberOfComponents=\"1\"/>\n");
	}
	if (view->Q) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Q\" NumberOfComponents=\"1\"/>\n");
	}
	if (view->alpha) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Alpha\" NumberOfComponents=\"1\"/>\n");
	}
	if (view->FK) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"FK\" NumberOfComponents=\"1\"/>\n");
	}
	*/

	fprintf( vtk_fp, "\t\t</PPointData>\n");

	/* write out the parallel information */
	DMViewVTK_write_PieceExtend(vtk_fp,2,da,local_file_prefix);


	/* close the file */
	fprintf( vtk_fp, "\t</PStructuredGrid>\n");
	fprintf( vtk_fp, "</VTKFile>\n");

	if(vtk_fp!=NULL){
		fclose( vtk_fp );
		vtk_fp = NULL;
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LaMEMViewQuadraturePoints_Meta3DPVTU"
PetscErrorCode LaMEMViewQuadraturePoints_Meta3DPVTU( LaMEMView_QuadratureFields *view, const char file_prefix[], const char local_file_prefix[])
{
	PetscMPIInt nproc,rank;
	char *vtk_filename, *name;
	FILE *vtk_fp;
	PetscInt i;

	PetscFunctionBegin;

	/* only master generates this file */
	MPI_Comm_size( PETSC_COMM_WORLD, &nproc );
	MPI_Comm_rank( PETSC_COMM_WORLD, &rank );

	if( rank != 0 ) { PetscFunctionReturn(0); }

	/* create file name */
	asprintf( &vtk_filename, "%s.pvtu", file_prefix );
	vtk_fp = fopen( vtk_filename, "w" );
	if( vtk_fp == NULL ) {
		printf("ERROR(%s): Cannot open file = %s on rank %lld \n", __FUNCT__, vtk_filename, (LLD)rank );
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"File opening error");
	}
	free( vtk_filename );

	/* (VTK) generate pvts header */
	fprintf( vtk_fp, "<?xml version=\"1.0\"?>\n");

#ifdef PETSC_WORDS_BIGENDIAN
	fprintf( vtk_fp, "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	fprintf( vtk_fp, "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif

	/* define size of the nodal mesh based on the cell DM */
	fprintf( vtk_fp, "\t<PUnstructuredGrid GhostLevel=\"0\">\n" ); /* note overlap = 0 */

	/* DUMP THE CELL REFERENCES */
	fprintf( vtk_fp, "\t\t<PCellData>\n");
	fprintf( vtk_fp, "\t\t</PCellData>\n");


	fprintf( vtk_fp, "\t\t<PPoints>\n");
	fprintf( vtk_fp, "\t\t\t<PDataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\"/>\n");
	fprintf( vtk_fp, "\t\t</PPoints>\n");

	fprintf( vtk_fp, "\t\t<PPointData>\n");

	if (view->mu) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"mu\" NumberOfComponents=\"1\"/>\n");
	}
	if (view->rho) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"rho\" NumberOfComponents=\"1\"/>\n");
	}
	if (view->G) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"G\" NumberOfComponents=\"1\"/>\n");
	}

	if (view->phase) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"phase\" NumberOfComponents=\"1\"/>\n");
	}
	if (view->numParticles) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"NumParticles\" NumberOfComponents=\"1\"/>\n");
	}
	if (view->pressure) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Pressure\" NumberOfComponents=\"1\"/>\n");
	}
	if (view->temperature) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Temperature\" NumberOfComponents=\"1\"/>\n");
	}
	if (view->devStrainrateInv) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"DevStrainrateInv\" NumberOfComponents=\"1\"/>\n");
	}
	if (view->devStressInv) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"DevStressInv\" NumberOfComponents=\"1\"/>\n");
	}

	if (view->strain) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Strain\" NumberOfComponents=\"1\"/>\n");
	}
	if (view->plasticStrain) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"PlasticStrain\" NumberOfComponents=\"1\"/>\n");
	}

	if (view->stress) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Txx\" NumberOfComponents=\"1\"/>\n");
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Tyy\" NumberOfComponents=\"1\"/>\n");
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Tzz\" NumberOfComponents=\"1\"/>\n");
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Txy\" NumberOfComponents=\"1\"/>\n");
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Txz\" NumberOfComponents=\"1\"/>\n");
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Tyz\" NumberOfComponents=\"1\"/>\n");
	}
	if (view->strainrate) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Exx\" NumberOfComponents=\"1\"/>\n");
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Eyy\" NumberOfComponents=\"1\"/>\n");
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Ezz\" NumberOfComponents=\"1\"/>\n");
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Exy\" NumberOfComponents=\"1\"/>\n");
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Exz\" NumberOfComponents=\"1\"/>\n");
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Eyz\" NumberOfComponents=\"1\"/>\n");
	}

	/* These are not supported */
	/*
	if (view->n) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"n\" NumberOfComponents=\"1\"/>\n");
	}
	if (view->C) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"C\" NumberOfComponents=\"1\"/>\n");
	}
	if (view->phi) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"phi\" NumberOfComponents=\"1\"/>\n");
	}
	if (view->k) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"k\" NumberOfComponents=\"1\"/>\n");
	}
	if (view->Cp) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Cp\" NumberOfComponents=\"1\"/>\n");
	}
	if (view->Q) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Q\" NumberOfComponents=\"1\"/>\n");
	}
	if (view->alpha) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Alpha\" NumberOfComponents=\"1\"/>\n");
	}
	if (view->FK) {
		fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"FK\" NumberOfComponents=\"1\"/>\n");
	}
	*/

	fprintf( vtk_fp, "\t\t</PPointData>\n");
	for(i=0;i<nproc;i++){
		asprintf( &name, "subdomain-%s-p%1.6lld.vtu", local_file_prefix, (LLD)i );

		fprintf( vtk_fp, "\t\t<Piece Source=\"%s\"/>\n",name);
		free(name);
	}






	/* close the file */
	fprintf( vtk_fp, "\t</PUnstructuredGrid>\n");
	fprintf( vtk_fp, "</VTKFile>\n");

	if(vtk_fp!=NULL){
		fclose( vtk_fp );
		vtk_fp = NULL;
	}
	PetscFunctionReturn(0);
}

//==========================================================================================================
/* LaMEMViewQuadraturePoints_3DPVTU
 * Prints the quadrature points as points (unstructured VTU file format)
 */
#undef __FUNCT__
#define __FUNCT__ "LaMEMViewQuadraturePoints_3DPVTU"
PetscErrorCode LaMEMViewQuadraturePoints_3DPVTU( LaMEMView_QuadratureFields *view,  UserContext *user, DM DA_Processors_fine, DM DA_Materials_fine, Vec Material_fine, PetscInt ngp_vel, const char NAME[] )
{
	char *vts_filename;
	char *pvts_filename;
	PetscInt xm,ym,zm;
	PetscErrorCode ierr;

	PetscFunctionBegin;
	if(view->vtk_ascii) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Quadrature point fields can only be printed in binary");
	}
	else {
		ierr = DMDAGetCorners(DA_Processors_fine,  PETSC_NULL,PETSC_NULL,PETSC_NULL, &xm, &ym, &zm); CHKERRQ(ierr);


		/* This prints the fields as unstructured point clouds */
		asprintf( &vts_filename, "%s-points", NAME );
		ierr = LaMEMViewQuadraturePoints_3DPVTU_appended( view, user, xm,ym,zm, DA_Materials_fine, Material_fine, ngp_vel, vts_filename );CHKERRQ(ierr);

		asprintf( &pvts_filename, "%s-points", NAME );
		ierr = LaMEMViewQuadraturePoints_Meta3DPVTU( view, pvts_filename, vts_filename );CHKERRQ(ierr);

		free(pvts_filename);
		free(vts_filename);


	}

	PetscFunctionReturn(0);
}


//==========================================================================================================
/* LaMEMViewQuadraturePoints_3DPVTS
 * Prints the quadrature point information as meshes (structured VTS file format)
 */
#undef __FUNCT__
#define __FUNCT__ "LaMEMViewQuadraturePoints_3DPVTS"
PetscErrorCode LaMEMViewQuadraturePoints_3DPVTS( LaMEMView_QuadratureFields *view,  UserContext *user, DM DA_Processors_fine, DM DA_Materials_fine, Vec Material_fine, PetscInt ngp_vel_1D, PetscInt itime, const char DirectoryName[] )
{
	char *vts_filename, *pfield_name, *mfield_name;
	char *pvts_filename, *NAME;
	PetscInt xm,ym,zm;
	PetscErrorCode ierr;


	PetscFunctionBegin;


	asprintf(&pfield_name, "%s_quadrature_timeseries.pvd",user->OutputFile);
	if (itime == 0) {
		/* create pvd file for temperature */
		ierr = ParaviewPVDOpen(pfield_name);CHKERRQ(ierr);
	}


	asprintf(&NAME, "%s_quadrature_step%1.6lld", user->OutputFile, (LLD)itime );		// Quadrature points

	if(view->vtk_ascii) {
		ierr = DMDAGetCorners(DA_Processors_fine,  PETSC_NULL,PETSC_NULL,PETSC_NULL, &xm, &ym, &zm); CHKERRQ(ierr);


		asprintf( &vts_filename, "%s-mesh", NAME );
		ierr = LaMEM_DMView_3DVTK_StructuredGrid_QuadPoints(view, DA_Materials_fine, DA_Processors_fine, user->DA_Quadrature, Material_fine,vts_filename, ngp_vel_1D, user, DirectoryName);				CHKERRQ(ierr);

		asprintf( &pvts_filename, "%s-mesh", NAME );
		ierr = LaMEMViewQuadraturePoints_Meta3DPVTS( view, user, pvts_filename, vts_filename, user->DA_Quadrature, DirectoryName  );CHKERRQ(ierr);

		free(pvts_filename);
		free(vts_filename);


	}
	else{
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Quadrature mesh field files can only be printed in ASCII");
	}
	PetscPrintf(PETSC_COMM_WORLD,"# Saved Paraview VTS point-wise quadrature output in file %s-mesh.pvts \n",NAME); CHKERRQ(ierr);


	/* add datafile to the time series */
	asprintf(&mfield_name, "%s-mesh.pvts", NAME );
	if (user->DimensionalUnits==1){
		ierr = ParaviewPVDAppend(pfield_name,user->time * user->Characteristic.Time/user->Characteristic.SecYear/1e6, mfield_name, DirectoryName );CHKERRQ(ierr);		// in Myrs
	}
	else{
		ierr = ParaviewPVDAppend(pfield_name,user->time , mfield_name, DirectoryName );CHKERRQ(ierr);
	}



	free(NAME);
	free(mfield_name);
	free(pfield_name);



	PetscFunctionReturn(0);
}

//==========================================================================================================
// Set the quadrature point coordinates to the Quadrature DM
#undef __FUNCT__
#define __FUNCT__ "LaMEM_SetQuadraturePointCoords_to_QuadratureDA"
PetscErrorCode LaMEM_SetQuadraturePointCoords_to_QuadratureDA(DM DA_Quadrature, PetscInt ngp_vel_1D, DM DA_Materials_fine, Vec Materials_fine )
{
	DM 						cda;
	Vec 					gc;
	DMDACoor3d 				***coors;
	PetscInt 				nx,ny,nz, i,j,k;
	PetscInt				xs,ys,zs,xm,ym,zm,ix,iy,iz, intp, intpx, intpy, intpz;
	PetscErrorCode 			ierr;
	MaterialsElementDynamic material_data;
	PetscScalar 			***materials_array;


	PetscFunctionBegin;

	/* retrieve data array from quad points */
	ierr = DMDAVecGetArray(DA_Materials_fine	,   Materials_fine, &materials_array); CHKERRQ(ierr);
	ierr = DMDAGetCorners(DA_Materials_fine,0,0,0,&nx,&ny,&nz); CHKERRQ(ierr);

	/* init all coordinate vectors - just in case no coordinates are currently set on the da */
	ierr = DMDASetUniformCoordinates(DA_Quadrature,-1.0,1.0, -1.0,1.0, -1.0,1.0 ); CHKERRQ(ierr);

	/* get coordinate vector and zero it out, we will ONLY set coords on the locally owned part, ghosts are not defined */
	ierr = DMGetCoordinateDM(DA_Quadrature,			&cda); CHKERRQ(ierr);
	ierr = DMGetCoordinates(DA_Quadrature,		&gc); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,gc,			&coors); CHKERRQ(ierr);

	ierr = DMDAGetCorners(DA_Quadrature,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
	for (k=0; k<nz; k++){
		for (j=0; j<ny; j++){
			for (i=0; i<nx; i++){
				LaMEMSetMaterialDataMemoryFromArray( &material_data, i,j,k, ngp_vel_1D*ngp_vel_1D*ngp_vel_1D, materials_array );

				for (intpz=0; intpz<ngp_vel_1D; intpz++){
					for (intpy=0; intpy<ngp_vel_1D; intpy++){
						for (intpx=0; intpx<ngp_vel_1D; intpx++){
							double x,y,z;

							intp = intpz + intpy*ngp_vel_1D +  intpx*ngp_vel_1D*ngp_vel_1D;		// the correct integration point

							/* Retrieve the coordinate of the integration point */
							x = material_data.Coord[0][intp];
							y = material_data.Coord[1][intp];
							z = material_data.Coord[2][intp];

							ix = xs + i*ngp_vel_1D + intpx;
							iy = ys + j*ngp_vel_1D + intpy;
							iz = zs + k*ngp_vel_1D + intpz;

							/* Update coordinates of the Quadrature DM */

							// BUGGY FOR CERTAIN COMBINATIONS OF PROCESSORS!!!!
							coors[iz][iy][ix].x = x;
							coors[iz][iy][ix].y = y;
							coors[iz][iy][ix].z = z;
				}}}
	}}}
	ierr = DMDAVecRestoreArray(cda,gc,&coors); CHKERRQ(ierr);
	//ierr = DMDestroy(cda); CHKERRQ(ierr);
	//ierr = VecDestroy(gc);CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(DA_Materials_fine	,   Materials_fine, &materials_array); CHKERRQ(ierr);

	/* Safe way to set the coordinates */

	/* define the values of the coordinates in the ghost locations */
	ierr = DAUpdatedGhostedCoordinates(DA_Quadrature);CHKERRQ(ierr);

	PetscFunctionReturn(0);

}

//==========================================================================================================
// Set Data from the materials array to the Quadrature DA, so we can write VTS files of the quad points
#undef __FUNCT__
#define __FUNCT__ "LaMEM_SetDataFromQuadraturePointsToQuadratureDA"
PetscErrorCode LaMEM_SetDataFromQuadraturePointsToQuadratureDA(DM DA_Quadrature, PetscInt ngp_vel_1D, DM DA_Materials_fine,
		Vec Materials_fine, Vec DataArray, UserContext *user, PetscInt DataType )
{
	Vec 					DataArray_local;
	PetscInt 				nx,ny,nz, i,j,k;
	PetscInt				xs,ys,zs,xm,ym,zm,ix,iy,iz, intp, intpx, intpy, intpz;
	PetscErrorCode 			ierr;
	MaterialsElementDynamic material_data;
	PetscScalar 			***materials_array, ***data_array;


	PetscFunctionBegin;

	/* retrieve data from quad points */
	ierr = DMGetLocalVector(DA_Quadrature,		&DataArray_local);						CHKERRQ(ierr);
	ierr = VecZeroEntries(DataArray_local);												CHKERRQ(ierr);

	ierr = DMDAVecGetArray(DA_Quadrature,			DataArray_local,	&data_array); 		CHKERRQ(ierr);
	ierr = DMDAGetCorners(DA_Quadrature,			&xs,&ys,&zs,&xm,&ym,&zm); 				CHKERRQ(ierr);
	ierr = DMDAGetCorners(DA_Materials_fine,		0,0,0,&nx,&ny,&nz); 					CHKERRQ(ierr);
	ierr = DMDAVecGetArray(DA_Materials_fine	,   Materials_fine, 	&materials_array); 	CHKERRQ(ierr);

	for (k=0; k<nz; k++){
		for (j=0; j<ny; j++){
			for (i=0; i<nx; i++){
				LaMEMSetMaterialDataMemoryFromArray( &material_data, i,j,k, ngp_vel_1D*ngp_vel_1D*ngp_vel_1D, materials_array );

				for (intpz=0; intpz<ngp_vel_1D; intpz++){
					for (intpy=0; intpy<ngp_vel_1D; intpy++){
						for (intpx=0; intpx<ngp_vel_1D; intpx++){
							double data;

							intp = intpz + intpy*ngp_vel_1D +  intpx*ngp_vel_1D*ngp_vel_1D;		// the correct integration point

							/* Retrieve the data at the integration point */
							if 		(DataType==0){
								data = material_data.Viscosity[intp];
								if (user->DimensionalUnits==1){
									data = data*user->Characteristic.Viscosity;	// in Pas
								}
							}
							else if (DataType==1){
								data = material_data.Phases[intp];
							}
							else if (DataType==2){
								data = material_data.Density[intp];
								if (user->DimensionalUnits==1){
									data = data*user->Characteristic.Density;	// in kg/m3
								}
							}
							else if (DataType==3){
								data = material_data.ElasticShearModule[intp];
								if (user->DimensionalUnits==1){
									data = data*user->Characteristic.Stress;	// in Pa;
								}
							}
							else if (DataType==4){
								data = material_data.NumParticles[intp];
							}
							else if (DataType==5){
								data = material_data.Pressure[intp];
								if (user->DimensionalUnits==1){
									data = data*user->Characteristic.Stress/1e6;	// in MPa;
								}
							}
							else if (DataType==6){
								data = material_data.SecondInvariantDevStrainrate[intp];
								if (user->DimensionalUnits==1){
									data = data*(1/user->Characteristic.Time);	// in 1/s
								}
							}
							else if (DataType==7){
								data = material_data.SecondInvariantDevStress[intp];
								if (user->DimensionalUnits==1){
									data = data*user->Characteristic.Stress/1e6;	// in MPa;
								}
							}
							else if (DataType==8){
								data = material_data.Strain[intp];
							}
							else if (DataType==9){
								data = material_data.PlasticStrain[intp];
							}
							else if (DataType==10){
								data = material_data.Temperature[intp];
								if (user->DimensionalUnits==1){
									data = data*user->Characteristic.Temperature;	// in K;
								}
							}
							else {
								SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Unknown data type");
								data = 0;

							}

							ix = xs + i*ngp_vel_1D + intpx;
							iy = ys + j*ngp_vel_1D + intpy;
							iz = zs + k*ngp_vel_1D + intpz;

							/* Update coordinates of the Quadrature DM */
							data_array[iz][iy][ix] = data;
			}}} /* end of quadrature points */
	}}} /* end on elements */


	ierr = DMDAVecRestoreArray(DA_Quadrature,			DataArray_local,&data_array); 		CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(DA_Materials_fine	,   Materials_fine, &materials_array); 	CHKERRQ(ierr);

	/* only scatter the locally owned part of the vector, ignore the ghosts as we DID not fill them in */
	ierr = DMLocalToGlobalBegin(DA_Quadrature,		DataArray_local,INSERT_VALUES,DataArray);	CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(  DA_Quadrature,		DataArray_local,INSERT_VALUES,DataArray);	CHKERRQ(ierr);

	ierr = DMRestoreLocalVector(DA_Quadrature,		&DataArray_local);			CHKERRQ(ierr);


	PetscFunctionReturn(0);

}



//==========================================================================================================
#undef __FUNCT__
#define __FUNCT__ "LaMEM_DMView_3DVTK_StructuredGrid_QuadPoints"
PetscErrorCode LaMEM_DMView_3DVTK_StructuredGrid_QuadPoints(LaMEMView_QuadratureFields *view, DM DA_Materials_fine, DM DA_Processors,
		DM DA_Quadrature, Vec Materials_fine, const char file_prefix[], PetscInt ngp_vel_1D, UserContext *user, const char DirectoryName[] )
{
	char 					*vtk_filename;
	PetscMPIInt 			rank;
	MPI_Comm 				comm;
	FILE 					*vtk_fp;
	PetscInt 				si,sj,sk, nx,ny,nz, i,j,k;
	DM 						cda;
	Vec 					coords, QuadraturePoint_DataArray_local;
	DMDACoor3d 				***_coords;
	Vec 					QuadraturePoint_DataArray;
	PetscErrorCode 			ierr;
	PetscScalar 			***materials_array, ***_data;

	PetscFunctionBegin;

	/* Update local quadrature point coordinates into global coordinate vector that is attached to the DA_Quadrature DA*/
	ierr = LaMEM_SetQuadraturePointCoords_to_QuadratureDA(DA_Quadrature, ngp_vel_1D, DA_Materials_fine, Materials_fine );	CHKERRQ(ierr);

	ierr = 	DMCreateGlobalVector(DA_Quadrature,&QuadraturePoint_DataArray); CHKERRQ(ierr);


	/* create file name */
	PetscObjectGetComm( (PetscObject)DA_Processors, &comm );
	MPI_Comm_rank( comm, &rank );
	asprintf( &vtk_filename, "./%s/subdomain-%s-p%1.6d.vts", DirectoryName, file_prefix, rank );

	/* Extract Material data at quadrature points */
	ierr = DMDAVecGetArray(DA_Materials_fine,Materials_fine, &materials_array); CHKERRQ(ierr);

	/* open file and write header */
	vtk_fp = fopen( vtk_filename, "w" );
	if( vtk_fp == NULL ) {
		printf("ERROR(DAView_3DVTK_StructuredGrid_QuadPoints): Cannot open file = %s \n", vtk_filename );
		exit(0);
	}

	fprintf( vtk_fp, "<?xml version=\"1.0\"?>\n");

	/* Coordinates */
	ierr = DMDAGetGhostCorners( DA_Quadrature, &si,&sj,&sk, &nx,&ny,&nz );CHKERRQ(ierr);

#ifdef PETSC_WORDS_BIGENDIAN
	fprintf( vtk_fp, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	fprintf( vtk_fp, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif

	fprintf( vtk_fp, "\t<StructuredGrid WholeExtent=\"%lld %lld %lld %lld %lld %lld\">\n", (LLD)si, (LLD)(si+nx-1), (LLD)sj,(LLD)(sj+ny-1), (LLD)sk,(LLD)(sk+nz-1));
	fprintf( vtk_fp, "\t\t<Piece Extent=\"%lld %lld %lld %lld %lld %lld\">\n", (LLD)si, (LLD)(si+nx-1), (LLD)sj,(LLD)(sj+ny-1), (LLD)sk,(LLD)(sk+nz-1));
	fprintf( vtk_fp, "\t\t\t<CellData></CellData>\n");

	fprintf( vtk_fp, "\t\t\t<Points>\n");


	/* Add coordinates to file */
	ierr = DMGetCoordinateDM( DA_Quadrature, &cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal( DA_Quadrature,&coords );CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coords,&_coords);CHKERRQ(ierr);
	fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
	/////////
	for( k=sk; k<sk+nz; k++ ) {
		for( j=sj; j<sj+ny; j++ ) {
			for( i=si; i<si+nx; i++ ) {
				double x,y,z;
				x = _coords[k][j][i].x;
				y = _coords[k][j][i].y;
				z = _coords[k][j][i].z;
				if (user->DimensionalUnits==1){
					x = x*user->Characteristic.Length/1e3; 	// in km
					y = y*user->Characteristic.Length/1e3; 	// in km
					z = z*user->Characteristic.Length/1e3; 	// in km
				}


				fprintf( vtk_fp, "\t\t\t\t\t%1.8f %1.8f %1.8f\n", x, y, z );

			}}}
	/////////
	fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
	fprintf( vtk_fp, "\t\t\t</Points>\n");
	ierr = DMDAVecRestoreArray(cda,coords,&_coords);CHKERRQ(ierr);
	//ierr = VecDestroy(&coords);CHKERRQ(ierr);
	//ierr = DMDestroy(cda);CHKERRQ(ierr);



	fprintf( vtk_fp, "\n");
	if (user->DimensionalUnits==1){
		fprintf( vtk_fp, "\t\t\t<PointData Scalars=\"mu[Pas]");
	}
	else{
		fprintf( vtk_fp, "\t\t\t<PointData Scalars=\"mu[]");
	}
	fprintf( vtk_fp, "\">\n");

	/* fields */
	if (view->mu) {
		ierr = VecZeroEntries(QuadraturePoint_DataArray);CHKERRQ(ierr);

		LaMEM_SetDataFromQuadraturePointsToQuadratureDA(DA_Quadrature, ngp_vel_1D, DA_Materials_fine, Materials_fine, QuadraturePoint_DataArray, user, 0 );

		/* Put global data into local vector */
		ierr = DMGetLocalVector		(DA_Quadrature,&QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin	(DA_Quadrature,QuadraturePoint_DataArray,INSERT_VALUES,QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd  	(DA_Quadrature,QuadraturePoint_DataArray,INSERT_VALUES,QuadraturePoint_DataArray_local);	CHKERRQ(ierr);


		ierr = DMDAGetGhostCorners( DA_Quadrature, &si,&sj,&sk, &nx,&ny,&nz );CHKERRQ(ierr);
		ierr = DMDAVecGetArray(DA_Quadrature,QuadraturePoint_DataArray_local,&_data);CHKERRQ(ierr);

		if (user->DimensionalUnits==1){
			fprintf(vtk_fp,"\t\t\t\t<DataArray type=\"Float64\" Name=\"mu [Pas]\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		}
		else{
			fprintf(vtk_fp,"\t\t\t\t<DataArray type=\"Float64\" Name=\"mu [ ]\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		}

		/////////
		for( k=sk; k<sk+nz; k++ ) {
			for( j=sj; j<sj+ny; j++ ) {
				for( i=si; i<si+nx; i++ ) {
					double quad_data;

					quad_data = _data[k][j][i];		// retrieve data

					fprintf( vtk_fp, "\t\t\t\t\t %5.8e\n", quad_data );	// write data to file

				}}}
		/////////
		fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");

		/* Cleanup */
		ierr = DMDAVecRestoreArray(DA_Quadrature,QuadraturePoint_DataArray_local,&_data);	CHKERRQ(ierr);
		ierr = DMRestoreLocalVector	(DA_Quadrature,&QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
	}
	if (view->phase) {
		ierr = VecZeroEntries(QuadraturePoint_DataArray);CHKERRQ(ierr);

		LaMEM_SetDataFromQuadraturePointsToQuadratureDA(DA_Quadrature, ngp_vel_1D, DA_Materials_fine, Materials_fine, QuadraturePoint_DataArray, user, 1 );

		/* Put global data into local vector */
		ierr = DMGetLocalVector		(DA_Quadrature,&QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin	(DA_Quadrature,QuadraturePoint_DataArray,INSERT_VALUES,QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd  	(DA_Quadrature,QuadraturePoint_DataArray,INSERT_VALUES,QuadraturePoint_DataArray_local);	CHKERRQ(ierr);


		ierr = DMDAGetGhostCorners( DA_Quadrature, &si,&sj,&sk, &nx,&ny,&nz );CHKERRQ(ierr);
		ierr = DMDAVecGetArray(DA_Quadrature,QuadraturePoint_DataArray_local,&_data);CHKERRQ(ierr);

		fprintf(vtk_fp,"\t\t\t\t<DataArray type=\"Float64\" Name=\"phase\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		/////////
		for( k=sk; k<sk+nz; k++ ) {
			for( j=sj; j<sj+ny; j++ ) {
				for( i=si; i<si+nx; i++ ) {
					double quad_data;

					quad_data = _data[k][j][i];		// retrieve data

					fprintf( vtk_fp, "\t\t\t\t\t %1.8f\n", quad_data );	// write data to file

				}}}
		/////////
		fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");

		/* Cleanup */
		ierr = DMDAVecRestoreArray(DA_Quadrature,QuadraturePoint_DataArray_local,&_data);	CHKERRQ(ierr);
		ierr = DMRestoreLocalVector	(DA_Quadrature,&QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
	}
	if (view->rho) {
		ierr = VecZeroEntries(QuadraturePoint_DataArray);CHKERRQ(ierr);

		LaMEM_SetDataFromQuadraturePointsToQuadratureDA(DA_Quadrature, ngp_vel_1D, DA_Materials_fine, Materials_fine, QuadraturePoint_DataArray, user, 2 );

		/* Put global data into local vector */
		ierr = DMGetLocalVector		(DA_Quadrature,&QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin	(DA_Quadrature,QuadraturePoint_DataArray,INSERT_VALUES,QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd  	(DA_Quadrature,QuadraturePoint_DataArray,INSERT_VALUES,QuadraturePoint_DataArray_local);	CHKERRQ(ierr);


		ierr = DMDAGetGhostCorners( DA_Quadrature, &si,&sj,&sk, &nx,&ny,&nz );CHKERRQ(ierr);
		ierr = DMDAVecGetArray(DA_Quadrature,QuadraturePoint_DataArray_local,&_data);CHKERRQ(ierr);

		if (user->DimensionalUnits==1){
			fprintf(vtk_fp,"\t\t\t\t<DataArray type=\"Float64\" Name=\"rho [kg/m3]\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		}
		else{
			fprintf(vtk_fp,"\t\t\t\t<DataArray type=\"Float64\" Name=\"rho [ ]\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		}
		/////////
		for( k=sk; k<sk+nz; k++ ) {
			for( j=sj; j<sj+ny; j++ ) {
				for( i=si; i<si+nx; i++ ) {
					double quad_data;

					quad_data = _data[k][j][i];		// retrieve data

					fprintf( vtk_fp, "\t\t\t\t\t %5.8e\n", quad_data );	// write data to file

				}}}
		/////////
		fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");

		/* Cleanup */
		ierr = DMDAVecRestoreArray(DA_Quadrature,QuadraturePoint_DataArray_local,&_data);	CHKERRQ(ierr);
		ierr = DMRestoreLocalVector	(DA_Quadrature,&QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
	}

	if (view->G) {
		ierr = VecZeroEntries(QuadraturePoint_DataArray);CHKERRQ(ierr);

		LaMEM_SetDataFromQuadraturePointsToQuadratureDA(DA_Quadrature, ngp_vel_1D, DA_Materials_fine, Materials_fine, QuadraturePoint_DataArray, user, 3 );

		/* Put global data into local vector */
		ierr = DMGetLocalVector		(DA_Quadrature,&QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin	(DA_Quadrature,QuadraturePoint_DataArray,INSERT_VALUES,QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd  	(DA_Quadrature,QuadraturePoint_DataArray,INSERT_VALUES,QuadraturePoint_DataArray_local);	CHKERRQ(ierr);


		ierr = DMDAGetGhostCorners( DA_Quadrature, &si,&sj,&sk, &nx,&ny,&nz );CHKERRQ(ierr);
		ierr = DMDAVecGetArray(DA_Quadrature,QuadraturePoint_DataArray_local,&_data);CHKERRQ(ierr);
		if (user->DimensionalUnits==1){
			fprintf(vtk_fp,"\t\t\t\t<DataArray type=\"Float64\" Name=\"G [Pa]\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		}
		else{
			fprintf(vtk_fp,"\t\t\t\t<DataArray type=\"Float64\" Name=\"G [ ]\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		}
		/////////
		for( k=sk; k<sk+nz; k++ ) {
			for( j=sj; j<sj+ny; j++ ) {
				for( i=si; i<si+nx; i++ ) {
					double quad_data;

					quad_data = _data[k][j][i];		// retrieve data

					fprintf( vtk_fp, "\t\t\t\t\t %5.8e\n", quad_data );	// write data to file

				}}}
		/////////
		fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");

		/* Cleanup */
		ierr = DMDAVecRestoreArray(DA_Quadrature,QuadraturePoint_DataArray_local,&_data);	CHKERRQ(ierr);
		ierr = DMRestoreLocalVector	(DA_Quadrature,&QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
	}


	if (view->numParticles) {
		ierr = VecZeroEntries(QuadraturePoint_DataArray);CHKERRQ(ierr);

		LaMEM_SetDataFromQuadraturePointsToQuadratureDA(DA_Quadrature, ngp_vel_1D, DA_Materials_fine, Materials_fine, QuadraturePoint_DataArray, user, 4 );

		/* Put global data into local vector */
		ierr = DMGetLocalVector		(DA_Quadrature,&QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin	(DA_Quadrature,QuadraturePoint_DataArray,INSERT_VALUES,QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd  	(DA_Quadrature,QuadraturePoint_DataArray,INSERT_VALUES,QuadraturePoint_DataArray_local);	CHKERRQ(ierr);


		ierr = DMDAGetGhostCorners( DA_Quadrature, &si,&sj,&sk, &nx,&ny,&nz );CHKERRQ(ierr);
		ierr = DMDAVecGetArray(DA_Quadrature,QuadraturePoint_DataArray_local,&_data);CHKERRQ(ierr);

		fprintf(vtk_fp,"\t\t\t\t<DataArray type=\"Float64\" Name=\"NumParticles\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		/////////
		for( k=sk; k<sk+nz; k++ ) {
			for( j=sj; j<sj+ny; j++ ) {
				for( i=si; i<si+nx; i++ ) {
					double quad_data;

					quad_data = _data[k][j][i];		// retrieve data

					fprintf( vtk_fp, "\t\t\t\t\t %1.8f\n", quad_data );	// write data to file

				}}}
		/////////
		fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");

		/* Cleanup */
		ierr = DMDAVecRestoreArray(DA_Quadrature,QuadraturePoint_DataArray_local,&_data);	CHKERRQ(ierr);
		ierr = DMRestoreLocalVector	(DA_Quadrature,&QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
	}
	if (view->pressure) {
		ierr = VecZeroEntries(QuadraturePoint_DataArray);CHKERRQ(ierr);

		LaMEM_SetDataFromQuadraturePointsToQuadratureDA(DA_Quadrature, ngp_vel_1D, DA_Materials_fine, Materials_fine, QuadraturePoint_DataArray, user, 5 );

		/* Put global data into local vector */
		ierr = DMGetLocalVector		(DA_Quadrature,&QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin	(DA_Quadrature,QuadraturePoint_DataArray,INSERT_VALUES,QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd  	(DA_Quadrature,QuadraturePoint_DataArray,INSERT_VALUES,QuadraturePoint_DataArray_local);	CHKERRQ(ierr);


		ierr = DMDAGetGhostCorners( DA_Quadrature, &si,&sj,&sk, &nx,&ny,&nz );CHKERRQ(ierr);
		ierr = DMDAVecGetArray(DA_Quadrature,QuadraturePoint_DataArray_local,&_data);CHKERRQ(ierr);
		if (user->DimensionalUnits==1){
			fprintf(vtk_fp,"\t\t\t\t<DataArray type=\"Float64\" Name=\"Pressure [MPa]\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		}
		else {
			fprintf(vtk_fp,"\t\t\t\t<DataArray type=\"Float64\" Name=\"Pressure [ ]\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		}
		/////////
		for( k=sk; k<sk+nz; k++ ) {
			for( j=sj; j<sj+ny; j++ ) {
				for( i=si; i<si+nx; i++ ) {
					double quad_data;

					quad_data = _data[k][j][i];		// retrieve data

					fprintf( vtk_fp, "\t\t\t\t\t %5.8e\n", quad_data );	// write data to file

				}}}
		/////////
		fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");

		/* Cleanup */
		ierr = DMDAVecRestoreArray(DA_Quadrature,QuadraturePoint_DataArray_local,&_data);	CHKERRQ(ierr);
		ierr = DMRestoreLocalVector	(DA_Quadrature,&QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
	}

	if (view->temperature) {
		ierr = VecZeroEntries(QuadraturePoint_DataArray);CHKERRQ(ierr);

		LaMEM_SetDataFromQuadraturePointsToQuadratureDA(DA_Quadrature, ngp_vel_1D, DA_Materials_fine, Materials_fine, QuadraturePoint_DataArray, user, 10 );

		/* Put global data into local vector */
		ierr = DMGetLocalVector		(DA_Quadrature,&QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin	(DA_Quadrature,QuadraturePoint_DataArray,INSERT_VALUES,QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd  	(DA_Quadrature,QuadraturePoint_DataArray,INSERT_VALUES,QuadraturePoint_DataArray_local);	CHKERRQ(ierr);


		ierr = DMDAGetGhostCorners( DA_Quadrature, &si,&sj,&sk, &nx,&ny,&nz );CHKERRQ(ierr);
		ierr = DMDAVecGetArray(DA_Quadrature,QuadraturePoint_DataArray_local,&_data);CHKERRQ(ierr);
		if (user->DimensionalUnits==1){
			fprintf(vtk_fp,"\t\t\t\t<DataArray type=\"Float64\" Name=\"Temperature [K]\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		}
		else {
			fprintf(vtk_fp,"\t\t\t\t<DataArray type=\"Float64\" Name=\"Temperature [ ]\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		}
		/////////
		for( k=sk; k<sk+nz; k++ ) {
			for( j=sj; j<sj+ny; j++ ) {
				for( i=si; i<si+nx; i++ ) {
					double quad_data;

					quad_data = _data[k][j][i];		// retrieve data

					fprintf( vtk_fp, "\t\t\t\t\t %5.8e\n", quad_data );	// write data to file

				}}}
		/////////
		fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");

		/* Cleanup */
		ierr = DMDAVecRestoreArray(DA_Quadrature,QuadraturePoint_DataArray_local,&_data);	CHKERRQ(ierr);
		ierr = DMRestoreLocalVector	(DA_Quadrature,&QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
	}

	if (view->devStrainrateInv) {
		ierr = VecZeroEntries(QuadraturePoint_DataArray);CHKERRQ(ierr);

		LaMEM_SetDataFromQuadraturePointsToQuadratureDA(DA_Quadrature, ngp_vel_1D, DA_Materials_fine, Materials_fine, QuadraturePoint_DataArray, user, 6 );

		/* Put global data into local vector */
		ierr = DMGetLocalVector		(DA_Quadrature,&QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin	(DA_Quadrature,QuadraturePoint_DataArray,INSERT_VALUES,QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd  	(DA_Quadrature,QuadraturePoint_DataArray,INSERT_VALUES,QuadraturePoint_DataArray_local);	CHKERRQ(ierr);


		ierr = DMDAGetGhostCorners( DA_Quadrature, &si,&sj,&sk, &nx,&ny,&nz );CHKERRQ(ierr);
		ierr = DMDAVecGetArray(DA_Quadrature,QuadraturePoint_DataArray_local,&_data);CHKERRQ(ierr);
		if (user->DimensionalUnits==1){
			fprintf(vtk_fp,"\t\t\t\t<DataArray type=\"Float64\" Name=\"DevStrainrateInv [1/s]\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		}
		else {
			fprintf(vtk_fp,"\t\t\t\t<DataArray type=\"Float64\" Name=\"DevStrainrateInv [ ]\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		}
		/////////
		for( k=sk; k<sk+nz; k++ ) {
			for( j=sj; j<sj+ny; j++ ) {
				for( i=si; i<si+nx; i++ ) {
					double quad_data;

					quad_data = _data[k][j][i];		// retrieve data

					fprintf( vtk_fp, "\t\t\t\t\t %5.8e\n", quad_data );	// write data to file

				}}}
		/////////
		fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");

		/* Cleanup */
		ierr = DMDAVecRestoreArray(DA_Quadrature,QuadraturePoint_DataArray_local,&_data);	CHKERRQ(ierr);
		ierr = DMRestoreLocalVector	(DA_Quadrature,&QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
	}


	if (view->devStressInv) {
		ierr = VecZeroEntries(QuadraturePoint_DataArray);CHKERRQ(ierr);

		LaMEM_SetDataFromQuadraturePointsToQuadratureDA(DA_Quadrature, ngp_vel_1D, DA_Materials_fine, Materials_fine, QuadraturePoint_DataArray, user, 7 );

		/* Put global data into local vector */
		ierr = DMGetLocalVector		(DA_Quadrature,&QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin	(DA_Quadrature,QuadraturePoint_DataArray,INSERT_VALUES,QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd  	(DA_Quadrature,QuadraturePoint_DataArray,INSERT_VALUES,QuadraturePoint_DataArray_local);	CHKERRQ(ierr);


		ierr = DMDAGetGhostCorners( DA_Quadrature, &si,&sj,&sk, &nx,&ny,&nz );CHKERRQ(ierr);
		ierr = DMDAVecGetArray(DA_Quadrature,QuadraturePoint_DataArray_local,&_data);CHKERRQ(ierr);
		if (user->DimensionalUnits==1){
			fprintf(vtk_fp,"\t\t\t\t<DataArray type=\"Float64\" Name=\"DevStressInv [MPa]\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		}
		else{
			fprintf(vtk_fp,"\t\t\t\t<DataArray type=\"Float64\" Name=\"DevStressInv [ ]\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		}
		/////////
		for( k=sk; k<sk+nz; k++ ) {
			for( j=sj; j<sj+ny; j++ ) {
				for( i=si; i<si+nx; i++ ) {
					double quad_data;

					quad_data = _data[k][j][i];		// retrieve data

					fprintf( vtk_fp, "\t\t\t\t\t %5.8e\n", quad_data );	// write data to file

				}}}
		/////////
		fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");

		/* Cleanup */
		ierr = DMDAVecRestoreArray(DA_Quadrature,QuadraturePoint_DataArray_local,&_data);	CHKERRQ(ierr);
		ierr = DMRestoreLocalVector	(DA_Quadrature,&QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
	}

	if (view->strain) {
		ierr = VecZeroEntries(QuadraturePoint_DataArray);CHKERRQ(ierr);

		LaMEM_SetDataFromQuadraturePointsToQuadratureDA(DA_Quadrature, ngp_vel_1D, DA_Materials_fine, Materials_fine, QuadraturePoint_DataArray, user, 8 );

		/* Put global data into local vector */
		ierr = DMGetLocalVector		(DA_Quadrature,&QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin	(DA_Quadrature,QuadraturePoint_DataArray,INSERT_VALUES,QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd  	(DA_Quadrature,QuadraturePoint_DataArray,INSERT_VALUES,QuadraturePoint_DataArray_local);	CHKERRQ(ierr);


		ierr = DMDAGetGhostCorners( DA_Quadrature, &si,&sj,&sk, &nx,&ny,&nz );CHKERRQ(ierr);
		ierr = DMDAVecGetArray(DA_Quadrature,QuadraturePoint_DataArray_local,&_data);CHKERRQ(ierr);

		fprintf(vtk_fp,"\t\t\t\t<DataArray type=\"Float64\" Name=\"Strain [ ]\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		/////////
		for( k=sk; k<sk+nz; k++ ) {
			for( j=sj; j<sj+ny; j++ ) {
				for( i=si; i<si+nx; i++ ) {
					double quad_data;

					quad_data = _data[k][j][i];		// retrieve data

					fprintf( vtk_fp, "\t\t\t\t\t %5.8e\n", quad_data );	// write data to file

				}}}
		/////////
		fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");

		/* Cleanup */
		ierr = DMDAVecRestoreArray(DA_Quadrature,QuadraturePoint_DataArray_local,&_data);	CHKERRQ(ierr);
		ierr = DMRestoreLocalVector	(DA_Quadrature,&QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
	}

	if (view->plasticStrain) {
		ierr = VecZeroEntries(QuadraturePoint_DataArray);CHKERRQ(ierr);

		LaMEM_SetDataFromQuadraturePointsToQuadratureDA(DA_Quadrature, ngp_vel_1D, DA_Materials_fine, Materials_fine, QuadraturePoint_DataArray, user, 9 );

		/* Put global data into local vector */
		ierr = DMGetLocalVector		(DA_Quadrature,&QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin	(DA_Quadrature,QuadraturePoint_DataArray,INSERT_VALUES,QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd  	(DA_Quadrature,QuadraturePoint_DataArray,INSERT_VALUES,QuadraturePoint_DataArray_local);	CHKERRQ(ierr);


		ierr = DMDAGetGhostCorners( DA_Quadrature, &si,&sj,&sk, &nx,&ny,&nz );CHKERRQ(ierr);
		ierr = DMDAVecGetArray(DA_Quadrature,QuadraturePoint_DataArray_local,&_data);CHKERRQ(ierr);

		fprintf(vtk_fp,"\t\t\t\t<DataArray type=\"Float64\" Name=\"PlasticStrain [ ]\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		/////////
		for( k=sk; k<sk+nz; k++ ) {
			for( j=sj; j<sj+ny; j++ ) {
				for( i=si; i<si+nx; i++ ) {
					double quad_data;

					quad_data = _data[k][j][i];		// retrieve data

					fprintf( vtk_fp, "\t\t\t\t\t %5.8e\n", quad_data );	// write data to file

				}}}
		/////////
		fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");

		/* Cleanup */
		ierr = DMDAVecRestoreArray(DA_Quadrature,QuadraturePoint_DataArray_local,&_data);	CHKERRQ(ierr);
		ierr = DMRestoreLocalVector	(DA_Quadrature,&QuadraturePoint_DataArray_local);	CHKERRQ(ierr);
	}


	fprintf( vtk_fp, "\t\t\t</PointData>\n");

	fprintf( vtk_fp, "\t\t</Piece>\n");
	fprintf( vtk_fp, "\t</StructuredGrid>\n");
	fprintf( vtk_fp, "</VTKFile>\n");

	// close file
	if( vtk_fp!= NULL ) {
		fclose( vtk_fp );
		vtk_fp = NULL;
	}

	// clean_up
	free(vtk_filename);


	ierr = DMDAVecRestoreArray(DA_Materials_fine,Materials_fine, &materials_array); CHKERRQ(ierr);
	ierr = VecDestroy(&QuadraturePoint_DataArray); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

//==========================================================================================================
// Create a structured grid for 3 component vectors (e.g. velocity)
#undef __FUNCT__
#define __FUNCT__ "DMView_3DVTK_StructuredGrid_3ComponentVector"
PetscErrorCode DMView_3DVTK_StructuredGrid_3ComponentVector( DM da, Vec FIELD, const char file_prefix[], const char DirectoryName[], const char vectorfield_name[], PetscScalar scaling_length)
{
	char 				*vtk_filename;
	PetscMPIInt 		rank;
	MPI_Comm 			comm;
	FILE 				*vtk_fp;
	PetscInt 			si,sj,sk, nx,ny,nz, i,j,k;
	PetscInt 			n_fields, N;
	DM 					cda;
	Vec 				coords;
	DMDACoor3d			***_coords;
	Vec 				l_FIELD;
	PetscScalar 		*_L_FIELD;
	PetscErrorCode 		ierr;

	PetscFunctionBegin;

	// create file name
	PetscObjectGetComm( (PetscObject)da, &comm );
	MPI_Comm_rank( comm, &rank );
	asprintf( &vtk_filename, "./%s/subdomain-%s-p%1.6d.vts", DirectoryName, file_prefix, rank );

	// open file and write header
	vtk_fp = fopen( vtk_filename, "w" );
	if( vtk_fp == NULL ) {
		printf("ERROR(DMView_3DVTK_StructuredGrid_3ComponentVectorView_3DVTK_StructuredGrid_3ComponentVector): Cannot open file = %s \n", vtk_filename );
		exit(0);
	}

	fprintf( vtk_fp, "<?xml version=\"1.0\"?>\n");

	// coords
	ierr = DMDAGetGhostCorners( da, &si,&sj,&sk, &nx,&ny,&nz );CHKERRQ(ierr);
	N = nx * ny * nz;

#ifdef PETSC_WORDS_BIGENDIAN
	fprintf( vtk_fp, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	fprintf( vtk_fp, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
	fprintf( vtk_fp, "\t<StructuredGrid WholeExtent=\"%lld %lld %lld %lld %lld %lld\">\n", (LLD)si,(LLD)(si+nx-1), (LLD)sj,(LLD)(sj+ny-1), (LLD)sk,(LLD)(sk+nz-1));
	fprintf( vtk_fp, "\t\t<Piece Extent=\"%lld %lld %lld %lld %lld %lld\">\n", (LLD)si,(LLD)(si+nx-1), (LLD)sj,(LLD)(sj+ny-1), (LLD)sk,(LLD)(sk+nz-1));
	fprintf( vtk_fp, "\t\t\t<CellData></CellData>\n");
	fprintf( vtk_fp, "\t\t\t<Points>\n");

	// copy coordinates
	ierr = DMGetCoordinateDM( da, &cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal( da,&coords );CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coords,&_coords);CHKERRQ(ierr);

	fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
	for( k=sk; k<sk+nz; k++ ) {
		for( j=sj; j<sj+ny; j++ ) {
			for( i=si; i<si+nx; i++ ) {
				double x,y,z;

				x = _coords[k][j][i].x * scaling_length;
				y = _coords[k][j][i].y * scaling_length;
				z = _coords[k][j][i].z * scaling_length;

				fprintf( vtk_fp, "\t\t\t\t\t%1.8f %1.8f %1.8f\n", x, y, z );
			}
		}
	}
	fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
	fprintf( vtk_fp, "\t\t\t</Points>\n");
	fprintf( vtk_fp, "\t\t\t<PointData Vectors=\"%s\">\n",vectorfield_name);

	ierr = DMDAVecRestoreArray(cda,coords,&_coords);CHKERRQ(ierr);
	//ierr = VecDestroy(&coords);CHKERRQ(ierr);
	//ierr = DMDestroy(cda);CHKERRQ(ierr);



	// check that the field indeed has 3 components
	ierr = DMDAGetInfo( da, 0, 0,0,0, 0,0,0, &n_fields, 0, 0, 0,0,0 );
	if (n_fields != 3){
		printf("ERROR(DMView_3DVTK_StructuredGrid_3ComponentVector): The DM has %lld fields and not 3 fields \n", (LLD)n_fields );
		exit(0);
	}

	// Write the 3-component vector field
	ierr = DMGetLocalVector( da, &l_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin( da, FIELD, INSERT_VALUES, l_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd( da, FIELD, INSERT_VALUES, l_FIELD );CHKERRQ(ierr);
	ierr = VecGetArray( l_FIELD, &_L_FIELD );CHKERRQ(ierr);

	fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"3\" format=\"ascii\">\n",vectorfield_name );
	for( i=0; i<N; i++ ) {
		double val0, val1, val2;

		val0 = _L_FIELD[ n_fields*i + 0 ];
		val1 = _L_FIELD[ n_fields*i + 1 ];
		val2 = _L_FIELD[ n_fields*i + 2 ];

		fprintf( vtk_fp, "\t\t\t\t\t%1.8f %1.8f %1.8f \n", val0, val1, val2 );
	}
	fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");

	ierr = VecRestoreArray( l_FIELD, &_L_FIELD );CHKERRQ(ierr);
	ierr = DMRestoreLocalVector( da, &l_FIELD );CHKERRQ(ierr);

	fprintf( vtk_fp, "\t\t\t</PointData>\n");
	fprintf( vtk_fp, "\t\t</Piece>\n");
	fprintf( vtk_fp, "\t</StructuredGrid>\n");
	fprintf( vtk_fp, "</VTKFile>\n");

	// close file
	if( vtk_fp!= NULL ) {
		fclose( vtk_fp );
		vtk_fp = NULL;
	}

	// clean up
	free(vtk_filename);

	PetscFunctionReturn(0);
}


//==========================================================================================================
// Create a 2D structured grid to write topography data
#undef __FUNCT__
#define __FUNCT__ "DMView_2DVTK_StructuredGrid_Topo"
PetscErrorCode DMView_2DVTK_StructuredGrid_Topo( DM da, Vec FIELD, const char file_prefix[],  const char DirectoryName[], Vec SurfaceVelocity_Vx, Vec SurfaceVelocity_Vy, Vec SurfaceVelocity_Vz,PetscScalar scaling_length )

{
	char 			*vtk_filename;
	PetscMPIInt 	rank;
	MPI_Comm 		comm;
	FILE 			*vtk_fp;
	PetscInt		si,sj,sk, nx,ny,nz, i,j,k;
	PetscInt		f, n_fields, size;
	DM				cda;
	Vec 			coords;
	DMDACoor3d 		***_coords;
	Vec				l_FIELD;
	PetscScalar		***ZCOORD;
	PetscScalar 	AverageHeight;
	PetscErrorCode 	ierr;


	PetscFunctionBegin;


	// create file name
	PetscObjectGetComm( (PetscObject)da, &comm );
	MPI_Comm_rank( comm, &rank );
	asprintf( &vtk_filename, "./%s/subdomain-%s-p%1.6d.vts", DirectoryName, file_prefix, rank );

//	printf("DMView_3DVTK_StructuredGrid: file = %s, rank = %lld \n", vtk_filename,rank );


	// open file and write header
	vtk_fp = fopen( vtk_filename, "w" );
	if( vtk_fp == NULL ) {
		printf("ERROR(DAView_2DVTK_StructuredGrid): Cannot open file = %s \n", vtk_filename );
		exit(0);
	}

	fprintf( vtk_fp, "<?xml version=\"1.0\"?>\n");

	// coords
	ierr = DMDAGetGhostCorners( da, &si,&sj,&sk, &nx,&ny,&nz );CHKERRQ(ierr);

	fprintf( vtk_fp, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	fprintf( vtk_fp, "\t<StructuredGrid WholeExtent=\"%lld %lld %lld %lld %lld %lld\">\n", (LLD)si,(LLD)(si+nx-1), (LLD)sj,(LLD)(sj+ny-1), 0LL,0LL);
	fprintf( vtk_fp, "\t\t<Piece Extent=\"%lld %lld %lld %lld %lld %lld\">\n", (LLD)si,(LLD)(si+nx-1), (LLD)sj,(LLD)(sj+ny-1), (LLD)sk,(LLD)sk);
	fprintf( vtk_fp, "\t\t\t<CellData></CellData>\n");
	fprintf( vtk_fp, "\t\t\t<Points>\n");

	//copy coordinates
	ierr = DMGetCoordinateDM( da, &cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal( da,&coords );CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coords,&_coords);CHKERRQ(ierr);

	ierr = DMGetLocalVector( da, &l_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin( da, FIELD, INSERT_VALUES, l_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd( da, FIELD, INSERT_VALUES, l_FIELD );CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, l_FIELD, &ZCOORD );CHKERRQ(ierr);

	//ierr = VecGetArray( l_FIELD, &_L_FIELD );CHKERRQ(ierr);

	fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
//	for( k=sk+nz-1; k<sk+nz; k++ ) {
		for( k=sk; k<(sk+nz); k++ ) {
			for( j=sj; j<sj+ny; j++ ) {
				for( i=si; i<si+nx; i++ ) {
					double x,y,z;

					x = _coords[k][j][i].x 	* scaling_length;
					y = _coords[k][j][i].y 	* scaling_length;
					z = ZCOORD[k][j][i]		* scaling_length;
					fprintf( vtk_fp, "\t\t\t\t\t%1.8f %1.8f %1.8f\n", x, y, z );
				}
			}
		}
	fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
	fprintf( vtk_fp, "\t\t\t</Points>\n");
	ierr = DMDAVecRestoreArray(cda,coords,&_coords);CHKERRQ(ierr);
	//ierr = VecDestroy(&coords);CHKERRQ(ierr);
	//ierr = DMDestroy(cda);CHKERRQ(ierr);

	//ierr = VecRestoreArray( l_FIELD, &_L_FIELD );CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da, l_FIELD, &ZCOORD );CHKERRQ(ierr);
	ierr = DMRestoreLocalVector( da, &l_FIELD );CHKERRQ(ierr);


	fprintf( vtk_fp, "\t\t\t<PointData Scalars=\" ");
	ierr = DMDAGetInfo( da, 0, 0,0,0, 0,0,0, &n_fields, 0, 0, 0,0,0 );
	for( f=0; f<n_fields; f++ ) {
		const char *field_name;
		ierr = DMDAGetFieldName( da, f, &field_name );CHKERRQ(ierr);
		fprintf( vtk_fp, "%s ", field_name );
	}
	fprintf( vtk_fp, "%s ", "SurfaceVelocity_Vx[]" );
	fprintf( vtk_fp, "%s ", "SurfaceVelocity_Vy[]" );
	fprintf( vtk_fp, "%s ", "SurfaceVelocity_Vz[]" );
	fprintf( vtk_fp, "%s ", "SurfaceAmplitude" );

	fprintf( vtk_fp, "\">\n");


	ierr = DMGetLocalVector( da, &l_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin( da, FIELD, INSERT_VALUES, l_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd( da, FIELD, INSERT_VALUES, l_FIELD );CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, l_FIELD, &ZCOORD );CHKERRQ(ierr);


	{
		const char *field_name;
		ierr = DMDAGetFieldName( da, 0, &field_name );CHKERRQ(ierr);
		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", field_name );
		for( k=sk+nz-1; k<sk+nz; k++ ) {
			for( j=sj; j<sj+ny; j++ ) {
				for( i=si; i<si+nx; i++ ) {
					double val;
					val = ZCOORD[k][j][i] 	* scaling_length;
					fprintf( vtk_fp, "\t\t\t\t\t%1.8f\n", val );
				}
			}
		}
	}
	fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
	ierr = DMDAVecRestoreArray(da, l_FIELD, &ZCOORD );CHKERRQ(ierr);
	ierr = DMRestoreLocalVector( da, &l_FIELD );CHKERRQ(ierr);


	/* Write Vx surface velocity in nondimensional units*/
	ierr = DMGetLocalVector( da, &l_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin( da, SurfaceVelocity_Vx, INSERT_VALUES, l_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd( da, SurfaceVelocity_Vx, INSERT_VALUES, l_FIELD );CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, l_FIELD, &ZCOORD );CHKERRQ(ierr);


	{
		const char *field_name;
		ierr = DMDAGetFieldName( da, 0, &field_name );CHKERRQ(ierr);
		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", "SurfaceVelocity_Vx[]" );
		for( k=sk+nz-1; k<sk+nz; k++ ) {
			for( j=sj; j<sj+ny; j++ ) {
				for( i=si; i<si+nx; i++ ) {
					double val;
					val = ZCOORD[k][j][i];
					fprintf( vtk_fp, "\t\t\t\t\t%1.8f\n", val );
				}
			}
		}

	}
	fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
	ierr = DMDAVecRestoreArray(da, l_FIELD, &ZCOORD );CHKERRQ(ierr);
	ierr = DMRestoreLocalVector( da, &l_FIELD );CHKERRQ(ierr);

	/* Write Vy surface velocity in nondimensional units*/
	ierr = DMGetLocalVector( da, &l_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin( da, SurfaceVelocity_Vy, INSERT_VALUES, l_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd( da, SurfaceVelocity_Vy, INSERT_VALUES, l_FIELD );CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, l_FIELD, &ZCOORD );CHKERRQ(ierr);


	{
		const char *field_name;
		ierr = DMDAGetFieldName( da, 0, &field_name );CHKERRQ(ierr);
		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", "SurfaceVelocity_Vy[]" );
		for( k=sk+nz-1; k<sk+nz; k++ ) {
			for( j=sj; j<sj+ny; j++ ) {
				for( i=si; i<si+nx; i++ ) {
					double val;
					val = ZCOORD[k][j][i];
					fprintf( vtk_fp, "\t\t\t\t\t%1.8f\n", val );
				}
			}
		}

	}
	fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
	ierr = DMDAVecRestoreArray(da, l_FIELD, &ZCOORD );CHKERRQ(ierr);
	ierr = DMRestoreLocalVector( da, &l_FIELD );CHKERRQ(ierr);

	/* Write Vz surface velocity in nondimensional units*/
	ierr = DMGetLocalVector( da, &l_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin( da, SurfaceVelocity_Vz, INSERT_VALUES, l_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd( da, SurfaceVelocity_Vz, INSERT_VALUES, l_FIELD );CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, l_FIELD, &ZCOORD );CHKERRQ(ierr);


	{
		const char *field_name;
		ierr = DMDAGetFieldName( da, 0, &field_name );CHKERRQ(ierr);
		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", "SurfaceVelocity_Vz[]" );
		for( k=sk+nz-1; k<sk+nz; k++ ) {
			for( j=sj; j<sj+ny; j++ ) {
				for( i=si; i<si+nx; i++ ) {
					double val;
					val = ZCOORD[k][j][i];
					fprintf( vtk_fp, "\t\t\t\t\t%1.8f\n", val );
				}
			}
		}

	}
	fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
	ierr = DMDAVecRestoreArray(da, l_FIELD, &ZCOORD );CHKERRQ(ierr);
	ierr = DMRestoreLocalVector( da, &l_FIELD );CHKERRQ(ierr);


	/* Write the amplitude of the surface deformation; usefull to look at folding */
	ierr = DMGetLocalVector( da, &l_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin( da, FIELD, INSERT_VALUES, l_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd( da, FIELD, INSERT_VALUES, l_FIELD );CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, l_FIELD, &ZCOORD );CHKERRQ(ierr);


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

	ierr = DMDAVecRestoreArray(da, l_FIELD, &ZCOORD );CHKERRQ(ierr);
	ierr = DMRestoreLocalVector( da, &l_FIELD );CHKERRQ(ierr);


	fprintf( vtk_fp, "\t\t\t</PointData>\n");
	fprintf( vtk_fp, "\t\t</Piece>\n");
	fprintf( vtk_fp, "\t</StructuredGrid>\n");
	fprintf( vtk_fp, "</VTKFile>\n");


	if( vtk_fp!= NULL ) {
		fclose( vtk_fp );
		vtk_fp = NULL;
	}

	free(vtk_filename);

	PetscFunctionReturn(0);
}


//==========================================================================================================
// Create a 2D structured grid to write topography data used for the 2D erosion code
#undef __FUNCT__
#define __FUNCT__ "DMView_2DVTK_StructuredGrid_Topo_Erosion"
PetscErrorCode DMView_2DVTK_StructuredGrid_Topo_Erosion( DM da, Vec FIELD, const char file_prefix[],  const char DirectoryName[], PetscScalar scaling_length )

{
	char 			*vtk_filename;
	PetscMPIInt 	rank;
	MPI_Comm 		comm;
	FILE 			*vtk_fp;
	PetscInt		si,sj, nx,ny, i,j;
	PetscInt		n_fields, size;
	DM				cda;
	Vec 			coords;
	DMDACoor2d 		**_coords;
	Vec				l_FIELD;
	PetscScalar		**ZCOORD;
	PetscScalar 	AverageHeight;
	PetscErrorCode 	ierr;


	PetscFunctionBegin;


	// create file name
	PetscObjectGetComm( (PetscObject)da, &comm );
	MPI_Comm_rank( comm, &rank );
	asprintf( &vtk_filename, "./%s/subdomain-%s-p%1.6d.vts", DirectoryName, file_prefix, rank );

//	printf("DMView_3DVTK_StructuredGrid: file = %s, rank = %lld \n", vtk_filename,rank );


	// open file and write header
	vtk_fp = fopen( vtk_filename, "w" );
	if( vtk_fp == NULL ) {
		printf("ERROR(DAView_2DVTK_StructuredGrid): Cannot open file = %s \n", vtk_filename );
		exit(0);
	}

	fprintf( vtk_fp, "<?xml version=\"1.0\"?>\n");

	// coords
	ierr = DMDAGetGhostCorners( da, &si,&sj,PETSC_NULL, &nx,&ny,PETSC_NULL);CHKERRQ(ierr);

	fprintf( vtk_fp, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	fprintf( vtk_fp, "\t<StructuredGrid WholeExtent=\"%lld %lld %lld %lld %lld %lld\">\n", (LLD)si,(LLD)(si+nx-1), (LLD)sj,(LLD)(sj+ny-1), 0LL,0LL);
	fprintf( vtk_fp, "\t\t<Piece Extent=\"%lld %lld %lld %lld %lld %lld\">\n", (LLD)si,(LLD)(si+nx-1), (LLD)sj,(LLD)(sj+ny-1), (LLD)0,(LLD)0);
	fprintf( vtk_fp, "\t\t\t<CellData></CellData>\n");
	fprintf( vtk_fp, "\t\t\t<Points>\n");


	//copy coordinates
	ierr = DMGetCoordinateDM( da, &cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal( da,&coords );CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coords,&_coords);CHKERRQ(ierr);

	ierr = DMGetLocalVector( da, &l_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin( da, FIELD, INSERT_VALUES, l_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd( da, FIELD, INSERT_VALUES, l_FIELD );CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, l_FIELD, &ZCOORD );CHKERRQ(ierr);

	//ierr = VecGetArray( l_FIELD, &_L_FIELD );CHKERRQ(ierr);

	fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
//	for( k=sk+nz-1; k<sk+nz; k++ ) {
			for( j=sj; j<sj+ny; j++ ) {
				for( i=si; i<si+nx; i++ ) {
					double x,y,z;

					x = _coords[j][i].x 	* scaling_length;
					y = _coords[j][i].y 	* scaling_length;
					z = ZCOORD[j][i]		* scaling_length;
					fprintf( vtk_fp, "\t\t\t\t\t%1.8f %1.8f %1.8f\n", x, y, z );
				}
			}
	fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
	fprintf( vtk_fp, "\t\t\t</Points>\n");
	ierr = DMDAVecRestoreArray(cda,coords,&_coords);CHKERRQ(ierr);
	//ierr = VecDestroy(&coords);CHKERRQ(ierr);
	//ierr = DMDestroy(cda);CHKERRQ(ierr);

	//ierr = VecRestoreArray( l_FIELD, &_L_FIELD );CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da, l_FIELD, &ZCOORD );CHKERRQ(ierr);
	ierr = DMRestoreLocalVector( da, &l_FIELD );CHKERRQ(ierr);


	fprintf( vtk_fp, "\t\t\t<PointData Scalars=\" ");
	ierr = DMDAGetInfo( da, 0, 0,0,0, 0,0,0, &n_fields, 0, 0, 0,0,0 );
	//for( f=0; f<n_fields; f++ ) {
	//	const char *field_name;
	//	ierr = DMDAGetFieldName( da, f, &field_name );CHKERRQ(ierr);
	//	fprintf( vtk_fp, "%s ", field_name );
	//}
	fprintf( vtk_fp, "%s ", "Topography[km]" );
	fprintf( vtk_fp, "%s ", "SurfaceAmplitude" );

	fprintf( vtk_fp, "\">\n");


	ierr = DMGetLocalVector( da, &l_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin( da, FIELD, INSERT_VALUES, l_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd( da, FIELD, INSERT_VALUES, l_FIELD );CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, l_FIELD, &ZCOORD );CHKERRQ(ierr);


	{
		const char *field_name;
		ierr = DMDAGetFieldName( da, 0, &field_name );CHKERRQ(ierr);
		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", "Topography[km]" );
			for( j=sj; j<sj+ny; j++ ) {
				for( i=si; i<si+nx; i++ ) {
					double val;
					val = ZCOORD[j][i] 	* scaling_length;
					fprintf( vtk_fp, "\t\t\t\t\t%1.8f\n", val );
				}
			}
	}
	fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
	ierr = DMDAVecRestoreArray(da, l_FIELD, &ZCOORD );CHKERRQ(ierr);
	ierr = DMRestoreLocalVector( da, &l_FIELD );CHKERRQ(ierr);



	/* Write the amplitude of the surface deformation; usefull to look at folding */
	ierr = DMGetLocalVector( da, &l_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin( da, FIELD, INSERT_VALUES, l_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd( da, FIELD, INSERT_VALUES, l_FIELD );CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, l_FIELD, &ZCOORD );CHKERRQ(ierr);


	// Compute average height & subtract it
	VecNorm(FIELD,NORM_1, &AverageHeight);
	VecGetSize(FIELD,&size);
	AverageHeight = AverageHeight/((PetscScalar) size);
	{
		const char *field_name;
		ierr = DMDAGetFieldName( da, 0, &field_name );CHKERRQ(ierr);
		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", "SurfaceAmplitude" );
			for( j=sj; j<sj+ny; j++ ) {
				for( i=si; i<si+nx; i++ ) {
					double val;
					val = (ZCOORD[j][i] - AverageHeight)	* scaling_length;
					fprintf( vtk_fp, "\t\t\t\t\t%1.8f\n", val );
				}
			}

	}
	fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");

	ierr = DMDAVecRestoreArray(da, l_FIELD, &ZCOORD );CHKERRQ(ierr);
	ierr = DMRestoreLocalVector( da, &l_FIELD );CHKERRQ(ierr);


	fprintf( vtk_fp, "\t\t\t</PointData>\n");
	fprintf( vtk_fp, "\t\t</Piece>\n");
	fprintf( vtk_fp, "\t</StructuredGrid>\n");
	fprintf( vtk_fp, "</VTKFile>\n");


	if( vtk_fp!= NULL ) {
		fclose( vtk_fp );
		vtk_fp = NULL;
	}

	free(vtk_filename);

	PetscFunctionReturn(0);
}


//==========================================================================================================
// Create a 3D structured grid
#undef __FUNCT__
#define __FUNCT__ "DMView_3DVTK_StructuredGrid"
PetscErrorCode DMView_3DVTK_StructuredGrid( DM da, Vec FIELD, const char file_prefix[],  const char DirectoryName[],PetscScalar scaling_length )
{
	char 			*vtk_filename;
	PetscMPIInt 	rank;
	MPI_Comm 		comm;
	FILE 			*vtk_fp;
	PetscInt 		si,sj,sk, nx,ny,nz, i,j,k;
	PetscInt 		f,n_fields, N;
	DM 				cda;
	Vec 			coords;
	DMDACoor3d 		***_coords;
	Vec 			l_FIELD;
	PetscScalar 	*_L_FIELD;
	PetscErrorCode 	ierr;


	PetscFunctionBegin;

	// create file name
	PetscObjectGetComm( (PetscObject)da, &comm );
	MPI_Comm_rank( comm, &rank );
	asprintf( &vtk_filename, "./%s/subdomain-%s-p%1.6d.vts", DirectoryName, file_prefix, rank );

//	printf("DMView_3DVTK_StructuredGrid: file = %s, rank = %lld \n", vtk_filename,rank );


	// open file and write header
	vtk_fp = fopen( vtk_filename, "w" );
	if( vtk_fp == NULL ) {
		printf("ERROR(DMView_3DVTK_StructuredGrid): Cannot open file = %s \n", vtk_filename );
		exit(0);
	}

	fprintf( vtk_fp, "<?xml version=\"1.0\"?>\n");

	// coords
	ierr = DMDAGetGhostCorners( da, &si,&sj,&sk, &nx,&ny,&nz );CHKERRQ(ierr);
	N = nx * ny * nz;

#ifdef PETSC_WORDS_BIGENDIAN
	fprintf( vtk_fp, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	fprintf( vtk_fp, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
	fprintf( vtk_fp, "\t<StructuredGrid WholeExtent=\"%lld %lld %lld %lld %lld %lld\">\n", (LLD)si,(LLD)(si+nx-1), (LLD)sj,(LLD)(sj+ny-1), (LLD)sk,(LLD)(sk+nz-1));
	fprintf( vtk_fp, "\t\t<Piece Extent=\"%lld %lld %lld %lld %lld %lld\">\n", (LLD)si,(LLD)(si+nx-1), (LLD)sj,(LLD)(sj+ny-1), (LLD)sk,(LLD)(sk+nz-1));
	fprintf( vtk_fp, "\t\t\t<CellData></CellData>\n");
	fprintf( vtk_fp, "\t\t\t<Points>\n");

	// copy coordinates
	ierr = DMGetCoordinateDM( da, &cda);CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal( da,&coords );CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,coords,&_coords);CHKERRQ(ierr);

	fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
	for( k=sk; k<sk+nz; k++ ) {
		for( j=sj; j<sj+ny; j++ ) {
			for( i=si; i<si+nx; i++ ) {
				double x,y,z;

				x = _coords[k][j][i].x	* scaling_length;
				y = _coords[k][j][i].y	* scaling_length;
				z = _coords[k][j][i].z	* scaling_length;

				fprintf( vtk_fp, "\t\t\t\t\t%1.8f %1.8f %1.8f\n", x, y, z );
			}
		}
	}
	fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
	fprintf( vtk_fp, "\t\t\t</Points>\n");
	ierr = DMDAVecRestoreArray(cda,coords,&_coords);CHKERRQ(ierr);
	//ierr = VecDestroy(&coords);CHKERRQ(ierr);
	//ierr = DMDestroy(cda);CHKERRQ(ierr);

	fprintf( vtk_fp, "\t\t\t<PointData Scalars=\" ");
	ierr = DMDAGetInfo( da, 0, 0,0,0, 0,0,0, &n_fields, 0, 0, 0,0,0 );
	for( f=0; f<n_fields; f++ ) {
		const char *field_name;
		ierr = DMDAGetFieldName( da, f, &field_name );CHKERRQ(ierr);
		fprintf( vtk_fp, "%s ", field_name );
	}
	fprintf( vtk_fp, "\">\n");

	ierr = DMGetLocalVector( da, &l_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin( da, FIELD, INSERT_VALUES, l_FIELD );CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd( da, FIELD, INSERT_VALUES, l_FIELD );CHKERRQ(ierr);
	ierr = VecGetArray( l_FIELD, &_L_FIELD );CHKERRQ(ierr);
	for( f=0; f<n_fields; f++ ) {
		const char *field_name;

		ierr = DMDAGetFieldName( da, f, &field_name );CHKERRQ(ierr);
		fprintf( vtk_fp, "\t\t\t\t<DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", field_name );
		////
		for( i=0; i<N; i++ ) {
			double val;

			val = _L_FIELD[ n_fields*i + f ];
			fprintf( vtk_fp, "\t\t\t\t\t%1.8f\n", val );
		}
		////
		fprintf( vtk_fp, "\t\t\t\t</DataArray>\n");
	}
	ierr = VecRestoreArray( l_FIELD, &_L_FIELD );CHKERRQ(ierr);
	ierr = DMRestoreLocalVector( da, &l_FIELD );CHKERRQ(ierr);


	fprintf( vtk_fp, "\t\t\t</PointData>\n");
	fprintf( vtk_fp, "\t\t</Piece>\n");
	fprintf( vtk_fp, "\t</StructuredGrid>\n");
	fprintf( vtk_fp, "</VTKFile>\n");


	if( vtk_fp!= NULL ) {
		fclose( vtk_fp );
		vtk_fp = NULL;
	}

	// clean up
	free(vtk_filename);

	PetscFunctionReturn(0);
}

//==========================================================================================================
#undef __FUNCT__
#define __FUNCT__ "DMViewVTK_write_PieceExtend_Topo"
PetscErrorCode DMViewVTK_write_PieceExtend_Topo( FILE *vtk_fp, PetscInt indent_level, DM da, const char local_file_prefix[] )
{
	PetscMPIInt nproc,rank;
	MPI_Comm comm;
	const PetscInt *lx,*ly,*lz;
	PetscInt M,N,P,pM,pN,pP,sum,*olx,*oly,*olz;
	PetscInt *osx,*osy,*osz,*oex,*oey,*oez;
	PetscInt i,j,k,II,stencil;
	PetscErrorCode ierr;

	PetscFunctionBegin;
	/* create file name */
	PetscObjectGetComm( (PetscObject)da, &comm );
	MPI_Comm_size( comm, &nproc );
	MPI_Comm_rank( comm, &rank );

	ierr = DMDAGetInfo( da, 0, &M,&N,&P, &pM,&pN,&pP, 0, &stencil, 0, 0, 0, 0 );CHKERRQ(ierr);
	ierr = DMDAGetOwnershipRanges( da,&lx,&ly,&lz );CHKERRQ(ierr);

	/*generate start,end list */
	ierr = PetscMalloc(sizeof(PetscInt)*(size_t)(pM+1), &olx );CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(size_t)(pN+1), &oly );CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(size_t)(pP+1), &olz );CHKERRQ(ierr);
	sum = 0;
	for( i=0; i<pM; i++ ) {

		olx[i] = sum;

		sum = sum + lx[i];

	}

	olx[pM] = sum;
	sum = 0;
	for( i=0; i<pN; i++ ) {

		oly[i] = sum;

		sum = sum + ly[i];

	}

	oly[pN] = sum;
	sum = 0;
	for( i=0; i<pP; i++ ) {

		olz[i] = sum;

		sum = sum + lz[i];

	}

	olz[pP] = sum;
	ierr = PetscMalloc(sizeof(PetscInt)*(size_t)(pM), &osx );CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(size_t)(pN), &osy );CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(size_t)(pP), &osz );CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(size_t)(pM), &oex );CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(size_t)(pN), &oey );CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(size_t)(pP), &oez );CHKERRQ(ierr);
	for( i=0; i<pM; i++ ) {
		osx[i] = olx[i]   - stencil;
		oex[i] = olx[i] + lx[i] + stencil;
		if(osx[i]<0)osx[i]=0;
		if(oex[i]>M)oex[i]=M;
	}

	for( i=0; i<pN; i++ ) {
		osy[i] = oly[i]   - stencil;
		oey[i] = oly[i] + ly[i] + stencil;
		if(osy[i]<0)osy[i]=0;
		if(oey[i]>M)oey[i]=N;
	}
	for( i=0; i<pP; i++ ) {
		osz[i] = olz[i]   - stencil;
		oez[i] = olz[i] + lz[i] + stencil;
		if(osz[i]<0)osz[i]=0;
		if(oez[i]>P)oez[i]=P;
	}




	for( k=0;k<pP;k++ ) {
		for( j=0;j<pN;j++ ) {
			for( i=0;i<pM;i++ ) {
				char *name;
				PetscInt procid = i + j*pM + k*pM*pN; /* convert proc(i,j,k) to pid */
				asprintf( &name, "subdomain-%s-p%1.6lld.vts", local_file_prefix, (LLD)procid );
				for( II=0; II<indent_level; II++ ) {
					fprintf(vtk_fp,"\t");
				}
				fprintf( vtk_fp, "<Piece Extent=\"%lld %lld %lld %lld %lld %lld\"      Source=\"%s\"/>\n",
						(LLD)osx[i],(LLD)oex[i]-1LL,
						(LLD)osy[j],(LLD)oey[j]-1LL,
								0LL,0LL, name);
				free(name);
			}
		}
	}
	ierr = PetscFree( olx );CHKERRQ(ierr);
	ierr = PetscFree( oly );CHKERRQ(ierr);
	ierr = PetscFree( olz );CHKERRQ(ierr);
	ierr = PetscFree( osx );CHKERRQ(ierr);
	ierr = PetscFree( osy );CHKERRQ(ierr);
	ierr = PetscFree( osz );CHKERRQ(ierr);
	ierr = PetscFree( oex );CHKERRQ(ierr);
	ierr = PetscFree( oey );CHKERRQ(ierr);
	ierr = PetscFree( oez );CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

//==========================================================================================================
#undef __FUNCT__
#define __FUNCT__ "DMViewVTK_write_PieceExtend"
PetscErrorCode DMViewVTK_write_PieceExtend( FILE *vtk_fp, PetscInt indent_level, DM da, const char local_file_prefix[] )
{
	PetscMPIInt nproc,rank;
	MPI_Comm comm;
	const PetscInt *lx,*ly,*lz;
	PetscInt M,N,P,pM,pN,pP,sum,*olx,*oly,*olz;
	PetscInt *osx,*osy,*osz,*oex,*oey,*oez;
	PetscInt i,j,k,II,stencil;
	PetscErrorCode ierr;

	PetscFunctionBegin;
	/* create file name */
	PetscObjectGetComm( (PetscObject)da, &comm );
	MPI_Comm_size( comm, &nproc );
	MPI_Comm_rank( comm, &rank );

	ierr = DMDAGetInfo( da, 0, &M,&N,&P, &pM,&pN,&pP, 0, &stencil, 0, 0, 0, 0 );CHKERRQ(ierr);
	ierr = DMDAGetOwnershipRanges( da,&lx,&ly,&lz );CHKERRQ(ierr);

	/*generate start,end list */
	ierr = PetscMalloc(sizeof(PetscInt)*(size_t)(pM+1), &olx );CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(size_t)(pN+1), &oly );CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(size_t)(pP+1), &olz );CHKERRQ(ierr);
	sum = 0;
	for( i=0; i<pM; i++ ) {

		olx[i] = sum;

		sum = sum + lx[i];

	}

	olx[pM] = sum;
	sum = 0;
	for( i=0; i<pN; i++ ) {

		oly[i] = sum;

		sum = sum + ly[i];

	}

	oly[pN] = sum;
	sum = 0;
	for( i=0; i<pP; i++ ) {

		olz[i] = sum;

		sum = sum + lz[i];

	}

	olz[pP] = sum;
	ierr = PetscMalloc(sizeof(PetscInt)*(size_t)(pM), &osx );CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(size_t)(pN), &osy );CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(size_t)(pP), &osz );CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(size_t)(pM), &oex );CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(size_t)(pN), &oey );CHKERRQ(ierr);
	ierr = PetscMalloc(sizeof(PetscInt)*(size_t)(pP), &oez );CHKERRQ(ierr);
	for( i=0; i<pM; i++ ) {
		osx[i] = olx[i]   - stencil;
		oex[i] = olx[i] + lx[i] + stencil;
		if(osx[i]<0)osx[i]=0;
		if(oex[i]>M)oex[i]=M;
	}

	for( i=0; i<pN; i++ ) {
		osy[i] = oly[i]   - stencil;
		oey[i] = oly[i] + ly[i] + stencil;
		if(osy[i]<0)osy[i]=0;
		if(oey[i]>M)oey[i]=N;
	}
	for( i=0; i<pP; i++ ) {
		osz[i] = olz[i]   - stencil;
		oez[i] = olz[i] + lz[i] + stencil;
		if(osz[i]<0)osz[i]=0;
		if(oez[i]>P)oez[i]=P;
	}




	for( k=0;k<pP;k++ ) {
		for( j=0;j<pN;j++ ) {
			for( i=0;i<pM;i++ ) {
				char *name;
				PetscInt procid = i + j*pM + k*pM*pN; /* convert proc(i,j,k) to pid */
				asprintf( &name, "subdomain-%s-p%1.6lld.vts", local_file_prefix, (LLD)procid );
				for( II=0; II<indent_level; II++ ) {
					fprintf(vtk_fp,"\t");
				}
				fprintf( vtk_fp, "<Piece Extent=\"%lld %lld %lld %lld %lld %lld\"      Source=\"%s\"/>\n",
						(LLD)osx[i],(LLD)oex[i]-1LL,
						(LLD)osy[j],(LLD)oey[j]-1LL,
						(LLD)osz[k],(LLD)oez[k]-1LL, name);
				free(name);
			}
		}
	}
	ierr = PetscFree( olx );CHKERRQ(ierr);
	ierr = PetscFree( oly );CHKERRQ(ierr);
	ierr = PetscFree( olz );CHKERRQ(ierr);
	ierr = PetscFree( osx );CHKERRQ(ierr);
	ierr = PetscFree( osy );CHKERRQ(ierr);
	ierr = PetscFree( osz );CHKERRQ(ierr);
	ierr = PetscFree( oex );CHKERRQ(ierr);
	ierr = PetscFree( oey );CHKERRQ(ierr);
	ierr = PetscFree( oez );CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

//==========================================================================================================
#undef __FUNCT__
#define __FUNCT__ "DMView_2DVTK_PStructuredGrid"
PetscErrorCode DMView_2DVTK_PStructuredGrid( DM da, const char file_prefix[], const char local_file_prefix[], UserContext *user )
{
	MPI_Comm comm;
	PetscMPIInt nproc,rank;
	char *vtk_filename;
	FILE *vtk_fp;
	PetscInt M,N,P, si,sj,sk, nx,ny,nz, dim;
	PetscInt i,dofs,swidth;
	PetscErrorCode ierr;

	PetscFunctionBegin;

	/* only master generates this file */
	PetscObjectGetComm( (PetscObject)da, &comm );
	MPI_Comm_size( comm, &nproc );
	MPI_Comm_rank( comm, &rank );

	if( rank != 0 ) { PetscFunctionReturn(0); }

	/* create file name */
	asprintf( &vtk_filename, "%s.pvts", file_prefix );
	vtk_fp = fopen( vtk_filename, "w" );
	if( vtk_fp == NULL ) {
		printf("ERROR(DMView_2DVTK_PStructuredGrid): Cannot open file = %s \n", vtk_filename );
		exit(0);
	}
	free( vtk_filename );

	/* (VTK) generate pvts header */
	fprintf( vtk_fp, "<?xml version=\"1.0\"?>\n");

	fprintf( vtk_fp, "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");

	/* define size of the nodal mesh based on the cell DM */
	DMDAGetInfo( da, &dim, &M,&N,&P, 0,0,0, &dofs,&swidth,0,0, 0, 0 );
	DMDAGetGhostCorners( da, &si,&sj,&sk, &nx,&ny,&nz );
	if (swidth==1) {
		fprintf( vtk_fp, "\t<PStructuredGrid GhostLevel=\"1\" WholeExtent=\"%lld %lld %lld %lld %lld %lld\">\n", 0LL,(LLD)M-1LL, 0LL,(LLD)N-1LL, 0LL,0LL ); /* note overlap = 1 for Q1 */
	}
	if (swidth==2) {
		fprintf( vtk_fp, "\t<PStructuredGrid GhostLevel=\"2\" WholeExtent=\"%lld %lld %lld %lld %lld %lld\">\n", 0LL,(LLD)M-1LL, 0LL,(LLD)N-1LL, 0LL,0LL ); /* note overlap = 2 for Q2 */
	}


	/* DUMP THE CELL REFERENCES */
	fprintf( vtk_fp, "\t\t<PCellData>\n");
	fprintf( vtk_fp, "\t\t</PCellData>\n");


	/* DUMP THE NODAL REFERENCES */
	/*

	 <PPoints>
	 <PDataArray type="Float64" Name="Points" NumberOfComponents="3"/>
	 </PPoints>
	 <PPointData>
	 <PDataArray type="Float64" Name="vx" NumberOfComponents="1"/>
	 <PDataArray type="Float64" Name="vy" NumberOfComponents="1"/>
	 <PDataArray type="Float64" Name="vz" NumberOfComponents="1"/>
	 </PPointData>
	 <Piece Extent="0 40 0 4 0 20"      Source="step-0018-mesh-p0000.vts"/>

	 */


	fprintf( vtk_fp, "\t\t<PPoints>\n");
	fprintf( vtk_fp, "\t\t\t<PDataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\"/>\n");
	fprintf( vtk_fp, "\t\t</PPoints>\n");

	if (dofs!=3){
		fprintf( vtk_fp, "\t\t<PPointData>\n");
		if (dim==3){
			for(i=0;i<dofs;i++){
				const char *fieldname;
				ierr = DMDAGetFieldName(da,i,&fieldname);CHKERRQ(ierr);
				fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"1\"/>\n",fieldname);

			}

			// assuming we deal with a surface topography file:
			fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"1\"/>\n","SurfaceVelocity_Vx[]");
			fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"1\"/>\n","SurfaceVelocity_Vy[]");
			fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"1\"/>\n","SurfaceVelocity_Vz[]");
			fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"1\"/>\n","SurfaceAmplitude");
		}
		else if (dim==2){
			// erosion surface topography
			fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"1\"/>\n","Topography[km]");
			fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"1\"/>\n","SurfaceAmplitude");

		}

		fprintf( vtk_fp, "\t\t</PPointData>\n");
	}
	else{
		fprintf( vtk_fp, "\t\t<PPointData>\n");
		if (user->DimensionalUnits==1){
			fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Velocity [cm/yr]\" NumberOfComponents=\"3\"/>\n");
		}
		else{
			fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Velocity [ ]\" NumberOfComponents=\"3\"/>\n");
		}
		fprintf( vtk_fp, "\t\t</PPointData>\n");

	}

	/* write out the parallel information */
	if (dim==3){
		ierr = DMViewVTK_write_PieceExtend_Topo(vtk_fp,2,da,local_file_prefix);CHKERRQ(ierr);
	}
	else if ((dim==2) && (rank==0)){
		char *name;
		PetscInt II;

//		ierr = DMViewVTK_write_PieceExtend_Topo(vtk_fp,2,da,local_file_prefix);CHKERRQ(ierr);

		asprintf( &name, "subdomain-%s-p%1.6lld.vts", local_file_prefix, (LLD)0 );
		for( II=0; II<2; II++ ) {
			fprintf(vtk_fp,"\t");
		}
		fprintf( vtk_fp, "<Piece Extent=\"%lld %lld %lld %lld %lld %lld\"      Source=\"%s\"/>\n",
				(LLD)0,(LLD)M-1,(LLD)0,(LLD)N-1,(LLD)0,(LLD)0, name);

	}
	/* close the file */
	fprintf( vtk_fp, "\t</PStructuredGrid>\n");
	fprintf( vtk_fp, "</VTKFile>\n");

	if(vtk_fp!=NULL){
		fclose( vtk_fp );
		vtk_fp = NULL;
	}

	PetscFunctionReturn(0);
}

//==========================================================================================================
#undef __FUNCT__
#define __FUNCT__ "DMView_3DVTK_PStructuredGrid"
PetscErrorCode DMView_3DVTK_PStructuredGrid( DM da, const char file_prefix[], const char local_file_prefix[], UserContext *user )
{
	MPI_Comm comm;
	PetscMPIInt nproc,rank;
	char *vtk_filename;
	FILE *vtk_fp;
	PetscInt M,N,P, si,sj,sk, nx,ny,nz;
	PetscInt i,dofs,swidth;
	PetscErrorCode ierr;

	PetscFunctionBegin;

	/* only master generates this file */
	PetscObjectGetComm( (PetscObject)da, &comm );
	MPI_Comm_size( comm, &nproc );
	MPI_Comm_rank( comm, &rank );

	if( rank != 0 ) { PetscFunctionReturn(0); }

	/* create file name */
	asprintf( &vtk_filename, "%s.pvts", file_prefix );
	vtk_fp = fopen( vtk_filename, "w" );
	if( vtk_fp == NULL ) {
		printf("ERROR(DMView_3DVTK_PStructuredGrid): Cannot open file = %s \n", vtk_filename );
		exit(0);
	}
	free( vtk_filename );

	/* (VTK) generate pvts header */
	fprintf( vtk_fp, "<?xml version=\"1.0\"?>\n");

#ifdef PETSC_WORDS_BIGENDIAN
	fprintf( vtk_fp, "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	fprintf( vtk_fp, "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif

	/* define size of the nodal mesh based on the cell DM */
	DMDAGetInfo( da, 0, &M,&N,&P, 0,0,0, &dofs,&swidth,0,0, 0, 0 );
	DMDAGetGhostCorners( da, &si,&sj,&sk, &nx,&ny,&nz );
	if (swidth==1) {
		fprintf( vtk_fp, "\t<PStructuredGrid GhostLevel=\"1\" WholeExtent=\"%lld %lld %lld %lld %lld %lld\">\n", 0LL,M-1LL, 0LL,N-1LL, 0LL,P-1LL ); /* note overlap = 1 for Q1 */
	}
	if (swidth==2) {
		fprintf( vtk_fp, "\t<PStructuredGrid GhostLevel=\"2\" WholeExtent=\"%lld %lld %lld %lld %lld %lld\">\n", 0LL,M-1LL, 0LL,N-1LL, 0LL,P-1LL ); /* note overlap = 2 for Q2 */
	}


	/* DUMP THE CELL REFERENCES */
	fprintf( vtk_fp, "\t\t<PCellData>\n");
	fprintf( vtk_fp, "\t\t</PCellData>\n");


	/* DUMP THE NODAL REFERENCES */
	/*

	 <PPoints>
	 <PDataArray type="Float64" Name="Points" NumberOfComponents="3"/>
	 </PPoints>
	 <PPointData>
	 <PDataArray type="Float64" Name="vx" NumberOfComponents="1"/>
	 <PDataArray type="Float64" Name="vy" NumberOfComponents="1"/>
	 <PDataArray type="Float64" Name="vz" NumberOfComponents="1"/>
	 </PPointData>
	 <Piece Extent="0 40 0 4 0 20"      Source="step-0018-mesh-p0000.vts"/>

	 */


	fprintf( vtk_fp, "\t\t<PPoints>\n");
	fprintf( vtk_fp, "\t\t\t<PDataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\"/>\n");
	fprintf( vtk_fp, "\t\t</PPoints>\n");

	if (dofs!=3){
		fprintf( vtk_fp, "\t\t<PPointData>\n");
		for(i=0;i<dofs;i++){
			const char *fieldname;
			ierr = DMDAGetFieldName(da,i,&fieldname);CHKERRQ(ierr);
			fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"1\"/>\n",fieldname);
		}
		fprintf( vtk_fp, "\t\t</PPointData>\n");
	}
	else{
		fprintf( vtk_fp, "\t\t<PPointData>\n");
		if (user->DimensionalUnits==1){
			fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Velocity [cm/yr]\" NumberOfComponents=\"3\"/>\n");
		}
		else{
			fprintf(vtk_fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Velocity [ ]\" NumberOfComponents=\"3\"/>\n");
		}
		fprintf( vtk_fp, "\t\t</PPointData>\n");

	}
	/* write out the parallel information */
	ierr = DMViewVTK_write_PieceExtend(vtk_fp,2,da,local_file_prefix);CHKERRQ(ierr);


	/* close the file */
	fprintf( vtk_fp, "\t</PStructuredGrid>\n");
	fprintf( vtk_fp, "</VTKFile>\n");

	if(vtk_fp!=NULL){
		fclose( vtk_fp );
		vtk_fp = NULL;
	}

	PetscFunctionReturn(0);
}

//==========================================================================================================
/*
 Generates the following

 On each processor you will get a file like
 subdomain-[NAME]-mesh-p0000.vts
 subdomain-[NAME]-mesh-p0001.vts
 ...

 For the set of all files, you will get a parallel metdata file
 [NAME]-mesh.pvts
 This will reference the files
 subdomain-[NAME]-mesh-p0000.vts
 subdomain-[NAME]-mesh-p0001.vts
 ...

 This routine is specifically designed to write topography files
 */
#undef __FUNCT__
#define __FUNCT__ "DMView2DPVTS_Topo"
PetscErrorCode DMView2DPVTS_Topo(DM da, Vec field, const char NAME[], UserContext *user, const char DirectoryName[],PetscScalar scaling_length)
{
	char 			*vts_filename, *pvts_filename;
	PetscInt 		n_fields, dim;
	PetscErrorCode 	ierr;


	PetscFunctionBegin;

	asprintf( &vts_filename, "%s-mesh", NAME );
//#ifdef __VTK_ASCII__
	//ierr = DMView_3DVTK_StructuredGrid(da,field,vts_filename);CHKERRQ(ierr);

	// get the number of fields
	ierr = DMDAGetInfo( da, &dim, 0,0,0, 0,0,0, &n_fields, 0, 0, 0, 0, 0 );CHKERRQ(ierr);

	/* Only writes topography */
	if (dim==3){
		/* We are dealing with the surface and bottom topographies which are distributed over the code */
		ierr = DMView_2DVTK_StructuredGrid_Topo(da,field,vts_filename,DirectoryName, user->SurfaceTopography_Vx, user->SurfaceTopography_Vy, user->SurfaceTopography_Vz,scaling_length);CHKERRQ(ierr); // other fields

		asprintf( &pvts_filename, "./%s/%s-mesh", DirectoryName, NAME );
		ierr = DMView_2DVTK_PStructuredGrid( da, pvts_filename, vts_filename, user );CHKERRQ(ierr);
		free(pvts_filename);
	}
	else if (dim==2){
		/* We are most likely dealing with the surface topography that is used for the erosion code.
		 * This is present on a 2D DA, and only on rank 0
		 */

		ierr = DMView_2DVTK_StructuredGrid_Topo_Erosion(da,field,vts_filename,DirectoryName, scaling_length);CHKERRQ(ierr); // other fields

		asprintf( &pvts_filename, "./%s/%s-mesh", DirectoryName, NAME );
		ierr = DMView_2DVTK_PStructuredGrid( da, pvts_filename, vts_filename, user );CHKERRQ(ierr);
		free(pvts_filename);

	}

		//#endif
		//#ifndef __VTK_ASCII__
		//	ierr = DAView_3DVTK_StructuredGrid_appended(da_nodes,velocity,temperature,vts_filename);CHKERRQ(ierr);
		//#endif
		// Binary files will be approx half the size, but will load into paraview in around 10 x times faster.
		// Binary output files will be generated 10 times faster than the ascii output files.






	// clean up
	free(vts_filename);

	PetscFunctionReturn(0);
}

//==========================================================================================================
#undef __FUNCT__
#define __FUNCT__ "DMView3DPVTS"
PetscErrorCode DMView3DPVTS(DM da, Vec field, const char NAME[], UserContext *user, const char DirectoryName[],PetscScalar scaling_length)
{
	char 			*vts_filename;
	char 			*pvts_filename;
	PetscInt 		n_fields;
	PetscErrorCode 	ierr;

	PetscFunctionBegin;

	asprintf( &vts_filename, "%s-mesh", NAME );
//#ifdef __VTK_ASCII__
	//ierr = DMView_3DVTK_StructuredGrid(da,field,vts_filename);CHKERRQ(ierr);

	// get the number of fields
	ierr = DMDAGetInfo( da, 0, 0,0,0, 0,0,0, &n_fields, 0, 0, 0,0,0 );

	// If we have 3 fields within LaMEM it is most likely the velocity field, otherwise all fields are treated as scalar
	if (n_fields==3){
		if (user->DimensionalUnits==1){
			ierr = DMView_3DVTK_StructuredGrid_3ComponentVector(da,field,vts_filename,DirectoryName, "Velocity [cm/yr]",scaling_length); CHKERRQ(ierr); // velocity
		}
		else{
			ierr = DMView_3DVTK_StructuredGrid_3ComponentVector(da,field,vts_filename,DirectoryName, "Velocity [  ]",scaling_length); CHKERRQ(ierr); // velocity
		}
	}
	else{
		ierr = DMView_3DVTK_StructuredGrid(da,field,vts_filename,DirectoryName,scaling_length);CHKERRQ(ierr); // other fields
	}


//#endif
//#ifndef __VTK_ASCII__
//	ierr = DAView_3DVTK_StructuredGrid_appended(da_nodes,velocity,temperature,vts_filename);CHKERRQ(ierr);
//#endif
	// Binary files will be approx half the size, but will load into paraview in around 10 x times faster.
	// Binary output files will be generated 10 times faster than the ascii output files.

	asprintf( &pvts_filename, "./%s/%s-mesh", DirectoryName, NAME );
	ierr = DMView_3DVTK_PStructuredGrid( da, pvts_filename, vts_filename, user );CHKERRQ(ierr);

	// clean up
	free(pvts_filename);
	free(vts_filename);

	PetscFunctionReturn(0);
}

//==========================================================================================================
#undef __FUNCT__
#define __FUNCT__ "ParaviewPVDOpen"
PetscErrorCode ParaviewPVDOpen(const char pvdfilename[])
{
	PetscMPIInt rank;
	FILE *fp;

	PetscFunctionBegin;
	/* only master generates this file */
	MPI_Comm_rank( PETSC_COMM_WORLD, &rank );
	if( rank != 0 ) { PetscFunctionReturn(0); }

	fp = fopen(pvdfilename,"w");
	fprintf(fp,"<?xml version=\"1.0\"?>\n");
#ifdef PETSC_WORDS_BIGENDIAN
	fprintf(fp,"<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
	fprintf(fp,"<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif

	fprintf(fp,"<Collection>\n");

	fprintf(fp,"</Collection>\n");
	fprintf(fp,"</VTKFile>\n");
	fclose(fp);
	PetscFunctionReturn(0);
}

//==========================================================================================================
#undef __FUNCT__
#define __FUNCT__ "ParaviewPVDAppend"
PetscErrorCode ParaviewPVDAppend(const char pvdfilename[],double time,const char datafile[], const char DirectoryName[])
{
	PetscMPIInt rank;
	FILE *fp;
	char line[10000];
	PetscInt key_L;
	char key[] = "</Collection>";
	char *copy,*tmp;

	PetscFunctionBegin;
	/* only master generates this file */
	MPI_Comm_rank( PETSC_COMM_WORLD, &rank );
	if( rank != 0 ) { PetscFunctionReturn(0); }


	fp = fopen(pvdfilename,"r");
	/* reset to start of file */
	rewind(fp);

	copy = NULL;
	key_L = (PetscInt) strlen( key );
	while( !feof(fp) ) {
		fgets( line, 10000-1, fp );
		if ( strncmp(key,line, (size_t)key_L)!=0 ) {

			/* copy line */
			if (copy!=NULL) {
			  asprintf(&tmp,"%s",copy);
			  free(copy);
                          asprintf(&copy,"%s%s",tmp,line);
                          free(tmp);
			}
			else {
			  asprintf(&copy,"%s",line);
			}
		}
		else {
			break;
		}
	}
	fclose(fp);

	/* open new file - clobbering the old */
	fp = fopen(pvdfilename,"w");

	/* write all copied chars */
	fprintf(fp,"%s",copy);

	/* write new data */
	fprintf(fp,"  <DataSet timestep=\"%1.6e\" file=\"./%s/%s\"/>\n",time, DirectoryName, datafile );

	/* close tag */
	fprintf(fp,"</Collection>\n");
	fprintf(fp,"</VTKFile>\n");

	fclose(fp);
	free(copy);

	PetscFunctionReturn(0);
}
//==========================================================================================================
// Write a velocity output file and perform scaling/unscaling of the variables
#undef __FUNCT__
#define __FUNCT__ "WriteTopographyOutputFile_VTS"
PetscErrorCode WriteTopographyOutputFile_VTS(UserContext *user,PetscInt itime, const char DirectoryName[])
{
	char 				*mfield_name, *pfield_name;
	PetscScalar 		scaling_length;
	PetscErrorCode 		ierr;

	PetscFunctionBegin;

	// --- Part 1: surface topography ---

	if (user->Output.surface_topography==1) {

	// create filename
	asprintf(&pfield_name, "%s_Surface_topography_timeseries.pvd",user->OutputFile);

	// create pvd file for topography
	if (itime == 0) {
		ierr = ParaviewPVDOpen(pfield_name);CHKERRQ(ierr);
	}

	// scale topography and lengths in km
	if (user->DimensionalUnits==1){

		scaling_length = user->Characteristic.Length/1000.0;
	}
	else{
		scaling_length = 1.0;
	}

	// create datafile
	asprintf(&mfield_name, "%s_Surface_topography_step%1.6lld",user->OutputFile, (LLD)itime );
	ierr = DMView2DPVTS_Topo(user->DA_SurfaceTopography, user->SurfaceTopography, mfield_name, user, DirectoryName,scaling_length);CHKERRQ(ierr);
	free(mfield_name);


	// add datafile to the time series
	asprintf(&mfield_name, "%s_Surface_topography_step%1.6lld-mesh.pvts", user->OutputFile, (LLD)itime );
	if (user->DimensionalUnits==1){
		ierr = ParaviewPVDAppend(pfield_name,user->time * user->Characteristic.Time/user->Characteristic.SecYear/1e6, mfield_name, DirectoryName );CHKERRQ(ierr);		// in Myrs
	}
	else{
		ierr = ParaviewPVDAppend(pfield_name,user->time , mfield_name, DirectoryName );CHKERRQ(ierr);
	}

	// Print info
	PetscPrintf(PETSC_COMM_WORLD,"# Saved Paraview VTS Surface Topography output in file %s \n",mfield_name); CHKERRQ(ierr);

	// clean up
	free(mfield_name);
	free(pfield_name);

	}

#if 1		// deactivated by default as it is typically flat (or does not change during the simulation)


	// --- Part 2: Bottom topography ---

	if (user->Output.bottom_topography==1) {

	// create filename
	asprintf(&pfield_name, "%s_Bottom_topography_timeseries.pvd",user->OutputFile);

	// create pvd file for topography
	if (itime == 0) {
		ierr = ParaviewPVDOpen(pfield_name);CHKERRQ(ierr);
	}

	// create datafile
	asprintf(&mfield_name, "%s_Bottom_topography_step%1.6lld",user->OutputFile, (LLD)itime );
	ierr = DMView2DPVTS_Topo(user->DA_BottomTopography, user->BottomTopography, mfield_name, user, DirectoryName,scaling_length);CHKERRQ(ierr);
	free(mfield_name);

	// add datafile to the time series
	asprintf(&mfield_name, "%s_Bottom_topography_step%1.6lld-mesh.pvts", user->OutputFile, (LLD)itime );
	if (user->DimensionalUnits==1){
		ierr = ParaviewPVDAppend(pfield_name,user->time * user->Characteristic.Time/user->Characteristic.SecYear/1e6, mfield_name, DirectoryName );CHKERRQ(ierr);		// in Myrs
	}
	else{
		ierr = ParaviewPVDAppend(pfield_name,user->time , mfield_name, DirectoryName );CHKERRQ(ierr);
	}

	// Print info
	PetscPrintf(PETSC_COMM_WORLD,"# Saved Paraview VTS Bottom Topography output in file %s \n",mfield_name); CHKERRQ(ierr);

	// clean up
	free(mfield_name);
	free(pfield_name);

	}

#endif


	if (user->ErosionParameters.ErosionModel==2){
		PetscMPIInt			rank;

		MPI_Comm_rank(PETSC_COMM_WORLD,&rank);


		if (rank==0){  // serial erosion code, so surface only present on rank 0

			// --- Part 3 - Write the internal erosion surface if requested --
			// create filename
			asprintf(&pfield_name, "%s_Erosion_topography_timeseries.pvd",user->OutputFile);

			// create pvd file for topography
			if (itime == 0) {
				ierr = ParaviewPVDOpen(pfield_name);CHKERRQ(ierr);
			}

			// create datafile
			asprintf(&mfield_name, "%s_Erosion_topography_step%1.6lld",user->OutputFile, (LLD)itime );

			ierr = DMView2DPVTS_Topo(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode, user->ErosionParameters.FE_ErosionCode.ErosionSurface, mfield_name, user, DirectoryName,scaling_length);CHKERRQ(ierr);

			free(mfield_name);

			// add datafile to the time series
			asprintf(&mfield_name, "%s_Erosion_topography_step%1.6lld-mesh.pvts", user->OutputFile, (LLD)itime );
			if (user->DimensionalUnits==1){
				ierr = ParaviewPVDAppend(pfield_name,user->time * user->Characteristic.Time/user->Characteristic.SecYear/1e6, mfield_name, DirectoryName );CHKERRQ(ierr);		// in Myrs
			}
			else{
				ierr = ParaviewPVDAppend(pfield_name,user->time , mfield_name, DirectoryName );CHKERRQ(ierr);
			}

			// Print info
			PetscPrintf(PETSC_COMM_WORLD,"# Saved Paraview VTS Erosion Topography output in file %s \n",mfield_name); CHKERRQ(ierr);

			// clean up
			free(mfield_name);
			free(pfield_name);
		}

	}





	PetscFunctionReturn(0);
}

//==========================================================================================================
// Write a velocity output file and perform scaling/unscaling of the variables
#undef __FUNCT__
#define __FUNCT__ "WriteVelocityOutputFile_VTS"
PetscErrorCode WriteVelocityOutputFile_VTS(UserContext *user,PetscInt itime, Vec sol, const char DirectoryName[])
{
	char 				*mfield_name, *pfield_name;
	PetscErrorCode 		ierr;
	PetscScalar 		scaling_length;


	PetscFunctionBegin;

	// create filename
	asprintf(&pfield_name, "%s_velocity_timeseries.pvd",user->OutputFile);

	// create pvd file for velocity
	if (itime == 0) {
		ierr = ParaviewPVDOpen(pfield_name);CHKERRQ(ierr);
	}

	// scale vel and create a scaling factor for writing coordinates in km
	if (user->DimensionalUnits==1){
		ierr = VecScale(sol,user->Characteristic.Velocity*user->Characteristic.cmYear);CHKERRQ(ierr);	// in cm/yr

		scaling_length = user->Characteristic.Length/1000.0;
	}
	else{
		scaling_length = 1.0;
	}


	// create datafile
	asprintf(&mfield_name, "%s_velocity_step%1.6lld",user->OutputFile, (LLD)itime );
	ierr = DMView3DPVTS(user->DA_Vel, sol, mfield_name, user, DirectoryName,scaling_length);CHKERRQ(ierr);
	free(mfield_name);


	// unscale velocity
	if (user->DimensionalUnits==1){
		ierr = VecScale(sol,1.0/(user->Characteristic.Velocity*user->Characteristic.cmYear));CHKERRQ(ierr);
	}

	// add datafile to the time series
	asprintf(&mfield_name, "%s_velocity_step%1.6lld-mesh.pvts", user->OutputFile, (LLD)itime );
	if (user->DimensionalUnits==1){
		ierr = ParaviewPVDAppend(pfield_name,user->time * user->Characteristic.Time/user->Characteristic.SecYear/1e6, mfield_name, DirectoryName );CHKERRQ(ierr);		// in Myrs
	}
	else{
		ierr = ParaviewPVDAppend(pfield_name,user->time , mfield_name, DirectoryName );CHKERRQ(ierr);
	}

	// print info
	PetscPrintf(PETSC_COMM_WORLD,"# Saved Paraview VTS Velocity output in file %s \n",mfield_name); CHKERRQ(ierr);


	// clean up
	free(mfield_name);
	free(pfield_name);


	PetscFunctionReturn(0);
}


//==========================================================================================================
// Write a temperature output file and perform scaling/unscaling of the variables
#undef __FUNCT__
#define __FUNCT__ "WriteTemperatureOutputFile_VTS"
PetscErrorCode WriteTemperatureOutputFile_VTS(UserContext *user,PetscInt itime, Vec Temp, const char DirectoryName[])
{
	char 			*mfield_name, *pfield_name;
	PetscScalar		scaling_length;
	PetscErrorCode 	ierr;

	PetscFunctionBegin;

	// create pvd file for temperature
	asprintf(&pfield_name, "%s_temperature_timeseries.pvd",user->OutputFile);
	if (itime == 0) {
		ierr = ParaviewPVDOpen(pfield_name);CHKERRQ(ierr);
	}

	// scale temperature and and create a scaling factor for writing coordinates in km
	if (user->DimensionalUnits==1){
		ierr = VecScale(Temp,user->Characteristic.Temperature);CHKERRQ(ierr);

		scaling_length = user->Characteristic.Length/1000.0;
	}
	else{
		scaling_length = 1.0;
	}


	// write data
	asprintf(&mfield_name, "%s_temperature_step%1.6lld", user->OutputFile, (LLD)itime );
	ierr = DMView3DPVTS(user->DA_Temp, Temp,mfield_name, user, DirectoryName,scaling_length);CHKERRQ(ierr);

	// get rid of memory leak
	free(mfield_name);

	// unscale temperature
	if (user->DimensionalUnits==1){
		ierr = VecScale(Temp,1.0/user->Characteristic.Temperature);		CHKERRQ(ierr);
	}

	// add datafile to the time series
	asprintf(&mfield_name, "%s_temperature_step%1.6lld-mesh.pvts", user->OutputFile, (LLD)itime );
	if (user->DimensionalUnits==1){
		ierr = ParaviewPVDAppend(pfield_name,user->time * user->Characteristic.Time/user->Characteristic.SecYear/1e6, mfield_name, DirectoryName );CHKERRQ(ierr);		// in Myrs
	}
	else{
		ierr = ParaviewPVDAppend(pfield_name,user->time , mfield_name, DirectoryName );CHKERRQ(ierr);
	}


	// print info
	PetscPrintf(PETSC_COMM_WORLD,"# Saved Paraview VTS Temperature output in file %s \n",mfield_name); CHKERRQ(ierr);

	// clean up
	free(mfield_name);
	free(pfield_name);

	PetscFunctionReturn(0);
}

//==========================================================================================================
// Create a 3D Voronoi diagram from particles with phase information, write the file to disk and perform scaling/unscaling of the variables
#undef __FUNCT__
#define __FUNCT__ "WritePhasesOutputFile_VTS"
PetscErrorCode WritePhasesOutputFile_VTS(LaMEMVelPressureDA C, UserContext *user,PetscInt itime, const char DirectoryName[])
{
	char *mfield_name, *pfield_name;
	float scaling;
	PetscErrorCode ierr;

	PetscFunctionBegin;

	asprintf(&pfield_name, "%s_phases_timeseries.pvd",user->OutputFile);
	if (itime == 0) {
		/* create pvd file for phases */
		ierr = ParaviewPVDOpen(pfield_name);CHKERRQ(ierr);
	}

	/* Are we computing in dimensional or dimensionless units? */
	scaling = 1.0;
	if (user->DimensionalUnits==1){
		scaling = (float)(user->Characteristic.Length/1e3);	// in km
	}

	/* write data */
	asprintf(&mfield_name, "%s_phases_step%1.6lld", user->OutputFile, (LLD)itime );
	ierr = AVDPhaseViewerExecute(C,user->DA_Vel,user->DA_Processors,user->num_particle_local,
			user->ParticlesLocal, mfield_name, DirectoryName, scaling);CHKERRQ(ierr);

	// get rid of memory leak
	free(mfield_name);

	/* add datafile to the time series */
	asprintf(&mfield_name, "%s_phases_step%1.6lld.pvtr", user->OutputFile, (LLD)itime );

	if (user->DimensionalUnits==1){  // time in  Myrs
		ierr = ParaviewPVDAppend(pfield_name,user->time * user->Characteristic.Time/user->Characteristic.SecYear/1e6, mfield_name, DirectoryName );CHKERRQ(ierr);		// in Myrs
	}
	else{
		ierr = ParaviewPVDAppend(pfield_name,user->time , mfield_name, DirectoryName );CHKERRQ(ierr);
	}

	/* Print info*/
	PetscPrintf(PETSC_COMM_WORLD,"# Saved Paraview VTS phases output in file %s \n",mfield_name); CHKERRQ(ierr);


	free(mfield_name);
	free(pfield_name);

	PetscFunctionReturn(0);
}

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

 							LaMEM_Particles.c

 Contains particle-related routines that are used in LaMEM

 $Id$

 Routines:

 		InitializeParticles				-	Initializes particle distribution
 		ElementInfo						-	Information such as min,max -coordinates of a given element
 		Inside_8NodeElement				-	Are we inside an 8-node (linear) element or not?
 		NaturalCoordinates				-	Get the natural coordinates of a particle in an element
 		GetParticleNodes				-	Find in which 8-node element a particle is, move particles to different CPUs
 											 if required, delete them if they're outside the computational domain and
 											 inject new particles if required.
 		AdvectParticles					-	Advect particles forward in time
 		InterpolateParticles			-	Interpolate particles after reading in modified mesh from file
 		SetInitialTracerProperties		-	Set initial tracer properties by hard-coding.
 		SetInitialTracerPhasesFromFile	-	Read particle-phases file from disk and set particle phases and temperature
 											accordingly.
 		MaterialPropertiesFromTracers	-	Compute material properties @ integration points from particles.
 		WriteParticlesToDisc			-	Write particle distributions to disk.
 		RemeshGrid						-	Remesh the computational grid (for use with ALE mode)
		ParticlePhaseTransitions		-	Perform phase transitions on particles
 		ComputePropertiesAtParticles	-	Computes material properties such as stress tensor at particles

This version is compatible with PETSc version 3.3-p2
The code cannot be distributed, used for commercial purposes or handed on,
without the explicit agreement of Boris Kaus.

 */

#include "LaMEM.h"
#include "LaMEM_Particles.h"
#include "Utils.h"
#include "Elements.h"
#include "Material.h"
#include "Mesh.h"
#include "LaMEM_Temperature.h"
#include "Quadrature.h"
#include "StokesOperators.h"


#ifdef PARTICLES
/* ==================================================================================================== */
/* Initialize particles on the various processors, based on the local element geometry (which works consistently
 * for low order elements as well) */
//DONE
#undef __FUNCT__
#define __FUNCT__ "InitializeParticles_ElementWise"
PetscErrorCode InitializeParticles_ElementWise( LaMEMVelPressureDA C, DM da, UserContext *user )
{
	PetscMPIInt     rank;
	PetscErrorCode  ierr;
	Particles 		ParticleLocal;
	PetscInt		nx,ny,nz,xs,ys,zs,xm,ym,zm,  num_particle_local, num_particle_start;
	PetscInt		xsp,ysp,zsp, xmp,ymp,zmp, iel_x, iel_y, iel_z, i, j, k, nel_x, nel_y, nel_z;
	PetscInt		iix, iiy, iiz, ip_x, ip_y, ip_z, ip_num, nnel;
	PetscInt		n_particle_props, AddRandomNoiseParticles;
	PetscScalar		dx_natural, dy_natural, dz_natural, Point[3];
	PetscRandom		rctx;
	PetscScalar		cf_rand;
	DAVPElementType element_type;
	DMDACoor3d		 ***coords, coord_elem[MAX_nnel], CoordPoint;
	DM			 	cda;
	Vec			 	gc;
	PetscScalar	 	ShapeVel[MAX_nnel],  **dhdsVel;

	element_type = C->type;
	nnel 		 = C->nnel;


	/* Do we add random noise or not ? */
	AddRandomNoiseParticles = 1;
	ierr = PetscOptionsGetInt(PETSC_NULL,"-AddRandomNoiseParticles",&AddRandomNoiseParticles	, PETSC_NULL); CHKERRQ(ierr);
	if (AddRandomNoiseParticles==0){
		PetscPrintf( PETSC_COMM_WORLD, "# Not adding random noise to initial particle distribution. \n");
	}
	else {
		PetscPrintf( PETSC_COMM_WORLD, "# Adding random noise to initial particle distribution. \n");
		AddRandomNoiseParticles = 1;
	}

	/* Check the number of particle properties is consistent */
	n_particle_props = sizeof(Particles)/sizeof(PetscScalar);
	if( n_particle_props != particle_props ) {
		PetscPrintf( PETSC_COMM_WORLD, "*** It appears that the number of particle properties is incorrect ***\n");
		PetscPrintf( PETSC_COMM_WORLD, "*** Please check that all members of typedef struct { } Particles; are of type PetscScalar in LaMEM.h ***\n");
		PetscPrintf( PETSC_COMM_WORLD, "*** Please check that #define particle_props in LaMEM.h matches the number of members in typedef struct { } Particles; ***\n");
		EMERGENCY_EXIT("Number of particle properties is inconsistent");
	}


	ierr = DMDAGetInfo(da,0,&nx,&ny,&nz,0,0,0,0,0,0,0,0,0);	 	CHKERRQ(ierr);  	// # of nodes in all directions
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); 		CHKERRQ(ierr);

	ierr = DMDAGetCorners(user->DA_Processors,&xsp,&ysp,&zsp,&xmp,&ymp,&zmp);		CHKERRQ(ierr);
	ierr = DMDAGetInfo(user->DA_Processors,0,&nel_x,&nel_y,&nel_z,0,0,0,0,0,0,0,0,0);	CHKERRQ(ierr);

	nx = user->finest_nnode_x;
	ny = user->finest_nnode_y;
	nz = user->finest_nnode_z;

	// correct for periodicity
	if ( (user->BC.LeftBound==3) && ((xs+xm)==user->finest_nnode_x-1) ) {
		xm = 	xm+1;
	}
	if ( (user->BC.FrontBound==3) && ((ys+ym)==user->finest_nnode_y-1) ) {
		ym =	ym+1;
	}


	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	// Initialize the Random number generator
	ierr = PetscRandomCreate(PETSC_COMM_WORLD, &rctx); CHKERRQ(ierr);
	ierr = PetscRandomSetFromOptions(rctx); CHKERRQ(ierr);

	/* Allocate the local particles array (which has a slight buffer) */
	user->MaxNumLocalParticles = ((PetscInt)(particle_buffer*(xmp*user->NumPartX)*(ymp*user->NumPartY)*(zmp*user->NumPartZ)));


	ierr = PetscMalloc( (size_t)user->MaxNumLocalParticles*sizeof(Particles),  	&user->ParticlesLocal	); CHKERRQ(ierr);
	ierr = PetscMemzero(user->ParticlesLocal,(size_t)user->MaxNumLocalParticles*sizeof(Particles)); CHKERRQ(ierr);

	num_particle_local = 0;
	num_particle_start = 0;
	if (zsp>0){	num_particle_start = num_particle_start +	(zsp-1)*user->NumPartZ*(nel_y*user->NumPartY)*(nel_x*user->NumPartX);	}
	if (ysp>0){	num_particle_start = num_particle_start +	(ysp-1)*user->NumPartY*(nel_x*user->NumPartX) ;    					}
	if (xsp>0){	num_particle_start = num_particle_start +	(xsp-1)*user->NumPartX*(nel_x); 				 					}

	ierr = LaMEMCreate2dArray( 3,nnel, &dhdsVel, PETSC_NULL ); CHKERRQ(ierr);

	/* Extract coordinate arrays */
	ierr = DMGetCoordinateDM(da,&cda); CHKERRQ(ierr);	//coordinates
	ierr = DMGetCoordinatesLocal(da,&gc); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,gc,&coords); CHKERRQ(ierr);

	/* Make a loop over all local elements */
	for(iel_x=xsp; iel_x<xsp+xmp; iel_x++){
		for (iel_y=ysp; iel_y<ysp+ymp; iel_y++){
			for (iel_z=zsp; iel_z<zsp+zmp; iel_z++){

				if( (element_type==DAVP_Q1P0) || (element_type==DAVP_Q1Q1) || (element_type==DAVP_FDSTAG)) {
					i = iel_x;
					j = iel_y;
					k = iel_z;
				}
				else if( (element_type==DAVP_Q2PM1L) || (element_type==DAVP_Q2PM1G) ){
					i = 2*iel_x;
					j = 2*iel_y;
					k = 2*iel_z;
				}
				else {
					EMERGENCY_EXIT("Unknown element type");
				}

				/* Extract coordinates of the local element in correct order */
				ierr = GetElementCoords(coord_elem, coords, i,j,k, 1); 					CHKERRQ(ierr);
				ierr = CorrectElementCoordsForPeriodicity( coord_elem, i, j, user,1 ); 	CHKERRQ(ierr);

				/* Add particles to current element */
				dx_natural = 2.0/((PetscScalar) (user->NumPartX));
				dy_natural = 2.0/((PetscScalar) (user->NumPartY));
				dz_natural = 2.0/((PetscScalar) (user->NumPartZ));
				for (iix=0; iix<user->NumPartX; iix++){
					for (iiy=0; iiy<user->NumPartY; iiy++){
						for (iiz=0; iiz<user->NumPartZ; iiz++){

							/* Set local coordinates of the particle */
							ParticleLocal.eta     = ((PetscScalar)(iix))*dx_natural - 1.0  +  dx_natural/2;
							ParticleLocal.zetha   = ((PetscScalar)(iiy))*dy_natural - 1.0  +  dy_natural/2;
							ParticleLocal.phi     = ((PetscScalar)(iiz))*dz_natural - 1.0  +  dz_natural/2;

							/* Add random noise if required 		*/
							if (AddRandomNoiseParticles==1	){

								ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
								cf_rand = (cf_rand-0.5);
								ParticleLocal.eta	=	ParticleLocal.eta 	+ cf_rand*dx_natural;

								ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
								cf_rand = (cf_rand-0.5);
								ParticleLocal.zetha	=	ParticleLocal.zetha + cf_rand*dy_natural;

								ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
								cf_rand = (cf_rand-0.5);
								ParticleLocal.phi	=	ParticleLocal.phi 	+ cf_rand*dz_natural;

							}

							/* Compute global coordinates of the particle */
							Point[0] = ParticleLocal.eta;
							Point[1] = ParticleLocal.zetha;
							Point[2] = ParticleLocal.phi;

							ComputeShapeFunctionVelocity(ShapeVel, dhdsVel, Point);				 		// Velocity shape function
							ComputeCoordIntp(nnel,ShapeVel, coord_elem,  &CoordPoint);				 	// Real coordinate of point

							/* Set coordinates to the particle */
							ParticleLocal.x = CoordPoint.x;
							ParticleLocal.y = CoordPoint.y;
							ParticleLocal.z = CoordPoint.z;

							/* helps finding the particle later - WILL BECOME OBSOLETE AT SOME POINT AS WE SHOULD SWITCH TO ELEMENT-BASED NUMBERING*/
							ParticleLocal.ix    = ((PetscScalar)(i));
							ParticleLocal.iy    = ((PetscScalar)(j));
							ParticleLocal.iz    = ((PetscScalar)(k));

							// initialize all variables to zero - important on some machines!
							ParticleLocal.phase			= 0.0;
							ParticleLocal.Plastic		= 0.0;
							ParticleLocal.Mu_eff		= 0.0;
							ParticleLocal.P				= 0.0;
							ParticleLocal.PlasticStrain = 0.0;
							ParticleLocal.Strain		= 0.0;
							ParticleLocal.Txx			= 0.0;
							ParticleLocal.Tyy			= 0.0;
							ParticleLocal.Tzz			= 0.0;
							ParticleLocal.Txy			= 0.0;
							ParticleLocal.Txz			= 0.0;
							ParticleLocal.Tyz			= 0.0;
							ParticleLocal.E2nd			= 0.0;
							ParticleLocal.T2nd			= 0.0;
							ParticleLocal.E2nd_pl		= 0.0;
							ParticleLocal.Mu_viscous	= 0.0;
							ParticleLocal.T				= 0.0;

							// Give every particle a unique number
							ip_x 				= 	iel_x*user->NumPartX + iix;
							ip_y 				= 	iel_y*user->NumPartY + iiy;
							ip_z 				= 	iel_z*user->NumPartZ + iiz;
							ip_num              =   ip_z*( ((nel_x)*user->NumPartX)*((nel_y)*user->NumPartY) ) + ip_y*((nel_x)*user->NumPartX) + ip_x;
							ParticleLocal.num   = 	((PetscScalar)(ip_num	));

							ParticleLocal.cpu   = ((PetscScalar)(rank 	));		// the cpu of the current particle

							// Add to particle array
							user->ParticlesLocal[num_particle_local] = ParticleLocal;

							num_particle_local = num_particle_local+1;
						}
					}
				}


			}
		}
	}
	user->num_particle_local = num_particle_local;

	/* Free coordinate arrays */
	ierr = DMDAVecRestoreArray(cda,gc,&coords); 						CHKERRQ(ierr);
	//	ierr = DMDestroy(cda); 											CHKERRQ(ierr);	//coordinates
	//	ierr = VecDestroy(gc); 											CHKERRQ(ierr);
	ierr = LaMEMDestroy2dArray(&dhdsVel, PETSC_NULL ); 	CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD,"Initialized particles \n");

	/* Cleaning up */
	ierr = PetscRandomDestroy(&rctx); CHKERRQ(ierr);	// Destroy random context

	PetscFunctionReturn(0);
}
/* ==================================================================================================== */




/* ==================================================================================================== */
/* Initialize particles on the various processors */
//DONE
#undef __FUNCT__
#define __FUNCT__ "InitializeParticles"
PetscErrorCode InitializeParticles( DM da, UserContext *user )
{
	PetscMPIInt     rank;
	PetscErrorCode  ierr;
	Particles 		ParticleLocal;
	PetscInt		nx,ny,nz,xs,ys,zs,xm,ym,zm,  num_particle_local, num_particle_start;
	PetscInt		ix,iy,iz, iix, iiy, iiz, ip_x, ip_y, ip_z, ip_num;
	PetscInt		n_particle_props, AddRandomNoiseParticles;
	PetscScalar		dx,dy,dz;
	PetscRandom		rctx;
	PetscScalar	    cf_rand;

	/* Do we add random noise or not ? */
	AddRandomNoiseParticles = 1;
	ierr = PetscOptionsGetInt(PETSC_NULL,"-AddRandomNoiseParticles",&AddRandomNoiseParticles	, PETSC_NULL); CHKERRQ(ierr);
	if (AddRandomNoiseParticles==0){
		PetscPrintf( PETSC_COMM_WORLD, "# Not adding random noise to initial particle distribution. \n");
	}
	else {
		PetscPrintf( PETSC_COMM_WORLD, "# Adding random noise to initial particle distribution. \n");
		AddRandomNoiseParticles = 1;
	}




	/* Check the number of particle properties is consistent */
	n_particle_props = sizeof(Particles)/sizeof(PetscScalar);
	if( n_particle_props != particle_props ) {
		PetscPrintf( PETSC_COMM_WORLD, "*** It appears that the number of particle properties is incorrect ***\n");
		PetscPrintf( PETSC_COMM_WORLD, "*** Please check that all members of typedef struct { } Particles; are of type PetscScalar in LaMEM.h ***\n");
		PetscPrintf( PETSC_COMM_WORLD, "*** Please check that #define particle_props in LaMEM.h matches the number of members in typedef struct { } Particles; ***\n");
		EMERGENCY_EXIT("Number of particle properties is inconsistent");
	}


	ierr = DMDAGetInfo(da,0,&nx,&ny,&nz,0,0,0,0,0,0,0,0,0);	 CHKERRQ(ierr);  	// # of nodes in all directions
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);

	nx = user->finest_nnode_x;
	ny = user->finest_nnode_y;
	nz = user->finest_nnode_z;

	// correct for periodicity
	if ( (user->BC.LeftBound==3) && ((xs+xm)==user->finest_nnode_x-1) ) {
		xm = 	xm+1;
	}
	if ( (user->BC.FrontBound==3) && ((ys+ym)==user->finest_nnode_y-1) ) {
		ym =	ym+1;
	}


	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	// Initialize the Random number generator
	ierr = PetscRandomCreate(PETSC_COMM_WORLD, &rctx); CHKERRQ(ierr);
	ierr = PetscRandomSetFromOptions(rctx); CHKERRQ(ierr);


	dx = user->W/((double) ((nx-1)*user->NumPartX));
	dy = user->L/((double) ((ny-1)*user->NumPartY));
	dz = user->H/((double) ((nz-1)*user->NumPartZ));

	/* Allocate the local particles array (which has a slight buffer) */
	user->MaxNumLocalParticles = ((PetscInt)(particle_buffer*(xm*user->NumPartX)*(ym*user->NumPartY)*(zm*user->NumPartZ)));
	ierr = PetscMalloc( (size_t)user->MaxNumLocalParticles*sizeof(Particles),  	&user->ParticlesLocal	); CHKERRQ(ierr);
	ierr = PetscMemzero(user->ParticlesLocal,(size_t)user->MaxNumLocalParticles*sizeof(Particles)); CHKERRQ(ierr);


	num_particle_local = 0;
	num_particle_start = 0;
	if (zs>0){	num_particle_start = num_particle_start +	(zs-1)*user->NumPartZ*((ny-1)*user->NumPartY)*((nx-1)*user->NumPartX);	}
	if (ys>0){	num_particle_start = num_particle_start +	(ys-1)*user->NumPartY*((nx-1)*user->NumPartX) ;    					}
	if (xs>0){	num_particle_start = num_particle_start +	(xs-1)*user->NumPartX*(nx-1); 				 					}

	for (ix=xs;ix<PetscMin(xs+xm, nx-1); ix++){
		for (iy=ys;iy<PetscMin(ys+ym, ny-1); iy++){
			for (iz=zs;iz<PetscMin(zs+zm, nz-1); iz++){

				for (iix=0; iix<user->NumPartX; iix++){
					for (iiy=0; iiy<user->NumPartY; iiy++){
						for (iiz=0; iiz<user->NumPartZ; iiz++){
							ParticleLocal.x     = ((PetscScalar)(ix*user->NumPartX + iix))*dx + user->x_left  +  0.5*dx;
							ParticleLocal.y     = ((PetscScalar)(iy*user->NumPartY + iiy))*dy + user->y_front +  0.5*dy;
							ParticleLocal.z     = ((PetscScalar)(iz*user->NumPartZ + iiz))*dz + user->z_bot   +  0.5*dz;

							ParticleLocal.ix    = ((PetscScalar)(ix));
							ParticleLocal.iy    = ((PetscScalar)(iy));
							ParticleLocal.iz    = ((PetscScalar)(iz));

							ParticleLocal.phase   = 0.0;
							ParticleLocal.Plastic = 0.0;
							ParticleLocal.Mu_eff  = 0.0;
							ParticleLocal.eta			= 0.0;
							ParticleLocal.zetha			= 0.0;
							ParticleLocal.phi			= 0.0;
							ParticleLocal.P				= 0.0;
							ParticleLocal.PlasticStrain = 0.0;
							ParticleLocal.Strain		= 0.0;
							ParticleLocal.Txx			= 0.0;
							ParticleLocal.Tyy			= 0.0;
							ParticleLocal.Tzz			= 0.0;
							ParticleLocal.Txy			= 0.0;
							ParticleLocal.Txz			= 0.0;
							ParticleLocal.Tyz			= 0.0;
							ParticleLocal.E2nd			= 0.0;
							ParticleLocal.T2nd			= 0.0;
							ParticleLocal.E2nd_pl		= 0.0;
							ParticleLocal.Mu_viscous	= 0.0;
							ParticleLocal.T				= 0.0;

							/* Give every particle a unique number */
							ip_x 				= 	ix*user->NumPartX + iix;
							ip_y 				= 	iy*user->NumPartY + iiy;
							ip_z 				= 	iz*user->NumPartZ + iiz;
							ip_num              =   ip_z*( ((nx-1)*user->NumPartX)*((ny-1)*user->NumPartY) ) + ip_y*((nx-1)*user->NumPartX) + ip_x;
							ParticleLocal.num   = 	((PetscScalar)(ip_num	));

							ParticleLocal.cpu   = ((PetscScalar)(rank 	));		// the cpu of the current particle

							// Add random noise if required
							if (AddRandomNoiseParticles==1	){
								// Add some random component to the particle
								ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
								cf_rand=cf_rand-0.5;
								ParticleLocal.x		=	ParticleLocal.x + cf_rand*dx;
								ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
								cf_rand=cf_rand-0.5;
								ParticleLocal.y		=	ParticleLocal.y + cf_rand*dy;
								ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
								cf_rand=cf_rand-0.5;
								ParticleLocal.z		=	ParticleLocal.z + cf_rand*dz;
							}


							// Add to particle array
							user->ParticlesLocal[num_particle_local] = ParticleLocal;

							num_particle_local = num_particle_local+1;

						}
					}
				}


			}
		}
	}
	user->num_particle_local = num_particle_local;

	PetscPrintf(PETSC_COMM_WORLD,"Initialized particles \n");

	/* Cleaning up */
	ierr = PetscRandomDestroy(&rctx); CHKERRQ(ierr);	// Destroy random context

	PetscFunctionReturn(0);
}
/* ==================================================================================================== */


/* ==================================================================================================== */
/* compute useful information about the element, such as minimum and maximum x-coordinates,
 *  maximum x-coordinate on left side etc
 */
#undef __FUNCT__
#define __FUNCT__ "ElementInfo"
PetscErrorCode ElementInfo( DMDACoor3d coord_elem[], PetscScalar *minx, PetscScalar *maxx, PetscScalar *miny,
		PetscScalar *maxy, PetscScalar *minz, PetscScalar *maxz, PetscScalar *max_leftx, PetscScalar *min_rightx,
		PetscScalar *max_fronty, PetscScalar *min_backy, PetscScalar *max_botz, PetscScalar *min_topz )
{
	PetscInt	inode;
	PetscScalar mix, max, miy,may, miz, maz, ma_leftx, mi_rightx;
	PetscScalar	ma_fronty, mi_backy, ma_botz, mi_topz;

	mix=1e8; max=-1e8; miy=1e8; may=-1e8; miz=1e8; maz=-1e8;
	for (inode=0; inode<8; inode++){
		mix = PetscMin(mix,coord_elem[inode].x);	max = PetscMax(max,coord_elem[inode].x);
		miy = PetscMin(miy,coord_elem[inode].y);	may = PetscMax(may,coord_elem[inode].y);
		miz = PetscMin(miz,coord_elem[inode].z);	maz = PetscMax(maz,coord_elem[inode].z);
	}
	*minx = mix;	*maxx = max;  *miny = miy; *maxy=may; *minz=miz; *maxz=maz;

	// Maximum x-coordinate on left side of element
	ma_leftx 	= PetscMax(coord_elem[0].x	, coord_elem[3].x);  ma_leftx 	= PetscMax(ma_leftx	, coord_elem[4].x);
	ma_leftx 	= PetscMax(ma_leftx			, coord_elem[7].x);	 *max_leftx = ma_leftx;

	// Minimum x-coordinate on right side of element
	mi_rightx 	= PetscMin(coord_elem[1].x	,coord_elem[2].x);   mi_rightx 	= PetscMin(mi_rightx,coord_elem[5].x);
	mi_rightx 	= PetscMin(mi_rightx		,coord_elem[6].x);   *min_rightx= mi_rightx;

	// Max. y-coordinate in front of element
	ma_fronty 	= PetscMax(coord_elem[0].y	, coord_elem[1].y);  ma_fronty 	= PetscMax(ma_fronty, coord_elem[4].y);
	ma_fronty 	= PetscMax(ma_fronty		, coord_elem[5].y);  *max_fronty= ma_fronty;

	// Min y-coord in back of element
	mi_backy 	= PetscMin(coord_elem[2].y	,coord_elem[3].y);   mi_backy 	= PetscMin(mi_backy	,coord_elem[6].y);
	mi_backy 	= PetscMin(mi_backy			,coord_elem[7].y);   *min_backy	= mi_backy;

	// Max. z-coordinate @ bottom of element
	ma_botz 	= PetscMax(coord_elem[0].z	, coord_elem[1].z);  ma_botz 	= PetscMax(ma_botz, coord_elem[2].z);
	ma_botz 	= PetscMax(ma_botz			, coord_elem[3].z);  *max_botz	= ma_botz;

	// Min. z-coordinate @ top of element
	mi_topz 	= PetscMin(coord_elem[5].z	, coord_elem[4].z);  mi_topz 	= PetscMin(mi_topz, coord_elem[6].z);
	mi_topz 	= PetscMin(mi_topz			, coord_elem[7].z);  *min_topz	= mi_topz;

	PetscFunctionReturn(0);
}
/* ==================================================================================================== */


/* ==================================================================================================== */
/* ---------------------------------------------------------------------
 * Finds the local coordinates of a particle inside a 3D, 8-node element
 * 	Algorithm is based on an idea of Thomas Kocher, and uses Newton Rhapson iterations
 * 		that converge rapidly if the particle is inside but diverge if the particle is outside.
 */
#undef __FUNCT__
#define __FUNCT__ "Inside_8NodeElement"
PetscErrorCode  Inside_8NodeElement( PetscScalar x_real, PetscScalar y_real, PetscScalar z_real,
		PetscScalar *eta_out, PetscScalar *zetha_out, PetscScalar *phi_out, DMDACoor3d coord_elem[],	PetscInt *inside, PetscInt DisplayIterations)
{
	PetscScalar 	a, b, c, d, e,f,g,h,i, coord[8][3], error_level;
	PetscScalar		x, y, z, rhs[3], dhdsVel[3][8], ShapeVel[8];
	PetscScalar		Point_real[3], diff[3], criterion;
	PetscInt		iter, stopped_iterations, jj, ii, kk;
	PetscScalar		Jacob[3][3],  InvJacob[3][3],	DeltaCoord[3], CriticalValue;

	/* preparations */
	for(ii=0;ii<8;ii++){
		coord[ii][0] = coord_elem[ii].x; coord[ii][1] = coord_elem[ii].y; coord[ii][2] = coord_elem[ii].z;
	}

	if (DisplayIterations==1){
		PetscPrintf(PETSC_COMM_SELF, "Inside_8NodeElement: coord_elem_x=[%f,%f,%f,%f,%f,%f,%f,%f] \n",coord[0][0],coord[1][0],coord[2][0],coord[3][0],
				coord[4][0],coord[5][0],coord[6][0],coord[7][0]);

		PetscPrintf(PETSC_COMM_SELF, "Inside_8NodeElement: coord_elem_y=[%f,%f,%f,%f,%f,%f,%f,%f] \n",coord[0][1],coord[1][1],coord[2][1],coord[3][1],
				coord[4][1],coord[5][1],coord[6][1],coord[7][1]);

		PetscPrintf(PETSC_COMM_SELF, "Inside_8NodeElement: coord_elem_z=[%f,%f,%f,%f,%f,%f,%f,%f] \n",coord[0][2],coord[1][2],coord[2][2],coord[3][2],
				coord[4][2],coord[5][2],coord[6][2],coord[7][2]);

	}

	/* Start with some initial guess of the natural coordinate of the point */
	x		=	*eta_out;
	y		=	*zetha_out;
	z		=	*phi_out;

	rhs[0]			=	x_real; rhs[1]		  =	y_real;	rhs[2]		  =	z_real;

	criterion           = 1;
	stopped_iterations  = 0;
	iter                = 0;
	error_level         = 1e-10;
	while ( criterion > error_level){

		/* Shape function and derivative @ the natural coordinate -----------------------*/
		ShapeVel[0]		=	0.125*(1.0-x)*(1.0-y)*(1.0-z);  ShapeVel[1]=0.125*(1.0+x)*(1.0-y)*(1.0-z);
		ShapeVel[2]		=	0.125*(1.0+x)*(1.0+y)*(1.0-z);  ShapeVel[3]=0.125*(1.0-x)*(1.0+y)*(1.0-z);
		ShapeVel[4]		=	0.125*(1.0-x)*(1.0-y)*(1.0+z);  ShapeVel[5]=0.125*(1.0+x)*(1.0-y)*(1.0+z);
		ShapeVel[6]		=	0.125*(1.0+x)*(1.0+y)*(1.0+z);  ShapeVel[7]=0.125*(1.0-x)*(1.0+y)*(1.0+z);

		dhdsVel[0][0]	=  -0.125*(1-y)*(1-z);	dhdsVel[0][1]=	0.125*(1-y)*(1-z);	dhdsVel[0][2]= 0.125*(1+y)*(1-z);
		dhdsVel[0][3]	=  -0.125*(1+y)*(1-z);	dhdsVel[0][4]= -0.125*(1-y)*(1+z);	dhdsVel[0][5]= 0.125*(1-y)*(1+z);
		dhdsVel[0][6]	=	0.125*(1+y)*(1+z);	dhdsVel[0][7]= -0.125*(1+y)*(1+z);
		dhdsVel[1][0]	=  -0.125*(1-x)*(1-z);	dhdsVel[1][1]= -0.125*(1+x)*(1-z);	dhdsVel[1][2]= 0.125*(1+x)*(1-z);
		dhdsVel[1][3]	=   0.125*(1-x)*(1-z);	dhdsVel[1][4]= -0.125*(1-x)*(1+z);	dhdsVel[1][5]=-0.125*(1+x)*(1+z);
		dhdsVel[1][6]	=   0.125*(1+x)*(1+z);	dhdsVel[1][7]=  0.125*(1-x)*(1+z);
		dhdsVel[2][0]	=  -0.125*(1-x)*(1-y);	dhdsVel[2][1]= -0.125*(1+x)*(1-y);	dhdsVel[2][2]=-0.125*(1+x)*(1+y);
		dhdsVel[2][3]	=  -0.125*(1-x)*(1+y); 	dhdsVel[2][4]=	0.125*(1-x)*(1-y); 	dhdsVel[2][5]= 0.125*(1+x)*(1-y);
		dhdsVel[2][6]	=	0.125*(1+x)*(1+y);	dhdsVel[2][7]=	0.125*(1-x)*(1+y);
		/*------------------------------------------------------------------------------*/

		/* Compute Jacobian and it's inverse -------------------------------------------*/
		for( ii = 0; ii<3;ii++){ for(jj=0;jj<3;jj++){ Jacob[ii][jj] = 0; }}
		for( ii = 0; ii < 3; ii++){
			for( jj = 0; jj < 3; jj++){
				for( kk = 0; kk < 8; kk++){
					Jacob[ii][jj] +=  dhdsVel[ii][kk]*coord[kk][jj];
				}
			}
		}
		a = Jacob[0][0]; b = Jacob[0][1]; c = Jacob[0][2]; d = Jacob[1][0];
		e = Jacob[1][1]; f = Jacob[1][2]; g = Jacob[2][0]; h = Jacob[2][1]; i = Jacob[2][2];
		InvJacob[0][0]	=	-(e*i-f*h)/(-a*e*i+a*f*h+d*b*i-d*c*h-g*b*f+g*c*e);
		InvJacob[0][1]	=	 (b*i-c*h)/(-a*e*i+a*f*h+d*b*i-d*c*h-g*b*f+g*c*e);
		InvJacob[0][2]	=	-(b*f-c*e)/(-a*e*i+a*f*h+d*b*i-d*c*h-g*b*f+g*c*e);
		InvJacob[1][0]	=	 (d*i-f*g)/(-a*e*i+a*f*h+d*b*i-d*c*h-g*b*f+g*c*e);
		InvJacob[1][1]	=	-(a*i-c*g)/(-a*e*i+a*f*h+d*b*i-d*c*h-g*b*f+g*c*e);
		InvJacob[1][2]	=	 (a*f-c*d)/(-a*e*i+a*f*h+d*b*i-d*c*h-g*b*f+g*c*e);
		InvJacob[2][0]	=	-(d*h-e*g)/(-a*e*i+a*f*h+d*b*i-d*c*h-g*b*f+g*c*e);
		InvJacob[2][1]	=	-(-a*h+b*g)/(-a*e*i+a*f*h+d*b*i-d*c*h-g*b*f+g*c*e);
		InvJacob[2][2]	=	 (-a*e+b*d)/(-a*e*i+a*f*h+d*b*i-d*c*h-g*b*f+g*c*e);
		/*------------------------------------------------------------------------------*/

		if (DisplayIterations==1){
			PetscPrintf(PETSC_COMM_SELF, "Inside_8NodeElement: iter=%lld, ShapeVel=[%f,%f,%f,%f,%f,%f,%f,%f] \n",(LLD)iter, ShapeVel[0], ShapeVel[1],ShapeVel[2],ShapeVel[3],ShapeVel[4],ShapeVel[5],ShapeVel[6],ShapeVel[7]);
		}

		/* Compute real coordinate for this natural coordinate */
		Point_real[0]=0; Point_real[1]=0; Point_real[2]=0;
		for (ii=0; ii<8; ii++){
			Point_real[0] = Point_real[0] + ShapeVel[ii]*coord_elem[ii].x;
			Point_real[1] = Point_real[1] + ShapeVel[ii]*coord_elem[ii].y;
			Point_real[2] = Point_real[2] + ShapeVel[ii]*coord_elem[ii].z;
		}

		for (ii=0; ii<3; ii++){diff[ii]	= rhs[ii]-	Point_real[ii];}							/* Compute difference */

		/* Compute error between natural and real point */
		criterion           = 	diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2];
		iter				=   iter + 1;

		if (criterion>error_level){
			/* Compute increment in coordinates from the jacobian	*/
			for (ii=0; ii<3; ii++){	DeltaCoord[ii]	=	0;	}
			for (ii=0; ii<3; ii++){
				for (jj=0; jj<3; jj++){
					DeltaCoord[ii] = DeltaCoord[ii] + InvJacob[ii][jj]*diff[jj];
				}
			}


			/* Update local coordinates */
			x = x + 0.2*DeltaCoord[0];
			y = y + 0.2*DeltaCoord[1];
			z = z + 0.2*DeltaCoord[2];


		}
		if (DisplayIterations==1){
			PetscPrintf(PETSC_COMM_SELF, "Inside_8NodeElement: Iteration %lld, error=%g [x,y,z]=[%f,%f,%f] x_real=[%f,%f,%f] \n",(LLD)iter, criterion,x,y,z,x_real, y_real, z_real);
		}


		if (iter>100 || PetscAbs(z)>5 || PetscAbs(x)>5	){
			criterion          = 0;
			stopped_iterations = 1;
		}

	}

	*eta_out   = x;
	*zetha_out = y;
	*phi_out   = z;

	*inside = 0;

	CriticalValue = 1.0 + 1e-8;  // if we're a tiny bit outside its also fine
	if ( (PetscAbs(x)<CriticalValue) && (PetscAbs(y)<CriticalValue) && (PetscAbs(z)<CriticalValue) && (criterion<error_level) && (stopped_iterations=0)){
		*inside = 1;
	}

	//MPI_Abort(PETSC_COMM_WORLD,1);
	PetscFunctionReturn(0);
}
/* ==================================================================================================== */


/* ---------------------------------------------------------------------
 * Finds the local coordinates of a particle inside a 3D, Q2P1 or Q1P0 element
 */
#undef __FUNCT__
#define __FUNCT__ "NaturalCoordinates"
PetscErrorCode  NaturalCoordinates( const PetscInt nnel, PetscScalar x_real, PetscScalar y_real, PetscScalar z_real,
		PetscScalar *eta_out, PetscScalar *zetha_out, PetscScalar *phi_out, DMDACoor3d coord_elem[], PetscInt *inside )
{
	double 			error_level;
	PetscScalar		x, y, z, rhs[3], **dhdsVel, ShapeVel[MAX_nnel];
	PetscScalar		Point_real[3], diff[3], criterion;
	PetscInt		iter, jj, ii;
	PetscScalar		**Jacob,  **InvJacob,	DeltaCoord[3], CriticalValue, Point[3], DetJacob;

	LaMEMCreate2dArray( 3,nnel, &dhdsVel, PETSC_NULL );
	LaMEMCreate2dArray( 3,3, &Jacob, PETSC_NULL );
	LaMEMCreate2dArray( 3,3, &InvJacob, PETSC_NULL );

	/* Start with some initial guess of the natural coordinate of the point */
	x		=	*eta_out;
	y		=	*zetha_out;
	z		=	*phi_out;

	rhs[0]			=	x_real; rhs[1]		  =	y_real;	rhs[2]		  =	z_real;

	criterion           = 1;
	iter                = 0;
	error_level         = 1e-10;
	while ( criterion > error_level){

		/* Shape function and derivative @ the natural coordinate -----------------------*/
		Point[0] = x; Point[1]=y; Point[2]=z;
		ComputeShapeFunctionVelocity(ShapeVel, dhdsVel, Point);
		/*------------------------------------------------------------------------------*/

		/* Compute Jacobian and it's inverse -------------------------------------------*/
		ComputeJacobianElement(nnel,dhdsVel,coord_elem,Jacob,InvJacob,&DetJacob);  // Jacobian etc. of element
		/*------------------------------------------------------------------------------*/


		if (DetJacob<0){
			//		PetscPrintf(PETSC_COMM_WORLD,"In routine NaturalCoordinates: Negative Jacobian, DetJacob=%g  \n",DetJacob);
			//			PetscPrintf(PETSC_COMM_WORLD," x-coords of local element: [%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g] \n",coord_elem[0].x,coord_elem[1].x,coord_elem[2].x,coord_elem[3].x,coord_elem[4].x,coord_elem[5].x,coord_elem[6].x,coord_elem[7].x,
			//					coord_elem[8].x,coord_elem[9].x,coord_elem[10].x,coord_elem[11].x,coord_elem[12].x,coord_elem[13].x,coord_elem[14].x,coord_elem[15].x);

			//		MPI_Abort(PETSC_COMM_WORLD,1);
		}

		/* Compute real coordinate for this natural coordinate */
		Point_real[0]=0; Point_real[1]=0; Point_real[2]=0;
		for (ii=0; ii<nnel; ii++){
			Point_real[0] = Point_real[0] + ShapeVel[ii]*coord_elem[ii].x;
			Point_real[1] = Point_real[1] + ShapeVel[ii]*coord_elem[ii].y;
			Point_real[2] = Point_real[2] + ShapeVel[ii]*coord_elem[ii].z;
		}

		for (ii=0; ii<3; ii++){diff[ii]	= rhs[ii]-	Point_real[ii];}							/* Compute difference */

		/* Compute error between natural and real point */
		criterion           = 	diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2];
		iter				=   iter + 1;

		if (criterion>error_level){
			/* Compute increment in coordinates from the jacobian	*/
			for (ii=0; ii<3; ii++){	DeltaCoord[ii]	=	0.0;	}
			for (ii=0; ii<3; ii++){
				for (jj=0; jj<3; jj++){
					DeltaCoord[ii] = DeltaCoord[ii] + InvJacob[ii][jj]*diff[jj];
				}
			}

			/* Update local coordinates */
			x = x + 0.2*DeltaCoord[0];
			y = y + 0.2*DeltaCoord[1];
			z = z + 0.2*DeltaCoord[2];
		}


		if (iter>100 || PetscAbs(z)>5 || PetscAbs(x)>5	|| PetscAbs(y)>5 ){
			criterion          = 0;
		}

	}

	*eta_out   = x;
	*zetha_out = y;
	*phi_out   = z;

	*inside = 0;
	CriticalValue = 1.0 + 1e-3;  // if we're a tiny bit outside its also fine
	if ( (PetscAbs(x)<CriticalValue) && (PetscAbs(y)<CriticalValue) && (PetscAbs(z)<CriticalValue) && (criterion<error_level) ){
		*inside = 1;
	}


	LaMEMDestroy2dArray(&dhdsVel, PETSC_NULL );
	LaMEMDestroy2dArray(&Jacob, PETSC_NULL );
	LaMEMDestroy2dArray(&InvJacob, PETSC_NULL );

	PetscFunctionReturn(0);

}
/*end of GetParticleNodes_old
 *
 *  ==================================================================================================== */




/* ==================================================================================================== */
/* Finds which the surrounding nodes are for given particles.
  for particles we only work with the surrounding 8 nodes

  The routine consists in the following parts:

  (1) Loop through the particles and determine whether particles that used to be close to a boundary are still on the current CPU
	  (thereby taking into account periodicity if required). If not, add them to a list of particles to be send over.
  (2) Send all particles to the correct CPU
  (3) Loop through all particles and determine the local coordinates of the particle.

 */
#undef __FUNCT__
#define __FUNCT__ "GetParticleNodes"
PetscErrorCode GetParticleNodes( LaMEMVelPressureDA C, DM da, UserContext *user)
{
	PetscMPIInt     rank, size;
	PetscErrorCode  ierr;
	PetscInt 		ipart, ix, iy, iz, inside, i, j, k, num, AddRandomNoiseParticles;
	PetscInt		TracersOutsideComputationalDomain, DeactivatedTracers, TracersSendToOtherCPU;
	PetscInt		ElementFound, TracerDeactivate,	InjectionPhase;
	PetscInt		nx,ny,nz, num_particle_local_new, num_backwards;
	PetscInt		xs,ys,zs,xm,ym,zm,xs_g, ys_g, zs_g, xm_g, ym_g, zm_g;
	PetscInt		iproc, cpu, NumParticlesMove[MaxNumCPU], LengthParticleCache;
	PetscInt		*NumberParticleElement, num_local_part, TotalNumberOfParticles;
	PetscInt		iix,iiy,iiz;
	PetscScalar		eta, zetha, phi, x,y,z, ShapeVel[8], dx_natural, dy_natural, dz_natural;
	Particles 		ParticleLocal, **ParticlesSendAway, *ParticlesLost, ParticleLost;
	DM				cda;
	DMDACoor3d		 ***coords, coord_elem[8];
	Vec			 	gc;
	PetscLogDouble	cputime_start, cputime_end, cputime_start0, cputime_end0;
	MPI_Datatype 	ParticleType, orgdatatype[2];
	MPI_Aint 		offsets[2], extent;
	PetscMPIInt     blockcounts[2];
	PetscRandom		rctx;
	PetscScalar		cf_rand;
	PetscInt		DisplayOutput, NumberLostParticles, LostParticle, LostParticle_Recv;
	PetscInt		mod, nnel;
	DMDACoor3d 		coord_elemQ2[MAX_nnel];
	DAVPElementType vpt_element_type;
	PetscInt 		NeighborCPU[3][3][3], StorageSpaceNeighborCPU[3][3][3], cpu_x, cpu_y, cpu_z, storage_space_cpu;
    
	/* Some comments:
	 *  The array MaxNumCPU will NOT scale well in parallel; as it is an integer array it's probably not much of an issue
	 *
	 */


	DisplayOutput = 0;	// 0-little output, 1-Lot's of output, 2-mostly timing output

	if (DisplayOutput == 1){
		PetscPrintf(PETSC_COMM_WORLD,"Starting GetParticleNodes routine \n");
	}

	vpt_element_type = C->type;
	nnel            = C->nnel;

	/* Detect potential bugs */
	if (sizeof(Particles)/sizeof(PetscScalar) != particle_props){
		PetscPrintf(PETSC_COMM_WORLD,"************************************************** \n");
		PetscPrintf(PETSC_COMM_WORLD,"Oops could it be that the number of particle props are defined wrong in LaMEM.h? \n");
		PetscPrintf(PETSC_COMM_WORLD,"Is specified as %lld, but I computed %lld \n", (LLD)particle_props, (LLD)(sizeof(Particles)/sizeof(PetscScalar)) );
		PetscPrintf(PETSC_COMM_WORLD,"************************************************** \n");

		MPI_Barrier(PETSC_COMM_WORLD);

		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"LaMEM_Particles.c: Are particle properties wrongly defined?");

		//ierr = MPI_Abort(PETSC_COMM_WORLD,1); CHKERRQ(ierr);
	}

	PetscTime(&cputime_start);
    if (DisplayOutput > 1){
        PetscTime(&cputime_start0);
	}
    
	// # of nodes in all directions [takes periodicity into account]
	nx = user->finest_nnode_x;
	ny = user->finest_nnode_y;
	nz = user->finest_nnode_z;

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&size);
	ComputeNeighbors(da, NeighborCPU);			/* Compute CPU neighbors to the current processor*/

	/* We temporarily store particles that will be send over to another CPU in a 2D temporary array called ParticlesSendAway
	 * As in 3D with structured meshes we can only have a maximum of 26 neighbours, the temporary array should be of size [LengthParticleCache][27]
	 * To be consistent with the 3D neighbor array defined above we here create a 3D numbering array (could probably be done in a nicer way).
	 *
	 */
	num=0;
	for (i=0; i<3; i++){
		for (j=0; j<3; j++){
			for (k=0; k<3; k++){
				StorageSpaceNeighborCPU[i][j][k] 	= 	num;
				num 								=	num+1;
			}
		}
	}

	// Initialize the Random number generator
	ierr = PetscRandomCreate(PETSC_COMM_WORLD, &rctx); 	CHKERRQ(ierr); // DAM-freed
	ierr = PetscRandomSetFromOptions(rctx); 			CHKERRQ(ierr);

	// Allocate 2D array for sending particles around

	//LengthParticleCache = 100000;
	LengthParticleCache = (PetscInt)(0.1*(PetscScalar)user->num_particle_local); // no more then 10% of the particles on this proc should be send away.
	if (DisplayOutput == 1){
		PetscPrintf(PETSC_COMM_WORLD,"num_particle_local=%lld LengthParticleCache=%lld \n", (LLD)(user->num_particle_local),(LLD)LengthParticleCache );
	}

	ierr = PetscMalloc( sizeof(Particles*) *(size_t)LengthParticleCache, &ParticlesSendAway ); CHKERRQ(ierr);
	for( i=0; i< LengthParticleCache; i++ ) {
		ierr = PetscMalloc( sizeof(Particles) * 27, &ParticlesSendAway[i] ); CHKERRQ(ierr);		// should be not size but 27
	}

	//	ParticlesLost =  (Particles *)malloc(sizeof(Particles )*LenghParticleCache);
	ierr = PetscMalloc(sizeof(Particles )*(size_t)LengthParticleCache, &ParticlesLost); CHKERRQ(ierr);
	ierr = PetscMemzero(ParticlesLost, sizeof(Particles )*(size_t)LengthParticleCache);CHKERRQ(ierr);


	/*----------------------------------------------------------------------------------------
	 * Define the MPI structured datatype for sending particle information around
	 *
	 * ParticleType is the derived data type with currently two sections (double and double)
	 *  this can however easily be changed in future by following the rules below
	 */

	/* Allocates the 3 MPI_DOUBLE for the double fields x,y,z */
	offsets[0]     = 0;				orgdatatype[0] = MPI_DOUBLE;	  	blockcounts[0] = 3;

	/* Allocates space for (PARTICLE_PROPS-3) MPI_DOUBLE fields after figuring  out the offset by getting size of MPIU_INT */
	MPI_Type_extent(MPI_DOUBLE, &extent);
	offsets[1] 	 = 3 * extent;	 	orgdatatype[1] = MPI_DOUBLE;		blockcounts[1] = particle_props-3;

	/* Now define structured type and commit it */
	MPI_Type_struct(2, blockcounts, offsets,   orgdatatype, &ParticleType);
	MPI_Type_commit(&ParticleType);

	/*
	 * End of definition of ParticleType
	 * ----------------------------------------------------------------------------------------*/



	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(da,&xs_g,&ys_g,&zs_g,&xm_g,&ym_g,&zm_g); CHKERRQ(ierr);

	ierr = PetscMalloc((size_t)(10*xm_g*ym_g*zm_g)*sizeof(PetscInt), &NumberParticleElement); CHKERRQ(ierr);
	///	PetscMalloc((1*(user->nnode_x)*(user->nnode_y)*(user->nnode_z))*sizeof(PetscInt), &NumberParticleElement);
	for (ix=xs_g; ix<(xs_g+xm_g-1); ix++){
		for (iy=ys_g; iy<(ys_g+ym_g-1); iy++){
			for (iz=zs_g; iz<(zs_g+zm_g-1); iz++){
				NumberParticleElement[ (iz-zs_g)*( (xm_g-1)*(ym_g-1)) + (iy-ys_g)*(xm_g-1) + (ix-xs_g) ] = 0;
			}
		}
	}


	ierr = DMGetCoordinateDM(da,&cda); 		CHKERRQ(ierr);			//coordinates
	ierr = DMGetCoordinatesLocal(da,&gc); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,gc,&coords); 	CHKERRQ(ierr);


	// Initialize stuff
	for (i=0; i<MaxNumCPU; i++){ NumParticlesMove[i] = 0;};   // DOES NOT SCALE !!
	DeactivatedTracers 				  = 0;
	TracersOutsideComputationalDomain = 0;
	TracersSendToOtherCPU			  = 0;
	if (DisplayOutput == 1){
		PetscPrintf(PETSC_COMM_WORLD,"Starting particle loop \n");
	}

	/* Loop over all particles and determine which ones should be send over to another CPU */
	for (ipart=0; ipart<user->num_particle_local; ipart++){

		ParticleLocal 	= 	user->ParticlesLocal[ipart];	// extract particle

		if (ParticleLocal.phase>=0){	// only do this for active particles


			ix 				= 	((PetscInt) ParticleLocal.ix);
			iy 				= 	((PetscInt) ParticleLocal.iy);
			iz 				= 	((PetscInt) ParticleLocal.iz);


			/* If a particle is close to a boundary of a current proc. verify that it really still is on the correct cpu */

			// Figure out whether the particle is on the current CPU and if yes, what ix,iy & iz are
			TracerDeactivate =	 0;
			ElementFound 	 =	 0;
			if ( (ix<xs_g) || (ix>(xs_g+xm_g-2)) || (iy<ys_g) || (iy>(ys_g+ym_g-2)) || (iz<zs_g) || (iz>(zs_g+zm_g-2))   ){
				TracerDeactivate = 1;	ElementFound      = 1;  	// emergency break  [THIS SHOULD NEVER HAPPEN, SINCE PARTICLES
				// SHOULD HAVE BEEN DISTRIBUTED CORRECTLY IN A PREVIOUS STEP!

				PetscPrintf(PETSC_COMM_SELF,"[rank %d]   %d %d %d    in [%d %d] [%d %d] [%d %d]  \n",rank,ix,iy,iz , xs_g,xs_g+xm_g-2,ys_g,ys_g+ym_g-2,zs_g,zs_g+zm_g-2);

				PetscPrintf(PETSC_COMM_SELF,"Something went horribly wrong on rank %lld, with Particle %g, since it should not be on this PROC \n",(LLD)rank,ParticleLocal.num);

				SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"GetParticleNodes: Something went horribly wrong on rank %lld, with Particle %g, since it should not be on this PROC");

				//ierr = MPI_Abort(PETSC_COMM_WORLD,1); CHKERRQ(ierr);

			}
			else {
				// Find the 8 nodes that surround this particle. Note that this algorithm may fail if
				//  the grid is too distorted and particles have to move several cells
				ElementFound 		= 0;
				TracerDeactivate 	= 0;
			}

			cpu_x = 1;	cpu_y = 1;	cpu_z = 1;
			GetParticleElementAndLocalCoordinates(&ix,&iy,&iz,xs_g,ys_g,zs_g,xm_g,ym_g,zm_g,xs,ys,zs,xm,ym,zm,
					coords,&ParticleLocal,user, vpt_element_type,	nx,ny,
					&TracerDeactivate, 	&ElementFound, &cpu_x, &cpu_y, &cpu_z, &DisplayOutput);

			cpu 				=	NeighborCPU[cpu_y][cpu_x][cpu_z];					// to which PROC do we send these particles?
			storage_space_cpu 	=	StorageSpaceNeighborCPU[cpu_y][cpu_x][cpu_z];		// where in the ParticlesSendAway array do we store these data?
			ParticleLocal.cpu 	= 	((double) cpu);



			if ( (ix>=(nx-1)) || (iy>=(ny-1)) || (iz>=(nz-1))   ){
				// Particle send outside the bounding box
				ParticleLocal.phase  = ((double) -4);		// mark it for deletion on the current processor

				// deactivate it locally
				TracerDeactivate = 1;
				cpu=-2;	// ensure that it is not send to a different cpu
			}


			/* for debugging ######################### what a crap!
		//	if ((ParticleLocal.num== 2622415) && (1==1)){
			if ((ipart==32202) && (DisplayOutput==1)){
				PetscPrintf(PETSC_COMM_SELF,"Particle %10.0f with coords = [%g,%g,%g] and [eta,zetha,phi]=[%f,%f,%f] is on rank %lld, with [cpu_x,cpu_y,cpu_z]=[%lld,%lld,%lld], [ix,iy,iz]=[%lld,%lld,%lld],"
						"[xs,ys,zs]=[%lld,%lld,%lld]-[%lld,%lld,%lld]  and goes to rank %lld; x_left = %g TracerDeactivate=%lld\n",ParticleLocal.num, ParticleLocal.x, ParticleLocal.y, ParticleLocal.z,
						ParticleLocal.eta, ParticleLocal.zetha,ParticleLocal.phi,
						rank,
						cpu_x,cpu_y,cpu_z,ix,iy,iz,xs,ys,zs,xs+xm,ys+ym,zs+zm, cpu, user->x_left, TracerDeactivate);
			}
			// ######################################  */


			if ( (cpu != rank) && (cpu != -2)) {
				/* We'll have to send this particle to a new cpu.
				 * Put them in a queue first and send them around in a next step */

				if (NumParticlesMove[cpu]<LengthParticleCache){
					ParticlesSendAway[NumParticlesMove[cpu]][storage_space_cpu] = ParticleLocal;
					NumParticlesMove[cpu] = NumParticlesMove[cpu] + 1;

				}
				else {
					PetscPrintf(PETSC_COMM_WORLD,"Oops, you try to send too many particles to another CPU!");
					PetscPrintf(PETSC_COMM_WORLD," Increase LengthParticleCache in GetParticleNodes (currently %lld) and recompile code!",(LLD)LengthParticleCache);

					SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER,"GetParticleNodes: Oops, you try to send too many particles (more than %lld) to another CPU! Increase LengthParticleCache in LaMEM_Particle.c", (LLD)LengthParticleCache);


					//	ierr = MPI_Abort(PETSC_COMM_WORLD,1); CHKERRQ(ierr);
				}
				TracersSendToOtherCPU = TracersSendToOtherCPU+1;
				ParticleLocal.phase  = ((double) -4);		// mark it for deletion on the current processor

				// deactivate it locally
				TracerDeactivate = 1;
			}

			/* Particle is outside the bounds of the computational domain; delete it */
			if (cpu == -2){
				// Particle send outside the bounds
				ParticleLocal.phase  = ((double) -4);		// mark it for deletion on the current processor

				// deactivate it locally
				TracerDeactivate = 1;

			}






			if (TracerDeactivate==1){
				ParticleLocal.phase  = ((double) -4);		// mark it for deletion on current cpu

			}
			else {
				// Keep track of how many particles are available in every cell on the current CPU
				//  this way we can add particle in almost empty cells
				NumberParticleElement[ (iz-zs_g)*( (xm_g-1)*(ym_g-1)) + (iy-ys_g)*(xm_g-1) + (ix-xs_g) ] = NumberParticleElement[ (iz-zs_g)*( (xm_g-1)*(ym_g-1)) + (iy-ys_g)*(xm_g-1) + (ix-xs_g) ] + 1;
			}

		}

		/* Put particle back in array */
		user->ParticlesLocal[ipart] = ParticleLocal;


	} // loop over particles


	if (DisplayOutput == 1){
		PetscPrintf(PETSC_COMM_WORLD,"rank=%lld Finished Getting Nodes %lld \n", (LLD)rank,(LLD)ipart);

		// information about tracers outside computational box
		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank=%lld, particles deleted because outside computational domain: %lld \n",(LLD)rank, (LLD)TracersOutsideComputationalDomain);
		PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);

		// information about how many tracers should be send to another CPU
		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank=%lld, particles to be send to another CPU: %lld \n",(LLD)rank, (LLD)TracersSendToOtherCPU);
		PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
	}



	/* At this stage, we have processed all particles
	 *  The next steps are:
	 * 		(1)  Reshuffle the list on the local CPU
	 *  	(2)  Receive and send particles to other CPUs 	(& add to local list)
	 *  	(3)  Check for empty elements and inject particles if empty cells are found required.
	 */

	/* Reorder the local list of particles */
	num_particle_local_new = 0;
	num_backwards          = 0;
	for (ipart=0; ipart<user->num_particle_local; ipart++){
		ParticleLocal 		   = user->ParticlesLocal[ipart];	// extract particle
		if (ParticleLocal.phase<0){
			num_backwards = num_backwards+1;
		}
		else{
			user->ParticlesLocal[ipart-num_backwards] 	= ParticleLocal;
			num_particle_local_new 						= num_particle_local_new+1;
		}
	}
	user->num_particle_local = num_particle_local_new;

	if (DisplayOutput == 1){
		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank=%lld, particles left after reordering but before redistribution: %lld \n",(LLD)rank, (LLD)(user->num_particle_local));
		PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
	}
    if (DisplayOutput > 1){
        PetscTime(&cputime_end0);
        PetscPrintf(PETSC_COMM_WORLD,"#  reordering but not yet redistributing took %g s \n",cputime_end0 - cputime_start0);
	}

	/* Send particles around */
	MPI_Barrier(PETSC_COMM_WORLD);
    
    if (DisplayOutput > 1){
        PetscTime(&cputime_start0);
	}

	SendParticleToCorrectCPU(NumParticlesMove, StorageSpaceNeighborCPU, NeighborCPU, ParticlesSendAway, user, DisplayOutput);
    if (DisplayOutput > 1){
        PetscTime(&cputime_end0);
        PetscPrintf(PETSC_COMM_WORLD,"#  Sending particles to other PROCS took %g s \n",cputime_end0 - cputime_start0);
	}


	if (DisplayOutput == 1){
		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank=%lld, Starting computation of local coordinates  \n",(LLD)rank);
		PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
	}

	/* Make a loop over all particles and compute its local coordinates */
    PetscTime(&cputime_start0);
    
	LostParticle 			=	0;
	for (ipart=0; ipart<user->num_particle_local; ipart++){
		ParticleLocal = user->ParticlesLocal[ipart];	// extract particle

		ix 					= ((PetscInt) ParticleLocal.ix);
		iy 					= ((PetscInt) ParticleLocal.iy);
		iz 					= ((PetscInt) ParticleLocal.iz);
		TracerDeactivate 	= 0;
		ElementFound 	 	= 0;



		cpu_x = 1;	cpu_y = 1;	cpu_z = 1;
		ierr = GetParticleElementAndLocalCoordinates(&ix,&iy,&iz,xs_g,ys_g,zs_g,xm_g,ym_g,zm_g,xs,ys,zs,xm,ym,zm,
				coords,&ParticleLocal, user, vpt_element_type,	nx,ny,
				&TracerDeactivate, 	&ElementFound, &cpu_x, &cpu_y, &cpu_z, &DisplayOutput);	CHKERRQ(ierr);


		cpu 				=	NeighborCPU[cpu_y][cpu_x][cpu_z];
		if (cpu != rank){
			/* Oops, it turns out that this particle is actually NOT on the correct CPU (even though it should have been with the
			 * 	redistribution algorithm at the beginning).
			 * Therefore, add the particle to an array with 'lost' particles
			 */
			PetscPrintf(PETSC_COMM_SELF," Oops; particle %lld on rank %lld was supposed to be on cpu =%lld \n",(LLD)ipart, (LLD)rank, (LLD)cpu);


			// Add it to array with 'lost' particles
			ParticlesLost[LostParticle] = 	ParticleLocal;
			if (LostParticle>LengthParticleCache){

				SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER, "You lost more particles than there is space in LenghParticleCache, which is %lld",(LLD)LengthParticleCache);

				SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"You lost more particles than there is space in LenghParticleCache");


				//ierr = MPI_Abort(PETSC_COMM_WORLD,1); CHKERRQ(ierr);
			}
			LostParticle 		= 	LostParticle + 1;

			ParticleLocal.phase = 	-4.0;

		}
		ParticleLocal.cpu 	= 	((double) cpu);

		/* Put particle back in array */
		user->ParticlesLocal[ipart] = ParticleLocal;

	}
	if (DisplayOutput == 1){
		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank=%lld, finished computation of local coordinates  \n",(LLD)rank);
		PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);

	}
	/* ----------------------------------------------------------------------------------------------------
	 * Deal with LOST particles - this should ofcourse never really happen in reality.
	 */

	/* (1) Send lost particles to every CPU 
      THIS IS NOT SCALABLE - TO BE FIXED. Doesn't seem to be a major bottleneck though.
     */
	NumberLostParticles 	=	0;
	for (iproc=0; iproc<size; iproc++){

		if (iproc==rank){	LostParticle_Recv = LostParticle;	}
		else{				LostParticle_Recv =	0;				}
        
        MPI_Bcast(&LostParticle_Recv, 1, MPIU_INT, iproc, PETSC_COMM_WORLD);

		for (i=0; i<LostParticle_Recv; i++){

			// send lost particle
			if (iproc==rank){	ParticleLost = ParticlesLost[i];	}
			MPI_Bcast(&ParticleLost, 1, ParticleType, iproc, PETSC_COMM_WORLD);

			// Add to list with local particles
			user->ParticlesLocal[user->num_particle_local + NumberLostParticles] = ParticleLost;

			NumberLostParticles 	=	NumberLostParticles + 1;
		}

	}
	MPI_Barrier(PETSC_COMM_WORLD);

	/* (2) Check whether the newly arrived particles are on the current CPU */
	for (ipart=user->num_particle_local; ipart<user->num_particle_local+NumberLostParticles; ipart++){
		ParticleLocal = user->ParticlesLocal[ipart];	// extract particle
		ix 					= ((PetscInt) ParticleLocal.ix);
		iy 					= ((PetscInt) ParticleLocal.iy);
		iz 					= ((PetscInt) ParticleLocal.iz);
		TracerDeactivate 	= 0;
		ElementFound 	 	= 0;

		if ( (ix>=xs_g) && (ix<(xs_g+xm_g)) && (iy>=ys_g) && (iy<(ys_g+ym_g)) && (iz>=zs_g) && (iz<(zs_g+zm_g))   ){
			cpu_x = 1;	cpu_y = 1;	cpu_z = 1;
			ierr = GetParticleElementAndLocalCoordinates(&ix,&iy,&iz,xs_g,ys_g,zs_g,xm_g,ym_g,zm_g,xs,ys,zs,xm,ym,zm,
					coords,&ParticleLocal, user, vpt_element_type,	nx,ny,
					&TracerDeactivate, 	&ElementFound, &cpu_x, &cpu_y, &cpu_z, &DisplayOutput); CHKERRQ(ierr);

			cpu 				=	NeighborCPU[cpu_y][cpu_x][cpu_z];
			if (cpu != rank){
				// not on current CPU; delete the particle
				ParticleLocal.phase = -4.0;
			}

		}
		else{
			ParticleLocal.phase = -4.0;
		}

		/* Put particle back in array */
		user->ParticlesLocal[ipart] = ParticleLocal;



	}
	user->num_particle_local = user->num_particle_local + NumberLostParticles;
    
    if (DisplayOutput > 1){
        PetscTime(&cputime_end0);
        PetscPrintf(PETSC_COMM_WORLD,"#  Dealing with lost particles took took %g s \n",cputime_end0 - cputime_start0);
	}
    
	/*
	 * All lost particles have been redistributed at this point
	 * -----------------------------------------------------------------------------------------------------*/


	/* Reshuffle the list of particles and get rid of potentially deleted particles */
	num_particle_local_new = 0;
	num_backwards          = 0;
	for (ipart=0; ipart<user->num_particle_local; ipart++){
		ParticleLocal 		   = user->ParticlesLocal[ipart];	// extract particle
		if (ParticleLocal.phase<0){
			num_backwards = num_backwards+1;
		}
		else{
			if (num_backwards>0){
				user->ParticlesLocal[ipart-num_backwards] 	= ParticleLocal;
			}
			num_particle_local_new 						= num_particle_local_new+1;

		}
	}
	user->num_particle_local = num_particle_local_new;


	if (1==1){
		// INJECT Particles
		/* The only thing left to do is to check whether there are empty cells in the local domain
		 *  if yes, we add some particles to this cell
		 *
		 */

		AddRandomNoiseParticles = 1;
		ierr = PetscOptionsGetInt(PETSC_NULL,"-AddRandomNoiseParticles",&AddRandomNoiseParticles	, PETSC_NULL); CHKERRQ(ierr);


		for (ix=xs; ix<(xs+xm-1); ix++){
			for (iy=ys; iy<(ys+ym-1); iy++){
				for (iz=zs; iz<(zs+zm-1); iz++){
					num_local_part = 	NumberParticleElement[ (iz-zs_g)*( (xm_g-1)*(ym_g-1)) + (iy-ys_g)*(xm_g-1) + (ix-xs_g) ];
					if (num_local_part<=user->NumParticlesToStartInjection){
						// inject particles

						// get shape function
						GetElementCoords(coord_elem, coords, ix,iy,iz,0);
						CorrectElementCoordsForPeriodicity( coord_elem, ix, iy, user, 0);

						if (1==1){
							// add new particles
							dx_natural = 2.0/((PetscScalar) (user->NumPartX));
							dy_natural = 2.0/((PetscScalar) (user->NumPartY));
							dz_natural = 2.0/((PetscScalar) (user->NumPartZ));
							for (iix=0; iix<user->NumPartX; iix++){
								for (iiy=0; iiy<user->NumPartY; iiy++){
									for (iiz=0; iiz<user->NumPartZ; iiz++){

										/* Generate random eta,zetha,phi coordinates */

										/* Set local coordinates of the particle */
										ParticleLocal.eta     = ((PetscScalar)(iix))*dx_natural - 1.0  +  dx_natural/2;
										ParticleLocal.zetha   = ((PetscScalar)(iiy))*dy_natural - 1.0  +  dy_natural/2;
										ParticleLocal.phi     = ((PetscScalar)(iiz))*dz_natural - 1.0  +  dz_natural/2;

										/* Add random noise if required 		*/
										if (AddRandomNoiseParticles==1	){

											ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
											cf_rand = (cf_rand-0.5);
											ParticleLocal.eta	=	ParticleLocal.eta 	+ cf_rand*dx_natural;

											ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
											cf_rand = (cf_rand-0.5);
											ParticleLocal.zetha	=	ParticleLocal.zetha + cf_rand*dy_natural;

											ierr = PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
											cf_rand = (cf_rand-0.5);
											ParticleLocal.phi	=	ParticleLocal.phi 	+ cf_rand*dz_natural;

										}

										ParticleLocal.ix  = ((double) ix  );
										ParticleLocal.iy  = ((double) iy  );
										ParticleLocal.iz  = ((double) iz  );
										ParticleLocal.cpu = ((double) rank);

										/* Compute real coordinates from these natural coordinates */
										x = ParticleLocal.eta; y = ParticleLocal.zetha; z = ParticleLocal.phi;
										ShapeVel[0]		=	0.125*(1.0-x)*(1.0-y)*(1.0-z);  ShapeVel[1]=0.125*(1.0+x)*(1.0-y)*(1.0-z);
										ShapeVel[2]		=	0.125*(1.0+x)*(1.0+y)*(1.0-z);  ShapeVel[3]=0.125*(1.0-x)*(1.0+y)*(1.0-z);
										ShapeVel[4]		=	0.125*(1.0-x)*(1.0-y)*(1.0+z);  ShapeVel[5]=0.125*(1.0+x)*(1.0-y)*(1.0+z);
										ShapeVel[6]		=	0.125*(1.0+x)*(1.0+y)*(1.0+z);  ShapeVel[7]=0.125*(1.0-x)*(1.0+y)*(1.0+z);
										ParticleLocal.x =   0;
										ParticleLocal.y =   0;
										ParticleLocal.z =   0;
										for (i=0; i<8; i++){
											ParticleLocal.x = ParticleLocal.x + ShapeVel[i]*coord_elem[i].x;
											ParticleLocal.y = ParticleLocal.y + ShapeVel[i]*coord_elem[i].y;
											ParticleLocal.z = ParticleLocal.z + ShapeVel[i]*coord_elem[i].z;
										}

										/* Give the particle a phase etc. */

										/* If we are using a FDSTAG formulation, we can compute the to-be-injected phase from the surrounding nodes.
										 * 	if not, we use a predefined injection phase
										 */


										InjectionPhase = user->ParticleInjectionPhase;
										if (vpt_element_type==DAVP_FDSTAG){
											PetscInt DominantPhase_WithAir, DominantPhase_WithoutAir, DominantPhase_WithoutAir_TopCell;


											ierr =  FDSTAG_DominantPhaseAtCorners(user, ix, iy,iz, &DominantPhase_WithAir, &DominantPhase_WithoutAir, &DominantPhase_WithoutAir_TopCell); CHKERRQ(ierr);
											if (DominantPhase_WithAir>-1){
												InjectionPhase = DominantPhase_WithAir;
											}
										}

										ParticleLocal.phase = InjectionPhase;
										ParticleLocal.num   = -1;

										if (user->num_particle_local< (user->MaxNumLocalParticles-1) ){ 		// if there is still space available
											user->ParticlesLocal[user->num_particle_local] 	= ParticleLocal;
											user->num_particle_local = user->num_particle_local+1;

										}

									}
								}
							}

						}
					}
				}
			}
		}
		// End of particle injection
		if (DisplayOutput == 1){
			PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank=%lld, particles after injection: %lld max # of particle that have space on current proc: %lld \n",(LLD)rank, (LLD)(user->num_particle_local), (LLD)(user->MaxNumLocalParticles));
			PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
		}
	}

	/* If we have Q2P1 elements, there is still some more work to do,
	 *  since the routines above computed natural coordinates for the case of Q1P0 elements
	 */
	if( __ELEMENT_TYPE__ == ELEMENT_Q2P1 ) {
		for (ipart=0; ipart<user->num_particle_local; ipart++){
			ParticleLocal 		   = user->ParticlesLocal[ipart];	// extract particle

			// re-construct element coordinates
			ix = ((PetscInt) ParticleLocal.ix); LaMEMMod(ix, 2, &mod); if (mod>0){ix = ix-1;};
			iy = ((PetscInt) ParticleLocal.iy); LaMEMMod(iy, 2, &mod); if (mod>0){iy = iy-1;};
			iz = ((PetscInt) ParticleLocal.iz); LaMEMMod(iz, 2, &mod); if (mod>0){iz = iz-1;};

			// local coordinates
			GetElementCoords(coord_elemQ2, coords, ix,iy,iz, 1);
			CorrectElementCoordsForPeriodicity( coord_elemQ2, ix, iy, user, 1);

			// natural coordinates
			eta 	= ParticleLocal.eta;	zetha 	= ParticleLocal.zetha;		phi 	= ParticleLocal.phi;
			NaturalCoordinates( nnel, ParticleLocal.x, ParticleLocal.y, ParticleLocal.z, &eta, &zetha, &phi, coord_elemQ2, &inside);
			ParticleLocal.eta 	= eta;	 	ParticleLocal.zetha	= zetha;		ParticleLocal.phi	= phi;

			user->ParticlesLocal[ipart]  = ParticleLocal;	// extract particle

		}
	}


	PetscTime(&cputime_end);
	if (DisplayOutput == 1){
		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"#rank=%lld, Attributing %lld particles to elements took %g s \n",(LLD)rank, (LLD)(user->num_particle_local),	cputime_end-cputime_start);
		PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
	} else if (DisplayOutput == 2){
		PetscPrintf(PETSC_COMM_WORLD,"#rank=%lld, Attributing %lld particles to elements took %g s \n",(LLD)rank, (LLD)(user->num_particle_local),	cputime_end-cputime_start);
    } else if (DisplayOutput == 0){
		PetscPrintf(PETSC_COMM_WORLD,"#rank=%lld, Attributing %lld particles to elements took %g s \n",(LLD)rank, (LLD)(user->num_particle_local),	cputime_end-cputime_start);
	}

	/* Sum total amount of particles */
	ierr = MPI_Allreduce(&user->num_particle_local,	&TotalNumberOfParticles,	1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"The total # of particles present is %lld \n",(LLD)TotalNumberOfParticles);




	/* Give warnings if required */
	if (DeactivatedTracers>0){

		PetscPrintf(PETSC_COMM_SELF,"#********************************************************** \n");
		PetscPrintf(PETSC_COMM_SELF,"# Warning: %lld particles could not be located on proc %lld!   \n",(LLD)DeactivatedTracers, (LLD)rank);
		PetscPrintf(PETSC_COMM_SELF,"# Try remeshing more frequently.   %lld %lld \n",(LLD)DeactivatedTracers, (LLD)rank);
		PetscPrintf(PETSC_COMM_SELF,"# ********************************************************** \n");
	}

	// Cleaning up
	ierr = DMDAVecRestoreArray(cda,gc,&coords); CHKERRQ(ierr);
	//	ierr = DMDestroy( cda ); CHKERRQ(ierr);
	//	ierr = VecDestroy( gc ); CHKERRQ(ierr);


	//	for (i=0; i<LenghParticleCache; i++){
	//		free(ParticlesSendAway[i]);
	//	}
	//	free(ParticlesSendAway);

	for (i=0; i<LengthParticleCache; i++){
		ierr = PetscFree(ParticlesSendAway[i]); CHKERRQ(ierr);
	}
	ierr = PetscFree(ParticlesSendAway); CHKERRQ(ierr);

	ierr =  PetscFree(ParticlesLost); 		   CHKERRQ(ierr);
	ierr = 	PetscFree(NumberParticleElement);  CHKERRQ(ierr);
	ierr = PetscRandomDestroy(&rctx); CHKERRQ(ierr);		// Destroy random context

	MPI_Type_free( &ParticleType );

	PetscFunctionReturn(0);

}
/* ==================================================================================================== */

/* ====================================================================================================
 * This finds the element in which the Particle resides (or determines that a particle is outside the current
 * CPU, in which case it'll tell the most likely direction the particle traveled.
 * */
#undef __FUNCT__
#define __FUNCT__ "SendParticleToCorrectCPU"
PetscErrorCode SendParticleToCorrectCPU(PetscInt NumParticlesMove[MaxNumCPU], PetscInt StorageSpaceNeighborCPU[3][3][3], PetscInt NeighborCPU[3][3][3], Particles **ParticlesSendAway,
		UserContext *user, PetscInt DisplayOutput)
{
	PetscMPIInt     rank, size, num_data_send, i, iproc, proc, ipart;
	PetscErrorCode  ierr;
	PetscInt 		ii, jj, kk, storage_cpu, NeighborCPU_vec[27], num, iii, jjj;
	MPI_Status		status;
	MPI_Datatype 	ParticleType, orgdatatype[2];
	MPI_Aint 		offsets[2], extent;
	PetscMPIInt     blockcounts[2];
	PetscInt 		tag=1;
	Particles 		*Particles_Send, *Particles_Recv;


	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&size);

	/*----------------------------------------------------------------------------------------
	 * Define the MPI structured datatype for sending particle information around
	 *
	 * ParticleType is the derived data type with currently two sections (double and double)
	 *  this can however easily be changed in future by following the rules below
	 */

	/* Allocates the 3 MPI_DOUBLE for the double fields x,y,z */
	offsets[0]     = 0;				orgdatatype[0] = MPI_DOUBLE;	  	blockcounts[0] = 3;

	/* Allocates space for (PARTICLE_PROPS-3) MPI_DOUBLE fields after figuring  out the offset by getting size of MPIU_INT */
	MPI_Type_extent(MPI_DOUBLE, &extent);
	offsets[1] 	 = 3 * extent;	 	orgdatatype[1] = MPI_DOUBLE;		blockcounts[1] = particle_props-3;

	/* Now define structured type and commit it */
	MPI_Type_struct(2, blockcounts, offsets,   orgdatatype, &ParticleType);
	MPI_Type_commit(&ParticleType);

	/*
	 * End of definition of ParticleType
	 * ----------------------------------------------------------------------------------------*/

    
    /* create a vector that contains the neighbouring PROCs numbers */
    num=0;
    for (ii=0; ii<3; ii++){
        for (jj=0; jj<3; jj++){
            for (kk=0; kk<3; kk++){
                NeighborCPU_vec[num] = NeighborCPU[ii][jj][kk];
                num                  = num+1;
            }
        }
    }
    
    

    for (iii=0; iii<27; iii++){  //loop over 27 neighbours including itself
    	iproc = NeighborCPU_vec[iii];       // negative numbers indicate non-existing neighbours (e.g. if the cell is at a boundary)

    	if ((iproc!=rank) && (iproc>=0)){

    		// Send data from rank to iproc: in this case info about how many particles will follow
    		num_data_send = NumParticlesMove[iproc];
    		ierr = MPI_Send(&num_data_send, 1, MPIU_INT, iproc,  9, PETSC_COMM_WORLD); CHKERRQ(ierr);  // this command is being send procs^2 times

    		// PetscPrintf(PETSC_COMM_SELF," I am proc %lld sending to proc %lld: %lld particles \n",(LLD)rank,(LLD)iproc,(LLD)num_data_send);

    		if (num_data_send>0){
    			// send the actual data

    			/* Find out where this is stored in the array ParticlesSendAway */
    			storage_cpu=-1;
    			for (ii=0; ii<3; ii++){
    				for (jj=0; jj<3; jj++){
    					for (kk=0; kk<3; kk++){
    						if (NeighborCPU[ii][jj][kk]==iproc){
    							storage_cpu =  StorageSpaceNeighborCPU[ii][jj][kk];		// that's where it was stored
    						}

    					}
    				}
    			}
    			if (storage_cpu==-1){
    				/* should not occur but catch error if it happens */
    				SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"SendParticleToCorrectCPU: dont find the correct neighbour storage space");
    			}


    			// Prepare the data for sending
    			PetscMalloc((size_t)num_data_send*sizeof(Particles),  	&Particles_Send	);
    			for (i=0; i<num_data_send; i++){
    				Particles_Send[i] = 	ParticlesSendAway[i][storage_cpu];
    			}

    			// Send the data
    			ierr = MPI_Send(Particles_Send, num_data_send, ParticleType, iproc,  tag, PETSC_COMM_WORLD); CHKERRQ(ierr);

    			PetscFree(Particles_Send);	// free array
    		}
    	}

    	else if (iproc==rank) {

    		// Receive data from neighbours
    		for (jjj=0; jjj<27; jjj++){  //loop over 27 neighbours including itself
    			proc = NeighborCPU_vec[jjj];

    			if ((proc!=rank) && (proc>=0)){
    				// Number of particles that are coming
    				ierr = MPI_Recv(&num_data_send, 1, MPIU_INT, proc,  MPI_ANY_TAG, PETSC_COMM_WORLD, &status); CHKERRQ(ierr);

    				// PetscPrintf(PETSC_COMM_SELF," I am proc %lld receiving from proc %lld: %lld particles \n",(LLD)proc,(LLD)rank,(LLD)num_data_send);

    				if (num_data_send>0){
    					// Allocate appropriate array
    					PetscMalloc((size_t)num_data_send*sizeof(Particles),  	&Particles_Recv	);


    					// Receive data
    					ierr = MPI_Recv(Particles_Recv, num_data_send, ParticleType, proc,  tag, PETSC_COMM_WORLD, &status); CHKERRQ(ierr);

    					// Add data to local list
    					for (ipart=0; ipart<num_data_send; ipart++){
    						if (user->num_particle_local< (user->MaxNumLocalParticles-1) ){ 		// if there is still space available
    							user->ParticlesLocal[user->num_particle_local] 	= Particles_Recv[ipart];

    							user->num_particle_local = user->num_particle_local+1;

    						}
    						else {
    							PetscPrintf(PETSC_COMM_SELF,"# ********************************************************** \n");
    							PetscPrintf(PETSC_COMM_SELF,"# Warning: insufficient space for particles on proc %lld!   \n",(LLD)rank);
    							PetscPrintf(PETSC_COMM_SELF,"# I currently have %lld particles and I can hold a maximum of %lld particles on this CPU !   \n",(LLD)(user->num_particle_local), (LLD)(user->MaxNumLocalParticles));
    							PetscPrintf(PETSC_COMM_SELF,"# Increase particle_buffer!   \n");
    							PetscPrintf(PETSC_COMM_SELF,"# ********************************************************** \n");
    						}
    					}
    					if (DisplayOutput == 1){
    						PetscPrintf(PETSC_COMM_SELF," Receiving %lld particles on proc %lld from proc %lld \n",(LLD)num_data_send,(LLD)rank,(LLD)proc);

    					}
    					PetscFree(Particles_Recv);	// free array
    				}

    			}

    		}
    	}
    }


	//PetscSynchronizedFlush(PETSC_COMM_WORLD);
	if (DisplayOutput == 1){
		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank=%lld, particles left after redistribution: %lld \n",(LLD)rank, (LLD)(user->num_particle_local));
		PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
	}

	MPI_Type_free(&ParticleType);

	PetscFunctionReturn(0);
}
/*
 * End of sending particles to correct CPU
 * ==================================================================================================== */


/* ====================================================================================================
 * This finds the element in which the Particle resides (or determines that a particle is outside the current
 * CPU, in which case it'll tell the most likely direction the particle traveled.
 * */
#undef __FUNCT__
#define __FUNCT__ "GetParticleElementAndLocalCoordinates"
PetscErrorCode GetParticleElementAndLocalCoordinates(
		PetscInt *ix_in, 	PetscInt *iy_in, 	PetscInt *iz_in,		PetscInt xs_g, 		PetscInt ys_g, 	PetscInt zs_g,
		PetscInt  xm_g, 	PetscInt ym_g, 		PetscInt zm_g, 			PetscInt xs, 		PetscInt ys, 	PetscInt zs,
		PetscInt  xm,		PetscInt ym, 		PetscInt zm,			DMDACoor3d ***coords, Particles 		*ParticleLocal_in,
		UserContext *user, 				DAVPElementType vpt_element_type,	PetscInt nx, 		PetscInt ny,
		PetscInt 	*TracerDeactivate_in, 		PetscInt *ElementFound_Outside,
		PetscInt *cpu_x_in, PetscInt *cpu_y_in, PetscInt *cpu_z_in, PetscInt *DisplayOutput)
{
	PetscMPIInt     rank;
	PetscErrorCode 	ierr;
	PetscInt 		n_iter;
	DMDACoor3d		coord_elem[8];
	PetscScalar 	eta, zetha, phi;
	PetscScalar		minx,maxx,miny,maxy,minz,maxz;
	PetscScalar		max_leftx, min_rightx, max_fronty, min_backy, max_botz, min_topz;
	PetscInt 		inside;
	PetscInt		ix_new=0, iy_new=0, iz_new=0, shift, ElementFound;
	PetscInt		ix,iy,iz,TracerDeactivate;
	PetscInt	 	cpu_x, cpu_y, cpu_z, PeriodicCorrection;
	Particles 		ParticleLocal;

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);



	// verify that starting location of particle is on current cpu
	n_iter 				= 	1;
	ElementFound 		= 	*ElementFound_Outside;
	ix 					=	*ix_in;
	iy 					=	*iy_in;
	iz 					=	*iz_in;
	TracerDeactivate 	= 	*TracerDeactivate_in;
	ParticleLocal 		=	*ParticleLocal_in;

	eta					=	ParticleLocal.eta;
	zetha				=	ParticleLocal.zetha;
	phi					=	ParticleLocal.phi;
	inside 				=	0;


	cpu_x 				=	*cpu_x_in;
	cpu_y 				=	*cpu_y_in;
	cpu_z 				=	*cpu_z_in;

	/* Check for periodicity */
	PeriodicCorrection = 0;
	if ((user->BC.LeftBound==3) || (user->BC.FrontBound==3)){
		if ((ParticleLocal.x> (user->x_left+user->W))  	& (user->BC.LeftBound==3)	)  {
			ParticleLocal.x 	=	ParticleLocal.x 	- 	user->W;
			ix 					=	ix 	-	(nx-2);
			cpu_x 				=	cpu_x+1;
			TracerDeactivate 	= 	0;
			ElementFound 		=	1;
			PeriodicCorrection 	=	1;
		}
		else if ((ParticleLocal.x< (user->x_left)) 		& (user->BC.LeftBound==3) 	)  {
			ParticleLocal.x 	=	ParticleLocal.x 	+ 	user->W;
			ix 					=	ix 	+	(nx-2);
			cpu_x 				=	cpu_x-1;
			TracerDeactivate 	= 	0;
			ElementFound 		=	1;
			PeriodicCorrection 	=	1;
		}
		if ((ParticleLocal.y> (user->y_front+user->L)) 	&  (user->BC.FrontBound==3) )  {
			ParticleLocal.y 	=	ParticleLocal.y 	- 	user->L;
			iy 					=	iy 	-	(ny-2);
			cpu_y 				=	cpu_y+1;
			TracerDeactivate 	= 	0;
			ElementFound 		=	1;
			PeriodicCorrection 	=	1;
		}
		else if ((ParticleLocal.y< (user->y_front)) 	&  (user->BC.FrontBound==3) )  {
			ParticleLocal.y 	=	ParticleLocal.y 	+ 	user->L;
			iy 					=	iy 	+	(ny-2);
			cpu_y 				=	cpu_y-1;
			TracerDeactivate 	= 	0;
			ElementFound 		=	1;
			PeriodicCorrection 	=	1;
		}
	}


	/* Determine whether ix,iy,iz are within the local domain */
	if ((iz==zs_g+zm_g-1)){
		// this is bad news:
		PetscPrintf(PETSC_COMM_SELF,"rank = %lld Particle:  [ix,iy,iz]=[%lld,%lld,%lld] domain=[%lld,%lld,%lld]-[%lld,%lld,%lld] has wrong iz coordinate and is actially outside the domain \n",(LLD)rank, (LLD)ix,(LLD)iy,(LLD)iz,(LLD)xs,(LLD)ys,(LLD)zs,(LLD)(xs+xm-1),(LLD)(zs+zm-1),(LLD)(zs+zm-1));

		//SETERRQ2(PETSC_ERR_MEMC,"GetParticleElementAndLocalCoordinates: Particle %f on rank %lld has the wrong iz coordinate and is actually outside the domain. It should have been deactivated", ParticleLocal.num, rank);

		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"GetParticleElementAndLocalCoordinates: Particle  has the wrong iz coordinate and is actually outside the domain");

		//ierr = MPI_Abort(PETSC_COMM_WORLD,1); CHKERRQ(ierr);

	}




	/* Determine the correct CPU of the particle */
	while ( (ElementFound<1) &&  (n_iter<18)){

		// Extract the suggested coordinates of the surrounding 8 nodes
		ierr = GetElementCoords(coord_elem, coords, ix,iy,iz, 0); 							CHKERRQ(ierr);

		ierr = CorrectElementCoordsForPeriodicity( coord_elem, ix, iy, user, 0);			CHKERRQ(ierr);

		ierr = ElementInfo(coord_elem,&minx,&maxx,&miny,&maxy,&minz,&maxz,
				&max_leftx, &min_rightx, &max_fronty, &min_backy, &max_botz, &min_topz);	CHKERRQ(ierr);

		// Check whether the current particle is in the element that is suggested
		if 	((ParticleLocal.x>=max_leftx) && (ParticleLocal.x<=min_rightx) &&
				(ParticleLocal.y>=max_fronty)&& (ParticleLocal.y<=min_backy)  &&
				(ParticleLocal.z>=max_botz)  && (ParticleLocal.z<=min_topz)    ){
			// the particle is clearly inside the element
			ElementFound 	= 	1;
			ix_new			= 	ix;
			iy_new			=	iy;
			iz_new			=	iz;

			// compute the particle's local coordinates in case of Q1P0 element
			eta					=	((ParticleLocal.x-minx)/(maxx-minx)-0.5)*2.0;	// initial guess
			zetha				=	((ParticleLocal.y-miny)/(maxy-miny)-0.5)*2.0; // initial guess
			phi					=	((ParticleLocal.z-minz)/(maxz-minz)-0.5)*2.0; // initial guess
			if( (vpt_element_type == DAVP_Q1P0)  || (vpt_element_type == DAVP_Q1Q1) || (vpt_element_type == DAVP_FDSTAG)) {
				PetscInt DisplayIterations;

				DisplayIterations=0;
				ierr = Inside_8NodeElement(ParticleLocal.x, ParticleLocal.y, ParticleLocal.z, &eta, &zetha, &phi, coord_elem,	&inside, DisplayIterations); CHKERRQ(ierr);

				// Ensure that the starting local coords make sense; if
				if ((eta<-1.0) || (eta>1.0) || (zetha<-1.0) || (zetha>1.0) || (phi<-1.0) || (phi>1.0)){
					// The particle 'should' be inside but the  Inside_8NodeElement algorithm sometimes fails.
					// in that case, deactivate the tracer altogether (doesn't happen much anyways)
					TracerDeactivate = 1;
				}

			}
			ParticleLocal.eta	=	eta;
			ParticleLocal.zetha	=	zetha;
			ParticleLocal.phi	=	phi;



		}


		if (	(ElementFound < 1)       &&
				(ParticleLocal.x>=minx) && 	(ParticleLocal.x<=maxx) &&
				(ParticleLocal.y>=miny) && 	(ParticleLocal.y<=maxy) &&
				(ParticleLocal.z>=minz) && 	(ParticleLocal.z<=maxz)     ){
			// could be in element, but requires more close checking

			if (n_iter>2){
				/* We're lost and trying to find  the location of the particle */
				eta 	= ((ParticleLocal.x-minx)/(maxx-minx)-0.5)*2.0;
				zetha 	= ((ParticleLocal.y-miny)/(maxy-miny)-0.5)*2.0;
				phi  	= ((ParticleLocal.z-minz)/(maxz-minz)-0.5)*2.0;
			}

			// check here with Newton-Rhapson iterations whether we're inside
			ierr = Inside_8NodeElement(ParticleLocal.x, ParticleLocal.y, ParticleLocal.z, &eta, &zetha, &phi, coord_elem,	&inside, 0); CHKERRQ(ierr);

			if 	(inside==1){
				// found!
				ElementFound 	= 	1;
				ix_new			= 	ix;
				iy_new			=	iy;
				iz_new			=	iz;
				if( (vpt_element_type == DAVP_Q1P0) || (vpt_element_type == DAVP_Q1Q1) || (vpt_element_type == DAVP_FDSTAG) ) {
					/* Note that FDSTAG is treated as if it would be a linear element */
					ParticleLocal.eta	=	eta;
					ParticleLocal.zetha	=	zetha;
					ParticleLocal.phi	=	phi;
				}

			}
			else if (inside==0){
				ElementFound 	= 	-1;

				// nope, let's check in which direction we should continue searching
				ix_new =ix; iy_new = iy; iz_new = iz;

				/* Our grid is typically mainly deformed in z-direction;
				 *  That's why this direction is given preference here */
				shift = 0;
				if (phi   <-1.0 && shift==0)	{ iz_new = iz-1;  shift=1;}
				if (phi   > 1.0 && shift==0)	{ iz_new = iz+1;  shift=1;}
				if (eta   <-1.0 && shift==0)	{ ix_new = ix-1;  shift=1;}
				if (eta   > 1.0 && shift==0)	{ ix_new = ix+1;  shift=1;}
				if (zetha <-1.0 && shift==0)	{ iy_new = iy-1;  shift=1;}
				if (zetha > 1.0 && shift==0)	{ iy_new = iy+1;  shift=1;}

			}

		}
		else if (	(ElementFound < 10) && (
				(ParticleLocal.x<minx) || 	(ParticleLocal.x>maxx) ||
				(ParticleLocal.y<miny) || 	(ParticleLocal.y>maxy) ||
				(ParticleLocal.z<minz) || 	(ParticleLocal.z>maxz) )    ){

			// must be outside element - find out towards which direction it most likely went
			ix_new =ix; iy_new = iy; iz_new = iz;
			if (ParticleLocal.x<minx){ ix_new = ix-1; }
			if (ParticleLocal.x>maxx){ ix_new = ix+1; }
			if (ParticleLocal.y<miny){ iy_new = iy-1; }
			if (ParticleLocal.y>maxy){ iy_new = iy+1; }
			if (ParticleLocal.z<minz){ iz_new = iz-1; }
			if (ParticleLocal.z>maxz){ iz_new = iz+1; }
			ElementFound      = -2;
		}


		// for debugging
		if ( (ParticleLocal.num==32202 )   && (*DisplayOutput==1) ){
			PetscPrintf(PETSC_COMM_SELF, "n_iter=%lld, rank=%lld, ElementFound=%lld coords=[%g,%g,%g] \n",(LLD)n_iter,(LLD)rank,(LLD)ElementFound,ParticleLocal.x,ParticleLocal.y,ParticleLocal.z );
			PetscPrintf(PETSC_COMM_SELF, "           ix,iy,iz            =[%lld,%lld,%lld] \n",(LLD)ix,(LLD)iy,(LLD)iz);
			PetscPrintf(PETSC_COMM_SELF, "           ix_new,iy_new,iz_new=[%lld,%lld,%lld] \n",(LLD)ix_new,(LLD)iy_new,(LLD)iz_new);
			PetscPrintf(PETSC_COMM_SELF, "           ind=[%10.0f] \n",ParticleLocal.num);
			PetscPrintf(PETSC_COMM_SELF, "           domain: [xs_g,ys_g,zs_g] = [%lld,%lld,%lld]-[%lld,%lld,%lld] \n",(LLD)xs_g,(LLD)ys_g,(LLD)zs_g,(LLD)(xs_g+xm_g-1),(LLD)(ys_g+ym_g-1),(LLD)(zs_g+zm_g-1));


			PetscPrintf(PETSC_COMM_SELF, "           minx,maxx=[%g,%g], miny,maxy=[%g,%g], minz,maxz=[%g,%g] \n",minx,maxx,miny,maxy,minz,maxz);
			PetscPrintf(PETSC_COMM_SELF, "           minx,maxx=[%g,%g], miny,maxy=[%g,%g], minz,maxz=[%g,%g] \n",max_leftx,min_rightx,max_fronty,min_backy,max_botz,min_topz);

			// perform the in-element test again
			PetscPrintf(PETSC_COMM_SELF, "           ParticleLocal: [ix,iy,iz] 	 = [%g,%g,%g] \n",ParticleLocal.ix,  ParticleLocal.iy, 		ParticleLocal.iz);
			PetscPrintf(PETSC_COMM_SELF, "           ParticleLocal: [x,y,z] 		 = [%g,%g,%g] \n",ParticleLocal.x,   ParticleLocal.y, 		ParticleLocal.z);
			PetscPrintf(PETSC_COMM_SELF, "           ParticleLocal: [eta,zetha,phi] = [%g,%g,%g] \n",ParticleLocal.eta, ParticleLocal.zetha, 	ParticleLocal.phi);

			//		eta = ParticleLocal.eta; zetha = ParticleLocal.zetha; phi = ParticleLocal.phi;
			Inside_8NodeElement(ParticleLocal.x, ParticleLocal.y, ParticleLocal.z, &eta, &zetha, &phi, coord_elem,	&inside, 0);
			PetscPrintf(PETSC_COMM_SELF, "           Computed: [eta,zetha,phi] = [%g,%g,%g] inside=%lld \n",eta, zetha, 	phi, (LLD)inside);

			PetscPrintf(PETSC_COMM_SELF," x-coords of local element: [%g,%g,%g,%g,%g,%g,%g,%g] \n",coord_elem[0].x,coord_elem[1].x,coord_elem[2].x,coord_elem[3].x,coord_elem[4].x,coord_elem[5].x,coord_elem[6].x,coord_elem[7].x);
			PetscPrintf(PETSC_COMM_SELF," y-coords of local element: [%g,%g,%g,%g,%g,%g,%g,%g] \n",coord_elem[0].y,coord_elem[1].y,coord_elem[2].y,coord_elem[3].y,coord_elem[4].y,coord_elem[5].y,coord_elem[6].y,coord_elem[7].y);
			PetscPrintf(PETSC_COMM_SELF," z-coords of local element: [%g,%g,%g,%g,%g,%g,%g,%g] \n",coord_elem[0].z,coord_elem[1].z,coord_elem[2].z,coord_elem[3].z,coord_elem[4].z,coord_elem[5].z,coord_elem[6].z,coord_elem[7].z);



		}


		ix = ix_new; iy=iy_new; iz=iz_new;

		// Check that we don't go outside the ghostpoints with the new element.
		// if this occurs make an emergency break and delete the tracer (which is
		// a hack that could probably be improved)
		if ( (ix<xs_g) || (ix>(xs_g+xm_g-2)) || (iy<ys_g) || (iy>(ys_g+ym_g-2)) || (iz<zs_g) || (iz>(zs_g+zm_g-2))   ){
			TracerDeactivate = 1;	ElementFound      = 1;  	// emergency break
			//PetscPrintf(PETSC_COMM_SELF," deactivated tracer %g \n",ParticleLocal.num);

		}

		n_iter = n_iter+1;


	}	// end



	if (PeriodicCorrection==0){
		/* Check to which CPU the particle should be send */
		if ((ix<xs))				{	cpu_x = cpu_x-1;	}
		else if ((ix>(xs+xm-1)))	{	cpu_x = cpu_x+1;	}
		if ((iy<ys) )				{	cpu_y = cpu_y-1;	}
		else if ((iy>(ys+ym-1)))	{	cpu_y = cpu_y+1;	}
	}
	if (iz<zs)					{	cpu_z = cpu_z-1;	}
	else if (iz>(zs+zm-1)) 		{	cpu_z = cpu_z+1;	}

	ParticleLocal.ix  = ((double) ix);
	ParticleLocal.iy  = ((double) iy);
	ParticleLocal.iz  = ((double) iz);


	*ElementFound_Outside	=	ElementFound;
	*ix_in 					=	ix;
	*iy_in 					=	iy;
	*iz_in 					=	iz;
	*cpu_x_in				=	cpu_x;
	*cpu_y_in				=	cpu_y;
	*cpu_z_in				=	cpu_z;

	*TracerDeactivate_in	=	TracerDeactivate;
	*ParticleLocal_in		=	ParticleLocal;


	PetscFunctionReturn(0);
}
/* ==================================================================================================== */




/* ====================================================================================================
 * Advect particles 																					*/
#undef __FUNCT__
#define __FUNCT__ "AdvectParticles"
PetscErrorCode AdvectParticles( LaMEMVelPressureDA C, DM da, UserContext *user, Vec Velocity, PetscScalar dt )
{
	PetscMPIInt	rank, size;
	PetscErrorCode ierr;
	PetscInt	ipart, ix, iy, iz, inode;
	PetscInt	ElementType_loc, nnode_z;
	PetscInt	xs,ys,zs,xm,ym,zm;
	DM			cda;
	Vec			local_Vel, gc, global;
	DMDACoor3d	***coords, coord_elem[MAX_nnel], coord_elem_new[MAX_nnel];
	DMDACoor3d 	CoordIntp;
	Particles	ParticleLocal;
	Field 		***velocity;
	PetscScalar	V_element[MAX_edof], Point[3], **dhdsVel;
	double      ShapeVel[MAX_nnel];
	PetscInt		mod;
	PetscInt nnel;
	DAVPElementType vpt_element_type;

	vpt_element_type = C->type;
	nnel = C->nnel;

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&size);


	ierr = LaMEMCreate2dArray( 3, nnel, &dhdsVel, PETSC_NULL ); CHKERRQ(ierr);

	/* coordinates */
	ierr = DMGetCoordinateDM(da,&cda); CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(da,&gc); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,gc,&coords); CHKERRQ(ierr);

	/* velocity */
	ierr = DMGetLocalVector(da,&local_Vel); CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da, Velocity, INSERT_VALUES, local_Vel); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(da,   Velocity, INSERT_VALUES, local_Vel); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, local_Vel,	&velocity); CHKERRQ(ierr);

	/* Number of elements */
	ierr = DMDAGetInfo(da, 0, 0, 0, &nnode_z, 0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);


	/* Loop over all particles */
	ix =0; iy=0; iz=0;
	for (ipart=0; ipart<user->num_particle_local; ipart++){

		/* Extract coordinates of element */
		ParticleLocal	=	user->ParticlesLocal[ipart];

		/* Only advect active particles */
		if (ParticleLocal.phase >= 0){		// only if particle is active


			if( (vpt_element_type == DAVP_Q1P0) || (vpt_element_type == DAVP_Q1Q1) || (vpt_element_type == DAVP_FDSTAG)) {
				ix				=	((PetscInt) ParticleLocal.ix);
				iy				=	((PetscInt) ParticleLocal.iy);
				iz				=	((PetscInt) ParticleLocal.iz);
			}
			else if ( (vpt_element_type==DAVP_Q2PM1L) || (vpt_element_type==DAVP_Q2PM1G) ) {
				ix = ((PetscInt) ParticleLocal.ix); LaMEMMod(ix, 2, &mod); if (mod>0){ix = ix-1;};
				iy = ((PetscInt) ParticleLocal.iy); LaMEMMod(iy, 2, &mod); if (mod>0){iy = iy-1;};
				iz = ((PetscInt) ParticleLocal.iz); LaMEMMod(iz, 2, &mod); if (mod>0){iz = iz-1;};
			}
			else if (vpt_element_type == DAVP_FDSTAG) {
				SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Routine does not yet work for FDSTAG!" );
				//	ierr = MPI_Abort(PETSC_COMM_WORLD,1);  CHKERRQ(ierr);
			}

			/* Get type of element */
			ElementType_loc=0;
			if (user->GridAdvectionMethod==1) 				{	ElementType_loc=2;}	// lagrangian element
			if 	(user->GridAdvectionMethod==2){		// ALE mode
				if (iz> (nnode_z-user->NumSurfaceNodes) )	{	ElementType_loc=2;}	// close-to-surface elements (lagrangian)
				if (iz==(nnode_z-user->NumSurfaceNodes) ) 	{	ElementType_loc=1;}	// transition zone (intermediate element)
			}

			if ( (ix>=xs) && (ix<(xs+xm)) &&  (iy>=ys) && (iy<(ys+ym)) &&  (iz>=zs) && (iz<(zs+zm))){

				/* Verify that element is on current CPU */
				// re-construct element coordinates
				ierr = GetElementCoords(coord_elem, coords, ix,iy,iz, 1);					CHKERRQ(ierr);
				ierr = CorrectElementCoordsForPeriodicity( coord_elem, ix, iy, user, 1);	CHKERRQ(ierr);


				/* Extract element velocity */
				ierr = GetVelocityElement(velocity, V_element, ix,iy,iz);	CHKERRQ(ierr);

				/* Modify velocity, if we have a lagrangian or a half-lagrangian element */
				if (ElementType_loc>0){
					if (ElementType_loc==2){
						for (inode=0; inode<nnel; inode++){
							V_element[3*inode + 0] = 0;
							V_element[3*inode + 1] = 0;
							V_element[3*inode + 2] = 0;
						}
					}
					if (ElementType_loc==1){
						// intermediate layer (only the bottom part of the element is advected,
						//	since the rest of the element has been advected already).

						if( (vpt_element_type == DAVP_Q1P0) || (vpt_element_type == DAVP_Q1Q1) || (vpt_element_type == DAVP_FDSTAG) ) {
							for (inode=4; inode<nnel; inode++){
								V_element[3*inode + 0] = 0;
								V_element[3*inode + 1] = 0;
								V_element[3*inode + 2] = 0;
							}
						}

						if( (vpt_element_type==DAVP_Q2PM1L) || (vpt_element_type==DAVP_Q2PM1G) ) {
							for (inode=4;  inode<8; inode++){
								V_element[3*inode + 0] = 0;	V_element[3*inode + 1] = 0;	V_element[3*inode + 2] = 0;
							}
							for (inode=12; inode<20; inode++){
								V_element[3*inode + 0] = 0;	V_element[3*inode + 1] = 0;	V_element[3*inode + 2] = 0;
							}
							for (inode=21; inode<27; inode++){
								V_element[3*inode + 0] = 0;	V_element[3*inode + 1] = 0;	V_element[3*inode + 2] = 0;
							}
						}

					}
				}


				/* Compute new element shape */
				for (inode=0; inode<nnel; inode++){
					coord_elem_new[inode].x = coord_elem[inode].x + dt*V_element[3*inode + 0];
					coord_elem_new[inode].y = coord_elem[inode].y + dt*V_element[3*inode + 1];
					coord_elem_new[inode].z = coord_elem[inode].z + dt*V_element[3*inode + 2];
				}

				/* Compute new tracer coordinate from new element shape */
				Point[0] = ParticleLocal.eta  ;
				Point[1] = ParticleLocal.zetha;
				Point[2] = ParticleLocal.phi  ;

				ierr = ComputeShapeFunctionVelocity(ShapeVel, dhdsVel, Point);			CHKERRQ(ierr);		// Velocity shape function
				ierr = ComputeCoordIntp(nnel,ShapeVel, coord_elem_new,  &CoordIntp);			CHKERRQ(ierr);// Real coordinate of particle

				/* Put back particle */
				user->ParticlesLocal[ipart].x	=	CoordIntp.x;
				user->ParticlesLocal[ipart].y	=	CoordIntp.y;
				user->ParticlesLocal[ipart].z	=	CoordIntp.z;
			}


		}	// end of only advect active particles

	} // end of particles loo[p
	ierr = DMDAVecRestoreArray(cda,gc,&coords); CHKERRQ(ierr);
	ierr = DMGetCoordinates(da,&global); CHKERRQ(ierr);

	ierr = DMDAVecRestoreArray(da,local_Vel,&velocity); CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(da,&local_Vel); CHKERRQ(ierr);


	ierr = LaMEMDestroy2dArray(&dhdsVel, PETSC_NULL ); CHKERRQ(ierr);

	//	ierr = DMDestroy( cda ); CHKERRQ(ierr);
	//	ierr = VecDestroy( gc ); CHKERRQ(ierr);
	//	ierr = VecDestroy( global ); CHKERRQ(ierr);

	PetscFunctionReturn(0);

}
/* ==================================================================================================== */

/* ====================================================================================================
 * Interpolate particles 																					*/
#undef __FUNCT__
#define __FUNCT__ "InterpolateParticles"
PetscErrorCode InterpolateParticles( LaMEMVelPressureDA C, DM da, UserContext *user)
{
	PetscMPIInt	rank, size;
	PetscInt	ipart, ix, iy, iz, inode;
	PetscInt	nnode_z;
	PetscInt	xs,ys,zs,xm,ym,zm;
	DM			cda;
	Vec			gc, global;
	DMDACoor3d	***coords, coord_elem[MAX_nnel], coord_elem_new[MAX_nnel];
	DMDACoor3d 	CoordIntp;
	Particles	ParticleLocal;
	PetscScalar	Point[3], **dhdsVel;
	double      ShapeVel[MAX_nnel];
	PetscInt		mod;
	PetscInt nnel;
	DAVPElementType vpt_element_type;

	vpt_element_type = C->type;
	nnel = C->nnel;

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&size);


	LaMEMCreate2dArray( 3, nnel, &dhdsVel, PETSC_NULL );

	/* coordinates */
	DMGetCoordinateDM(da,&cda);
	DMGetCoordinatesLocal(da,&gc);
	DMDAVecGetArray(cda,gc,&coords);

	/* Number of elements */
	DMDAGetInfo(da, 0, 0, 0, &nnode_z, 0,0,0,0,0,0,0,0,0);
	DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);


	/* Loop over all particles */
	ix=0; iy=0; iz=0;
	for (ipart=0; ipart<user->num_particle_local; ipart++){

		/* Extract coordinates of element */
		ParticleLocal	=	user->ParticlesLocal[ipart];

		/* Only advect active particles */
		if (ParticleLocal.phase>-4){		// only if particle is active

			if( (vpt_element_type == DAVP_Q1P0) || (vpt_element_type == DAVP_Q1Q1) || (vpt_element_type == DAVP_FDSTAG)) {
				ix				=	((PetscInt) ParticleLocal.ix);
				iy				=	((PetscInt) ParticleLocal.iy);
				iz				=	((PetscInt) ParticleLocal.iz);
			}
			else if( (vpt_element_type==DAVP_Q2PM1L) || (vpt_element_type==DAVP_Q2PM1G) ) {
				ix = ((PetscInt) ParticleLocal.ix); LaMEMMod(ix, 2, &mod); if (mod>0){ix = ix-1;};
				iy = ((PetscInt) ParticleLocal.iy); LaMEMMod(iy, 2, &mod); if (mod>0){iy = iy-1;};
				iz = ((PetscInt) ParticleLocal.iz); LaMEMMod(iz, 2, &mod); if (mod>0){iz = iz-1;};
			}

			if ( (ix>=xs) && (ix<(xs+xm)) &&  (iy>=ys) && (iy<(ys+ym)) &&  (iz>=zs) && (iz<(zs+zm))){

				/* Verify that element is on current CPU */
				// re-construct element coordinates
				GetElementCoords(coord_elem, coords, ix,iy,iz, 1);
				CorrectElementCoordsForPeriodicity( coord_elem, ix, iy, user, 1);

				/* Compute new element shape */
				for (inode=0; inode<nnel; inode++){
					coord_elem_new[inode].x = coord_elem[inode].x; // + dt*V_element[3*inode + 0];
					coord_elem_new[inode].y = coord_elem[inode].y; // + dt*V_element[3*inode + 1];
					coord_elem_new[inode].z = coord_elem[inode].z; // + dt*V_element[3*inode + 2];
				}

				/* Compute new tracer coordinate from new element shape */
				Point[0] = ParticleLocal.eta  ;
				Point[1] = ParticleLocal.zetha;
				Point[2] = ParticleLocal.phi  ;

				ComputeShapeFunctionVelocity(ShapeVel, dhdsVel, Point);				// Velocity shape function
				ComputeCoordIntp(nnel,ShapeVel, coord_elem_new,  &CoordIntp);		// Real coordinate of particle

				/* Put back particle */
				user->ParticlesLocal[ipart].x	=	CoordIntp.x;
				user->ParticlesLocal[ipart].y	=	CoordIntp.y;
				user->ParticlesLocal[ipart].z	=	CoordIntp.z;
			}

		}	// end of only advect active particles

	} // end of particles loo[p
	DMDAVecRestoreArray(cda,gc,&coords);
	DMGetCoordinates(da,&global);

	LaMEMDestroy2dArray(&dhdsVel, PETSC_NULL );

	//	DMDestroy( cda );
	//	VecDestroy( gc );
	//	VecDestroy( global );

	PetscFunctionReturn(0);

}
/* ==================================================================================================== */

/* ====================================================================================================
 * Set initial tracer phases 					*/
//DONE
#undef __FUNCT__
#define __FUNCT__ "SetInitialTracerProperties"
PetscErrorCode SetInitialTracerProperties( DM da, UserContext *user )
{
	PetscInt		ipart;
	Particles		ParticleLocal;
	PetscErrorCode 	ierr;

	/* Loop over all particles */




	/* Set Phases depending on the case of Model Setup */
	if 		(user->Setup.Model==0){				// diapir setup

		for (ipart=0; ipart<user->num_particle_local; ipart++){
			ParticleLocal 		= user->ParticlesLocal[ipart];
			ParticleLocal.phase = 0;

			//PetscPrintf(PETSC_COMM_WORLD," *** DIAPIR SETUP: Setup.Diapir_Hi=%g ParticleLocal.z=%g *** \n",user->Setup.Diapir_Hi,ParticleLocal.z);

			if ( ParticleLocal.z>=user->Setup.Diapir_Hi){
				ParticleLocal.phase = 1;
			}



			user->ParticlesLocal[ipart]	=	ParticleLocal;
		}


	}
	else if (user->Setup.Model==1){ 			// single layer folding setup

	}
	else if (user->Setup.Model==2){ 			//
		// Falling block setup
		PetscInt	 	nx,ny,nz, nel_x, nel_y, nel_z;
		PetscScalar 	dx,dy,dz;
		PetscScalar 	BlockLeft, BlockRight, BlockFront, BlockBack, BlockBottom, BlockTop;
		PetscScalar		BlockWidthX,BlockWidthY,BlockWidthZ;
		PetscBool 		Block_2D, Block_2Dy;

		// number of elements on finest resolution
		nel_x 				= 	user->finest_nelx;
		nel_y 				= 	user->finest_nely;
		nel_z 				= 	user->finest_nelz;


		ierr 				= 	DMDAGetInfo(da,0,&nx,&ny,&nz,0,0,0,0,0,0,0,0,0);	CHKERRQ(ierr);  	// # of nodes in all directions
		dx 					= 	user->W/((double) (nel_x));
		dy 					= 	user->L/((double) (nel_y));
		dz 					= 	user->H/((double) (nel_z));

		BlockWidthX 		=	((double) (nel_x/2))*dx;
		BlockWidthY 		=	((double) (nel_y/2))*dy;
		BlockWidthZ 		=	((double) (nel_z/2))*dz;


		// DEBUGGING!!!
		/*
		BlockWidthX = user->L;
		BlockWidthZ = user->H;
		BlockWidthY = user->W;

		BlockFront   =  user->y_front;
		//BlockFront = user->W/2;
		BlockLeft = user->L/2;
		//BlockLeft = user->x_left;

		BlockBottom  = user->z_bot;
		//BlockBottom  = user->H/2;
		 */


		BlockLeft 			=	((double) (1*nel_x/4))*dx + user->x_left;		//	left side of block
		BlockRight 			=	BlockLeft + BlockWidthX;					//	right side of block
		BlockFront 			=	((double) (1*nel_y/4))*dy + user->y_front;		//	front side of block
		BlockBack 			=	BlockFront + BlockWidthY;					//	back side of block
		BlockBottom 		=	((double) (1*nel_z/4))*dz + user->z_bot;			//	bottom side of block
		BlockTop 			=	BlockBottom + BlockWidthZ;					//	top side of block



		// Print information about the falling block setup:
		PetscPrintf(PETSC_COMM_WORLD," *** Using a falling block setup *** \n");
		PetscPrintf(PETSC_COMM_WORLD," Falling block coordinates:  [left,right]=[%f,%f]; [front,back]=[%f,%f]; [bottom,top]=[%f,%f] \n",
				BlockLeft,BlockRight,BlockFront,BlockBack,BlockBottom,BlockTop,nel_x);

		Block_2D  = PETSC_FALSE;
		Block_2Dy = PETSC_FALSE;
		ierr 	 = PetscOptionsGetBool( PETSC_NULL, "-Block_2D",   &Block_2D, PETSC_NULL); CHKERRQ(ierr);
		ierr 	 = PetscOptionsGetBool( PETSC_NULL, "-Block_2Dy",  &Block_2Dy, PETSC_NULL); CHKERRQ(ierr);

		for (ipart=0; ipart<user->num_particle_local; ipart++){

			ParticleLocal 		= user->ParticlesLocal[ipart];
			ParticleLocal.phase = 0;
			ParticleLocal.T 	= ParticleLocal.z;




			if ( 	(ParticleLocal.z>BlockBottom) 	&&
					(ParticleLocal.z<BlockTop) 		&&
					(ParticleLocal.x>BlockLeft) 	&&
					(ParticleLocal.x<BlockRight)	&& (Block_2Dy==PETSC_FALSE) ){


				if (Block_2D == PETSC_TRUE){
					// 2D block
					ParticleLocal.phase = 1;	// block
					ParticleLocal.T 	= 1.0;
				}
				else{
					if ((ParticleLocal.y>BlockFront) 	&&
							(ParticleLocal.y<BlockBack) ){

						// 3D block
						ParticleLocal.phase 	= 1;	// block
						ParticleLocal.T 	= 1.0;

					}
				}

			}

			if (Block_2Dy == PETSC_TRUE){
				// Create 2D block in y-direction
				if ( 	(ParticleLocal.z>BlockBottom) 	&&
						(ParticleLocal.z<BlockTop) 		&&
						(ParticleLocal.y>BlockFront) 	&&
						(ParticleLocal.y<BlockBack) ){

					// 2D block
					ParticleLocal.phase 	= 1;	// block
					ParticleLocal.T 	= 1.0;


				}
			}



			user->ParticlesLocal[ipart]	=	ParticleLocal;
		}

	}
	else if (user->Setup.Model==6){ 			// subduction setup with air
		PetscScalar 	H_air, SlabThickness, SlabWidth, SlabLength;
		PetscScalar		DistanceFromLeft, SlabMaxSubdDepth,Slab_Xright, dz;
		PetscScalar		Air_thickness, Slab_ThicknessFactor, Air_ThicknessFactor, Slab_WidthFactor;
		PetscInt		nx,ny,nz;

		ierr 				= DMDAGetInfo(da,0,&nx,&ny,&nz,0,0,0,0,0,0,0,0,0);	CHKERRQ(ierr);  	// # of nodes in all directions
		dz 					= user->H/((double) (nz-1));


		Slab_ThicknessFactor = 0.1;	// in % of the total H
		PetscOptionsGetReal(PETSC_NULL ,"-Slab_ThicknessFactor",	&Slab_ThicknessFactor	, PETSC_NULL);

		Air_ThicknessFactor = 0.1;	// in % of the total H
		PetscOptionsGetReal(PETSC_NULL ,"-Air_ThicknessFactor",		&Air_ThicknessFactor		, PETSC_NULL);

		Slab_WidthFactor = 1.0;
		PetscOptionsGetReal(PETSC_NULL ,"-Slab_WidthFactor",		&Slab_WidthFactor		, PETSC_NULL);


		//Air_thickness 		=	round((user->H/15)/dz)*dz*Slab_AirThicknessFactor;

		//SlabThickness 		= 	round((user->H/10)/dz)*dz;
		SlabThickness 		= 	round(user->H/dz)*dz*Slab_ThicknessFactor;
		Air_thickness 		= 	round(user->H/dz)*dz*Air_ThicknessFactor;
		H_air 	      		= 	user->z_bot + user->H - Air_thickness;


		SlabWidth 	  		= 	Slab_WidthFactor*user->L;
		SlabLength    		=	0.4*(user->W);
		DistanceFromLeft 	=	(user->W)/10;
		SlabMaxSubdDepth 	=	2*SlabThickness;
		Slab_Xright 		=	DistanceFromLeft + SlabLength-(user->H-H_air);


		PetscPrintf(PETSC_COMM_WORLD,"SlabThickness [-Slab_ThicknessFactor]: %f AirThickness [-Air_ThicknessFactor]: %f H: %f SlabWidthFactor [-Slab_WidthFactor]: %f \n",SlabThickness, Air_thickness, user->H, Slab_WidthFactor);

		//PetscPrintf(PETSC_COMM_WORLD,"Air: %f %f \n",Air_thickness, dz);


		for (ipart=0; ipart<user->num_particle_local; ipart++){
			ParticleLocal 		= user->ParticlesLocal[ipart];
			ParticleLocal.phase = 0;

			if ( 	(ParticleLocal.z>(H_air-SlabThickness)) &&
					(ParticleLocal.z<H_air) &&
					(ParticleLocal.x>(user->x_left+DistanceFromLeft)) &&
					(ParticleLocal.x<(user->x_left+DistanceFromLeft+SlabLength)) &&
					(ParticleLocal.y<(user->y_front+SlabWidth)) ){

				ParticleLocal.phase 	= 1;	// Slab top part

			}

			if ( 	(ParticleLocal.z<(user->H - (ParticleLocal.x-Slab_Xright))) &&
					(ParticleLocal.z>(user->H-2*SlabThickness - (ParticleLocal.x-Slab_Xright))) &&
					(ParticleLocal.y<(user->y_front+SlabWidth)) &&
					(ParticleLocal.z>(H_air-SlabMaxSubdDepth)) ){

				ParticleLocal.phase 	= 1;	// Slab inclined part

			}
			if ( ParticleLocal.z>H_air){
				ParticleLocal.phase = 2;	// Air
			}

			user->ParticlesLocal[ipart]	=	ParticleLocal;
		}

	}
	else if (user->Setup.Model==5){

		for (ipart=0; ipart<user->num_particle_local; ipart++){

			ParticleLocal 		= user->ParticlesLocal[ipart];
			ParticleLocal.phase = 0;
			if ( (((PetscInt) ParticleLocal.iz)<=(user->Setup.ind_fold_top)) && (((PetscInt) ParticleLocal.iz)>=(user->Setup.ind_fold_bot)) ){
				ParticleLocal.phase = 1;
			}
			user->ParticlesLocal[ipart]	=	ParticleLocal;
		}
	}
	else if (user->Setup.Model==8){
		// a square heterogeneity is present with specified dimensions
		PetscScalar x_start, x_end, y_start, y_end, z_start, z_end;

		user->num_phases = 2;

		// (1) Extract the correct lengthscales from the commandline
		x_start = 0.0; x_end = 0.0;
		y_start = 0.0; y_end = 0.0;
		z_start = 0.0; z_end = 0.0;

		PetscOptionsGetReal(PETSC_NULL,"-Heterogeneity_x_start"      ,	&x_start 		, PETSC_NULL);
		PetscOptionsGetReal(PETSC_NULL,"-Heterogeneity_x_end"        ,	&x_end 			, PETSC_NULL);
		PetscOptionsGetReal(PETSC_NULL,"-Heterogeneity_y_start"      ,	&y_start 		, PETSC_NULL);
		PetscOptionsGetReal(PETSC_NULL,"-Heterogeneity_y_end"        ,	&y_end 			, PETSC_NULL);
		PetscOptionsGetReal(PETSC_NULL,"-Heterogeneity_z_start"      ,	&z_start 		, PETSC_NULL);
		PetscOptionsGetReal(PETSC_NULL,"-Heterogeneity_z_end"      ,	&z_end 			, PETSC_NULL);


		// (2) Non-dimensionalize the values
		x_start = x_start/user->Characteristic.Length;
		x_end   = x_end/user->Characteristic.Length;
		y_start = y_start/user->Characteristic.Length;
		y_end   = y_end/user->Characteristic.Length;
		z_start = z_start/user->Characteristic.Length;
		z_end   = z_end/user->Characteristic.Length;

		for (ipart=0; ipart<user->num_particle_local; ipart++){
			ParticleLocal 		= user->ParticlesLocal[ipart];
			ParticleLocal.phase = 0;

			if ( 	(ParticleLocal.x>=x_start) &&
					(ParticleLocal.x<=x_end  ) &&
					(ParticleLocal.y>=y_start) &&
					(ParticleLocal.y<=y_end  ) &&
					(ParticleLocal.z>=z_start) &&
					(ParticleLocal.z<=z_end  ) ) {

				ParticleLocal.phase = 1;
			}

			user->ParticlesLocal[ipart]	=	ParticleLocal;

		}
	}

	else if (user->Setup.Model==10){
			//A one layer setup over a detachment with two linear shaped perturbation (as in Grasemann & Schmalholz 2012)
			//One perturbation is fixed in position, the other one's position is varied with -Heterogeneity_Offset

			PetscScalar 		z_bot, z_top, dz_nodes, Het_L, Het_W, Het_A, Het_Off;
			PetscBool     	flg, DisplayInfo;
			DisplayInfo = PETSC_FALSE;


			PetscOptionsGetReal(PETSC_NULL,"-Heterogeneity_L"      ,	&Het_L 		, PETSC_NULL);
			PetscOptionsGetReal(PETSC_NULL,"-Heterogeneity_W"      ,	&Het_W 			, PETSC_NULL);
			PetscOptionsGetReal(PETSC_NULL,"-Heterogeneity_A"      ,	&Het_A 		, PETSC_NULL);
			PetscOptionsGetReal(PETSC_NULL,"-Heterogeneity_Offset" ,	&Het_Off 			, PETSC_NULL);

			// Non-dimensionalize the values
			Het_L 		= Het_L/user->Characteristic.Length;
			Het_W 		= Het_W/user->Characteristic.Length;
			Het_A 		= Het_A/user->Characteristic.Length;
			Het_Off 	= Het_Off/user->Characteristic.Length;

			dz_nodes = user->H/((PetscScalar)(user->nel_z) * 1.0);
			for (ipart=0; ipart<user->num_particle_local; ipart++){
							ParticleLocal 		= user->ParticlesLocal[ipart];
							ParticleLocal.phase = 0;

							PetscOptionsGetReal(PETSC_NULL,"-Layer1_bottom"      ,	&z_bot 		, PETSC_NULL);
							PetscOptionsGetReal(PETSC_NULL,"-Layer1_top"         ,	&z_top 		, &flg);

							//One layer properties
							if ( (ParticleLocal.iz>=(z_bot/dz_nodes) ) && (ParticleLocal.iz<(z_top/dz_nodes)  ) && flg ) {
								ParticleLocal.phase = 1;
								if (DisplayInfo){
									PetscPrintf(PETSC_COMM_WORLD,"Adding Layer with Phase=1 from z=[%g,%g]\n",z_bot, z_top);
								}
							}

							//Seed Fixed
							if ( 	(ParticleLocal.y>= (user->L/2 - Het_L/2) ) &&
									(ParticleLocal.y<= (user->L/2 + Het_L/2)  ) &&
									(ParticleLocal.x<= (Het_W)) &&
									(ParticleLocal.iz<=(z_bot + Het_A)/dz_nodes  ) 	) {
									ParticleLocal.phase = 0;
							}

							//Seed with Offset
							if ( 	(ParticleLocal.y>= (Het_Off + user->L/2 - Het_L/2) ) &&
									(ParticleLocal.y<= (Het_Off + user->L/2 + Het_L/2)  ) &&
									(ParticleLocal.x>= (user->W - Het_W)) &&
									(ParticleLocal.iz<=(z_bot + Het_A)/dz_nodes  ) 	) {
									ParticleLocal.phase = 0;
							}


			user->ParticlesLocal[ipart]	=	ParticleLocal;
			}


	}



	else if (user->Setup.Model==9){
		// A multilayer folding setup that is relevant for the Zagros setup
		PetscScalar 		zbot[10], ztop[10], dz_nodes;
		PetscBool     	flg, DisplayInfo;
		PetscInt		PertType,i;
		//		PetscScalar 	dz_nodes;


		DisplayInfo = PETSC_FALSE;

		//Adapt the particles to the initial mesh deformed with a cos perturbation (ampl2D)
		//if((user->ampl2D)!=(0.0)){

		PertType = 0;	// defauly value
		PetscOptionsGetInt(PETSC_NULL,"-Perturbation"		, &PertType	, PETSC_NULL);

		// initialize arrays for fractions
		for (i=0;i<10;i++) { zbot[i] = 0.0; ztop[i] = 0.0;}

		// get input data
		PetscOptionsGetReal(PETSC_NULL,"-Layer1_bottom"      ,&zbot[0], PETSC_NULL);
		PetscOptionsGetReal(PETSC_NULL,"-Layer1_top"         ,&ztop[0], &flg);
		if (DisplayInfo && flg) PetscPrintf(PETSC_COMM_WORLD,"# - Adding Layer with Phase=1 from z=[%g,%g]\n",zbot[0],ztop[0]);

		PetscOptionsGetReal(PETSC_NULL,"-Layer2_bottom"      ,&zbot[1], PETSC_NULL);
		PetscOptionsGetReal(PETSC_NULL,"-Layer2_top"         ,&ztop[1], &flg);
		if (DisplayInfo && flg) PetscPrintf(PETSC_COMM_WORLD,"# - Adding Layer with Phase=2 from z=[%g,%g]\n",zbot[1],ztop[1]);

		PetscOptionsGetReal(PETSC_NULL,"-Layer3_bottom"      ,&zbot[2], PETSC_NULL);
		PetscOptionsGetReal(PETSC_NULL,"-Layer3_top"         ,&ztop[2], &flg);
		if (DisplayInfo && flg) PetscPrintf(PETSC_COMM_WORLD,"# - Adding Layer with Phase=3 from z=[%g,%g]\n",zbot[2],ztop[2]);

		PetscOptionsGetReal(PETSC_NULL,"-Layer4_bottom"      ,&zbot[3], PETSC_NULL);
		PetscOptionsGetReal(PETSC_NULL,"-Layer4_top"         ,&ztop[3], &flg);
		if (DisplayInfo && flg) PetscPrintf(PETSC_COMM_WORLD,"# - Adding Layer with Phase=4 from z=[%g,%g]\n",zbot[3],ztop[3]);

		PetscOptionsGetReal(PETSC_NULL,"-Layer5_bottom"      ,&zbot[4], PETSC_NULL);
		PetscOptionsGetReal(PETSC_NULL,"-Layer5_top"         ,&ztop[4], &flg);
		if (DisplayInfo && flg) PetscPrintf(PETSC_COMM_WORLD,"# - Adding Layer with Phase=5 from z=[%g,%g]\n",zbot[4],ztop[4]);

		PetscOptionsGetReal(PETSC_NULL,"-Layer6_bottom"      ,&zbot[5], PETSC_NULL);
		PetscOptionsGetReal(PETSC_NULL,"-Layer6_top"         ,&ztop[5], &flg);
		if (DisplayInfo && flg) PetscPrintf(PETSC_COMM_WORLD,"# - Adding Layer with Phase=6 from z=[%g,%g]\n",zbot[5],ztop[5]);

		PetscOptionsGetReal(PETSC_NULL,"-Layer7_bottom"      ,&zbot[6], PETSC_NULL);
		PetscOptionsGetReal(PETSC_NULL,"-Layer7_top"         ,&ztop[6], &flg);
		if (DisplayInfo && flg) PetscPrintf(PETSC_COMM_WORLD,"# - Adding Layer with Phase=7 from z=[%g,%g]\n",zbot[6],ztop[6]);

		PetscOptionsGetReal(PETSC_NULL,"-Layer8_bottom"      ,&zbot[7], PETSC_NULL);
		PetscOptionsGetReal(PETSC_NULL,"-Layer8_top"         ,&ztop[7], &flg);
		if (DisplayInfo && flg) PetscPrintf(PETSC_COMM_WORLD,"# - Adding Layer with Phase=8 from z=[%g,%g]\n",zbot[7],ztop[7]);

		PetscOptionsGetReal(PETSC_NULL,"-Layer9_bottom"      ,&zbot[8], PETSC_NULL);
		PetscOptionsGetReal(PETSC_NULL,"-Layer9_top"         ,&ztop[8], &flg);
		if (DisplayInfo && flg) PetscPrintf(PETSC_COMM_WORLD,"# - Adding Layer with Phase=9 from z=[%g,%g]\n",zbot[8],ztop[8]);

		PetscOptionsGetReal(PETSC_NULL,"-Layer10_bottom"      ,&zbot[9], PETSC_NULL);
		PetscOptionsGetReal(PETSC_NULL,"-Layer10_top"         ,&ztop[9], &flg);
		if (DisplayInfo && flg) PetscPrintf(PETSC_COMM_WORLD,"# - Adding Layer with Phase=10 from z=[%g,%g]\n",zbot[9],ztop[9]);

		dz_nodes = user->H/((PetscScalar)(user->nel_z) * 2.0);

		if (PertType == 1){

			for (ipart=0; ipart<user->num_particle_local; ipart++){
				ParticleLocal 		= user->ParticlesLocal[ipart];
				ParticleLocal.phase = 0;

				if ( (ParticleLocal.iz>=(zbot[0]/dz_nodes) ) && (ParticleLocal.iz<(ztop[0]/dz_nodes)) ) {
					ParticleLocal.phase = 1;
				}

				if ( (ParticleLocal.iz>=(zbot[1]/dz_nodes) ) && (ParticleLocal.iz<(ztop[1]/dz_nodes)) ) {
					ParticleLocal.phase = 2;
				}

				if ( (ParticleLocal.iz>=(zbot[2]/dz_nodes) ) && (ParticleLocal.iz<(ztop[2]/dz_nodes)) ) {
					ParticleLocal.phase = 3;
				}

				if ( (ParticleLocal.iz>=(zbot[3]/dz_nodes) ) && (ParticleLocal.iz<(ztop[3]/dz_nodes)) ) {
					ParticleLocal.phase =	4;
				}

				if ( (ParticleLocal.iz>=(zbot[4]/dz_nodes) ) && (ParticleLocal.iz<(ztop[4]/dz_nodes)) ) {
					ParticleLocal.phase =	5;
				}

				if ( (ParticleLocal.iz>=(zbot[5]/dz_nodes) ) && (ParticleLocal.iz<(ztop[5]/dz_nodes)) ) {
					ParticleLocal.phase =	6;
				}

				if ( (ParticleLocal.iz>=(zbot[6]/dz_nodes) ) && (ParticleLocal.iz<(ztop[6]/dz_nodes)) ) {
					ParticleLocal.phase =	7;
				}

				if ( (ParticleLocal.iz>=(zbot[7]/dz_nodes) ) && (ParticleLocal.iz<(ztop[7]/dz_nodes)) ) {
					ParticleLocal.phase =	8;
				}

				if ( (ParticleLocal.iz>=(zbot[8]/dz_nodes) ) && (ParticleLocal.iz<(ztop[8]/dz_nodes)) ) {
					ParticleLocal.phase =	9;
				}

				if ( (ParticleLocal.iz>=(zbot[9]/dz_nodes) ) && (ParticleLocal.iz<(ztop[9]/dz_nodes)) ) {
					ParticleLocal.phase =	10;
				}

				user->ParticlesLocal[ipart]	=	ParticleLocal;

			}


		}

		else {

			for (ipart=0; ipart<user->num_particle_local; ipart++){
				ParticleLocal 		= user->ParticlesLocal[ipart];
				ParticleLocal.phase = 0;

				if ( (ParticleLocal.z>=zbot[0]  ) && (ParticleLocal.z<=ztop[0]  )) {
					ParticleLocal.phase = 1;
				}

				if ( (ParticleLocal.z>=zbot[1]  ) && (ParticleLocal.z<=ztop[1]  )) {
					ParticleLocal.phase = 2;
				}

				if ( (ParticleLocal.z>=zbot[2]  ) && (ParticleLocal.z<=ztop[2]  )) {
					ParticleLocal.phase = 3;
				}

				if ( (ParticleLocal.z>=zbot[3]  ) && (ParticleLocal.z<=ztop[3]  )) {
					ParticleLocal.phase =4;
				}

				if ( (ParticleLocal.z>=zbot[4]  ) && (ParticleLocal.z<=ztop[4]  )) {
					ParticleLocal.phase =5;
				}

				if ( (ParticleLocal.z>=zbot[5]  ) && (ParticleLocal.z<=ztop[5]  )) {
					ParticleLocal.phase = 6;
				}

				if ( (ParticleLocal.z>=zbot[6]  ) && (ParticleLocal.z<=ztop[6]  )) {
					ParticleLocal.phase =7;
				}

				if ( (ParticleLocal.z>=zbot[7]  ) && (ParticleLocal.z<=ztop[7]  )) {
					ParticleLocal.phase =8;
				}

				if ( (ParticleLocal.z>=zbot[8]  ) && (ParticleLocal.z<=ztop[8]  )) {
					ParticleLocal.phase =9;
				}

				if ( (ParticleLocal.z>=zbot[9]  ) && (ParticleLocal.z<=ztop[9]  )) {
					ParticleLocal.phase =10;
				}

				user->ParticlesLocal[ipart]	=	ParticleLocal;


			}
		}
	}
	PetscFunctionReturn(0);
}
/* ==================================================================================================== */

/* ====================================================================================================
 * Set initial tracer phases from a file read from the disk				 								*/
#undef __FUNCT__
#define __FUNCT__ "SetInitialTracerPhasesFromFile"
PetscErrorCode SetInitialTracerPhasesFromFile( DM da, UserContext *user )
{
	PetscMPIInt		rank, size;
	PetscErrorCode  ierr;
	PetscInt		ipart, nz;
	Particles		ParticleLocal;
	PetscViewer		view_in;
	Vec				PhaseVec;
	PetscInt		vector_length, NumBgPz,NumBgPy,NumBgPx, maxvalue;
	PetscInt		***BackgroundParticles, i, j, px, py, pz, num, ix, iy, iz;
	PetscScalar		*data, ***Temperature, maxtemp, dBGPart_x, dBGPart_y, dBGPart_z;

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&size);


	/* Read the file from disk on proc 0 */
	PetscPrintf(PETSC_COMM_WORLD," Reading phase distribution from binary file: %s   \n", user->ParticleFilename);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, user->ParticleFilename,FILE_MODE_READ, &view_in); CHKERRQ(ierr);

	ierr = VecCreate(PETSC_COMM_SELF,&PhaseVec);	CHKERRQ(ierr);	// create vector [from petsc 3.2 onwards]	
	ierr = VecLoad(PhaseVec, view_in);				CHKERRQ(ierr);	// load vector
	ierr = PetscViewerDestroy(&view_in);			CHKERRQ(ierr);	// close binary
	MPI_Barrier(PETSC_COMM_WORLD);

	// Copy vector to local array
	ierr = VecGetArray(PhaseVec,&data); CHKERRQ(ierr);

	ierr = VecGetSize(PhaseVec,&vector_length); CHKERRQ(ierr);

	NumBgPz	=	(PetscInt)(data[0]);
	NumBgPy	=	(PetscInt)(data[1]);
	NumBgPx	=	(PetscInt)(data[2]);

	/* Allocate temporary 3D array that holds the particles */
	//BackgroundParticles = (PetscInt ***)malloc(sizeof(PetscInt **)*NumBgPz);
	//for (i=0; i<NumBgPz; i++){
	//	BackgroundParticles[i] =  (PetscInt **)malloc(sizeof(PetscInt *)*NumBgPy);
	//	for (j=0; j<NumBgPy; j++){
	//		BackgroundParticles[i][j] =  (PetscInt *)malloc(sizeof(PetscInt )*NumBgPx);
	//	}
	//}

	/* Allocate temporary 3D array that holds the particles */
	ierr = PetscMalloc(sizeof(PetscInt **)*(size_t)NumBgPz, &BackgroundParticles); CHKERRQ(ierr);
	for (i=0; i<NumBgPz; i++){
		ierr = PetscMalloc(sizeof(PetscInt *)*(size_t)NumBgPy, &BackgroundParticles[i]); CHKERRQ(ierr);
		for (j=0; j<NumBgPy; j++){
			ierr = PetscMalloc(sizeof(PetscInt )*(size_t)NumBgPx, &BackgroundParticles[i][j]); CHKERRQ(ierr);
		}
	}




	/* Allocate temporary 3D array that holds the temperature */
	ierr = PetscMalloc(sizeof(PetscScalar **)*(size_t)NumBgPz, &Temperature); CHKERRQ(ierr);
	for (i=0; i<NumBgPz; i++){
		ierr = PetscMalloc(sizeof(PetscScalar *)*(size_t)NumBgPy, &Temperature[i]); CHKERRQ(ierr);
		for (j=0; j<NumBgPy; j++){
			ierr = PetscMalloc(sizeof(PetscScalar )*(size_t)NumBgPx, &Temperature[i][j]); CHKERRQ(ierr);
		}
	}

	/* Put particle phases into 3D array */
	num = 3;
	maxvalue = 0;
	for (pz=0; pz<NumBgPz; pz++){
		for (py=0; py<NumBgPy; py++){
			for (px=0; px<NumBgPx; px++){
				BackgroundParticles[pz][py][px] = (PetscInt)(data[num]);
				if (BackgroundParticles[pz][py][px]>maxvalue){
					maxvalue = BackgroundParticles[pz][py][px];
				}
				num = num+1;
			}
		}
	}

	/* Put temperature into 3D array */
	num = 3+NumBgPx*NumBgPy*NumBgPz;
	maxtemp = 0;
	for (pz=0; pz<NumBgPz; pz++){
		for (py=0; py<NumBgPy; py++){
			for (px=0; px<NumBgPx; px++){

				/* read in temperature and non-dimensionalize it if necessary
				 * Make sure that the Matlab file produces temperature in the correct
				 * 	dimensions!!!!
				 */
				Temperature[pz][py][px] =   data[num]/user->Characteristic.Temperature; ///Characteristic.Temperature;
				if (Temperature[pz][py][px]>maxtemp){
					maxtemp = Temperature[pz][py][px];
				}

				num = num+1;
			}
		}
	}

	PetscPrintf(PETSC_COMM_WORLD," Maximum phase-number detected: %lld   maxtemp = %g, num_phases = %lld \n", (LLD)maxvalue, maxtemp, (LLD)(user->num_phases+1));

	/* Error checking */
	if ( (maxvalue+1) > (user->num_phases)){
		PetscPrintf(PETSC_COMM_WORLD," ***** Detected %lld as maximum phase in particle file   *****\n", (LLD)maxvalue);
		PetscPrintf(PETSC_COMM_WORLD," ***** You however only specified material properties for phases [0-%lld] in input file   *****\n", (LLD)(user->num_phases));
		PetscPrintf(PETSC_COMM_WORLD," ***** Aborting. Change input file!   *****\n");
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Wrong maximum phase selected!");


		//ierr = MPI_Abort(PETSC_COMM_WORLD,1); CHKERRQ(ierr);
	}


	// Compute spacing of BG particle grid
	dBGPart_x = user->W/( (PetscScalar)(NumBgPx-1));
	dBGPart_y = user->L/( (PetscScalar)(NumBgPy-1));
	dBGPart_z = user->H/( (PetscScalar)(NumBgPz-1));

	ierr = VecRestoreArray(PhaseVec,&data);
	ierr = VecDestroy(&PhaseVec);



	/* Make a loop over all local particles and set the phases of those local particles
	 * 	based on the closest particle read from grid 									*/
	ierr = DMDAGetInfo(da,0,0,0,&nz,0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);
	for (ipart=0; ipart<user->num_particle_local; ipart++){
		ParticleLocal = user->ParticlesLocal[ipart];


		/* Find the BG particle phases that are closest to the current particle */
		ix = ((PetscInt)((ParticleLocal.x-user->x_left + 0.5*dBGPart_x)/dBGPart_x));
		iy = ((PetscInt)((ParticleLocal.y-user->y_front+ 0.5*dBGPart_y)/dBGPart_y));
		iz = ((PetscInt)((ParticleLocal.z-user->z_bot  + 0.5*dBGPart_z)/dBGPart_z));

		ParticleLocal.phase	=	((PetscScalar) BackgroundParticles[iz][iy][ix]);
		ParticleLocal.T		=	((PetscScalar) Temperature[iz][iy][ix]);

		/* Hack in case we use a lithosphere-following grid */
		if ( (user->InitialMeshFromFile == 1) && ( ((PetscInt) ParticleLocal.iz) >= (nz-user->InitialMantleLevel-1) ) ){
			ParticleLocal.phase	=	((PetscScalar) BackgroundParticles[NumBgPz-1][iy][ix]);
		}


		/* Set back into array */
		user->ParticlesLocal[ipart] = ParticleLocal;

	}


	/* Free memory */
	for (i=0; i<NumBgPz; i++){
		for (j=0; j<NumBgPy; j++){
			ierr = PetscFree(BackgroundParticles[i][j]); CHKERRQ(ierr);
		}
		ierr = PetscFree(BackgroundParticles[i]); CHKERRQ(ierr);
	}
	ierr = PetscFree(BackgroundParticles); CHKERRQ(ierr);

	for (i=0; i<NumBgPz; i++){
		for (j=0; j<NumBgPy; j++){
			ierr = PetscFree(Temperature[i][j]); CHKERRQ(ierr);
		}
		ierr = PetscFree(Temperature[i]); CHKERRQ(ierr);
	}
	ierr = PetscFree(Temperature); CHKERRQ(ierr);

	PetscFunctionReturn(0);
}
/* ===================================================================================================== */


/* Material properties from tracers @ finest grid ====================================================== */
#undef __FUNCT__
#define __FUNCT__ "MaterialPropertiesFromTracers"
PetscErrorCode MaterialPropertiesFromTracers( LaMEMVelPressureDA C, UserContext *user)
{
	PetscMPIInt             rank;
	PetscErrorCode          ierr;
	//MaterialsElement		***materials_fine;
	PetscInt				ipart, ielx, iely, ielz, i;
	PetscInt				xs  ,ys  ,zs  ,xm  ,ym  ,zm;
	PetscInt				xs_g,ys_g,zs_g,xm_g,ym_g,zm_g, el_number;
	PetscInt				intp, phase, MostAbundantPhase, iphase, MaxNumberTracers, TotalNumberOfParticlesInElement;
	PetscInt				**NumberParticlesOfEachPhase, DisplayOutput;
	PetscScalar				IntPoint[3][MAX_ngp_vel],  IntWeight[MAX_ngp_vel], distance, min_dist;
	PetscScalar				eta, zetha, phi, MeanValue;
	PetscScalar				MeanStrain,		MeanPlasticStrain, 	MeanE2nd,	MeanT2nd;
	PetscScalar				MeanPressure,	MeanTemperature, 	MeanPlastic;
	PetscScalar				MeanTxx, 		MeanTyy,			MeanTzz;
	PetscScalar				MeanTxz, 		MeanTyz,			MeanTxy,	MeanMuViscous;
	Particles				ParticleLocal;
	PetscLogDouble			cputime_start;
	PetscInt				mod;
	DAVPElementType 		vpt_element_type;
	PetscInt 					ngp_vel, nnel, nintp_1D;
	PetscScalar 			***materials_fine_array;
	MaterialsElementDynamic material_fine_data;

	vpt_element_type = C->type;
	ngp_vel         = C->ngp_vel;
	nnel            = C->nnel;
	nintp_1D        = C->nintp_1D;


	DisplayOutput = 0;	// 0-little output, 1-Lot's of output (for debugging)
	if (DisplayOutput==1){
		PetscPrintf(PETSC_COMM_WORLD,"# Computing material properties from tracers... \n");
	}

	// This routine computes the phase of an element by looping through all local particles.
	// A number for now - it should probably be changed to real properties such as viscosity @ a later stage

	PetscTime(&cputime_start);
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	/* Initialize some stuff */
	IntegrationPoints( nnel, nintp_1D, IntPoint,  IntWeight);		// integration points and weight

	/* Initialize all properties (set to zero) */
	ierr = DMDAVecGetArray(user->DA_Materials,   user->Materials, &materials_fine_array); CHKERRQ(ierr);
	ierr = DMDAGetCorners(user->DA_Processors,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
	for (ielz=zs; ielz<zs+zm; ielz++){
		for (iely=ys; iely<ys+ym; iely++){
			for (ielx=xs; ielx<xs+xm; ielx++){
				for (intp=0; intp<ngp_vel; intp++){
					LaMEMSetMaterialDataMemoryFromArray( &material_fine_data, ielx-xs,iely-ys,ielz-zs, C->ngp_vel, materials_fine_array );

					// We ONLY interpolate memory-variables from particles to integration points!!!
					material_fine_data.NumParticles[intp]					=	0;

					material_fine_data.Strain[intp]							=	0;
					material_fine_data.PlasticStrain[intp]					=	0;
					material_fine_data.SecondInvariantDevStrainrate[intp]	=	0;
					material_fine_data.SecondInvariantDevStress[intp]		=	0;
					material_fine_data.Pressure[intp]						=	0;
					material_fine_data.Temperature[intp]					=	0;

					material_fine_data.PlasticViscosity[intp]				=	0;
					material_fine_data.Plastic[intp]						=	0;


					// also stresses should be interpolated (for viscoelastic simulations)
					material_fine_data.DevStress[0][intp]					=	0;
					material_fine_data.DevStress[1][intp]					=	0;
					material_fine_data.DevStress[2][intp]					=	0;
					material_fine_data.DevStress[3][intp]					=	0;
					material_fine_data.DevStress[4][intp]					=	0;
					material_fine_data.DevStress[5][intp]					=	0;

				}
			}
		}
	}
	ierr = DMDAVecRestoreArray(user->DA_Materials,   user->Materials, &materials_fine_array); CHKERRQ(ierr);

	if (DisplayOutput==1){
		PetscPrintf(PETSC_COMM_WORLD,"Initialized material properties to zero... \n");
	}
	MPI_Barrier(PETSC_COMM_WORLD);


	/* Allocate an array that is used to track the dominant phase at each integration point */
	ierr = DMDAGetCorners(user->DA_Processors,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(user->DA_Processors,&xs_g,&ys_g,&zs_g,&xm_g,&ym_g,&zm_g); CHKERRQ(ierr);

	ierr = PetscMalloc(sizeof(PetscInt *)*(size_t)(xm*ym*zm), &NumberParticlesOfEachPhase);	CHKERRQ(ierr);
	for (i=0; i<(xm*ym*zm); i++){
		ierr = PetscMalloc(sizeof(PetscInt )*(size_t)(user->num_phases), &NumberParticlesOfEachPhase[i]); CHKERRQ(ierr);
	}


	for (i=0; i<(xm*ym*zm); i++){
		for (iphase=0; iphase<user->num_phases; iphase++){
			NumberParticlesOfEachPhase[i][iphase] = 0;
		}
	}
	if (DisplayOutput==1){
		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank=%lld Allocated arrays and starting particle loop... xm,ym,zm=[%lld,%lld,%lld] ngp_vel=%lld num_phases=%lld \n",(LLD)rank, (LLD)xm,(LLD)ym,(LLD)zm,(LLD)ngp_vel, (LLD)(user->num_phases));
		PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);

	}
	MPI_Barrier(PETSC_COMM_WORLD);


	if (DisplayOutput==1){
		PetscPrintf(PETSC_COMM_WORLD,"# Starting particle loop. \n");
		MPI_Barrier(PETSC_COMM_WORLD);
	}

	/* Loop over all particles */
	ierr = DMDAVecGetArray(user->DA_Materials,   user->Materials, &materials_fine_array); CHKERRQ(ierr);
	for (ipart=0; ipart<user->num_particle_local; ipart++){

		/* Compute the element it belongs to  */
		ParticleLocal	=	user->ParticlesLocal[ipart];
		eta             =   ParticleLocal.eta;
		zetha           =   ParticleLocal.zetha;
		phi             =   ParticleLocal.phi;

		if( (vpt_element_type==DAVP_Q1P0) ||   (vpt_element_type==DAVP_Q1Q1) || (vpt_element_type == DAVP_FDSTAG)) {
			ielx = ((PetscInt) ParticleLocal.ix);
			iely = ((PetscInt) ParticleLocal.iy);
			ielz = ((PetscInt) ParticleLocal.iz);
		}
		else if( (vpt_element_type==DAVP_Q2PM1L) || (vpt_element_type==DAVP_Q2PM1G) ) {
			ielx = ((PetscInt) ParticleLocal.ix); LaMEMMod(ielx, 2, &mod); if (mod>0){ielx = ielx-1;}; 	ielx=ielx/2;
			iely = ((PetscInt) ParticleLocal.iy); LaMEMMod(iely, 2, &mod); if (mod>0){iely = iely-1;};	iely=iely/2;
			ielz = ((PetscInt) ParticleLocal.iz); LaMEMMod(ielz, 2, &mod); if (mod>0){ielz = ielz-1;};	ielz=ielz/2;
		}
		else {
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Element type not implemented" );
		}


		/* Check to which integration point we'll put the particle
		 *
		 * This routine can be speed up, by checking in x,y & z directions, which
		 * 	requires 3+3+3=9 computations in Q2P1 case, versus 27 now  (6 versus 8).
		 *
		 * It is however slightly more complicated to compute which is why I leave it as is
		 * 	for now (unless this routine becomes too  time-consuming )
		 *
		 * Further savings can be made for the Q1P0 element, since in this case we only
		 * 	have to check the signs of eta, zetha & phi
		 */
		intp = -1;
		min_dist = 1e4;
		for (i=0; i<ngp_vel; i++){
			distance = 	(IntPoint[0][i]-eta  )*(IntPoint[0][i]-eta  ) +
					(IntPoint[1][i]-zetha)*(IntPoint[1][i]-zetha) +
					(IntPoint[2][i]-phi  )*(IntPoint[2][i]-phi  );
			if (distance<min_dist){
				intp    	 = i;
				min_dist 	 = distance;
			}
		}
		if( intp==-1 ) {  EMERGENCY_EXIT("Did not find a minimum"); }

		// Catch possible errors (it happens almost never, but it crashes the code if, so just in case):
		if ( (ielx<xs) || (ielx>=(xs+xm)) || (iely<ys) || (iely>=(ys+ym)) || (ielz<zs) || (ielz>=(zs+zm)) ){
			PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Particle %lld with number %g on rank %lld with [ix,iy,iz]=[%lld,%lld,%lld] outside domain which goes: [%lld,%lld,%lld] - [%lld,%lld,%lld]. \n",
					(LLD)ipart,ParticleLocal.num,(LLD)rank,(LLD)ielx,(LLD)iely,(LLD)ielz,(LLD)xs,(LLD)ys,(LLD)zs,(LLD)(xm+xs),(LLD)(ym+ys),(LLD)(zm+zs));
			PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
			PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Check NumMaterialPropsElem in LaMEM.h if this occurs frequently! \n");
			PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);


			//
			ParticleLocal.phase = -4;   //deactivate particle
		}

		if (isnan(ParticleLocal.phase)){
			PetscPrintf(PETSC_COMM_SELF,"Detected a particle with phase==NaN! particle # %lld \n",(LLD)ipart);
		}


		if ( (ielx>=xs) && (ielx<(xs+xm)) && (iely>=ys) && (iely<(ys+ym)) && (ielz>=zs) && (ielz<(zs+zm)) ){

			/* Initialize array */
			LaMEMSetMaterialDataMemoryFromArray( &material_fine_data, ielx-xs,iely-ys,ielz-zs, C->ngp_vel, materials_fine_array );

			/* Add properties to the integration point */
			phase 	= 	((PetscInt) ParticleLocal.phase);

			if ( phase>=0  && phase<user->num_phases){


				material_fine_data.NumParticles[intp] 					= 	material_fine_data.NumParticles[intp] 		+ 1.0;

				material_fine_data.Strain[intp]							=	material_fine_data.Strain[intp]								+ ParticleLocal.Strain;
				material_fine_data.PlasticStrain[intp]					=	material_fine_data.PlasticStrain[intp]						+ ParticleLocal.PlasticStrain;
				material_fine_data.SecondInvariantDevStrainrate[intp]	=	material_fine_data.SecondInvariantDevStrainrate[intp]	+ ParticleLocal.E2nd;
				material_fine_data.SecondInvariantDevStress[intp]		=	material_fine_data.SecondInvariantDevStress[intp]			+ ParticleLocal.T2nd;
				material_fine_data.Pressure[intp]						=	material_fine_data.Pressure[intp]							+ ParticleLocal.P;
				material_fine_data.Temperature[intp]					=	material_fine_data.Temperature[intp]						+ ParticleLocal.T;
				material_fine_data.Plastic[intp]						=	material_fine_data.Plastic[intp]							+ ParticleLocal.Plastic;


				/* 3D deviatoric stress tensor */
				material_fine_data.DevStress[0][intp]	 				=	material_fine_data.DevStress[0][intp]						+	ParticleLocal.Txx;
				material_fine_data.DevStress[1][intp]	 				=	material_fine_data.DevStress[1][intp]						+	ParticleLocal.Tyy;
				material_fine_data.DevStress[2][intp]	 				=	material_fine_data.DevStress[2][intp]						+	ParticleLocal.Tzz;
				material_fine_data.DevStress[3][intp]	 				=	material_fine_data.DevStress[3][intp]						+	ParticleLocal.Txy;
				material_fine_data.DevStress[4][intp]	 				=	material_fine_data.DevStress[4][intp]						+	ParticleLocal.Txz;
				material_fine_data.DevStress[5][intp]	 				=	material_fine_data.DevStress[5][intp]						+	ParticleLocal.Tyz;


				/* Keep track of which phases are being set at this element */
				el_number = 	(ielz-zs_g)*( (xm_g)*(ym_g) ) + (iely-ys_g)*(xm_g) + (ielx-xs_g);

				NumberParticlesOfEachPhase[el_number][phase] = NumberParticlesOfEachPhase[el_number][phase] + 1;

			}
		}

		user->ParticlesLocal[ipart] =	ParticleLocal;	// put back particle, in case it gets deactivated along the way


	}
	MPI_Barrier(PETSC_COMM_WORLD);
	if (DisplayOutput==1){
		PetscPrintf(PETSC_COMM_WORLD,"# Finished particle loop... \n");
	}

	/* Set the most abundant phase to every integration point in a cell*/
	ierr = DMDAGetCorners(user->DA_Processors,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
	for (ielz=zs; ielz<zs+zm; ielz++){
		for (iely=ys; iely<ys+ym; iely++){
			for (ielx=xs; ielx<xs+xm; ielx++){

				LaMEMSetMaterialDataMemoryFromArray( &material_fine_data, ielx-xs, iely-ys, ielz-zs, C->ngp_vel, materials_fine_array );
				TotalNumberOfParticlesInElement = 0;
				for (intp=0; intp<ngp_vel; intp++){
					TotalNumberOfParticlesInElement = TotalNumberOfParticlesInElement  + (PetscInt) material_fine_data.NumParticles[intp];

					el_number = 	(ielz-zs_g)*( (xm_g)*(ym_g) ) + (iely-ys_g)*(xm_g) + (ielx-xs_g);

					MostAbundantPhase = 0;
					MaxNumberTracers  = 0;
					for (iphase=0; iphase<user->num_phases; iphase++){
						if (NumberParticlesOfEachPhase[el_number][iphase]>MaxNumberTracers){
							MaxNumberTracers  = NumberParticlesOfEachPhase[el_number][iphase];
							MostAbundantPhase = iphase;
						}
					}
				}

				for (intp=0; intp<ngp_vel; intp++){
					material_fine_data.Phases[intp] = ((double) MostAbundantPhase);  //set phase of element to the most abundant phase in the element
				}

			}
		}
	}




	/* Most material properties are computed by using the dominating phase at each integration point (as computed above).
	 * An exception is made, however, for purely history variables such as strain, plastic strain or stresses.
	 *
	 * These are computed here by arithmetic averaging over the part surrounding an integration point.
	 *
	 * Loop over local elements and finish the averaging in every intp.
	 * 		Compute mean (element-wise) properties
	 * 		If an integration point has no particles attributed to it, set the mean value
	 */
	if (DisplayOutput==1){
		PetscPrintf(PETSC_COMM_WORLD,"Starting averaging of material properties at integration points... \n");
	}
	ierr = DMDAGetCorners(user->DA_Processors,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
	for (ielz=zs; ielz<zs+zm; ielz++){
		for (iely=ys; iely<ys+ym; iely++){
			for (ielx=xs; ielx<xs+xm; ielx++){
				LaMEMSetMaterialDataMemoryFromArray( &material_fine_data, ielx-xs,iely-ys,ielz-zs, C->ngp_vel, materials_fine_array );

				// Make a loop over each integration point
				MeanStrain 		= 	MeanPlasticStrain 	= 	MeanE2nd = MeanT2nd = 0.0;
				MeanPressure 	= 	MeanTemperature 	= 	MeanPlastic = 0.0;
				MeanTxx 		=	MeanTyy				=	MeanTzz = 0.0;
				MeanTxz 		=	MeanTyz				=	MeanTxy = 0.0;
				MeanMuViscous 	=	MeanValue 			=			0.0;

				for (intp=0; intp<ngp_vel; intp++){

					if (material_fine_data.NumParticles[intp]>0){

						material_fine_data.Strain[intp]							=	material_fine_data.Strain[intp]/material_fine_data.NumParticles[intp];
						material_fine_data.PlasticStrain[intp]					=	material_fine_data.PlasticStrain[intp]/material_fine_data.NumParticles[intp];
						material_fine_data.SecondInvariantDevStrainrate[intp]	=	material_fine_data.SecondInvariantDevStrainrate[intp]/material_fine_data.NumParticles[intp];
						material_fine_data.SecondInvariantDevStress[intp]		=	material_fine_data.SecondInvariantDevStress[intp]/material_fine_data.NumParticles[intp];
						material_fine_data.Pressure[intp]						=	material_fine_data.Pressure[intp]/material_fine_data.NumParticles[intp];
						material_fine_data.Temperature[intp]					=	material_fine_data.Temperature[intp]/material_fine_data.NumParticles[intp];

						material_fine_data.Plastic[intp]						=	material_fine_data.Plastic[intp]/material_fine_data.NumParticles[intp];

						MeanStrain 				= 	MeanStrain 			+ 	material_fine_data.Strain[intp]	;
						MeanPlasticStrain 		= 	MeanPlasticStrain 	+ 	material_fine_data.PlasticStrain[intp]	;
						MeanE2nd 				=	MeanE2nd			+	material_fine_data.SecondInvariantDevStrainrate[intp];
						MeanT2nd 				=	MeanT2nd			+	material_fine_data.SecondInvariantDevStress[intp];
						MeanPressure 			=	MeanPressure		+	material_fine_data.Pressure[intp];
						MeanTemperature 		=	MeanTemperature		+	material_fine_data.Temperature[intp];
						MeanPlastic 			=	MeanPlastic 		+	material_fine_data.Plastic[intp];

						MeanTxx 				=	MeanTxx 			+	material_fine_data.DevStress[0][intp];
						MeanTyy 				=	MeanTyy 			+	material_fine_data.DevStress[1][intp];
						MeanTzz 				=	MeanTzz 			+	material_fine_data.DevStress[2][intp];
						MeanTxy 				=	MeanTxy 			+	material_fine_data.DevStress[3][intp];
						MeanTxz 				=	MeanTxz 			+	material_fine_data.DevStress[4][intp];
						MeanTyz 				=	MeanTyz 			+	material_fine_data.DevStress[5][intp];

						/* 3D deviatoric stress tensor */
						material_fine_data.DevStress[0][intp]	 				=	material_fine_data.DevStress[0][intp]/material_fine_data.NumParticles[intp];
						material_fine_data.DevStress[1][intp]	 				=	material_fine_data.DevStress[1][intp]/material_fine_data.NumParticles[intp];
						material_fine_data.DevStress[2][intp]	 				=	material_fine_data.DevStress[2][intp]/material_fine_data.NumParticles[intp];
						material_fine_data.DevStress[3][intp]	 				=	material_fine_data.DevStress[3][intp]/material_fine_data.NumParticles[intp];
						material_fine_data.DevStress[4][intp]	 				=	material_fine_data.DevStress[4][intp]/material_fine_data.NumParticles[intp];
						material_fine_data.DevStress[5][intp]	 				=	material_fine_data.DevStress[5][intp]/material_fine_data.NumParticles[intp];



						MeanValue     		= MeanValue + 1;
					}


				}

				if (MeanValue>0){
					MeanStrain 			=	MeanStrain/MeanValue;
					MeanPlasticStrain	=	MeanPlasticStrain/MeanValue;
					MeanPlastic			=	MeanPlastic/MeanValue;
					MeanMuViscous 		=	MeanMuViscous/MeanValue;

					MeanE2nd			=	MeanE2nd/MeanValue;
					MeanT2nd			=	MeanT2nd/MeanValue;
					MeanPressure		=	MeanPressure/MeanValue;
					MeanTemperature		=	MeanTemperature/MeanValue;

					MeanTxx 			=	MeanTxx/MeanValue;
					MeanTyy 			=	MeanTyy/MeanValue;
					MeanTzz 			=	MeanTzz/MeanValue;
					MeanTxy 			=	MeanTxy/MeanValue;
					MeanTxz 			=	MeanTxz/MeanValue;
					MeanTyz 			=	MeanTyz/MeanValue;

				}
				else	{
					/* If the elements is really empty, inject a phase */
					PetscPrintf(PETSC_COMM_WORLD,"***** Possible error detected: mean viscosity =0 in [ielx,iely,ielz]=[%lld, %lld, %lld]: using values of phase %lld ***** \n",(LLD)ielx, (LLD)iely, (LLD)ielz, (LLD)(user->ParticleInjectionPhase) );
					MeanStrain 			=	0.0;
					MeanPlasticStrain 	=	0.0;
					MeanPlastic			=	0.0;
					MeanE2nd 			=	0.0;
					MeanT2nd 			=	0.0;
					MeanPressure		=	0.0;
					MeanTemperature 	=	0.0;

					MeanTxx 			=	0.0;
					MeanTyy 			=	0.0;
					MeanTzz 			=	0.0;
					MeanTxy 			=	0.0;
					MeanTxz 			=	0.0;
					MeanTyz 			=	0.0;

				}


				// Check for empty integration points
				for (intp=0; intp<ngp_vel; intp++){
					if (material_fine_data.NumParticles[intp]==0){

						material_fine_data.Strain[intp]								=	MeanStrain;
						material_fine_data.PlasticStrain[intp]						=	MeanPlasticStrain;
						material_fine_data.Plastic[intp]							=	MeanPlastic;
						material_fine_data.SecondInvariantDevStrainrate[intp]		=	MeanE2nd;
						material_fine_data.SecondInvariantDevStress[intp]			=	MeanT2nd;
						material_fine_data.Pressure[intp]							=	MeanPressure;
						material_fine_data.Temperature[intp]						=	MeanTemperature;

						material_fine_data.DevStress[0][intp]	 					=	MeanTxx;
						material_fine_data.DevStress[1][intp]	 					=	MeanTyy;
						material_fine_data.DevStress[2][intp]	 					=	MeanTzz;
						material_fine_data.DevStress[3][intp]	 					=	MeanTxy;
						material_fine_data.DevStress[4][intp]	 					=	MeanTxz;
						material_fine_data.DevStress[5][intp]	 					=	MeanTyz;

					}

					// Uncomment this, if you want only mean props per element
					//material_fine_data.Strain[intp]			= 	MeanStrain;
					//material_fine_data.PlasticStrain[intp] 	=	MeanPlasticStrain;

				}

			}
		}
	}


	// Make one more loop and determine whether an integration point is plastic or not. If it is plastic, compute the plastic viscosity
	for (ielz=zs; ielz<zs+zm; ielz++){
		for (iely=ys; iely<ys+ym; iely++){
			for (ielx=xs; ielx<xs+xm; ielx++){
				LaMEMSetMaterialDataMemoryFromArray( &material_fine_data, ielx-xs,iely-ys,ielz-zs, C->ngp_vel, materials_fine_array );
				for (intp=0; intp<ngp_vel; intp++){
					if (material_fine_data.Plastic[intp]>0.5){
						material_fine_data.Plastic[intp] 			= 	1.0;
						material_fine_data.PlasticViscosity[intp] 	= 	material_fine_data.SecondInvariantDevStress[intp]/(2.0*material_fine_data.SecondInvariantDevStrainrate[intp]);

					}
					else {
						material_fine_data.Plastic[intp] 			=	0.0;
						material_fine_data.PlasticViscosity[intp] 	= 	0.0;
					}

				}
			}
		}
	}


	ierr = DMDAVecRestoreArray(user->DA_Materials, user->Materials, &materials_fine_array); CHKERRQ(ierr);
	if (DisplayOutput==1){
		PetscPrintf(PETSC_COMM_WORLD,"Finished averaging of material properties at integration points... \n");
	}

	/* Cleaning up */
	ierr = DMDAGetCorners(user->DA_Processors,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
	for (i=0; i<(xm*ym*zm); i++){
		ierr =  PetscFree(NumberParticlesOfEachPhase[i]);  CHKERRQ(ierr);
	}
	ierr= PetscFree(NumberParticlesOfEachPhase); CHKERRQ(ierr);

	PetscFunctionReturn(0);

}
/* ==================================================================================================== */



/* ====================================================================================================
 * Compute material properties and history variables (temperature, stress, strain etc.) at all points where
 * the FDSTAG stencil needs them:
 * 	- center points
 *	- Sxz points
 *	- Sxy points
 *	- Syz points
 *
 *	We employ a distance-based weighting of properties. At every point we compute the fraction of each phase that is present
 *	such that we -at a later stage- we can determine for example the effective viscosity or effective density at this point in
 *	a straightforward manner (even if, for example, phase 0 has a linear viscosity, and say phase 2 a powerlaw T- and P-dependent viscosity.
 *	Yet, at the same time, the code also needs to know some history-dependent variables such as strain, stress-tensor, temperature, pressure
 *	at the finite difference points. These variables are also computed in this routine and are stored as global vectors.
 *
 * */
#undef __FUNCT__
#define __FUNCT__ "ComputeHistoryPropertiesFromParticles_FDSTAG"
PetscErrorCode ComputeHistoryPropertiesFromParticles_FDSTAG(LaMEMVelPressureDA C, UserContext *user)
{
	PetscMPIInt   			rank;
	PetscErrorCode			ierr;
	PetscLogDouble			t0, t1;
	DAVPElementType 		vpt_element_type;

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	PetscTime(&t0);
	PetscPrintf( PETSC_COMM_WORLD, "#  Computing material properties from particles for FDSTAG ... ", (LLD) user->num_phases );


	vpt_element_type = C->type;
	if (vpt_element_type != DAVP_FDSTAG){
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "This routine only works for FDSTAG ! " );
	}



	/* Evaluate properties at center points ----------------------------------*/
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(PETSC_NULL, 							user, "PhaseProportions","CenterNodes"); CHKERRQ(ierr);			// compute phase proportions
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(user->FDSTAG.Center_Temperature, 		user, "Temperature"		,"CenterNodes"); CHKERRQ(ierr);			// average temperature
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(user->FDSTAG.Center_Pressure, 			user, "Pressure"		,"CenterNodes"); CHKERRQ(ierr);			// average pressure
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(user->FDSTAG.Center_Strain, 			user, "Strain"			,"CenterNodes"); CHKERRQ(ierr);			// average strain
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(user->FDSTAG.Center_PlasticStrain,		user, "PlasticStrain"	,"CenterNodes"); CHKERRQ(ierr);			// average plastic strain
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(user->FDSTAG.Center_T2nd, 				user, "T2nd"			,"CenterNodes"); CHKERRQ(ierr);			// average 2nd invarian of deviatoric stress tensor
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(user->FDSTAG.Center_E2nd, 				user, "E2nd"			,"CenterNodes"); CHKERRQ(ierr);			// average 2nd invariant of strain rate tensor
	/* End of evaluating properties at center points --------------------------*/

	/* Evaluate properties at nodal points ------------------------------------*/
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(PETSC_NULL, 							user, "PhaseProportions","CornerNodes"); CHKERRQ(ierr);			// compute phase proportions
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(user->FDSTAG.Corner_Temperature, 		user, "Temperature"		,"CornerNodes"); CHKERRQ(ierr);			// average temperature
	/* End of evaluating properties at nodal points ---------------------------*/

	/* Evaluate properties at XY points ---------------------------------------*/
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(PETSC_NULL, 							user, "PhaseProportions","XYNodes"); CHKERRQ(ierr);				// compute phase proportions
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(user->FDSTAG.XYPoints_Temperature, 	user, "Temperature"		,"XYNodes"); CHKERRQ(ierr);				// average temperature
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(user->FDSTAG.XYPoints_Pressure, 		user, "Pressure"		,"XYNodes"); CHKERRQ(ierr);				// average pressure
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(user->FDSTAG.XYPoints_Strain, 			user, "Strain"			,"XYNodes"); CHKERRQ(ierr);				// average strain
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(user->FDSTAG.XYPoints_PlasticStrain,	user, "PlasticStrain"	,"XYNodes"); CHKERRQ(ierr);				// average plastic strain
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(user->FDSTAG.XYPoints_T2nd, 			user, "T2nd"			,"XYNodes"); CHKERRQ(ierr);				// average 2nd invarian of deviatoric stress tensor
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(user->FDSTAG.XYPoints_E2nd, 			user, "E2nd"			,"XYNodes"); CHKERRQ(ierr);				// average 2nd invariant of strain rate tensor
	/* End of evaluating properties at XY points ------------------------------*/

	/* Evaluate properties at XZ points ---------------------------------------*/
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(PETSC_NULL, 							user, "PhaseProportions","XZNodes"); CHKERRQ(ierr);				// compute phase proportions
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(user->FDSTAG.XZPoints_Temperature, 	user, "Temperature"		,"XZNodes"); CHKERRQ(ierr);				// average temperature
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(user->FDSTAG.XZPoints_Pressure, 		user, "Pressure"		,"XZNodes"); CHKERRQ(ierr);				// average pressure
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(user->FDSTAG.XZPoints_Strain, 			user, "Strain"			,"XZNodes"); CHKERRQ(ierr);				// average strain
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(user->FDSTAG.XZPoints_PlasticStrain,	user, "PlasticStrain"	,"XZNodes"); CHKERRQ(ierr);				// average plastic strain
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(user->FDSTAG.XZPoints_T2nd, 			user, "T2nd"			,"XZNodes"); CHKERRQ(ierr);				// average 2nd invarian of deviatoric stress tensor
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(user->FDSTAG.XZPoints_E2nd, 			user, "E2nd"			,"XZNodes"); CHKERRQ(ierr);				// average 2nd invariant of strain rate tensor
	/* End of evaluating properties at XZ points ------------------------------*/

	/* Evaluate properties at YZ points ---------------------------------------*/
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(PETSC_NULL, 							user, "PhaseProportions","YZNodes"); CHKERRQ(ierr);				// compute phase proportions
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(user->FDSTAG.YZPoints_Temperature, 	user, "Temperature"		,"YZNodes"); CHKERRQ(ierr);				// average temperature
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(user->FDSTAG.YZPoints_Pressure, 		user, "Pressure"		,"YZNodes"); CHKERRQ(ierr);				// average pressure
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(user->FDSTAG.YZPoints_Strain, 			user, "Strain"			,"YZNodes"); CHKERRQ(ierr);				// average strain
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(user->FDSTAG.YZPoints_PlasticStrain,	user, "PlasticStrain"	,"YZNodes"); CHKERRQ(ierr);				// average plastic strain
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(user->FDSTAG.YZPoints_T2nd, 			user, "T2nd"			,"YZNodes"); CHKERRQ(ierr);				// average 2nd invarian of deviatoric stress tensor
	ierr =  AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(user->FDSTAG.YZPoints_E2nd, 			user, "E2nd"			,"YZNodes"); CHKERRQ(ierr);				// average 2nd invariant of strain rate tensor
	/* End of evaluating properties at YZ points ------------------------------*/

	PetscTime(&t1);
	PetscPrintf( PETSC_COMM_WORLD, "  [%g s] \n",t1-t0 );



	/* For debugging: write some stuff to disk */
	{
		PetscBool  OutputVTK_PhaseProportions;
		PetscInt iphase;

		OutputVTK_PhaseProportions = PETSC_FALSE;
		ierr = PetscOptionsHasName(PETSC_NULL,"-OutputVTK_PhaseProportions",&OutputVTK_PhaseProportions); CHKERRQ(ierr);


		if (OutputVTK_PhaseProportions){
			PetscViewer vv;
			char 						*PhaseName = NULL;



			ierr = UpdateMaterialProperties_FDSTAG(user); CHKERRQ(ierr); 	// UPDATE MATERIAL PROPERTIES @ THE NODES - FOR DEBUGGING ONLY [SHOULD BE IN MAIN CODE]


			ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "MaterialProperties_Center_FDSTAG.vtk", &vv);CHKERRQ(ierr);
			ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
			ierr = DMView(user->FDSTAG.DA_CENTER, vv);CHKERRQ(ierr);
			for (iphase=0;  iphase<user->num_phases;  iphase++){
				asprintf(&PhaseName, "Phase_%1lld",(LLD)iphase);

				ierr = PetscObjectSetName( (PetscObject)user->FDSTAG.Center_PhaseProportions[iphase], PhaseName );CHKERRQ(ierr);
				ierr = VecView(user->FDSTAG.Center_PhaseProportions[iphase], vv);CHKERRQ(ierr);
			}
			ierr = PetscObjectSetName( (PetscObject)user->FDSTAG.Center_Temperature, "Temperature" );CHKERRQ(ierr);
			ierr = VecView(user->FDSTAG.Center_Temperature, vv);CHKERRQ(ierr);
			ierr = PetscObjectSetName( (PetscObject)user->FDSTAG.Center_EffectiveViscosity, "Effective_Viscosity" );CHKERRQ(ierr);
			ierr = VecView(user->FDSTAG.Center_EffectiveViscosity, vv);CHKERRQ(ierr);
			ierr = PetscViewerDestroy(&vv);CHKERRQ(ierr);
			PetscPrintf(PETSC_COMM_WORLD,"FDSTAG DEBUG: wrote center material properties to  MaterialProperties_Center_FDSTAG.vtk \n");

			ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "MaterialProperties_Corner_FDSTAG.vtk", &vv);CHKERRQ(ierr);
			ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
			ierr = DMView(user->FDSTAG.DA_CORNER, vv);CHKERRQ(ierr);
			for (iphase=0;  iphase<user->num_phases;  iphase++){
				asprintf(&PhaseName, "Phase_%1lld",(LLD)iphase);

				ierr = PetscObjectSetName( (PetscObject)user->FDSTAG.Corner_PhaseProportions[iphase], PhaseName );CHKERRQ(ierr);
				ierr = VecView(user->FDSTAG.Corner_PhaseProportions[iphase], vv);CHKERRQ(ierr);
			}
			ierr = PetscObjectSetName( (PetscObject)user->FDSTAG.Corner_Temperature, "Temperature" );CHKERRQ(ierr);
			ierr = VecView(user->FDSTAG.Corner_Temperature, vv);CHKERRQ(ierr);
			ierr = PetscViewerDestroy(&vv);CHKERRQ(ierr);
			PetscPrintf(PETSC_COMM_WORLD,"FDSTAG DEBUG: wrote center material properties to  MaterialProperties_Corner_FDSTAG.vtk \n");

			ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "MaterialProperties_XY_FDSTAG.vtk", &vv);CHKERRQ(ierr);
			ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
			ierr = DMView(user->FDSTAG.DA_XY_POINTS, vv);CHKERRQ(ierr);
			for (iphase=0;  iphase<user->num_phases;  iphase++){
				asprintf(&PhaseName, "Phase_%1lld",(LLD)iphase);
				ierr = PetscObjectSetName( (PetscObject)user->FDSTAG.XYPoints_PhaseProportions[iphase], PhaseName );CHKERRQ(ierr);
				ierr = VecView(user->FDSTAG.XYPoints_PhaseProportions[iphase], vv);CHKERRQ(ierr);
			}
			ierr = PetscObjectSetName( (PetscObject)user->FDSTAG.XYPoints_EffectiveViscosity, "Effective_Viscosity" );CHKERRQ(ierr);
			ierr = VecView(user->FDSTAG.XYPoints_EffectiveViscosity, vv);CHKERRQ(ierr);
			ierr = PetscViewerDestroy(&vv);CHKERRQ(ierr);
			PetscPrintf(PETSC_COMM_WORLD,"FDSTAG DEBUG: wrote XY points to  MaterialProperties_XY_FDSTAG.vtk \n");

			ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "MaterialProperties_XZ_FDSTAG.vtk", &vv);CHKERRQ(ierr);
			ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
			ierr = DMView(user->FDSTAG.DA_XZ_POINTS, vv);CHKERRQ(ierr);
			for (iphase=0;  iphase<user->num_phases;  iphase++){
				asprintf(&PhaseName, "Phase_%1lld",(LLD)iphase);
				ierr = PetscObjectSetName( (PetscObject)user->FDSTAG.XZPoints_PhaseProportions[iphase], PhaseName );CHKERRQ(ierr);
				ierr = VecView(user->FDSTAG.XZPoints_PhaseProportions[iphase], vv);CHKERRQ(ierr);
			}
			ierr = PetscObjectSetName( (PetscObject)user->FDSTAG.XZPoints_EffectiveViscosity, "Effective_Viscosity" );CHKERRQ(ierr);
			ierr = VecView(user->FDSTAG.XZPoints_EffectiveViscosity, vv);CHKERRQ(ierr);
			ierr = PetscViewerDestroy(&vv);CHKERRQ(ierr);
			PetscPrintf(PETSC_COMM_WORLD,"FDSTAG DEBUG: wrote XZ points to  MaterialProperties_XZ_FDSTAG.vtk \n");

			ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "MaterialProperties_YZ_FDSTAG.vtk", &vv);CHKERRQ(ierr);
			ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
			ierr = DMView(user->FDSTAG.DA_YZ_POINTS, vv);CHKERRQ(ierr);
			for (iphase=0;  iphase<user->num_phases;  iphase++){
				asprintf(&PhaseName, "Phase_%1lld",(LLD)iphase);
				ierr = PetscObjectSetName( (PetscObject)user->FDSTAG.YZPoints_PhaseProportions[iphase], PhaseName );CHKERRQ(ierr);
				ierr = VecView(user->FDSTAG.YZPoints_PhaseProportions[iphase], vv);CHKERRQ(ierr);
			}
			ierr = PetscObjectSetName( (PetscObject)user->FDSTAG.YZPoints_EffectiveViscosity, "Effective_Viscosity" );CHKERRQ(ierr);
			ierr = VecView(user->FDSTAG.YZPoints_EffectiveViscosity, vv);CHKERRQ(ierr);

			ierr = PetscViewerDestroy(&vv);CHKERRQ(ierr);
			PetscPrintf(PETSC_COMM_WORLD,"FDSTAG DEBUG: wrote YZ points to  MaterialProperties_YZ_FDSTAG.vtk \n");
		}
	}

	PetscFunctionReturn(0);
}
/* ==================================================================================================== */


/* ====================================================================================================
 * Computes (history) properties from particles at center nodes using a distance-based averaging scheme.
 * A string is given as input to indicate which properties are being averaged.
 *
 */
#undef __FUNCT__
#define __FUNCT__ "AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG"
PetscErrorCode AverageHistoryPropertiesFromParticles_AroundNodes_FDSTAG(Vec Global_MaterialProperties, UserContext *user, const char PropertyToBeAveraged[], const char NodeAtWhichAveragingIsPerformed[])
{
	PetscMPIInt   			rank;
	PetscErrorCode			ierr;
	PetscBool 				ComputingPhaseProportions=PETSC_FALSE, flag=PETSC_FALSE;
	PetscInt				i, ipart, phase, PropertyNumberToBeAveraged, AveragingNode;
	PetscInt 				ielx, iely, ielz;
	Particles				ParticleLocal;
	DM 						da;
	Vec 					Vec_PhaseProportions_local[max_num_phases], NumberParticles_local, ShapeFunctions_local;
	Vec 					NumberParticles_global, ShapeFunctions_global, Local_MaterialProperties_Vec;
	PetscScalar 			***Numberparticles, ***PhaseProportions_local[max_num_phases], ***ShapeFunctions, ***MaterialProperties_local;
	PetscScalar 			eta, zetha, phi, x, y, z, PropertyValue;
	PetscScalar				Point[3], ShapeFunction[MAX_nnel], **dhds;
	PetscInt				xs,ys,zs,xm,ym,zm;

	LaMEMCreate2dArray( 3,8, &dhds, PETSC_NULL );


	/* Around which nodes do we average properties? */
	AveragingNode=-1; da = user->FDSTAG.DA_CENTER;
	PetscStrcmp(NodeAtWhichAveragingIsPerformed,"CenterNodes",&flag); 	if (flag){	AveragingNode=0;	da = user->FDSTAG.DA_CENTER; }
	PetscStrcmp(NodeAtWhichAveragingIsPerformed,"CornerNodes",&flag); 	if (flag){	AveragingNode=1;	da = user->FDSTAG.DA_CORNER; }
	PetscStrcmp(NodeAtWhichAveragingIsPerformed,"XYNodes",&flag); 		if (flag){	AveragingNode=2;	da = user->FDSTAG.DA_XY_POINTS; }
	PetscStrcmp(NodeAtWhichAveragingIsPerformed,"XZNodes",&flag); 		if (flag){	AveragingNode=3;	da = user->FDSTAG.DA_XZ_POINTS; }
	PetscStrcmp(NodeAtWhichAveragingIsPerformed,"YZNodes",&flag); 		if (flag){	AveragingNode=4;	da = user->FDSTAG.DA_YZ_POINTS; }


	if (AveragingNode==-1){
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_SUP, "We did not yet implement averaging for %s",NodeAtWhichAveragingIsPerformed);
	}


	/* Are we dealing with PhaseProportions or are we dealing with something else?
	 * PhaseProportions are dealt with seperately in this routine
	 */
	PetscStrcmp(PropertyToBeAveraged,"PhaseProportions",&ComputingPhaseProportions);
	/* Which property are we averaging? */
	PropertyNumberToBeAveraged=-1;
	if (ComputingPhaseProportions){
		PropertyNumberToBeAveraged=0;
	}
	else{
		PetscStrcmp(PropertyToBeAveraged,"Temperature",&flag); 		if (flag){	PropertyNumberToBeAveraged=1;	}
		PetscStrcmp(PropertyToBeAveraged,"Pressure",&flag); 		if (flag){	PropertyNumberToBeAveraged=2;	}
		PetscStrcmp(PropertyToBeAveraged,"Strain",&flag); 			if (flag){	PropertyNumberToBeAveraged=3;	}
		PetscStrcmp(PropertyToBeAveraged,"PlasticStrain",&flag); 	if (flag){	PropertyNumberToBeAveraged=4;	}
		PetscStrcmp(PropertyToBeAveraged,"T2nd",&flag); 			if (flag){	PropertyNumberToBeAveraged=5;	}
		PetscStrcmp(PropertyToBeAveraged,"E2nd",&flag); 			if (flag){	PropertyNumberToBeAveraged=6;	}
	}


	/* Retrieve vectors for storing phase proportions - phase proportions are dealt with in a special manner*/
	if (ComputingPhaseProportions){
		for (i=0; i<user->num_phases; i++){
			ierr = DMGetLocalVector(da,	&Vec_PhaseProportions_local[i]);	CHKERRQ(ierr);		// local vector that will store local phase proportions.
			ierr = VecZeroEntries(Vec_PhaseProportions_local[i]);		CHKERRQ(ierr);
			ierr = DMDAVecGetArray(da,   	Vec_PhaseProportions_local[i], 	&PhaseProportions_local[i]	 ); CHKERRQ(ierr); // Phase proportions
		}
		PropertyNumberToBeAveraged = 0;
	}
	else{
		ierr = DMGetLocalVector(da,	&Local_MaterialProperties_Vec);	CHKERRQ(ierr);		// local vector that will store local material properties
		ierr = VecZeroEntries(Local_MaterialProperties_Vec);		CHKERRQ(ierr);
		ierr = DMDAVecGetArray(da,   	Local_MaterialProperties_Vec, 	&MaterialProperties_local	 ); CHKERRQ(ierr); 	// will contain averaged material properties
	}


	/* These vectors & arrays have to be allocated in any case */
	ierr = DMGetLocalVector(da,	&NumberParticles_local);	CHKERRQ(ierr);
	ierr = VecZeroEntries(NumberParticles_local);			CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da,   	NumberParticles_local, 	&Numberparticles	 ); CHKERRQ(ierr); // Viscosity @ nodes
	ierr = DMGetLocalVector(da,	&ShapeFunctions_local);	CHKERRQ(ierr);
	ierr = VecZeroEntries(ShapeFunctions_local);			CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da,   	ShapeFunctions_local, 	&ShapeFunctions	 ); CHKERRQ(ierr); // Viscosity @ nodes


	ierr = DMCreateGlobalVector(da,	&NumberParticles_global); 	CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da,	&ShapeFunctions_global); 	CHKERRQ(ierr);


	DMDAGetGhostCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);



	/* PART 1: Make a loop over all local particles, and compute viscosity and density @ nodes */
	for (ipart=0; ipart<user->num_particle_local; ipart++){

		/* Compute the element the particle belongs to */
		ParticleLocal	=	user->ParticlesLocal[ipart];
		phase			=	((PetscInt) ParticleLocal.phase);  // phase of current particle



		if (phase>=0){ // only active particles count
			eta             =   ParticleLocal.eta;
			zetha           =   ParticleLocal.zetha;
			phi             =   ParticleLocal.phi;

			ielx 			= 	((PetscInt) ParticleLocal.ix);
			iely 			= 	((PetscInt) ParticleLocal.iy);
			ielz 			= 	((PetscInt) ParticleLocal.iz);

			/* Which local material properties are we averaging? */
			PropertyValue =	0;
			if 		(PropertyNumberToBeAveraged==0){	// PhaseProportions - dealt with separately
			}
			else if (PropertyNumberToBeAveraged==1){ 	// Temperature
				PropertyValue = ParticleLocal.T;
			}
			else if (PropertyNumberToBeAveraged==2){ 	// Pressure
				PropertyValue = ParticleLocal.P;
			}
			else if (PropertyNumberToBeAveraged==3){ 	// Strain
				PropertyValue = ParticleLocal.Strain;
			}
			else if (PropertyNumberToBeAveraged==4){ 	// PlasticStrain
				PropertyValue = ParticleLocal.PlasticStrain;
			}
			else if (PropertyNumberToBeAveraged==5){ 	// T2nd
				PropertyValue = ParticleLocal.T2nd;
			}
			else if (PropertyNumberToBeAveraged==6){ 	// E2nd
				PropertyValue = ParticleLocal.E2nd;
			}
			else{
				SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Unknown property to be averaged from Particles" );
			}

			/* Add phase to the correct vector */
			if 		(AveragingNode==0){				/* Average center nodes */

				// shift shape function
				x 				=	1.0-PetscAbs(eta);
				y 				=	1.0-PetscAbs(zetha);
				z 				=	1.0-PetscAbs(phi);

				ShapeFunction[0] 			=	x*y*z;

				if (phase>=user->num_phases){
					SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Particle has phase %i, but the maximum number of phases in the code is %i. Check whether num_phases is specified correctly in the input file.",phase,user->num_phases);
				}

				if (ComputingPhaseProportions){

					PhaseProportions_local[phase][ielz  ][iely  ][ielx  ] += ShapeFunction[0];		// Add linear shape function
				}
				else{
					MaterialProperties_local[ielz  ][iely  ][ielx  ]  +=	PropertyValue*ShapeFunction[0];
				}
				ShapeFunctions[ielz  ][iely  ][ielx  ] += ShapeFunction[0];				//	collect shape function contributions

				/* Keep track of the number of particles added */
				Numberparticles[ielz  ][iely  ][ielx  ] += 1;

			}										/* End of averaging center nodes */
			else if (AveragingNode==1){				/* Average corner nodes */

				/* Evaluate shape function at (eta,zetha,phi) */
				Point[0] 		= 	eta;
				Point[1] 		= 	zetha;
				Point[2] 		=	phi;
				ierr = ComputeShapeFunctionVelocity(ShapeFunction, dhds, Point); CHKERRQ(ierr);


				/* Add contributions to the 8 surrounding nodes */
				Numberparticles[ielz  ][iely  ][ielx  ] += 1;
				Numberparticles[ielz  ][iely  ][ielx+1] += 1;
				Numberparticles[ielz  ][iely+1][ielx+1] += 1;
				Numberparticles[ielz  ][iely+1][ielx  ] += 1;
				Numberparticles[ielz+1][iely  ][ielx  ] += 1;
				Numberparticles[ielz+1][iely  ][ielx+1] += 1;
				Numberparticles[ielz+1][iely+1][ielx+1] += 1;
				Numberparticles[ielz+1][iely+1][ielx  ] += 1;

				ShapeFunctions[ielz  ][iely  ][ielx  ] += ShapeFunction[0];
				ShapeFunctions[ielz  ][iely  ][ielx+1] += ShapeFunction[1];
				ShapeFunctions[ielz  ][iely+1][ielx+1] += ShapeFunction[2];
				ShapeFunctions[ielz  ][iely+1][ielx  ] += ShapeFunction[3];
				ShapeFunctions[ielz+1][iely  ][ielx  ] += ShapeFunction[4];
				ShapeFunctions[ielz+1][iely  ][ielx+1] += ShapeFunction[5];
				ShapeFunctions[ielz+1][iely+1][ielx+1] += ShapeFunction[6];
				ShapeFunctions[ielz+1][iely+1][ielx  ] += ShapeFunction[7];

				if (ComputingPhaseProportions){
					PhaseProportions_local[phase][ielz  ][iely  ][ielx  ] += ShapeFunction[0];
					PhaseProportions_local[phase][ielz  ][iely  ][ielx+1] += ShapeFunction[1];
					PhaseProportions_local[phase][ielz  ][iely+1][ielx+1] += ShapeFunction[2];
					PhaseProportions_local[phase][ielz  ][iely+1][ielx  ] += ShapeFunction[3];
					PhaseProportions_local[phase][ielz+1][iely  ][ielx  ] += ShapeFunction[4];
					PhaseProportions_local[phase][ielz+1][iely  ][ielx+1] += ShapeFunction[5];
					PhaseProportions_local[phase][ielz+1][iely+1][ielx+1] += ShapeFunction[6];
					PhaseProportions_local[phase][ielz+1][iely+1][ielx  ] += ShapeFunction[7];
				}
				else{
					MaterialProperties_local[ielz  ][iely  ][ielx  ] += PropertyValue*ShapeFunction[0];
					MaterialProperties_local[ielz  ][iely  ][ielx+1] += PropertyValue*ShapeFunction[1];
					MaterialProperties_local[ielz  ][iely+1][ielx+1] += PropertyValue*ShapeFunction[2];
					MaterialProperties_local[ielz  ][iely+1][ielx  ] += PropertyValue*ShapeFunction[3];
					MaterialProperties_local[ielz+1][iely  ][ielx  ] += PropertyValue*ShapeFunction[4];
					MaterialProperties_local[ielz+1][iely  ][ielx+1] += PropertyValue*ShapeFunction[5];
					MaterialProperties_local[ielz+1][iely+1][ielx+1] += PropertyValue*ShapeFunction[6];
					MaterialProperties_local[ielz+1][iely+1][ielx  ] += PropertyValue*ShapeFunction[7];
				}
			}										/* End of averaging corner nodes */
			else if (AveragingNode==2){				/* Average around XY nodes */

				/* Evaluate shape function at (eta,zetha,phi) */
				x = eta;
				y = zetha;
				z = 1.0-PetscAbs(phi);

				//ShapeFunction[0] =
				ShapeFunction[0]		=	0.25*(1.0-x)*(1.0-y)*z;
				ShapeFunction[1]		=	0.25*(1.0+x)*(1.0-y)*z;
				ShapeFunction[2]		=	0.25*(1.0+x)*(1.0+y)*z;
				ShapeFunction[3]		=	0.25*(1.0-x)*(1.0+y)*z;


				/* Add contributions to the 8 surrounding nodes */
				Numberparticles[ielz  ][iely  ][ielx  ] += 1;
				Numberparticles[ielz  ][iely  ][ielx+1] += 1;
				Numberparticles[ielz  ][iely+1][ielx+1] += 1;
				Numberparticles[ielz  ][iely+1][ielx  ] += 1;


				ShapeFunctions[ielz  ][iely  ][ielx  ] += ShapeFunction[0];
				ShapeFunctions[ielz  ][iely  ][ielx+1] += ShapeFunction[1];
				ShapeFunctions[ielz  ][iely+1][ielx+1] += ShapeFunction[2];
				ShapeFunctions[ielz  ][iely+1][ielx  ] += ShapeFunction[3];


				if (ComputingPhaseProportions){
					PhaseProportions_local[phase][ielz  ][iely  ][ielx  ] += ShapeFunction[0];
					PhaseProportions_local[phase][ielz  ][iely  ][ielx+1] += ShapeFunction[1];
					PhaseProportions_local[phase][ielz  ][iely+1][ielx+1] += ShapeFunction[2];
					PhaseProportions_local[phase][ielz  ][iely+1][ielx  ] += ShapeFunction[3];
				}
				else{
					MaterialProperties_local[ielz  ][iely  ][ielx  ] += PropertyValue*ShapeFunction[0];
					MaterialProperties_local[ielz  ][iely  ][ielx+1] += PropertyValue*ShapeFunction[1];
					MaterialProperties_local[ielz  ][iely+1][ielx+1] += PropertyValue*ShapeFunction[2];
					MaterialProperties_local[ielz  ][iely+1][ielx  ] += PropertyValue*ShapeFunction[3];
				}

			}
			else if (AveragingNode==3){				/* Average around XZ nodes */

				/* Evaluate shape function at (eta,zetha,phi) */
				x = eta;
				y = 1.0-PetscAbs(zetha);
				z = phi;

				ShapeFunction[0]		=	0.25*(1.0-x)*y*(1.0-z);
				ShapeFunction[1]		=	0.25*(1.0+x)*y*(1.0-z);
				ShapeFunction[2]		=	0.25*(1.0-x)*y*(1.0+z);
				ShapeFunction[3]		=	0.25*(1.0+x)*y*(1.0+z);


				/* Add contributions to the 4 surrounding nodes */
				Numberparticles[ielz  ][iely  ][ielx  ] += 1;
				Numberparticles[ielz  ][iely  ][ielx+1] += 1;
				Numberparticles[ielz+1][iely  ][ielx  ] += 1;
				Numberparticles[ielz+1][iely  ][ielx+1] += 1;


				ShapeFunctions[ielz  ][iely  ][ielx  ] += ShapeFunction[0];
				ShapeFunctions[ielz  ][iely  ][ielx+1] += ShapeFunction[1];
				ShapeFunctions[ielz+1][iely  ][ielx  ] += ShapeFunction[2];
				ShapeFunctions[ielz+1][iely  ][ielx+1] += ShapeFunction[3];


				if (ComputingPhaseProportions){
					PhaseProportions_local[phase][ielz  ][iely  ][ielx  ] += ShapeFunction[0];
					PhaseProportions_local[phase][ielz  ][iely  ][ielx+1] += ShapeFunction[1];
					PhaseProportions_local[phase][ielz+1][iely  ][ielx  ] += ShapeFunction[2];
					PhaseProportions_local[phase][ielz+1][iely  ][ielx+1] += ShapeFunction[3];
				}
				else{
					MaterialProperties_local[ielz  ][iely  ][ielx  ] += PropertyValue*ShapeFunction[0];
					MaterialProperties_local[ielz  ][iely  ][ielx+1] += PropertyValue*ShapeFunction[1];
					MaterialProperties_local[ielz+1][iely  ][ielx  ] += PropertyValue*ShapeFunction[2];
					MaterialProperties_local[ielz+1][iely  ][ielx+1] += PropertyValue*ShapeFunction[3];
				}

			}
			else if (AveragingNode==4){				/* Average around YZ nodes */

				/* Evaluate shape function at (eta,zetha,phi) */
				x = 1.0-PetscAbs(eta);
				y = zetha;
				z = phi;

				ShapeFunction[0]		=	0.25*x*(1.0-y)*(1.0-z);
				ShapeFunction[1]		=	0.25*x*(1.0+y)*(1.0-z);
				ShapeFunction[2]		=	0.25*x*(1.0-y)*(1.0+z);
				ShapeFunction[3]		=	0.25*x*(1.0+y)*(1.0+z);

				/* Add contributions to the 4 surrounding nodes */
				Numberparticles[ielz  ][iely  ][ielx  ] += 1;
				Numberparticles[ielz  ][iely+1][ielx  ] += 1;
				Numberparticles[ielz+1][iely  ][ielx  ] += 1;
				Numberparticles[ielz+1][iely+1][ielx  ] += 1;


				ShapeFunctions[ielz  ][iely  ][ielx  ] += ShapeFunction[0];
				ShapeFunctions[ielz  ][iely+1][ielx  ] += ShapeFunction[1];
				ShapeFunctions[ielz+1][iely  ][ielx  ] += ShapeFunction[2];
				ShapeFunctions[ielz+1][iely+1][ielx  ] += ShapeFunction[3];


				if (ComputingPhaseProportions){
					PhaseProportions_local[phase][ielz  ][iely  ][ielx  ] += ShapeFunction[0];
					PhaseProportions_local[phase][ielz  ][iely+1][ielx  ] += ShapeFunction[1];
					PhaseProportions_local[phase][ielz+1][iely  ][ielx  ] += ShapeFunction[2];
					PhaseProportions_local[phase][ielz+1][iely+1][ielx  ] += ShapeFunction[3];
				}
				else{
					MaterialProperties_local[ielz  ][iely  ][ielx  ] += PropertyValue*ShapeFunction[0];
					MaterialProperties_local[ielz  ][iely+1][ielx  ] += PropertyValue*ShapeFunction[1];
					MaterialProperties_local[ielz+1][iely  ][ielx  ] += PropertyValue*ShapeFunction[2];
					MaterialProperties_local[ielz+1][iely+1][ielx  ] += PropertyValue*ShapeFunction[3];
				}

			}

		}


	}
	ierr = DMDAVecRestoreArray(da,   	NumberParticles_local, 	&Numberparticles	 	); CHKERRQ(ierr); // Viscosity @ nodes
	ierr = DMDAVecRestoreArray(da,    ShapeFunctions_local,	&ShapeFunctions		); CHKERRQ(ierr);

	if (ComputingPhaseProportions){
		for (i=0; i<user->num_phases; i++ ){
			ierr = DMDAVecRestoreArray(da,   	Vec_PhaseProportions_local[i], 	&PhaseProportions_local[i]	 ); CHKERRQ(ierr);
		}
	}
	else{
		ierr = DMDAVecRestoreArray(da,   	Local_MaterialProperties_Vec, 	&MaterialProperties_local	 ); CHKERRQ(ierr); 	// will contain averaged material properties
	}


	// Send local->global vectors
	if (ComputingPhaseProportions){
		for (i=0;  i<user->num_phases; i++){
			if 		(AveragingNode==0){		// Center points
				ierr = VecZeroEntries(user->FDSTAG.Center_PhaseProportions[i]); CHKERRQ(ierr);
				ierr = DMLocalToGlobalBegin(da, Vec_PhaseProportions_local[i], ADD_VALUES, user->FDSTAG.Center_PhaseProportions[i]);	CHKERRQ(ierr);
				ierr = DMLocalToGlobalEnd  (da, Vec_PhaseProportions_local[i], ADD_VALUES, user->FDSTAG.Center_PhaseProportions[i]);	CHKERRQ(ierr);
			}
			else if (AveragingNode==1){		// Corner points
				ierr = VecZeroEntries(user->FDSTAG.Corner_PhaseProportions[i]); CHKERRQ(ierr);
				ierr = DMLocalToGlobalBegin(da, Vec_PhaseProportions_local[i], ADD_VALUES,  user->FDSTAG.Corner_PhaseProportions[i]);	CHKERRQ(ierr);
				ierr = DMLocalToGlobalEnd  (da, Vec_PhaseProportions_local[i], ADD_VALUES, user->FDSTAG.Corner_PhaseProportions[i]);	CHKERRQ(ierr);
			}
			else if (AveragingNode==2){		// XY points
				ierr = VecZeroEntries(user->FDSTAG.XYPoints_PhaseProportions[i]); CHKERRQ(ierr);
				ierr = DMLocalToGlobalBegin(da, Vec_PhaseProportions_local[i], ADD_VALUES, user->FDSTAG.XYPoints_PhaseProportions[i]);	CHKERRQ(ierr);
				ierr = DMLocalToGlobalEnd  (da, Vec_PhaseProportions_local[i], ADD_VALUES, user->FDSTAG.XYPoints_PhaseProportions[i]);	CHKERRQ(ierr);
			}
			else if (AveragingNode==3){		// XZ points
				ierr = VecZeroEntries(user->FDSTAG.XZPoints_PhaseProportions[i]); CHKERRQ(ierr);
				ierr = DMLocalToGlobalBegin(da, Vec_PhaseProportions_local[i], ADD_VALUES, user->FDSTAG.XZPoints_PhaseProportions[i]);	CHKERRQ(ierr);
				ierr = DMLocalToGlobalEnd  (da, Vec_PhaseProportions_local[i], ADD_VALUES, user->FDSTAG.XZPoints_PhaseProportions[i]);	CHKERRQ(ierr);
			}
			else if (AveragingNode==4){		// YZ points
				ierr = VecZeroEntries(user->FDSTAG.YZPoints_PhaseProportions[i]); CHKERRQ(ierr);
				ierr = DMLocalToGlobalBegin(da, Vec_PhaseProportions_local[i], ADD_VALUES, user->FDSTAG.YZPoints_PhaseProportions[i]);	CHKERRQ(ierr);
				ierr = DMLocalToGlobalEnd  (da, Vec_PhaseProportions_local[i], ADD_VALUES, user->FDSTAG.YZPoints_PhaseProportions[i]);	CHKERRQ(ierr);
			}

		}
	}
	else{
		ierr = VecZeroEntries(Global_MaterialProperties); CHKERRQ(ierr);
		ierr = DMLocalToGlobalBegin(da,Local_MaterialProperties_Vec, ADD_VALUES, Global_MaterialProperties);		CHKERRQ(ierr);
		ierr = DMLocalToGlobalEnd  (da, Local_MaterialProperties_Vec, ADD_VALUES, Global_MaterialProperties);	CHKERRQ(ierr);
	}
	ierr = VecZeroEntries(ShapeFunctions_global); CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da,ShapeFunctions_local, ADD_VALUES, ShapeFunctions_global);	CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd (da, ShapeFunctions_local, ADD_VALUES, ShapeFunctions_global);	CHKERRQ(ierr);

	ierr = VecZeroEntries(NumberParticles_global); CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da,NumberParticles_local, ADD_VALUES, NumberParticles_global);	CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd (da, NumberParticles_local, ADD_VALUES, NumberParticles_global);	CHKERRQ(ierr);

	if (AveragingNode==0){
		/* Store # of particles in center node [mainly for debugging purposes */
		ierr = VecZeroEntries(user->FDSTAG.Center_NumParticles); CHKERRQ(ierr);
		ierr = DMLocalToGlobalBegin(da, NumberParticles_local, ADD_VALUES, user->FDSTAG.Center_NumParticles);	CHKERRQ(ierr);
		ierr = DMLocalToGlobalEnd  (da, NumberParticles_local, ADD_VALUES, user->FDSTAG.Center_NumParticles);	CHKERRQ(ierr);
	}



	// Pointwise division
	if (ComputingPhaseProportions){
		for (i=0;  i<user->num_phases; i++){
			PetscScalar 	minValue, maxValue;

			if 		(AveragingNode==0){		// Center points
				ierr = VecPointwiseDivide(user->FDSTAG.Center_PhaseProportions[i], user->FDSTAG.Center_PhaseProportions[i],   ShapeFunctions_global); CHKERRQ(ierr);


				/* Catch errors */
				VecMin(user->FDSTAG.Center_PhaseProportions[i],PETSC_NULL, &minValue);
				VecMax(user->FDSTAG.Center_PhaseProportions[i],PETSC_NULL, &maxValue);

				if ((minValue<-1e-6) || (maxValue>1.1)){
					MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
					PetscPrintf(PETSC_COMM_WORLD,"rank %i, user->FDSTAG.Center_PhaseProportions[%i] outside bounds: min,max=[%f,%f] \n",rank,i,minValue,maxValue);
					SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_FP,"Outside bounds");
				}


			}
			else if (AveragingNode==1){		// Corner points
				ierr = VecPointwiseDivide(user->FDSTAG.Corner_PhaseProportions[i], user->FDSTAG.Corner_PhaseProportions[i],   ShapeFunctions_global); CHKERRQ(ierr);


				/* Catch errors */
				VecMin(user->FDSTAG.Corner_PhaseProportions[i],PETSC_NULL, &minValue);
				VecMax(user->FDSTAG.Corner_PhaseProportions[i],PETSC_NULL, &maxValue);

				if ((minValue<-1e-6) || (maxValue>1.1)){
					PetscPrintf(PETSC_COMM_WORLD,"rank %i, user->FDSTAG.Corner_PhaseProportions[%i] outside bounds: min,max=[%f,%f] \n",rank,i,minValue,maxValue);
					SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_FP,"Outside bounds");
				}

			}
			else if (AveragingNode==2){		// XY points
				ierr = VecPointwiseDivide(user->FDSTAG.XYPoints_PhaseProportions[i], user->FDSTAG.XYPoints_PhaseProportions[i],   ShapeFunctions_global); CHKERRQ(ierr);

				/* Catch errors */
				VecMin(user->FDSTAG.XYPoints_PhaseProportions[i],PETSC_NULL, &minValue);
				VecMax(user->FDSTAG.XYPoints_PhaseProportions[i],PETSC_NULL, &maxValue);

				if ((minValue<-1e-6) || (maxValue>1.1)){
					PetscPrintf(PETSC_COMM_WORLD,"rank %i, user->FDSTAG.XYPoints_PhaseProportions[%i] outside bounds: min,max=[%f,%f] \n",rank,i,minValue,maxValue);
					SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_FP,"Outside bounds");
				}

			}
			else if (AveragingNode==3){		// XZ points
				ierr = VecPointwiseDivide(user->FDSTAG.XZPoints_PhaseProportions[i], user->FDSTAG.XZPoints_PhaseProportions[i],   ShapeFunctions_global); CHKERRQ(ierr);

				/* Catch errors */
				VecMin(user->FDSTAG.XZPoints_PhaseProportions[i],PETSC_NULL, &minValue);
				VecMax(user->FDSTAG.XZPoints_PhaseProportions[i],PETSC_NULL, &maxValue);

				if ((minValue<-1e-6) || (maxValue>1.1)){
					PetscPrintf(PETSC_COMM_WORLD,"rank %i, user->FDSTAG.XZPoints_PhaseProportions[%i] outside bounds: min,max=[%f,%f] \n",rank,i,minValue,maxValue);
					SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_FP,"Outside bounds");
				}


			}
			else if (AveragingNode==4){		// YZ points
				ierr = VecPointwiseDivide(user->FDSTAG.YZPoints_PhaseProportions[i], user->FDSTAG.YZPoints_PhaseProportions[i],   ShapeFunctions_global); CHKERRQ(ierr);

				/* Catch errors */
				VecMin(user->FDSTAG.YZPoints_PhaseProportions[i],PETSC_NULL, &minValue);
				VecMax(user->FDSTAG.YZPoints_PhaseProportions[i],PETSC_NULL, &maxValue);

				if ((minValue<-1e-6) || (maxValue>1.1)){
					PetscPrintf(PETSC_COMM_WORLD,"rank %i, user->FDSTAG.YZPoints_PhaseProportions[%i] outside bounds: min,max=[%f,%f] \n",rank,i,minValue,maxValue);
					SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_FP,"Outside bounds");
				}

			}

		}
	}
	else{
		ierr = VecPointwiseDivide(Global_MaterialProperties, Global_MaterialProperties,   ShapeFunctions_global); CHKERRQ(ierr);
	}
	ierr = VecPointwiseDivide(NumberParticles_global, 		NumberParticles_global,  		ShapeFunctions_global); CHKERRQ(ierr);


	/* Correct arrays for prevent roundoff errors */
	/* Catch floating point errors, which occur if ShapeFunctions_global has a zero entry because of too few particles */
	{
		PetscScalar 	*ShapeFunctions_global_array, *DataArray, *DataPhaseProportionsArray[max_num_phases];
		PetscInt		i_el, nlocal;

		VecGetArray(ShapeFunctions_global,&ShapeFunctions_global_array);
		//VecGetOwnershipRange(ShapeFunctions_global,&start,&end);
		VecGetLocalSize(ShapeFunctions_global, &nlocal);


		/* retrieve arrays */
		if (ComputingPhaseProportions){
			for (i=0;  i<user->num_phases; i++){
				if 		(AveragingNode==0){		// Center points
					VecGetArray(user->FDSTAG.Center_PhaseProportions[i],&DataPhaseProportionsArray[i]);
				}
				else if (AveragingNode==1){		// Corner points
					VecGetArray(user->FDSTAG.Corner_PhaseProportions[i],&DataPhaseProportionsArray[i]);
				}
				else if (AveragingNode==2){		// XY points
					VecGetArray(user->FDSTAG.XYPoints_PhaseProportions[i],&DataPhaseProportionsArray[i]);
				}
				else if (AveragingNode==3){		// XZ points
					VecGetArray(user->FDSTAG.XZPoints_PhaseProportions[i],&DataPhaseProportionsArray[i]);
				}
				else if (AveragingNode==4){		// YZ points
					VecGetArray(user->FDSTAG.YZPoints_PhaseProportions[i],&DataPhaseProportionsArray[i]);
				}

			}
		}
		else{
			VecGetArray(Global_MaterialProperties,&DataArray);
		}

		/* Set values manually to zero */
		for (i_el=0; i_el<nlocal; i_el++){
			if (PetscAbs(ShapeFunctions_global_array[i_el])<1e-13){
				/* oops, we had division by zero */

				/* Depending on which array we were dealing with we have to manually set the value to zero */
				if (ComputingPhaseProportions){
					for (i=0;  i<user->num_phases; i++){
						DataPhaseProportionsArray[i][i_el] = 0;
					}
				}
				else{
					DataArray[i_el] = 0;
				}
			}
		}


		/* restore arrays */
		if (ComputingPhaseProportions){
			for (i=0;  i<user->num_phases; i++){
				if 		(AveragingNode==0){		// Center points
					VecRestoreArray(user->FDSTAG.Center_PhaseProportions[i],&DataPhaseProportionsArray[i]);
				}
				else if (AveragingNode==1){		// Corner points
					VecRestoreArray(user->FDSTAG.Corner_PhaseProportions[i],&DataPhaseProportionsArray[i]);
				}
				else if (AveragingNode==2){		// XY points
					VecRestoreArray(user->FDSTAG.XYPoints_PhaseProportions[i],&DataPhaseProportionsArray[i]);
				}
				else if (AveragingNode==3){		// XZ points
					VecRestoreArray(user->FDSTAG.XZPoints_PhaseProportions[i],&DataPhaseProportionsArray[i]);
				}
				else if (AveragingNode==4){		// YZ points
					VecRestoreArray(user->FDSTAG.YZPoints_PhaseProportions[i],&DataPhaseProportionsArray[i]);
				}

			}
		}
		else{
			VecRestoreArray(Global_MaterialProperties,&DataArray);
		}

		VecRestoreArray(ShapeFunctions_global,&ShapeFunctions_global_array);

	}


	/* Clean up */
	if (ComputingPhaseProportions){
		for (i=0;  i<user->num_phases; i++){
			ierr = DMRestoreLocalVector(da,	&Vec_PhaseProportions_local[i]);	CHKERRQ(ierr);									// local vector that will store local phase proportions.
		}
	}
	else{
		ierr = DMRestoreLocalVector(da,	&Local_MaterialProperties_Vec);	CHKERRQ(ierr);

	}

	ierr = DMRestoreLocalVector(da,	&NumberParticles_local);	CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(da,	&ShapeFunctions_local);	CHKERRQ(ierr);
	ierr = VecDestroy(&NumberParticles_global);	CHKERRQ(ierr);
	ierr = VecDestroy(&ShapeFunctions_global);	CHKERRQ(ierr);
	LaMEMDestroy2dArray(&dhds, PETSC_NULL );


	PetscFunctionReturn(0);
}
/* ==================================================================================================== */




/* Write viscosity and density fields to disc in case of FDSTAG ======================================= */
#undef __FUNCT__
#define __FUNCT__ "OutputVTK_ViscosityDensityFields_FDSTAG"
PetscErrorCode OutputVTK_ViscosityDensityFields_FDSTAG( UserContext *user, Vec Temp)
{
	PetscErrorCode 	ierr;
	PetscViewer		vv;
	PetscInt 		iphase;
	PetscBool		MATLAB;
	char 			*PhaseName = NULL;

	/* VTK or MATLAB viewer? */
	MATLAB=PETSC_FALSE;
	//	MATLAB=PETSC_TRUE;

	//PetscPrintf(PETSC_COMM_WORLD,"OutputVTK_ViscosityDensityFields_FDSTAG: writing VTK output files \n");

	/* Write center viscosity and density */
	if (MATLAB){
		PetscPrintf(PETSC_COMM_WORLD,"FDSTAG DEBUG: Writing output files in MATLAB format - Read with ReadViscosity_FDSTAG.m \n");

//		ierr = PetscViewerBinaryMatlabOpen(PETSC_COMM_WORLD,"Viscosity_Center_FDSTAG.dat",&vv); CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Viscosity_Center_FDSTAG.dat",FILE_MODE_WRITE,&vv); CHKERRQ(ierr);
		ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_BINARY_MATLAB); CHKERRQ(ierr);

		ierr = DMDASetFieldName(user->FDSTAG.DA_CENTER,0,"data");																CHKERRQ(ierr);

//		ierr = PetscViewerBinaryMatlabOutputVecDA(vv,"ViscosityCenter",user->FDSTAG.Center_EffectiveViscosity,user->FDSTAG.DA_CENTER);	CHKERRQ(ierr);
		ierr = PetscObjectSetName((PetscObject)user->FDSTAG.Center_EffectiveViscosity,"ViscosityCenter"); CHKERRQ(ierr);
		ierr = VecView(user->FDSTAG.Center_EffectiveViscosity, vv); CHKERRQ(ierr);

//		ierr = PetscViewerBinaryMatlabOutputVecDA(vv,"DensityCenter",  user->FDSTAG.Center_Density,user->FDSTAG.DA_CENTER);	CHKERRQ(ierr);
		ierr = PetscObjectSetName((PetscObject)user->FDSTAG.Center_Density,"DensityCenter"); CHKERRQ(ierr);
		ierr = VecView(user->FDSTAG.Center_Density, vv); CHKERRQ(ierr);

		for (iphase=0;  iphase<user->num_phases;  iphase++){
			asprintf(&PhaseName, "Phase_%1lld",(LLD)iphase);
//			ierr = PetscViewerBinaryMatlabOutputVecDA(vv,PhaseName,  user->FDSTAG.Center_PhaseProportions[iphase],user->FDSTAG.DA_CENTER);	CHKERRQ(ierr);
			ierr = PetscObjectSetName((PetscObject)user->FDSTAG.Center_PhaseProportions[iphase],PhaseName); CHKERRQ(ierr);
			ierr = VecView(user->FDSTAG.Center_PhaseProportions[iphase], vv); CHKERRQ(ierr);
		}

		ierr = PetscViewerDestroy(&vv); 																			CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"FDSTAG DEBUG: wrote center viscosity & density info to  Viscosity_Center_FDSTAG.dat \n");
	}
	else{
		ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "Viscosity_Center_FDSTAG.vtk", &vv);CHKERRQ(ierr);
		ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"FDSTAG DEBUG: Writing output files in VTK format - Read with Paraview \n");

		ierr = PetscObjectSetName( (PetscObject)user->FDSTAG.Center_EffectiveViscosity, "Center_Viscosity" );CHKERRQ(ierr);
		ierr = PetscObjectSetName( (PetscObject)user->FDSTAG.Center_Density, 			"Center_Density" );CHKERRQ(ierr);
		for (iphase=0;  iphase<user->num_phases;  iphase++){
			asprintf(&PhaseName, "Phase_%1lld",(LLD)iphase);
			ierr = PetscObjectSetName( (PetscObject)user->FDSTAG.Center_PhaseProportions[iphase], PhaseName );CHKERRQ(ierr);
		}


		ierr = DMView(user->FDSTAG.DA_CENTER, vv);CHKERRQ(ierr);


		ierr = VecScale(user->FDSTAG.Center_EffectiveViscosity, user->Characteristic.Viscosity);			// scale  [dimensional units]
		ierr = VecView(user->FDSTAG.Center_EffectiveViscosity, vv);CHKERRQ(ierr);
		ierr = VecScale(user->FDSTAG.Center_EffectiveViscosity, 1/user->Characteristic.Viscosity);			// unscale  [dimensional units]

		ierr = VecScale(user->FDSTAG.Center_Density, user->Characteristic.Density);			// scale  [dimensional units]
		ierr = VecView(user->FDSTAG.Center_Density, vv);CHKERRQ(ierr);
		ierr = VecScale(user->FDSTAG.Center_Density, 1/user->Characteristic.Density);		// unscale [non-dimensional units]

		for (iphase=0;  iphase<user->num_phases;  iphase++){
			ierr = VecView(user->FDSTAG.Center_PhaseProportions[iphase], vv);CHKERRQ(ierr);
		}

		ierr = PetscViewerDestroy(&vv);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"FDSTAG DEBUG: wrote center viscosity & density & phase info to  Viscosity_Center_FDSTAG.vtk \n");

	}

#if 0
	/* Write nodal viscosity  */
	if (MATLAB){

//		ierr = PetscViewerBinaryMatlabOpen(PETSC_COMM_WORLD,"Viscosity_Nodes_FDSTAG.dat",&vv); CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Viscosity_Nodes_FDSTAG.dat",FILE_MODE_WRITE,&vv); CHKERRQ(ierr);
		ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_BINARY_MATLAB); CHKERRQ(ierr);

		ierr = DMDASetFieldName(user->DA_FDSTAG_NODES,0,"data");
																		CHKERRQ(ierr);
//		ierr = PetscViewerBinaryMatlabOutputVecDA(vv,"ViscosityNodes",user->FDSTAG_Viscosity_Nodes,user->DA_FDSTAG_NODES);	CHKERRQ(ierr);
		ierr = PetscObjectSetName((PetscObject)user->FDSTAG.FDSTAG_Viscosity_Nodes,"ViscosityNodes"); CHKERRQ(ierr);
		ierr = VecView(user->FDSTAG.FDSTAG_Viscosity_Nodes, vv); CHKERRQ(ierr);

//		ierr = PetscViewerBinaryMatlabOutputVecDA(vv,"DensityNodes",user->FDSTAG_Density_Nodes,user->DA_FDSTAG_NODES);	CHKERRQ(ierr);
		ierr = PetscObjectSetName((PetscObject)user->FDSTAG.FDSTAG_Density_Nodes,"DensityNodes"); CHKERRQ(ierr);
		ierr = VecView(user->FDSTAG.FDSTAG_Density_Nodes, vv); CHKERRQ(ierr);

		ierr = PetscViewerDestroy(&vv); 																			CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"FDSTAG DEBUG: wrote nodal viscosity & density info to  Viscosity_Nodes_FDSTAG.dat \n");
	}
	else{
		ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "Viscosity_Nodes_FDSTAG.vtk", &vv);CHKERRQ(ierr);
		ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);

		ierr = PetscObjectSetName( (PetscObject)user->FDSTAG_Viscosity_Nodes, "Nodal_Viscosity" );CHKERRQ(ierr);
		ierr = DMView(user->DA_FDSTAG_NODES, vv);CHKERRQ(ierr);

		ierr = VecScale(user->FDSTAG_Viscosity_Nodes, user->Characteristic.Viscosity);			// scale  [dimensional units]
		ierr = VecView(user->FDSTAG_Viscosity_Nodes, vv);CHKERRQ(ierr);
		ierr = VecScale(user->FDSTAG_Viscosity_Nodes, 1/user->Characteristic.Viscosity);		// unscale  [dimensional units]

		ierr = VecScale(user->FDSTAG_Density_Nodes, user->Characteristic.Density);				// scale  [dimensional units]
		ierr = VecView(user->FDSTAG_Density_Nodes, vv);CHKERRQ(ierr);
		ierr = VecScale(user->FDSTAG_Density_Nodes, 1/user->Characteristic.Density);			// unscale

		ierr = PetscViewerDestroy(&vv);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"FDSTAG DEBUG: wrote nodal viscosity , density info to  Viscosity_Nodes_FDSTAG.vtk \n");

	}


	/* Write density @ Vz points */
	if (MATLAB){
//		ierr = PetscViewerBinaryMatlabOpen(PETSC_COMM_WORLD,"Density_Vz_FDSTAG.dat",&vv); CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Density_Vz_FDSTAG.dat",FILE_MODE_WRITE,&vv); CHKERRQ(ierr);
		ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_BINARY_MATLAB); CHKERRQ(ierr);

		ierr = DMDASetFieldName(user->DA_FDSTAG_RHO_VZ,0,"data");																CHKERRQ(ierr);
//		ierr = PetscViewerBinaryMatlabOutputVecDA(vv,"VZ_DENSITY",  user->FDSTAG_Density_Vz,user->DA_FDSTAG_RHO_VZ);	CHKERRQ(ierr);
		ierr = PetscObjectSetName((PetscObject)user->FDSTAG.FDSTAG_Density_Vz,"VZ_DENSITY"); CHKERRQ(ierr);
		ierr = VecView(user->FDSTAG.FDSTAG_Density_Vz, vv); CHKERRQ(ierr);

		ierr = PetscViewerDestroy(&vv); 																			CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"FDSTAG DEBUG: wrote center viscosity & density info to  Density_Vz_FDSTAG.dat \n");
	}
	else{
		ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "Density_Vz_FDSTAG.vtk", &vv);CHKERRQ(ierr);
		ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
		ierr = PetscObjectSetName( (PetscObject)user->FDSTAG_Density_Vz, "VZ_DENSITY" );CHKERRQ(ierr);
		ierr = DMView(user->DA_FDSTAG_RHO_VZ, vv);CHKERRQ(ierr);

		ierr = VecScale(user->FDSTAG_Density_Vz, user->Characteristic.Density);				// scale  [dimensional units]
		ierr = VecView(user->FDSTAG_Density_Vz, vv);CHKERRQ(ierr);
		ierr = VecScale(user->FDSTAG_Density_Vz, 1/user->Characteristic.Density);				// unscale  [dimensional units]

		ierr = PetscViewerDestroy(&vv);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"FDSTAG DEBUG: wrote Vz density info to  Density_Vz_FDSTAG.vtk \n");
	}
#endif


	/* Write viscosity at Sxy points */
	if (MATLAB){
//		ierr = PetscViewerBinaryMatlabOpen(PETSC_COMM_WORLD,"Viscosity_XY_FDSTAG.dat",&vv); CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Viscosity_XY_FDSTAG.dat",FILE_MODE_WRITE,&vv); CHKERRQ(ierr);
		ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_BINARY_MATLAB); CHKERRQ(ierr);

		ierr = DMDASetFieldName(user->FDSTAG.DA_XY_POINTS,0,"data");																CHKERRQ(ierr);
		
//		ierr = PetscViewerBinaryMatlabOutputVecDA(vv,"Viscosity_XY",  user->FDSTAG.XYPoints_EffectiveViscosity,user->FDSTAG.DA_XY_POINTS);	CHKERRQ(ierr);
		ierr = PetscObjectSetName((PetscObject)user->FDSTAG.XYPoints_EffectiveViscosity,"Viscosity_XY"); CHKERRQ(ierr);
		ierr = VecView(user->FDSTAG.XYPoints_EffectiveViscosity, vv); CHKERRQ(ierr);

		for (iphase=0;  iphase<user->num_phases;  iphase++){
			asprintf(&PhaseName, "Phase_%1lld",(LLD)iphase);
			
//			ierr = PetscViewerBinaryMatlabOutputVecDA(vv,PhaseName,  user->FDSTAG.XYPoints_PhaseProportions[iphase],user->FDSTAG.DA_XY_POINTS);	CHKERRQ(ierr);
			ierr = PetscObjectSetName((PetscObject)user->FDSTAG.XYPoints_PhaseProportions[iphase],"PhaseName"); CHKERRQ(ierr);
			ierr = VecView(user->FDSTAG.XYPoints_PhaseProportions[iphase], vv); CHKERRQ(ierr);

		}

		ierr = PetscViewerDestroy(&vv); 																			CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"FDSTAG DEBUG: wrote center viscosity & density info to  Viscosity_XY_FDSTAG.dat \n");
	}
	else{
		ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "Viscosity_XY_FDSTAG.vtk", &vv);CHKERRQ(ierr);
		ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
		ierr = PetscObjectSetName( (PetscObject)user->FDSTAG.XYPoints_EffectiveViscosity, "Viscosity_XY" );CHKERRQ(ierr);
		for (iphase=0;  iphase<user->num_phases;  iphase++){
			asprintf(&PhaseName, "Phase_%1lld",(LLD)iphase);
			ierr = PetscObjectSetName( (PetscObject)user->FDSTAG.XYPoints_PhaseProportions[iphase], PhaseName );CHKERRQ(ierr);
		}

		ierr = DMView(user->FDSTAG.DA_XY_POINTS, vv);CHKERRQ(ierr);

		ierr = VecScale(user->FDSTAG.XYPoints_EffectiveViscosity, user->Characteristic.Viscosity);		// scale
		ierr = VecView(user->FDSTAG.XYPoints_EffectiveViscosity, vv);CHKERRQ(ierr);
		ierr = VecScale(user->FDSTAG.XYPoints_EffectiveViscosity, 1/user->Characteristic.Viscosity);	// unscale
		for (iphase=0;  iphase<user->num_phases;  iphase++){
			ierr = VecView(user->FDSTAG.XYPoints_PhaseProportions[iphase], vv);CHKERRQ(ierr);
		}

		ierr = PetscViewerDestroy(&vv);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"FDSTAG DEBUG: wrote corner viscosity and phase info to  Viscosity_XY_FDSTAG.vtk \n");
	}

	/* Write viscosity at Sxz points */
	if (MATLAB){
//		ierr = PetscViewerBinaryMatlabOpen(PETSC_COMM_WORLD,"Viscosity_XZ_FDSTAG.dat",&vv); CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Viscosity_XZ_FDSTAG.dat",FILE_MODE_WRITE,&vv); CHKERRQ(ierr);
		ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_BINARY_MATLAB); CHKERRQ(ierr);

		ierr = DMDASetFieldName(user->FDSTAG.DA_XZ_POINTS,0,"data");
																			CHKERRQ(ierr);
//		ierr = PetscViewerBinaryMatlabOutputVecDA(vv,"Viscosity_XZ",  user->FDSTAG.XZPoints_EffectiveViscosity,user->FDSTAG.DA_XZ_POINTS);	CHKERRQ(ierr);
		ierr = PetscObjectSetName((PetscObject)user->FDSTAG.XZPoints_EffectiveViscosity,"Viscosity_XZ"); CHKERRQ(ierr);
		ierr = VecView(user->FDSTAG.XZPoints_EffectiveViscosity, vv); CHKERRQ(ierr);

		for (iphase=0;  iphase<user->num_phases;  iphase++){
			asprintf(&PhaseName, "Phase_%1lld",(LLD)iphase);
			
//			ierr = PetscViewerBinaryMatlabOutputVecDA(vv,PhaseName,  user->FDSTAG.XZPoints_PhaseProportions[iphase],user->FDSTAG.DA_XZ_POINTS);	CHKERRQ(ierr);
			ierr = PetscObjectSetName((PetscObject)user->FDSTAG.XZPoints_PhaseProportions[iphase],PhaseName); CHKERRQ(ierr);
			ierr = VecView(user->FDSTAG.XZPoints_PhaseProportions[iphase], vv); CHKERRQ(ierr);

		}
		ierr = PetscViewerDestroy(&vv); 																			CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"FDSTAG DEBUG: wrote center viscosity & density info to  Viscosity_XZ_FDSTAG.dat \n");
	}
	else{
		ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "Viscosity_XZ_FDSTAG.vtk", &vv);CHKERRQ(ierr);
		ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
		ierr = PetscObjectSetName( (PetscObject)user->FDSTAG.XZPoints_EffectiveViscosity, "XZ_Viscosity" );CHKERRQ(ierr);
		for (iphase=0;  iphase<user->num_phases;  iphase++){
			asprintf(&PhaseName, "Phase_%1lld",(LLD)iphase);
			ierr = PetscObjectSetName( (PetscObject)user->FDSTAG.XZPoints_PhaseProportions[iphase], PhaseName );CHKERRQ(ierr);
		}

		ierr = DMView(user->FDSTAG.DA_XZ_POINTS, vv);CHKERRQ(ierr);

		ierr = VecScale(user->FDSTAG.XZPoints_EffectiveViscosity, user->Characteristic.Viscosity);		// scale
		ierr = VecView(user->FDSTAG.XZPoints_EffectiveViscosity, vv);CHKERRQ(ierr);
		ierr = VecScale(user->FDSTAG.XZPoints_EffectiveViscosity, 1/user->Characteristic.Viscosity);	// unscale
		for (iphase=0;  iphase<user->num_phases;  iphase++){
			ierr = VecView(user->FDSTAG.XZPoints_PhaseProportions[iphase], vv);CHKERRQ(ierr);
		}

		ierr = PetscViewerDestroy(&vv);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"FDSTAG DEBUG: wrote corner viscosity and phase info to  Viscosity_XZ_FDSTAG.vtk \n");
	}


	/* Write viscosity at Syz points */
	if (MATLAB){
//		ierr = PetscViewerBinaryMatlabOpen(PETSC_COMM_WORLD,"Viscosity_YZ_FDSTAG.dat",&vv); CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Viscosity_YZ_FDSTAG.dat",FILE_MODE_WRITE,&vv); CHKERRQ(ierr);
		ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_BINARY_MATLAB); CHKERRQ(ierr);
		
		ierr = DMDASetFieldName(user->FDSTAG.DA_YZ_POINTS,0,"data");																CHKERRQ(ierr);

//		ierr = PetscViewerBinaryMatlabOutputVecDA(vv,"Viscosity_YZ",  user->FDSTAG.YZPoints_EffectiveViscosity,user->FDSTAG.DA_YZ_POINTS);	CHKERRQ(ierr);
		ierr = PetscObjectSetName((PetscObject)user->FDSTAG.YZPoints_EffectiveViscosity,"Viscosity_YZ"); CHKERRQ(ierr);
		ierr = VecView(user->FDSTAG.YZPoints_EffectiveViscosity, vv); CHKERRQ(ierr);

		for (iphase=0;  iphase<user->num_phases;  iphase++){
			asprintf(&PhaseName, "Phase_%1lld",(LLD)iphase);
//			ierr = PetscViewerBinaryMatlabOutputVecDA(vv,PhaseName,  user->FDSTAG.YZPoints_PhaseProportions[iphase],user->FDSTAG.DA_YZ_POINTS);	CHKERRQ(ierr);
			ierr = PetscObjectSetName((PetscObject)user->FDSTAG.YZPoints_PhaseProportions[iphase],PhaseName); CHKERRQ(ierr);
			ierr = VecView(user->FDSTAG.YZPoints_PhaseProportions[iphase], vv); CHKERRQ(ierr);

		}
		ierr = PetscViewerDestroy(&vv); 																			CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"FDSTAG DEBUG: wrote center viscosity & density info to  Viscosity_YZ_FDSTAG.dat \n");
	}
	else{
		ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "Viscosity_YZ_FDSTAG.vtk", &vv);CHKERRQ(ierr);
		ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
		ierr = PetscObjectSetName( (PetscObject)user->FDSTAG.YZPoints_EffectiveViscosity, "YZ_Viscosity" );CHKERRQ(ierr);
		for (iphase=0;  iphase<user->num_phases;  iphase++){
			asprintf(&PhaseName, "Phase_%1lld",(LLD)iphase);
			ierr = PetscObjectSetName( (PetscObject)user->FDSTAG.YZPoints_PhaseProportions[iphase], PhaseName );CHKERRQ(ierr);
		}

		ierr = DMView(user->FDSTAG.DA_YZ_POINTS, vv);CHKERRQ(ierr);

		ierr = VecScale(user->FDSTAG.YZPoints_EffectiveViscosity, user->Characteristic.Viscosity);		// scale
		ierr = VecView(user->FDSTAG.YZPoints_EffectiveViscosity, vv);CHKERRQ(ierr);
		ierr = VecScale(user->FDSTAG.YZPoints_EffectiveViscosity, 1/user->Characteristic.Viscosity);	// unscale
		for (iphase=0;  iphase<user->num_phases;  iphase++){
			ierr = VecView(user->FDSTAG.YZPoints_PhaseProportions[iphase], vv);CHKERRQ(ierr);
		}

		ierr = PetscViewerDestroy(&vv);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"FDSTAG DEBUG: wrote corner viscosity and phase info to  Viscosity_YZ_FDSTAG.vtk \n");
	}




#ifdef TEMPERATURE
	if (MATLAB){
//		ierr = PetscViewerBinaryMatlabOpen(PETSC_COMM_WORLD,"EnergyProperties_FDSTAG.dat",&vv); CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"EnergyProperties_FDSTAG.dat",FILE_MODE_WRITE,&vv); CHKERRQ(ierr);
		ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_BINARY_MATLAB); CHKERRQ(ierr);

		ierr = DMDASetFieldName(user->FDSTAG.DA_CORNER,0,"data");																				CHKERRQ(ierr);
		
//		ierr = PetscViewerBinaryMatlabOutputVecDA(vv,"HeatCapacity",	user->FDSTAG.Corner_HeatCapacity,user->FDSTAG.DA_CORNER);				CHKERRQ(ierr);
		ierr = PetscObjectSetName((PetscObject)user->FDSTAG.Corner_HeatCapacity,"HeatCapacity"); CHKERRQ(ierr);
		ierr = VecView(user->FDSTAG.Corner_HeatCapacity, vv); CHKERRQ(ierr);

//		ierr = PetscViewerBinaryMatlabOutputVecDA(vv,"Conductivity",	user->FDSTAG.Corner_Conductivity,user->FDSTAG.DA_CORNER);				CHKERRQ(ierr);
		ierr = PetscObjectSetName((PetscObject)user->FDSTAG.Corner_Conductivity,"Conductivity"); CHKERRQ(ierr);
		ierr = VecView(user->FDSTAG.Corner_Conductivity, vv); CHKERRQ(ierr);

//		ierr = PetscViewerBinaryMatlabOutputVecDA(vv,"Radioactiveheat",	user->FDSTAG.Corner_RadioactiveHeat,user->FDSTAG.DA_CORNER);	CHKERRQ(ierr);
		ierr = PetscObjectSetName((PetscObject)user->FDSTAG.Corner_RadioactiveHeat,"Radioactiveheat"); CHKERRQ(ierr);
		ierr = VecView(user->FDSTAG.Corner_RadioactiveHeat, vv); CHKERRQ(ierr);

//		ierr = PetscViewerBinaryMatlabOutputVecDA(vv,"Temperature",	Temp,user->DA_Temp);	CHKERRQ(ierr);
		ierr = PetscObjectSetName((PetscObject)Temp,"Temperature"); CHKERRQ(ierr);
		ierr = VecView(Temp, vv); CHKERRQ(ierr);

		for (iphase=0;  iphase<user->num_phases;  iphase++){
			asprintf(&PhaseName, "Phase_%1lld",(LLD)iphase);
//			ierr = PetscViewerBinaryMatlabOutputVecDA(vv,PhaseName,  user->FDSTAG.Corner_PhaseProportions[iphase],user->FDSTAG.DA_CORNER);	CHKERRQ(ierr);
			ierr = PetscObjectSetName((PetscObject)user->FDSTAG.Corner_PhaseProportions[iphase],PhaseName); CHKERRQ(ierr);
			ierr = VecView(user->FDSTAG.Corner_PhaseProportions[iphase], vv); CHKERRQ(ierr);

		}
		ierr = PetscViewerDestroy(&vv); 																			CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"FDSTAG DEBUG: wrote nodal heat capacity, conductivity and radioactive heatinfo to  EnergyProperties_FDSTAG.dat \n");
	}
	else{
		ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "EnergyProperties_FDSTAG.vtk", &vv);CHKERRQ(ierr);
		ierr = PetscViewerSetFormat(vv, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);

		ierr = PetscObjectSetName( (PetscObject)user->FDSTAG.Corner_HeatCapacity, "HeatCapacity" );CHKERRQ(ierr);
		ierr = PetscObjectSetName( (PetscObject)user->FDSTAG.Corner_Conductivity, "Conductivity" );CHKERRQ(ierr);
		ierr = PetscObjectSetName( (PetscObject)user->FDSTAG.Corner_RadioactiveHeat, "Radioactiveheat" );CHKERRQ(ierr);
		ierr = PetscObjectSetName( (PetscObject)Temp, "Temperature" );CHKERRQ(ierr);
		for (iphase=0;  iphase<user->num_phases;  iphase++){
			asprintf(&PhaseName, "Phase_%1lld",(LLD)iphase);
			ierr = PetscObjectSetName( (PetscObject)user->FDSTAG.Corner_PhaseProportions[iphase], PhaseName );CHKERRQ(ierr);
		}

		ierr = DMView(user->FDSTAG.DA_CORNER, vv);CHKERRQ(ierr);
		ierr = VecView(user->FDSTAG.Corner_HeatCapacity, vv);CHKERRQ(ierr);
		ierr = VecView(user->FDSTAG.Corner_Conductivity, vv);CHKERRQ(ierr);
		ierr = VecView(user->FDSTAG.Corner_RadioactiveHeat, vv);CHKERRQ(ierr);
		ierr = VecView(Temp, vv);CHKERRQ(ierr);

		for (iphase=0;  iphase<user->num_phases;  iphase++){
			ierr = VecView(user->FDSTAG.Corner_PhaseProportions[iphase], vv);CHKERRQ(ierr);
		}

		ierr = PetscViewerDestroy(&vv);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"FDSTAG DEBUG: wrote nodal heat capacity, conductivity and radioactive heat, T and phases info to  EnergyProperties_FDSTAG.vtk \n");


	}

#endif


	PetscFunctionReturn(0);

}
/* ==================================================================================================== */

/* Write initial Particles to disc =================================================================== */
#undef __FUNCT__
#define __FUNCT__ "WriteInitialParticlesToDisc"
PetscErrorCode WriteInitialParticlesToDisc( UserContext *user)
{
	PetscMPIInt			rank, size;
	PetscErrorCode 		ierr;
	char				SaveFileNameParticles[PETSC_MAX_PATH_LEN];
	PetscViewer			view_out;
	Vec 				Particle_Vec, information;


	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&size);

	/* Create directory */
	if (rank==0){
		mkdir("InitialParticles", S_IRWXU);			// only rank 0 does this
	}
	MPI_Barrier(PETSC_COMM_WORLD);					// All other procs should wait


	sprintf(SaveFileNameParticles  ,"%s.%lld.out",	"./InitialParticles/Initial_Particles"				, (LLD)rank	); 	// construct the filename with Particles


	ierr = VecCreate(PETSC_COMM_SELF,&information); CHKERRQ(ierr);
	ierr = VecSetSizes(information,PETSC_DECIDE,40); CHKERRQ(ierr);
	ierr = VecSetFromOptions(information); CHKERRQ(ierr);
	ierr = VecSetValue(information,0 ,	(PetscScalar)(size), 						INSERT_VALUES); CHKERRQ(ierr);	// # of processors
	ierr = VecSetValue(information,1 ,	(PetscScalar)(user->num_particle_local), 	INSERT_VALUES); CHKERRQ(ierr);	// # of local particles
	ierr = VecSetValue(information,2 ,	(PetscScalar)(particle_props), 			INSERT_VALUES); CHKERRQ(ierr);	// # of particle properties


	ierr = VecAssemblyBegin(information);  CHKERRQ(ierr);
	ierr = VecAssemblyEnd(information); CHKERRQ(ierr);

	/* Put all particles (and corresponding data such as T, P, location etc.) in a single vector */

	// Ad-hoc solution; shoud be replaced
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,
			1,
			particle_props*user->num_particle_local,
			(PetscScalar*)(user->ParticlesLocal),
			&Particle_Vec); CHKERRQ(ierr);


	ierr = VecAssemblyBegin(Particle_Vec); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(Particle_Vec); CHKERRQ(ierr);

	/* write file */
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,SaveFileNameParticles,FILE_MODE_WRITE, &view_out); CHKERRQ(ierr);
	ierr = VecView(information, 		view_out); 				CHKERRQ(ierr);
	ierr = VecView(Particle_Vec, 		view_out); 				CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&view_out); 						CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD,"Saved initial particle distribution to directory ./InitialParticles \n");


	/* Clean up */
	ierr = VecDestroy(&information); 	 	CHKERRQ(ierr);
	ierr = VecDestroy(&Particle_Vec); 		CHKERRQ(ierr);

	PetscFunctionReturn(0);

}
/* ==================================================================================================== */


/* Load initial Particles from disc =================================================================== */

#undef __FUNCT__
#define __FUNCT__ "LoadInitialParticlesFromDisc_FDSTAG"
PetscErrorCode LoadInitialParticlesFromDisc_FDSTAG( UserContext *user)
{
	PetscMPIInt			rank, size;
	PetscErrorCode 		ierr;
	PetscInt			ipart,ix[1];
	PetscInt	 		num0, num_total,AddRandomNoiseParticles;
	PetscScalar 		num[1];
	char				LoadFileNameParticles[PETSC_MAX_PATH_LEN];
	PetscViewer			view_in;
	Vec 				lvec_prtcls,lvec_info;
	Vec 				lvec_load;
	PetscScalar         *larray_prtcls;
	PetscScalar			*data,*larray_info;
	PetscInt            num_props=0,num_entries,num_procs=0;
	PetscInt			num_prtcls,Nprocx,Nprocy,Nprocz;
	PetscInt			nump_x,nump_y,nump_z,npartcls_x,npartcls_y,npartcls_z,nx=0,ny=0,nz=0;
	PetscInt            xs,xm,ys,ym,zs,zm,k;
	PetscRandom			rctx;
	PetscScalar	    	cf_rand;
	PetscScalar			dx,dy,dz;
	DM                  cda;
	Vec                 gc;
	PetscScalar         *xc,*yc,*zc;
	DMDACoor3d       ***coors;

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&size);


	// [0]. Load particles from InitialParticles directory for FDSTAG
	if (user->LoadInitialParticlesFromDisc==1){
		PetscPrintf(PETSC_COMM_WORLD,"Load initial particles \n");
		// load particles
		sprintf(LoadFileNameParticles  ,"./%s/%s.%lld.out",user->LoadInitialParticlesDirectory,	"Initial_Particles",(LLD)rank);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,LoadFileNameParticles,FILE_MODE_READ, &view_in); CHKERRQ(ierr);
		ierr = VecCreate(PETSC_COMM_SELF,&lvec_info);                 CHKERRQ(ierr);
		ierr = VecCreate(PETSC_COMM_SELF,&lvec_prtcls);               CHKERRQ(ierr);
		ierr = VecLoad(lvec_info,  view_in);                          CHKERRQ(ierr);
		ierr = VecLoad(lvec_prtcls, view_in);                         CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&view_in);                          CHKERRQ(ierr);


		// Error checking
		ix[0] = 0; VecGetValues(lvec_info,1,ix,num); num_procs = (PetscInt)(num[0]);
		if (size != num_procs){
			SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_SUP,"InitialParticles were created on %lld procs; yet you currently run with %lld procs!",(LLD)num_procs, (LLD)size);
		}
		ix[0] = 1; VecGetValues(lvec_info,1,ix,num);
		ix[0] = 2; VecGetValues(lvec_info,1,ix,num); num_props  = (PetscInt)(num[0]);


		ierr = VecGetSize(lvec_prtcls, &num_entries);                 CHKERRQ(ierr);
		user->num_particle_local = num_entries/particle_props;

		ierr = VecGetArray(lvec_prtcls,&larray_prtcls);               CHKERRQ(ierr);
	}

	// [1]. Load particles from Matlab generated files
	if (user->LoadInitialParticlesFromDisc==2){
		// Loading particle information from Matlab generated files - recommended for parallel runs
		// You can add random noise to particles (after they were loaded from Matlab) by setting AddRandomNoiseParticles = 1 in the ParamFile;
		// load files
		PetscPrintf(PETSC_COMM_WORLD,"# Loading particle information from MATLAB files... \n");

		sprintf(LoadFileNameParticles  ,"%s.%lld.out",	"./MatlabInputParticles/Input_Particles",(LLD)rank);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,LoadFileNameParticles,FILE_MODE_READ, &view_in); CHKERRQ(ierr);
		ierr = VecCreate(PETSC_COMM_SELF,&lvec_load);                 CHKERRQ(ierr);
		ierr = VecLoad(lvec_load,  view_in);                          CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&view_in);                          CHKERRQ(ierr);

		ierr = VecGetArray(lvec_load,&data); CHKERRQ(ierr);
		ierr = VecGetSize(lvec_load, &num_total);                 CHKERRQ(ierr);

		// split loaded data into 2 arrays larray_info and larray_prtcls
		ierr = PetscMalloc((size_t)12*sizeof(PetscScalar), &larray_info); CHKERRQ(ierr);
		ierr = PetscMalloc((size_t)(num_total-12)*sizeof(PetscScalar), &larray_prtcls); CHKERRQ(ierr);

		for(k=0;k<12;k++){
			larray_info[k]=data[k];
		}

		num0=0;
		for(k=12;k<num_total;k++){
			larray_prtcls[num0]=data[k];
			num0++;
		}

		ierr = VecRestoreArray(lvec_load,&data); CHKERRQ(ierr);
		ierr = VecDestroy(&lvec_load); 		                          CHKERRQ(ierr);

		/* Note: this method is consistent with creating initial particles files. The 2 arrays larray_info
		   and larray_prtcls contain the needed information. */

		// ERROR CHECKING - load information from larray_info
		num_procs  = (PetscInt)(larray_info[0]); 	//no of processors
		num_prtcls = (PetscInt)(larray_info[1]);	//no of particles per partition
		num_props  = (PetscInt)(larray_info[2]);  	//this might be different from the num_prop=28 that is used for storing particle properties
		Nprocx     = (PetscInt)(larray_info[3]);	//no of proc x dir
		Nprocy     = (PetscInt)(larray_info[4]);	//no of proc y dir
		Nprocz     = (PetscInt)(larray_info[5]);	//no of proc z dir
		nump_x     = (PetscInt)(larray_info[6]);	//no of particles x dir
		nump_y     = (PetscInt)(larray_info[7]);	//no of particles y dir
		nump_z     = (PetscInt)(larray_info[8]);	//no of particles z dir
		npartcls_x = (PetscInt)(larray_info[9]);	//no of particles/cell x dir
		npartcls_y = (PetscInt)(larray_info[10]);	//no of particles/cell y dir
		npartcls_z = (PetscInt)(larray_info[11]);	//no of particles/cell z dir
		nx		   = nump_x/npartcls_x;
		ny		   = nump_y/npartcls_y;
		nz		   = nump_z/npartcls_z;

		if (size != num_procs){
			SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_SUP,"Input particles were created on %lld procs in Matlab; yet you currently run with %lld procs!",(LLD)num_procs, (LLD)size);
		}
		if (user->nel_x!= nx){
			SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_SUP,"Model setup was created with %lld cells in the X-dir, while here you require %lld cells in the X-dir!",(LLD)nx, (LLD)user->nel_x);
		}
		if (user->nel_y!= ny){
			SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_SUP,"Model setup was created with %lld cells in the Y-dir, while here you require %lld cells in the Y-dir!",(LLD)ny, (LLD)user->nel_y);
		}
		if (user->nel_z!= nz){
			SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_SUP,"Model setup was created with %lld cells in the Z-dir, while here you require %lld cells in the Z-dir!",(LLD)nz, (LLD)user->nel_z);
		}
		if (user->cpu_x!= Nprocx){
			SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_SUP,"Processor partitioning error: particles created on %lld procs in X-dir; yet you currently run on %lld procs in X-dir!",(LLD)Nprocx, (LLD)user->cpu_x);
		}
		if (user->cpu_y!= Nprocy){
			SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_SUP,"Processor partitioning error: particles created on %lld procs in Y-dir; yet you currently run on %lld procs in Y-dir!",(LLD)Nprocy, (LLD)user->cpu_y);
		}
		if (user->cpu_z!= Nprocz){
			SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_SUP,"Processor partitioning error: particles created on %lld procs in Z-dir; yet you currently run on %lld procs in Z-dir!",(LLD)Nprocz, (LLD)user->cpu_z);
		}
		if (num0/num_props!=num_prtcls){
			SETERRQ3(PETSC_COMM_WORLD, PETSC_ERR_SUP,"Matlab error: No of particles present (%lld) is inconsistent with the no of particles required (%lld) on processor [%lld]",(LLD)(num0/num_props),(LLD)num_prtcls,(LLD)rank);
		}
		if (user->NumPartX!= npartcls_x){
			SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_SUP,"Model setup was created with %lld particles/cell in the X-dir, while here you require %lld particles/cell in the X-dir!",(LLD)npartcls_x, (LLD)user->NumPartX);
		}
		if (user->NumPartY!= npartcls_y){
			SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_SUP,"Model setup was created with %lld particles/cell in the Y-dir, while here you require %lld particles/cell in the Y-dir!",(LLD)npartcls_y, (LLD)user->NumPartY);
		}
		if (user->NumPartZ!= npartcls_z){
			SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_SUP,"Model setup was created with %lld particles/cell in the Z-dir, while here you require %lld particles/cell in the Z-dir!",(LLD)npartcls_z, (LLD)user->NumPartZ);
		}

		user->num_particle_local = num0/num_props;
	}

	// Get coordinates axes of finest nodes
	ierr = DMDAGetCorners(user->DA_Vel,&xs,&ys,&zs,&xm,&ym,&zm);  CHKERRQ(ierr);
	ierr = DMGetCoordinates(user->DA_Vel , &gc);                CHKERRQ(ierr);
	ierr = DMGetCoordinateDM(user->DA_Vel, &cda);               CHKERRQ(ierr);

	ierr = PetscMalloc((size_t)xm*sizeof(PetscScalar),&xc);               CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)ym*sizeof(PetscScalar),&yc);               CHKERRQ(ierr);
	ierr = PetscMalloc((size_t)zm*sizeof(PetscScalar),&zc);               CHKERRQ(ierr);

	ierr = DMDAVecGetArray(cda, gc, &coors);                      CHKERRQ(ierr);  
	for(k=xs;k<xs+xm;k++) xc[k-xs] = coors[zs][ys][k].x;
	for(k=ys;k<ys+ym;k++) yc[k-ys] = coors[zs][k][xs].y;
	for(k=zs;k<zs+zm;k++) zc[k-zs] = coors[k][ys][xs].z;
	ierr = DMDAVecRestoreArray(cda, gc, &coors);                  CHKERRQ(ierr);


	//Allocate the local particles array (which has a slight buffer)
	user->MaxNumLocalParticles = (PetscInt) (1.20*(PetscScalar)user->num_particle_local);
	ierr = PetscMalloc((size_t)user->MaxNumLocalParticles*sizeof(Particles),&user->ParticlesLocal); CHKERRQ(ierr);
	ierr = PetscMemzero(user->ParticlesLocal,(size_t)user->MaxNumLocalParticles*sizeof(Particles)); CHKERRQ(ierr);

	for (ipart=0; ipart<user->num_particle_local; ipart++){
		user->ParticlesLocal[ipart].x             = larray_prtcls[ipart*num_props+0];
		user->ParticlesLocal[ipart].y             = larray_prtcls[ipart*num_props+1];
		user->ParticlesLocal[ipart].z             = larray_prtcls[ipart*num_props+2];
		user->ParticlesLocal[ipart].num           = larray_prtcls[ipart*num_props+3];
		user->ParticlesLocal[ipart].phase         = larray_prtcls[ipart*num_props+4];


		// why only from 0 to xm-2? Is this correct? Or do we have a mistake in GetParticleNodes (line 965) ?



		if (num_props==28){
			// if we load initial particles that were dumped to disk before by LaMEM, we should read them all here in the correct
			// order.
			//PetscScalar x,y,z,num,phase,cpu,ix,iy,iz,eta,zetha,phi,T,P,PlasticStrain,Strain, Txx,Tyy,Tzz,Txy,Txz,Tyz,E2nd,T2nd,E2nd_pl, Mu_eff, Plastic, Mu_viscous;

			// Read the other parameters in the order that they are stored within the structure:
			user->ParticlesLocal[ipart].cpu           = larray_prtcls[ipart*num_props+5];
			user->ParticlesLocal[ipart].ix            = larray_prtcls[ipart*num_props+6];
			user->ParticlesLocal[ipart].iy            = larray_prtcls[ipart*num_props+7];
			user->ParticlesLocal[ipart].iz            = larray_prtcls[ipart*num_props+8];
			user->ParticlesLocal[ipart].eta           = larray_prtcls[ipart*num_props+9];
			user->ParticlesLocal[ipart].zetha         = larray_prtcls[ipart*num_props+10];
			user->ParticlesLocal[ipart].phi           = larray_prtcls[ipart*num_props+11];
			user->ParticlesLocal[ipart].T             = larray_prtcls[ipart*num_props+12];
			user->ParticlesLocal[ipart].P             = larray_prtcls[ipart*num_props+13];
			user->ParticlesLocal[ipart].PlasticStrain = larray_prtcls[ipart*num_props+14];
			user->ParticlesLocal[ipart].Strain        = larray_prtcls[ipart*num_props+15];
			user->ParticlesLocal[ipart].Txx           = larray_prtcls[ipart*num_props+16];
			user->ParticlesLocal[ipart].Tyy           = larray_prtcls[ipart*num_props+17];
			user->ParticlesLocal[ipart].Tzz           = larray_prtcls[ipart*num_props+18];
			user->ParticlesLocal[ipart].Txy           = larray_prtcls[ipart*num_props+19];
			user->ParticlesLocal[ipart].Txz           = larray_prtcls[ipart*num_props+20];
			user->ParticlesLocal[ipart].Tyz           = larray_prtcls[ipart*num_props+21];
			user->ParticlesLocal[ipart].E2nd          = larray_prtcls[ipart*num_props+22];
			user->ParticlesLocal[ipart].T2nd          = larray_prtcls[ipart*num_props+23];
			user->ParticlesLocal[ipart].E2nd_pl       = larray_prtcls[ipart*num_props+24];
			user->ParticlesLocal[ipart].Mu_eff        = larray_prtcls[ipart*num_props+25];
			user->ParticlesLocal[ipart].Plastic       = larray_prtcls[ipart*num_props+26];
			user->ParticlesLocal[ipart].Mu_viscous    = larray_prtcls[ipart*num_props+27];
		}
		else if (num_props==5 || num_props==6){
			// If, on the other hand, the initial particles were generated by MATLAB, only 5 properties were saved to disk, and the ix,iy,iz can be reconstructed
			user->ParticlesLocal[ipart].ix            = (PetscScalar) (xs+Bisection(xc,0,xm-2,user->ParticlesLocal[ipart].x));
			user->ParticlesLocal[ipart].iy            = (PetscScalar) (ys+Bisection(yc,0,ym-2,user->ParticlesLocal[ipart].y));
			user->ParticlesLocal[ipart].iz            = (PetscScalar) (zs+Bisection(zc,0,zm-2,user->ParticlesLocal[ipart].z));
			user->ParticlesLocal[ipart].cpu           = rank;

			user->ParticlesLocal[ipart].eta           = 0.0;
			user->ParticlesLocal[ipart].zetha         = 0.0;
			user->ParticlesLocal[ipart].phi           = 0.0;

			user->ParticlesLocal[ipart].Plastic       = 0.0;
			user->ParticlesLocal[ipart].Mu_eff        = 0.0;
			user->ParticlesLocal[ipart].P			  = 0.0;
			user->ParticlesLocal[ipart].PlasticStrain = 0.0;
			user->ParticlesLocal[ipart].Strain	      = 0.0;
			user->ParticlesLocal[ipart].Txx		      = 0.0;
			user->ParticlesLocal[ipart].Tyy		      = 0.0;
			user->ParticlesLocal[ipart].Tzz		      = 0.0;
			user->ParticlesLocal[ipart].Txy		      = 0.0;
			user->ParticlesLocal[ipart].Txz		      = 0.0;
			user->ParticlesLocal[ipart].Tyz		      = 0.0;
			user->ParticlesLocal[ipart].E2nd	      = 0.0;
			user->ParticlesLocal[ipart].T2nd	      = 0.0;
			user->ParticlesLocal[ipart].E2nd_pl		  = 0.0;
			user->ParticlesLocal[ipart].Mu_viscous	  = 0.0;

			if(num_props==6)
			{
				user->ParticlesLocal[ipart].T         = larray_prtcls[ipart*num_props+5];
			}
			else
			{
				user->ParticlesLocal[ipart].T         = 0.0;
			}

		}
		else {

			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP,"Oops; the file we are loading has neither 5 nor 28 properties on the particles, so I don't know how to deal with it");

		}





	}

	if (user->LoadInitialParticlesFromDisc==1){
		// Free specific allocated memory
		ierr = VecRestoreArray(lvec_prtcls,&larray_prtcls);            CHKERRQ(ierr);
		ierr = VecDestroy(&lvec_prtcls); 		                      CHKERRQ(ierr);
		ierr = VecDestroy(&lvec_info); 		                          CHKERRQ(ierr);

		PetscPrintf(PETSC_COMM_WORLD,"[%lld] Loaded initial particle distribution (%lld prtcls) from directory ./%s on %lld procs \n",(LLD)rank,(LLD)user->num_particle_local,user->LoadInitialParticlesDirectory, (LLD)num_procs );
	}

	if (user->LoadInitialParticlesFromDisc==2){
		PetscPrintf(PETSC_COMM_WORLD,"[%lld] Loaded particle distribution (%lld prtcls) from directory ./MatlabInputParticles on %lld procs \n",(LLD)rank,(LLD)user->num_particle_local,(LLD)num_procs);

		/* Check for AddRandomNoiseParticles */
		AddRandomNoiseParticles = 0;
		ierr = PetscOptionsGetInt(PETSC_NULL,"-AddRandomNoiseParticles",&AddRandomNoiseParticles	, PETSC_NULL); CHKERRQ(ierr);

		dx = user->W/((double) ((nx-1)*user->NumPartX));
		dy = user->L/((double) ((ny-1)*user->NumPartY));
		dz = user->H/((double) ((nz-1)*user->NumPartZ));

		// Add random noise if required
		if (AddRandomNoiseParticles==0){
			PetscPrintf( PETSC_COMM_WORLD, "# Not adding random noise to particle distribution... \n");
		}
		if (AddRandomNoiseParticles==1	){
			PetscPrintf( PETSC_COMM_WORLD, "# Adding random noise to particle distribution... \n");

			// Initialize the Random number generator
			ierr = PetscRandomCreate(PETSC_COMM_WORLD, &rctx); CHKERRQ(ierr);
			ierr = PetscRandomSetFromOptions(rctx); CHKERRQ(ierr);
			//PetscRandomSetInterval(rctx,0.0,0.1); //option to control/reduce the noise in the particles default =(0,1)

			// Add some random component to the particle
			for (ipart=0; ipart<user->num_particle_local; ipart++){
				ierr 	= PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
				cf_rand = cf_rand-0.5;
				user->ParticlesLocal[ipart].x =	user->ParticlesLocal[ipart].x + cf_rand*dx;
				ierr 	= PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
				cf_rand = cf_rand-0.5;
				user->ParticlesLocal[ipart].y =	user->ParticlesLocal[ipart].y + cf_rand*dy;
				ierr 	= PetscRandomGetValueReal(rctx, &cf_rand); CHKERRQ(ierr);
				cf_rand = cf_rand-0.5;
				user->ParticlesLocal[ipart].z =	user->ParticlesLocal[ipart].z + cf_rand*dz;
			}

			ierr = PetscRandomDestroy(&rctx); 							  CHKERRQ(ierr);	// Destroy random context
		}

		// Free specific allocated memory
		ierr = PetscFree(larray_info); 								  CHKERRQ(ierr);
		ierr = PetscFree(larray_prtcls); 							  CHKERRQ(ierr);

	}

	// Free allocated memory
	ierr = PetscFree(xc);                                         CHKERRQ(ierr);
	ierr = PetscFree(yc);                                         CHKERRQ(ierr);
	ierr = PetscFree(zc);                                         CHKERRQ(ierr);

	MPI_Barrier(PETSC_COMM_WORLD);
	PetscFunctionReturn(0);
}



/* Load initial Particles from disc =================================================================== */
#undef __FUNCT__
#define __FUNCT__ "LoadInitialParticlesFromDisc"
PetscErrorCode LoadInitialParticlesFromDisc( UserContext *user)
{
	PetscMPIInt			rank, size;
	PetscErrorCode 		ierr;
	PetscInt			num_procs, numParticles, ipart,ix[1];
	PetscScalar 		num[1];
	char				LoadFileNameParticles[PETSC_MAX_PATH_LEN];
	PetscViewer			view_in;
	Vec 				Particle_Vec, information;
	Particles			*Particle_array;

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&size);


	sprintf(LoadFileNameParticles  ,"%s.%lld.out",	"./InitialParticles/Initial_Particles"				, (LLD)rank	); 	// construct the filename with Particles


	/* load particles */
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,LoadFileNameParticles,FILE_MODE_READ, &view_in); CHKERRQ(ierr);
	ierr = VecCreate(PETSC_COMM_SELF,&information);		CHKERRQ(ierr);
	ierr = VecCreate(PETSC_COMM_SELF,&Particle_Vec);	CHKERRQ(ierr);
	ierr = VecLoad(information, view_in);				CHKERRQ(ierr);
	ierr = VecLoad(Particle_Vec, view_in);				CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&view_in);				CHKERRQ(ierr);


	/* Error checking */
	ix[0] = 0; VecGetValues(information,1,ix,num);	num_procs = (PetscInt)(num[0]);
	if (size != num_procs){
		SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_SUP,"InitialParticles were created on %lld procs; yet you currently run with %lld procs!",(LLD)num_procs, (LLD)size);
	}



	/* The number of particles loaded on the current PROC is equal to the length of the vector */
	ierr = VecGetSize(Particle_Vec, &numParticles);				CHKERRQ(ierr);
	user->num_particle_local = numParticles/particle_props;
	PetscPrintf(PETSC_COMM_SELF,"[%lld] Loaded initial particle distribution (%lld prtcls) from directory ./InitialParticles on %lld procs \n",(LLD)rank,(LLD)user->num_particle_local,(LLD)num_procs);


	// Put particles at the correct location
	ierr = VecGetArray(Particle_Vec,(PetscScalar**)&Particle_array); CHKERRQ(ierr);
	for (ipart=0; ipart<user->num_particle_local; ipart++){

		//PetscPrintf(PETSC_COMM_S,"Loaded particles number %i \n",(LLD)ipart);
		user->ParticlesLocal[ipart] = Particle_array[ipart];

	}
	ierr = VecRestoreArray(Particle_Vec,(PetscScalar**)&Particle_array); CHKERRQ(ierr);

	ierr = VecDestroy(&Particle_Vec); 		CHKERRQ(ierr);
	ierr = VecDestroy(&information); 		CHKERRQ(ierr);


	PetscPrintf(PETSC_COMM_WORLD,"Loaded initial particle distribution from directory ./InitialParticles on %lld procs \n",(LLD)num_procs);

	PetscFunctionReturn(0);

}
/* ==================================================================================================== */


/* Write particles to disc - Could be improved memory-wise! =========================================== */
#undef __FUNCT__
#define __FUNCT__ "WriteParticlesToDisc"
PetscErrorCode WriteParticlesToDisc( UserContext *user, PetscInt itime )
{
	PetscErrorCode ierr;
	Vec				Particle_Vec, information;
	PetscMPIInt		rank, size;
	PetscViewer		view_out;
	char			SaveFileName[PETSC_MAX_PATH_LEN];

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&size);


	/* Create a vector and store particle in this vector - this is potentially memory-wasteful  */

	// What is done here is terrible in terms of common sense and language standards. Durty ad-hoc programming.
	// We are lucky that particles consist of "PetscScalar"s only
	// This should be replaced

	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,
			1,
			user->num_particle_local*particle_props,
			(PetscScalar*)(user->ParticlesLocal),
			&Particle_Vec); CHKERRQ(ierr);


	ierr = VecAssemblyBegin(Particle_Vec); 	VecAssemblyEnd(Particle_Vec); CHKERRQ(ierr);

	/* Create a vector with information that's helpfull in matlab */
	ierr = VecCreateSeq(PETSC_COMM_SELF,10, &information); CHKERRQ(ierr);
	ierr = VecSetValue(information,0, ((PetscScalar)(user->num_particle_local)), INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,1, ((PetscScalar)(size)), 					INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(information,2, ((PetscScalar)(particle_props)), 			INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(information); 	VecAssemblyEnd(information); CHKERRQ(ierr);


	// Write info to disc --------------------------------------------------------
	sprintf(SaveFileName,"Particles.%lld.%lld.out",(LLD)rank,(LLD)itime+1000000LL);  	 	// Construct the filename
//	ierr = PetscViewerBinaryMatlabOpen(PETSC_COMM_SELF,SaveFileName,&view_out); CHKERRQ(ierr);	// Open file
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,SaveFileName,FILE_MODE_WRITE,&view_out); CHKERRQ(ierr);	// Open file
	ierr = PetscViewerSetFormat(view_out, PETSC_VIEWER_BINARY_MATLAB); CHKERRQ(ierr);

	ierr = VecView(information, 				view_out); CHKERRQ(ierr);
	ierr = VecView(Particle_Vec, 				view_out); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&view_out); CHKERRQ(ierr);											// Close file
	// ----------------------------------------------------------------------------

	PetscPrintf(PETSC_COMM_WORLD,"# Saved Particles to file %s \n",SaveFileName);

	ierr = VecDestroy(&Particle_Vec); CHKERRQ(ierr);
	ierr = VecDestroy(&information); CHKERRQ(ierr);
	//	}

	PetscFunctionReturn(0);
}
/* ==================================================================================================== */

#endif


/* ====================================================================================================
 * Perform phase transitions on particles 					*/
#undef __FUNCT__
#define __FUNCT__ "ParticlePhaseTransitions"
PetscErrorCode ParticlePhaseTransitions( UserContext *user )
{
	PetscInt		ipart,	phase, PhaseTransformed, iphase, Num_PT;
	Particles		ParticleLocal;

	Num_PT = 0;
	if (user->num_phase_transitions>0){

		/* Loop over all particles */
		for (ipart=0; ipart<user->num_particle_local; ipart++){
			ParticleLocal = user->ParticlesLocal[ipart];
			phase 		  = ((PetscInt) ParticleLocal.phase);

			PhaseTransformed		  = 0;
			for (iphase=0; iphase<user->num_phase_transitions; iphase++){

				if (phase == user->PhaseTransitions[iphase].InitialPhase){
					/* Particle could be affected */

					if (user->PhaseTransitions[iphase].TransitionType==1){
						/* Depth-based transition */
						if (user->PhaseTransitions[iphase].TransitionBelow==1){
							if (ParticleLocal.z < user->PhaseTransitions[iphase].TransitionDepth){
								phase = user->PhaseTransitions[iphase].TransformedPhase;
								PhaseTransformed = 1;
							}
						} else {
							if (ParticleLocal.z > user->PhaseTransitions[iphase].TransitionDepth ){
								phase = user->PhaseTransitions[iphase].TransformedPhase;
								PhaseTransformed = 1;
							}
						}
					}
					else if (user->PhaseTransitions[iphase].TransitionType==2){
						SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER," This type of phase transition has not yet been implemented! \n");

						//	ierr = MPI_Abort(PETSC_COMM_WORLD,1); CHKERRQ(ierr);
					}

				}

			}


			if (PhaseTransformed==1){
				ParticleLocal.phase 		= ((double) phase);
				user->ParticlesLocal[ipart]	=	ParticleLocal;
				Num_PT 						= Num_PT+1;
			}
		}
		PetscPrintf(PETSC_COMM_WORLD," # of phase transitions = %lld \n", (LLD)Num_PT);
	}

	PetscFunctionReturn(0);
}
/* ==================================================================================================== */

/* ====================================================================================================
 *  Correct the phase of particles depending on whether they are above or below the internal free surface.
 *
 *  If they are above, rock particles are transformed into air.
 *  If they are below, air particles are transformed into rock.
 *
 */
#undef __FUNCT__
#define __FUNCT__ "CorrectPhasesForInternalFreeSurface"
PetscErrorCode CorrectPhasesForInternalFreeSurface( UserContext *user )
{
	PetscErrorCode 	ierr;
	DM				cda_SurfaceTopo;
	Vec 			gc_SurfaceTopo, LocalSurfaceTopography_vec;
	PetscMPIInt     rank;
	PetscBool		InjectSpecificPhase=PETSC_FALSE,flg;
	DMDACoor3d		***coors_SurfaceTopo;
	PetscInt		PhaseAir, ipart, Num_PT, phase;
	PetscInt		xs,ys,zs,xm,ym;
	PetscScalar		***LocalSurfaceTopography, minSurfaceTopo, maxSurfaceTopo;
	Particles		ParticleLocal;
	PetscLogDouble	cputime_start, cputime_end;

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	if (user->ErosionParameters.UseInternalFreeSurface != 0){
		PetscPrintf(PETSC_COMM_WORLD,"# Correcting particle phase due to internal free surface motion...  ");

		PetscTime(&cputime_start);

		/* Retrieve internal free surface */
		ierr = DMGetCoordinateDM(user->DA_SurfaceTopography,						&cda_SurfaceTopo		); 	CHKERRQ(ierr);
		ierr = DMGetCoordinatesLocal(user->DA_SurfaceTopography,					&gc_SurfaceTopo			); 	CHKERRQ(ierr);
		ierr = DMDAVecGetArray(cda_SurfaceTopo,gc_SurfaceTopo,						&coors_SurfaceTopo		); 	CHKERRQ(ierr);

		ierr = DMGetLocalVector(user->DA_SurfaceTopography,	&LocalSurfaceTopography_vec); CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin(user->DA_SurfaceTopography, user->SurfaceTopography, INSERT_VALUES, LocalSurfaceTopography_vec); CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd(user->DA_SurfaceTopography,   user->SurfaceTopography, INSERT_VALUES, LocalSurfaceTopography_vec); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(user->DA_SurfaceTopography,LocalSurfaceTopography_vec, 	&LocalSurfaceTopography	);	CHKERRQ(ierr);

		/* Minimum & maximum of surface topograpy [to speed up code later]			*/
		VecMin(user->SurfaceTopography,		PETSC_NULL, 	&minSurfaceTopo);
		VecMax(user->SurfaceTopography,  	PETSC_NULL, 	&maxSurfaceTopo);

		PhaseAir 	=	user->ErosionParameters.StickyAirPhase;

		ierr = PetscOptionsGetBool( PETSC_NULL, "-InjectSpecificPhaseBelowFreeSurface", &InjectSpecificPhase, &flg );	// do we want to inject a specific phase?

		/* Size of free surface grid */
		ierr = DMDAGetCorners(user->DA_SurfaceTopography,&xs,&ys,&zs,&xm,&ym,PETSC_NULL); CHKERRQ(ierr);

		/* Loop over all particles */
		Num_PT = 0;
		for (ipart=0; ipart<user->num_particle_local; ipart++){
			PetscScalar		part_z, dx, dy, eta, xsi, N[4], z_FreeSurface;
			PetscInt 		ix,iy,iz;

			ParticleLocal	= user->ParticlesLocal[ipart];
			phase 		  	= ((PetscInt) ParticleLocal.phase);

			part_z 			=	ParticleLocal.z;

			if 	((phase>=0)	&& ((( (phase==PhaseAir) && (part_z<maxSurfaceTopo)) ||
					((phase!=PhaseAir) && (part_z>minSurfaceTopo)) ))) {

				/* Air particle that might be below the internal free surface, or a rock particle that is now in the air
				 *
				 * we have to evaluate these cases more carefully by interpolating the precise free surface height above the particle */

				/* In which element is this particle? */
				ix 				=	(PetscInt) ParticleLocal.ix;
				iy 				=	(PetscInt) ParticleLocal.iy;
				iz 				=	(PetscInt) ParticleLocal.iz;

				if ((ix<xs) || (ix>=(xs+xm)) || (iy<ys) || (iy>=(ys+ym))){
					SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER," Particle seems to be on a the wrong processor. \n");
				}

				/* Interpolate free surface height using linear shape */
				dx 				=	coors_SurfaceTopo[zs][iy  ][ix+1].x-coors_SurfaceTopo[zs][iy  ][ix  ].x;
				dy 				=	coors_SurfaceTopo[zs][iy+1][ix  ].y-coors_SurfaceTopo[zs][iy  ][ix  ].y;

				/* natural coordinates [straightforward as it is an undeformed element */
				eta 			=	(ParticleLocal.x - coors_SurfaceTopo[zs][iy][ix].x)/dx*2  - 1.0;
				xsi 			=	(ParticleLocal.y - coors_SurfaceTopo[zs][iy][ix].y)/dy*2  - 1.0;

				/* Shape function	*/
				N[0] 			=	0.25*(1-eta)*(1-xsi);
				N[1] 			=	0.25*(1+eta)*(1-xsi);
				N[2] 			=	0.25*(1+eta)*(1+xsi);
				N[3] 			=	0.25*(1-eta)*(1+xsi);

				/* Interpolate height of the free surface at the location of the particle*/
				z_FreeSurface 	=	N[0]*LocalSurfaceTopography[zs][iy  ][ix  ] + N[1]*LocalSurfaceTopography[zs][iy  ][ix+1] +
						N[2]*LocalSurfaceTopography[zs][iy+1][ix+1] + N[3]*LocalSurfaceTopography[zs][iy+1][ix  ];



				/* Rock particle that is now in the air -> transform into air */
				if ((phase!=PhaseAir) && (part_z>z_FreeSurface)){
					ParticleLocal.phase 		= 	((double) PhaseAir);
					Num_PT 						= 	Num_PT+1;

				}


				/* Air particle is below the free surface
				 *
				 * There are a number of ways to deal with this.
				 * 1) Use the dominant, non-air phase	[default]
				 * 2) Inject a specific phase.	[select -InjectSpecificPhaseBelowFreeSurface on command line. The phase is set with -ParticleInjectionPhase]
				 *
				 * 3) If we use	a sedimentation algorithm, select a specific phase that might change with time.
				 *    This is not implemented here, but rather in the routine ApplySedimentationAndFastErosion
				 *
				 * 4) Use the dominant phase at the TOP of the cell without air (or without other selected phases; select those with e.g., -PhasesToBeExcludedInDominantPhaseCalculation 1,2,3 from the command-line))
				 *	 [select -UseDominantNonAirPhase_TopCell at the command line]
				 *
				 */
				if ((phase==PhaseAir) &&  (part_z<z_FreeSurface)){
					PetscInt 	DominantPhase_WithAir, DominantPhase_WithoutAir, DominantPhase_WithoutAir_TopCell;
					PetscBool 	found;



					DominantPhase_WithoutAir = user->ParticleInjectionPhase;
					/* compute dominant, non-air phase */
					ierr =  FDSTAG_DominantPhaseAtCorners(user, ix, iy,iz, &DominantPhase_WithAir, &DominantPhase_WithoutAir, &DominantPhase_WithoutAir_TopCell); CHKERRQ(ierr);
					//PetscPrintf(PETSC_COMM_SELF,"Air->Rock; DominantPhase_WithAir=%i, DominantPhase_WithoutAir=%i DominantPhase_WithoutAir_TopCell=%i\n",DominantPhase_WithAir,DominantPhase_WithoutAir, DominantPhase_WithoutAir_TopCell);

					if (DominantPhase_WithoutAir<0){	// no dominant phase found, use default
						DominantPhase_WithoutAir = user->ParticleInjectionPhase;
					}
					if (DominantPhase_WithoutAir_TopCell<0){
						DominantPhase_WithoutAir_TopCell = user->ParticleInjectionPhase;
					}

					/* In case we want to inject a very specific phase. Add command-line options
					 * -InjectSpecificPhaseBelowFreeSurface	-ParticleInjectionPhase PhaseNumberYouWantToInject
					 */
					if (InjectSpecificPhase){
						DominantPhase_WithoutAir = user->ParticleInjectionPhase;
					}
					if  (PetscAbs(user->ErosionParameters.SedimentationRate)>0.0){
						// the sedimentation algorithm is active and tells us which sediment phase should be added
						DominantPhase_WithoutAir = user->ErosionParameters.PhaseSedimented;
					}


					PetscOptionsName("-UseDominantNonAirPhase_TopCell",PETSC_NULL,PETSC_NULL,&found);
					if (found){
						ParticleLocal.phase 		= 	((double) DominantPhase_WithoutAir_TopCell);
					}
					else{
						ParticleLocal.phase 		= 	((double) DominantPhase_WithoutAir);
					}
					Num_PT 						= 	Num_PT+1;
				}

				user->ParticlesLocal[ipart]	=	ParticleLocal;

			}



		}	// end of loop over particles
		PetscTime(&cputime_end);

		PetscPrintf(PETSC_COMM_WORLD,"Corrected %i Particles... [%f s]\n", Num_PT, cputime_end-cputime_start);

		/* Cleaning up */
		ierr 	= 	DMDAVecRestoreArray(cda_SurfaceTopo,gc_SurfaceTopo,	&coors_SurfaceTopo); 							CHKERRQ(ierr);
		//ierr 	= 	DMDAVecRestoreArray(user->DA_SurfaceTopography,user->SurfaceTopography, 	&LocalSurfaceTopography	);	CHKERRQ(ierr);
		ierr 	= 	DMDAVecRestoreArray(user->DA_SurfaceTopography,LocalSurfaceTopography_vec, 	&LocalSurfaceTopography	);	CHKERRQ(ierr);
		ierr 	= 	DMRestoreLocalVector(user->DA_SurfaceTopography,	&LocalSurfaceTopography_vec); CHKERRQ(ierr);

	}


	PetscFunctionReturn(0);
}




/*==========================================================================================================*/
/* Computes properties such as pressure, strainrate, deviatoric stress @ integration points
 *
 * Note that temperature is added to the particles @ a different location.
 *
 */
#if 1
#undef __FUNCT__
#define __FUNCT__ "ComputePropertiesAtParticles"
PetscErrorCode ComputePropertiesAtParticles(LaMEMVelPressureDA C, DM da, Vec Velocity, DM da_pres, Vec Pressure,
		DM da_temp, Vec Temp, Vec Temp_old,  UserContext *user, PetscScalar dt )
{
	PetscMPIInt             rank;
	PetscErrorCode 			ierr;
	PetscInt       			ipart, mx,my,mz, ipres, iel_x, iel_y, iel_z, ii, i,j,k;
	PetscInt		 		xmp,ymp,zmp,xsp,ysp,zsp, ComputeFull;
	PetscInt				xm,ym,zm,xs,ys,zs;
	DAVPElementType 		element_type;
	PetscInt				npres;
	PetscInt		 		phase;
	DMDACoor3d		 		***coords, coord_elem[MAX_nnel];
	DM				 		cda;
	Vec			 			gc;
	Field			 		***velocity;
	PetscScalar	 			mu, mu_plastic, mu_viscous, P_lithos, sizeparticles, TotalMemoryOfParticles, MaxMemoryOfParticles, MinMemoryOfParticles;
	PetscScalar				V_element[MAX_edof], P_element[MAX_npres], Point[3], DevStress[6];
	PetscScalar 			***temperature, ***temperature_old, T_element[MAX_edof_temp],  T_old_element[MAX_edof_temp];
	Particles				ParticleLocal;
	Vec						local_Vel, Temp_local, Temp_old_local, Pressure_local, Pressure_local_duplicate;
	MaterialsElement		***materials;
	PointWiseInformation   	PointInformation;
	PetscInt				mod;
	P_array 				***PressureContinuous_array;
	PressureElemDynamic 	PressureNew_data;
	PetscScalar 			***PressureNew_array;
	PressureElem	 		***PressureNew;


	PetscFunctionBegin;


	element_type = C->type;
	npres    = C->npres;

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	//coordinates
	ierr = DMGetCoordinateDM(da,&cda); CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(da,&gc); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(cda,gc,&coords); CHKERRQ(ierr);



	/* Copy pressure to local array, including ghostpoints */
	ierr = DMCreateLocalVector(da_pres,&Pressure_local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da_pres, Pressure, INSERT_VALUES, Pressure_local); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(da_pres,   Pressure, INSERT_VALUES, Pressure_local); CHKERRQ(ierr);

	ierr = VecDuplicate(Pressure_local, &Pressure_local_duplicate); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da_pres, 			Pressure_local, 		&PressureNew); CHKERRQ(ierr);
	if ((element_type==DAVP_Q1Q1) || (element_type == DAVP_FDSTAG)){
		ierr = DMDAVecGetArray(da_pres, 		Pressure_local_duplicate, 		&PressureContinuous_array); CHKERRQ(ierr);
	}
	else {
		ierr = DMDAVecGetArray(da_pres, 		Pressure_local_duplicate, 		&PressureNew_array); CHKERRQ(ierr);
	}


	/* Copy global velocity solution to local processor, including ghostpoints */
	ierr = DMGetLocalVector(da,&local_Vel); CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da, Velocity, INSERT_VALUES, local_Vel); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(da,   Velocity, INSERT_VALUES, local_Vel); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, local_Vel,	&velocity);		CHKERRQ(ierr);

	/* Copy temperature to local processor */
	if (da_temp != PETSC_NULL){
		ierr = DMCreateLocalVector(da_temp,&Temp_local); CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin(da_temp, Temp, INSERT_VALUES, Temp_local); CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd(da_temp,   Temp, INSERT_VALUES, Temp_local); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(da_temp, Temp_local,	&temperature); CHKERRQ(ierr);

		ierr = DMCreateLocalVector(da_temp,&Temp_old_local); CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin(da_temp, Temp_old, INSERT_VALUES, Temp_old_local); CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd(da_temp,   Temp_old, INSERT_VALUES, Temp_old_local); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(da_temp, Temp_old_local,	&temperature_old); CHKERRQ(ierr);
	}

	// get material properties
	ierr = DMDAVecGetArray(user->DA_Materials, user->Materials, &materials); CHKERRQ(ierr);

	ierr = DMDAGetInfo(user->DA_Processors,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0);	CHKERRQ(ierr);  	// # of nodes in all directions
	ierr = DMDAGetCorners(user->DA_Processors,&xsp,&ysp,&zsp,&xmp,&ymp,&zmp);CHKERRQ(ierr);	// the local node portion

	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);	//


	/* Loop over all particles */
	iel_x =0; iel_y=0; iel_z=0;
	for (ipart=0; ipart<user->num_particle_local; ipart++){
		ParticleLocal = user->ParticlesLocal[ipart];

		/* Compute the element it belongs to  */
		ParticleLocal	=	user->ParticlesLocal[ipart];
		if (ParticleLocal.phase>-4){		// only if particle is active

			/* Detect in which element the particle is in */
			if(  (__ELEMENT_TYPE__ == ELEMENT_Q1P0) ||  (element_type==DAVP_Q1Q1) || (element_type == DAVP_FDSTAG) ) {
				iel_x	=	((PetscInt) ParticleLocal.ix);
				iel_y	=	((PetscInt) ParticleLocal.iy);
				iel_z	=	((PetscInt) ParticleLocal.iz);
				i 		= 	iel_x;
				j 		= 	iel_y;
				k 		= 	iel_z;

			}
			else if( __ELEMENT_TYPE__ == ELEMENT_Q2P1 ){
				iel_x 	= ((PetscInt) ParticleLocal.ix); LaMEMMod(iel_x, 2, &mod); if (mod>0){iel_x = iel_x-1;}; iel_x = iel_x/2;
				iel_y 	= ((PetscInt) ParticleLocal.iy); LaMEMMod(iel_y, 2, &mod); if (mod>0){iel_y = iel_y-1;}; iel_y = iel_y/2;
				iel_z 	= ((PetscInt) ParticleLocal.iz); LaMEMMod(iel_z, 2, &mod); if (mod>0){iel_z = iel_z-1;}; iel_z = iel_z/2;
				i 		= 2*iel_x;
				j 		= 2*iel_y;
				k 		= 2*iel_z;
			}
			else {
				SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Unknown element type" );
			}


			// debugging !what a crap, not going to update this!
			if (rank==1){
				//	PetscPrintf(PETSC_COMM_SELF,"Particle %lld on rank %lld has iel=[%lld,%lld,%lld] and [ix,iy,iz]=[%lld,%lld,%lld]   but the local velocity grid goes from x=[%lld-%lld] y=[%lld-%lld] z=[%lld-%lld]\n",ipart, rank, iel_x, iel_y, iel_z, ((PetscInt) ParticleLocal.ix), ((PetscInt) ParticleLocal.iy), ((PetscInt) ParticleLocal.iz), xsp,xsp+xmp-1,ysp,ysp+ymp-1,zsp,zsp+zmp-1);
			}


			/* Extract coordinates of the local element in correct order */
			ierr = GetElementCoords(coord_elem, coords, i,j,k,1);CHKERRQ(ierr);
			ierr = CorrectElementCoordsForPeriodicity( coord_elem, i, j, user, 1);CHKERRQ(ierr);


			/* Extract Velocity of the current element */
			ierr = GetVelocityElement(velocity, V_element, i,j,k); CHKERRQ(ierr);

			Point[0] 	= ParticleLocal.eta;
			Point[1]	= ParticleLocal.zetha;
			Point[2]	= ParticleLocal.phi;

			/* Extract Pressure shape function of the current element */
			if ((element_type==DAVP_Q1Q1 ) || (element_type == DAVP_FDSTAG)){
				/* Continuous pressure */
				ierr = GetPressureElement(PressureContinuous_array, P_element, iel_x,iel_y,iel_z); CHKERRQ(ierr);
			}
			else {
				/* Discontinuous pressure */
				LaMEMSetPressueElemDataMemoryFromArray( &PressureNew_data, iel_x,iel_y,iel_z, npres, PressureNew_array );
				for (ipres=0; ipres<npres; ipres++){
					P_element[ipres]  = PressureNew_data.P[ipres];
				}
			}


			/* Extract Temperature nodal values of the current element */
			if (da_temp != PETSC_NULL){
				ierr = GetTemperatureElement(temperature    , T_element	, i,j,k);CHKERRQ(ierr);
				ierr = GetTemperatureElement(temperature_old, T_old_element, i,j,k);CHKERRQ(ierr);

			}
			else {
				for (ii=0; ii<MAX_edof_temp; ii++){	T_element[ii] = T_old_element[ii] = 0.0;	}
			}

			/* Extract old stress from tracer */
			//==========================================================
			DevStress[0] = 	ParticleLocal.Txx;		DevStress[1] = 	ParticleLocal.Tyy;
			DevStress[2] = 	ParticleLocal.Tzz; 		DevStress[3] = 	ParticleLocal.Txy;
			DevStress[4] = 	ParticleLocal.Txz;		DevStress[5] = 	ParticleLocal.Tyz;
			//==========================================================


			/* Extract material parameters from tracer */
			phase 	= 	((PetscInt) ParticleLocal.phase);							/* Phase */
			//			G		=	user->PhaseProperties.ElasticShearModule[phase];		/* Elastic shear module */
			P_lithos=	 (user->z_bot + user->H - ParticleLocal.z)*user->PhaseProperties.rho[phase]*user->Gravity;

			mu 			= 	ParticleLocal.Mu_eff;									/* Effective viscosity of LAST timestep */
			mu_viscous  =	ParticleLocal.Mu_viscous;								/* Real viscosity of LAST timestep */
			if (PetscAbsScalar(mu)<1e-3){
				mu		=	user->PhaseProperties.mu[phase];
				mu_viscous = mu;
			}



			/* Compute properties at the tracer location */
			ComputeFull = 1;
			ierr = ComputePointWiseProperties( C, &PointInformation, Point, V_element, P_element, T_element, T_old_element,
					DevStress, coord_elem, mu, mu_viscous, ComputeFull ); CHKERRQ(ierr);

			// if the particle failed plastically before, compute the plastic viscosity based on last stresses and strain rates
			if ( ( (PetscInt) ParticleLocal.Plastic)==1 ){
				mu_plastic = ParticleLocal.T2nd/(2.0*ParticleLocal.E2nd);
			}
			else {
				mu_plastic = 0.0;
			}

			// Update the effective viscosity of the particle
			ierr = ComputeEffectiveViscosity(phase, user, PointInformation.Temperature, PointInformation.SecondInvariantDeviatoricStrainrate,
					PointInformation.SecondInvariantDeviatoricStress,  PointInformation.Pressure,
					P_lithos, ParticleLocal.PlasticStrain, 0,
					&mu_viscous, &mu_plastic, &mu); CHKERRQ(ierr);

			// if it was NOT plastic before, but it is now, stresses should be re-estimated
			if ( ( ( (PetscInt) ParticleLocal.Plastic)==0 ) && (PetscAbsScalar(mu_plastic)>0.0) ){
				ierr =	ComputePointWiseProperties( C, &PointInformation, Point, V_element, P_element, T_element, T_old_element,
						DevStress, coord_elem, mu, mu_viscous, ComputeFull ); CHKERRQ(ierr);
			}

			ParticleLocal.Mu_eff  		= mu;
			ParticleLocal.Mu_viscous 	= mu_viscous;
			if (PetscAbsScalar(mu_plastic)>0){
				ParticleLocal.Plastic = 1;
			}
			else{
				ParticleLocal.Plastic = 0;
			}


			//==========================================================
			/* Store information on particle */
			ParticleLocal.Txx  = PointInformation.DeviatoricStress[0];
			ParticleLocal.Tyy  = PointInformation.DeviatoricStress[1];
			ParticleLocal.Tzz  = PointInformation.DeviatoricStress[2];
			ParticleLocal.Txy  = PointInformation.DeviatoricStress[3];
			ParticleLocal.Txz  = PointInformation.DeviatoricStress[4];
			ParticleLocal.Tyz  = PointInformation.DeviatoricStress[5];

			//==========================================================

			ParticleLocal.T2nd 			= 	PointInformation.SecondInvariantDeviatoricStress;
			ParticleLocal.E2nd 			= 	PointInformation.SecondInvariantDeviatoricStrainrate;
			ParticleLocal.E2nd_pl 		= 	PointInformation.PlasticStrainrate2ndInvariant;

			// Temperature
			if (!user->ArtTemp)
			{
				ParticleLocal.T 			= 	PointInformation.Temperature;		// update complete temperature
			}
			//ParticleLocal.T 			= 	ParticleLocal.T + PointInformation.Temperature_diff;		// update temperature difference

			ParticleLocal.P    			= 	PointInformation.Pressure;
			ParticleLocal.Strain 		=	ParticleLocal.Strain		+ 	dt*ParticleLocal.E2nd;
			ParticleLocal.PlasticStrain	=	ParticleLocal.PlasticStrain	+ 	dt*ParticleLocal.E2nd_pl;


			/* Put back particle */
			user->ParticlesLocal[ipart] = 	ParticleLocal;

		} // end of if phase>-4
	} // end of loop over particles



	/* restore materials */
	ierr = DMDAVecRestoreArray(user->DA_Materials,   user->Materials, &materials); CHKERRQ(ierr);

	/* Restore temperature */
	if (da_temp != PETSC_NULL){
		ierr = DMDAVecRestoreArray(da_temp, Temp_local,		&temperature	); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(da_temp, Temp_old_local,	&temperature_old); CHKERRQ(ierr);

		ierr = VecDestroy(&Temp_local); CHKERRQ(ierr);
		ierr = VecDestroy(&Temp_old_local); CHKERRQ(ierr);
	}

	/* velocity */
	ierr = DMDAVecRestoreArray(da,local_Vel,&velocity);	CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(da,&local_Vel); CHKERRQ(ierr);

	/* restore pressure */
	ierr = DMDAVecRestoreArray(da_pres,			Pressure_local,			&PressureNew); CHKERRQ(ierr);
	if ((element_type==DAVP_Q1Q1) || (element_type == DAVP_FDSTAG)){
		ierr = DMDAVecRestoreArray(da_pres, 		Pressure_local_duplicate, 		&PressureContinuous_array); CHKERRQ(ierr);
	}
	else {
		ierr = DMDAVecRestoreArray(da_pres, 		Pressure_local_duplicate, 		&PressureNew_array); CHKERRQ(ierr);
	}
	ierr = VecDestroy(&Pressure_local); CHKERRQ(ierr);
	ierr = VecDestroy(&Pressure_local_duplicate); CHKERRQ(ierr);



	ierr = DMDAVecRestoreArray(cda,gc,&coords);	CHKERRQ(ierr);
	//	ierr = DMDestroy(cda); CHKERRQ(ierr);
	//	ierr = VecDestroy(gc);

	sizeparticles = (double)((size_t)user->MaxNumLocalParticles*sizeof(Particles))*1e-6;
	ierr = MPI_Allreduce(&sizeparticles, &TotalMemoryOfParticles, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"The total memory used by particles is %e MB\n",TotalMemoryOfParticles);
	ierr = MPI_Allreduce(&sizeparticles, &MaxMemoryOfParticles, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"The max memory used by particles is %e MB\n",MaxMemoryOfParticles);
	ierr = MPI_Allreduce(&sizeparticles, &MinMemoryOfParticles, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"The min memory used by particles is %e MB\n",MinMemoryOfParticles);
	//MPI_Abort(PETSC_COMM_WORLD,1);

	/* Clean up */
	PetscFunctionReturn(0);
}
/*==========================================================================================================*/
#endif



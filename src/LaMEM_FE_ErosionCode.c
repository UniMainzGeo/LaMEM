/*
 LaMEM_FE_ErosionCode
 
 Includes routines to model erosion using a finite difference technique, which can be coupled to LaMEM.
 
 
 
 Code developed by Marine Collignon and Boris Kaus
 
 $Id$
 
 */

#include "LaMEM.h"
#include "LaMEM_FE_ErosionCode.h"


/*==============================================================================================================================
 * Solve one tectonic timestep of the erosion code, which is typically split into a number of sub-timesteps
 *
 */
#undef __FUNCT__
#define __FUNCT__ "FE_ErosionCode_TectonicTimestep"
PetscErrorCode FE_ErosionCode_TectonicTimestep( UserContext *user)
{
	PetscErrorCode 			ierr;
	PetscInt				num_timesteps, itime, M, N;
	PetscScalar				dt_tectonics, dt_erosion, dt, time;
	FE_Erosion_UserData 	Erosion_UserData;
    PetscInt                *bcdof;
    PetscScalar             *bcval, SecYear;
    PetscScalar		        L2_Topo;
    PetscLogDouble     		cputime_start, cputime_end,cputime_start0, cputime_end0;
    
    SecYear = 3600*24*365.25;
    
    
    PetscTime(&cputime_start);
    
    PetscTime(&cputime_start0);
	/* Compute number of timesteps that should be performed with the erosion code */
	dt_tectonics 	= 	user->dt * user->Characteristic.Time;					// tectonic timestep in seconds
	dt_erosion 		=	user->ErosionParameters.FE_ErosionCode.dt;				// ideal erosion timestep
	num_timesteps 	=	((PetscInt) (dt_tectonics/dt_erosion)) + 1;				// how many timesteps do we take?
    
	dt 				=	dt_tectonics/((PetscScalar) num_timesteps);				// dt employed in Erosion code
    
	/* Initialise local vectors and DMDA's */
    
	/* LaMEM computes in nondimensional units - FE_ErosionCode in dimensional units
	 * Transfer the erosion grid from nondimensional to dimensional units.
	 * NOTE: this must be change in the future
	 */
	ierr 	=	TransferLaMEMDataToAndInitialize_FE_ErosionCode( user, &Erosion_UserData, dt);          CHKERRQ(ierr);
    PetscTime(&cputime_end0);
    PetscPrintf(PETSC_COMM_SELF,"  FE Erosion code: TransferLaMEMDataToAndInitialize_FE_ErosionCode took %f s \n",cputime_end0-cputime_start0);
    
    /* Define Boundary Conditions for the FE code: 4-options: BC=0, BC=1, BC=2, or BC=something else */
    PetscTime(&cputime_start0);
    
	if 	(Erosion_UserData.BC==0){           // left-right fixed, front/back zero-flux
		N = Erosion_UserData.N;
		ierr = PetscMalloc(sizeof(PetscInt)*2*(size_t)N, &bcdof);                                                   CHKERRQ(ierr);
		ierr = PetscMalloc(sizeof(PetscScalar)*2*(size_t)N, &bcval);                                                CHKERRQ(ierr);
	}
	else if 	(Erosion_UserData.BC==1){   // left-right zero-flux, front/back fixed
		M = Erosion_UserData.M;
		ierr = PetscMalloc(sizeof(PetscInt)*2*(size_t)M, &bcdof);                                                   CHKERRQ(ierr);
		ierr = PetscMalloc(sizeof(PetscScalar)*2*(size_t)M, &bcval);                                                CHKERRQ(ierr);
	}
    
    else if (Erosion_UserData.BC==2){       // fixed topography on the 4 boundaries
        N = Erosion_UserData.N;
        M = Erosion_UserData.M;
        ierr = PetscMalloc(sizeof(PetscInt)*2*(size_t)(M+N-2), &bcdof);                                                   CHKERRQ(ierr);
		ierr = PetscMalloc(sizeof(PetscScalar)*2*(size_t)(M+N-2), &bcval);                                                CHKERRQ(ierr);
    }
    
    else if (Erosion_UserData.BC==3){           // fixed one side bottom side correspond to front of orogen. no flux on the back boundary 
        M = Erosion_UserData.M; 
        ierr = PetscMalloc(sizeof(PetscInt)*(size_t)M, &bcdof);                                                   CHKERRQ(ierr);
		ierr = PetscMalloc(sizeof(PetscScalar)*(size_t)M, &bcval);                                                CHKERRQ(ierr);
        
    }
    
    
    if ((Erosion_UserData.BC==0) || (Erosion_UserData.BC==1) || (Erosion_UserData.BC==2) || (Erosion_UserData.BC==3)){
        ierr = DefineBC(&Erosion_UserData, bcdof, bcval);                                                   CHKERRQ(ierr);
    }
    
    PetscTime(&cputime_end0);
    PetscPrintf(PETSC_COMM_SELF,"  FE Erosion code: Initial part took %f s \n",cputime_end0-cputime_start0);
    
    
    /* time loop */
    time = 0.0;
    for (itime=0; itime<num_timesteps; itime++){
        
        /* call the function for streamflow and set coefficient for diffusion */
    	PetscTime(&cputime_start0);
        ierr = Set_Diffusion_Coefficient(&Erosion_UserData);                                            CHKERRQ(ierr);
        PetscTime(&cputime_end0);
        PetscPrintf(PETSC_COMM_SELF,"  FE Erosion code: Set_Diffusion_Coefficient took %f s \n",cputime_end0-cputime_start0);
        
        /* Non Linear Diffusion */
        PetscTime(&cputime_start0);
        ierr = Solve_NonlinearDiffusion(&Erosion_UserData, bcdof, bcval);                               CHKERRQ(ierr);
        PetscTime(&cputime_end0);
        PetscPrintf(PETSC_COMM_SELF,"  FE Erosion code: Solve nonlinear diffusion took %f s \n",cputime_end0-cputime_start0);
        
        /* Solve non linear diffusion - smoother */
        /*    ierr = Solve_linearDiffusion( &Erosion_UserData, bcdof, bcval);  CHKERRQ(ierr); */
        
		VecNorm(Erosion_UserData.Elevation,NORM_2,&L2_Topo);
		time = time + Erosion_UserData.dt;
        
		PetscPrintf(PETSC_COMM_SELF,"FE Erosion code: sub-timestep %i/%i, time = %1.2f [kyrs] || Elevation ||_2 = %f \n",itime+1, num_timesteps, time/SecYear/1e3, L2_Topo);
        
        
		ierr =  FE_Erosion_UpdateTopography_DA(&Erosion_UserData);                                      CHKERRQ(ierr);	// Update topography DA for visualization
        
        
#if 1
		{
			PetscViewer viewer;
			char 		name[PETSC_MAX_PATH_LEN];
			PetscBool 	flag=PETSC_FALSE;

        
        


			PetscOptionsHasName(PETSC_NULL,"-FE_ErosionCode_View",&flag);

        
			if (flag){

				sprintf(name,"ErosionCode_%d.vtk",itime);
				ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF, name, &viewer);                CHKERRQ(ierr);
				ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_VTK);                CHKERRQ(ierr);
				ierr = DMView(Erosion_UserData.da_elevation_output, viewer);                CHKERRQ(ierr);
				ierr = VecView(Erosion_UserData.Elevation, viewer);                         CHKERRQ(ierr);
				ierr = PetscViewerDestroy( &viewer );


				sprintf(name,"ErosionCode_Properties_%d.vtk",itime);
				ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF, name, &viewer);                CHKERRQ(ierr);
				ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_VTK);                CHKERRQ(ierr);
				ierr = DMView(Erosion_UserData.da_element, viewer);                CHKERRQ(ierr);
				ierr = VecView(Erosion_UserData.properties, viewer);                   CHKERRQ(ierr);
				ierr = PetscViewerDestroy( &viewer );

			}
		}

#endif
        
        
        
        
    }
    
	/* Transfer data back from dimensional to dimensionless units, for use in LaMEM */
    ierr = VecScale(Erosion_UserData.Elevation, 1/user->Characteristic.Length); 						CHKERRQ(ierr);
    ierr = VecCopy(Erosion_UserData.Elevation, user->ErosionParameters.FE_ErosionCode.ErosionSurface); 	CHKERRQ(ierr);
    
    /* Record time of Erosion solver */
    PetscTime(&cputime_end);
	PetscPrintf(PETSC_COMM_SELF,"Total time taken for the Erosion Solver %f \n",cputime_end-cputime_start);
	PetscPrintf(PETSC_COMM_SELF,"****************************************\n");
    
    
    /* Cleanup */
    ierr = PetscFree(bcdof);                                    CHKERRQ(ierr);
    ierr = PetscFree(bcval);                                    CHKERRQ(ierr);
    ierr = VecDestroy(&Erosion_UserData.Elevation);				CHKERRQ(ierr);
    ierr = DMDestroy(&Erosion_UserData.da_elevation);			CHKERRQ(ierr);
    ierr = DMDestroy(&Erosion_UserData.da_elevation_output);	CHKERRQ(ierr);
    ierr = DMDestroy(&Erosion_UserData.da_element);             CHKERRQ(ierr);
    ierr = VecDestroy(&Erosion_UserData.properties);			CHKERRQ(ierr);
    ierr = VecDestroy(&Erosion_UserData.Elevation_init);		CHKERRQ(ierr);
    
	PetscFunctionReturn(0);
    
}
/*==============================================================================================================================*/


/*==============================================================================================================================
 * Transfer the non-dimensional LaMEM data to the (dimensional) data used in FE_ErosionCode
 *
 * Obviously, this needs to be changed soon and the erosion code should also compute in dimensionless units
 *
 */
#undef __FUNCT__
#define __FUNCT__ "TransferLaMEMDataToAndInitialize_FE_ErosionCode"
PetscErrorCode TransferLaMEMDataToAndInitialize_FE_ErosionCode( UserContext *user, FE_Erosion_UserData *Erosion_UserData,
                                                               PetscScalar dt_erosion)
{
	PetscErrorCode  ierr;
	PetscInt		M,N, nelx, nely,i;
	PetscScalar 	day, year, rain_m_year, k0, c, hack, rain;
    PetscScalar     rain_river, rain_river_year;
	Vec 			gc1, gc;
	PetscScalar		xymin[2], xymax[2];
    PetscInt        BC, fill_lake, mode_river, nbre_river;
    PetscScalar     location_river[100];
    
	// copy elevation  data
	ierr 	= 	VecDuplicate(user->ErosionParameters.FE_ErosionCode.ErosionSurface, &Erosion_UserData->Elevation);              CHKERRQ(ierr);		//
    PetscObjectSetName((PetscObject)Erosion_UserData->Elevation,"Elevation");
    
	ierr 	= 	DMDAGetInfo(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode,PETSC_NULL,&M,&N,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);       CHKERRQ(ierr);
	ierr 	= 	DMDACreate2d(PETSC_COMM_SELF,DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,M,N,PETSC_DECIDE,PETSC_DECIDE,1,1,PETSC_NULL,PETSC_NULL,&Erosion_UserData->da_elevation);                            CHKERRQ(ierr);
    
    
	/* Within LaMEM, everything is in ND units - the erosion code currently works with dimensional units */
	ierr 	=	VecCopy(user->ErosionParameters.FE_ErosionCode.ErosionSurface,Erosion_UserData->Elevation);						CHKERRQ(ierr);		// copy coordinate values
	ierr 	= 	VecScale(Erosion_UserData->Elevation, user->Characteristic.Length);                                             CHKERRQ(ierr);		// in meters
    
	// copy coordinates from LaMEM erosion surface to the local erosion coordinates
	ierr 	= 	DMDASetUniformCoordinates(Erosion_UserData->da_elevation,0,1,0.0,1,0,0);                      CHKERRQ(ierr);		// set coordinates
	ierr	=	DMGetCoordinates(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode,&gc1);                              CHKERRQ(ierr);		// coordinates
	ierr	=	DMGetCoordinates(Erosion_UserData->da_elevation,&gc);                                                         CHKERRQ(ierr);		// coordinates
	ierr 	=	VecCopy(gc1, gc);                                                                                               CHKERRQ(ierr);		// copy coordinate values
	ierr 	= 	VecScale(gc, user->Characteristic.Length);                                                                      CHKERRQ(ierr);		// in meters
	ierr	=	DMSetCoordinates(Erosion_UserData->da_elevation,gc);                                                          CHKERRQ(ierr);		// set back coordinates
    
	//VecView(Erosion_UserData->Elevation, PETSC_VIEWER_STDOUT_SELF);
    
    
	ierr 	= 	DMDAGetInfo(Erosion_UserData->da_elevation,PETSC_NULL,&Erosion_UserData->M,&Erosion_UserData->N,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);     CHKERRQ(ierr);
    
	ierr					=	DMDAGetBoundingBox(Erosion_UserData->da_elevation,xymin,xymax);                                     CHKERRQ(ierr);
	Erosion_UserData->lx 	= 	xymax[0]-xymin[0];
	Erosion_UserData->ly 	= 	xymax[1]-xymin[1];
	Erosion_UserData->M 	=	M;
	Erosion_UserData->N 	=	N;
	nelx 					=	M-1;
	nely 					=	N-1;
	Erosion_UserData->nelx 	= 	nelx;
    Erosion_UserData->nely 	= 	nely;
    
    
	/* For visualization only */
	ierr = DMDACreate3d(PETSC_COMM_SELF,DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR,M,N,	1,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,1,PETSC_NULL,PETSC_NULL,PETSC_NULL,&Erosion_UserData->da_elevation_output); CHKERRQ(ierr);
    ierr = DMDASetUniformCoordinates(Erosion_UserData->da_elevation_output,0,xymax[0],0,xymax[1],0,1);                          CHKERRQ(ierr);		// set coordinates
    
    
    ierr = DMDACreate2d(PETSC_COMM_SELF,DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,nelx,nely,PETSC_DECIDE,PETSC_DECIDE,4,0,PETSC_NULL,PETSC_NULL,&Erosion_UserData->da_element); CHKERRQ(ierr);
    ierr = DMDASetUniformCoordinates(Erosion_UserData->da_element,0,xymax[0],0,xymax[1],0,0); CHKERRQ(ierr);		// set coordinates
    
    ierr = DMDASetFieldName(Erosion_UserData->da_element,0, "H_el" );           CHKERRQ(ierr);
    ierr = DMDASetFieldName(Erosion_UserData->da_element,1, "k_el" );           CHKERRQ(ierr);
    ierr = DMDASetFieldName(Erosion_UserData->da_element,2, "Q_el" );           CHKERRQ(ierr);
    ierr = DMDASetFieldName(Erosion_UserData->da_element,3, "s_el" );           CHKERRQ(ierr);
    
    
    /* Create global vectors for properties */
    ierr = DMCreateGlobalVector(Erosion_UserData->da_element, &Erosion_UserData->properties);  CHKERRQ(ierr);        // global vector that contains properties data for each element (H_el, k_el, Q_el, s_el)
    PetscObjectSetName((PetscObject)Erosion_UserData->properties,"properties_element");
    
    /* Create global vector for elevation_initial */
    ierr = VecDuplicate(Erosion_UserData->Elevation, &Erosion_UserData->Elevation_init);                                       CHKERRQ(ierr);
    PetscObjectSetName((PetscObject)Erosion_UserData->Elevation_init,"Elevation_initial");
    
    ierr = VecCopy(Erosion_UserData->Elevation, Erosion_UserData->Elevation_init);                                   CHKERRQ(ierr);
    
    
    
	/* time parameters */
	day		= 3600*24;
	year	= 365.25 * day;
    
	/* Landscape evolution parameters */

	rain_m_year             =	user->ErosionParameters.FE_ErosionCode.rain_m_year;
    rain_river_year         =	user->ErosionParameters.FE_ErosionCode.rain_river_year;
	k0                      =	user->ErosionParameters.FE_ErosionCode.k0;
	c                       =	user->ErosionParameters.FE_ErosionCode.c;
	hack                    =	user->ErosionParameters.FE_ErosionCode.n;
    BC                      =   user->ErosionParameters.FE_ErosionCode.BC;
    fill_lake               =   user->ErosionParameters.FE_ErosionCode.fill_lake; 
    mode_river              =   user->ErosionParameters.FE_ErosionCode.mode_river;
    nbre_river              =   user->ErosionParameters.FE_ErosionCode.nbre_river;
    
    for (i=0; i<nbre_river; i++){
        location_river[i] = user->ErosionParameters.FE_ErosionCode.location_river[i];
        Erosion_UserData->location_river[i] = location_river[i];
    }
    
    rain                    = rain_m_year/year;
    rain_river              = rain_river_year/year;
    
	Erosion_UserData->rain            = rain;	// in m-s
    Erosion_UserData->rain_river      = rain_river;
	Erosion_UserData->k0              = k0;
	Erosion_UserData->c               = c;
	Erosion_UserData->hack            = hack;
    Erosion_UserData->dt              = dt_erosion;				
	Erosion_UserData->BC              = BC;		                            
    Erosion_UserData->fill_lake       = fill_lake;
    Erosion_UserData->mode_river      = mode_river;
    Erosion_UserData->nbre_river      = nbre_river; 
    
	PetscFunctionReturn(0);
}
/*==============================================================================================================================*/





/*===============================================================================================================*/
/*
 FUNCTIONS FOR FINITE ELEMENT EROSION CODE AND OTHER ROUTINE
 
 */


/*--------------------------------------------------------------------------------------------------------------------------------------------------------------*/
/*   Boundary conditions: fixed topography in x-directions and no flux in y-directions
 */

#undef __FUNCT__
#define __FUNCT__ "DefineBC"
PetscErrorCode DefineBC( FE_Erosion_UserData *Erosion_UserData, PetscInt bcdof[], PetscScalar bcval[])
{
    PetscErrorCode ierr;
    PetscInt        N, M;
    PetscInt        i,j,id,mstart,m,nstart,n;
    Vec             Elevation;
    Vec             Elevation_local;
    PetscScalar     *Hvec, **H;
    DM              da_elevation;
    
    da_elevation = Erosion_UserData->da_elevation;
    Elevation  =  Erosion_UserData->Elevation;
    M = Erosion_UserData-> M;
    N = Erosion_UserData-> N;
    
    
    /* Define Hvec global vecteur */
    ierr = PetscMalloc((size_t)(M*N)*sizeof(PetscScalar), &Hvec);                                 CHKERRQ(ierr);
    
    /* Copy elevation vector to local processor */
    ierr = DMGetLocalVector(da_elevation,&Elevation_local);                             CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(da_elevation,Elevation,INSERT_VALUES,Elevation_local);	CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(da_elevation,Elevation,INSERT_VALUES,Elevation_local);	CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_elevation,Elevation_local,&H);							CHKERRQ(ierr);
    
    DMDAGetCorners(Erosion_UserData->da_elevation,&mstart,&nstart,PETSC_NULL,&m,&n,PETSC_NULL);
    
    //num = 0;
    for (j=nstart; j<nstart+n; j++) {
        for (i=mstart; i<mstart+m; i++) {
            
            id = i +j*M;
            Hvec[id] = H[j][i];
            //      num = num+1;
        }
    }
    
	/* send elevation data back to global vector */
    ierr = DMDAVecRestoreArray(da_elevation,Elevation_local,&H);											    CHKERRQ(ierr);
    ierr = DMLocalToGlobalBegin(da_elevation,Elevation_local,INSERT_VALUES,Elevation);      CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(da_elevation,Elevation_local,INSERT_VALUES,Elevation);        CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(da_elevation,&Elevation_local);                                               CHKERRQ(ierr);
    
    computeBC (N, M, Hvec, bcdof, bcval, Erosion_UserData->BC);
    
    
    PetscFree(Hvec);
    
    PetscFunctionReturn(0);
}
/*--------------------------------------------------------------------------------------------------------------------------------------------------------------*/

#undef __FUNCT__
#define __FUNCT__ "computeBC"
PetscErrorCode computeBC(PetscInt N, PetscInt M, PetscScalar H[], PetscInt bcdof[], PetscScalar bcval[], PetscInt BC)
{
    PetscInt i,j, id;
    
    
    for (j=0; j<N; j++){
        for (i=0; i<M; i++){
            id = i + j*M;
            
            if (BC==0){
            	//left
            	if (i==0){
            		bcdof[j]=id;
            		bcval[j]=H[id];
            	}
                
            	// right
            	if (i==M-1){
            		bcdof[N+j] = id;
            		bcval[N+j]=H[id];
            	}
            }
            else if (BC==1){
            	// front
            	if (j==0){
            		bcdof[i]=id;
            		bcval[i]=H[id];
            	}
                
            	//back
            	if (j==N-1){
            		bcdof[M+i] = id;
            		bcval[M+i]=H[id];
            	}
            }
            
            else if (BC==2){
                // front 
                if (j==0){
            		bcdof[i]=id;
            		bcval[i]=H[id];
            	}
                // back 
                if (j==N-1){
            		bcdof[M+i] = id;
            		bcval[M+i]=H[id];
            	}
                // left    
                if ((i==0) && (j>0) && (j<N-1)){
            		bcdof[2*M-1+j]=id;
            		bcval[2*M-1+j]=H[id];
            	}
                // right
                if ((i==M-1) && (j>0) && (j<N-1)){
            		bcdof[2*M+N-3+j] = id;
            		bcval[2*M+N-3+j]=H[id];
            	}                
            }
            
            else if (BC==3){
                if (j==0){
            		bcdof[i]=id;
            		bcval[i]=H[id];
            	}
                
            }
            
            else {
                continue;
            }
            
        }
    }
    PetscFunctionReturn(0);
}

/*------------------------------------------------------------------------------------------------------------------------------------------------------------------*/


/*------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
/* Set default parameters & receive input parameters from the command-line if specified
 */
#undef __FUNCT__
#define __FUNCT__ "GetInputParameters"
PetscErrorCode GetInputParameters( FE_Erosion_UserData *Erosion_UserData)
{
	PetscErrorCode ierr;
	PetscScalar    lx=100e3, ly=50e3;				//	length and width of the model
	PetscScalar    Hmax = 400, noise = 1;		//	maximum elevation and random noise
	PetscScalar    day,year,dt, dt_years = 1000; //  needed to compute time in s
	PetscScalar	   rain, k0, c, hack, rain_m_year; // erosion parameters
	PetscInt       M = 101 , N = 51;            // M = # of nodes in x-direction, N = # of nodes in y-direction
    PetscInt       nt = 1, output_freq=1; // # time steps and BC
	PetscBool	   flg;
    
    
	/* Mesh size */
	ierr = PetscOptionsGetInt(PETSC_NULL,"-M",&M,PETSC_NULL);                   CHKERRQ(ierr);		/* number of nodes in x-direction */
	Erosion_UserData->M = M;
    Erosion_UserData->nelx = M-1;
    
    
	ierr = PetscOptionsGetInt(PETSC_NULL,"-N",&N,PETSC_NULL);                   CHKERRQ(ierr); 		/* number of nodes in y-direction */
	Erosion_UserData->N = N;
    Erosion_UserData->nely = N-1;
    
	ierr = PetscOptionsGetScalar(PETSC_NULL,"-lx",&lx,PETSC_NULL);              CHKERRQ(ierr); 		/* length in x- direction [m] */
	Erosion_UserData->lx = lx;
    
    ierr = PetscOptionsGetScalar(PETSC_NULL,"-ly",&ly,PETSC_NULL);              CHKERRQ(ierr); 		/* length in y- direction [m] */
	Erosion_UserData->ly = ly;
    
	/* Other parameters that depend on the initial noise & surface */
    
	/* topological parameters */
    ierr = PetscOptionsGetScalar(PETSC_NULL,"-Hmax",&Hmax,PETSC_NULL);          CHKERRQ(ierr);      /* Maximum elevation of the initial topography */
    Erosion_UserData->Hmax = Hmax;
    
    ierr = PetscOptionsGetScalar(PETSC_NULL,"-noise",&noise,PETSC_NULL);        CHKERRQ(ierr);      /* initial roughness of the topography*/
	Erosion_UserData->noise = noise;
    
	Erosion_UserData->nt	=	nt;
	ierr = PetscOptionsGetInt(PETSC_NULL,"-nt",&nt,&flg);        CHKERRQ(ierr);      /* number of timesteps */
	if (flg){
		Erosion_UserData->nt	=	nt;
	}
    
    
    ierr = PetscOptionsGetInt(PETSC_NULL,"-output_freq",&output_freq,PETSC_NULL);                   CHKERRQ(ierr);		/* number of nodes in x-direction */
	Erosion_UserData->output_freq = output_freq;
    
    
    
	ierr = PetscOptionsGetScalar(PETSC_NULL,"-dt_years",&dt_years,&flg);        CHKERRQ(ierr);      /* number of timesteps */
    
    
	/* time parameters */
    day		= 3600*24;
    year	= 365.25 * day;
    dt		= dt_years * year;		// in seconds
    
    ierr = PetscOptionsGetScalar(PETSC_NULL,"-dt",&dt,PETSC_NULL);              CHKERRQ(ierr);      /* Time Step [s] */
    Erosion_UserData->dt = dt;
    
    //ierr = PetscOptionsGetInt(PETSC_NULL,"-nt",&nt,PETSC_NULL);              CHKERRQ(ierr);      /* Number of integration */
    // Erosion_UserData->nt = nt;
    
    
	/* Landscape evolution parameters */
    rain_m_year	= 	0.3;		// in m/year
    k0 			= 	3.2e-12;
    c 			= 	1;
    hack        = 	2;
    
    ierr = PetscOptionsGetScalar(PETSC_NULL,"-rain_m_year",&rain_m_year,PETSC_NULL);          CHKERRQ(ierr);      /* Effective rain [m/year] average annual precipitation */
	rain = rain_m_year/year;
    
    Erosion_UserData->rain 	= rain;	// in m-s
    
    ierr = PetscOptionsGetScalar(PETSC_NULL,"-k0",&k0,PETSC_NULL);              CHKERRQ(ierr);      /* coeffiscient of diffusivity */
    Erosion_UserData->k0 	= k0;
    
    ierr = PetscOptionsGetScalar(PETSC_NULL,"-c",&c,PETSC_NULL);                CHKERRQ(ierr);      /* fluvial coeffiscient (reflect in part the erodibility of rocks) */
    Erosion_UserData->c 	= c;
    
    ierr = PetscOptionsGetScalar(PETSC_NULL,"-hack",&hack,PETSC_NULL);           CHKERRQ(ierr);        /* power exponent --> related to the Hack law's */
    Erosion_UserData->hack 	= hack;
    
    
	/* Print some useful info to the screen */
	PetscPrintf(PETSC_COMM_SELF," --              FD_ErosionModel            -- \n");
	PetscPrintf(PETSC_COMM_SELF," Dimensions: %1.0fx%1.0f [m   ]	(change with -lx,-ly)   \n", lx,ly);
	PetscPrintf(PETSC_COMM_SELF," Gridpoints: %ix%i [   ]	(change with -M,-N)   \n", M,N);
	PetscPrintf(PETSC_COMM_SELF," Timesteps:  %i           (change with -nt) \n", nt);
	PetscPrintf(PETSC_COMM_SELF," Timestep:   %f [yrs]     (change with -dt_years) \n \n", dt_years);
    
    
    
	PetscPrintf(PETSC_COMM_SELF," Erodability k: k=k0 + c*Q^n, with Q=streamflow [m2/s] \n");
	PetscPrintf(PETSC_COMM_SELF," Rain :                      %1.3f [m/yr]   (change with -rain_m_year) \n", rain_m_year);
	PetscPrintf(PETSC_COMM_SELF," k0 :                        %1.3e [m2/s]   (change with -k0) \n", k0);
	PetscPrintf(PETSC_COMM_SELF," c :                         %1.3e [((m^2/s)^(1-n)]   (change with -c) \n", c);
	PetscPrintf(PETSC_COMM_SELF," hack :                      %1.3f [    ]   (change with -hack) \n", hack);
    
	PetscPrintf(PETSC_COMM_SELF," --------------------------------------------- \n");
    
    
	PetscFunctionReturn(0);
}
/*----------------------------------------------------------------------------------------------------------------------------------------------------------------*/

/*
 *-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 *		Updates the z-coordinates of the topography DA [for nicer visualization]
 */
#undef __FUNCT__
#define __FUNCT__ "UpdateTopography_DA"
PetscErrorCode FE_Erosion_UpdateTopography_DA( FE_Erosion_UserData *Erosion_UserData)
{
	PetscErrorCode	ierr;
	PetscInt		mstart, nstart, m,n,i,j;
	PetscScalar		**H;
	DM 				cda;
	Vec				gc, Elevation_local, global;
	DMDACoor3d       ***coors;
    
    
	/* Set the initial elevation ============================================================*/
    ierr = DMGetCoordinateDM(Erosion_UserData->da_elevation_output,&cda); 				CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(Erosion_UserData->da_elevation_output,&gc);			CHKERRQ(ierr);
    ierr = DMDAVecGetArray(cda,gc,&coors);										    CHKERRQ(ierr);
    
    /* Copy elevation vector to local processor */
	ierr = DMGetLocalVector(Erosion_UserData->da_elevation,&Elevation_local);                                               CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(Erosion_UserData->da_elevation,Erosion_UserData->Elevation,INSERT_VALUES,Elevation_local);	CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(Erosion_UserData->da_elevation,Erosion_UserData->Elevation,INSERT_VALUES,Elevation_local);	CHKERRQ(ierr);
    ierr = DMDAVecGetArray(Erosion_UserData->da_elevation,Elevation_local,&H);
    
    
    DMDAGetCorners(cda,&mstart,&nstart,PETSC_NULL,&m,&n,PETSC_NULL);
    for (i=mstart; i<mstart+m; i++) {
    	for (j=nstart; j<nstart+n; j++) {
            
    		coors[0][j][i].z = H[j][i];
            
            ///	PetscPrintf(PETSC_COMM_WORLD,"coors[0][%i][%i] = %f \n", coors[0][j][i].z,j,i);
            
            
    	}
    }
    ierr = DMDAVecRestoreArray(cda,gc,&coors);										 CHKERRQ(ierr);
    
    
	DMGetCoordinates(Erosion_UserData->da_elevation_output,&global);
	DMLocalToGlobalBegin(cda,gc,INSERT_VALUES,global);
	DMLocalToGlobalEnd(cda,gc,INSERT_VALUES,global);
    
    
    /* send elevation data back to global vector */
    ierr = DMDAVecRestoreArray(Erosion_UserData->da_elevation,Elevation_local,&H);											    CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(Erosion_UserData->da_elevation,&Elevation_local);                                              CHKERRQ(ierr);
    //ierr = VecDestroy(&Elevation_local);
    
    
    
	PetscFunctionReturn(0);
    
}
/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/

#undef __FUNCT__
#define __FUNCT__ "Solve_NonlinearDiffusion"
PetscErrorCode Solve_NonlinearDiffusion( FE_Erosion_UserData *Erosion_UserData, PetscInt bcdof[], PetscScalar bcval[])
{
	PetscErrorCode ierr;
    Vec            b, FG;                /*  RHS */
    Mat            KLG, KRG;             /* linear system matrix */
    KSP            ksp;                  /* linear solver context */
	PC			   pc;                   /* preconditioner context */
    // PetscLogStage  stage;
    
    Vec             Elevation;
    DM              da_elevation;
    PetscInt        M, N;
    
    da_elevation	=  Erosion_UserData->da_elevation;
    Elevation		=  Erosion_UserData->Elevation;
    M               =  Erosion_UserData->M;
    N               =  Erosion_UserData->N;
    
    
    /* Generate a matrix with the correct non-zero pattern of type AIJ. This will work in parallel and serial */
	ierr =  DMCreateMatrix(da_elevation,&KRG);                         CHKERRQ(ierr);
    ierr =  DMCreateMatrix(da_elevation,&KLG);                         CHKERRQ(ierr);
    ierr =  DMCreateGlobalVector(da_elevation,&FG);                           CHKERRQ(ierr);
    
    ierr = Assemble_matrixes_vector( Erosion_UserData, KLG, KRG, FG);        CHKERRQ(ierr);
    
    /* form vector b to solve the system */
    ierr = DMCreateGlobalVector(da_elevation,&b);                   CHKERRQ(ierr);
    ierr = MatMult(KRG,Elevation,b);                                CHKERRQ(ierr);
    ierr = VecAXPY(b,1.0,FG);                                      CHKERRQ(ierr);
    
    
    /* Apply Boundary conditions */
    
    // in case we fix left & right boundaries
    if (Erosion_UserData->BC==0){
    	ierr = MatZeroRows(KLG, 2*N,bcdof, 1, PETSC_NULL, PETSC_NULL);  CHKERRQ(ierr);
    	ierr = VecSetValues(b, 2*N, bcdof, bcval, INSERT_VALUES);       CHKERRQ(ierr);
    }
    // in case we fix front & back boundaries
    else if (Erosion_UserData->BC==1){
    	ierr = MatZeroRows(KLG, 2*M,bcdof, 1, PETSC_NULL, PETSC_NULL);  CHKERRQ(ierr);
    	ierr = VecSetValues(b, 2*M, bcdof, bcval, INSERT_VALUES);       CHKERRQ(ierr);
    }
    // in case we fix front, back, left and right boundaries
    else if (Erosion_UserData->BC==2){
    	ierr = MatZeroRows(KLG, 2*(M+N-2),bcdof, 1, PETSC_NULL, PETSC_NULL);  CHKERRQ(ierr);
    	ierr = VecSetValues(b, 2*(M+N-2), bcdof, bcval, INSERT_VALUES);       CHKERRQ(ierr);
    }
    
    // in case we fix front and let the backside with no flux
    
    else if (Erosion_UserData->BC==3){
        ierr = MatZeroRows(KLG, M,bcdof, 1, PETSC_NULL, PETSC_NULL);  CHKERRQ(ierr);
    	ierr = VecSetValues(b, M, bcdof, bcval, INSERT_VALUES);       CHKERRQ(ierr);
    }
    
    
    
    ierr = VecAssemblyBegin(b);                                     CHKERRQ(ierr);
    ierr = VecAssemblyEnd(b);                                       CHKERRQ(ierr);
    
    
    /*  printf("matrix KLG after applying BC");
     ierr = MatView(KLG, PETSC_VIEWER_STDOUT_WORLD);                 CHKERRQ(ierr);
     printf("matrix KRG after applying BC");
     ierr = MatView(KRG, PETSC_VIEWER_STDOUT_WORLD);                 CHKERRQ(ierr);
     printf("vector b after BC");
     ierr = VecView(b, PETSC_VIEWER_STDOUT_WORLD);                   CHKERRQ(ierr);*/
    
    /*  sprintf(name,"matrix.mat");
     ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF, name, &viewer);                CHKERRQ(ierr);
     ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);                CHKERRQ(ierr);
     ierr = MatView(KLG, viewer);                CHKERRQ(ierr);
     ierr = PetscViewerDestroy( &viewer );             CHKERRQ(ierr);         */
    
    /*ierr  = PetscViewerBinaryOpen(PETSC_COMM_SELF,"kr.pmat",FILE_MODE_WRITE,&viewer); CHKERRQ(ierr);
     ierr = MatView(KRG, viewer);  CHKERRQ(ierr);
     ierr = PetscViewerDestroy( &viewer );             CHKERRQ(ierr);
     
     ierr  = PetscViewerBinaryOpen(PETSC_COMM_SELF,"kl.pmat",FILE_MODE_WRITE,&viewer); CHKERRQ(ierr);
     ierr = MatView(KLG, viewer);  CHKERRQ(ierr);
     ierr = PetscViewerDestroy( &viewer );             CHKERRQ(ierr);     */
    
    /*------------------------------------------------------
     Create the linear solver and set various options
     ------------------------------------------------------*/
    
    /* Create linear solver context */
    ierr = KSPCreate(PETSC_COMM_SELF,&ksp);                         CHKERRQ(ierr);
    ierr = KSPSetOptionsPrefix(ksp,"ErosionModel_");				CHKERRQ(ierr);		// set prefix such that we can specify solver options from the command line
    
    
    /* Set operators */
    ierr =  KSPSetOperators(ksp,KLG,KLG);      CHKERRQ(ierr);
    
    /* Set runtime option */
    ierr = KSPSetType(ksp, "preonly");								CHKERRQ(ierr);		// default solver
    ierr = KSPSetFromOptions(ksp);                                  CHKERRQ(ierr);		// allow command-line override
    
    ierr = KSPGetPC(ksp,&pc);										CHKERRQ(ierr);
    ierr = PCSetType(pc,PCLU);										CHKERRQ(ierr);
    ierr = PCSetFromOptions(pc);									CHKERRQ(ierr);		// allow command-line override
    
    
    /*------------------------------------------------------
     Solve the linear system
     ------------------------------------------------------*/
    
    ierr =  KSPSolve(ksp,b,Elevation);                              CHKERRQ(ierr);
    
    
    
    
    ierr = MatDestroy(&KLG);                                       CHKERRQ(ierr);
    ierr = MatDestroy(&KRG);                                       CHKERRQ(ierr);
    ierr = VecDestroy(&FG);                                        CHKERRQ(ierr);
    ierr = VecDestroy(&b);                                         CHKERRQ(ierr);
    ierr = KSPDestroy(&ksp);                                       CHKERRQ(ierr);
    
    
    PetscFunctionReturn(0);
    
}
/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/

#undef __FUNCT__
#define __FUNCT__ "Assemble_matrixes_vector"
PetscErrorCode Assemble_matrixes_vector(FE_Erosion_UserData *Erosion_UserData, Mat KLG, Mat KRG, Vec FG)
{
    
    PetscErrorCode ierr;
    PetscScalar	    dt;
    Vec             properties, l_properties;
    Vec             gc;
    DM              da_elevation, cda, da_element;
    DMDACoor2d      **coors;
    PetscInt        nx,ny;
    
    const PetscInt *iglobal;
    PetscInt        sex,sey, nelx, nely, ei, ej, i, j, id;
    PetscInt        n_int = 4;
    PetscScalar     k_el, s_el;
    PetscScalar     MM[16], KM[16], F[4], KL[16], KR[16];
    PetscScalar     coord_el[8];
    PetscInt        row[4];
    FE_ElementInfo  **element_info;
    ISLocalToGlobalMapping ltogm;

    // variables associated to the structure Erosion_UserData
    
    dt				=	Erosion_UserData->dt;
    da_elevation	=   Erosion_UserData->da_elevation;
    da_element      =   Erosion_UserData->da_element;
    properties      =   Erosion_UserData->properties;
    
    
    
    /* Access the da_elevation (nodes) to get the coordinates */
    ierr = DMGetCoordinateDM(da_elevation,&cda);                                                              CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(da_elevation,&gc);                                                         CHKERRQ(ierr);
    ierr = DMDAVecGetArray(cda,gc,&coors);                                                                      CHKERRQ(ierr);
    
    /* Access the da_element to get the element properties */
	ierr = DMGetLocalVector(da_element, &l_properties );                                                  CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da_element, properties, INSERT_VALUES, l_properties );                    CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(da_element, properties, INSERT_VALUES, l_properties );                      CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da_element, l_properties, &element_info );                                       CHKERRQ(ierr);
    
    
    DMDAGetGhostCorners(da_elevation,0,0,0,&nx,&ny,0);
    DAGetElementCorners( da_elevation, &sex, &sey, &nelx, &nely, 0 );
	for( ej=sey; ej<sey+nely; ej++ ) {
		for( ei=sex; ei<sex+nelx; ei++ ) {
            
            /* get coordinates for elements */
            
            ierr = GetElementCoord(coord_el, ei, ej, coors);                                                    CHKERRQ(ierr);
            k_el = element_info[ej][ei].k_el;
            s_el = 0.0;
            
            //  printf("(%d,%d) k_el: %1.15e, s_el: %1.10e\n", ei,ej,k_el, s_el);
            
            /*  printf("element\n") */
            
            
            /* Initialize vectors */
            ierr = PetscMemzero( F, sizeof(PetscScalar)*4);                         CHKERRQ(ierr);
            ierr = PetscMemzero( KM, sizeof(PetscScalar)*4*4 );                     CHKERRQ(ierr);
            ierr = PetscMemzero( MM, sizeof(PetscScalar)*4*4 );                     CHKERRQ(ierr);
            ierr = PetscMemzero( KL, sizeof(PetscScalar)*4*4 );                     CHKERRQ(ierr);
            ierr = PetscMemzero( KR, sizeof(PetscScalar)*4*4 );                     CHKERRQ(ierr);
            
            
            ierr = Evaluate_Integrale( n_int, k_el, s_el, coord_el, MM, KM, F);     CHKERRQ(ierr);
            
            /* Compute KL and KR to insert in the global matrixes */
            for (j=0; j<4; j++){
                for (i=0; i<4; i++){
                    id = i*4 + j;
                    KL[id] = MM[id] * 1.0/dt + KM[id];
                    KR[id]= MM[id] * 1.0/dt;
                    //  printf("KL: %.15f, KR: %.15f\n",KL[id], KR[id]);
                    
                }
            }
            
            
            /* insert element matrix into global matrix */
//          ierr = DMDAGetGlobalIndices(da_elevation, PETSC_NULL, &iglobal);        CHKERRQ(ierr);

        	ierr = DMGetLocalToGlobalMapping(da_elevation, &ltogm);   CHKERRQ(ierr);
        	ierr = ISLocalToGlobalMappingGetIndices(ltogm, &iglobal); CHKERRQ(ierr);

            ierr = ConstructRow( nx, sex, sey, ei, ej, iglobal, row);               CHKERRQ(ierr);
            
            // printf(" ei,ej %d,%d : (%d  , %d , %d , %d)\n",ei,ej,row[0], row[1], row[2], row[3]);
            ierr = MatSetValues( KLG, 4, row, 4, row, KL, ADD_VALUES);              CHKERRQ(ierr);
            ierr = MatSetValues(KRG, 4, row, 4, row, KR, ADD_VALUES);               CHKERRQ(ierr);
            ierr = VecSetValues(FG, 4, row, F, ADD_VALUES);                         CHKERRQ(ierr);
        }
    }
    
    
    ierr = MatAssemblyBegin(KLG, MAT_FINAL_ASSEMBLY);                               CHKERRQ(ierr);
	ierr = MatAssemblyEnd(KLG, MAT_FINAL_ASSEMBLY);                                 CHKERRQ(ierr);
    
    ierr = MatAssemblyBegin(KRG, MAT_FINAL_ASSEMBLY);                               CHKERRQ(ierr);
	ierr = MatAssemblyEnd(KRG, MAT_FINAL_ASSEMBLY);                                 CHKERRQ(ierr);
    
    ierr = VecAssemblyBegin(FG);                                                    CHKERRQ(ierr);
	ierr = VecAssemblyEnd(FG);                                                      CHKERRQ(ierr);
    
    /* printf("matrix KLG before applying BC");
     ierr = MatView(KLG, PETSC_VIEWER_STDOUT_WORLD);                                 CHKERRQ(ierr);
     printf("matrix KRG before applying BC");
     ierr = MatView(KRG, PETSC_VIEWER_STDOUT_WORLD);                                 CHKERRQ(ierr);
     printf("vector FG before applying BC");
     ierr = VecView(FG, PETSC_VIEWER_STDOUT_WORLD);                                  CHKERRQ(ierr);
     */
    /* send elevation and coordinates data back to global vector */
    
    ierr = DMDAVecRestoreArray(cda,gc,&coors);                                      CHKERRQ(ierr);
    
    
    ierr = DMDAVecRestoreArray(da_element,l_properties,&element_info);				CHKERRQ(ierr);
    ierr = DMLocalToGlobalBegin(da_element,l_properties,INSERT_VALUES,properties); 	CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(da_element,l_properties,INSERT_VALUES,properties);    CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(da_element,&l_properties);                          CHKERRQ(ierr);
    
    
    PetscFunctionReturn(0);
}

/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/

#undef __FUNCT__
#define __FUNCT__ "Solve_linearDiffusion"
PetscErrorCode Solve_linearDiffusion( FE_Erosion_UserData *Erosion_UserData, PetscInt bcdof[], PetscScalar bcval[])
{
	PetscErrorCode ierr;
    Vec            b, FG;                /*  RHS */
    Mat            KLG, KRG;             /* linear system matrix */
    KSP            ksp;                  /* linear solver context */
	PC			   pc;                   /* preconditioner context */
    // PetscLogStage  stage;
    
    Vec             Elevation;
    DM              da_elevation; 
    PetscInt        M, N;
    da_elevation	=  Erosion_UserData->da_elevation;
    Elevation		=  Erosion_UserData->Elevation;
    M               =  Erosion_UserData->M;
    N               =  Erosion_UserData->N; 
    
    
    
    
    
    /* Generate a matrix with the correct non-zero pattern of type AIJ. This will work in parallel and serial */
	ierr =  DMCreateMatrix(da_elevation,&KRG);                         CHKERRQ(ierr);
    ierr =  DMCreateMatrix(da_elevation,&KLG);                         CHKERRQ(ierr);
    ierr =  DMCreateGlobalVector(da_elevation,&FG);                           CHKERRQ(ierr); 
    
    ierr = Assemble_matrixes_vector_smoother( Erosion_UserData, KLG, KRG, FG);        CHKERRQ(ierr);
    
    /* form vector b to solve the system */
    ierr = DMCreateGlobalVector(da_elevation,&b);                   CHKERRQ(ierr);
    ierr = MatMult(KRG,Elevation,b);                                CHKERRQ(ierr);
    ierr = VecAXPY(b,1.0,FG);                                      CHKERRQ(ierr);
    
    /* Apply Boundary conditions */
    
    // in case we fix left & right boundaries
    if (Erosion_UserData->BC==0){
    	ierr = MatZeroRows(KLG, 2*N,bcdof, 1, PETSC_NULL, PETSC_NULL);  CHKERRQ(ierr);
    	ierr = VecSetValues(b, 2*N, bcdof, bcval, INSERT_VALUES);       CHKERRQ(ierr);
    }
    // in case we fix front & back boundaries
    else if (Erosion_UserData->BC==1){
    	ierr = MatZeroRows(KLG, 2*M,bcdof, 1, PETSC_NULL, PETSC_NULL);  CHKERRQ(ierr);
    	ierr = VecSetValues(b, 2*M, bcdof, bcval, INSERT_VALUES);       CHKERRQ(ierr);
    }
    // in case we fix front, back, left and right boundaries
    else if (Erosion_UserData->BC==2){
    	ierr = MatZeroRows(KLG, 2*(M+N-2),bcdof, 1, PETSC_NULL, PETSC_NULL);  CHKERRQ(ierr);
    	ierr = VecSetValues(b, 2*(M+N-2), bcdof, bcval, INSERT_VALUES);       CHKERRQ(ierr);
    }
    
    // in case we fix front and let the backside with no flux
    
    else if (Erosion_UserData->BC==3){
        ierr = MatZeroRows(KLG, M,bcdof, 1, PETSC_NULL, PETSC_NULL);  CHKERRQ(ierr);
    	ierr = VecSetValues(b, M, bcdof, bcval, INSERT_VALUES);       CHKERRQ(ierr);
    }
    
    
    ierr = VecAssemblyBegin(b);                                     CHKERRQ(ierr);
    ierr = VecAssemblyEnd(b);                                       CHKERRQ(ierr);
    
    
    /*  printf("matrix KLG after applying BC");
     ierr = MatView(KLG, PETSC_VIEWER_STDOUT_WORLD);                 CHKERRQ(ierr); 
     printf("matrix KRG after applying BC");
     ierr = MatView(KRG, PETSC_VIEWER_STDOUT_WORLD);                 CHKERRQ(ierr); 
     printf("vector b after BC");
     ierr = VecView(b, PETSC_VIEWER_STDOUT_WORLD);                   CHKERRQ(ierr);*/
    
    /*  sprintf(name,"matrix.mat");
     ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF, name, &viewer);                CHKERRQ(ierr);
     ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);                CHKERRQ(ierr);
     ierr = MatView(KLG, viewer);                CHKERRQ(ierr);
     ierr = PetscViewerDestroy( &viewer );             CHKERRQ(ierr);         */
    
    /*ierr  = PetscViewerBinaryOpen(PETSC_COMM_SELF,"kr.pmat",FILE_MODE_WRITE,&viewer); CHKERRQ(ierr);
     ierr = MatView(KRG, viewer);  CHKERRQ(ierr);
     ierr = PetscViewerDestroy( &viewer );             CHKERRQ(ierr);     
     
     ierr  = PetscViewerBinaryOpen(PETSC_COMM_SELF,"kl.pmat",FILE_MODE_WRITE,&viewer); CHKERRQ(ierr);
     ierr = MatView(KLG, viewer);  CHKERRQ(ierr);
     ierr = PetscViewerDestroy( &viewer );             CHKERRQ(ierr);     */
    
    /*------------------------------------------------------
     Create the linear solver and set various options 
     ------------------------------------------------------*/
    
    /* Create linear solver context */
    ierr = KSPCreate(PETSC_COMM_SELF,&ksp);                         CHKERRQ(ierr);
    ierr = KSPSetOptionsPrefix(ksp,"ErosionModel_");				CHKERRQ(ierr);		// set prefix such that we can specify solver options from the command line
    
    
    /* Set operators */
    ierr =  KSPSetOperators(ksp,KLG,KLG);      CHKERRQ(ierr);
    
    /* Set runtime option */
    ierr = KSPSetType(ksp, "preonly");								CHKERRQ(ierr);		// default solver
    ierr = KSPSetFromOptions(ksp);                                  CHKERRQ(ierr);		// allow command-line override
    
    ierr = KSPGetPC(ksp,&pc);										CHKERRQ(ierr);	
    ierr = PCSetType(pc,PCLU);										CHKERRQ(ierr);		
    ierr = PCSetFromOptions(pc);									CHKERRQ(ierr);		// allow command-line override
    
    
    /*------------------------------------------------------
     Solve the linear system
     ------------------------------------------------------*/
    
    ierr =  KSPSolve(ksp,b,Elevation);                              CHKERRQ(ierr);  
    
    
    
    
    ierr = MatDestroy(&KLG);                                       CHKERRQ(ierr);
    ierr = MatDestroy(&KRG);                                       CHKERRQ(ierr);
    ierr = VecDestroy(&FG);                                        CHKERRQ(ierr);
    ierr = VecDestroy(&b);                                         CHKERRQ(ierr);
    ierr = KSPDestroy(&ksp);                                       CHKERRQ(ierr);
    
    
    PetscFunctionReturn(0);
    
}

/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/


#undef __FUNCT__
#define __FUNCT__ "Assemble_matrixes_vector_smoother"
PetscErrorCode Assemble_matrixes_vector_smoother(FE_Erosion_UserData *Erosion_UserData, Mat KLG, Mat KRG, Vec FG)    
{
    
    PetscErrorCode ierr;    
    PetscScalar	    dt, rain;
    Vec             gc; 
    DM              da_elevation, cda;
    DMDACoor2d      **coors;
    PetscInt        M, nx,ny;
    
    const PetscInt  *iglobal;
    PetscInt        sex,sey, nelx, nely, ei, ej, i, j, id; 
    PetscInt        n_int = 4;
    PetscScalar     k_el, s_el, lx, dx;
    PetscScalar     MM[16], KM[16], F[4], KL[16], KR[16];
    PetscScalar     coord_el[8];
    PetscInt        row[4];
	ISLocalToGlobalMapping ltogm;

    // variables associated to the structure Erosion_UserData 
    
    M				=	Erosion_UserData->M;
    dt				=	Erosion_UserData->dt;
    da_elevation	=   Erosion_UserData->da_elevation;
    lx              =   Erosion_UserData->lx;
    rain            =   Erosion_UserData->rain; 
    
    dx = lx/(M-1);
    
    /* Access the da_elevation (nodes) to get the coordinates */     
    ierr = DMGetCoordinateDM(da_elevation,&cda);                                                              CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(da_elevation,&gc);                                                         CHKERRQ(ierr);
    ierr = DMDAVecGetArray(cda,gc,&coors);                                                                      CHKERRQ(ierr);
    
    DMDAGetGhostCorners(da_elevation,0,0,0,&nx,&ny,0);    
    DAGetElementCorners(da_elevation, &sex, &sey, &nelx, &nely, 0);
	for( ej=sey; ej<sey+nely; ej++ ) {
		for( ei=sex; ei<sex+nelx; ei++ ) {
            
            /* get coordinates for elements */
            
            ierr = GetElementCoord(coord_el, ei, ej, coors);                                                    CHKERRQ(ierr);
            k_el = 0.0001*dx*rain;
            s_el = 0;
            
            //  printf("(%d,%d) k_el: %1.15e, s_el: %1.10e\n", ei,ej,k_el, s_el);
            
            /*  printf("element\n") */
            
            
            /* Initialize vectors */
            ierr = PetscMemzero( F, sizeof(PetscScalar)*4);                         CHKERRQ(ierr);
            ierr = PetscMemzero( KM, sizeof(PetscScalar)*4*4 );                     CHKERRQ(ierr);  
            ierr = PetscMemzero( MM, sizeof(PetscScalar)*4*4 );                     CHKERRQ(ierr);
            ierr = PetscMemzero( KL, sizeof(PetscScalar)*4*4 );                     CHKERRQ(ierr);
            ierr = PetscMemzero( KR, sizeof(PetscScalar)*4*4 );                     CHKERRQ(ierr);
            
            
            ierr = Evaluate_Integrale( n_int, k_el, s_el, coord_el, MM, KM, F);     CHKERRQ(ierr);
            
            /* Compute KL and KR to insert in the global matrixes */
            for (j=0; j<4; j++){
                for (i=0; i<4; i++){
                    id = i*4 + j;
                    KL[id] = MM[id] * 1.0/dt + KM[id];
                    KR[id]= MM[id] * 1.0/dt;
                    //  printf("KL: %.15f, KR: %.15f\n",KL[id], KR[id]);
                    
                }
            }
            
            
            /* insert element matrix into global matrix */
//          ierr = DMDAGetGlobalIndices(da_elevation, PETSC_NULL, &iglobal);        CHKERRQ(ierr);

        	ierr = DMGetLocalToGlobalMapping(da_elevation, &ltogm);   CHKERRQ(ierr);
        	ierr = ISLocalToGlobalMappingGetIndices(ltogm, &iglobal); CHKERRQ(ierr);

            ierr = ConstructRow( nx, sex, sey, ei, ej, iglobal, row);               CHKERRQ(ierr);
            
            // printf(" ei,ej %d,%d : (%d  , %d , %d , %d)\n",ei,ej,row[0], row[1], row[2], row[3]); 
            ierr = MatSetValues( KLG, 4, row, 4, row, KL, ADD_VALUES);              CHKERRQ(ierr);
            ierr = MatSetValues(KRG, 4, row, 4, row, KR, ADD_VALUES);               CHKERRQ(ierr);
            ierr = VecSetValues(FG, 4, row, F, ADD_VALUES);                         CHKERRQ(ierr);
        }
    }
    
    
    ierr = MatAssemblyBegin(KLG, MAT_FINAL_ASSEMBLY);                               CHKERRQ(ierr);
	ierr = MatAssemblyEnd(KLG, MAT_FINAL_ASSEMBLY);                                 CHKERRQ(ierr);
    
    ierr = MatAssemblyBegin(KRG, MAT_FINAL_ASSEMBLY);                               CHKERRQ(ierr);
	ierr = MatAssemblyEnd(KRG, MAT_FINAL_ASSEMBLY);                                 CHKERRQ(ierr);
    
    ierr = VecAssemblyBegin(FG);                                                    CHKERRQ(ierr);
	ierr = VecAssemblyEnd(FG);                                                      CHKERRQ(ierr);
    
    /* printf("matrix KLG before applying BC");
     ierr = MatView(KLG, PETSC_VIEWER_STDOUT_WORLD);                                 CHKERRQ(ierr); 
     printf("matrix KRG before applying BC");
     ierr = MatView(KRG, PETSC_VIEWER_STDOUT_WORLD);                                 CHKERRQ(ierr); 
     printf("vector FG before applying BC");
     ierr = VecView(FG, PETSC_VIEWER_STDOUT_WORLD);                                  CHKERRQ(ierr);
     */
    /* send elevation and coordinates data back to global vector */
    
    ierr = DMDAVecRestoreArray(cda,gc,&coors);                                      CHKERRQ(ierr);
    
    
    PetscFunctionReturn(0);
}

/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/


/*=========== COMPUTE FEM ROUTINE ===================*/

/* Compute FEM routine */
/*
 Local numbering
 
 3-------2
 |       |
 |       |
 |       |
 0-------1
 
 */



#undef __FUNCT__
#define __FUNCT__ "GetShapeFunction"
PetscErrorCode GetShapeFunction( PetscScalar coord_int0 ,PetscScalar coord_int1 , PetscScalar Ni[], PetscScalar dNi[2][4])
{
    PetscScalar xi = coord_int0;
    PetscScalar eta = coord_int1;
    
    Ni[0] = 0.25*(1-xi)*(1-eta);
    Ni[1] = 0.25*(1+xi)*(1-eta);
    Ni[2] = 0.25*(1+xi)*(1+eta);
    Ni[3] = 0.25*(1-xi)*(1+eta);
    
    
    dNi[0][0] = -0.25 * (1-eta);
    dNi[0][1] = 0.25 * (1-eta);
    dNi[0][2] = 0.25 * (1+eta);
    dNi[0][3] = -0.25 * (1+eta);
    
    dNi[1][0] = -0.25 * (1-xi);
    dNi[1][1] = -0.25 * (1+xi);
    dNi[1][2] = 0.25 * (1+xi);
    dNi[1][3] = 0.25 * (1-xi);
    
    PetscFunctionReturn(0);
}
/*---------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/

#undef __FUNCT__
#define __FUNCT__ "ConstructDetJacobian"
PetscErrorCode ConstructDetJacobian( PetscScalar dNi[2][4], PetscScalar coords[], PetscScalar *det_J)
{
    PetscScalar J00,J01,J10,J11,J;
    PetscInt i, npe=4;
    J00=J01=J10=J11=J=0;
    for (i=0; i<npe; i++){
        PetscScalar cx = coords [2*i];
        PetscScalar cy = coords [2*i+1];
        
        J00 = J00 + dNi[0][i] * cx;
        J01 = J01 + dNi[1][i] * cx;
        J10 = J10 + dNi[0][i] * cy;
        J11 = J11 + dNi[1][i] * cy;
    }
    
    J = (J00*J11) - (J01*J10);
    
    *det_J = J;
    
    PetscFunctionReturn(0);
}
/*-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/

#undef __FUNCT__
#define __FUNCT__ "ConstructdNdx"
PetscErrorCode ConstructdNdx( PetscScalar dNi[2][4], PetscScalar coords[], PetscScalar dNdx[2][4])
{
    PetscScalar J00,J01,J10,J11,J;
    PetscScalar iJ00,iJ01,iJ10,iJ11;
    PetscInt i, npe=4;
    J00=J01=J10=J11=J=0;
    for (i=0; i<npe; i++){
        PetscScalar cx = coords [2*i];
        PetscScalar cy = coords [2*i+1];
        
        J00 = J00 + dNi[0][i] * cx;
        J01 = J01 + dNi[1][i] * cx;
        J10 = J10 + dNi[0][i] * cy;
        J11 = J11 + dNi[1][i] * cy;
    }
    
    J = (J00*J11)-(J01*J10);
    
    iJ00 = J00/J;
    iJ01 = J01/J;
    iJ10 = J10/J;
    iJ11 = J11/J;
    
    for (i=0; i<npe; i++){
        dNdx[0][i] = dNi[0][i] * iJ00 + dNi[1][i] * iJ01;
        dNdx[1][i] = dNi[0][i] * iJ10 + dNi[1][i] * iJ11;
    }
    
    PetscFunctionReturn(0);
}
/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/

/* Compute poisson coeffiscient - loop over all integration points - */

#undef __FUNCT__
#define __FUNCT__ "Evaluate_Integrale"
PetscErrorCode Evaluate_Integrale (PetscInt n_int, PetscScalar k_el, PetscScalar s_el, PetscScalar coord_el[], PetscScalar MM[], PetscScalar KM[], PetscScalar F[])
{
    PetscErrorCode ierr;
    PetscInt i_int, i,j, npe=4;
    PetscScalar Ni[4],weight[4];
    
    //  static double **dNi,**dNdx;
    // 	just use local arrays instead
	PetscScalar dNi[2][4], dNdx[2][4];
    
    PetscScalar coord_int[4][2];
    PetscScalar det_J;
    
    //  don't need this here anymore
    //  static PetscInt been_here = 0;
    
    
    /* Allocate space if 1st time 
     
     // this is a crapppppy memory leak!
     
     if( been_here == 0 ) {
     dNi = (double**)malloc( sizeof(double*) * 2 );
     dNdx = (double**)malloc( sizeof(double*) * 2 );
     for( d=0; d<2; d++ ) {
     dNi[d] = (double*)malloc( sizeof(double) * 4 );
     dNdx[d] = (double*)malloc( sizeof(double) * 4 );
     }
     been_here = 1;
     }
     
     */
    
    /* Define quadrature rule for FEM */
    
    coord_int[0][0]= -sqrt(1.0/3.0);       coord_int[0][1]=-sqrt(1.0/3.0);
    coord_int[1][0]= sqrt(1.0/3.0);        coord_int[1][1]=-sqrt(1.0/3.0);
    coord_int[2][0]= sqrt(1.0/3.0);        coord_int[2][1]=sqrt(1.0/3.0);
    coord_int[3][0]= -sqrt(1.0/3.0);       coord_int[3][1]=sqrt(1.0/3.0);
    
    weight[0]=1.0;
    weight[1]=1.0;
    weight[2]=1.0;
    weight[3]=1.0;
    
    
    
    for (i_int=0; i_int<n_int; i_int++){
        
        
        ierr = GetShapeFunction(coord_int[i_int][0] , coord_int[i_int][1] , Ni, dNi);                   CHKERRQ(ierr);
        ierr = ConstructDetJacobian( dNi, coord_el, &det_J);                                       CHKERRQ(ierr);
        ierr = ConstructdNdx( dNi, coord_el, dNdx);                                                 CHKERRQ(ierr);
        
        
        
        /*  for (i=0; i<4; i++){
         printf(" Ni[%d]: %f\n",i,Ni[i]);
         }
         
         for (j=0; j<2; j++){
         for(i=0; i<4; i++){
         printf("dNi[%d][%d]: %f\n",j,i,dNi[j][i]);
         }
         }
         
         for (j=0; j<2; j++){
         for(i=0; i<4; i++){
         printf("dNdx[%d][%d]: %f\n",j,i,dNdx[j][i]);
         }
         }
         */
        
        /* Compute MM Vector */
        
        for (j=0; j<4; j++){
            for (i=0; i<4; i++){
                MM[i*4+j] = MM[i*4+j] + Ni[i]*Ni[j]*det_J*weight[i_int];
                //printf("MM : %.20f\n ",MM[i*4+j]);
            }
        }
        
        
        /* Compute KM Vector */
        
        for (j=0; j<4; j++){
            for (i=0; i<4; i++){
                KM[i*4+j] = KM[i*4+j]  + k_el*(dNdx[0][j]*dNdx[0][i]+dNdx[1][j]*dNdx[1][i])*det_J*weight[i_int];
                //  printf("KM : %.20f\n ",KM[i*4+j]);
            }
        }
        
        
        /* compute vector F */
        
        for (i=0; i<npe; i++){
            F[i]=F[i]+s_el*Ni[i]*det_J*weight[i_int];
            // printf("F : %.8f\n ",F[i]);
        }
        
        /* visualization matrixes : debbuging*/
        
        
    }
    /*for (j=0; j<4; j++){
     for (i=0; i<4; i++){
     printf("MM : %.20f, KM : %.20f\n ",MM[i*4+j], KM[i*4+j]);
     }
     }*/
    
    PetscFunctionReturn(0);
}
/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/

#undef __FUNCT__
#define __FUNCT__ "GetElementCoord"
PetscErrorCode GetElementCoord(PetscScalar coord_el[], PetscInt i, PetscInt j, DMDACoor2d **coors)
{
    
    
    coord_el[0]= coors[j][i].x;                     coord_el[1]= coors[j][i].y;
    coord_el[2]= coors[j][i+1].x;                   coord_el[3]= coors[j][i+1].y;
    coord_el[4]= coors[j+1][i+1].x;                 coord_el[5]= coors[j+1][i+1].y;
    coord_el[6]= coors[j+1][i].x;                   coord_el[7]= coors[j+1][i].y;
    
    PetscFunctionReturn(0);
}



/*=========== FUNCTIONS THAT DEALS WITH ELEMENTS IN PARALLEL ===================*/

#undef __FUNCT__
#define __FUNCT__ "DAGetElementCorners"
PetscErrorCode DAGetElementCorners( DM da, PetscInt *sx, PetscInt *sy,PetscInt *mx, PetscInt *my,PetscInt *mz)
{
    PetscErrorCode ierr;
	PetscInt si,sj,sk;
    
    
	ierr = DMDAGetGhostCorners(da,&si,&sj,&sk,0,0,0);  CHKERRQ(ierr);
    
	*sx = si;
	if( si != 0 ) {
		*sx = si + 1;
	}
    
	*sy = sj;
	if( sj != 0 ) {
		*sy = sj + 1;
	}
    
	ierr = DAGetLocalElementSize( da, mx, my ,mz);  CHKERRQ(ierr);
    
	PetscFunctionReturn(0);
}
/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/



/* procs to the left claim the ghost node as their element */
#undef __FUNCT__
#define __FUNCT__ "DAGetLocalElementSize"
PetscErrorCode DAGetLocalElementSize( DM da, PetscInt *mxl, PetscInt *myl,PetscInt *mzl)
{
    PetscErrorCode ierr;
	PetscInt m,n,p, M,N,P;
    /*	PetscInt MM,NN,PP; */
	PetscInt sx,sy,sz;
    
    
	ierr = DMDAGetInfo( da, 0, &M,&N,&P, 0,0,0, 0,0,0,0,0,0 );      CHKERRQ(ierr);
	ierr = DMDAGetCorners( da, &sx,&sy,&sz, &m,&n,&p );         CHKERRQ(ierr);
    
	if( mxl != PETSC_NULL ) {
		*mxl = m;
		if( (sx+m) == M ) { /* last proc */
			*mxl = m - 1;
		}
	}
	if( myl != PETSC_NULL ) {
		*myl = n;
		if( (sy+n) == N ) { /* last proc */
			*myl = n - 1;
		}
	}
	if( mzl != PETSC_NULL ) {
		*mzl = p;
		if( (sz+p) == P ) { /* last proc */
			*myl = p - 1;
		}
	}
    
	/* check */
	/*
     DAGetInfo( da, 0, &M,&N,&P, 0,0,0, 0,0,0,0 );
     
     MPI_Allreduce( &ml, &MM, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD );
     MPI_Allreduce( &nl, &NN, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD );
     MPI_Allreduce( &pl, &PP, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD );
     
     if( (mxl!=PETSC_NULL) && (MM != (M-1) ) ) {		SETERRQ( PETSC_ERR_USER, "Number of elements in x is wrong" );		}
     if( (myl!=PETSC_NULL) && (NN != (N-1) ) ) {		SETERRQ( PETSC_ERR_USER, "Number of elements in y is wrong" );		}
     if( (mzl!=PETSC_NULL) && (PP != (P-1) ) ) {		SETERRQ( PETSC_ERR_USER, "Number of elements in z is wrong" );		}
     */
    
	PetscFunctionReturn(0);
}
/*--------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/

#undef __FUNCT__
#define __FUNCT__ "ConstructRow"
PetscErrorCode ConstructRow(const PetscInt m, const PetscInt istart, const PetscInt jstart, const PetscInt ei, const PetscInt ej, const PetscInt iglobal[], PetscInt row[])
{
    PetscInt i,j;
    
    i = ei-istart;
    j = ej-jstart;
    
    row[0] = iglobal[i+j*m];
    row[1] = iglobal[i+1+j*m];
    row[2] = iglobal[i+1+(j+1)*m];
    row[3] = iglobal[i+(j+1)*m];
    
    PetscFunctionReturn(0);
}
/*---------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/


PetscErrorCode DAGetElementOwnershipRanges( DM da, PetscInt **_lx, PetscInt **_ly )
{
    PetscErrorCode ierr;
	PetscMPIInt rank;
	PetscInt proc_I, proc_J;
	PetscInt cpu_x, cpu_y;
	PetscInt local_mx, local_my;
	Vec vlx, vly;
	PetscInt *LX, *LY, i;
	PetscScalar *_a;
	Vec V_SEQ;
	VecScatter ctx;
    
    
	MPI_Comm_rank( PETSC_COMM_WORLD, &rank );
    
	DMDAGetInfo( da, 0, 0,0,0, &cpu_x,&cpu_y,0, 0,0,0,0,0,0);
    
	proc_J = rank/cpu_x;
	proc_I = rank - cpu_x * proc_J;
    
	PetscMalloc( sizeof(PetscInt)*(size_t)cpu_x, &LX );
	PetscMalloc( sizeof(PetscInt)*(size_t)cpu_y, &LY );
    
	ierr = DAGetLocalElementSize( da, &local_mx,&local_my,PETSC_NULL);                  CHKERRQ(ierr);
	ierr = VecCreate( PETSC_COMM_WORLD, &vlx );                                         CHKERRQ(ierr);
	ierr = VecSetSizes( vlx, PETSC_DECIDE, cpu_x );                                     CHKERRQ(ierr);
	ierr = VecSetFromOptions( vlx );                                                    CHKERRQ(ierr);
    
	ierr = VecCreate( PETSC_COMM_WORLD, &vly );                                         CHKERRQ(ierr);
	ierr = VecSetSizes( vly, PETSC_DECIDE, cpu_y );                                     CHKERRQ(ierr);
	ierr = VecSetFromOptions( vly );                                                    CHKERRQ(ierr);
    
	ierr = VecSetValue( vlx, proc_I, (PetscScalar)(local_mx+1.0e-9), INSERT_VALUES );   CHKERRQ(ierr);
	ierr = VecSetValue( vly, proc_J, (PetscScalar)(local_my+1.0e-9), INSERT_VALUES );   CHKERRQ(ierr);
	ierr = VecAssemblyBegin(vlx);                                                       CHKERRQ(ierr);
    ierr = VecAssemblyEnd(vlx);                                                         CHKERRQ(ierr);
	ierr = VecAssemblyBegin(vly);                                                       CHKERRQ(ierr);
    ierr = VecAssemblyEnd(vly);                                                         CHKERRQ(ierr);
    
    
    
	ierr = VecScatterCreateToAll(vlx,&ctx,&V_SEQ);                                      CHKERRQ(ierr);
	ierr = VecScatterBegin(ctx,vlx,V_SEQ,INSERT_VALUES,SCATTER_FORWARD);                CHKERRQ(ierr);
	ierr = VecScatterEnd(ctx,vlx,V_SEQ,INSERT_VALUES,SCATTER_FORWARD);                  CHKERRQ(ierr);
	ierr = VecGetArray( V_SEQ, &_a );                                                   CHKERRQ(ierr);
	for( i=0; i<cpu_x; i++ ){
		LX[i] = (PetscInt)_a[i];
	}
	ierr = VecRestoreArray(V_SEQ,&_a);                                                  CHKERRQ(ierr);
	ierr = VecScatterDestroy(&ctx);                                                      CHKERRQ(ierr);
	ierr = VecDestroy(&V_SEQ);                                                           CHKERRQ(ierr);
    
	ierr = VecScatterCreateToAll(vly,&ctx,&V_SEQ);                                      CHKERRQ(ierr);
	ierr = VecScatterBegin(ctx,vly,V_SEQ,INSERT_VALUES,SCATTER_FORWARD);                CHKERRQ(ierr);
	ierr = VecScatterEnd(ctx,vly,V_SEQ,INSERT_VALUES,SCATTER_FORWARD);                  CHKERRQ(ierr);
	ierr = VecGetArray( V_SEQ, &_a );                                                   CHKERRQ(ierr);
	for( i=0; i<cpu_y; i++ ){
		LY[i] = (PetscInt)_a[i];
	}
	ierr = VecRestoreArray(V_SEQ,&_a);                                                  CHKERRQ(ierr);
	ierr = VecScatterDestroy(&ctx);                                                      CHKERRQ(ierr);
	ierr = VecDestroy(&V_SEQ);                                                           CHKERRQ(ierr);
    
    
    
	/*
     VecGetArray( vlx, &_a );
     for( i=0; i<proc_I; i++ ){
     LX[i] = (PetscInt)_a[i];
     }
     VecRestoreArray(vlx,&_a);
     VecGetArray( vly, &_a );
     for( i=0; i<proc_J; i++ ){
     LY[i] = (PetscInt)_a[i];
     }
     VecRestoreArray(vlx,&_a);
     */
	*_lx = LX;
	*_ly = LY;
    
	ierr = VecDestroy(&vlx);                                                             CHKERRQ(ierr);
	ierr = VecDestroy(&vly);                                                             CHKERRQ(ierr);
    
	PetscFunctionReturn(0);
}
/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/

#undef __FUNCT__
#define __FUNCT__ "Set_Diffusion_Coefficient"
PetscErrorCode Set_Diffusion_Coefficient( FE_Erosion_UserData *Erosion_UserData)
{
    
    /* Get variables from the structure Erosion_UserData */
    
    PetscErrorCode ierr;
    PetscScalar     lx, ly;
    PetscScalar     rain, k0, c, k, hack;
    PetscScalar     H_mean[4];
    PetscScalar     mean_high;
    PetscScalar     rain_river;
    PetscInt        nelx, nely, fill_lake;
//    PetscInt		mode_river, nbre_river;
    PetscInt        i, j, id, m, n, nstart, mstart;
    Vec             Elevation;
    Vec             Elevation_local;
    Vec             properties, l_properties;
    PetscScalar     *Hvec, *Qvec, **H;
    DM              da_elevation, da_element;
    FE_ElementInfo  **element_info;
    
    
    nelx = Erosion_UserData->nelx;
    nely = Erosion_UserData->nely;
	lx   = Erosion_UserData->lx;
    ly   = Erosion_UserData->ly; 
    rain = Erosion_UserData->rain;
    c    = Erosion_UserData->c;
    k0   = Erosion_UserData->k0;
    hack = Erosion_UserData->hack;
    da_elevation    =  Erosion_UserData->da_elevation;
    da_element      =  Erosion_UserData->da_element;
    Elevation       =  Erosion_UserData->Elevation;
    properties      =  Erosion_UserData->properties;
    fill_lake       =  Erosion_UserData->fill_lake;
//    mode_river      =  Erosion_UserData->mode_river;
//    nbre_river      =  Erosion_UserData->nbre_river;
    rain_river      =  Erosion_UserData->rain_river;
    rain            =  Erosion_UserData->rain;
    
    /* Define Hvec and Qvec, vectors used for the computing the stream flow */
    ierr = PetscMalloc((size_t)(nelx*nely)*sizeof(PetscScalar), &Hvec);            CHKERRQ(ierr);
    ierr = PetscMalloc((size_t)(nelx*nely)*sizeof(PetscScalar), &Qvec);            CHKERRQ(ierr);
    
    
    /* get the elevation per nodes */
    ierr = DMGetLocalVector(da_elevation,&Elevation_local);                                      CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(da_elevation,Elevation,INSERT_VALUES,Elevation_local);           CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(da_elevation,Elevation,INSERT_VALUES,Elevation_local);             CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_elevation,Elevation_local,&H);                                     CHKERRQ(ierr);
    
    
    /* Copy properties vector to local processor */
	ierr = DMGetLocalVector(da_element,&l_properties);                                          CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(da_element,properties,INSERT_VALUES,l_properties);              CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(da_element,properties,INSERT_VALUES,l_properties);                CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_element,l_properties,&element_info);                                           CHKERRQ(ierr);
    
    DMDAGetCorners(da_element,&mstart,&nstart,PETSC_NULL,&m,&n,PETSC_NULL);
    for (j=nstart; j<nstart+n; j++) {
        for (i=mstart; i<mstart+m; i++) {
            
            id = i+j*nelx;
            ierr = GetAverageNodes(H_mean, i, j, H, &mean_high);
            element_info[j][i].H_el = mean_high;
            Hvec[id] = mean_high;
            
        }
	}
    
    /* send elevation data back to global vector */
    ierr = DMDAVecRestoreArray(da_elevation,Elevation_local,&H);                                CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(da_elevation,&Elevation_local);                                 CHKERRQ(ierr);
    
    
    
    ierr = DMDAVecRestoreArray(da_element,l_properties,&element_info);                           CHKERRQ(ierr);
    
    
    /* Call the function stream_flow_D8 */

     ierr = stream_flow(nely, nelx, rain, rain_river, Qvec, lx, ly, Hvec, Erosion_UserData->location_river, Erosion_UserData->nbre_river, fill_lake);      CHKERRQ(ierr);
    
    
    /*================== Define the coefficients for diffusion ==========================*/
    
    /* Copy properties vector to local processor */
	ierr = DMDAVecGetArray(da_element,l_properties,&element_info);                              CHKERRQ(ierr);
    
    DMDAGetCorners(da_element,&mstart,&nstart,PETSC_NULL,&m,&n,PETSC_NULL);
    for (j=nstart; j<nstart+n; j++) {
        for (i=mstart; i<mstart+m; i++) {
            
            id = i+j*nelx;
            
            element_info[j][i].Q_el = Qvec[id];
            
            //            if (Qvec[id]<2.0*threshold){
            //                Qvec[id]=10.0*threshold;
            //            }
            
            k = k0 + c*pow(Qvec[id],hack);
            
            element_info[j][i].k_el = k;
            
            //  printf("H_element[%d][%d]: %.10e, k_el: %.20e, Q_element[%d][%d]: %.10e, s_el: %.10e\n", j,i, element_info[j][i].H_el, k, j, i, element_info[j][i].Q_el,s);
        }
    }
    
    
    /* Send properties data back to global vector */
    
    
    ierr = DMDAVecRestoreArray(da_element, l_properties, &element_info);                              CHKERRQ(ierr);
    ierr = DMLocalToGlobalBegin(da_element,l_properties,INSERT_VALUES,properties);                    CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(da_element,l_properties,INSERT_VALUES,properties);                      CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(da_element,&l_properties);                                            CHKERRQ(ierr);
    
    PetscFree(Qvec);
    PetscFree(Hvec);
    
    
    PetscFunctionReturn(0);
}

/*==================================FUNCTIONS FOR THE STREAMFLOW AND OTHERS==========================================*/

#undef __FUNCT__
#define __FUNCT__ "GetAverageNodes"
PetscErrorCode GetAverageNodes(PetscScalar table[], PetscInt i, PetscInt j, PetscScalar **matrix, PetscScalar *mean){
    
    PetscScalar average;
    table[0] = matrix[j][i];
    table[1] = matrix[j][i+1];
    table[2] = matrix[j+1][i];
    table[3] = matrix[j+1][i+1];
    
    average = 0.25*(table[0]+table[1]+table[2]+table[3]);
    *mean = average;
    
    // printf("matrix[j][i], matrix[j][i], matrix[j][i+1], matrix[j+1][i+1], matrix[j+1][i]: %.10f, %.10f, %.10f, %.10f\n", matrix[j][i],matrix[j][i+1],matrix[j+1][i+1],matrix[j+1][i]);
    PetscFunctionReturn(0);
}




/*--------------------------------------------------------------------------------------------------------------------
 Compare and sort Data. the data stored: elevation (H), number, and lowest neighbor of each nodes of the grid*/

#undef __FUNCT__
#define __FUNCT__ "NodeDataComparison"
PetscErrorCode NodeDataComparison(const void *a, const void *b){
    
    const NodeData *nd_a;
    const NodeData *nd_b;
    
    nd_a = (const NodeData*)a;
    nd_b = (const NodeData*)b;
    
    if ( nd_a->H < nd_b->H) {
        return 1;
    } else if ( nd_a->H > nd_b->H) {
        return -1;
    } else {
        return 0;
    }
    PetscFunctionReturn(0);
    
}

#undef __FUNCT__
#define __FUNCT__ "NodeDataSort"
PetscErrorCode NodeDataSort(PetscInt L,PetscScalar H[],PetscInt low[],PetscInt indices[]){
    
    NodeData *list;
    PetscInt i;
    PetscErrorCode ierr;
    
    ierr = PetscMalloc((size_t)L*sizeof( NodeData),&list); CHKERRQ(ierr);
    
    memset( list, 0, sizeof(NodeData)*(size_t)L );
    
    for (i=0; i<L; i++) {
        list[i].H = H[i];
        list[i].low = low[i];
        list[i].num = indices[i];
        //  printf("input : %.4d H,low,num = %1.4f,%.4d,%.4d\n",i,H[i],low[i],indices[i]);
        
    }
    
    
    // sort
    qsort( &list[0], (size_t)L, sizeof(NodeData), &NodeDataComparison );
    
    // permute original
    for (i=0; i<L; i++) {
        H[i]      = list[i].H;
        low[i]    = list[i].low;
        indices[i] =list[i].num;
        //printf("output: %.4d H,low,num = %1.4f,%.4d,%.4d\n",i,H[i],low[i],indices[i]);
    }
    
    // clean up
    PetscFree(list);
    PetscFunctionReturn(0);
    
}
/*---------------------------------------------------------------------------------------------------------------------*/



/* find the minimum in an array*/

#undef __FUNCT__
#define __FUNCT__ "lamem_min"
PetscErrorCode lamem_min(PetscScalar table[], PetscInt sizetable, PetscInt *location_minvalue, PetscScalar *minvalue){
    
    PetscScalar _min = table[0];
    PetscInt i, c;
    
    c = 0;
    for (i = 0; i < sizetable; i++){
        if (table[i]<_min){
            _min = table[i];
            c = i;
        }
    }
    
    *location_minvalue = c;
    *minvalue = _min;
    
    PetscFunctionReturn(0);
}
/*---------------------------------------------------------------------------------------------------------------------*/

/* Choose between ceil and floor to locate the river */
#undef __FUNCT__
#define __FUNCT__ "choise_int"
PetscErrorCode choise_int(PetscScalar arrondi_sup, PetscScalar arrondi_inf, PetscScalar result, PetscInt *result_int)
{
    PetscScalar diff_ceil, diff_floor,min;
    PetscScalar diff[2];
    PetscInt  c;
    
    
    diff_ceil  = PetscAbsScalar(arrondi_sup-result);
    diff_floor = PetscAbsScalar(result-arrondi_inf);
    
    diff[0] = diff_ceil;
    diff[1] = diff_floor;
    
    min = diff[0];
    if (diff[1]<diff[0]){
        min = diff[1]; 
    }
    if (min==diff[0]){
        c = (PetscInt)arrondi_sup;
    }
    else {
        c= (PetscInt)arrondi_inf;
    }
    
    *result_int = c; 
    
    PetscFunctionReturn(0); 
}
/*---------------------------------------------------------------------------------------------------------------------*/



/* Check the lowest neighbor among nine nodes*/

#undef __FUNCT__
#define __FUNCT__ "low_neigh"
PetscErrorCode low_neigh(PetscScalar H[], PetscInt low [], PetscInt M, PetscInt N){
    PetscInt i,j, id;
    
    for (j =0; j<N; j++){
        for (i =0; i<M; i++){
            
            id = i + j*M;
            
            if (i==0 && j==0)
            {
                PetscScalar Hloc[] = {H[M], H[1], H[M+1], H[0]};
                PetscInt  numloc[] = {M, 1, M+1, 0};
                PetscInt location;
                PetscScalar minvalue;
                
                lamem_min(Hloc,4, &location,&minvalue);
                low[id] = numloc[location];
            }
            
            else if (i==M-1 && j==0)
            {
                PetscScalar Hloc[] = {H[2*M-1], H[M-2], H[2*M-2], H[M-1]};
                PetscInt numloc[]  = {2*M-1, M-2, 2*M-2, M-1};
                PetscInt location;
                PetscScalar minvalue;
                
                lamem_min(Hloc,4, &location,&minvalue);
                low[id] = numloc[location];
            }
            
            else if (i==0 && j==N-1)
            {
                PetscScalar Hloc[] = {H[id-M], H[id+1], H[id-M+1], H[id]};
                PetscInt numloc[]  = {id-M, id+1, id-M+1, id};
                PetscInt location;
                PetscScalar minvalue;
                
                lamem_min(Hloc,4, &location,&minvalue);
                low[id] = numloc[location];
            }
            
            else if (i==M-1 && j==N-1)
            {
                PetscScalar Hloc[] = {H[id-M], H[id-1], H[id-M-1], H[id]};
                PetscInt numloc[]  = {id-M, id-1, id-M-1, id};
                PetscInt location;
                PetscScalar minvalue;
                
                lamem_min(Hloc,4, &location,&minvalue);
                low[id] = numloc[location];
                
            }
            
            else if (i==0)
            {
                PetscScalar Hloc[] = {H[id-M], H[id+M], H[id+1], H[id+M+1], H[id-M+1], H[id]};
                PetscInt numloc[] = {id-M, id+M, id+1, id+M+1, id-M+1, id};
                PetscInt location;
                PetscScalar minvalue;
                
                lamem_min(Hloc,6, &location,&minvalue);
                low[id] = numloc[location];
            }
            
            else if (i==M-1)
            {
                PetscScalar Hloc[] = {H[id-M], H[id+M], H[id-1], H[id-1-M], H[id-1+M], H[id]};
                PetscInt numloc[] = {id-M, id+M, id-1, id-1-M, id-1+M, id};
                PetscInt location;
                PetscScalar minvalue;
                
                lamem_min(Hloc,6, &location,&minvalue);
                low[id] = numloc[location];
                
            }
            
            else if (j==0)
            {
                PetscScalar Hloc[] = {H[id+M], H[id+1], H[id-1], H[id+M+1], H[id+M-1], H[id]};
                PetscInt numloc[] = {id+M, id+1, id-1, id+M+1, id+M-1, id};
                PetscInt location;
                PetscScalar minvalue;
                
                lamem_min(Hloc,6, &location,&minvalue);
                low[id] = numloc[location];
            }
            
            else if (j==N-1)
            {
                PetscScalar Hloc[] = {H[id-M], H[id+1], H[id-1], H[id-M-1], H[id-M+1], H[id]};
                PetscInt numloc[] = {id-M, id+1, id-1, id-M-1, id-M+1, id};
                PetscInt location;
                PetscScalar minvalue;
                
                lamem_min(Hloc,6, &location,&minvalue);
                low[id] = numloc[location];
                
            }
            
            else
            {
                PetscScalar Hloc[] ={H[id-M], H[id+M], H[id-1], H[id+1], H[id+1+M], H[id+1-M], H[id-1+M], H[id-1-M], H[id]};
                PetscInt numloc[] = {id-M, id+M, id-1, id+1, id+1+M, id+1-M, id-1+M, id-1-M, id};
                PetscInt location;
                PetscScalar minvalue;
                
                lamem_min(Hloc,9, &location,&minvalue);
                low[id] = numloc[location];
            }
            
        }
    }
    
    
    PetscFunctionReturn(0);
}
/*---------------------------------------------------------------------------------------------------------------*/

#undef __FUNCT__
#define __FUNCT__ "NodeNeighbour_LowestHeight"
PetscErrorCode NodeNeighbour_LowestHeight(PetscInt Mx,PetscInt My,PetscScalar *H,PetscInt low[])
{
    PetscInt i,j,ni,nj,II,JJ,n;
    PetscInt min_n_idx,min_idx;
    PetscScalar min_H;
    PetscInt n_neighbours,n_list[9];
    
    PetscFunctionBegin;
    
    for (j=0; j<My; j++) {
        for (i=0; i<Mx; i++) {
            PetscInt idx = i + j * Mx;
            
            n_neighbours = 0;
            for (nj=-1; nj<2; nj++) {
                for (ni=-1; ni<2; ni++) {
                    
                    II = i+ni;
                    JJ = j+nj;
                    
                    if ( (II<0) || (II>=Mx) ) { continue; }
                    if ( (JJ<0) || (JJ>=My) ) { continue; }
                    
                    n_list[n_neighbours] = II + JJ * Mx;
                    
                    n_neighbours++;
                }
            }                      
            /* sanity check */
            if (n_neighbours == 0) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"num neighbours = 0");
            }
            if (n_neighbours > 9) {
                SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"num neighbours > 9");
            }
            
            min_n_idx = 0; /* index of neighbour which is lowest (0..8) */
            min_H   = H[ n_list[0] ];
            for (n=1; n<n_neighbours; n++) {
                if (H[ n_list[n] ] < min_H) {
                    min_H = H[ n_list[n] ];
                    min_n_idx = n;
                }
            }
            min_idx = n_list[min_n_idx];
            
            low[idx] = min_idx;
            
        }
    }
    
    PetscFunctionReturn(0);
}
/*---------------------------------------------------------------------------------------------------------------*/


/* check the lowest neighbourg according to gradient */


#undef __FUNCT__
#define __FUNCT__ "low_neigh_grad"
PetscErrorCode low_neigh_grad(PetscScalar H[], PetscInt low [], PetscInt M, PetscInt N, PetscScalar lx, PetscScalar ly){
    
    PetscInt i,j, id;
    PetscScalar dx, dy, dX, dY, hyp;
    
    
    dx = lx/M;
    dy = ly/N;
    dX = pow(dx,2.0);
    dY = pow(dy,2.0);
    hyp = sqrt(dX+dY);
    
    for (j =0; j<N; j++){
        for (i =0; i<M; i++){
            
            id = i + j*M;
            
            if (i==0 && j==0)
            {
                PetscScalar Hloc[] = {(H[M]-H[0])/dy, (H[1]-H[0])/dx, (H[M+1]-H[0])/hyp, (H[0]-H[0])};
                PetscInt  numloc[] = {M, 1, M+1, 0};
                PetscInt location;
                PetscScalar minvalue;
                
                lamem_min(Hloc,4, &location,&minvalue);
                low[id] = numloc[location];
            }
            
            else if (i==M-1 && j==0)
            {
                PetscScalar Hloc[] = {(H[2*M-1]-H[M-1])/dy, (H[M-2]-H[M-1])/dx, (H[2*M-2]-H[M-1])/hyp, (H[M-1]-H[M-1])};
                PetscInt numloc[]  = {2*M-1, M-2, 2*M-2, M-1};
                PetscInt location;
                PetscScalar minvalue;
                
                lamem_min(Hloc,4, &location,&minvalue);
                low[id] = numloc[location];
            }
            
            else if (i==0 && j==N-1)
            {
                PetscScalar Hloc[] = {(H[id-M]-H[id])/dy, (H[id+1]-H[id])/dx, (H[id-M+1]-H[id])/hyp, (H[id]-H[id])};
                PetscInt numloc[]  = {id-M, id+1, id-M+1, id};
                PetscInt location;
                PetscScalar minvalue;
                
                lamem_min(Hloc,4, &location,&minvalue);
                low[id] = numloc[location];
            }
            
            else if (i==M-1 && j==N-1)
            {
                PetscScalar Hloc[] = {(H[id-M]-H[id])/dy, (H[id-1]-H[id])/dx, (H[id-M-1]-H[id])/hyp, (H[id]-H[id])};
                PetscInt numloc[]  = {id-M, id-1, id-M-1, id};
                PetscInt location;
                PetscScalar minvalue;
                
                lamem_min(Hloc,4, &location,&minvalue);
                low[id] = numloc[location];
                
            }
            
            else if (i==0)
            {
                PetscScalar Hloc[] = {(H[id-M]-H[id])/dy, (H[id+M]-H[id])/dy, (H[id+1]-H[id])/dy, (H[id+M+1]-H[id])/hyp, (H[id-M+1]-H[id])/hyp, (H[id]-H[id])};
                PetscInt numloc[] = {id-M, id+M, id+1, id+M+1, id-M+1, id};
                PetscInt location;
                PetscScalar minvalue;
                
                lamem_min(Hloc,6, &location,&minvalue);
                low[id] = numloc[location];
            }
            
            else if (i==M-1)
            {
                PetscScalar Hloc[] = {(H[id-M]-H[id])/dy, (H[id+M]-H[id])/dy, (H[id-1]-H[id])/dx, (H[id-1-M]-H[id])/hyp, (H[id-1+M]-H[id])/hyp, (H[id]-H[id])};
                PetscInt numloc[] = {id-M, id+M, id-1, id-1-M, id-1+M, id};
                PetscInt location;
                PetscScalar minvalue;
                
                lamem_min(Hloc,6, &location,&minvalue);
                low[id] = numloc[location];
                
            }
            
            else if (j==0)
            {
                PetscScalar Hloc[] = {(H[id+M]-H[id])/dy, (H[id+1]-H[id])/dx, (H[id-1]-H[id])/dx, (H[id+M+1]-H[id])/hyp, (H[id+M-1]-H[id])/hyp, (H[id]-H[id])};
                PetscInt numloc[] = {id+M, id+1, id-1, id+M+1, id+M-1, id};
                PetscInt location;
                PetscScalar minvalue;
                
                lamem_min(Hloc,6, &location,&minvalue);
                low[id] = numloc[location];
            }
            
            else if (j==N-1)
            {
                PetscScalar Hloc[] = {(H[id-M]-H[id])/dy, (H[id+1]-H[id])/dx, (H[id-1]-H[id])/dx, (H[id-M-1]-H[id])/hyp, (H[id-M+1]-H[id])/hyp, (H[id]-H[id])};
                PetscInt numloc[] = {id-M, id+1, id-1, id-M-1, id-M+1, id};
                PetscInt location;
                PetscScalar minvalue;
                
                lamem_min(Hloc,6, &location,&minvalue);
                low[id] = numloc[location];
                
            }
            
            else
            {
                PetscScalar Hloc[] ={(H[id-M]-H[id])/dy, (H[id+M]-H[id])/dy, (H[id-1]-H[id])/dx, (H[id+1]-H[id])/dx, (H[id+1+M]-H[id])/hyp, (H[id+1-M]-H[id])/hyp, (H[id-1+M]-H[id])/hyp, (H[id-1-M]-H[id])/hyp, (H[id]-H[id])};
                PetscInt numloc[] = {id-M, id+M, id-1, id+1, id+1+M, id+1-M, id-1+M, id-1-M, id};
                PetscInt location;
                PetscScalar minvalue;
                
                lamem_min(Hloc,9, &location,&minvalue);
                low[id] = numloc[location];
            }
            
        }
    }
    
    
    PetscFunctionReturn(0);
}







/*---------------------------------------------------------------------------------
 Compute the water surface discharge (Streamflow = Q[m2/s]) using the D8 algorithm
 D8 algorithm modified from O'Callaghan & Mark (1984)
 ----------------------------------------------------------------------------------*/

#undef __FUNCT__
#define __FUNCT__ "stream_flow_D8"
PetscErrorCode stream_flow_D8(PetscInt N, PetscInt M, PetscScalar rain, PetscScalar Qvec[], PetscScalar ly, PetscScalar Hvec[], PetscInt fill_lake){
    
    PetscErrorCode ierr;
    PetscInt i, j, id;
    
    
    /*-----------------Initialization of the q matrix---------------------------*/
    PetscScalar *length_stream;
    PetscScalar dy;
    
    ierr = PetscMalloc((size_t)(M*N)*sizeof(PetscScalar),&length_stream);  CHKERRQ(ierr);
    
    dy = ly/N;
    
    for(j=0; j<N; j++){
        for (i=0; i<M; i++){
            
            id = i + j*M;
            length_stream[id]=dy;
        }
    }
    
    /*-----------------Indices matrix from 0 to N*M-1----------------------------*/
    for (i=0; i<M*N; i++){
        
        Qvec[i]= rain * length_stream[i];
        
    }
    
    /*-----------------Indices matrix from 0 to N*M-1----------------------------*/
    PetscInt *num;
    ierr = PetscMalloc((size_t)(M*N)*sizeof(PetscInt),&num);   CHKERRQ(ierr);
    
    for (j= 0; j<N; j++){
        for (i = 0; i<M; i++){
            
            id = i + j*M;
            num[id] = id;
        }
    }
    
    /*------------------Check for the lowest neighbour-----------------------------*/
    PetscInt *low;
    
    ierr = PetscMalloc((size_t)(M*N)*sizeof(PetscInt),&low);                    CHKERRQ(ierr);
    // ierr = low_neigh(Hvec, low, M, N);                                       CHKERRQ(ierr);
    ierr = NodeNeighbour_LowestHeight(M,N,Hvec,low);                            CHKERRQ(ierr);
    

    /*------------------Sort the cell from highest to lowest------------------------*/
    ierr = NodeDataSort(M*N, Hvec, low, num);                                   CHKERRQ(ierr);
    
    
    /*------------------Define q for each cell--------------------------------------*/    
    if (fill_lake==1){ 
        for (j=0; j<N; j++){
            for (i=0; i<M; i++){
                id = i + j*M;
                if (Hvec[num[id]]==Hvec[low[id]] && i!=0 && i!=M-1 && j!=0 && j!=N-1){
                    Qvec[low[id]] = 0.0;     
                }
                else {
                    Qvec[low[id]] = Qvec[low[id]] + Qvec[num[id]];   
                }
            }
        }
    }
    
    else {
        for (i=0; i<M*N; i++){
            if (Hvec[num[i]]!=Hvec[low[i]]){
                Qvec[low[i]] = Qvec[low[i]] + Qvec[num[i]];
            }
        }
    }
    

    /*------------------Free memory-------------------------------------------------*/
    ierr = PetscFree(low);                                                      CHKERRQ(ierr);
    ierr = PetscFree(num);                                                      CHKERRQ(ierr);
    ierr = PetscFree(length_stream);                                            CHKERRQ(ierr);


PetscFunctionReturn(0);
}
/*---------------------------------------------------------------------------------------------*/


#undef __FUNCT__
#define __FUNCT__ "stream_flow_D8_location_river"
PetscErrorCode stream_flow_D8_location_river(PetscInt N, PetscInt M, PetscScalar rain, PetscScalar Qvec[], PetscScalar lx,PetscScalar ly, PetscScalar Hvec[], PetscScalar location_river[], PetscInt nbre_river)
{
    
    PetscErrorCode ierr;
    PetscInt i, j, id;
    
    
    /*----------------Length associated to each node------------------*/
    PetscScalar *length_stream;
    PetscScalar dy;
    
    ierr = PetscMalloc((size_t)(M*N)*sizeof(PetscScalar),&length_stream);                 CHKERRQ(ierr);
    
    dy = ly/N;
    
    for(j=0; j<N; j++){
        for (i=0; i<M; i++){
            id = i + j*M;
            length_stream[id]=dy;
        }
    }
    
    /*-----------------Initialization of the rain matrix----------------*/
    PetscScalar *rain_vec; 
    PetscInt nodes_x[nbre_river], indices_river[nbre_river];
    PetscScalar result; 
    PetscScalar arrondi_inf, arrondi_sup;
    PetscInt	result_int;
    
    ierr = PetscMalloc((size_t)(M*N)*sizeof(PetscScalar),&rain_vec);                      CHKERRQ(ierr);
    
    for (i=0; i<nbre_river; i++){
        result  = location_river[i]*(M/lx);
      //  printf("location_river[%d]: %f, nbre_river:%d\n",i,location_river[i], nbre_river);
        
        arrondi_inf = floor(result);
        arrondi_sup = ceil(result);  
        
        ierr = choise_int(arrondi_sup, arrondi_inf, result, &result_int);
        nodes_x[i] = result_int; 
        indices_river[i] = nodes_x[i] + (N-1)*M;
    }
    
    for (j= 0; j<N; j++){
        for (i = 0; i<M; i++){ 
            id = i+j*M;
            rain_vec[id]=0.0;
        }
    }
    
    for (i=0; i<nbre_river; i++){
        id = indices_river[i]; 
        rain_vec[id]=rain;
    }
    
    
    /*-----------------Initialization of the q matrix---------------------------*/
    for (i=0; i<M*N; i++){
        Qvec[i]= rain_vec[i]*length_stream[i];
    }
    
    
    /*-----------------Indices matrix from 0 to N*M-1----------------------------*/
    PetscInt *num; 
    ierr = PetscMalloc((size_t)(M*N)*sizeof(PetscInt),&num);                              CHKERRQ(ierr);
    
    for (j= 0; j<N; j++){
        for (i = 0; i<M; i++){  
            id = i + j*M;
            num[id] = id;
        }
    }
    
    
    /*------------------Check for the lowest neighbour-----------------------------*/
    PetscInt *low;
    ierr = PetscMalloc((size_t)(M*N)*sizeof(PetscInt),&low);                              CHKERRQ(ierr);
    ierr = low_neigh(Hvec, low, M, N);                                          CHKERRQ(ierr);
    
    
    /*------------------Sort the cell from highest to lowest------------------------*/
    ierr = NodeDataSort(M*N, Hvec, low, num);                                   CHKERRQ(ierr);
    
    
    /*------------------Define q for each cell--------------------------------------*/
    PetscScalar cell; 
    cell = rain*dy; 
    
    for (i=0; i<M*N; i++){  
        if (Qvec[num[i]]!=0.0){
            Qvec[low[i]] = Qvec[low[i]] + Qvec[num[i]]+cell;
        }
    }
    
    
    /*------------------Free memory-------------------------------------------------*/
    ierr = PetscFree(low);CHKERRQ(ierr);
    ierr = PetscFree(num);CHKERRQ(ierr);
    ierr = PetscFree(length_stream); CHKERRQ(ierr); 
    ierr = PetscFree(rain_vec); CHKERRQ(ierr); 
    
    
    PetscFunctionReturn(0);    
}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/


#undef __FUNCT__
#define __FUNCT__ "stream_flow"
PetscErrorCode stream_flow(PetscInt N, PetscInt M, PetscScalar rain, PetscScalar rain_river, PetscScalar Qvec[], PetscScalar lx,PetscScalar ly, PetscScalar Hvec[], PetscScalar location_river[], PetscInt nbre_river, PetscInt fill_lake)
{
    
    PetscErrorCode ierr;
    PetscInt i, j;
    PetscInt id;
    
    
    /*----------------Length associated to each node------------------*/
    PetscScalar *length_stream;
    PetscScalar dy;
    
    ierr = PetscMalloc((size_t)(M*N)*sizeof(PetscScalar),&length_stream);                 CHKERRQ(ierr);
    
    dy = ly/N;
    
    for(j=0; j<N; j++){
        for (i=0; i<M; i++){
            id = i + j*M;
            length_stream[id]=dy;
        }
    }
    
    /*-----------------Initialization of the rain matrix----------------*/
    PetscScalar *rain_vec; 
    PetscInt nodes_x[nbre_river], indices_river[nbre_river];
    PetscScalar result; 
    PetscScalar arrondi_inf, arrondi_sup;
    PetscInt    result_int;
    
    ierr = PetscMalloc((size_t)(M*N)*sizeof(PetscScalar),&rain_vec);                      CHKERRQ(ierr);
    
    for (i=0; i<nbre_river; i++){
        result  = location_river[i]*(M/lx);
        arrondi_inf = floor(result);
        arrondi_sup = ceil(result);
        
        ierr = choise_int(arrondi_sup, arrondi_inf, result, &result_int);
        nodes_x[i] = result_int; 
        indices_river[i] = nodes_x[i] + (N-1)*M;
    }
    
    for (j= 0; j<N; j++){
        for (i = 0; i<M; i++){ 
            id = i+j*M;
            rain_vec[id]=rain;
        }
    }
    
    for (i=0; i<nbre_river; i++){
        id = indices_river[i]; 
        rain_vec[id]=rain+rain_river;
    }
    
    
    /*-----------------Initialization of the q matrix---------------------------*/
    for (i=0; i<M*N; i++){
        Qvec[i]= rain_vec[i]*length_stream[i];
    }
    
    
    /*-----------------Indices matrix from 0 to N*M-1----------------------------*/
    PetscInt *num;
    ierr = PetscMalloc((size_t)(M*N)*sizeof(PetscInt),&num);                              CHKERRQ(ierr);
    
    for (j= 0; j<N; j++){
        for (i = 0; i<M; i++){  
            id = i + j*M;
            num[id] = id;
        }
    }
    
    
    /*------------------Check for the lowest neighbour-----------------------------*/
    PetscInt *low;
    ierr = PetscMalloc((size_t)(M*N)*sizeof(PetscInt),&low);                              CHKERRQ(ierr);
    ierr = low_neigh(Hvec, low, M, N);                                          CHKERRQ(ierr);
    
    
    /*------------------Sort the cell from highest to lowest------------------------*/
    ierr = NodeDataSort(M*N, Hvec, low, num);                                   CHKERRQ(ierr);
    
    
    /*------------------Define q for each cell--------------------------------------*/
    PetscScalar cell, background; 
    cell = rain_river*dy; 
    background = rain*dy;
    
    printf("fill_lake = %d",fill_lake);
    if(rain_river!=0.0){
        for (i=0; i<M*N; i++){  
            if (Qvec[num[i]]>background){
                Qvec[low[i]] = Qvec[low[i]] + Qvec[num[i]]+cell;
            }
        }
    }
    
    if (rain!=0.0){
        if (fill_lake==1){ 
            for (j=0; j<N; j++){
                for (i=0; i<M; i++){
                    id = i + j*M;
                    if (Hvec[num[id]]==Hvec[low[id]] && i!=0 && i!=M-1 && j!=0 && j!=N-1){
                        Qvec[low[id]] = 0.0;     
                    }
                    else {
                        Qvec[low[id]] = Qvec[low[id]] + Qvec[num[id]];   
                    }
                }
            }
        }
        else {
            for (i=0; i<M*N; i++){
                if (Hvec[num[i]]!=Hvec[low[i]]){
                    Qvec[low[i]] = Qvec[low[i]] + Qvec[num[i]];
                }
            }
        }
    }
    
    
    /*------------------Free memory-------------------------------------------------*/
    ierr = PetscFree(low);CHKERRQ(ierr);
    ierr = PetscFree(num);CHKERRQ(ierr);
    ierr = PetscFree(length_stream); CHKERRQ(ierr); 
    ierr = PetscFree(rain_vec); CHKERRQ(ierr); 
    
    
    PetscFunctionReturn(0);    
}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/



/*==========================================================================================================*/
/* Save initial erosion surface that contains the random noise */
#undef __FUNCT__
#define __FUNCT__ "SaveInitialErosionSurface"
PetscErrorCode SaveInitialErosionSurface( UserContext *user, const char *FileName)
{
    PetscMPIInt rank, size;
    PetscErrorCode 		ierr;
    char				SaveFileName[PETSC_MAX_PATH_LEN];
    PetscViewer			view_out;
    
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    
    if (rank == 0){
        
        /* Create output directory */
        mkdir("InitialErosionSurface", S_IRWXU);
        
        sprintf(SaveFileName,"./InitialErosionSurface/%s.%lld.out",FileName,(LLD)rank);  // construct the filename
        
        /* Write the actual file */
        // Save erosion surface information
        ierr = PetscViewerCreate(PETSC_COMM_SELF,&view_out);				CHKERRQ(ierr);
        ierr = PetscViewerSetType(view_out,PETSCVIEWERBINARY);				CHKERRQ(ierr);
        ierr = PetscViewerFileSetMode(view_out, FILE_MODE_WRITE);			CHKERRQ(ierr);
        ierr = PetscViewerFileSetName(view_out,SaveFileName);		CHKERRQ(ierr);
        ierr = VecView(user->ErosionParameters.FE_ErosionCode.ErosionSurface, view_out); 	CHKERRQ(ierr);		// erosion surface
        
        ierr = PetscViewerDestroy(&view_out); CHKERRQ(ierr);
        
        PetscPrintf(PETSC_COMM_WORLD,"  Saved Initial Erosion Surface to file %s \n",SaveFileName);
        
    }
    
    PetscFunctionReturn(0);
}
/*==========================================================================================================*/

/* Load initial erosion surface that contains the random noise */
#undef __FUNCT__
#define __FUNCT__ "LoadInitialErosionSurface"
PetscErrorCode LoadInitialErosionSurface(UserContext *user)
{
    PetscMPIInt 		rank, size;
    PetscInt 			Nx,Ny,intVecEntries;
    PetscErrorCode 		ierr;
    char				LoadFileName[PETSC_MAX_PATH_LEN];
    PetscViewer			view_in1;
    Vec					tempErosionSurface;
    
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    
    if (rank == 0) {
        
        sprintf(LoadFileName  ,"%s.out",	"./InitialErosionSurface/InitialErosionSurface.0");
        
        /* Load vector that contains erosion surface from the specified folder */
        ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,LoadFileName,FILE_MODE_READ, &view_in1); CHKERRQ(ierr);
        ierr = VecCreate(PETSC_COMM_SELF,&tempErosionSurface);CHKERRQ(ierr);
        ierr = VecLoad(tempErosionSurface,view_in1); 		CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&view_in1); 						CHKERRQ(ierr);
        
        /* Get some info about the initialized DA_FE_ErosionCode array and about the loaded vector*/
        ierr = DMDAGetInfo(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode,PETSC_NULL,&Nx,&Ny,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);
        ierr = VecGetSize(tempErosionSurface, &intVecEntries);
        
        /*check if the loaded InitialErosionSurface can be fit to the Initialized DA_FE_ErosionCode array*/
        if((Nx*Ny)==intVecEntries){
            ierr = VecCopy(tempErosionSurface, user->ErosionParameters.FE_ErosionCode.ErosionSurface); 		CHKERRQ(ierr);
            ierr = VecDestroy(&tempErosionSurface); 		CHKERRQ(ierr);
            /*Set coordinates of the Erosion_Code DM*/
            ierr = DMDASetUniformCoordinates(user->ErosionParameters.FE_ErosionCode.DA_FE_ErosionCode,user->x_left,user->x_left+user->W,user->y_front,user->y_front + user->L,0,0); CHKERRQ(ierr);
            
            PetscPrintf(PETSC_COMM_WORLD,"  Loaded initial erosion surface from file %s \n", LoadFileName);
        }
        else{
            PetscPrintf(PETSC_COMM_WORLD,"  Warning!!! Size of InitialErosionSurface and the Initialized DA_FE_ErosionCode do not fit.The initialized one will be used \n ");
        }
    }
    
    PetscFunctionReturn(0);
}







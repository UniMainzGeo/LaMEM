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

$Id:$
$Date:$

This version is compatible with PETSc version 3.3-p2
The code cannot be distributed, used for commercial purposes or handed on,
without the explicit agreement of Boris Kaus.
*/

#include "LaMEM.h"

/* declare prototypes for constructors here
extern PetscErrorCode LaMEMVelPressureDACreate_Q1P0( LaMEMVelPressureDA C );
extern PetscErrorCode LaMEMVelPressureDACreate_Q2PM1( LaMEMVelPressureDA C );
extern PetscErrorCode LaMEMVelPressureDACreate_Q1Q1( LaMEMVelPressureDA C );
extern PetscErrorCode LaMEMVelPressureDACreate_FDSTAG( LaMEMVelPressureDA C );
*/

typedef PetscErrorCode (*DAVP_Constructor)( LaMEMVelPressureDA);


DAVP_Constructor DAVP_constructor_list[] = {
		LaMEMVelPressureDACreate_Q2PM1,
		LaMEMVelPressureDACreate_Q2PM1,
		LaMEMVelPressureDACreate_Q1P0,
		LaMEMVelPressureDACreate_Q1Q1,
		LaMEMVelPressureDACreate_FDSTAG,
		0 };

const char *DAVPElementType_NAME[] = {
	"Q2Pm1_local",
	"Q2Pm1_global",
	"Q1P0",
	"Q1Q1",
	"FDSTAG",
	0 };

const char *DAVPElementType_DESCRIPTION[] = {
	"Triquadratic velocity basis function (xi,eta,zeta) + discontinuous, linear pressure basis function (xi,eta,zeta)",
	"Triquadratic velocity basis function (xi,eta,zeta) + discontinuous, linear pressure basis function (x,y,z)",
	"Trilinear velocity basis function (xi,eta,zeta) + discontinuous, constant pressure basis function (xi,eta,zeta)",
	"Trilinear velocity basis function (xi,eta,zeta) + Trilinear pressure basis function (xi,eta,zeta)",
	"Finite difference formulation with a staggered velocity and pressure discretization",
	0 };


/*
Change the mesh from the command using the argument
  -vpt_element
*/
#undef __FUNCT__
#define __FUNCT__ "LaMEMVelPressureDACreate"
PetscErrorCode LaMEMVelPressureDACreate( const DAVPElementType type, LaMEMVelPressureDA *_C )
{
	PetscErrorCode ierr;
	DAVP_Constructor fp;
	LaMEMVelPressureDA C;
	PetscBool 	flg;
	char element_name[PETSC_MAX_PATH_LEN];
	DAVPElementType my_type;

	my_type = type;

	/*
	Check command line args for an element type
	If we find a valid option, we will override that specified in the argument, "type".
	*/
	flg = PETSC_FALSE;
	PetscOptionsGetString( PETSC_NULL, "-vpt_element", element_name, PETSC_MAX_PATH_LEN-1, &flg );
	if( flg == PETSC_TRUE ) {
		PetscInt k, type_k;
		PetscBool found;
		const char *type_name;

		// check type is valid //
		found = PETSC_FALSE;
		type_name = DAVPElementType_NAME[0];
		k = 0;
		while( type_name != 0 ) {
			if( strcmp( type_name, element_name ) == 0 ) {
				found = PETSC_TRUE;
				type_k = k;
				break;
			}

			type_name = DAVPElementType_NAME[k+1];
			k++;
		}


		// debug
		if( found == PETSC_FALSE ) {
			PetscPrintf( PETSC_COMM_WORLD, "LaMEMVelPressureDACreate: Unrecognised element type set on command line. \n");
			PetscPrintf( PETSC_COMM_WORLD, "Choose one of { Q2Pm1_local, Q2Pm1_global, Q1P0, Q1Q1, FDSTAG } \n");
		}
		else {
			my_type = (DAVPElementType)type_k;
		}
	}

	/* get the function pointer */
	fp = PETSC_NULL;
	fp = DAVP_constructor_list[ my_type ];


	/* allocate and fill out the base class info */
	ierr = PetscMalloc( sizeof(struct _p_LaMEMVelPressureDA), &C ); CHKERRQ(ierr);
	ierr = PetscMemzero( C, sizeof(struct _p_LaMEMVelPressureDA) ); CHKERRQ(ierr);

	/* textual info */
	C->type = my_type;
	asprintf( &C->type_name, "%s", DAVPElementType_NAME[ my_type ] );
	asprintf( &C->type_description, "%s", DAVPElementType_DESCRIPTION[ my_type ] );

	/* implmenentation info */
	if( fp ) {
		fp( C );
	}
	else {
		PetscPrintf( PETSC_COMM_WORLD, "LaMEMVelPressureDACreate: Element type %s does not appear to have constructor implemented. \n", C->type_name);
		PetscPrintf( PETSC_COMM_WORLD, "LaMEMVelPressureDACreate: Unrecognised type. Choose one of { DAVP_Q2PM1L, DAVP_Q2PM1G, DAVP_Q1P0, DAVP_Q1Q1, DAVP_FDSTAG } \n");
		SETERRQ( PETSC_COMM_WORLD, PETSC_ERR_SUP, "Element type either not supported or not implemented");
	}


	/* check buffers */
	if( C->nnel > MAX_nnel ) {
		SETERRQ( PETSC_COMM_WORLD, PETSC_ERR_MEM, "MAX_nnel is not large enough");
	}
	if( C->ngp_vel > MAX_ngp_vel ) {
		SETERRQ( PETSC_COMM_WORLD, PETSC_ERR_MEM, "MAX_ngp_vel is not large enough");
	}
	if( C->nintp_1D > MAX_nintp_1D ) {
		SETERRQ( PETSC_COMM_WORLD, PETSC_ERR_MEM, "MAX_nintp_1D is not large enough");
	}
	if( C->npres > MAX_npres ) {
		SETERRQ( PETSC_COMM_WORLD, PETSC_ERR_MEM, "MAX_npres is not large enough");
	}
	if( C->edof > MAX_edof ) {
		SETERRQ( PETSC_COMM_WORLD, PETSC_ERR_MEM, "MAX_edof is not large enough");
	}
	if( C->edof_temp > MAX_edof_temp ) {
		SETERRQ( PETSC_COMM_WORLD, PETSC_ERR_MEM, "MAX_edof_temp is not large enough");
	}
	if( C->nnode_el_1D > MAX_nnode_el_1D ) {
		SETERRQ( PETSC_COMM_WORLD, PETSC_ERR_MEM, "MAX_nnode_el_1D is not large enough");
	}
	if( C->nnel_1D > MAX_nnel_1D ) {
		SETERRQ( PETSC_COMM_WORLD, PETSC_ERR_MEM, "MAX_nnel_1D is not large enough");
	}
	if( C->ndim > MAX_ndim ) {
		SETERRQ( PETSC_COMM_WORLD, PETSC_ERR_MEM, "MAX_ndim is not large enough");
	}


	/* set pointer */
	*_C = C;

	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "LaMEMVelPressureDADestroy"
PetscErrorCode LaMEMVelPressureDADestroy( LaMEMVelPressureDA *C )
{
	PetscErrorCode ierr;

	free( (*C)->type_name );
	free( (*C)->type_description );

	ierr = PetscFree( *C ); CHKERRQ(ierr);
	*C = PETSC_NULL;

	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "LaMEMVelPressureDAView"
PetscErrorCode LaMEMVelPressureDAView( LaMEMVelPressureDA C )
{
	PetscPrintf( PETSC_COMM_WORLD, "LaMEMVelPressureDA Object:\n");
	PetscPrintf( PETSC_COMM_WORLD, "\ttype: %s (index=%lld)\n",C->type_name, (LLD)(C->type));
	PetscPrintf( PETSC_COMM_WORLD, "\tdescription: %s\n",C->type_description);
	PetscPrintf( PETSC_COMM_WORLD, "\tParameters:\n");
	PetscPrintf( PETSC_COMM_WORLD, "\t\telement type = %lld\n", (LLD)(C->ElementType));
	PetscPrintf( PETSC_COMM_WORLD, "\t\tnodes per element (velocity)      = %lld\n",(LLD)(C->nnel));
	PetscPrintf( PETSC_COMM_WORLD, "\t\tdofs per element (velocity)       = %lld\n",(LLD)(C->edof));
	PetscPrintf( PETSC_COMM_WORLD, "\t\tnodes per element [1d] (velocity) = %lld\n",(LLD)(C->nnel_1D));
	PetscPrintf( PETSC_COMM_WORLD, "\t\tgauss points per element          = %lld\n",(LLD)(C->ngp_vel));
	PetscPrintf( PETSC_COMM_WORLD, "\t\tgauss points per element [1d]     = %lld\n",(LLD)(C->nintp_1D));
	PetscPrintf( PETSC_COMM_WORLD, "\t\tnode per element (pressure)       = %lld\n",(LLD)(C->npres));
	PetscPrintf( PETSC_COMM_WORLD, "\t\tnodes per element (temperature)      = %lld\n",(LLD)(C->nnel));
	PetscPrintf( PETSC_COMM_WORLD, "\t\tdofs per element (temperature)       = %lld\n",(LLD)(C->edof_temp));
	PetscPrintf( PETSC_COMM_WORLD, "\t\tnodes per element [1d] (temperature) = %lld\n",(LLD)(C->nnel_1D));
	PetscPrintf( PETSC_COMM_WORLD, "\t\tspatial dimension                    = %lld\n", (LLD)(C->ndim));



	PetscFunctionReturn(0);
}

/* interfaces */
#undef __FUNCT__
#define __FUNCT__ "LaMEMVelPressureDAGetInfo"
PetscErrorCode LaMEMVelPressureDAGetInfo( LaMEMVelPressureDA C,
		DAVPElementType *type,
		char **type_name,
		char **type_description,
		PetscInt *nnel,
		PetscInt *ngp_vel,
		PetscInt *nintp_1D,
		PetscInt *ElementType,
		PetscInt *npres,
		PetscInt *edof,
		PetscInt *edof_temp,
		PetscInt *nnode_el_1D,
		PetscInt *nnel_1D,
		PetscInt *ndim )
{
	if( type 				!= PETSC_NULL ) {		*type 				= C->type;					}
	if( type_name			!= PETSC_NULL ) {		*type_name 			= C->type_name;				}
	if( type_description 	!= PETSC_NULL ) {		*type_description 	= C->type_description;		}

	if( nnel 		!= PETSC_NULL ) {	*nnel 		= C->nnel;			}
	if( ngp_vel 	!= PETSC_NULL ) {	*ngp_vel 	= C->ngp_vel;		}
	if( nintp_1D 	!= PETSC_NULL ) {	*nintp_1D 	= C->nintp_1D;		}

	if( ElementType 	!= PETSC_NULL ) {	*ElementType 	= C->ElementType;	}
	if( npres 			!= PETSC_NULL ) {	*npres 			= C->npres;		}
	if( edof 			!= PETSC_NULL ) {	*edof 			= C->edof;			}

	if( edof_temp 		!= PETSC_NULL ) {	*edof_temp 		= C->edof_temp;	}
	if( nnode_el_1D 	!= PETSC_NULL ) {	*nnode_el_1D 	= C->nnode_el_1D;	}
	if( nnel_1D 		!= PETSC_NULL ) {	*nnel_1D 		= C->nnel_1D;		}
	if( ndim 			!= PETSC_NULL ) {	*ndim 			= C->ndim;			}

	PetscFunctionReturn(0);
}

PetscErrorCode LaMEMVelPressureDAGetBCInfo( LaMEMVelPressureDA C, PetscInt **left, PetscInt **right, PetscInt **front, PetscInt **back, PetscInt **lower, PetscInt **upper )
{
	if( left != PETSC_NULL ) {		*left = C->boundary_nodes_left;		}
	if( right != PETSC_NULL ) {		*right = C->boundary_nodes_right;	}

	if( front != PETSC_NULL ) {		*front = C->boundary_nodes_front;	}
	if( back != PETSC_NULL ) {		*back = C->boundary_nodes_back;		}

	if( lower != PETSC_NULL ) {		*lower = C->boundary_nodes_lower;	}
	if( upper != PETSC_NULL ) {		*upper = C->boundary_nodes_upper;	}

	PetscFunctionReturn(0);
}


void VPT_ComputeElementIndex(
		LaMEMVelPressureDA C,
		Particles ParticleLocal, PetscInt *_ix, PetscInt *_iy, PetscInt *_iz )
{

	C->fp_ComputeElementIndex( C, ParticleLocal, _ix,_iy,_iz );
}

void VPT_CreateStencilInGlobalStiffnessTemp(
		LaMEMVelPressureDA C,
		MatStencil *row,MatStencil *col,const PetscInt i,const PetscInt j, const PetscInt k )
{
	C->fp_CreateStencilInGlobalStiffnessTemp( C, row, col, i,j,k );
}

void VPT_GetTemperatureElement(
		LaMEMVelPressureDA C,
		PetscScalar ***temperature, PetscScalar Temp_element[], PetscInt i, PetscInt j, PetscInt k )
{
	C->fp_GetTemperatureElement( C, temperature, Temp_element, i,j,k );
}

void VPT_SetValuesRHS_Temp( LaMEMVelPressureDA C, PetscScalar ***rhs, PetscScalar T_RHS[], PetscInt i, PetscInt j, PetscInt k )
{
	C->fp_SetValuesRHS_Temp( C, rhs, T_RHS, i,j,k );
}

void VPT_FindNearestNode(
		LaMEMVelPressureDA C,
		PetscScalar eta, PetscScalar zetha, PetscScalar phi, PetscInt *ix_add_out, PetscInt *iy_add_out, PetscInt *iz_add_out )
{
	C->fp_FindNearestNode( C, eta, zetha, phi, ix_add_out, iy_add_out, iz_add_out );
}


void VPT_ComputeVelocityLocal2Global(
		LaMEMVelPressureDA C,
		MatStencil *row, MatStencil *col,
		PetscInt local2global[],
		const PetscInt i ,const PetscInt j , const PetscInt k,
		Mat MATRIX )
{
	C->fp_ComputeVelocityLocal2Global( C, row, col, local2global, i,j,k, MATRIX );
}


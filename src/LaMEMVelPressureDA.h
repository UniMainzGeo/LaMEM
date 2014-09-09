
#ifndef __LaMEMVelPressureDA_h__
#define __LaMEMVelPressureDA_h__


/*
NOTE:
If you add an element, you must
  1) Edit the typedef enum DAVPElementType.
  2) Add you element type constructors protoype to LaMEMVelPressureDA.c
  3) Add your constructor function to the list, DAVP_Constructor DAVP_list[].
  4) Add a name for your element type to the list, DAVPElementType_NAME[].
     This textual name will be searched for on the command line.
  5) Add a description of your element into the list, DAVPElementType_DESCRIPTION[].

NOTE:
If I want to add an element, I do it in a much more clear way.

*** IMPORTANT ***
*** IT IS IMPERATIVE THAT THE ORDER OF THE ENTRIES IN THESE LISTS ARE CONSISTENT.
*** DO NOT JUMBLE THE ENTRIES UP.
*** THE ARRAYS' DAVP_list[], DAVPElementType_NAME[], DAVPElementType_DESCRIPTION[] MUST TERMINATE WITH NULL (IE 0)
*** THE ENUM DAVPElementType SHOULD TERMINATE WITH  __DAVP_LIST_TERMINATOR__.
*/


typedef enum {
			DAVP_Q2PM1L=0,
			DAVP_Q2PM1G,
			DAVP_Q1P0,
			DAVP_Q1Q1,
			DAVP_FDSTAG,
			__DAVP_LIST_TERMINATOR__
		} DAVPElementType;

typedef struct _p_LaMEMVelPressureDA* LaMEMVelPressureDA;


PetscErrorCode LaMEMVelPressureDACreate( DAVPElementType type, LaMEMVelPressureDA *_C );
PetscErrorCode LaMEMVelPressureDADestroy( LaMEMVelPressureDA *C );
PetscErrorCode LaMEMVelPressureDAView( LaMEMVelPressureDA C );

void VPT_ComputeElementIndex( LaMEMVelPressureDA C, Particles ParticleLocal, PetscInt *_ix, PetscInt *_iy, PetscInt *_iz );
void VPT_CreateStencilInGlobalStiffnessTemp( LaMEMVelPressureDA C, MatStencil *row,MatStencil *col,const PetscInt i,const PetscInt j, const PetscInt k );
void VPT_GetTemperatureElement( LaMEMVelPressureDA C, PetscScalar ***temperature, PetscScalar Temp_element[], PetscInt i, PetscInt j, PetscInt k );
void VPT_SetValuesRHS_Temp( LaMEMVelPressureDA C, PetscScalar ***rhs, PetscScalar T_RHS[], PetscInt i, PetscInt j, PetscInt k );
void VPT_FindNearestNode( LaMEMVelPressureDA C, PetscScalar eta, PetscScalar zetha, PetscScalar phi, PetscInt *ix_add_out, PetscInt *iy_add_out, PetscInt *iz_add_out );

PetscErrorCode LaMEMVelPressureDAGetInfo( LaMEMVelPressureDA C,
		DAVPElementType *type, char **type_name, char **type_description,
		PetscInt *_nnel, PetscInt *_ngp_vel, PetscInt *_nintp_1D, PetscInt *_ElementType,
		PetscInt *_npres, PetscInt *_edof, PetscInt *_edof_temp, PetscInt *_nnode_el_1D,
		PetscInt *_nnel_1D, PetscInt *_ndim );
PetscErrorCode LaMEMVelPressureDAGetBCInfo( LaMEMVelPressureDA C, PetscInt **left, PetscInt **right, PetscInt **front, PetscInt **back, PetscInt **lower, PetscInt **upper );


void VPT_ComputeVelocityLocal2Global(
		LaMEMVelPressureDA C,
		MatStencil *row, MatStencil *col,
		PetscInt local2global[],
		const PetscInt i ,const PetscInt j , const PetscInt k,
		Mat MATRIX );

#endif

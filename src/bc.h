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
 **    filename:   bc.h
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
//........................... BOUNDARY CONDITIONS ...........................
//---------------------------------------------------------------------------
#ifndef __bc_h__
#define __bc_h__
//---------------------------------------------------------------------------

struct FB;
struct Scaling;
struct TSSol;
struct FDSTAG;
struct Marker;
struct DBMat;
struct JacRes;

//---------------------------------------------------------------------------
// index shift type
enum ShiftType
{
	_LOCAL_TO_GLOBAL_,
	_GLOBAL_TO_LOCAL_

};

//---------------------------------------------------------------------------
// Bezier block (rotating polygon moving along Bezier curve)
//---------------------------------------------------------------------------

struct BCBlock
{
	// path description
	PetscInt    npath;                        // number of path points of Bezier curve
	PetscScalar theta[  _max_path_points_  ]; // orientation angles at path points
	PetscScalar time [  _max_path_points_  ]; // times at path points
	PetscScalar path [6*_max_path_points_-4]; // Bezier curve path & control points (3*n-2)

	// block description
	PetscInt    npoly;                      // number of polygon vertices
	PetscScalar poly [2*_max_poly_points_]; // polygon coordinates
	PetscScalar bot, top;                   // bottom & top coordinates of the block

	// WARNING bottom coordinate should be advected (how? average?)
	// Top of the box can be assumed to be the free surface
	// sticky air nodes should never be constrained (this is easy to check)

};

//---------------------------------------------------------------------------

// setup data structures
PetscErrorCode BCBlockCreate(BCBlock *bcb, Scaling *scal, FB *fb);

// compute position along the path and rotation angle as a function of time
PetscErrorCode BCBlockGetPosition(BCBlock *bcb, PetscScalar t, PetscInt *f, PetscScalar x[]);

// compute current polygon coordinates
PetscErrorCode BCBlockGetPolygon(BCBlock *bcb, PetscScalar Xb[], PetscScalar *cpoly);

//---------------------------------------------------------------------------
// Internal velocity boxes (rectangular boxes with constant prescribed velocity that are either fixed or move)
//---------------------------------------------------------------------------

struct VelBox
{
	PetscInt 	advect;  // box advection flag
	PetscScalar cenX;    // x-coordinates of center
	PetscScalar cenY;    // y-coordinates of center
	PetscScalar cenZ;    // z-coordinates of center
	PetscScalar widthX;  // Width in x
	PetscScalar widthY;  // Width in y
	PetscScalar widthZ;  // Width in z
	PetscScalar vx;      // Vx-velocity within box
	PetscScalar vy;      // Vy-velocity within box
	PetscScalar vz;      // Vz-velocity within box
};

//---------------------------------------------------------------------------
// Internal velocity cylinders (cylinders with constant prescribed velocity that are either fixed or move)
//---------------------------------------------------------------------------

struct VelCylinder
{
	PetscInt 	advect;  // cylinder advection flag
	PetscScalar baseX;   // x-coordinates of base
	PetscScalar baseY;   // y-coordinates of base
	PetscScalar baseZ;   // z-coordinates of base
	PetscScalar capX;    // x-coordinates of cap
	PetscScalar capY;    // y-coordinates of cap
	PetscScalar capZ;    // z-coordinates of cap
	PetscScalar rad;     // radius of cylidner
	PetscScalar vx;      // Vx-velocity within box
	PetscScalar vy;      // Vy-velocity within box
	PetscScalar vz;      // Vz-velocity within box
	PetscScalar vmag;    // velocity magnitude
	PetscInt    type;    // velocity profile type
};

//---------------------------------------------------------------------------

PetscErrorCode VelBoxCreate(VelBox *velbox, Scaling *scal, FB *fb);
PetscErrorCode VelBoxPrint (VelBox *velbox, Scaling *scal, PetscInt cnt);
PetscErrorCode VelCylinderCreate(VelCylinder *velcyl, Scaling *scal, FB *fb);
PetscErrorCode VelCylinderPrint (VelCylinder *velcyl, Scaling *scal, PetscInt cnt);

//---------------------------------------------------------------------------

// boundary condition context
struct BCCtx
{
	//=====================================================================
	//
	// Boundary condition vectors contain prescribed DOF values:
	//
	//    *Internal points (marked with positive number in the index arrays)
	//        DBL_MAX   - active DOF flag
	//        otherwise - single-point constraint (SPC) value
	//
	//    *Boundary ghost point (marked with -1 in the index arrays)
	//        DBL_MAX   - free-slip (zero-flux) condition flag
	//        otherwise - two-point constraint (TPC) value
	//
	// Boundary ghost points require consistent setting of constraints
	// on the processor boundaries (since PETSc doesn't exchange boundary
	// ghost point values). Internal ghost points should be synchronized
	// after initializing the single-point constraints. Synchronization
	// can be skipped if all ghost points are initialized redundantly
	// on all the processes (DO THIS!).
	//
	// Single point constraints are additionally stored as lists
	// for constraining matrices and vectors. Matrices require global
	// index space, vectors require local index space.
	//
	// Global v-p index space can be either monolithic or block
	// Local v-p index space is ALWAYS coupled (since all solvers are coupled)
	//
	// NOTE! It may be worth storing TPC also as lists (for speedup).
	//=====================================================================


	FDSTAG   *fs;   // staggered grid
	TSSol    *ts;   // time stepping parameters
	Scaling  *scal; // scaling parameters
	DBMat    *dbm;  // material database
	JacRes   *jr;   // Jacobian-residual context (CROSS-REFERENCE!)


	Vec bcvx, bcvy, bcvz, bcp, bcT; // local (ghosted)

	// single-point constraints
	ShiftType    stype;   // current index shift type
	PetscInt     numSPC;  // total number of constraints
	PetscInt    *SPCList; // local indices of SPC
	PetscScalar *SPCVals; // values of SPC

	// velocity
	PetscInt     vNumSPC;
	PetscInt    *vSPCList;
	PetscScalar *vSPCVals;

	// pressure
	PetscInt     pNumSPC;
	PetscInt    *pSPCList;
	PetscScalar *pSPCVals;

	// temperature
	PetscInt     tNumSPC;
	PetscInt    *tSPCList;
	PetscScalar *tSPCVals;

	// two-point constraints
//	PetscInt     numTPC;       // number of two-point constraints (TPC)
//	PetscInt    *TPCList;      // local indices of TPC (ghosted layout)
//	PetscInt    *TPCPrimeDOF;  // local indices of primary DOF (ghosted layout)
//	PetscScalar *TPCVals;      // values of TPC
//	PetscScalar *TPCLinComPar; // linear combination parameters

	//=====================
	// VELOCITY CONSTRAINTS
	//=====================

	// horizontal background strain-rate parameters
	PetscInt     ExxNumPeriods;
	PetscScalar  ExxTimeDelims [_max_periods_-1];
	PetscScalar  ExxStrainRates[_max_periods_  ];

	PetscInt     EyyNumPeriods;
	PetscScalar  EyyTimeDelims [_max_periods_-1];
	PetscScalar  EyyStrainRates[_max_periods_  ];

	// simple shear background strain-rate parameters
	PetscInt     ExyNumPeriods;
	PetscScalar  ExyTimeDelims [_max_periods_-1];
	PetscScalar  ExyStrainRates[_max_periods_  ];
	
	PetscInt     EyzNumPeriods;
	PetscScalar  EyzTimeDelims [_max_periods_-1];
	PetscScalar  EyzStrainRates[_max_periods_  ];

	PetscInt     ExzNumPeriods;
	PetscScalar  ExzTimeDelims [_max_periods_-1];
	PetscScalar  ExzStrainRates[_max_periods_  ];

	// background strain rate reference point
	PetscScalar  BGRefPoint[3];

	// Bezier block
	PetscInt 	 nblocks;             // number of Bezier blocks
	BCBlock      blocks[_max_boxes_]; // BC block

    // internal velocity boxes
	PetscInt 	 nboxes;              // number of velocity boxes
    VelBox       vboxes[_max_boxes_]; // velocity boxes

	// internal velocity cylinders
	PetscInt 	 ncylinders;                            // number of velocity boxes
    VelCylinder  vcylinders[_max_boxes_];               // velocity boxes

	// velocity inflow & outflow boundary condition
	PetscInt     face,face_out,num_phase_bc,phase[5];   // face (1-left 2-right 3-front 4-back) & phase identifiers
	PetscScalar  bot, top,relax_dist,phase_interval[6]; // bottom & top coordinates of the plate
	PetscScalar  velin, velout; 			            // inflow & outflow velocities
	PetscScalar  velbot, veltop; 		                // bottom/top inflow velocities
	PetscInt     bvel_temperature_inflow;
	PetscScalar  bvel_thermal_age,bvel_potential_temperature, bvel_temperature_top;
	PetscScalar  bvel_constant_temperature;
	PetscInt     VelNumPeriods; 			// number of periods when boundary inflow velocity will change , must be less than _max_periods_
    PetscScalar  VelTimeDelims [_max_periods_-1];
    PetscScalar  velin_array [_max_periods_];

	// Plume inflow bottom boundary condition
	PetscInt        Plume_Inflow;
	PetscInt		Plume_Type;                 // Do we have a plume-like inflow boundary?
	PetscInt 		Plume_flux_ctr;             // Plume flux is constrained or not?
	PetscInt 		Plume_Dimension;            // Type [2D=1, or 3D=2]
	PetscInt		Plume_Phase;                // Phase of plume
	PetscInt        Plume_Phase_Mantle;         // Mantle phase (astenosphere)
	PetscScalar     Plume_Depth;                // Column plume height
	PetscScalar		Plume_Temperature;          // Temperature
	PetscScalar		Plume_Center[2];            // center [x,y] coordinates (for 3D plume)		
	PetscScalar		Plume_Radius;               // radius of plume (for 3D plume)
	PetscScalar		Plume_Inflow_Velocity;      // inflow velocity
	PetscInt 		Plume_VelocityType;         // type of inflow [Gaussian=0=default or Poiseuille=1]
	PetscScalar     Plume_areaFrac;             // how much of the plume area is actually in the model
	PetscScalar     Plume_Pressure;             // Plume Pressure at the bottom of the model (i.e. the bottom pressure boundary condition)
	// open boundary flag
	PetscInt     	top_open;
	PetscInt 		bot_open;
	PetscInt        phase_inflow_bot;

	// no-slip boundary condition mask
	PetscInt     	noslip[6];

	// fixed phase (no-flow condition)
	PetscInt     	fixPhase;

	// fixed cells (no-flow condition)
	PetscInt       	fixCell;
	unsigned char 	*fixCellFlag;

	//========================
	// TEMPERATURE CONSTRAINTS
	//========================

	// temperature on top and bottom boundaries and initial guess activation flag
	// bottom T can change with time
	PetscInt     	TbotNumPeriods;
	PetscScalar  	TbotTimeDelims [_max_periods_-1];
	PetscScalar  	Tbot[_max_periods_  ];

	PetscScalar  	Ttop;
	PetscInt     	initTemp;

	//=====================
	// PRESSURE CONSTRAINTS
	//=====================

	// pressure on top and bottom boundaries and initial guess activation flag
	PetscScalar  	pbot, ptop;
	PetscInt     	initPres;

};
//---------------------------------------------------------------------------

// create boundary condition context
PetscErrorCode BCCreate(BCCtx *bc, FB *fb);

// read boundary condition context from restart database
PetscErrorCode BCReadRestart(BCCtx *bc, FILE *fp);

// write boundary condition context to restart database
PetscErrorCode BCWriteRestart(BCCtx *bc, FILE *fp);

// allocate internal vectors and arrays
PetscErrorCode BCCreateData(BCCtx *bc);

// destroy boundary condition context
PetscErrorCode BCDestroy(BCCtx *bc);

// read fixed cells from files in parallel
PetscErrorCode BCReadFixCell(BCCtx *bc, FB *fb);

// apply ALL boundary conditions
PetscErrorCode BCApply(BCCtx *bc);

// apply SPC to global solution vector
PetscErrorCode BCApplySPC(BCCtx *bc);

// shift indices of constrained nodes
PetscErrorCode BCShiftIndices(BCCtx *bc, ShiftType stype);

// plume pressure like boundary condition
PetscErrorCode BCApplyPres_Plume_Pressure(BCCtx *bc);


//---------------------------------------------------------------------------
// Specific constraints
//---------------------------------------------------------------------------

// apply pressure constraints
PetscErrorCode BCApplyPres(BCCtx *bc);

// apply temperature constraints
PetscErrorCode BCApplyTemp(BCCtx *bc);

// apply default velocity constraints on the boundaries
PetscErrorCode BCApplyVelDefault(BCCtx *bc);

// apply Bezier blocks
PetscErrorCode BCApplyBezier(BCCtx *bc);

// apply inflow/outflow boundary velocities
PetscErrorCode BCApplyBoundVel(BCCtx *bc);

// apply internal velocity boxes
PetscErrorCode BCApplyVelBox(BCCtx *bc);

// apply internal velocity cylinders
PetscErrorCode BCApplyVelCylinder(BCCtx *bc);

// constraint all cells containing phase
PetscErrorCode BCApplyPhase(BCCtx *bc);

// constrain cells marked in files in parallel
PetscErrorCode BCApplyCells(BCCtx *bc);

// create SPC constraint lists
PetscErrorCode BCListSPC(BCCtx *bc);

// apply two-point constraints on the boundaries
PetscErrorCode BCApplyVelTPC(BCCtx *bc);

// apply plume_open_boundary condition
PetscErrorCode BC_Plume_inflow(BCCtx *bc);

// apply pressure permeable boundary condition
PetscErrorCode BCApply_Permeable_Pressure(BCCtx *bc);

// Get the average lithostatic pressure at the bottom
PetscErrorCode GetAverageLithostatic(BCCtx *bc);

// Get the densities of the external material
PetscScalar GetDensity(BCCtx *bc,PetscInt Phase, PetscScalar T, PetscScalar p );



//---------------------------------------------------------------------------
// Service functions
//---------------------------------------------------------------------------

// get current background strain rates & reference point coordinates
PetscErrorCode BCGetBGStrainRates(
		BCCtx       *bc,
		PetscScalar *Exx_,
		PetscScalar *Eyy_,
		PetscScalar *Ezz_,
		PetscScalar *Exy_,
		PetscScalar *Eyz_,
		PetscScalar *Exz_,
		PetscScalar *Rxx_,
		PetscScalar *Ryy_,
		PetscScalar *Rzz_);
//change velin in accordance with given time intervals
PetscErrorCode BCGetVelins(
		BCCtx       *bc);

// Get current bottom temperature
PetscErrorCode BCGetTempBound(
		BCCtx       *bc,
		PetscScalar *Tbot);

// stretch staggered grid if background strain rates are defined
PetscErrorCode BCStretchGrid(BCCtx *bc);

// change phase of inflow markers if velocity boundary condition is defined
PetscErrorCode BCOverridePhase(BCCtx *bc, PetscInt cellID, Marker *P);

//---------------------------------------------------------------------------
// MACROS
//---------------------------------------------------------------------------

#define LIST_SPC(bc, list, vals, cnt, iter)\
	if(bc[k][j][i] != DBL_MAX) { list[cnt] = iter; vals[cnt] = bc[k][j][i]; cnt++; }

//---------------------------------------------------------------------------
#endif

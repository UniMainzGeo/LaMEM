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
//........................... BOUNDARY CONDITIONS ...........................
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "bc.h"
#include "JacRes.h"
#include "parsing.h"
#include "scaling.h"
#include "tssolve.h"
#include "fdstag.h"
#include "tools.h"
#include "Tensor.h"
#include "advect.h"
#include "phase.h"
#include "constEq.h"
#include "surf.h"
//---------------------------------------------------------------------------
// Bezier block functions
//---------------------------------------------------------------------------
PetscErrorCode BCBlockCreate(BCBlock *bcb, Scaling *scal, FB *fb)
{
	//	-npath    - Number of path points of Bezier curve (end-points only!)
	//	-path_dim - Path dimension: 2 = x-y plane (default), 3 = full 3D
	//	-theta    - Orientation angles at path points (counter-clockwise positive)
	//	-time     - Times at path points
	//	-path     - path points coordinates (x-y for 2D, x-y-z for 3D)
	//	-npoly    - Number of polygon vertices
	//	-poly     - Polygon x-y coordinates at initial time (absolute, moves with path)
	//	-bot      - Polygon bottom z-coordinate at initial time (absolute, moves with path for 3D)
	//	-top      - Polygon top z-coordinate at initial time (absolute, moves with path for 3D)

	
	PetscInt       numPathCoords;
	PetscFunctionBeginUser;

	// set defaults
	bcb->npath   = 2;
	bcb->pathDim = 2;
	bcb->npoly   = 4;

	PetscCall(getIntParam   (fb, _OPTIONAL_, "npath",    &bcb->npath,   1, _max_path_points_));
	PetscCall(getIntParam   (fb, _OPTIONAL_, "path_dim", &bcb->pathDim, 1, 3));

	// validate path_dim
	if(bcb->pathDim != 2 && bcb->pathDim != 3)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "path_dim must be 2 or 3, got: %" PetscInt_FMT "", bcb->pathDim);
	}

	// compute number of path coordinates based on dimension
	numPathCoords = bcb->pathDim * bcb->npath;

	PetscCall(getScalarParam(fb, _OPTIONAL_, "theta",  bcb->theta, bcb->npath,    scal->angle ));
	PetscCall(getScalarParam(fb, _REQUIRED_, "time",   bcb->time,  bcb->npath,    scal->time  ));
	PetscCall(getScalarParam(fb, _REQUIRED_, "path",   bcb->path,  numPathCoords, scal->length));

	PetscCall(getIntParam   (fb, _OPTIONAL_, "npoly", &bcb->npoly, 1,              _max_poly_points_));
	PetscCall(getScalarParam(fb, _REQUIRED_, "poly",   bcb->poly,  2*bcb->npoly,   scal->length     ));
	PetscCall(getScalarParam(fb, _REQUIRED_, "bot",   &bcb->bot,   1,              scal->length     ));
	PetscCall(getScalarParam(fb, _REQUIRED_, "top",   &bcb->top,   1,              scal->length     ));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode BCBlockPrint(BCBlock *bcb, Scaling *scal, PetscInt cnt)
{
	PetscFunctionBeginUser;

	PetscPrintf(PETSC_COMM_WORLD, "      Bezier block #                          : %" PetscInt_FMT " \n", cnt);
	PetscPrintf(PETSC_COMM_WORLD, "      Path dimension                          : %" PetscInt_FMT " \n", bcb->pathDim);
	PetscPrintf(PETSC_COMM_WORLD, "      Number of path points                   : %" PetscInt_FMT " \n", bcb->npath);
	PetscPrintf(PETSC_COMM_WORLD, "      Number of polygon vertices              : %" PetscInt_FMT " \n", bcb->npoly);

	PetscPrintf(PETSC_COMM_WORLD, "      Bot/Top initial z-coordinates          : %g / %g %s \n",
		bcb->bot*scal->length, bcb->top*scal->length, scal->lbl_length);

	if(bcb->pathDim == 3)
	{
		PetscPrintf(PETSC_COMM_WORLD, "      (bot/top move with 3D path)            @ \n");
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode BCBlockGetPosition(BCBlock *bcb, PetscScalar t, PetscInt *f, PetscScalar X[])
{
	// compute position along the path and rotation angle as a function of time
	// For 2D path (pathDim=2): X[0]=x, X[1]=y, X[2]=theta
	// For 3D path (pathDim=3): X[0]=x, X[1]=y, X[2]=z, X[3]=theta

	PetscInt      i, n, dim;
	PetscScalar   r, s;
	PetscScalar  *p1, *p2;
	PetscScalar  *path, *theta, *time;

	PetscFunctionBeginUser;

	n     = bcb->npath;
	dim   = bcb->pathDim;
	path  = bcb->path;
	theta = bcb->theta;
	time  = bcb->time;

	// set flag
	(*f) = 1; if(t < time[0] || t > time[n-1]) { (*f) = 0; PetscFunctionReturn(0); }

	// find time interval
	for(i = 1; i < n-1; i++) { if(t < time[i]) break; } i--;

	// get path and control points (stride depends on dimension)
	p1 = path + dim*i;
	p2 = p1   + dim;

	// compute interpolation parameters
	r  = (t - time[i])/(time[i+1] - time[i]);
	s  = 1.0 - r;

	// interpolate path and rotation angle
	X[0] = s*p1[0] + r*p2[0];  // x
	X[1] = s*p1[1] + r*p2[1];  // y

	if(dim == 3)
	{
		// 3D path: X[0]=x, X[1]=y, X[2]=z, X[3]=theta
		X[2] = s*p1[2]    + r*p2[2];      // z
		X[3] = s*theta[i] + r*theta[i+1]; // theta
	}
	else
	{
		// 2D path: X[0]=x, X[1]=y, X[2]=theta
		X[2] = s*theta[i] + r*theta[i+1]; // theta
	}

	//   [A] Bezier curves can be input directly.
	//   Bezier curve requires 4 points per segment (see e.g. wikipedia):
	//   path point P0 - control point P1 - control point P2 - path point P3.
	//   The last path point (P3) of every, but the last, interval is omitted due to continuity.
	//   Altogether, "path" variable should provide 3*npath-2 points.
	//   Every point has x and y coordinates, so total number of entries should be 6*npath-4.
	//   Bezier curves can be most easily generated using Inkscape software.
	//   Continuity of tangent lines can be imposed by the tool "make selected nodes symmetric"
	//   Coordinates of the curve points can be accessed using the XML editor in Inkscape.
	//   Alternatively one can process .svg files by geomIO software.

	//   [B] Alternative is to create smooth B-spline curves passing through the basic path points.
	//   Example (5 path points (S0 - S4), 4 Bezier segments):
	//   1) Solve for 3 B-control points (tri-diagonal system with 2 rhs & solution vectors one for x and one for y):
	//   | 4 1 0 |   | B1 |    | 6S1-S0 |
	//   | 1 4 1 | * | B2 | =  | 6S2    |
	//   | 0 1 4 |   | B3 |    | 6S3-S4 |
	//   End-points:
	//   B0 = S0
	//   B4 = S4
	//   2) Compute two Bezier control points for each segment form B-points:
	//   Example: Segment S1-S2
	//   Control points:
	//   P1=2/3*B1 + 1/3*B2
	//   P2=2/3*B2 + 1/3*B1

	//   [C] In any case Bezier curves and B-splines can not be used directly,
	//   since their curve parameter (t) maps nonlinearly on curve length, i.e:
	//   l(t=1/3) != L/3, where L in the total length of curve segment.
	//   This will lead to artificial "accelerations" along the curve path.
	//   Instead Bezier curves must be approximated by linear segments.
	//   This can be done adaptively by increasing number of subdivisions until approximate
	//   curve length converges to a loose relative tolerance (say 5-10%).

	//   [D] Code snippet:
	//   // get path and control points
	//   p1 = path + 6*i;
	//   p2 = p1 + 2;
	//   p3 = p2 + 2;
	//   p4 = p3 + 2;
	//   // compute interpolation parameters
	//   r  = (t - time[i])/(time[i+1] - time[i]);
	//   r2 = r*r;
	//   r3 = r2*r;
	//   s  = 1.0 - r;
	//   s2 = s*s;
	//   s3 = s2*s;
	//   // interpolate Bezier path
	//   X[0] = s3*p1[0] + 3.0*s2*r*p2[0] + 3.0*s*r2*p3[0] + r3*p4[0];
	//   X[1] = s3*p1[1] + 3.0*s2*r*p2[1] + 3.0*s*r2*p3[1] + r3*p4[1];

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode BCBlockGetPolygon(BCBlock *bcb, PetscScalar Xb[], PetscScalar *cpoly)
{
	// compute current polygon coordinates (2D x-y polygon)
	// For 2D path: Xb[0]=x, Xb[1]=y, Xb[2]=theta
	// For 3D path: Xb[0]=x, Xb[1]=y, Xb[2]=z, Xb[3]=theta

	PetscInt     i, dim;
	PetscScalar *xa, *xb;
	PetscScalar  Xa[4], thetaCur, thetaInit, theta, costh, sinth;

	PetscFunctionBeginUser;

	dim = bcb->pathDim;

	// get initial polygon position (x, y only for 2D rotation)
	Xa[0] = bcb->path[0];
	Xa[1] = bcb->path[1];
	thetaInit = bcb->theta[0];

	// get current rotation angle (theta is at different index for 2D vs 3D)
	if(dim == 3)
	{
		thetaCur = Xb[3];
	}
	else
	{
		thetaCur = Xb[2];
	}

	// get rotation matrix (rotation around z-axis)
	theta = thetaCur - thetaInit;
	costh = cos(theta);
	sinth = sin(theta);

	// compute current polygon coordinates
	for(i = 0; i < bcb->npoly; i++)
	{
		// get reference and current points
		xa = bcb->poly + 2*i;
		xb = cpoly     + 2*i;

		// rotate & displace (2D rotation in x-y plane)
		RotDispPoint2D(Xa, Xb, costh, sinth, xa, xb);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// Velocity box functions
//---------------------------------------------------------------------------
PetscErrorCode VelBoxCreate(VelBox *velbox, Scaling *scal, FB *fb)
{
	
	PetscFunctionBeginUser;

	//========================
	// velocity box parameters
	//========================

	velbox->vx = DBL_MAX;
	velbox->vy = DBL_MAX;
	velbox->vz = DBL_MAX;

	PetscCall(getScalarParam(fb, _REQUIRED_, "cenX",   &velbox->cenX,   1,  scal->length));
	PetscCall(getScalarParam(fb, _REQUIRED_, "cenY",   &velbox->cenY,   1,  scal->length));
	PetscCall(getScalarParam(fb, _REQUIRED_, "cenZ",   &velbox->cenZ,   1,  scal->length));
	PetscCall(getScalarParam(fb, _REQUIRED_, "widthX", &velbox->widthX, 1,  scal->length));
	PetscCall(getScalarParam(fb, _REQUIRED_, "widthY", &velbox->widthY, 1,  scal->length));
	PetscCall(getScalarParam(fb, _REQUIRED_, "widthZ", &velbox->widthZ, 1,  scal->length));
	PetscCall(getScalarParam(fb, _OPTIONAL_, "vx",     &velbox->vx,     1,  scal->velocity));
	PetscCall(getScalarParam(fb, _OPTIONAL_, "vy",     &velbox->vy,     1,  scal->velocity));
	PetscCall(getScalarParam(fb, _OPTIONAL_, "vz",     &velbox->vz,     1,  scal->velocity));
	PetscCall(getIntParam   (fb, _REQUIRED_, "advect", &velbox->advect, 1,  1));

	if(velbox->vx == DBL_MAX && velbox->vy == DBL_MAX && velbox->vz == DBL_MAX)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Velocity box should specify at least one velocity component");
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode VelBoxPrint(VelBox *velbox, Scaling *scal, PetscInt cnt)
{
	PetscFunctionBeginUser;

	PetscPrintf(PETSC_COMM_WORLD, "      Velocity box #                          : %" PetscInt_FMT " \n",  cnt);
	PetscPrintf(PETSC_COMM_WORLD, "      Box center                              : %g, %g, %g %s \n", velbox->cenX  *scal->length, velbox->cenY  *scal->length, velbox->cenZ  *scal->length, scal->lbl_length);
	PetscPrintf(PETSC_COMM_WORLD, "      Box width                               : %g, %g, %g %s \n", velbox->widthX*scal->length, velbox->widthY*scal->length, velbox->widthZ*scal->length, scal->lbl_length);
	if(velbox->vx != DBL_MAX)
	{
		PetscPrintf(PETSC_COMM_WORLD, "      X-velocity                              : %g %s \n", velbox->vx*scal->velocity, scal->lbl_velocity);
	}
	if(velbox->vy != DBL_MAX)
	{
		PetscPrintf(PETSC_COMM_WORLD, "      Y-velocity                              : %g %s \n", velbox->vy*scal->velocity, scal->lbl_velocity);
	}
	if(velbox->vz != DBL_MAX)
	{
		PetscPrintf(PETSC_COMM_WORLD, "      Z-velocity                              : %g %s \n", velbox->vz*scal->velocity, scal->lbl_velocity);
	}
	if(velbox->advect)
	{
		PetscPrintf(PETSC_COMM_WORLD, "      Advect velocity with flow               @  \n");
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// Velocity cylinder functions
//---------------------------------------------------------------------------
PetscErrorCode VelCylinderCreate(VelCylinder *velcyl, Scaling *scal, FB *fb)
{
	char           str_type[_str_len_];

	
	PetscFunctionBeginUser;

	//========================
	// velocity cylinder parameters
	//========================

	velcyl->vx   = DBL_MAX;
	velcyl->vy   = DBL_MAX;
	velcyl->vz   = DBL_MAX;
	velcyl->vmag = DBL_MAX;

	PetscCall(getScalarParam(fb, _REQUIRED_, "baseX",  &velcyl->baseX,  1,  scal->length));
	PetscCall(getScalarParam(fb, _REQUIRED_, "baseY",  &velcyl->baseY,  1,  scal->length));
	PetscCall(getScalarParam(fb, _REQUIRED_, "baseZ",  &velcyl->baseZ,  1,  scal->length));
	PetscCall(getScalarParam(fb, _REQUIRED_, "capX",   &velcyl->capX,   1,  scal->length));
	PetscCall(getScalarParam(fb, _REQUIRED_, "capY",   &velcyl->capY,   1,  scal->length));
	PetscCall(getScalarParam(fb, _REQUIRED_, "capZ",   &velcyl->capZ,   1,  scal->length));
	PetscCall(getScalarParam(fb, _REQUIRED_, "radius", &velcyl->rad,    1,  scal->length));
	PetscCall(getScalarParam(fb, _OPTIONAL_, "vx",     &velcyl->vx,     1,  scal->velocity));
	PetscCall(getScalarParam(fb, _OPTIONAL_, "vy",     &velcyl->vy,     1,  scal->velocity));
	PetscCall(getScalarParam(fb, _OPTIONAL_, "vz",     &velcyl->vz,     1,  scal->velocity));
	PetscCall(getScalarParam(fb, _OPTIONAL_, "vmag",   &velcyl->vmag,   1,  scal->velocity));
	PetscCall(getStringParam(fb, _OPTIONAL_, "type",    str_type,       "uniform"));
	PetscCall(getIntParam   (fb, _REQUIRED_, "advect", &velcyl->advect, 1,  1));

	if(!strcmp(str_type, "uniform"))
	{
		velcyl->type = 0;
	}
	else if(!strcmp(str_type, "parabolic"))
	{
		velcyl->type = 1;
	}
	else
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Velocity cylinder type must be uniform or parabolic");
	}

	if((velcyl->vx != DBL_MAX || velcyl->vy != DBL_MAX || velcyl->vz != DBL_MAX) && velcyl->vmag != DBL_MAX)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "For velocity cylinder, specify vmag or vx/vy/vz");
	}

	if(velcyl->vx == DBL_MAX && velcyl->vy == DBL_MAX && velcyl->vz == DBL_MAX && velcyl->vmag == DBL_MAX)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Velocity cylinder should specify at least one velocity component");
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode VelCylinderPrint(VelCylinder *velcyl, Scaling *scal, PetscInt cnt)
{
	PetscFunctionBeginUser;

	PetscPrintf(PETSC_COMM_WORLD, "      Velocity cylinder #                     : %" PetscInt_FMT " \n",  cnt);
	PetscPrintf(PETSC_COMM_WORLD, "      Cylinder base                           : %g, %g, %g %s \n", velcyl->baseX  *scal->length, velcyl->baseY  *scal->length, velcyl->baseZ  *scal->length, scal->lbl_length);
	PetscPrintf(PETSC_COMM_WORLD, "      Cylinder cap                            : %g, %g, %g %s \n", velcyl->capX*scal->length, velcyl->capY*scal->length, velcyl->capZ*scal->length, scal->lbl_length);
	PetscPrintf(PETSC_COMM_WORLD, "      Cylinder radius                         : %g %s \n", velcyl->rad*scal->length, scal->lbl_length);
	if(velcyl->vx != DBL_MAX)
	{
		PetscPrintf(PETSC_COMM_WORLD, "      X-velocity                              : %g %s \n", velcyl->vx*scal->velocity, scal->lbl_velocity);
	}
	if(velcyl->vy != DBL_MAX)
	{
		PetscPrintf(PETSC_COMM_WORLD, "      Y-velocity                              : %g %s \n", velcyl->vy*scal->velocity, scal->lbl_velocity);
	}
	if(velcyl->vz != DBL_MAX)
	{
		PetscPrintf(PETSC_COMM_WORLD, "      Z-velocity                              : %g %s \n", velcyl->vz*scal->velocity, scal->lbl_velocity);
	}
	if(velcyl->vmag != DBL_MAX)
	{
		PetscPrintf(PETSC_COMM_WORLD, "      velocity magnitude                      : %g %s \n", velcyl->vmag*scal->velocity, scal->lbl_velocity);
	}
	if(velcyl->type == 0)
	{
		PetscPrintf(PETSC_COMM_WORLD, "      velocity profile                        : uniform \n");
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "      velocity profile                        : parabolic \n");
	}
	if(velcyl->advect)
	{
		PetscPrintf(PETSC_COMM_WORLD, "      Advect velocity with flow               @  \n");
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// BCCtx functions
//---------------------------------------------------------------------------
PetscErrorCode BCCreate(BCCtx *bc, FB *fb)
{
	Scaling     *scal;
	FDSTAG      *fs;
	PetscInt     periodic;
	PetscInt     jj, mID;
	PetscScalar  bz;
	char         inflow_temp[_str_len_],str_inflow[_str_len_];

	
	PetscFunctionBeginUser;

	// access context
	scal = bc->scal;
	fs   = bc->fs;
	mID  = bc->dbm->numPhases-1;

	// set periodic flag
	periodic = fs->periodic;

	// initialize
	bc->Tbot[0]                 = -1.0;
	bc->Ttop                    = -1.0;
	bc->pbot                    = -1.0;
	bc->ptop                    = -1.0;
	bc->fixPhase                = -1;
	bc->num_phase_bc            =   -1;
	bc->velout                  =  DBL_MAX;
	bc->Plume_Inflow            = 0;
	bc->bvel_temperature_inflow = -1;

	//=====================
	// VELOCITY CONSTRAINTS
	//=====================

	// horizontal background strain-rate parameters
	PetscCall(getIntParam   (fb, _OPTIONAL_, "exx_num_periods",  &bc->ExxNumPeriods,  1,                   _max_periods_    ));
	PetscCall(getScalarParam(fb, _REQUIRED_, "exx_time_delims",   bc->ExxTimeDelims,  bc->ExxNumPeriods-1, scal->time       ));
	PetscCall(getScalarParam(fb, _REQUIRED_, "exx_strain_rates",  bc->ExxStrainRates, bc->ExxNumPeriods,   scal->strain_rate));
	PetscCall(getIntParam   (fb, _OPTIONAL_, "eyy_num_periods",  &bc->EyyNumPeriods,  1,                   _max_periods_    ));
	PetscCall(getScalarParam(fb, _REQUIRED_, "eyy_time_delims",   bc->EyyTimeDelims,  bc->EyyNumPeriods-1, scal->time       ));
	PetscCall(getScalarParam(fb, _REQUIRED_, "eyy_strain_rates",  bc->EyyStrainRates, bc->EyyNumPeriods,   scal->strain_rate));
	PetscCall(getScalarParam(fb, _OPTIONAL_, "bg_ref_point",      bc->BGRefPoint,     3,                   scal->length));

	// simple shear background strain-rate parameters
	PetscCall(getIntParam   (fb, _OPTIONAL_, "exy_num_periods",   &bc->ExyNumPeriods,  1,                   _max_periods_    ));
	PetscCall(getScalarParam(fb, _REQUIRED_, "exy_time_delims",    bc->ExyTimeDelims,  bc->ExyNumPeriods-1, scal->time       ));
	PetscCall(getScalarParam(fb, _REQUIRED_, "exy_strain_rates",   bc->ExyStrainRates, bc->ExyNumPeriods,   scal->strain_rate));

	// Bezier blocks
	PetscCall(FBFindBlocks(fb, _OPTIONAL_, "<BCBlockStart>", "<BCBlockEnd>"));

	if(fb->nblocks)
	{
		// error checking
		if(fb->nblocks > _max_boxes_)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many Bezier blocks! found: %" PetscInt_FMT ", max allowed: %" PetscInt_FMT "", fb->nblocks, _max_boxes_);
		}

		// store actual number of Bezier blocks
		bc->nblocks = fb->nblocks;

		// read Bezier blocks
		for(jj = 0; jj < fb->nblocks; jj++)
		{
			PetscCall(BCBlockCreate(bc->blocks + jj, scal, fb));

			fb->blockID++;
		}
	}

	PetscCall(FBFreeBlocks(fb));

	// velocity boxes
	PetscCall(FBFindBlocks(fb, _OPTIONAL_, "<VelBoxStart>", "<VelBoxEnd>"));

	if(fb->nblocks)
	{
		// error checking
		if(fb->nblocks > _max_boxes_)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many velocity boxes! found: %" PetscInt_FMT ", max allowed: %" PetscInt_FMT "", fb->nblocks, _max_boxes_);
		}

		// store actual number of velocity blocks
		bc->nboxes = fb->nblocks;

		// read velocity boxes
		for(jj = 0; jj < fb->nblocks; jj++)
		{
			PetscCall(VelBoxCreate(bc->vboxes + jj, scal, fb));

			fb->blockID++;
		}
	}

	PetscCall(FBFreeBlocks(fb));

	// velocity cylinders
	PetscCall(FBFindBlocks(fb, _OPTIONAL_, "<VelCylinderStart>", "<VelCylinderEnd>"));

	if(fb->nblocks)
	{
		// error checking
		if(fb->nblocks > _max_boxes_)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Too many velocity cylinders! found: %" PetscInt_FMT ", max allowed: %" PetscInt_FMT "", fb->nblocks, _max_boxes_);
		}

		// store actual number of velocity blocks
		bc->ncylinders = fb->nblocks;

		// read velocity boxes
		for(jj = 0; jj < fb->nblocks; jj++)
		{
			PetscCall(VelCylinderCreate(bc->vcylinders + jj, scal, fb));

			fb->blockID++;
		}
	}

	PetscCall(FBFreeBlocks(fb));

	// boundary inflow/outflow velocities
	PetscCall(getStringParam(fb, _OPTIONAL_, "bvel_face", str_inflow, NULL));  // must have component
	if      (!strcmp(str_inflow, "Left"))               bc->face=1;
	else if (!strcmp(str_inflow, "Right"))              bc->face=2;
	else if (!strcmp(str_inflow, "Front"))              bc->face=3;
	else if (!strcmp(str_inflow, "Back"))               bc->face=4;
	else if (!strcmp(str_inflow, "CompensatingInflow")) bc->face=5;

	PetscCall(getIntParam(fb, _OPTIONAL_, "bvel_face_out", &bc->face_out, 	1, -1));

	if(bc->face)
	{
		PetscCall(getIntParam   (fb, _OPTIONAL_, "bvel_num_phase",        &bc->num_phase_bc,     1, 5));
		PetscCall(getIntParam   (fb, _OPTIONAL_, "bvel_phase",             bc->phase,            bc->num_phase_bc, mID));
		PetscCall(getScalarParam(fb, _REQUIRED_, "bvel_bot",              &bc->bot,              1, scal->length  ));
		PetscCall(getScalarParam(fb, _REQUIRED_, "bvel_top",              &bc->top,              1, scal->length  ));
		PetscCall(getScalarParam(fb, _REQUIRED_, "bvel_velin",            &bc->velin,            1, scal->velocity));
		PetscCall(getScalarParam(fb, _OPTIONAL_, "bvel_velout",           &bc->velout,           1, scal->velocity));
		PetscCall(getIntParam   (fb, _OPTIONAL_, "velin_num_periods",     &bc->VelNumPeriods,    1, _max_periods_));
		PetscCall(getIntParam   (fb, _OPTIONAL_, "velin_net_num_periods", &bc->VelNetNumPeriods, 1, _max_periods_));
		PetscCall(getScalarParam(fb, _OPTIONAL_, "bvel_relax_d",          &bc->relax_dist,       1, scal->length ));

		if(bc->VelNumPeriods > 1)
		{
			PetscCall(getScalarParam(fb, _REQUIRED_, "velin_time_delims", bc->VelTimeDelims, bc->VelNumPeriods-1, scal->time    ));
			PetscCall(getScalarParam(fb, _REQUIRED_, "bvel_velin",        bc->velin_array,   bc->VelNumPeriods,   scal->velocity));
			PetscCall(BCGetVelins(bc));
		}

		if(bc->VelNetNumPeriods > 1)
		{
			PetscCall(getScalarParam(fb, _REQUIRED_, "velin_net_time_delims", bc->VelNetTimeDelims,  bc->VelNetNumPeriods-1, scal->time    ));
			PetscCall(getScalarParam(fb, _REQUIRED_, "bvel_velin_net",        bc->velin_net_array,   bc->VelNetNumPeriods,   scal->velocity));
			PetscCall(BCGetVelins(bc));
		}

		PetscCall(getScalarParam(fb, _OPTIONAL_, "bvel_phase_interval",     bc->phase_interval, bc->num_phase_bc+1, scal->length));
		PetscCall(getStringParam(fb, _OPTIONAL_, "bvel_temperature_inflow", inflow_temp , NULL));

		if(!strcmp(inflow_temp, "Constant_T_inflow")) { bc->bvel_temperature_inflow = 1; }
		if(!strcmp(inflow_temp, "Fixed_thermal_age")) { bc->bvel_temperature_inflow = 2; }

		if( bc->bvel_temperature_inflow == 2)
		{
			PetscCall(getScalarParam(fb, _REQUIRED_, "bvel_temperature_mantle", &bc->bvel_potential_temperature, 1, 1.0));
			PetscCall(getScalarParam(fb, _REQUIRED_, "bvel_temperature_top",    &bc->bvel_temperature_top,       1, 1.0));
			PetscCall(getScalarParam(fb, _REQUIRED_, "bvel_thermal_age",        &bc->bvel_thermal_age,           1, scal->time));
		}
		else if(bc->bvel_temperature_inflow == 1)
		{
			PetscCall(getScalarParam(fb, _REQUIRED_, "bvel_temperature_constant", &bc->bvel_constant_temperature, 1, scal->time));
		}

		PetscCall(getScalarParam(fb, _OPTIONAL_, "bvel_velbot", &bc->velbot, 1, scal->velocity));
		PetscCall(getScalarParam(fb, _OPTIONAL_, "bvel_veltop", &bc->veltop, 1, scal->velocity));

		PetscCall(FDSTAGGetGlobalBox(bc->fs, NULL, NULL, &bz, NULL, NULL, NULL));

		// compute outflow velocity (if required)
		if(bc->velout == DBL_MAX)
		{
			// INTRODUCE CORRECTION FOR CELL SIZES
			// MUST BE MASS CONSERVATIVE IN DISCRETE SENSE
			bc->velout = -bc->velin*(bc->top - bc->bot)/(bc->bot - bz);
		}
	}

	// open boundary flag
	PetscCall(getIntParam(fb, _OPTIONAL_, "open_top_bound", &bc->top_open, 1, -1));

	//open bottom boundary flag
	PetscCall(getIntParam(fb, _OPTIONAL_, "open_bot_bound", &bc->bot_open, 1, -1));

	if(bc->bot_open)
	{
		PetscCall(getIntParam(fb, _OPTIONAL_, "permeable_phase_inflow", &bc->phase_inflow_bot, 1, -1));
	}

	// no-slip boundary condition mask
	PetscCall(getIntParam(fb, _OPTIONAL_, "noslip",  bc->noslip, 6, -1));

	// fixed phase (no-flow condition)
	PetscCall(getIntParam(fb, _OPTIONAL_, "fix_phase", &bc->fixPhase, 1, mID));

	// fixed cells (no-flow condition)
	PetscCall(getIntParam(fb, _OPTIONAL_, "fix_cell",  &bc->fixCell, 1, mID));

	// plume-like inflow boundary condition @ bottom
	PetscCall(getIntParam(fb, _OPTIONAL_, "Plume_InflowBoundary", &bc->Plume_Inflow, 1, -1));

	if(bc->Plume_Inflow)
	{
		char str[_str_len_];

		// type of plume (2D or 3D)
		PetscCall(getStringParam(fb, _REQUIRED_, "Plume_Type", str, NULL));  // must have component

		// type of boundary conditions
		if(!strcmp(str, "Inflow_Type"))
		{
			bc->Plume_Type = 1; // velocity flux
		}
		else if(!strcmp(str, "Permeable_Type"))
		{
			bc->Plume_Type = 2; // activate open_bot boundary condition
			bc->bot_open   = 1; // open the bottom boundary
		}
		else
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Choose either [Influx_type; Permeable_Type] as parameter for Plume_Type, not %s",str);
		}

		if(bc->Plume_Type == 1)
		{
			bc->Plume_areaFrac = 1.0;

			PetscCall(getScalarParam(fb,_REQUIRED_,  "Plume_Inflow_Velocity", &bc->Plume_Inflow_Velocity, 1, scal->velocity));
			PetscCall(getStringParam(fb, _REQUIRED_, "Plume_VelocityType",    str, "Gaussian"));  // must have component
			PetscCall(getScalarParam(fb, _OPTIONAL_, "Plume_areaFrac",         &bc->Plume_areaFrac,       1,  1.0));

			if(!strcmp(str, "Poiseuille"))
			{
				bc->Plume_VelocityType = 0; // Poiseuille
			}
			else if(!strcmp(str, "Gaussian"))
			{
				bc->Plume_VelocityType = 1; // Gaussian perturbation (smoother)
			}
			else
			{
				SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Choose either [Poiseuille; Gaussian] as parameter for Plume_VelocityType, not %s",str);
			}
		}
		if(bc->Plume_Type == 2)
		{
			// bc->Plume_Pressure = -1;
			// PetscCall(getScalarParam(fb,_REQUIRED_,"Plume_Depth",    &bc->Plume_Depth,    1, scal->length));
			// PetscCall(getScalarParam(fb,_OPTIONAL_,"Plume_Pressure", &bc->Plume_Pressure, 1, scal->stress));
			PetscCall(getIntParam(fb, _REQUIRED_, "Plume_Phase_Mantle"  , &bc->phase_inflow_bot, 1, mID));
		}

		// 2D or 3D
		PetscCall(getStringParam(fb, _REQUIRED_, "Plume_Dimension", str, NULL));  // must have component

		if(!strcmp(str, "2D"))
		{
			bc->Plume_Dimension = 1; // 2D setup
		}
		else if(!strcmp(str, "3D"))
		{
			bc->Plume_Dimension = 2; // 3D (circular)
		}
		else
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Choose either [2D; 3D] as parameter for Plume_Type, not %s",str);
		}

		if(bc->Plume_Dimension == 1)
		{
			// 2D perturbation in x-direction
			PetscCall(getScalarParam(fb,_REQUIRED_, "Plume_Center",  bc->Plume_Center, 1, scal->length));
		}
		else if(bc->Plume_Dimension == 2)
		{
			// 3D circular inflow a given [X,Y] coordinates
			PetscCall(getScalarParam(fb,_REQUIRED_, "Plume_Center", bc->Plume_Center, 2, scal->length));
		}

		// other options
		PetscCall(getIntParam	 (fb, _REQUIRED_, "Plume_Phase"       , &bc->Plume_Phase,       1, mID));
		PetscCall(getScalarParam(fb, _REQUIRED_, "Plume_Temperature" , &bc->Plume_Temperature, 1, 1));
		PetscCall(getScalarParam(fb, _REQUIRED_, "Plume_Radius",       &bc->Plume_Radius,      1, scal->length));

	}

	if((bc->bot_open || bc->Plume_Type == 2) && !bc->phase_inflow_bot )
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "The permeable inflow phase or the mantle plume phase must be defined\n");
	}

	//========================
	// TEMPERATURE CONSTRAINTS
	//========================

	bc->TbotNumPeriods = 1;
	PetscCall(getIntParam (fb, _OPTIONAL_, "temp_bot_num_periods",  &bc->TbotNumPeriods,  1, _max_periods_ ));

	if(bc->TbotNumPeriods > 1)
	{
		PetscCall(getScalarParam(fb, _REQUIRED_, "temp_bot_time_delim", bc->TbotTimeDelims, bc->TbotNumPeriods-1, scal->time	));
		PetscCall(getScalarParam(fb, _REQUIRED_, "temp_bot",            bc->Tbot,           bc->TbotNumPeriods,   1.0));
	}
	else
	{
		PetscCall(getScalarParam(fb, _OPTIONAL_, "temp_bot",  bc->Tbot, 1,   1.0));
	}
	PetscCall(getScalarParam(fb, _OPTIONAL_, "temp_top",   &bc->Ttop,     1, 1.0));
	PetscCall(getIntParam   (fb, _OPTIONAL_, "init_temp",  &bc->initTemp, 1, -1));

	//=====================
	// PRESSURE CONSTRAINTS
	//=====================

	PetscCall(getScalarParam(fb, _OPTIONAL_, "pres_bot",  &bc->pbot,     1, 1.0));
	PetscCall(getScalarParam(fb, _OPTIONAL_, "pres_top",  &bc->ptop,     1, 1.0));
	PetscCall(getIntParam   (fb, _OPTIONAL_, "init_pres", &bc->initPres, 1, -1));

	//======
	// CHECK
	//======

	if((bc->Tbot[0] == bc->Ttop) && bc->initTemp)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Top and bottom temperatures give zero initial gradient (Tbot, Ttop, initTemp) \n");
	}

	if(bc->top_open && bc->noslip[5])
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "No-slip condition is incompatible with open boundary (open_top_bound, noslip) \n");
	}

	if(periodic)
	{
		if(!bc->noslip[4])
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Periodic condition requires no-slip on bottom boundary (periodic, noslip) \n");
		}

		if(bc->noslip[0] || bc->noslip[1])
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Periodic condition is incompatible with no-slip on left and right boundaries (periodic, noslip) \n");
		}

		if(bc->ExxNumPeriods || bc->ExyNumPeriods)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Periodic condition is incompatible with xx and xy background strain rates (periodic, exx_num_periods, exy_num_periods) \n");
		}

		if(bc->face)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Periodic condition is incompatible with boundary velocity (periodic, bvel_face) \n");
		}

		if(bc->fixPhase != -1)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Periodic condition is incompatible with fixed phase (periodic, fix_phase) \n");
		}

		if(bc->fixCell)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Periodic condition is incompatible with fixed cells (periodic, fix_cell) \n");
		}
	}

	//==============
	// print summary
	//==============

	PetscPrintf(PETSC_COMM_WORLD, "Boundary condition parameters: \n");

	PetscPrintf(PETSC_COMM_WORLD, "   No-slip boundary mask [lt rt ft bk bm tp]  : ");

	if(periodic) { PetscPrintf(PETSC_COMM_WORLD, "   Periodic bc in x-direction                 @ "); }

	for(jj = 0; jj < 6; jj++)
	{
		PetscPrintf(PETSC_COMM_WORLD, "%" PetscInt_FMT " ", bc->noslip[jj]);
	}

	PetscPrintf(PETSC_COMM_WORLD, "\n");

	if(bc->ExxNumPeriods) { PetscPrintf(PETSC_COMM_WORLD, "   Number of x-background strain rate periods : %" PetscInt_FMT " \n",  bc->ExxNumPeriods); }
	if(bc->EyyNumPeriods) { PetscPrintf(PETSC_COMM_WORLD, "   Number of y-background strain rate periods : %" PetscInt_FMT " \n",  bc->EyyNumPeriods); }
	if(bc->nblocks)       { PetscPrintf(PETSC_COMM_WORLD, "   Number of Bezier blocks                    : %" PetscInt_FMT " \n",  bc->nblocks);       }

	for(jj = 0; jj < bc->nblocks; jj++)
	{
		PetscCall(BCBlockPrint(bc->blocks + jj, scal, jj));
	}

	if(bc->nboxes) { PetscPrintf(PETSC_COMM_WORLD, "   Number of velocity boxes                   : %" PetscInt_FMT " \n",  bc->nboxes); }

	for(jj = 0; jj < bc->nboxes; jj++)
	{
		PetscCall(VelBoxPrint(bc->vboxes + jj, scal, jj));
	}

	for(jj = 0; jj < bc->ncylinders; jj++)
	{
		PetscCall(VelCylinderPrint(bc->vcylinders + jj, scal, jj));
	}
	if(bc->top_open) { PetscPrintf(PETSC_COMM_WORLD, "   Open top boundary                          @ \n"); }
	if(bc->bot_open) { PetscPrintf(PETSC_COMM_WORLD, "   Open bottom boundary                       @ \n"); }
	if(bc->bot_open && bc->Plume_Type == 2) {;}
	if(bc->fixPhase != -1)   { PetscPrintf(PETSC_COMM_WORLD, "   Fixed phase                                : %" PetscInt_FMT "  \n", bc->fixPhase); }
	if(bc->Ttop     != -1.0) { PetscPrintf(PETSC_COMM_WORLD, "   Top boundary temperature                   : %g %s \n", bc->Ttop, scal->lbl_temperature); }
	if(bc->TbotNumPeriods == 1)
	{
		if(bc->Tbot[0] != -1.0) { PetscPrintf(PETSC_COMM_WORLD, "   Bottom boundary temperature                : %g %s \n", bc->Tbot[0], scal->lbl_temperature); }
	}
	else
	{
		// we have a Tbot that changes with time
		PetscPrintf(PETSC_COMM_WORLD, "   Number of bottom boundary temp periods     : %" PetscInt_FMT "  \n", bc->TbotNumPeriods);
		PetscPrintf(PETSC_COMM_WORLD, "   Bottom boundary temperatures               : ");
		for (jj=0; jj<bc->TbotNumPeriods; jj++)
		{
			PetscPrintf(PETSC_COMM_WORLD, "%g ", bc->Tbot[jj]);
		}
		PetscPrintf(PETSC_COMM_WORLD, " %s \n", scal->lbl_temperature);
		PetscPrintf(PETSC_COMM_WORLD, "   Bottom boundary temp time periods          :     ");
		for (jj=0; jj<bc->TbotNumPeriods-1; jj++)
		{
			PetscPrintf(PETSC_COMM_WORLD, "%g ", bc->TbotTimeDelims[jj]*scal->time);
		}
		PetscPrintf(PETSC_COMM_WORLD, " %s \n", scal->lbl_time);
	}

	if(bc->Plume_Inflow == 1)
	{
		PetscPrintf(PETSC_COMM_WORLD, "   Adding plume inflow bottom condition       @ \n");
		if(bc->Plume_Type == 1)
		{
			PetscPrintf(PETSC_COMM_WORLD, "      Type of plume                           : Inflow \n");
		}
		else
		{
			PetscPrintf(PETSC_COMM_WORLD, "      Type of plume                           : Open Bottom \n");
		}
		if(bc->Plume_VelocityType == 0)
		{
			PetscPrintf(PETSC_COMM_WORLD, "      Type of velocity perturbation           : Poiseuille flow (and constant outflow) \n");
		}
		else
		{
			PetscPrintf(PETSC_COMM_WORLD, "      Type of velocity perturbation           : Gaussian in/out flow \n");
		}

		PetscPrintf(PETSC_COMM_WORLD, "      Temperature of plume                    : %g %s \n", bc->Plume_Temperature, 	 				scal->lbl_temperature);
		PetscPrintf(PETSC_COMM_WORLD, "      Phase of plume                          : %" PetscInt_FMT " \n",  bc->Plume_Phase);
		PetscPrintf(PETSC_COMM_WORLD, "      Inflow velocity                         : %g %s \n", bc->Plume_Inflow_Velocity*scal->velocity, scal->lbl_velocity);
		PetscPrintf(PETSC_COMM_WORLD, "      Area fraction of plume                  : %g \n", bc->Plume_areaFrac);

		if(bc->Plume_Dimension == 1)
		{
			PetscPrintf(PETSC_COMM_WORLD, "      Location of center                      : [%g] %s \n", bc->Plume_Center[0]*scal->length,       scal->lbl_length);
		}
		else
		{
			PetscPrintf(PETSC_COMM_WORLD, "      Location of center                      : [%g, %g] %s \n", bc->Plume_Center[0]*scal->length, bc->Plume_Center[1]*scal->length, scal->lbl_length);
		}

		PetscPrintf(PETSC_COMM_WORLD, "      Radius of plume                         : %g %s \n", bc->Plume_Radius*scal->length, scal->lbl_length);
	}

	if(bc->face)
	{
		PetscPrintf(PETSC_COMM_WORLD, "   Adding inflow velocity at boundary         @ \n");

		if(bc->VelNumPeriods>1)
		{
			PetscPrintf(PETSC_COMM_WORLD, "      Number of inflow periods                : %" PetscInt_FMT "   \n",  bc->VelNumPeriods);
		}
		else
		{
			PetscPrintf(PETSC_COMM_WORLD, "      Number of inflow periods                : 1   \n");
		}

		if(bc->VelNetNumPeriods > 1)
		{
			PetscPrintf(PETSC_COMM_WORLD, "      Number of net inflow periods                : %" PetscInt_FMT "   \n",  bc->VelNetNumPeriods);
		}
		else
		{
			PetscPrintf(PETSC_COMM_WORLD, "      Number of net inflow periods                : 1   \n");
		}

		PetscPrintf(PETSC_COMM_WORLD, "      Inflow velocity boundary                : %s \n", str_inflow);

		if(bc->face_out == 1)
		{
			PetscPrintf(PETSC_COMM_WORLD, "      Outflow at opposite boundary            @ \n");
		}

		if(bc->num_phase_bc >= 0)
		{
			PetscPrintf(PETSC_COMM_WORLD, "      Inflow phase                            : %" PetscInt_FMT " \n", bc->phase[0]);
		}
		else
		{
			PetscPrintf(PETSC_COMM_WORLD, "      Inflow phase from next to boundary      @ \n");
		}

		PetscPrintf(PETSC_COMM_WORLD, "      Inflow window [bottom, top]             : [%3.2f,%3.2f] %s \n", bc->bot*scal->length, bc->top*scal->length, scal->lbl_length);
		PetscPrintf(PETSC_COMM_WORLD, "      Inflow velocity                         : %1.2f %s \n", bc->velin*scal->velocity, scal->lbl_velocity);

		if(bc->velout > 0.0)
		{
			PetscPrintf(PETSC_COMM_WORLD, "      Outflow velocity                        : %1.2f %s \n", bc->velout*scal->velocity, scal->lbl_velocity);
		}
		else if (!bc->face_out)
		{
			PetscPrintf(PETSC_COMM_WORLD, "       Outflow velocity from mass balance     @ \n");
		}
		if(bc->face == 5)                  { PetscPrintf(PETSC_COMM_WORLD, "      Bottom flow velocity                    : %1.2f %s \n", bc->velbot*scal->velocity, scal->lbl_velocity); }
		if(bc->face == 5 && !bc->top_open) { PetscPrintf(PETSC_COMM_WORLD, "      Top flow velocity                       : %1.2f %s \n", bc->veltop*scal->velocity, scal->lbl_velocity); }

		if(bc->relax_dist > 0){ PetscPrintf(PETSC_COMM_WORLD, "      Velocity smoothening distance           : %1.2f %s \n", bc->relax_dist*scal->length, scal->lbl_length); }

		if(bc->bvel_temperature_inflow)
		{
			if(bc->bvel_temperature_inflow == 1)
			{
				PetscPrintf(PETSC_COMM_WORLD, "      Temperature type of inflow material     : Constant \n");
				PetscPrintf(PETSC_COMM_WORLD, "         Temperature                          : %g %s  \n",bc->bvel_constant_temperature*scal->time,    scal->lbl_temperature);
			}
			if(bc->bvel_temperature_inflow == 2)
			{
				PetscPrintf(PETSC_COMM_WORLD, "      Temperature type of inflow material     : Halfspace cooling \n");
				PetscPrintf(PETSC_COMM_WORLD, "         Thermal Age                          : %1.0f %s  \n",bc->bvel_thermal_age*scal->time,    scal->lbl_time);
				PetscPrintf(PETSC_COMM_WORLD, "         Temperature @ top                    : %1.1f %s  \n",bc->bvel_temperature_top,           scal->lbl_temperature );
				PetscPrintf(PETSC_COMM_WORLD, "         Temperature @ bottom                 : %1.1f %s  \n",bc->bvel_potential_temperature,     scal->lbl_temperature );
			}
		}
		else
		{
			 PetscPrintf(PETSC_COMM_WORLD, "      Inflow temperature from closest marker  @ \n");
		}
	}

	if(bc->ptop != -1.0) PetscPrintf(PETSC_COMM_WORLD, "   Top boundary pressure                      : %g %s \n", bc->ptop, scal->lbl_stress);
	if(bc->pbot != -1.0) PetscPrintf(PETSC_COMM_WORLD, "   Bottom boundary pressure                   : %g %s \n", bc->pbot, scal->lbl_stress);

	PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------\n");

	// nondimensionalize temperature & pressure
	if(bc->Ttop    != -1.0) { bc->Ttop  = (bc->Ttop + scal->Tshift)/scal->temperature; }
	if(bc->Tbot[0] != -1.0)
	{
		for (jj = 0; jj < bc->TbotNumPeriods; jj++)
		{
			bc->Tbot[jj] = (bc->Tbot[jj]  + scal->Tshift)/scal->temperature;
		}
	}

	if(bc->ptop != -1.0)            { bc->ptop /= scal->stress; }
	if(bc->pbot != -1.0)            { bc->pbot /= scal->stress; }
	bc->Plume_Temperature          = (bc->Plume_Temperature          + scal->Tshift)/scal->temperature;
	bc->bvel_potential_temperature = (bc->bvel_potential_temperature + scal->Tshift)/scal->temperature;
	bc->bvel_temperature_top       = (bc->bvel_temperature_top       + scal->Tshift)/scal->temperature;
	bc->bvel_constant_temperature  = (bc->bvel_constant_temperature  + scal->Tshift)/scal->temperature;

	// allocate vectors and arrays
	PetscCall(BCCreateData(bc));

	// read fixed cells from files in parallel
	PetscCall(BCReadFixCell(bc, fb));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode BCReadRestart(BCCtx *bc, FILE *fp)
{
	PetscInt nCells;

	
	PetscFunctionBeginUser;

	nCells = bc->fs->nCells;

	// allocate memory
	PetscCall(BCCreateData(bc));

	// read fixed cell IDs
	if(bc->fixCell)
	{
		fread(bc->fixCellFlag, (size_t)nCells, 1, fp);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode BCWriteRestart(BCCtx *bc, FILE *fp)
{
	PetscInt nCells;

	PetscFunctionBeginUser;

	nCells = bc->fs->nCells;

	// write fixed cell IDs
	if(bc->fixCell)
	{
		fwrite(bc->fixCellFlag, (size_t)nCells, 1, fp);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode BCCreateData(BCCtx *bc)
{
	FDSTAG   *fs;
	DOFIndex *dof;

	
	PetscFunctionBeginUser;

	fs  =  bc->fs;
	dof = &fs->dof;

	// create boundary conditions vectors (velocity, pressure, temperature)
	PetscCall(DMCreateLocalVector(fs->DA_X,   &bc->bcvx));
	PetscCall(DMCreateLocalVector(fs->DA_Y,   &bc->bcvy));
	PetscCall(DMCreateLocalVector(fs->DA_Z,   &bc->bcvz));
	PetscCall(DMCreateLocalVector(fs->DA_CEN, &bc->bcp));
	PetscCall(DMCreateLocalVector(fs->DA_CEN, &bc->bcT));

	// SPC velocity-pressure
	PetscCall(makeIntArray (&bc->SPCList, NULL, dof->ln));
	PetscCall(makeScalArray(&bc->SPCVals, NULL, dof->ln));

	// SPC (temperature)
	PetscCall(makeIntArray (&bc->tSPCList, NULL, dof->lnp));
	PetscCall(makeScalArray(&bc->tSPCVals, NULL, dof->lnp));

	if(bc->fixCell)
	{
		PetscCall(PetscMalloc((size_t)fs->nCells, &bc->fixCellFlag));
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode BCDestroy(BCCtx *bc)
{
	
	PetscFunctionBeginUser;

	// destroy boundary conditions vectors (velocity, pressure, temperature)
	PetscCall(VecDestroy(&bc->bcvx));
	PetscCall(VecDestroy(&bc->bcvy));
	PetscCall(VecDestroy(&bc->bcvz));
	PetscCall(VecDestroy(&bc->bcp));
	PetscCall(VecDestroy(&bc->bcT));

	// SPC velocity-pressure
	PetscCall(PetscFree(bc->SPCList));
	PetscCall(PetscFree(bc->SPCVals));

	// SPC temperature
	PetscCall(PetscFree(bc->tSPCList));
	PetscCall(PetscFree(bc->tSPCVals));

	// fixed cell IDs
	PetscCall(PetscFree(bc->fixCellFlag));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode BCReadFixCell(BCCtx *bc, FB *fb)
{
	FILE           *fp;
	PetscLogDouble  t;
	PetscMPIInt     rank;
	struct          stat sb;
	char           *filename, file[_str_len_];

	
	PetscFunctionBeginUser;

	// check activation
	if(!bc->fixCell) PetscFunctionReturn(0);

	// get file name
	PetscCall(getStringParam(fb, _OPTIONAL_, "fix_cell_file", file, "./bc/cdb"));

	PrintStart(&t, "Loading fixed cell flags in parallel from", file);

	// compile input file name with extension
	PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));

	asprintf(&filename, "%s.%1.8" PetscMPIInt_FMT ".dat", file, rank);

	// open file
	fp = fopen(filename, "rb");

	if(fp == NULL)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Cannot open input file %s\n", filename);
	}

	// check file size
	stat(filename, &sb);

	if((PetscInt)sb.st_size != bc->fs->nCells)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Wrong fixed cell file size %s\n", filename);
	}

	// read flags
	fread(bc->fixCellFlag, (size_t)bc->fs->nCells, 1, fp);

	// close file
	fclose(fp);
	free(filename);

	PrintDone(t);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode BCApply(BCCtx *bc)
{
	FDSTAG *fs;

	
	PetscFunctionBeginUser;

	// access context
	fs = bc->fs;

	// mark all variables unconstrained
	PetscCall(VecSet(bc->bcvx, DBL_MAX));
	PetscCall(VecSet(bc->bcvy, DBL_MAX));
	PetscCall(VecSet(bc->bcvz, DBL_MAX));
	PetscCall(VecSet(bc->bcp,  DBL_MAX));
	PetscCall(VecSet(bc->bcT,  DBL_MAX));

	//============
	// TEMPERATURE
	//============

	// WARNING! Synchronization is necessary if SPC constraints are active
	// LOCAL_TO_LOCAL(fs->DA_CEN, bc->bcT)

	PetscCall(BCApplyTemp(bc));

	//==========================================
	// PRESSURE (must be called before velocity)
	//==========================================

	// WARNING! Synchronization is necessary if SPC constraints are active
	// LOCAL_TO_LOCAL(fs->DA_CEN, bc->bcp)

	PetscCall(BCApplyPres(bc));

	//=============================
	// VELOCITY (RESTRUCTURE THIS!)
	//=============================

	// apply default velocity constraints
	PetscCall(BCApplyVelDefault(bc));

	// apply Bezier block constraints
	PetscCall(BCApplyBezier(bc));

	// apply prescribed boundary velocity
	PetscCall(BCApplyBoundVel(bc));

	// apply velocity boxes
	PetscCall(BCApplyVelBox(bc));

	// apply velocity cylinders
	PetscCall(BCApplyVelCylinder(bc));

	// fix all cells occupied by phase
	PetscCall(BCApplyPhase(bc));

	// fix specific cells
	PetscCall(BCApplyCells(bc));

	// plume like boundary condition
	if(bc->Plume_Type == 1)
	{
		PetscCall(BC_Plume_inflow(bc));
	}

	// synchronize SPC constraints in the internal ghost points
	LOCAL_TO_LOCAL(fs->DA_X,   bc->bcvx)
	LOCAL_TO_LOCAL(fs->DA_Y,   bc->bcvy)
	LOCAL_TO_LOCAL(fs->DA_Z,   bc->bcvz)

	// apply two-point constraints
	PetscCall(BCApplyVelTPC(bc));

	// form SPC constraint lists
	PetscCall(BCListSPC(bc));

	// apply SPC to global solution vector
	PetscCall(BCApplySPC(bc));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode BCApplySPC(BCCtx *bc)
{
	// apply SPC to global solution vector

	PetscScalar *sol, *vals;
	PetscInt    i, num, *list;

	
	PetscFunctionBeginUser;

	PetscCall(VecGetArray(bc->jr->gsol, &sol));

	//============================================
	// enforce single point constraints (velocity)
	//============================================

	num   = bc->vNumSPC;
	list  = bc->vSPCList;
	vals  = bc->vSPCVals;

	for(i = 0; i < num; i++) sol[list[i]] = vals[i];

	//============================================
	// enforce single point constraints (pressure)
	//============================================

	num   = bc->pNumSPC;
	list  = bc->pSPCList;
	vals  = bc->pSPCVals;

	for(i = 0; i < num; i++) sol[list[i]] = vals[i];

	PetscCall(VecRestoreArray(bc->jr->gsol, &sol));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// Specific constraints
//---------------------------------------------------------------------------
PetscErrorCode BCApplyPres(BCCtx *bc)
{
	// apply pressure constraints

	FDSTAG      *fs;
	PetscScalar pbot, ptop;
	PetscInt    mcz;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***bcp;

	
	PetscFunctionBeginUser;

	// access context
	fs = bc->fs;

	// get boundary pressure
	pbot = bc->pbot;
	ptop = bc->ptop;

	// initialize index bounds
	mcz = fs->dsz.tcels - 1;

	PetscCall(DMDAVecGetArray(fs->DA_CEN, bc->bcp, &bcp));

	//-----------------------------------------------------
	// P points (TPC only, hence looping over ghost points)
	//-----------------------------------------------------
	if(pbot >= 0.0 || ptop >= 0.0)
	{
		GET_CELL_RANGE_GHOST_INT(nx, sx, fs->dsx)
		GET_CELL_RANGE_GHOST_INT(ny, sy, fs->dsy)
		GET_CELL_RANGE_GHOST_INT(nz, sz, fs->dsz)

		START_STD_LOOP
		{
			// only positive pressure!
			// negative will set normal velocity BC automatically
			if(pbot >= 0.0 && k == 0)   { bcp[k-1][j][i] = pbot; }
			if(ptop >= 0.0 && k == mcz) { bcp[k+1][j][i] = ptop; }
		}
		END_STD_LOOP
	}

	// restore access
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, bc->bcp, &bcp));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode BCApplyTemp(BCCtx *bc)
{
	// apply temperature constraints

	FDSTAG      *fs;
	PetscScalar Tbot, Ttop;
	PetscInt    mcz;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***bcT;
	PetscScalar x,y;
	PetscScalar xmin, xmax;
	PetscScalar rad_plume_squared, rad_squared;

	
	PetscFunctionBeginUser;

	// access context
	fs = bc->fs;

	// get boundary temperatures
	PetscCall(BCGetTempBound(bc, &Tbot));

	Ttop = bc->Ttop;

	// initialize index bounds
	mcz = fs->dsz.tcels - 1;

	PetscCall(DMDAVecGetArray(fs->DA_CEN, bc->bcT, &bcT));

	//-----------------------------------------------------
	// T points (TPC only, hence looping over ghost points)
	//-----------------------------------------------------
	if(Tbot >= 0.0 || Ttop >= 0.0)
	{
		GET_CELL_RANGE_GHOST_INT(nx, sx, fs->dsx)
		GET_CELL_RANGE_GHOST_INT(ny, sy, fs->dsy)
		GET_CELL_RANGE_GHOST_INT(nz, sz, fs->dsz)

		START_STD_LOOP
		{
			// only positive temperature!
			// negative will set zero-flux BC automatically
			if(Tbot >= 0.0 && k == 0)   { bcT[k-1][j][i] = Tbot; }
			if(Ttop >= 0.0 && k == mcz) { bcT[k+1][j][i] = Ttop; }

			// in case we have a plume-like inflow boundary condition:
			if(bc->Plume_Inflow == 1 && k==0)
			{
				x = COORD_CELL_GHOST(i, fs->dsx);
				y = COORD_CELL_GHOST(j, fs->dsy);

				// 2D plume
				if(bc->Plume_Dimension==1)
				{
					xmin =  bc->Plume_Center[0] - 3.0*bc->Plume_Radius;
					xmax =  bc->Plume_Center[0] + 3.0*bc->Plume_Radius;

					if((x >= xmin) && (x <= xmax))
					{
						bcT[k-1][j][i] = Tbot + (bc->Plume_Temperature-Tbot)*PetscExpScalar(-PetscPowScalar(x-bc->Plume_Center[0], 2.0)/(PetscPowScalar(bc->Plume_Radius, 2.0)));
					}

				}
				// 3D plume
				else
				{
					rad_plume_squared = PetscPowScalar(bc->Plume_Radius,2.0);
					rad_squared       = PetscPowScalar(x - bc->Plume_Center[0],2.0) + PetscPowScalar(y - bc->Plume_Center[1],2.0);

					if(rad_squared <= 9.0*rad_plume_squared)
					{
						bcT[k-1][j][i] = Tbot + (bc->Plume_Temperature-Tbot)*PetscExpScalar(-(rad_squared/(2.0*rad_plume_squared)));
					}
				}
			}
		}
		END_STD_LOOP
	}

	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, bc->bcT, &bcT));

PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode BCApplyVelDefault(BCCtx *bc)
{
	// apply default velocity constraints on the boundaries

	FDSTAG      *fs;
	PetscScalar Exx, Eyy, Ezz, Exy;
	PetscScalar Rxx, Ryy, Rzz;
	PetscScalar bx,  by,  bz;
	PetscScalar ex,  ey,  ez;
	PetscScalar vbx, vby, vbz;
	PetscScalar vex, vey, vez;
	PetscScalar y;
	PetscInt    periodic;
	PetscInt    mnx, mny, mnz;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, top_open;
	PetscScalar ***bcvx,  ***bcvy,  ***bcvz, ***bcp;

	
	PetscFunctionBeginUser;

	// access context
	fs = bc->fs;

	// set periodic flag
	periodic = fs->periodic;

	// set open boundary flag
	top_open = bc->top_open;

	// initialize index bounds
	mnx = fs->dsx.tnods - 1;
	mny = fs->dsy.tnods - 1;
	mnz = fs->dsz.tnods - 1;

	// get current coordinates of the mesh boundaries
	PetscCall(FDSTAGGetGlobalBox(fs, &bx, &by, &bz, &ex, &ey, &ez));

	// get background strain rates
	PetscCall(BCGetBGStrainRates(bc, &Exx, &Eyy, &Ezz, &Exy, &Rxx, &Ryy, &Rzz));

	// get boundary velocities
	// reference point is assumed to be fixed
	// velocity is a product of strain rate and coordinate w.r.t. reference point
	vbx = (bx - Rxx)*Exx;   vex = (ex - Rxx)*Exx;
	vby = (by - Ryy)*Eyy;   vey = (ey - Ryy)*Eyy;
	vbz = (bz - Rzz)*Ezz;   vez = (ez - Rzz)*Ezz;

	if(top_open)
	{
		vez = 0.0;
		vbz = 0.0;
	}

	// access constraint vectors
	PetscCall(DMDAVecGetArray(fs->DA_X,   bc->bcvx, &bcvx));
	PetscCall(DMDAVecGetArray(fs->DA_Y,   bc->bcvy, &bcvy));
	PetscCall(DMDAVecGetArray(fs->DA_Z,   bc->bcvz, &bcvz));
	PetscCall(DMDAVecGetArray(fs->DA_CEN, bc->bcp,  &bcp));

	//=========================================================================
	// SPC (normal velocities)
	//=========================================================================

	//------------------
	// X points SPC only
	//------------------
	if(!periodic)
	{
		PetscCall(DMDAGetCorners(fs->DA_X, &sx, &sy, &sz, &nx, &ny, &nz));

		START_STD_LOOP
		{
			// get coordinate
			y = COORD_CELL(j, sy, fs->dsy);

			if(i == 0   && bcp[k][j][-1 ] == DBL_MAX) { bcvx[k][j][i] = vbx + (y-Ryy)*Exy; }
			if(i == mnx && bcp[k][j][mnx] == DBL_MAX) { bcvx[k][j][i] = vex + (y-Ryy)*Exy; }

		}
		END_STD_LOOP
	}

	//------------------
	// Y points SPC only
	//------------------
	PetscCall(DMDAGetCorners(fs->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		if(j == 0   && bcp[k][-1 ][i] == DBL_MAX) { bcvy[k][j][i] = vby; }
		if(j == mny && bcp[k][mny][i] == DBL_MAX) { bcvy[k][j][i] = vey; }
	}
	END_STD_LOOP

	//------------------
	// Z points SPC only
	//------------------
	PetscCall(DMDAGetCorners(fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		if(k == 0                && bcp[-1 ][j][i] == DBL_MAX) { bcvz[k][j][i] = vbz; }
		if(k == mnz && !top_open && bcp[mnz][j][i] == DBL_MAX) { bcvz[k][j][i] = vez; }
	}
	END_STD_LOOP

	// restore access
	PetscCall(DMDAVecRestoreArray(fs->DA_X,   bc->bcvx, &bcvx));
	PetscCall(DMDAVecRestoreArray(fs->DA_Y,   bc->bcvy, &bcvy));
	PetscCall(DMDAVecRestoreArray(fs->DA_Z,   bc->bcvz, &bcvz));
	PetscCall(DMDAVecRestoreArray(fs->DA_CEN, bc->bcp,  &bcp));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode BCGetVelins(BCCtx *bc)
{
	Scaling *scal = bc->scal;
	PetscScalar  bz;
	PetscInt    jj, kk;
	PetscScalar time;

	
	PetscFunctionBegin;

	// initialize
	time = bc->ts->time;

	if (bc->VelNumPeriods > 1 && bc->VelNetNumPeriods > 1)
	{
		// both velin and velin_net have period arrays
		for (jj = 0; jj < bc->VelNumPeriods-1; jj++)
		{
			if (time < bc->VelTimeDelims[jj]) break;
		}
		for (kk = 0; kk < bc->VelNetNumPeriods-1; kk++)
		{
			if (time < bc->VelNetTimeDelims[kk]) break;
		}

		PetscCall(FDSTAGGetGlobalBox(bc->fs, NULL, NULL, &bz, NULL, NULL, NULL));

		bc->velin  = bc->velin_array[jj] + bc->velin_net_array[kk];
		bc->velout = -bc->velin*(bc->top - bc->bot)/(bc->bot - bz);

		PetscPrintf(PETSC_COMM_WORLD,
				"BCGetVelins BOTH: time=%g (Myr) jj=%" PetscInt_FMT " kk=%" PetscInt_FMT " velin_base=%g velin_net=%g velin_total=%g\n",
				(PetscScalar)(time*scal->time),
				jj, kk,
				(PetscScalar)(bc->velin_array[jj]*scal->velocity),
				(PetscScalar)(bc->velin_net_array[kk]*scal->velocity),
				(PetscScalar)((bc->velin_array[jj] + bc->velin_net_array[kk])*scal->velocity));
	}
	else if (bc->VelNumPeriods)
	{
		// only velin array provided
		for(jj = 0; jj < bc->VelNumPeriods-1; jj++)
		{
			if (time < bc->VelTimeDelims[jj]) break;
		}
		PetscCall(FDSTAGGetGlobalBox(bc->fs, NULL, NULL, &bz, NULL, NULL, NULL));
		bc->velin  = bc->velin_array[jj];
		bc->velout = -bc->velin*(bc->top - bc->bot)/(bc->bot - bz);
	}
	else if (bc->VelNetNumPeriods > 1)
	{
		// only velin_net array provided; velin scalar should have been read earlier */
		static PetscBool   baseVelinInitialized = PETSC_FALSE;
		static PetscScalar baseVelin            = 0.0;
		if(!baseVelinInitialized)
		{
			baseVelin            = bc->velin;  // store the original scalar inflow velocity
			baseVelinInitialized = PETSC_TRUE;
		}
		for(kk = 0; kk < bc->VelNetNumPeriods-1; kk++)
		{
			if (time < bc->VelNetTimeDelims[kk]) break;
		}
		PetscCall(FDSTAGGetGlobalBox(bc->fs, NULL, NULL, &bz, NULL, NULL, NULL));
		// piecewise constant: base scalar inflow plus the current net offset;
		// no accumulation over timesteps
		bc->velin  = baseVelin + bc->velin_net_array[kk];
		bc->velout = -bc->velin*(bc->top - bc->bot)/(bc->bot - bz);
	}
	else
	{
		// neither periods array provided; use scalar velin (already read) and compute velout
		PetscCall(FDSTAGGetGlobalBox(bc->fs, NULL, NULL, &bz, NULL, NULL, NULL));
		bc->velout = -bc->velin*(bc->top - bc->bot)/(bc->bot - bz);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode BCApplyVelTPC(BCCtx *bc)
{
	// apply two-point constraints on the boundaries

	FDSTAG      *fs;
	PetscInt    mcx, mcy, mcz;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscInt    nsLeft, nsRight, nsFront, nsBack, nsBottom, nsTop;
	PetscScalar Exy, Ryy, y, dy, vx;
	PetscScalar ***bcvx,  ***bcvy,  ***bcvz;

	
	PetscFunctionBeginUser;

	// access context
	fs = bc->fs;

	// initialize index bounds
	mcx = fs->dsx.tcels - 1;
	mcy = fs->dsy.tcels - 1;
	mcz = fs->dsz.tcels - 1;

	// initialize no-slip flags
	nsLeft   = bc->noslip[0];
	nsRight  = bc->noslip[1];
	nsFront  = bc->noslip[2];
	nsBack   = bc->noslip[3];
	nsBottom = bc->noslip[4];
	nsTop    = bc->noslip[5];

	// get background shear strain rate
	PetscCall(BCGetBGStrainRates(bc, NULL, NULL, NULL, &Exy, NULL, &Ryy, NULL));

	//=========================================================================
	// TPC (no-slip boundary conditions)
	//=========================================================================

	// access constraint vectors
	PetscCall(DMDAVecGetArray(fs->DA_X,   bc->bcvx, &bcvx));
	PetscCall(DMDAVecGetArray(fs->DA_Y,   bc->bcvy, &bcvy));
	PetscCall(DMDAVecGetArray(fs->DA_Z,   bc->bcvz, &bcvz));

	//-----------------------------------------------------
	// X points (TPC only, hence looping over ghost points)
	//-----------------------------------------------------
	if(nsFront || nsBack || nsBottom || nsTop || Exy)
	{
		GET_NODE_RANGE_GHOST_INT(nx, sx, fs->dsx)
		GET_CELL_RANGE_GHOST_INT(ny, sy, fs->dsy)
		GET_CELL_RANGE_GHOST_INT(nz, sz, fs->dsz)

		START_STD_LOOP
		{
			// set tangential velocity
			if(Exy)
			{
				// get coordinate and cell size
				y  = COORD_CELL_GHOST(j, fs->dsy);
				dy = SIZE_CELL(j, sy, fs->dsy);

				// compute velocity
				vx = (y-Ryy)*Exy;

				if(j == 0)   { bcvx[k][j-1][i] = vx - dy*Exy/2.0; }
				if(j == mcy) { bcvx[k][j+1][i] = vx + dy*Exy/2.0; }
			}
			else
			{
				if(nsFront  && j == 0)   { bcvx[k][j-1][i] = 0.0; }
				if(nsBack   && j == mcy) { bcvx[k][j+1][i] = 0.0; }
				if(nsBottom && k == 0)   { bcvx[k-1][j][i] = 0.0; }
				if(nsTop    && k == mcz) { bcvx[k+1][j][i] = 0.0; }
			}
		}
		END_STD_LOOP
	}

	//-----------------------------------------------------
	// Y points (TPC only, hence looping over ghost points)
	//-----------------------------------------------------
	if(nsLeft || nsRight || nsBottom || nsTop)
	{
		GET_CELL_RANGE_GHOST_INT(nx, sx, fs->dsx)
		GET_NODE_RANGE_GHOST_INT(ny, sy, fs->dsy)
		GET_CELL_RANGE_GHOST_INT(nz, sz, fs->dsz)

		START_STD_LOOP
		{
			if(nsLeft   && i == 0)   { bcvy[k][j][i-1] = 0.0; }
			if(nsRight  && i == mcx) { bcvy[k][j][i+1] = 0.0; }
			if(nsBottom && k == 0)   { bcvy[k-1][j][i] = 0.0; }
			if(nsTop    && k == mcz) { bcvy[k+1][j][i] = 0.0; }
		}
		END_STD_LOOP
	}

	//-----------------------------------------------------
	// Z points (TPC only, hence looping over ghost points)
	//-----------------------------------------------------
	if(nsLeft || nsRight || nsFront || nsBack)
	{
		GET_CELL_RANGE_GHOST_INT(nx, sx, fs->dsx)
		GET_CELL_RANGE_GHOST_INT(ny, sy, fs->dsy)
		GET_NODE_RANGE_GHOST_INT(nz, sz, fs->dsz)

		START_STD_LOOP
		{
			if(nsLeft  && i == 0)   { bcvz[k][j][i-1] = 0.0; }
			if(nsRight && i == mcx) { bcvz[k][j][i+1] = 0.0; }
			if(nsFront && j == 0)   { bcvz[k][j-1][i] = 0.0; }
			if(nsBack  && j == mcy) { bcvz[k][j+1][i] = 0.0; }
		}
		END_STD_LOOP
	}

	// restore access
	PetscCall(DMDAVecRestoreArray(fs->DA_X,   bc->bcvx, &bcvx));
	PetscCall(DMDAVecRestoreArray(fs->DA_Y,   bc->bcvy, &bcvy));
	PetscCall(DMDAVecRestoreArray(fs->DA_Z,   bc->bcvz, &bcvz));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode BCApplyBezier(BCCtx *bc)
{
	FDSTAG      *fs;
	BCBlock     *bcb;
	PetscInt    fbeg, fend, npoly, in, dim;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, ib;
	PetscScalar ***bcvx, ***bcvy, ***bcvz;
	PetscScalar t, dt, theta, costh, sinth, atol, bot, top, vel, zOffset, velz;
	PetscScalar Xbeg[4], Xend[4], xbeg[3], xend[3], box[4], cpoly[2*_max_poly_points_];

	
	PetscFunctionBeginUser;

	// check whether constraint is activated
	if(!bc->nblocks) PetscFunctionReturn(0);

	// access context
	fs    =  bc->fs;
	t     =  bc->ts->time;
	dt    =  bc->ts->dt;

	// access velocity constraint vectors
	PetscCall(DMDAVecGetArray(fs->DA_X, bc->bcvx, &bcvx));
	PetscCall(DMDAVecGetArray(fs->DA_Y, bc->bcvy, &bcvy));
	PetscCall(DMDAVecGetArray(fs->DA_Z, bc->bcvz, &bcvz));

	// loop over all bezier blocks
	for(ib = 0; ib < bc->nblocks; ib++)
	{
		bcb   =  bc->blocks + ib;
		dim   =  bcb->pathDim;
		npoly =  bcb->npoly;

		// get polygon positions in the beginning & end of the time step
		PetscCall(BCBlockGetPosition(bcb, t,    &fbeg, Xbeg));
		PetscCall(BCBlockGetPosition(bcb, t+dt, &fend, Xend));

		// check whether constraint applies to the current time step
		if(!fbeg || !fend) continue;

		// compute bot/top z-coordinates
		// For 2D path: bot/top are absolute z-coordinates (no z-movement)
		// For 3D path: bot/top are initial absolute z-coordinates that move with the path
		//              (same as polygon x-y vertices which are also initial absolute)
		if(dim == 3)
		{
			// compute z-displacement from initial path position (same logic as polygon x-y)
			zOffset = Xbeg[2] - bcb->path[2];  // current_path_z - initial_path_z
			bot     = bcb->bot + zOffset;
			top     = bcb->top + zOffset;
			velz    = (Xend[2] - Xbeg[2])/dt;  // z-velocity
		}
		else
		{
			bot  = bcb->bot;
			top  = bcb->top;
			velz = 0.0;
		}

		// get current polygon geometry
		PetscCall(BCBlockGetPolygon(bcb, Xbeg, cpoly));

		// get bounding box
		polygon_box(&npoly, cpoly, 1e-12, &atol, box);

		// get time step rotation matrix
		// theta is at index 2 for 2D path, index 3 for 3D path
		if(dim == 3)
		{
			theta = Xend[3] - Xbeg[3];
		}
		else
		{
			theta = Xend[2] - Xbeg[2];
		}
		costh = cos(theta);
		sinth = sin(theta);

		//---------
		// X points
		//---------
		PetscCall(DMDAGetCorners(fs->DA_X, &sx, &sy, &sz, &nx, &ny, &nz));

		START_STD_LOOP
		{
			// get node coordinates in the beginning of time step
			xbeg[0] = COORD_NODE(i, sx, fs->dsx);
			xbeg[1] = COORD_CELL(j, sy, fs->dsy);
			xbeg[2] = COORD_CELL(k, sz, fs->dsz);

			// perform point test
			if(xbeg[2] >= bot && xbeg[2] <= top)
			{
				in_polygon(1, xbeg, npoly, cpoly, box, atol, &in);

				// check whether point is inside polygon
				if(in)
				{
					// compute point position in the end of time step
					RotDispPoint2D(Xbeg, Xend, costh, sinth, xbeg, xend);

					// compute & set x-velocity
					vel = (xend[0] - xbeg[0])/dt;

					bcvx[k][j][i] = vel;
				}
			}
		}
		END_STD_LOOP

		//---------
		// Y points
		//---------
		PetscCall(DMDAGetCorners(fs->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz));

		START_STD_LOOP
		{
			// get node coordinates in the beginning of time step
			xbeg[0] = COORD_CELL(i, sx, fs->dsx);
			xbeg[1] = COORD_NODE(j, sy, fs->dsy);
			xbeg[2] = COORD_CELL(k, sz, fs->dsz);

			// perform point test
			if(xbeg[2] >= bot && xbeg[2] <= top)
			{
				in_polygon(1, xbeg, npoly, cpoly, box, atol, &in);

				// check whether point is inside polygon
				if(in)
				{
					// compute point position in the end of time step
					RotDispPoint2D(Xbeg, Xend, costh, sinth, xbeg, xend);

					// compute & set y-velocity
					vel = (xend[1] - xbeg[1])/dt;

					bcvy[k][j][i] = vel;
				}
			}
		}
		END_STD_LOOP

		//---------
		// Z points (only for 3D path)
		//---------
		if(dim == 3)
		{
			PetscCall(DMDAGetCorners(fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz));

			START_STD_LOOP
			{
				// get node coordinates in the beginning of time step
				xbeg[0] = COORD_CELL(i, sx, fs->dsx);
				xbeg[1] = COORD_CELL(j, sy, fs->dsy);
				xbeg[2] = COORD_NODE(k, sz, fs->dsz);

				// perform point test
				if(xbeg[2] >= bot && xbeg[2] <= top)
				{
					in_polygon(1, xbeg, npoly, cpoly, box, atol, &in);

					// check whether point is inside polygon
					if(in)
					{
						// set z-velocity (uniform for the block)
						bcvz[k][j][i] = velz;
					}
				}
			}
			END_STD_LOOP
		}
	}
	// restore access
	PetscCall(DMDAVecRestoreArray(fs->DA_X, bc->bcvx, &bcvx));
	PetscCall(DMDAVecRestoreArray(fs->DA_Y, bc->bcvy, &bcvy));
	PetscCall(DMDAVecRestoreArray(fs->DA_Z, bc->bcvz, &bcvz));

PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode BCApplyBoundVel(BCCtx *bc)
{
	FDSTAG      *fs;
	PetscInt    mnz, mnx, mny;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, kk;
	PetscInt    top_open, bot_open;
	PetscScalar ***bcvx,  ***bcvy, ***bcvz;
	PetscScalar z, bot, top, vel, velin, velout,relax_dist, velbot, veltop, time;

	
	PetscFunctionBeginUser;

	// check whether constraint is activated
	if(!bc->face) PetscFunctionReturn(0);

	// update inflow velocity value for current timestep
	PetscCall(BCGetVelins(bc));

	// access context
	fs         = bc->fs;
	bot        = bc->bot;
	top        = bc->top;
	velin      = bc->velin;
	velout     = bc->velout;
	relax_dist = bc->relax_dist;
	velbot     = bc->velbot;
	veltop     = bc->veltop;
	time       = bc->ts->time;
	top_open   = bc->top_open;
	bot_open   = bc->bot_open;

	// precompute net inflow value once for this timestep to avoid repeated inner-loop work */ // pkongpet 11/12/25
	PetscScalar velin_net = 0.0;
	if (bc->VelNetNumPeriods > 1)
	{
		for (kk = 0; kk < bc->VelNetNumPeriods-1; kk++)
		{
			if (time < bc->VelNetTimeDelims[kk]) break;
		}
		velin_net = bc->velin_net_array[kk]; // nondimensionalize
	}

	// initialize maximal index in all directions
	mnx = fs->dsx.tnods - 1;
	mny = fs->dsy.tnods - 1;
	mnz = fs->dsz.tnods - 1;

	// access velocity constraint vectors
	PetscCall(DMDAVecGetArray(fs->DA_X, bc->bcvx, &bcvx));
	PetscCall(DMDAVecGetArray(fs->DA_Y, bc->bcvy, &bcvy));
	PetscCall(DMDAVecGetArray(fs->DA_Z, bc->bcvz, &bcvz));

	//---------
	// X points
	//---------
	PetscCall(DMDAGetCorners(fs->DA_X, &sx, &sy, &sz, &nx, &ny, &nz));
	
	if(bc->face == 1 || bc->face == 2)
	{
		START_STD_LOOP
		{
			z   = COORD_CELL(k, sz, fs->dsz);
			vel = 0.0;
			if(bc->face_out)
			{
				if(z <= top && z >= bot)           vel = velin;
				if(z >= top && z<= top+relax_dist) vel = velin-(velin/(relax_dist))*(z-top);
				if(z <= bot && z>= bot-relax_dist) vel = velin+(velin/(relax_dist))*(z-bot);
				if(bc->face_out != 1)
				{
					if(z < bot-relax_dist)             vel = velout;
				}

				if((bc->face  == 1) && i == 0)                         { bcvx[k][j][i] =  vel; }
				if((bc->face  == 1) && i == mnx && bc->face_out == 1)  { bcvx[k][j][i] =  vel; }

				if((bc->face  == 2) && i == 0   && bc->face_out == 1)  { bcvx[k][j][i] = -vel; }
				if((bc->face  == 2) && i == mnx)                       { bcvx[k][j][i] = -vel; }

				if((bc->face  == 1) && i == 0   && bc->face_out == -1) { bcvx[k][j][i] =  vel; }
				if((bc->face  == 1) && i == mnx && bc->face_out == -1) { bcvx[k][j][i] = -vel; }

				if((bc->face  == 2) && i == 0   && bc->face_out == -1) { bcvx[k][j][i] =  vel; }
				if((bc->face  == 2) && i == mnx && bc->face_out == -1) { bcvx[k][j][i] = -vel; }

			}
			else
			{
				if(z <= top && z >= bot) vel = velin;
				if(z < bot)              vel = velout;

				if((bc->face == 1)  && i == 0 )   { bcvx[k][j][i] = vel; }
				if((bc->face == 2) && i == mnx) { bcvx[k][j][i] = vel; }
			}
		}
		END_STD_LOOP
	}

	if(bc->face == 5)
	{
		START_STD_LOOP
		{
			z   = COORD_CELL(k, sz, fs->dsz);
			vel = 0.0;
			if(z <= top && z >= bot) vel = velin;

			if(i == 0)   { bcvx[k][j][i] = vel; }
			if(i == mnx)
			{
				bcvx[k][j][i] = -vel + (2.0*velin_net); // right boundary (compensating inflow)
			}
		}
		END_STD_LOOP
	}

	//---------
	// Y points
	//---------
	PetscCall(DMDAGetCorners(fs->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz));

	if(bc->face == 3 || bc->face == 4)
	{
		START_STD_LOOP
		{
			z   = COORD_CELL(k, sz, fs->dsz);
			vel = 0.0;
			vel = 0.0;
			if(bc->face_out)
			{
				if(z <= top && z >= bot) vel = velin;
				if(z >= top && z<= top+relax_dist) vel = velin-(velin/(relax_dist))*(z-top);
				if(z <= bot && z>= bot-relax_dist) vel = velin+(velin/(relax_dist))*(z-bot);
				if(bc->face_out != 1)
				{
					if(z < bot-relax_dist)         vel = velout;
				}

				if((bc->face  == 3) && j == 0)                         { bcvy[k][j][i] =  vel; }
				if((bc->face  == 3) && j == mny && bc->face_out == 1)  { bcvy[k][j][i] =  vel; }

				if((bc->face  == 4) && j == 0   && bc->face_out == 1)  { bcvy[k][j][i] = -vel; }
				if((bc->face  == 4) && j == mny)                       { bcvy[k][j][i] = -vel; }


				if((bc->face  == 3) && j == 0   && bc->face_out == -1) { bcvy[k][j][i] =  vel; }
				if((bc->face  == 3) && j == mny && bc->face_out == -1) { bcvy[k][j][i] = -vel; }

				if((bc->face  == 4) && j == 0   && bc->face_out == -1) { bcvy[k][j][i] =  vel; }
				if((bc->face  == 4) && j == mny && bc->face_out == -1) { bcvy[k][j][i] = -vel; }
			}
			else
			{
				if(z <= top && z >= bot) vel = velin;
				if(z < bot)              vel = velout;

				if(bc->face == 3 && j == 0)   { bcvy[k][j][i] = vel; }
				if(bc->face == 4 && j == mny) { bcvy[k][j][i] = vel; }
			}
		}
		END_STD_LOOP
	}

	//---------
	// Z points
	//---------
	PetscCall(DMDAGetCorners(fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz));

	if(bc->face == 5 )
	{
		START_STD_LOOP
		{
			vel = 0.0;

			if(k == 0   && !bot_open)  vel = velbot;
			if(k == mnz && !top_open)  vel = veltop;

			if(k == 0   && !bot_open) { bcvz[k][j][i] = vel; }
			if(k == mnz && !top_open) { bcvz[k][j][i] = vel; }
		}
		END_STD_LOOP
	}

	// restore access
	PetscCall(DMDAVecRestoreArray(fs->DA_X, bc->bcvx, &bcvx));
	PetscCall(DMDAVecRestoreArray(fs->DA_Y, bc->bcvy, &bcvy));
	PetscCall(DMDAVecRestoreArray(fs->DA_Z, bc->bcvz, &bcvz));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode BCApplyVelBox(BCCtx *bc)
{
	FDSTAG      *fs;
	VelBox      *velbox;
	PetscScalar ***bcvx, ***bcvy, ***bcvz;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, ib;
	PetscScalar x, y, z, cx, cy, cz, dx, dy, dz, t, vx, vy, vz;
	PetscScalar xmin, xmax, ymin, ymax, zmin, zmax;

	
	PetscFunctionBeginUser;

	// skip initial guess
	if(bc->jr->ctrl.initGuess) PetscFunctionReturn(0);

	// check whether internal velocity box condition is activated
	if(!bc->nboxes) PetscFunctionReturn(0);

	// access context
	fs    =  bc->fs;
	t     =  bc->ts->time;

	// access velocity constraint vectors
	PetscCall(DMDAVecGetArray(fs->DA_X, bc->bcvx, &bcvx));
	PetscCall(DMDAVecGetArray(fs->DA_Y, bc->bcvy, &bcvy));
	PetscCall(DMDAVecGetArray(fs->DA_Z, bc->bcvz, &bcvz));

	// loop over all boxes
	for(ib = 0; ib < bc->nboxes; ib++)
	{
		// get current box
		velbox = bc->vboxes + ib;

		vx = velbox->vx; cx = velbox->cenX; dx = velbox->widthX;
		vy = velbox->vy; cy = velbox->cenY; dy = velbox->widthY;
		vz = velbox->vz; cz = velbox->cenZ; dz = velbox->widthZ;

		// advect box (if requested)
		if(velbox->advect)
		{
			if(vx != DBL_MAX) cx += vx*t;
			if(vy != DBL_MAX) cy += vy*t;
			if(vz != DBL_MAX) cz += vz*t;
		}

		// get bounds
		xmin = cx - dx/2.0; xmax = cx + dx/2.0;
		ymin = cy - dy/2.0; ymax = cy + dy/2.0;
		zmin = cz - dz/2.0; zmax = cz + dz/2.0;

		//---------
		// X points
		//----------
		if(vx != DBL_MAX)
		{
			PetscCall(DMDAGetCorners(fs->DA_X, &sx, &sy, &sz, &nx, &ny, &nz));

			START_STD_LOOP
			{
				// get coordinates
				x = COORD_NODE(i, sx, fs->dsx);
				y = COORD_CELL(j, sy, fs->dsy);
				z = COORD_CELL(k, sz, fs->dsz);

				if(x >= xmin && x <= xmax
						&& y >= ymin && y <= ymax
						&& z >= zmin && z <= zmax)
				{
					bcvx[k][j][i] = vx;
				}
			}
			END_STD_LOOP
		}

		//---------
		// Y points
		//----------
		if(vy != DBL_MAX)
		{
			PetscCall(DMDAGetCorners(fs->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz));

			START_STD_LOOP
			{
				// get coordinates
				x = COORD_CELL(i, sx, fs->dsx);
				y = COORD_NODE(j, sy, fs->dsy);
				z = COORD_CELL(k, sz, fs->dsz);
				if(x >= xmin && x <= xmax
						&& y >= ymin && y <= ymax
						&& z >= zmin && z <= zmax)
				{
					bcvy[k][j][i] = vy;
				}
			}
			END_STD_LOOP
		}

		//---------
		// Z points
		//----------
		if(vz != DBL_MAX)
		{
			PetscCall(DMDAGetCorners(fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz));

			START_STD_LOOP
			{
				// get coordinates
				x = COORD_CELL(i, sx, fs->dsx);
				y = COORD_CELL(j, sy, fs->dsy);
				z = COORD_NODE(k, sz, fs->dsz);
				if(x >= xmin && x <= xmax
						&& y >= ymin && y <= ymax
						&& z >= zmin && z <= zmax)
				{
					bcvz[k][j][i] = vz;
				}
			}
			END_STD_LOOP
		}
	}

	// restore access
	PetscCall(DMDAVecRestoreArray(fs->DA_X, bc->bcvx, &bcvx));
	PetscCall(DMDAVecRestoreArray(fs->DA_Y, bc->bcvy, &bcvy));
	PetscCall(DMDAVecRestoreArray(fs->DA_Z, bc->bcvz, &bcvz));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode BCApplyVelCylinder(BCCtx *bc)
{
	FDSTAG      *fs;
	VelCylinder *velcyl;
	PetscScalar ***bcvx, ***bcvy, ***bcvz;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, ic;
	PetscScalar x, y, z, cx, cy, cz, bx, by, bz, t, r, vx, vy, vz, vmag;
	PetscScalar a, ax, ay, az, px, py, pz, npc, dx, dy, dz, dr, rr;
	PetscScalar velType;

	
	PetscFunctionBeginUser;

	// skip initial guess
	if(bc->jr->ctrl.initGuess) PetscFunctionReturn(0);

	// check whether internal velocity cylinder condition is activated
	if(!bc->ncylinders) PetscFunctionReturn(0);

	// access context
	fs    =  bc->fs;
	t     =  bc->ts->time;

	// access velocity constraint vectors
	PetscCall(DMDAVecGetArray(fs->DA_X, bc->bcvx, &bcvx));
	PetscCall(DMDAVecGetArray(fs->DA_Y, bc->bcvy, &bcvy));
	PetscCall(DMDAVecGetArray(fs->DA_Z, bc->bcvz, &bcvz));

	// loop over all cylinders
	for(ic = 0; ic < bc->ncylinders; ic++)
	{
		// get current cylinder
		velcyl = bc->vcylinders + ic;

		// get coordinates
		by = velcyl->baseY; cy = velcyl->capY;
		bx = velcyl->baseX; cx = velcyl->capX;
		bz = velcyl->baseZ; cz = velcyl->capZ;
		r  = velcyl->rad;

		// get type of velocity profile
		velType = (PetscScalar)velcyl->type;

		// get velocity components
		vmag = velcyl->vmag;
		if(vmag != DBL_MAX)
		{
			// get cylinder axis vector
			ax = cx - bx;
			ay = cy - by;
			az = cz - bz;
			a  = sqrt(ax*ax + ay*ay + az*az);

			// partition velocities
			vx = vmag * ax / a;
			vy = vmag * ay / a;
			vz = vmag * az / a;

		}
		else
		{
			vy = velcyl->vy;
			vx = velcyl->vx;
			vz = velcyl->vz;
		}

		// advect cylinder (if requested)
		if(velcyl->advect)
		{
			if(vx != DBL_MAX) {bx += vx*t; cx += vx*t;}
			if(vy != DBL_MAX) {by += vy*t; cy += vy*t;}
			if(vz != DBL_MAX) {bz += vz*t; cz += vz*t;}
		}

		// get cylinder axis vector
		ax = cx - bx;
		ay = cy - by;
		az = cz - bz;

		//---------
		// X points
		//----------
		if(vx != DBL_MAX)
		{
			PetscCall(DMDAGetCorners(fs->DA_X, &sx, &sy, &sz, &nx, &ny, &nz));

			START_STD_LOOP
			{
				// get coordinates
				x = COORD_NODE(i, sx, fs->dsx);
				y = COORD_CELL(j, sy, fs->dsy);
				z = COORD_CELL(k, sz, fs->dsz);

				// get vector between a test point and cylinder base
				px = x - bx;
				py = y - by;
				pz = z - bz;

				// find normalized parametric coordinate of a point-axis projection
				npc = (ax*px + ay*py + az*pz)/(ax*ax + ay*ay + az*az);

				// find distance vector between point and axis
				dx = px - npc*ax;
				dy = py - npc*ay;
				dz = pz - npc*az;

				// compare position to radius
				dr = sqrt(dx*dx + dy*dy + dz*dz);
				rr = dr / r;

				// check cylinder
				if(npc >= 0.0 && npc <= 1.0 && rr <= 1.0)
				{
					bcvx[k][j][i] = vx * (1 - rr*rr*velType);
				}
			}
			END_STD_LOOP
		}

		//---------
		// Y points
		//----------
		if(vy != DBL_MAX)
		{
			PetscCall(DMDAGetCorners(fs->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz));

			START_STD_LOOP
			{
				// get coordinates
				x = COORD_CELL(i, sx, fs->dsx);
				y = COORD_NODE(j, sy, fs->dsy);
				z = COORD_CELL(k, sz, fs->dsz);

				// get vector between a test point and cylinder base
				px = x - bx;
				py = y - by;
				pz = z - bz;

				// find normalized parametric coordinate of a point-axis projection
				npc = (ax*px + ay*py + az*pz)/(ax*ax + ay*ay + az*az);

				// find distance vector between point and axis
				dx = px - npc*ax;
				dy = py - npc*ay;
				dz = pz - npc*az;

				// compare position to radius
				dr = sqrt(dx*dx + dy*dy + dz*dz);
				rr = dr / r;

				// check cylinder
				if(npc >= 0.0 && npc <= 1.0 && rr <= 1.0)
				{
					bcvy[k][j][i] = vy * (1 - rr*rr*velType);
				}
			}
			END_STD_LOOP
		}

		//---------
		// Z points
		//----------
		if(vz != DBL_MAX)
		{
			PetscCall(DMDAGetCorners(fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz));

			START_STD_LOOP
			{
				// get coordinates
				x = COORD_CELL(i, sx, fs->dsx);
				y = COORD_CELL(j, sy, fs->dsy);
				z = COORD_NODE(k, sz, fs->dsz);

				// get vector between a test point and cylinder base
				px = x - bx;
				py = y - by;
				pz = z - bz;

				// find normalized parametric coordinate of a point-axis projection
				npc = (ax*px + ay*py + az*pz)/(ax*ax + ay*ay + az*az);

				// find distance vector between point and axis
				dx = px - npc*ax;
				dy = py - npc*ay;
				dz = pz - npc*az;

				// compare position to radius
				dr = sqrt(dx*dx + dy*dy + dz*dz);
				rr = dr / r;

				// check cylinder
				if(npc >= 0.0 && npc <= 1.0 && rr <= 1)
				{
					bcvz[k][j][i] = vz * (1 - rr*rr*velType);
				}
			}
			END_STD_LOOP
		}
	}

	// restore access
	PetscCall(DMDAVecRestoreArray(fs->DA_X, bc->bcvx, &bcvx));
	PetscCall(DMDAVecRestoreArray(fs->DA_Y, bc->bcvy, &bcvy));
	PetscCall(DMDAVecRestoreArray(fs->DA_Z, bc->bcvz, &bcvz));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode BCApplyPhase(BCCtx *bc)
{
	// apply default velocity constraints on the boundaries

	FDSTAG      *fs;
	SolVarCell  *svCell;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz, iter, fixPhase;
	PetscScalar ***bcvx,  ***bcvy,  ***bcvz;

	
	PetscFunctionBeginUser;

	// access context
	fs       = bc->fs;
	fixPhase = bc->fixPhase;
	svCell   = bc->jr->svCell;

	// check constraint activation
	if(fixPhase == -1) PetscFunctionReturn(0);

	// access constraint vectors
	PetscCall(DMDAVecGetArray(fs->DA_X, bc->bcvx, &bcvx));
	PetscCall(DMDAVecGetArray(fs->DA_Y, bc->bcvy, &bcvy));
	PetscCall(DMDAVecGetArray(fs->DA_Z, bc->bcvz, &bcvz));

	// get local grid sizes
	PetscCall(DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz));

	iter = 0;

	START_STD_LOOP
	{
		// check for constrained cell
		if(svCell[iter++].phRat[fixPhase] == 1.0)
		{
			bcvx[k][j][i]   = 0.0;
			bcvx[k][j][i+1] = 0.0;

			bcvy[k][j][i]   = 0.0;
			bcvy[k][j+1][i] = 0.0;

			bcvz[k][j][i]   = 0.0;
			bcvz[k+1][j][i] = 0.0;
		}
	}
	END_STD_LOOP

	// restore access
	PetscCall(DMDAVecRestoreArray(fs->DA_X, bc->bcvx, &bcvx));
	PetscCall(DMDAVecRestoreArray(fs->DA_Y, bc->bcvy, &bcvy));
	PetscCall(DMDAVecRestoreArray(fs->DA_Z, bc->bcvz, &bcvz));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode BCApplyCells(BCCtx *bc)
{
	// apply default velocity constraints on the boundaries

	FDSTAG        *fs;
	unsigned char *fixCellFlag;
	PetscInt      i, j, k, nx, ny, nz, sx, sy, sz, iter;
	PetscScalar   ***bcvx,  ***bcvy,  ***bcvz;

	
	PetscFunctionBeginUser;

	// check activation
	if(!bc->fixCell) PetscFunctionReturn(0);

	// access context
	fs          = bc->fs;
	fixCellFlag = bc->fixCellFlag;

	// access constraint vectors
	PetscCall(DMDAVecGetArray(fs->DA_X, bc->bcvx, &bcvx));
	PetscCall(DMDAVecGetArray(fs->DA_Y, bc->bcvy, &bcvy));
	PetscCall(DMDAVecGetArray(fs->DA_Z, bc->bcvz, &bcvz));

	// get local grid sizes
	PetscCall(DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz));

	iter = 0;

	START_STD_LOOP
	{
		// check for constrained cell
		if(fixCellFlag[iter++])
		{
			bcvx[k][j][i]   = 0.0;
			bcvx[k][j][i+1] = 0.0;

			bcvy[k][j][i]   = 0.0;
			bcvy[k][j+1][i] = 0.0;

			bcvz[k][j][i]   = 0.0;
			bcvz[k+1][j][i] = 0.0;
		}
	}
	END_STD_LOOP

	// restore access
	PetscCall(DMDAVecRestoreArray(fs->DA_X, bc->bcvx, &bcvx));
	PetscCall(DMDAVecRestoreArray(fs->DA_Y, bc->bcvy, &bcvy));
	PetscCall(DMDAVecRestoreArray(fs->DA_Z, bc->bcvz, &bcvz));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode BCListSPC(BCCtx *bc)
{
	// create SPC constraint lists

	FDSTAG      *fs;
	DOFIndex    *dof;
	PetscInt    iter, numSPC, *SPCList;
	PetscInt    i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar ***bcvx,  ***bcvy,  ***bcvz, *SPCVals;

	
	PetscFunctionBeginUser;

	// access context
	fs      = bc->fs;
	dof     = &fs->dof;
	SPCVals = bc->SPCVals;
	SPCList = bc->SPCList;

	// clear constraints
	PetscCall(PetscMemzero(SPCVals, sizeof(PetscScalar)*(size_t)dof->ln));
	PetscCall(PetscMemzero(SPCList, sizeof(PetscInt)   *(size_t)dof->ln));

	// access vectors
	PetscCall(DMDAVecGetArray(fs->DA_X, bc->bcvx, &bcvx));
	PetscCall(DMDAVecGetArray(fs->DA_Y, bc->bcvy, &bcvy));
	PetscCall(DMDAVecGetArray(fs->DA_Z, bc->bcvz, &bcvz));

	iter   = 0;
	numSPC = 0;

	//---------
	// X points
	//---------

	PetscCall(DMDAGetCorners(fs->DA_X, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		LIST_SPC(bcvx, SPCList, SPCVals, numSPC, iter)

		iter++;
	}
	END_STD_LOOP

	//---------
	// Y points
	//---------

	PetscCall(DMDAGetCorners(fs->DA_Y, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		LIST_SPC(bcvy, SPCList, SPCVals, numSPC, iter)

		iter++;
	}
	END_STD_LOOP

	//---------
	// Z points
	//---------

	PetscCall(DMDAGetCorners(fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		LIST_SPC(bcvz, SPCList, SPCVals, numSPC, iter)

		iter++;
	}
	END_STD_LOOP

	// store velocity list
	bc->vNumSPC  = numSPC;
	bc->vSPCList = SPCList;
	bc->vSPCVals = SPCVals;

	// WARNING! primary pressure constraints are not implemented, otherwise compute here
	bc->pNumSPC = 0;

	// WARNING! primary temperature constraints are not implemented, otherwise compute here
	bc->tNumSPC = 0;

	// store total number of SPC constraints
	bc->numSPC = numSPC;

	// restore access
	PetscCall(DMDAVecRestoreArray(fs->DA_X, bc->bcvx, &bcvx));
	PetscCall(DMDAVecRestoreArray(fs->DA_Y, bc->bcvy, &bcvy));
	PetscCall(DMDAVecRestoreArray(fs->DA_Z, bc->bcvz, &bcvz));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// Service functions
//---------------------------------------------------------------------------
PetscErrorCode BCGetBGStrainRates(
		BCCtx       *bc,
		PetscScalar *Exx_,
		PetscScalar *Eyy_,
		PetscScalar *Ezz_,
		PetscScalar *Exy_,
		PetscScalar *Rxx_,
		PetscScalar *Ryy_,
		PetscScalar *Rzz_)
{
	// get current background strain rates & reference point coordinates

	PetscInt    jj;
	PetscScalar time, Exx, Eyy, Ezz, Exy;

	PetscFunctionBeginUser;

	// initialize
	time = bc->ts->time;
	Exx  = 0.0;
	Eyy  = 0.0;
	Ezz  = 0.0;
	Exy  = 0.0;

	// x-direction background strain rate
	if(bc->ExxNumPeriods)
	{
		for(jj = 0; jj < bc->ExxNumPeriods-1; jj++)
		{
			if(time < bc->ExxTimeDelims[jj]) break;
		}

		Exx = bc->ExxStrainRates[jj];
	}

	// y-direction background strain rate
	if(bc->EyyNumPeriods)
	{
		for(jj = 0; jj < bc->EyyNumPeriods-1; jj++)
		{
			if(time < bc->EyyTimeDelims[jj]) break;
		}

		Eyy = bc->EyyStrainRates[jj];
	}

	// z-direction background strain rate
	Ezz = -(Exx+Eyy);

	// xy-direction background strain rate
	if(bc->ExyNumPeriods)
	{
		for(jj = 0; jj < bc->ExyNumPeriods-1; jj++)
		{
			if(time < bc->ExyTimeDelims[jj]) break;
		}

		// note: we add the factor 2 here, such the the second invariant gives the specified value
		Exy = bc->ExyStrainRates[jj]*2.0;
	}

	// store result
	if(Exx_) (*Exx_) = Exx;
	if(Eyy_) (*Eyy_) = Eyy;
	if(Ezz_) (*Ezz_) = Ezz;
	if(Exy_) (*Exy_) = Exy;
	if(Rxx_) (*Rxx_) = bc->BGRefPoint[0];
	if(Ryy_) (*Ryy_) = bc->BGRefPoint[1];
	if(Rzz_) (*Rzz_) = bc->BGRefPoint[2];

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode BCGetTempBound(
		BCCtx       *bc,
		PetscScalar *Tbot)
{
	// get current bottom temperature

	PetscInt    jj;
	PetscScalar time, Tbot_val;

	PetscFunctionBeginUser;

	// initialize
	time     = bc->ts->time;
	Tbot_val = 0.0;

	//
	if(bc->TbotNumPeriods)
	{
		for(jj = 0; jj < bc->TbotNumPeriods-1; jj++)
		{
			if(time < bc->TbotTimeDelims[jj]) break;
		}

		Tbot_val = bc->Tbot[jj];
	}

	// store result
	*Tbot = Tbot_val;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode BCStretchGrid(BCCtx *bc)
{
	// apply background strain-rate "DWINDLAR" BC (Bob Shaw "Ship of Strangers")

	// Stretch grid with constant stretch factor about reference point.
	// The reference point remains fixed, and the displacements of all points are
	// proportional to the distance from the reference point.
	// Stretch factor is positive at extension, i.e.:
	// eps   = (L_new - L_old)/L_old
	// L_new = L_old + eps*L_old
	// x_new = x_old + eps*(x_old - x_ref)

	TSSol       *ts;
	FDSTAG      *fs;
	PetscScalar Exx, Eyy, Ezz;
	PetscScalar Rxx, Ryy, Rzz;

	
	PetscFunctionBeginUser;

	// access context
	fs = bc->fs;
	ts = bc->ts;

	// get background strain rates
	PetscCall(BCGetBGStrainRates(bc, &Exx, &Eyy, &Ezz, NULL, &Rxx, &Ryy, &Rzz));

	// stretch grid
	if(Exx) { PetscCall(Discret1DStretch(&fs->dsx, Exx*ts->dt, Rxx)); }
	if(Eyy) { PetscCall(Discret1DStretch(&fs->dsy, Eyy*ts->dt, Ryy)); }
	if(Ezz) { PetscCall(Discret1DStretch(&fs->dsz, Ezz*ts->dt, Rzz)); }

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode BCOverridePhase(BCCtx *bc, PetscInt cellID, Marker *P)
{
	FDSTAG     *fs;
	PetscInt    i, j, k, M, N, mx, my, sx, sy,sz,ip;
	PetscScalar z,x, y, cmax,cmin,z_plate;
	PetscScalar Temp_age,k_thermal,dT_adiabatic,Z_Top,Tbot;
	PetscInt phase_inflow;
	PetscScalar T_inflow;
	
	PetscFunctionBeginUser;

	// get time-dependent Tbot
	PetscCall(BCGetTempBound(bc, &Tbot));

	if((bc->face) || bc->Plume_Inflow || bc->bot_open)
	{
		fs = bc->fs;
		M  = fs->dsx.ncels;
		N  = fs->dsy.ncels;
		sx = fs->dsx.pstart;
		sy = fs->dsy.pstart;
		sz = fs->dsz.pstart;
		mx = fs->dsx.tcels-1;
		my = fs->dsy.tcels-1;
		z  = P->X[2];
		x = P->X[0];
		y = P->X[1];

		GET_CELL_IJK(cellID, i, j, k, M, N);

		if(((bc->face == 1 && i + sx == 0)
				||  (bc->face == 2 && i + sx == mx)
				||  (bc->face == 3 && j + sy == 0)
				||  (bc->face == 4 && j + sy == my))
				&&  (z >= bc->bot && z <= bc->top) && (bc->bvel_temperature_inflow>0))
		{
			if(bc->jr->ctrl.Adiabatic_gr > 0.0)
			{
				if(bc->jr->surf->UseFreeSurf)
				{
					Z_Top = bc->jr->surf->InitLevel;
				}
				else
				{
					Z_Top = bc->fs->dsz.gcrdend;
				}

				dT_adiabatic= bc->jr->ctrl.Adiabatic_gr*PetscAbs(z-Z_Top);
			}
			else
			{
				dT_adiabatic = 0.0;
			}

			if(bc->bvel_temperature_inflow==2)
			{
				k_thermal= 1e-6/( (bc->scal->length_si)*(bc->scal->length_si)/(bc->scal->time_si));
				z_plate = PetscAbs(z-bc->top);
				Temp_age = (bc->bvel_potential_temperature-bc->bvel_temperature_top)*erf(z_plate/2.0/sqrt(k_thermal*bc->bvel_thermal_age)) + bc->bvel_temperature_top;
				P->T = Temp_age+dT_adiabatic;
			}
			else if(bc->bvel_temperature_inflow == 1)
			{
				P->T=bc->bvel_constant_temperature;
			}
		}

		if(bc->num_phase_bc >= 0)
		{
			// expand i, j, k cell indices
			if(((bc->face == 1 && i + sx == 0)
					||  (bc->face == 2 && i + sx == mx)
					||  (bc->face == 3 && j + sy == 0)
					||  (bc->face == 4 && j + sy == my))
					&&  (z >= bc->bot-bc->relax_dist && z <= bc->top+bc->relax_dist))
			{
				for(ip=0;ip<bc->num_phase_bc;ip++)
				{
					if(z>=bc->phase_interval[ip] && z<bc->phase_interval[ip+1])
					{
						P->phase = bc->phase[ip];
					}
				}
			}
		}

		// if we have have a inflow condition @ the lower boundary, we change the phase of the particles within the zone
		if(k+sz == 0)
		{
			if(bc->Plume_Inflow == 1)
			{
				// This routine handle the inflow and outflow. If the Plume boundary is "permeable" type, within the plume radius the phase that are
				// injected is the one prescribed for the plume. The particle injected has the same temperature of the TBot (i.e. according to a gaussian thermal
				// perturbation). Otherwise has the phase and temperature of the background mantle.

				phase_inflow = bc->phase_inflow_bot;

				if(bc->Plume_Dimension==1)
				{
					T_inflow = Tbot + (bc->Plume_Temperature-Tbot)*PetscExpScalar( - PetscPowScalar(x-bc->Plume_Center[0],2.0 ) /(PetscPowScalar(bc->Plume_Radius,2.0))) ;

					cmin = bc->Plume_Center[0] - bc->Plume_Radius;
					cmax = bc->Plume_Center[0] + bc->Plume_Radius;

					if(x>=cmin && x<=cmax)
					{
						phase_inflow = bc->Plume_Phase;
					}
				}
				else
				{
					T_inflow = Tbot + (bc->Plume_Temperature-Tbot)*PetscExpScalar( - ( PetscPowScalar(x-bc->Plume_Center[0],2.0 ) + PetscPowScalar(y-bc->Plume_Center[1],2.0 ) )/(PetscPowScalar(bc->Plume_Radius,2.0)));

					if (PetscPowScalar((x - bc->Plume_Center[0]),2.0) +
						PetscPowScalar((y - bc->Plume_Center[1]),2.0) <= PetscPowScalar( bc->Plume_Radius,2.0) )
					{
						phase_inflow = bc->Plume_Phase;
					}
				}

				P->phase  = phase_inflow;
				P->T      = T_inflow;
				//	PetscPrintf(PETSC_COMM_WORLD,"Plume Temperature P->T=%6f \n",P->T*bc->scal->temperature-bc->scal->Tshift);

			}
			else if(bc->bot_open)
			{
				P->phase = bc->phase_inflow_bot;
				P->T     = Tbot;
			}
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode BC_Plume_inflow(BCCtx *bc)
{
	FDSTAG          *fs;
	PetscInt        i, j, k, nx, ny, nz, sx, sy, sz;
	PetscScalar     ***bcvz;
	PetscScalar     vel, x_min,x_max,y_min,y_max,x,y;
	PetscScalar     Area_Bottom, Area_Inflow, Area_Outflow, V_avg, V_in, V_out, Qin, areaFrac;
	PetscScalar     radius2, R;

	
	PetscFunctionBeginUser;

	if(!bc->Plume_Inflow) PetscFunctionReturn(0);

	fs = bc->fs;

	PetscCall(FDSTAGGetGlobalBox(bc->fs, &x_min, &y_min,NULL, &x_max, &y_max, NULL));

	V_in      = bc->Plume_Inflow_Velocity; // max. inflow velocity
	areaFrac  = bc->Plume_areaFrac;

	if(bc->Plume_Dimension == 1)
	{
		// 2D
		Area_Bottom  = (x_max-x_min);
		Area_Inflow  =  2.0*bc->Plume_Radius;    // inflow length
		Area_Outflow = Area_Bottom-Area_Inflow;  // outflow length
	}
	else
	{
		Area_Bottom  = (x_max-x_min)*(y_max-y_min);
		Area_Inflow  = PETSC_PI*bc->Plume_Radius*bc->Plume_Radius; // inflow
		Area_Outflow = Area_Bottom-Area_Inflow;
	}

	if(bc->Plume_VelocityType==0 )
	{
		// Poiseuille-type inflow condition
		//  Note that this results in a velocity discontinuity at the border

		// We assume Poiseuille flow between plates (2D) or in a pipe (3D):
		if(bc->Plume_Dimension == 1) { V_avg = V_in*2.0/3.0; } // 2D
		else                         { V_avg = V_in*1.0/2.0; } // 3D

		// outflow velocity is based on mass conservation (i.e.: Qin+Qout=0)
		Qin   =  V_avg * Area_Inflow * areaFrac; // volume influx
		V_out = -Qin / Area_Outflow;             // outflow velocity
	}
	else
	{
		// Gaussian-like inflow perturbation
		if(bc->Plume_Dimension == 1)
		{
			// 2D
			PetscScalar a,b,c, xc;

			// Gaussian perturbation velocity - anything that creates a rigid plug is a problem
			// we integrate the velocity profile over the full domain as:
			// V = V_out + (V_in-V_out)*exp(-((x-xc)^2)/c^2) from x=xmin..xmax
			// We can do this with sympy, which gives:
			// V_avg = V_out + (V_in-V_out)*(sqrt(pi)*c*erf((-xc + x_max)/c)/2 - sqrt(pi)*c*erf((-xc + x_min)/c)/2))/(x_max-x_min)
			// V_avg = V_out + (V_in-V_out)*(a-b)  ->  V_out*(1-(a-b)) = -V_in*(a-b), so V_out =  -V_in*(a-b)/(1-(a-b))

			xc      =   bc->Plume_Center[0];
			c       =   bc->Plume_Radius;
			a       =   PetscSqrtScalar(PETSC_PI)*c*erf((-xc + x_max)/c)/2.0/(x_max-x_min); // dV
			b       =   PetscSqrtScalar(PETSC_PI)*c*erf((-xc + x_min)/c)/2.0/(x_max-x_min); // dV

			// average velocity should be zero
			V_out   =   -V_in*(a-b)/(1-(a-b))*areaFrac;

		}
		else
		{
			// 3D
			// In 3D, the expression for the velocity is:
			// V = V_out + (V_in-V_out)*exp(-((x-xc)^2 + (y-yc)^2)/c^2) from  x = xmin..xmax and y=y_min .. y_max

			PetscScalar a, b, d, e, c, xc, yc;

			xc      =   bc->Plume_Center[0];
			yc      =   bc->Plume_Center[1];
			c       =   bc->Plume_Radius;

			a       =   1.0/4.0 * PETSC_PI*erf((-xc + x_max)/c)*erf((-yc + y_max)/c)/Area_Bottom;
			b       =   1.0/4.0 * PETSC_PI*erf((-xc + x_min)/c)*erf((-yc + y_max)/c)/Area_Bottom;

			d       =   1.0/4.0 * PETSC_PI*erf((-xc + x_min)/c)*erf((-yc + y_min)/c)/Area_Bottom;
			e       =   1.0/4.0 * PETSC_PI*erf((-xc + x_max)/c)*erf((-yc + y_min)/c)/Area_Bottom;

			//      so V_avg = V_out + (V_in-V_out)*((a-b)/Area + (d-e)/Area)
			// since we want V_avg = 0, we can compute V_out as:

			V_out   =   -V_in*(a-b + d-e)/(1-(a-b + d-e))*areaFrac;                     // average velocity should be zero
		}
	}

	// access constraint vectors
	PetscCall(DMDAVecGetArray(fs->DA_Z,   bc->bcvz, &bcvz));

	//=========================================================================
	// SPC (normal velocities)
	//=========================================================================

	PetscCall(DMDAGetCorners(fs->DA_Z, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		x       = COORD_CELL(i, sx, fs->dsx);
		radius2 = PetscPowScalar(bc->Plume_Radius,2.0);

		if(bc->Plume_VelocityType==0 )
		{
			// Poiseuille type inflow
			if(bc->Plume_Dimension == 1)
			{
				// 2D
				R       =   PetscPowScalar((x-bc->Plume_Center[0]),2.0);
			}
			else
			{
				// 3D
				y       =   COORD_CELL(j, sy, fs->dsy);
				R       =   PetscPowScalar((x-bc->Plume_Center[0]),2.0) + PetscPowScalar((y-bc->Plume_Center[1]),2.0);
			}
			if  ( R <=  radius2 ) {  vel = V_in*(1.0 - R/radius2);   }
			else                  {  vel = V_out;                    }
		}

		else
		{
			// Gaussian plume inflow condition

			PetscScalar xc;

			xc = bc->Plume_Center[0];

			// gaussian type perturbation
			if(bc->Plume_Dimension == 1)
			{
				// 2D
				// Gaussian velocity perturbation
				vel = V_out + (V_in-V_out)*PetscExpScalar( - PetscPowScalar(x-xc,2.0 ) /radius2 ) ;

			}
			else
			{
				PetscScalar yc;
				yc  =   bc->Plume_Center[1];
				y   =   COORD_CELL(j, sy, fs->dsy);

				// Gaussian velocity perturbation
				vel = V_out+(V_in-V_out)*PetscExpScalar( - ( PetscPowScalar(x-xc,2.0 ) + PetscPowScalar(y-yc,2.0 ) )/radius2 );

			}
		}

		if(k == 0)
		{
			bcvz[k][j][i] = vel;
		}
	}
	END_STD_LOOP

	// restore access
	PetscCall(DMDAVecRestoreArray(fs->DA_Z,   bc->bcvz, &bcvz));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

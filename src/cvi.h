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
 **    filename:   cvi.h
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
//..........   Conservative Velocity Interpolation Routines (CVI)  .........
//---------------------------------------------------------------------------
#ifndef __cvi_h__
#define __cvi_h__
//-----------------------------------------------------------------------------
// structures
//-----------------------------------------------------------------------------
typedef struct
{
	PetscScalar      x0[3];    // initial position
	PetscScalar      x[3];     // position to interpolate
	PetscScalar      v[3];     // velocity interpolated
	PetscScalar      v_eff[3]; // effective velocity
	PetscInt         ind;      // global index of markers

} VelInterp;
//-----------------------------------------------------------------------------
typedef enum
{
	EULER,          // euler explicit in time
	RUNGE_KUTTA_2,  // runge-kutta 2nd order in space
	RUNGE_KUTTA_4   // runge-kutta 4th order in space

} AdvectionType;
//-----------------------------------------------------------------------------
typedef enum
{
	STAG,      // trilinear interp from fdstag points
	Q2,        // quadratic interpolation
	SQ2,       // spline quadratic interpolation
	NODES,     // bilinear interp to nodes, trilinear interp to markers + correction
	MINMOD,    // minmod interp to nodes, trilinear interp to markers + correction
	CorrQ2,    // quadratic interpolation + correction
	CorrSQ2,   // spline quadratic interpolation + correction
	STAG_P     // empirical approach (T. Gerya)

} VelInterpType;
//-----------------------------------------------------------------------------
typedef struct
{
	AdvectionType    advection;    // advection scheme
	VelInterpType    velinterp;    // velocity interpolation
	VelInterp        *interp;
	PetscInt         nmark;        // number of markers to interpolate
	PetscInt         nbuff;        // buffer size

	// FDSTAG context
	FDSTAG           *fs;
	JacRes           *jr;

	// point-cell interaction
	PetscInt         *cellnum;
	PetscInt         *markind;
	PetscInt         *markstart;

	// exchange
	MPI_Comm         icomm;   // distinct communicator for communicating markers
	PetscInt         nproc;   // total number of processors
	PetscInt         iproc;   // processor rank

	VelInterp        *sendbuf;
	VelInterp        *recvbuf;

	PetscInt         nsend;                // total number of markers to be sent (local)
	PetscInt         nsendm[_num_neighb_]; // number of markers to be sent to each process
	PetscInt         ptsend[_num_neighb_]; // send buffer pointers

	PetscInt         nrecv;
	PetscInt         nrecvm[_num_neighb_]; // number of markers to be received from each process
	PetscInt         ptrecv[_num_neighb_]; // receive buffer pointers

	PetscInt         ndel;
	PetscInt         *idel;

} AdvVelCtx;

//-----------------------------------------------------------------------------
// main routines
//-----------------------------------------------------------------------------

// major advection routines
PetscErrorCode ADVelAdvectMain     (AdvCtx *actx);
PetscErrorCode ADVelReadOptions    (AdvVelCtx *vi);
PetscErrorCode ADVelInterpPT       (AdvCtx *actx);
PetscErrorCode ADVelAdvectScheme   (AdvCtx *actx, AdvVelCtx *vi);
PetscErrorCode ADVelCollectIndices (AdvCtx *actx, AdvVelCtx *vi);

// subroutines
PetscErrorCode ADVelCreate         (AdvCtx *actx, AdvVelCtx *vi);
PetscErrorCode ADVelDestroy        (AdvVelCtx *vi);

// Runge-Kutta step
PetscErrorCode ADVelRungeKuttaStep (AdvVelCtx *vi, PetscScalar dt, PetscScalar a, PetscInt type);

// coordinate manipulation
PetscErrorCode ADVelInitCoord      (AdvCtx *actx, VelInterp *interp, PetscInt n);
PetscErrorCode ADVelRetrieveCoord  (AdvCtx *actx, VelInterp *interp, PetscInt n);
PetscErrorCode ADVelAdvectCoord    (VelInterp *interp, PetscInt n, PetscScalar dt, PetscInt type);
PetscErrorCode ADVelResetCoord     (VelInterp *interp, PetscInt n);

// effective velocity
PetscErrorCode ADVelCalcEffVel     (VelInterp *interp, PetscInt n, PetscScalar a);

// used for exchange
PetscErrorCode ADVelExchange       (AdvVelCtx *vi);
PetscErrorCode ADVelDeleteOutflow  (AdvVelCtx *vi);
PetscErrorCode ADVelMapToDomains   (AdvVelCtx *vi);
PetscErrorCode ADVelExchangeNMark  (AdvVelCtx *vi);
PetscErrorCode ADVelCreateMPIBuff  (AdvVelCtx *vi);
PetscErrorCode ADVelExchangeMark   (AdvVelCtx *vi);
PetscErrorCode ADVelDestroyMPIBuff (AdvVelCtx *vi);
PetscErrorCode ADVelCollectGarbage (AdvVelCtx *vi);
PetscErrorCode ADVelReAllocStorage (AdvVelCtx *vi, PetscInt nmark);
PetscErrorCode ADVelMapMarkToCells (AdvVelCtx *vi);

// velocity interpolation
PetscErrorCode ADVelInterpMain       (AdvVelCtx *vi);
PetscErrorCode ADVelInterpSTAG       (AdvVelCtx *vi);
PetscErrorCode ADVelInterpNODES      (AdvVelCtx *vi);
PetscErrorCode ADVelInterpMINMOD     (AdvVelCtx *vi);
PetscErrorCode ADVelInterpQ2         (AdvVelCtx *vi);
PetscErrorCode ADVelInterpSQ2        (AdvVelCtx *vi);
PetscErrorCode ADVelInterpCorrQ2     (AdvVelCtx *vi);
PetscErrorCode ADVelInterpCorrSQ2    (AdvVelCtx *vi);
PetscErrorCode ADVelInterpSTAGP      (AdvVelCtx *vi);

//-----------------------------------------------------------------------------
// service functions
//-----------------------------------------------------------------------------

static inline PetscScalar GenInterpLin3D(
	PetscScalar A[8],
	PetscScalar xe,
	PetscScalar ye,
	PetscScalar ze)
{
	PetscScalar v, xb, yb, zb;

	xb = 1.0 - xe;
	yb = 1.0 - ye;
	zb = 1.0 - ze;

	// interpolate & return result
	v =
	A[0]*xb*yb*zb +
	A[1]*xe*yb*zb +
	A[2]*xb*ye*zb +
	A[3]*xe*ye*zb +
	A[4]*xb*yb*ze +
	A[5]*xe*yb*ze +
	A[6]*xb*ye*ze +
	A[7]*xe*ye*ze;

	return v;
}
//-----------------------------------------------------------------------------
static inline PetscScalar GenInterpLin2D(
	PetscScalar A[4],
	PetscScalar xe,
	PetscScalar ye)
{
	PetscScalar v, xb, yb;

	xb = 1.0 - xe;
	yb = 1.0 - ye;

	// interpolate & return result
	v =
	A[0]*xb*yb +
	A[1]*xe*yb +
	A[2]*xb*ye +
	A[4]*xe*ye;

	return v;
}
//-----------------------------------------------------------------------------
static inline PetscScalar minmod(
	PetscScalar a,
	PetscScalar b)
{
	// minmod slope limiter
	PetscScalar X,c;

	X = 0.0;
	c = a*b;

	if ((c > 0) && fabs(a)<=fabs(b)) X = a;
	if ((c > 0) && fabs(a)> fabs(b)) X = b;

	return X;
}
//-----------------------------------------------------------------------------
static inline PetscScalar InterpLinMinmod(
	PetscScalar U1,
	PetscScalar U2,
	PetscScalar U3,
	PetscScalar dx,
	PetscScalar x1,
	PetscScalar x2,
	PetscScalar x3,
	PetscInt    d)
{
	PetscScalar X;

	X = 0.0;

	if (d==0) X = U2 - dx*0.5*minmod((U3-U2)/(x3-x2),(U2-U1)/(x2-x1));
	else      X = U2 + dx*0.5*minmod((U3-U2)/(x3-x2),(U2-U1)/(x2-x1));

	return X;
}
//-----------------------------------------------------------------------------
static inline PetscScalar QuadraticCoeff(
	PetscScalar v0, PetscScalar v1, PetscScalar v2,
	PetscScalar x0, PetscScalar x1, PetscScalar x2,
	PetscScalar x,
	PetscInt type)
{
	PetscScalar v;

	// type 0 - quadratic interpolation
	if (type == 0)
	{
		PetscScalar a, b, c, dx0, dx1, dx2;

		dx0 = x0-x1;
		dx1 = x0-x2;
		dx2 = x1-x2;

		// coefficients
		a =
				v0/  dx0 /  dx1 +
				v1/(-dx0)/  dx2 +
				v2/(-dx1)/(-dx2);

		b = (v0-v1)/dx0 - a*(x0+x1);
		c = v0 - a*x0*x0 - b*x0;

		// interpolate & return result
		v = a*x*x+b*x+c;
	}

	// type 1 - spline quadratic interpolation
	else if (type == 1)
	{
		PetscScalar a1, b0, b1, c0, c1;

		// coefficients
		b0 = (v1-v0)/(x1-x0);
		c0 = v0-b0*x0;

		a1 = ((v1-v2)/(x1-x2)-b0)/(x2-x1);
		b1 = b0-2*a1*x1;
		c1 = v1-b1*x1-a1*x1*x1;

		// interpolate & return result
		if (x < x1) v = b0*x+c0;
		else        v = a1*x*x+b1*x+c1;
	}

	return v;
}
//-----------------------------------------------------------------------------
static inline PetscScalar GenInterpQuadratic2D(
	PetscScalar vi[9],
	PetscScalar xi1[3],
	PetscScalar xi2[3],
	PetscScalar x1,
	PetscScalar x2,
	PetscInt type)
{
	PetscScalar v, v1, v2, v3;

	// intermediate velocities
	v1 = QuadraticCoeff(vi[0], vi[3], vi[6], xi2[0], xi2[1], xi2[2], x2, type);
	v2 = QuadraticCoeff(vi[1], vi[4], vi[7], xi2[0], xi2[1], xi2[2], x2, type);
	v3 = QuadraticCoeff(vi[2], vi[5], vi[8], xi2[0], xi2[1], xi2[2], x2, type);

	// new velocity
	v = QuadraticCoeff(v1, v2, v3, xi1[0], xi1[1], xi1[2], x1, type);

	return v;
}
//-----------------------------------------------------------------------------
static inline PetscScalar GenInterpQuadratic3D(
	PetscScalar vi[27],
	PetscScalar xi1[3],
	PetscScalar xi2[3],
	PetscScalar xi3[3],
	PetscScalar x1,
	PetscScalar x2,
	PetscScalar x3,
	PetscInt type)
{
	PetscScalar v, v1[9], v2[3];

	// intermediate velocities 1
	v1[0] = QuadraticCoeff(vi[0], vi[ 9], vi[18], xi3[0], xi3[1], xi3[2], x3, type);
	v1[1] = QuadraticCoeff(vi[1], vi[10], vi[19], xi3[0], xi3[1], xi3[2], x3, type);
	v1[2] = QuadraticCoeff(vi[2], vi[11], vi[20], xi3[0], xi3[1], xi3[2], x3, type);
	v1[3] = QuadraticCoeff(vi[3], vi[12], vi[21], xi3[0], xi3[1], xi3[2], x3, type);
	v1[4] = QuadraticCoeff(vi[4], vi[13], vi[22], xi3[0], xi3[1], xi3[2], x3, type);
	v1[5] = QuadraticCoeff(vi[5], vi[14], vi[23], xi3[0], xi3[1], xi3[2], x3, type);
	v1[6] = QuadraticCoeff(vi[6], vi[15], vi[24], xi3[0], xi3[1], xi3[2], x3, type);
	v1[7] = QuadraticCoeff(vi[7], vi[16], vi[25], xi3[0], xi3[1], xi3[2], x3, type);
	v1[8] = QuadraticCoeff(vi[8], vi[17], vi[26], xi3[0], xi3[1], xi3[2], x3, type);

	// intermediate velocities 2
	v2[0] = QuadraticCoeff(v1[0], v1[3], v1[6], xi2[0], xi2[1], xi2[2], x2, type);
	v2[1] = QuadraticCoeff(v1[1], v1[4], v1[7], xi2[0], xi2[1], xi2[2], x2, type);
	v2[2] = QuadraticCoeff(v1[2], v1[5], v1[8], xi2[0], xi2[1], xi2[2], x2, type);

	// new velocity
	v = QuadraticCoeff(v2[0], v2[1], v2[2], xi1[0], xi1[1], xi1[2], x1, type);

	return v;
}

#endif

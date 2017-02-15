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

struct VelInterp
{
	PetscScalar      x0[3];    // initial position
	PetscScalar      x[3];     // position to interpolate
	PetscScalar      v[3];     // velocity interpolated
	PetscScalar      v_eff[3]; // effective velocity
	PetscInt         ind;      // global index of markers

};

//-----------------------------------------------------------------------------
struct AdvVelCtx
{

	VelInterp        *interp;
	PetscInt         nmark;        // number of markers to interpolate
	PetscInt         nbuff;        // buffer size

	// FDSTAG context
	FDSTAG           *fs;
	JacRes           *jr;
	AdvCtx           *actx;

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

};

//-----------------------------------------------------------------------------
// main routines
//-----------------------------------------------------------------------------

// major advection routines
PetscErrorCode ADVelAdvectMain     (AdvCtx *actx);
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
PetscErrorCode ADVelInterpMINMOD     (AdvVelCtx *vi);
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
#endif

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
//..... LOCAL LEVEL INTEGRATION ALGORITHMS FOR CONSTITUTIVE EQUATIONS  ......
//---------------------------------------------------------------------------
#if defined(_WIN32) && !defined(_USE_MATH_DEFINES)
// MinGW's strict-ISO mode hides M_PI (a GNU/POSIX extension) from <math.h>/
// <cmath> unless this is defined before the first include of either.
#define _USE_MATH_DEFINES
#endif
#include "LaMEM.h"
#include "Tensor.h"
#include "tools.h"

//---------------------------------------------------------------------------
// Elastic stress rotation functions
//---------------------------------------------------------------------------
void GetRotationMatrix(
    Tensor2RN   *R,  // rotation matrix
    PetscScalar  dt, // time step
    PetscScalar  wx, // vorticity vector components
    PetscScalar  wy, // ...
    PetscScalar  wz) // ...
{
	// compute rotation matrix from axis & angle (Euler-Rodrigues formula)
	// WARNING! Courant criterion for rotation angle should be implemented
	// angle tolerance probably should be also larger than machine epsilon

	PetscScalar w, theta, ct, st, cf;

	// get length of the instantaneous rotation axis (vorticity intensity)
	w = sqrt(wx*wx + wy*wy + wz*wz);

	// get finite rotation angle (vorticity intensity is twice the average angular velocity)
	theta = dt*(w/2.0);

	// round-off
	if(theta <  2.0*DBL_EPSILON)
	{
		Tensor2RNUnit(R);

		return;
	}

	// get unit rotation axis
	wx /= w;
	wy /= w;
	wz /= w;

	// compute rotation operator using Euler-Rodrigues formula
	ct = cos(theta);
	st = sin(theta);
	cf = 1.0 - ct;

	R->xx =  ct    + cf*wx*wx;   R->xy = -st*wz + cf*wx*wy;   R->xz =  st*wy + cf*wx*wz;
	R->yx =  st*wz + cf*wy*wx;   R->yy =  ct    + cf*wy*wy;   R->yz = -st*wx + cf*wy*wz;
	R->zx = -st*wy + cf*wz*wx;   R->zy =  st*wx + cf*wz*wy;   R->zz =  ct    + cf*wz*wz;
}
//---------------------------------------------------------------------------
void RotateStress(Tensor2RN *R, Tensor2RS *S, Tensor2RS *SR)
{
	// rotate stress tensor
	// [SR] = [R] * [S] * [R]^t
	// SRij = Rik * Skl * Rjl

	PetscScalar dx, dy, dz;

	dx = R->xx*S->xx + R->xy*S->xy + R->xz*S->xz;
	dy = R->xx*S->xy + R->xy*S->yy + R->xz*S->yz;
	dz = R->xx*S->xz + R->xy*S->yz + R->xz*S->zz;

	SR->xx = dx*R->xx + dy*R->xy + dz*R->xz;
	SR->xy = dx*R->yx + dy*R->yy + dz*R->yz;
	SR->xz = dx*R->zx + dy*R->zy + dz*R->zz;

	dx = R->yx*S->xx + R->yy*S->xy + R->yz*S->xz;
	dy = R->yx*S->xy + R->yy*S->yy + R->yz*S->yz;
	dz = R->yx*S->xz + R->yy*S->yz + R->yz*S->zz;

	SR->yy = dx*R->yx + dy*R->yy + dz*R->yz;
	SR->yz = dx*R->zx + dy*R->zy + dz*R->zz;

	dx = R->zx*S->xx + R->zy*S->xy + R->zz*S->xz;
	dy = R->zx*S->xy + R->zy*S->yy + R->zz*S->yz;
	dz = R->zx*S->xz + R->zy*S->yz + R->zz*S->zz;

	SR->zz = dx*R->zx + dy*R->zy + dz*R->zz;
}
//---------------------------------------------------------------------------
void Tensor2RSCopy(Tensor2RS *A, Tensor2RS *B)
{
	// copy symmetric second order tensor B = A
	B->xx = A->xx;
	B->xy = A->xy; B->yy = A->yy;
	B->xz = A->xz; B->yz = A->yz; B->zz = A->zz;
}
//---------------------------------------------------------------------------
void Tensor2RNClear(Tensor2RN *A)
{
	A->xx = 0.0; A->xy = 0.0; A->xz = 0.0;
	A->yx = 0.0; A->yy = 0.0; A->yz = 0.0;
	A->zx = 0.0; A->zy = 0.0; A->zz = 0.0;
}
//---------------------------------------------------------------------------
void Tensor2RSClear(Tensor2RS *A)
{
	A->xx = 0.0;
	A->xy = 0.0; A->yy = 0.0;
	A->xz = 0.0; A->yz = 0.0; A->zz = 0.0;
}
//---------------------------------------------------------------------------
PetscInt Tensor2RNCheckEq(Tensor2RN *A, Tensor2RN *B, PetscScalar tol)
{
	if(!CHECKEQ(A->xx, B->xx, tol, DBL_EPSILON)) return 0;
	if(!CHECKEQ(A->xy, B->xy, tol, DBL_EPSILON)) return 0;
	if(!CHECKEQ(A->xz, B->xz, tol, DBL_EPSILON)) return 0;
	if(!CHECKEQ(A->yx, B->yx, tol, DBL_EPSILON)) return 0;
	if(!CHECKEQ(A->yy, B->yy, tol, DBL_EPSILON)) return 0;
	if(!CHECKEQ(A->yz, B->yz, tol, DBL_EPSILON)) return 0;
	if(!CHECKEQ(A->zx, B->zx, tol, DBL_EPSILON)) return 0;
	if(!CHECKEQ(A->zy, B->zy, tol, DBL_EPSILON)) return 0;
	if(!CHECKEQ(A->zz, B->zz, tol, DBL_EPSILON)) return 0;

	return 1;
}
//---------------------------------------------------------------------------
void Tensor2RNNorm(Tensor2RN *A, PetscScalar *pk)
{
	// k = |A|

	PetscScalar s, k;

	s = fabs(A->xx) + fabs(A->xy) + fabs(A->xz);           k = s;
	s = fabs(A->yx) + fabs(A->yy) + fabs(A->yz); if(s > k) k = s;
	s = fabs(A->zx) + fabs(A->zy) + fabs(A->zz); if(s > k) k = s;

	(*pk) = k;
}
//---------------------------------------------------------------------------
void Tensor2RSNorm(Tensor2RS *A, PetscScalar *pk)
{
	// k = |A|

	PetscScalar s, k;

	s = fabs(A->xx) + fabs(A->xy) + fabs(A->xz);           k = s;
	s = fabs(A->xy) + fabs(A->yy) + fabs(A->yz); if(s > k) k = s;
	s = fabs(A->xz) + fabs(A->yz) + fabs(A->zz); if(s > k) k = s;

	(*pk) = k;
}
//---------------------------------------------------------------------------
void Tensor2RNDivide(Tensor2RN *A, PetscScalar k)
{
	// A = A/k

	A->xx /= k; A->xy /= k; A->xz /= k;
	A->yx /= k; A->yy /= k; A->yz /= k;
	A->zx /= k; A->zy /= k; A->zz /= k;
}
//---------------------------------------------------------------------------
void Tensor2RNTrace(Tensor2RN *A)
{
	// A = A - tr[A]/3*I

	PetscScalar tr;

	tr = (A->xx + A->yy + A->zz)/3.0;

	A->xx -= tr;
	A->yy -= tr;
	A->zz -= tr;
}
//---------------------------------------------------------------------------
void Tensor2RNSym(Tensor2RN *A, Tensor2RN *B)
{
	// B = (A + A^t)/2

	B->xx =  A->xx;                B->xy = (A->yx + A->xy)/2.0;   B->xz = (A->zx + A->xz)/2.0;
	B->yx = (A->xy + A->yx)/2.0;   B->yy =  A->yy;                B->yz = (A->zy + A->yz)/2.0;
	B->zx = (A->xz + A->zx)/2.0;   B->zy = (A->yz + A->zy)/2.0;   B->zz =  A->zz;
}
//---------------------------------------------------------------------------
void Tensor2RNProduct(Tensor2RN *A, Tensor2RN *B, Tensor2RN *C)
{
	// C = A*B

	C->xx = A->xx*B->xx + A->xy*B->yx + A->xz*B->zx;
	C->xy = A->xx*B->xy + A->xy*B->yy + A->xz*B->zy;
	C->xz = A->xx*B->xz + A->xy*B->yz + A->xz*B->zz;
	C->yx = A->yx*B->xx + A->yy*B->yx + A->yz*B->zx;
	C->yy = A->yx*B->xy + A->yy*B->yy + A->yz*B->zy;
	C->yz = A->yx*B->xz + A->yy*B->yz + A->yz*B->zz;
	C->zx = A->zx*B->xx + A->zy*B->yx + A->zz*B->zx;
	C->zy = A->zx*B->xy + A->zy*B->yy + A->zz*B->zy;
	C->zz = A->zx*B->xz + A->zy*B->yz + A->zz*B->zz;
}
//---------------------------------------------------------------------------
void Tensor2RNTranspose(Tensor2RN *A, Tensor2RN *B)
{
	// B = A^t

	B->xx = A->xx; B->xy = A->yx; B->xz = A->zx;
	B->yx = A->xy; B->yy = A->yy; B->yz = A->zy;
	B->zx = A->xz; B->zy = A->yz; B->zz = A->zz;
}
//---------------------------------------------------------------------------
void Tensor2RNCopy(Tensor2RN *A, Tensor2RN *B)
{
	// B = A

	B->xx = A->xx; B->xy = A->xy; B->xz = A->xz;
	B->yx = A->yx; B->yy = A->yy; B->yz = A->yz;
	B->zx = A->zx; B->zy = A->zy; B->zz = A->zz;
}
//---------------------------------------------------------------------------
void Tensor2RNCopySym(Tensor2RN *A, Tensor2RS *B)
{
	// B = A

	B->xx = A->xx;
	B->xy = A->xy; B->yy = A->yy;
	B->xz = A->xz; B->yz = A->yz; B->zz = A->zz;
}

//---------------------------------------------------------------------------
void Tensor2RNUnit(Tensor2RN *A)
{
	// A = I

	A->xx = 1.0; A->xy = 0.0; A->xz = 0.0;
	A->yx = 0.0; A->yy = 1.0; A->yz = 0.0;
	A->zx = 0.0; A->zy = 0.0; A->zz = 1.0;
}
//---------------------------------------------------------------------------
void Tensor2RNSum3(
    Tensor2RN *A, PetscScalar ka,
    Tensor2RN *B, PetscScalar kb,
    Tensor2RN *C, PetscScalar kc,
    Tensor2RN *R)
{
	// R = ka*A + kb*B + kc*C

	R->xx = ka*A->xx + kb*B->xx + kc*C->xx;
	R->xy = ka*A->xy + kb*B->xy + kc*C->xy;
	R->xz = ka*A->xz + kb*B->xz + kc*C->xz;
	R->yx = ka*A->yx + kb*B->yx + kc*C->yx;
	R->yy = ka*A->yy + kb*B->yy + kc*C->yy;
	R->yz = ka*A->yz + kb*B->yz + kc*C->yz;
	R->zx = ka*A->zx + kb*B->zx + kc*C->zx;
	R->zy = ka*A->zy + kb*B->zy + kc*C->zy;
	R->zz = ka*A->zz + kb*B->zz + kc*C->zz;
}
//---------------------------------------------------------------------------
void Tensor2RNView(Tensor2RN *A, const char *msg)
{
	printf("%s: \n\n", msg);
	printf("%g %g %g \n",   A->xx, A->xy, A->xz);
	printf("%g %g %g \n",   A->yx, A->yy, A->yz);
	printf("%g %g %g \n\n", A->zx, A->zy, A->zz);
}
//---------------------------------------------------------------------------
void Tensor2RSView(Tensor2RS *A, const char *msg)
{
	printf("%s: \n\n", msg);
	printf("%g %g %g \n",   A->xx, A->xy, A->xz);
	printf("%g %g %g \n",   A->xy, A->yy, A->yz);
	printf("%g %g %g \n\n", A->xz, A->yz, A->zz);
}
//---------------------------------------------------------------------------
PetscErrorCode Tensor2RS2DSpectral(
    PetscScalar  axx,
    PetscScalar  ayy,
    PetscScalar  axy,
    PetscScalar *pa1,
    PetscScalar *pa2,
    PetscScalar  v1[],
    PetscScalar  v2[],
    PetscScalar  tol)
{
	PetscScalar theta, t, c, s, tau, nrm, sum, a1, a2, a, v[2];

	PetscFunctionBeginUser;

	// get stress norm
	sum = fabs(axx) + fabs(axy);               nrm = sum;
	sum = fabs(axy) + fabs(ayy); if(sum > nrm) nrm = sum;

	// initialize eigenvalues & eigenvectors
	a1 = axx;
	a2 = ayy;

	v1[0] = 1.0; v2[0] = 0.0;
	v1[1] = 0.0; v2[1] = 1.0;

	// compute eigenvectors & eigenvalues
	if(fabs(axy) > tol*nrm)
	{
		theta = 0.5*(ayy - axx)/axy;
		t     = 1.0/(fabs(theta) + sqrt(theta*theta + 1.0));
		if(theta < 0.0) t = -t;

		a1 -= t*axy;
		a2 += t*axy;

		c   = 1.0/sqrt(t*t + 1.0);
		s   = t*c;
		tau = s/(1.0 + c);

		v1[0] -= s*tau; v2[0] += s;
		v1[1] -= s;     v2[1] -= s*tau;
	}

	// sort principal values in descending order & permute principal directions
	if(a2 > a1)
	{
		a    = a1;    a1    = a2;    a2    = a;
		v[0] = v1[0]; v1[0] = v2[0]; v2[0] = v[0];
		v[1] = v1[1]; v1[1] = v2[1]; v2[1] = v[1];
	}

	(*pa1) = a1;
	(*pa2) = a2;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------


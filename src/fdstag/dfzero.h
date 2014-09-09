//---------------------------------------------------------------------------
//..................... DFZERO ROOT FINDING ALGORITHM .......................
//---------------------------------------------------------------------------
#ifndef __dfzero_h__
#define __dfzero_h__
//---------------------------------------------------------------------------
static inline PetscScalar FDMIN(PetscScalar a, PetscScalar b)
{
    if(a <= b) return a;
    else       return b;
}
//---------------------------------------------------------------------------
static inline PetscScalar FDMAX(PetscScalar a, PetscScalar b)
{
    if(a >= b) return a;
    else       return b;
}
//---------------------------------------------------------------------------
static inline PetscScalar FDSIGN(PetscScalar x, PetscScalar y)
{
    if(y >= 0.0) return  PetscAbsScalar(x);
    else         return -PetscAbsScalar(x);
}
//---------------------------------------------------------------------------
void DFZERO(
	PetscScalar (*F)(PetscScalar, void *), // nonlinear function with parameter and context
	void         *FCTX,                    // pointer to a function evaluation context
	PetscScalar  *_B,                      // left bound of root interval (output)
	PetscScalar  *_C,                      // right bound of root interval (output)
	PetscScalar    R,                      // initial guess
	PetscScalar    RE,                     // relative tolerance
	PetscScalar    AE,                     // absolute tolerance
	PetscInt     *_IFLAG);                 // error code (output)
//---------------------------------------------------------------------------
#endif

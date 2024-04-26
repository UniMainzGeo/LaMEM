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
//....................   FDSTAG VECTOR OUTPUT ROUTINES   ....................
//---------------------------------------------------------------------------
#ifndef __outFunct_h__
#define __outFunct_h__
//---------------------------------------------------------------------------

struct OutBuf;
struct JacRes;

//---------------------------------------------------------------------------
//...........  Multi-component output vector data structure .................
//---------------------------------------------------------------------------

typedef struct OutVec OutVec;

struct OutVec
{
	JacRes   *jr;
	OutBuf   *outbuf;
	PetscInt  ncomp;                        // number of components
	char      name      [_str_len_];        // output vector name
	PetscInt  phase_mask[_max_num_phases_]; // phase mask for phase aggregate
	PetscErrorCode (*OutVecWrite)(OutVec*); // output function pointer
};

void OutVecCreate(
	OutVec         *outvec,
	JacRes         *jr,
	OutBuf         *outbuf,
	const char     *name,
	const char     *label,
	PetscErrorCode (*OutVecWrite)(OutVec*),
	PetscInt        num,       // number of vector components or phases to aggregate
	PetscInt       *phase_ID); // phase IDs to aggregate

//---------------------------------------------------------------------------

PetscErrorCode PVOutWritePhase       (OutVec*);
PetscErrorCode PVOutWritePhaseAgg    (OutVec*);
PetscErrorCode PVOutWriteDensity     (OutVec*);
PetscErrorCode PVOutWriteViscTotal   (OutVec*);
PetscErrorCode PVOutWriteViscCreep   (OutVec*);
PetscErrorCode PVOutWriteVelocity    (OutVec*);
PetscErrorCode PVOutWritePressure    (OutVec*);
PetscErrorCode PVOutWriteGradient    (OutVec*);
PetscErrorCode PVOutWriteTotalPress  (OutVec*);
PetscErrorCode PVOutWriteEffPress    (OutVec*);
PetscErrorCode PVOutWriteOverPress   (OutVec*);
PetscErrorCode PVOutWriteLithoPress  (OutVec*);
PetscErrorCode PVOutWritePorePress   (OutVec*);
PetscErrorCode PVOutWriteTemperature (OutVec*);
PetscErrorCode PVOutWriteConductivity(OutVec*);
PetscErrorCode PVOutWriteDevStress   (OutVec*);
PetscErrorCode PVOutWriteJ2DevStress (OutVec*);
PetscErrorCode PVOutWriteStrainRate  (OutVec*);
PetscErrorCode PVOutWriteJ2StrainRate(OutVec*);
PetscErrorCode PVOutWriteMeltFraction(OutVec*);
PetscErrorCode PVOutWriteFluidDensity(OutVec*);
PetscErrorCode PVOutWriteVolRate     (OutVec*);
PetscErrorCode PVOutWriteVorticity   (OutVec*);
PetscErrorCode PVOutWriteAngVelMag   (OutVec*);
PetscErrorCode PVOutWriteTotStrain   (OutVec*);
PetscErrorCode PVOutWritePlastStrain (OutVec*);
PetscErrorCode PVOutWritePlastDissip (OutVec*);
PetscErrorCode PVOutWriteTotDispl    (OutVec*);
PetscErrorCode PVOutWriteSHmax       (OutVec*);
PetscErrorCode PVOutWriteStAngle     (OutVec*);
PetscErrorCode PVOutWriteEHmax       (OutVec*);
PetscErrorCode PVOutWriteYield       (OutVec*);
PetscErrorCode PVOutWriteRelDIIdif   (OutVec*);
PetscErrorCode PVOutWriteRelDIIdis   (OutVec*);
PetscErrorCode PVOutWriteRelDIIprl   (OutVec*);
PetscErrorCode PVOutWriteRelDIIpl    (OutVec*);
// === debug vectors ===============================================
PetscErrorCode PVOutWriteMomentRes   (OutVec*);
PetscErrorCode PVOutWriteContRes     (OutVec*);
PetscErrorCode PVOutWritEnergRes     (OutVec*);
PetscErrorCode PVOutWriteVelocityGr  (OutVec*);
PetscErrorCode PVOutWriteHeatSource  (OutVec*); // *djking
// ... add more output functions here

//---------------------------------------------------------------------------
#endif

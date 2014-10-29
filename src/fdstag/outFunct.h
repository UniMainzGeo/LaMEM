//---------------------------------------------------------------------------
//....................   FDSTAG VECTOR OUTPUT ROUTINES   ....................
//---------------------------------------------------------------------------
#ifndef __outFunct_h__
#define __outFunct_h__
//---------------------------------------------------------------------------

PetscErrorCode PVOutWritePhase       (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteDensity     (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteViscosity   (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteVelocity    (JacRes*, OutBuf*);
PetscErrorCode PVOutWritePressure    (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteTemperature (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteDevStress   (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteJ2DevStress (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteStrainRate  (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteJ2StrainRate(JacRes*, OutBuf*);
PetscErrorCode PVOutWriteVolRate     (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteVorticity   (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteAngVelMag   (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteTotStrain   (JacRes*, OutBuf*);
PetscErrorCode PVOutWritePlastStrain (JacRes*, OutBuf*);
PetscErrorCode PVOutWritePlastDissip (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteTotDispl    (JacRes*, OutBuf*);
// === debug	 vectors ===============================================
PetscErrorCode PVOutWriteMomentRes   (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteContRes     (JacRes*, OutBuf*);
PetscErrorCode PVOutWritEnergRes     (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteDII_CEN     (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteDII_XY      (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteDII_XZ      (JacRes*, OutBuf*);
PetscErrorCode PVOutWriteDII_YZ      (JacRes*, OutBuf*);

// ... add more output functions here

//---------------------------------------------------------------------------
#endif

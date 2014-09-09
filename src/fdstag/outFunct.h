//---------------------------------------------------------------------------
//....................   FDSTAG VECTOR OUTPUT ROUTINES   ....................
//---------------------------------------------------------------------------
#ifndef __outFunct_h__
#define __outFunct_h__
//---------------------------------------------------------------------------

PetscErrorCode PVOutWritePhase       (JacResCtx*, OutBuf*);
PetscErrorCode PVOutWriteDensity     (JacResCtx*, OutBuf*);
PetscErrorCode PVOutWriteViscosity   (JacResCtx*, OutBuf*);
PetscErrorCode PVOutWriteVelocity    (JacResCtx*, OutBuf*);
PetscErrorCode PVOutWritePressure    (JacResCtx*, OutBuf*);
PetscErrorCode PVOutWriteTemperature (JacResCtx*, OutBuf*);
PetscErrorCode PVOutWriteMomentRes   (JacResCtx*, OutBuf*);
PetscErrorCode PVOutWriteContRes     (JacResCtx*, OutBuf*);
PetscErrorCode PVOutWritEnergRes     (JacResCtx*, OutBuf*);
PetscErrorCode PVOutWriteDevStress   (JacResCtx*, OutBuf*);
PetscErrorCode PVOutWriteJ2DevStress (JacResCtx*, OutBuf*);
PetscErrorCode PVOutWriteStrainRate  (JacResCtx*, OutBuf*);
PetscErrorCode PVOutWriteJ2StrainRate(JacResCtx*, OutBuf*);
PetscErrorCode PVOutWriteVolRate     (JacResCtx*, OutBuf*);
PetscErrorCode PVOutWriteVorticity   (JacResCtx*, OutBuf*);
PetscErrorCode PVOutWriteAngVelMag   (JacResCtx*, OutBuf*);
PetscErrorCode PVOutWriteTotStrain   (JacResCtx*, OutBuf*);
PetscErrorCode PVOutWritePlastStrain (JacResCtx*, OutBuf*);
PetscErrorCode PVOutWritePlastDissip (JacResCtx*, OutBuf*);
PetscErrorCode PVOutWriteTotDispl    (JacResCtx*, OutBuf*);
PetscErrorCode PVOutWriteDII_CEN     (JacResCtx*, OutBuf*);
PetscErrorCode PVOutWriteDII_XY      (JacResCtx*, OutBuf*);
PetscErrorCode PVOutWriteDII_XZ      (JacResCtx*, OutBuf*);
PetscErrorCode PVOutWriteDII_YZ      (JacResCtx*, OutBuf*);

// ... add more output functions here

//---------------------------------------------------------------------------
#endif

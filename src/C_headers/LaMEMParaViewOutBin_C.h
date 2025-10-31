#ifndef LAMEMPARAVIEWOUTBIN_C_H
#define LAMEMPARAVIEWOUTBIN_C_H

#include <petsc.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

// Output mask structure for C interface
typedef struct {
    PetscInt phase, density, visc_total, visc_creep, velocity, pressure;
    PetscInt tot_pressure, gradient, eff_press, over_press, litho_press, pore_press;
    PetscInt temperature, conductivity, dev_stress, j2_dev_stress;
    PetscInt strain_rate, j2_strain_rate, vol_rate, vorticity, ang_vel_mag;
    PetscInt tot_strain, plast_strain, plast_dissip, tot_displ;
    PetscInt SHmax, StAngle, EHmax, yield;
    PetscInt DIIdif, DIIdis, DIIprl, DIIpl;
    PetscInt melt_fraction, fluid_density;
    PetscInt moment_res, cont_res, energ_res, vel_gr_tensor;
    PetscInt num_agg;
    char agg_name[16][256];  // max 16 aggregates, 256 chars each
    PetscInt agg_num_phase[16];
    PetscInt agg_phase_ID[16][32];  // max 32 phases per aggregate
} LaMEMOutMask;

// Output buffer functions
PetscErrorCode LaMEMParaViewOutBin_OutBufCreate(void *outbuf, void *jr);
PetscErrorCode LaMEMParaViewOutBin_OutBufDestroy(void *outbuf);
void LaMEMParaViewOutBin_OutBufConnectToFile(void *outbuf, FILE *fp);
void LaMEMParaViewOutBin_OutBufDump(void *outbuf);
void LaMEMParaViewOutBin_OutBufPutCoordVec(void *outbuf, void *ds, PetscScalar cf);
PetscErrorCode LaMEMParaViewOutBin_OutBufPut3DVecComp(void *outbuf, PetscInt ncomp, PetscInt dir, 
                                                      PetscScalar cf, PetscScalar shift);
PetscErrorCode LaMEMParaViewOutBin_OutBufZero3DVecComp(void *outbuf, PetscInt ncomp, PetscInt dir);

// Output mask functions
void LaMEMParaViewOutBin_OutMaskSetDefault(LaMEMOutMask *omask);
PetscInt LaMEMParaViewOutBin_OutMaskCountActive(LaMEMOutMask *omask);

// Main ParaView output functions
PetscErrorCode LaMEMParaViewOutBin_PVOutCreate(void *pvout, void *fb);
PetscErrorCode LaMEMParaViewOutBin_PVOutCreateData(void *pvout);
PetscErrorCode LaMEMParaViewOutBin_PVOutDestroy(void *pvout);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteTimeStep(void *pvout, const char *dirName, PetscScalar ttime);
PetscErrorCode LaMEMParaViewOutBin_PVOutWritePVTR(void *pvout, const char *dirName);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteVTR(void *pvout, const char *dirName);

// Output vector writing functions (large set of specialized functions)
PetscErrorCode LaMEMParaViewOutBin_PVOutWritePhase(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteDensity(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteViscTotal(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteViscCreep(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteVelocity(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWritePressure(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteTotalPress(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteGradient(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteEffPress(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteOverPress(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteLithoPress(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWritePorePress(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteTemperature(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteConductivity(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteDevStress(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteJ2DevStress(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteStrainRate(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteJ2StrainRate(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteVolRate(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteVorticity(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteAngVelMag(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteTotStrain(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWritePlastStrain(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWritePlastDissip(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteTotDispl(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteSHmax(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteStAngle(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteEHmax(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteYield(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteRelDIIdif(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteRelDIIdis(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteRelDIIprl(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteRelDIIpl(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteMeltFraction(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteFluidDensity(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteMomentRes(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteContRes(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWritEnergRes(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWriteVelocityGr(void *outvec);
PetscErrorCode LaMEMParaViewOutBin_PVOutWritePhaseAgg(void *outvec);

// Utility functions
void LaMEMParaViewOutBin_WriteXMLHeader(FILE *fp, const char *file_type);
PetscErrorCode LaMEMParaViewOutBin_UpdatePVDFile(const char *dirName, const char *outfile, 
                                                  const char *ext, long int *offset, 
                                                  PetscScalar ttime, PetscInt outpvd);

#ifdef __cplusplus
}
#endif

#endif
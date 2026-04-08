#include "LaMEMParaViewOutBin_C.h"
#include "../LaMEM.h"
#include "../paraViewOutBin.h"
#include "../scaling.h"
#include "../parsing.h"
#include "../fdstag.h"
#include "../JacRes.h"
#include "../phase.h"
#include "../outFunct.h"
#include "../tools.h"

extern "C" {

// Output buffer functions
PetscErrorCode LaMEMParaViewOutBin_OutBufCreate(void *outbuf, void *jr) {
    return OutBufCreate((OutBuf*)outbuf, (JacRes*)jr);
}

PetscErrorCode LaMEMParaViewOutBin_OutBufDestroy(void *outbuf) {
    return OutBufDestroy((OutBuf*)outbuf);
}

void LaMEMParaViewOutBin_OutBufConnectToFile(void *outbuf, FILE *fp) {
    OutBufConnectToFile((OutBuf*)outbuf, fp);
}

void LaMEMParaViewOutBin_OutBufDump(void *outbuf) {
    OutBufDump((OutBuf*)outbuf);
}

void LaMEMParaViewOutBin_OutBufPutCoordVec(void *outbuf, void *ds, PetscScalar cf) {
    OutBufPutCoordVec((OutBuf*)outbuf, (Discret1D*)ds, cf);
}

PetscErrorCode LaMEMParaViewOutBin_OutBufPut3DVecComp(void *outbuf, PetscInt ncomp, PetscInt dir, 
                                                      PetscScalar cf, PetscScalar shift) {
    return OutBufPut3DVecComp((OutBuf*)outbuf, ncomp, dir, cf, shift);
}

PetscErrorCode LaMEMParaViewOutBin_OutBufZero3DVecComp(void *outbuf, PetscInt ncomp, PetscInt dir) {
    return OutBufZero3DVecComp((OutBuf*)outbuf, ncomp, dir);
}

// Output mask functions
void LaMEMParaViewOutBin_OutMaskSetDefault(LaMEMOutMask *omask) {
    // Set C struct defaults matching OutMaskSetDefault
    omask->phase = 1;
    omask->visc_total = 1;
    omask->visc_creep = 1;
    omask->velocity = 1;
    omask->pressure = 1;
    // All other fields default to 0
    omask->density = 0;
    omask->tot_pressure = 0;
    omask->gradient = 0;
    omask->eff_press = 0;
    omask->over_press = 0;
    omask->litho_press = 0;
    omask->pore_press = 0;
    omask->temperature = 0;
    omask->conductivity = 0;
    omask->dev_stress = 0;
    omask->j2_dev_stress = 0;
    omask->strain_rate = 0;
    omask->j2_strain_rate = 0;
    omask->vol_rate = 0;
    omask->vorticity = 0;
    omask->ang_vel_mag = 0;
    omask->tot_strain = 0;
    omask->plast_strain = 0;
    omask->plast_dissip = 0;
    omask->tot_displ = 0;
    omask->SHmax = 0;
    omask->StAngle = 0;
    omask->EHmax = 0;
    omask->yield = 0;
    omask->DIIdif = 0;
    omask->DIIdis = 0;
    omask->DIIprl = 0;
    omask->DIIpl = 0;
    omask->melt_fraction = 0;
    omask->fluid_density = 0;
    omask->moment_res = 0;
    omask->cont_res = 0;
    omask->energ_res = 0;
    omask->vel_gr_tensor = 0;
    omask->num_agg = 0;
}

PetscInt LaMEMParaViewOutBin_OutMaskCountActive(LaMEMOutMask *omask) {
    PetscInt cnt = 0;
    if(omask->phase) cnt++;
    if(omask->density) cnt++;
    if(omask->visc_total) cnt++;
    if(omask->visc_creep) cnt++;
    if(omask->velocity) cnt++;
    if(omask->pressure) cnt++;
    if(omask->tot_pressure) cnt++;
    if(omask->gradient) cnt++;
    if(omask->eff_press) cnt++;
    if(omask->over_press) cnt++;
    if(omask->litho_press) cnt++;
    if(omask->pore_press) cnt++;
    if(omask->temperature) cnt++;
    if(omask->conductivity) cnt++;
    if(omask->dev_stress) cnt++;
    if(omask->j2_dev_stress) cnt++;
    if(omask->strain_rate) cnt++;
    if(omask->j2_strain_rate) cnt++;
    if(omask->vol_rate) cnt++;
    if(omask->vorticity) cnt++;
    if(omask->ang_vel_mag) cnt++;
    if(omask->tot_strain) cnt++;
    if(omask->plast_strain) cnt++;
    if(omask->plast_dissip) cnt++;
    if(omask->tot_displ) cnt++;
    if(omask->SHmax) cnt++;
    if(omask->StAngle) cnt++;
    if(omask->EHmax) cnt++;
    if(omask->yield) cnt++;
    if(omask->DIIdif) cnt++;
    if(omask->DIIdis) cnt++;
    if(omask->DIIprl) cnt++;
    if(omask->DIIpl) cnt++;
    if(omask->melt_fraction) cnt++;
    if(omask->fluid_density) cnt++;
    if(omask->moment_res) cnt++;
    if(omask->cont_res) cnt++;
    if(omask->energ_res) cnt++;
    if(omask->vel_gr_tensor) cnt++;
    cnt += omask->num_agg;
    return cnt;
}

// Main ParaView output functions
PetscErrorCode LaMEMParaViewOutBin_PVOutCreate(void *pvout, void *fb) {
    return PVOutCreate((PVOut*)pvout, (FB*)fb);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutCreateData(void *pvout) {
    return PVOutCreateData((PVOut*)pvout);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutDestroy(void *pvout) {
    return PVOutDestroy((PVOut*)pvout);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteTimeStep(void *pvout, const char *dirName, PetscScalar ttime) {
    return PVOutWriteTimeStep((PVOut*)pvout, dirName, ttime);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWritePVTR(void *pvout, const char *dirName) {
    return PVOutWritePVTR((PVOut*)pvout, dirName);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteVTR(void *pvout, const char *dirName) {
    return PVOutWriteVTR((PVOut*)pvout, dirName);
}

// Output vector writing functions (all the specialized field writers)
PetscErrorCode LaMEMParaViewOutBin_PVOutWritePhase(void *outvec) {
    return PVOutWritePhase((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteDensity(void *outvec) {
    return PVOutWriteDensity((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteViscTotal(void *outvec) {
    return PVOutWriteViscTotal((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteViscCreep(void *outvec) {
    return PVOutWriteViscCreep((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteVelocity(void *outvec) {
    return PVOutWriteVelocity((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWritePressure(void *outvec) {
    return PVOutWritePressure((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteTotalPress(void *outvec) {
    return PVOutWriteTotalPress((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteGradient(void *outvec) {
    return PVOutWriteGradient((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteEffPress(void *outvec) {
    return PVOutWriteEffPress((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteOverPress(void *outvec) {
    return PVOutWriteOverPress((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteLithoPress(void *outvec) {
    return PVOutWriteLithoPress((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWritePorePress(void *outvec) {
    return PVOutWritePorePress((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteTemperature(void *outvec) {
    return PVOutWriteTemperature((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteConductivity(void *outvec) {
    return PVOutWriteConductivity((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteDevStress(void *outvec) {
    return PVOutWriteDevStress((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteJ2DevStress(void *outvec) {
    return PVOutWriteJ2DevStress((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteStrainRate(void *outvec) {
    return PVOutWriteStrainRate((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteJ2StrainRate(void *outvec) {
    return PVOutWriteJ2StrainRate((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteVolRate(void *outvec) {
    return PVOutWriteVolRate((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteVorticity(void *outvec) {
    return PVOutWriteVorticity((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteAngVelMag(void *outvec) {
    return PVOutWriteAngVelMag((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteTotStrain(void *outvec) {
    return PVOutWriteTotStrain((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWritePlastStrain(void *outvec) {
    return PVOutWritePlastStrain((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWritePlastDissip(void *outvec) {
    return PVOutWritePlastDissip((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteTotDispl(void *outvec) {
    return PVOutWriteTotDispl((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteSHmax(void *outvec) {
    return PVOutWriteSHmax((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteStAngle(void *outvec) {
    return PVOutWriteStAngle((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteEHmax(void *outvec) {
    return PVOutWriteEHmax((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteYield(void *outvec) {
    return PVOutWriteYield((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteRelDIIdif(void *outvec) {
    return PVOutWriteRelDIIdif((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteRelDIIdis(void *outvec) {
    return PVOutWriteRelDIIdis((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteRelDIIprl(void *outvec) {
    return PVOutWriteRelDIIprl((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteRelDIIpl(void *outvec) {
    return PVOutWriteRelDIIpl((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteMeltFraction(void *outvec) {
    return PVOutWriteMeltFraction((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteFluidDensity(void *outvec) {
    return PVOutWriteFluidDensity((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteMomentRes(void *outvec) {
    return PVOutWriteMomentRes((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteContRes(void *outvec) {
    return PVOutWriteContRes((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWritEnergRes(void *outvec) {
    return PVOutWritEnergRes((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWriteVelocityGr(void *outvec) {
    return PVOutWriteVelocityGr((OutVec*)outvec);
}

PetscErrorCode LaMEMParaViewOutBin_PVOutWritePhaseAgg(void *outvec) {
    return PVOutWritePhaseAgg((OutVec*)outvec);
}

// Utility functions
void LaMEMParaViewOutBin_WriteXMLHeader(FILE *fp, const char *file_type) {
    WriteXMLHeader(fp, file_type);
}

PetscErrorCode LaMEMParaViewOutBin_UpdatePVDFile(const char *dirName, const char *outfile, 
                                                  const char *ext, long int *offset, 
                                                  PetscScalar ttime, PetscInt outpvd) {
    return UpdatePVDFile(dirName, outfile, ext, offset, ttime, outpvd);
}

}
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
 **    filename:   paraViewOutBin.h
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
//.................   FDSTAG PARAVIEW XML OUTPUT ROUTINES   .................
//---------------------------------------------------------------------------
// All output fields are interpolated onto the corner nodes before output.
// Each processor includes local number of nodes in every spatial direction
// plus one overlapping ghost node from the next processor. Last processor
// doesn't have a ghost node. Every output field is copied into single precision
// buffer defined on the local output nodes. Every components of the vector
// and tensor filed is scaled (if necessary) and arranged accordingly in the
// buffer: x, y, z components for the vector fields, and xx, yy, zz, xy, yz, xz
// components for the tensor filed (diagonal format). When buffer is arranged
// it is written to the output file. Prerequisite for the copying to the buffer
// is to have every component in the LOCAL corner vector, which is obtained
// by usual scattering from the GLOBAL corner vector. Interpolation to the GLOBAL
// corner vector is done from the LOCAL source vectors (center, edges, or faces).
// These vectors are also obtained by global-to-local scattering. Some vectors
// (velocity and momentum residual) are assumed to be available in LOCAL format.
// The overall scheme is as follows:
//    * loop over components (e.g. xx, yy, ... for stress tensor)
//       - copy field from context to global vector (only for center or edge)
//       - global-to-local scatter                  (only for center or edge)
//       - interpolate from local source vector to global corner vector
//       - scatter from global-to-local corner vectors
//       - scale and copy component to the buffer from local corner vector
//    * and of loop
//    * dump buffer to output file
//---------------------------------------------------------------------------
#ifndef __paraViewOutBin_h__
#define __paraViewOutBin_h__
//---------------------------------------------------------------------------

// maximum number of components in the output vector
#define _max_num_comp_ 9

//---------------------------------------------------------------------------

struct FB;
struct FDSTAG;
struct JacRes;
struct Discret1D;

//---------------------------------------------------------------------------
//............................. Output buffer ...............................
//---------------------------------------------------------------------------
struct OutBuf
{
	FDSTAG   *fs;    // staggered grid layout
	FILE     *fp;    // output file handler
	float    *buff;  // direct output buffer
	PetscInt  cn;    // current number of elements in the buffer

	// grid buffer vectors
	Vec lbcen, lbcor, lbxy, lbxz, lbyz; // local (ghosted)

};
//---------------------------------------------------------------------------
PetscErrorCode OutBufCreate(OutBuf *outbuf, JacRes *jr);

PetscErrorCode OutBufDestroy(OutBuf *outbuf);

void OutBufConnectToFile(OutBuf  *outbuf, FILE *fp);

// dump output buffer contents to disk
void OutBufDump(OutBuf  *outbuf);

// put FDSTAG coordinate vector to output buffer
void OutBufPutCoordVec(
	OutBuf      *outbuf,
	Discret1D   *ds,
	PetscScalar  cf); // scaling coefficient

// put component of 3D vector to output buffer
PetscErrorCode OutBufPut3DVecComp(
	OutBuf      *outbuf,
	PetscInt     ncomp,  // number of components
	PetscInt     dir,    // component identifier
	PetscScalar  cf,     // scaling coefficient
	PetscScalar  shift); // shift parameter (subtracted from scaled values)

PetscErrorCode OutBufZero3DVecComp(
	OutBuf      *outbuf,
	PetscInt     ncomp,  // number of components
	PetscInt     dir);   // component identifier

//---------------------------------------------------------------------------
// ...................  Vector output function pointer ......................
//---------------------------------------------------------------------------

typedef PetscErrorCode (*OutVecFunctPtr)(JacRes*, OutBuf*);

//---------------------------------------------------------------------------
//...........  Multi-component output vector data structure .................
//---------------------------------------------------------------------------
struct OutVec
{
	char          *name;        // output vector name
	OutVecFunctPtr OutVecFunct; // pointer to vector output function
	PetscInt       ncomp;       // number of components

};
//---------------------------------------------------------------------------
void OutVecCreate(
	OutVec         *outvec,
	const char     *name,
	const char     *label,
	OutVecFunctPtr  OutVecFunct,
	PetscInt        ncomp);

void OutVecDestroy(OutVec *outvec);

//---------------------------------------------------------------------------
//.......................... Vector output mask .............................
//---------------------------------------------------------------------------
struct OutMask
{
	PetscInt phase;          // phase
	PetscInt density;        // density
	PetscInt visc_total;     // total effective viscosity
	PetscInt visc_creep;     // creep effective viscosity
	PetscInt visc_plast;     // viscoplastic viscosity
	PetscInt velocity;       // velocity
	PetscInt pressure;       // pressure
	PetscInt eff_press;      // effective pressure
	PetscInt over_press;     // overpressure
	PetscInt litho_press;    // lithostatic pressure
	PetscInt pore_press;     // pore pressure
	PetscInt temperature;    // temperature
	PetscInt dev_stress;     // deviatoric stress tensor
	PetscInt j2_dev_stress;  // deviatoric stress second invariant
	PetscInt strain_rate;    // deviatoric strain rate tensor
	PetscInt j2_strain_rate; // deviatoric strain rate second invariant
	PetscInt vol_rate;       // volumetric strain rate
	PetscInt vorticity;      // vorticity vector
	PetscInt ang_vel_mag;    // average angular velocity magnitude
	PetscInt tot_strain;     // total strain
	PetscInt plast_strain;   // accumulated plastic strain
	PetscInt plast_dissip;   // plastic dissipation
	PetscInt tot_displ;      // total displacements
	PetscInt SHmax;          // maximum horizontal stress
	PetscInt EHmax;          // maximum horizontal extension
	PetscInt ISA;            // Infinite Strain Axis
	PetscInt GOL;            // Grain Orientation Lag
	PetscInt yield;          // yield stress
	// === debugging vectors ===============================================
	PetscInt moment_res;     // momentum residual
	PetscInt cont_res;       // continuity residual
	PetscInt energ_res;      // energy residual

};
//---------------------------------------------------------------------------

void OutMaskSetDefault(OutMask *omask);

PetscInt OutMaskCountActive(OutMask *omask);

//---------------------------------------------------------------------------
//...................... ParaView output driver object ......................
//---------------------------------------------------------------------------
struct PVOut
{
	JacRes   *jr;
	char      outfile[_STR_LEN_]; // output file name
	OutMask   omask;              // output vector mask
	PetscInt  nvec;               // number of output vectors
	OutVec   *outvecs;            // output vectors
	OutBuf    outbuf;             // output buffer
	long int  offset;             // pvd file offset
	PetscInt  outpvd;             // pvd file output flag

};
//---------------------------------------------------------------------------

// create ParaView output driver
PetscErrorCode PVOutCreate(PVOut *pvout, FB *fb);

// create output buffer and vectors
PetscErrorCode PVOutCreateData(PVOut *pvout);

// destroy ParaView output driver
PetscErrorCode PVOutDestroy(PVOut *pvout);

// write all time-step output files to disk (PVD, PVTR, VTR)
PetscErrorCode PVOutWriteTimeStep(PVOut *pvout, const char *dirName, PetscScalar ttime);

// write parallel PVTR file (called every time step on first processor)
// WARNING! this is potential bottleneck, get rid of writing every time-step
PetscErrorCode PVOutWritePVTR(PVOut *pvout, const char *dirName);

// write sequential VTR files on every processor (called every time step)
PetscErrorCode PVOutWriteVTR(PVOut *pvout, const char *dirName);

//---------------------------------------------------------------------------
//........................... Service Functions .............................
//---------------------------------------------------------------------------

// Add standard header to output file
void WriteXMLHeader(FILE *fp, const char *file_type);

// update PVD file (called every time step on first processor)
// WARNING! this is potential bottleneck, get rid of writing every time-step
PetscErrorCode UpdatePVDFile(
		const char *dirName, const char *outfile, const char *ext,
		long int *offset, PetscScalar ttime, PetscInt *outpvd);

//---------------------------------------------------------------------------
#endif

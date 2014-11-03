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
#define _timestep_buff_size_ 4096
// maximum number of components in the output vector
#define _max_num_comp_ 6
//---------------------------------------------------------------------------
//............................. Output buffer ...............................
//---------------------------------------------------------------------------
typedef struct
{
	FDSTAG   *fs;    // staggered grid layout
	FILE     *fp;    // output file handler
	float    *buff;  // direct output buffer
	PetscInt  cn;    // current number of elements in the buffer
	// grid buffer vectors
	Vec gbcen, gbcor, gbxy, gbxz, gbyz; // global
	Vec lbcen, lbcor, lbxy, lbxz, lbyz; // local (ghosted)

} OutBuf;
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
	PetscInt     ncomp, // number of components
	PetscInt     dir,   // component identifier
	PetscScalar  cf);   // scaling coefficient

//---------------------------------------------------------------------------
// ...................  Vector output function pointer ......................
//---------------------------------------------------------------------------

typedef PetscErrorCode (*OutVecFunctPtr)(JacRes*, OutBuf*);

//---------------------------------------------------------------------------
//...........  Multi-component output vector data structure .................
//---------------------------------------------------------------------------
typedef struct
{
	char          *name;        // output vector name
	OutVecFunctPtr OutVecFunct; // pointer to vector output function
	PetscInt       ncomp;       // number of components

} OutVec;
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
typedef struct
{
	PetscInt phase;          // phase
	PetscInt density;        // density
	PetscInt viscosity;      // effective viscosity
	PetscInt velocity;       // velocity
	PetscInt pressure;       // pressure
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
	PetscInt phrat[max_num_phases]; // phase ratios
	// === debugging vectors ===============================================
	PetscInt moment_res;     // momentum residual
	PetscInt cont_res;       // continuity residual
	PetscInt energ_res;      // energy residual
	PetscInt DII_CEN;        // effective strain rate invariant in center
	PetscInt DII_XY;         // effective strain rate invariant on xy-edge
	PetscInt DII_XZ;         // effective strain rate invariant on xz-edge
	PetscInt DII_YZ;         // effective strain rate invariant on yz-edge

	// ... add more output vector identifiers here

} OutMask;
//---------------------------------------------------------------------------

void OutMaskSetDefault(OutMask *omask);

PetscInt OutMaskCountActive(OutMask *omask);

//---------------------------------------------------------------------------
//...................... ParaView output driver object ......................
//---------------------------------------------------------------------------
typedef struct
{
	char        *outfile; // output file name
	PetscScalar  crdScal; // output scaling for coordinates
	OutMask      omask;   // output vector mask
	PetscInt     nvec;    // number of output vectors
	OutVec      *outvecs; // output vectors
	OutBuf       outbuf;  // output buffer
	//================================
	// time step buffer
	// to be replaced by python script
	//================================
	PetscInt     nstep;   // number of collected time steps
	PetscInt     bcnt;    // buffer counter
	char        *dbname;  // time step database name
	PetscScalar *bstamps; // buffer for output time stamps
	PetscInt    *bindexs; // buffer for indices of saved time steps

} PVOut;
//---------------------------------------------------------------------------

// create ParaView output driver
PetscErrorCode PVOutCreate(PVOut *pvout, JacRes *jr, const char *filename);

// destroy ParaView output driver
PetscErrorCode PVOutDestroy(PVOut *pvout);

// Add standard header to output file
void PVOutWriteXMLHeader(FILE *fp, const char *file_type);

// write all time-step output files to disk (PVTR, VTR) and update time step buffer data
PetscErrorCode PVOutWriteTimeStep(PVOut *pvout, JacRes *jr, PetscScalar ttime, PetscInt tindx);

// store time stamp and step index to the buffer
PetscErrorCode PVOutUpdateTimeStepBuffer(PVOut *pvout, PetscScalar ttime, PetscInt tindx);

// dump time-step buffer to disk
PetscErrorCode PVOutDumpTimeStepBuffer(PVOut *pvout);

// write final PVD file (time step collection) (called once per simulation)
PetscErrorCode PVOutWritePVD(PVOut *pvout);

// write parallel PVTR file (called every time step on first processor)
// WARNING! this is potential bottleneck, figure out how to get rid of writing every time-step
PetscErrorCode PVOutWritePVTR(PVOut *pvout, PetscInt tindx);

// write sequential VTR files on every processor (called every time step)
PetscErrorCode PVOutWriteVTR(PVOut *pvout, JacRes *jr, PetscInt tindx);

//---------------------------------------------------------------------------

#endif

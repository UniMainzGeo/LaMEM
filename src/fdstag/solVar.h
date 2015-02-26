//---------------------------------------------------------------------------
//.................   SOLUTION AND HISTORY VARIABLES   ......................
//---------------------------------------------------------------------------
#ifndef __solVar_h__
#define __solVar_h__

//---------------------------------------------------------------------------
//....    Non-Symmetric second rank tensor (gradient & rotation tensors) ....
//---------------------------------------------------------------------------

typedef struct
{
	PetscScalar xx, xy, xz;
	PetscScalar yx, yy, yz;
	PetscScalar zx, zy, zz;

} Tensor2RN;

//---------------------------------------------------------------------------
//.......   Symmetric second rank tensor (stress & strain tensors)   ........
//---------------------------------------------------------------------------

typedef struct
{
	PetscScalar xx, xy, xz;
	PetscScalar     yy, yz;
	PetscScalar         zz;

} Tensor2RS;

//---------------------------------------------------------------------------
//............   Material marker (history variables advection)   ............
//---------------------------------------------------------------------------

typedef struct
{
	PetscInt    phase; // phase identifier
	PetscScalar X[3];  // global coordinates
	PetscScalar p;     // pressure
	PetscScalar T;     // temperature
	PetscScalar APS;   // accumulated plastic strain
	Tensor2RS   S;     // deviatoric stress

} Marker;

//---------------------------------------------------------------------------
//.....................   Deviatoric solution variables   ...................
//---------------------------------------------------------------------------

typedef struct
{
	PetscScalar  DII;   // effective strain rate
	PetscScalar  eta;   // effective tangent viscosity
	PetscScalar  I2Gdt; // inverse elastic viscosity (1/2G/dt)
	PetscScalar  Hr;    // shear heating term (partial)
	PetscScalar  DIIpl; // plastic strain rate
	PetscScalar  APS;   // accumulated plastic strain
	PetscScalar  PSR;   // plastic strain-rate contribution

} SolVarDev;

//---------------------------------------------------------------------------
//.....................   Volumetric solution variables   ...................
//---------------------------------------------------------------------------

typedef struct
{
	PetscScalar  theta; // volumetric strain rate
	PetscScalar  rho;   // strain- & temperature-dependent density
	PetscScalar  IKdt;  // inverse bulk elastic viscosity (1/K/dt)
	PetscScalar  alpha; // effective thermal expansion
	PetscScalar  Tn;    // history temperature
	PetscScalar  pn;    // history pressure

} SolVarBulk;

//---------------------------------------------------------------------------
//........................   Cell solution variables   ......................
//---------------------------------------------------------------------------

typedef struct
{
	SolVarDev    svDev;         // deviatoric variables
	SolVarBulk   svBulk;        // volumetric variables
	PetscScalar  sxx, syy, szz; // deviatoric stress
	PetscScalar  hxx, hyy, hzz; // history stress (elastic)
	PetscScalar  dxx, dyy, dzz; // total deviatoric strain rate
	PetscScalar *phRat;         // phase ratios in the control volume

} SolVarCell;

//---------------------------------------------------------------------------
//........................   Edge solution variables   ......................
//---------------------------------------------------------------------------

typedef struct
{
	SolVarDev    svDev; // deviatoric variables
	PetscScalar  s;     // xy, xz, yz deviatoric stress components
	PetscScalar  h;     // xy, xz, yz history stress components (elastic)
	PetscScalar  d;     // xy, xz, yz total deviatoric strain rate components
	PetscScalar  ws;    // normalization for distance-dependent interpolation
	PetscScalar *phRat; // phase ratios in the control volume

} SolVarEdge;

//---------------------------------------------------------------------------
//.......................   Softening Law Parameters  .......................
//---------------------------------------------------------------------------

typedef struct
{
	PetscInt    ID;
	PetscScalar APS1; // begin of softening APS
	PetscScalar APS2; // end of softening APS
	PetscScalar A;    // reduction ratio

} Soft_t;

//---------------------------------------------------------------------------
//......................   Material parameter table   .......................
//---------------------------------------------------------------------------

typedef struct
{
	PetscInt     ID;      // material ID
	// physical parameters
	PetscScalar  rho;     // reference density
	// elasticity parameters
	PetscScalar  K;       // bulk modulus
	PetscScalar  Kp;      // pressure dependence parameter
	PetscScalar  G;       // shear modulus
	// diffusion creep parameters
	PetscScalar  Bd;      // pre-exponential constant
	PetscScalar  Ed;      // activation energy
	PetscScalar  Vd;      // activation volume
	// dislocation creep parameters
	PetscScalar  Bn;      // pre-exponential constant
	PetscScalar  n;       // power law exponent
	PetscScalar  En;      // activation energy
	PetscScalar  Vn;      // activation volume
	// Peierls creep parameters
	PetscScalar  Bp;      // pre-exponential constant
	PetscScalar  Ep;      // activation energy
	PetscScalar  Vp;      // activation volume
	PetscScalar  taup;    // scaling stress
	PetscScalar  gamma;   // approximation parameter
	PetscScalar  q;       // stress-dependence parameter
	// plasticity parameters
	PetscScalar  fr;      // friction coefficient
	PetscScalar  ch;      // cohesion
	PetscInt     chSoftID;// id of softening law
	PetscInt     frSoftID;// id of softening law
	Soft_t      *frSoft;  // friction softening law parameters
	Soft_t      *chSoft;  // cohesion softening law parameters
	// thermal parameters
	PetscScalar  alpha;   // thermal expansivity
	PetscScalar  Cp;      // cpecific heat (capacity)
	PetscScalar  k;       // thermal conductivity
	PetscScalar  A;       // radiogenic heat production

} Material_t;

//---------------------------------------------------------------------------
//.....................   Material parameter limiters .......................
//---------------------------------------------------------------------------

typedef struct
{
	// viscosity limits (and inverses)
	PetscScalar eta_min;
	PetscScalar eta_max;
	// reference viscosity (initial guess)
	PetscScalar eta_ref;
	// reference temperature
	PetscScalar TRef;
	// universal gas constant
	PetscScalar Rugc;
	// viscosity & strain-rate tolerances
	PetscScalar eta_atol; // viscosity absolute tolerance
	PetscScalar eta_rtol; // viscosity relative tolerance
	PetscScalar DII_atol; // strain rate absolute tolerance
	PetscScalar DII_rtol; // strain rate relative tolerance
	// background (reference) strain-rate
	PetscScalar DII_ref;
	// plasticity parameters limits
	PetscScalar minCh;  // minimum cohesion
	PetscScalar minFr;  // maximum friction
	PetscScalar tauUlt; // ultimate yield stress
	// thermo-mechanical coupling controls
	PetscScalar shearHeatEff; // shear heating efficiency parameter [0 - 1]
	// rheology controls
	PetscBool   quasiHarmAvg; // plasticity quasi-harmonic averaging flag
	PetscBool   initGuessFlg; // initial guess computation flag

} MatParLim;

//---------------------------------------------------------------------------
#endif

/*




static inline void Tensor2RNClear(Tensor2RN *A)
{
	// set all components to a constant (A_ij = 0.0)
	A->_11 = 0.0; A->_12 = 0.0; A->_13 = 0.0;
	A->_21 = 0.0; A->_22 = 0.0; A->_23 = 0.0;
	A->_31 = 0.0; A->_32 = 0.0; A->_33 = 0.0;
}

static inline void Tensor2RNUnit(Tensor2RN *A)
{
	// initialize Kronecker's delta (A_ij = delta_ij)
	A->_11 = 1.0; A->_12 = 0.0; A->_13 = 0.0;
	A->_21 = 0.0; A->_22 = 1.0; A->_23 = 0.0;
	A->_31 = 0.0; A->_32 = 0.0; A->_33 = 1.0;
}

static inline void Tensor2RNCopy(Tensor2RN *A, Tensor2RN *B)
{
	// copy one tensor to another (A = B)
	A->_11 = B->_11;  A->_12 = B->_12;  A->_13 = B->_13;
	A->_21 = B->_21;  A->_22 = B->_22;  A->_23 = B->_23;
	A->_31 = B->_31;  A->_32 = B->_32;  A->_33 = B->_33;
}
//---------------------------------------------------------------------------
inline void Tensor2RNScale(Tensor2RN *A, PetscScalar k)
{
	// set all components with a constant (A_ij = A_ij*k)
	A->_11 *= k; A->_12 *= k; A->_13 *= k;
	A->_21 *= k; A->_22 *= k; A->_23 *= k;
	A->_31 *= k; A->_32 *= k; A->_33 *= k;
}

//---------------------------------------------------------------------------
static inline void Tensor2RNUpdate(Tensor2RN *A, Tensor2RN *B, PetscScalar kb)
{
	// copy scaled tensor to another (A = A + kb*B)
	A->_11 += kb*B->_11;  A->_12 += kb*B->_12;  A->_13 += kb*B->_13;
	A->_21 += kb*B->_21;  A->_22 += kb*B->_22;  A->_23 += kb*B->_23;
	A->_31 += kb*B->_31;  A->_32 += kb*B->_32;  A->_33 += kb*B->_33;
}
//---------------------------------------------------------------------------
static inline void Tensor2RNCombine(Tensor2RN *A, Tensor2RN *B, Tensor2RN *C, PetscScalar kb, PetscScalar kc)
{
	// form a linear combination of two tensors (A = kb*B + kc*C)
	A->_11 = kb*B->_11 + kc*C->_11;  A->_12 = kb*B->_12 + kc*C->_12;  A->_13 = kb*B->_13 + kc*C->_13;
	A->_21 = kb*B->_21 + kc*C->_21;  A->_22 = kb*B->_22 + kc*C->_22;  A->_23 = kb*B->_23 + kc*C->_23;
	A->_31 = kb*B->_31 + kc*C->_31;  A->_32 = kb*B->_32 + kc*C->_32;  A->_33 = kb*B->_33 + kc*C->_33;
}

//---------------------------------------------------------------------------
static inline void Tensor2RNCompGrad(Tensor2RN *A,
	PetscScalar *ax, PetscScalar *ay, PetscScalar *az,
	PetscScalar  dx, PetscScalar  dy, PetscScalar  dz)
{
	// Compute gradient of a vector field given by component vectors ax, ay, az.
	// Each component vector is defined on a 6-point stencil with arrangement:
	// points 0 - 1 are aligned along x axis
	// points 2 - 3 are aligned along y axis
	// points 4 - 5 are aligned along z axis
	A->_11 = (ax[1] - ax[0])/dx;  A->_12 = (ax[3] - ax[2])/dy;  A->_13 = (ax[5] - ax[4])/dz;
	A->_21 = (ay[1] - ay[0])/dx;  A->_22 = (ay[3] - ay[2])/dy;  A->_23 = (ay[5] - ay[4])/dz;
	A->_31 = (az[1] - az[0])/dx;  A->_32 = (az[3] - az[2])/dy;  A->_33 = (az[5] - az[4])/dz;
}
//---------------------------------------------------------------------------
static inline PetscScalar Tensor2RNDet(Tensor2RN *A)
{
	// compute determinant of a tensor
	return A->_11*(A->_22*A->_33 - A->_23*A->_32)
	+      A->_12*(A->_23*A->_31 - A->_21*A->_33)
	+      A->_13*(A->_21*A->_32 - A->_22*A->_31);
}

static inline void Tensor2RSClear(Tensor2RS *A)
{
	// set all components to a constant (A_ij = 0.0)
	A->_11 = 0.0;
	A->_12 = 0.0; A->_22 = 0.0;
	A->_13 = 0.0; A->_23 = 0.0; A->_33 = 0.0;
}

static inline void Tensor2RSUpdate(Tensor2RS *A, Tensor2RS *B, PetscScalar kb)
{
	// copy scaled tensor to another (A = A + kb*B)
	A->_11 += kb*B->_11;
	A->_12 += kb*B->_12;  A->_22 += kb*B->_22;
	A->_13 += kb*B->_13;  A->_23 += kb*B->_23;  A->_33 += kb*B->_33;
}


*/


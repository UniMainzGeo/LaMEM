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
//..............   USER-DEFINED PHASE TRANSITION PLUGIN (Julia)   ...........
//---------------------------------------------------------------------------
#include "LaMEM.h"
#include "Tensor.h"
#include "advect.h"
#include "scaling.h"
#include "JacRes.h"
#include "fdstag.h"
#include "bc.h"
#include "tssolve.h"
#include "phase.h"
#include "phase_transition_plugin.h"
#include <cstddef>
#include <cstdint>
//---------------------------------------------------------------------------
// C ABI of the plugin function (see phase_transition_plugin.h for the
// full, documented signature). Note: phase_in/phase_out are `int` (Cint,
// 32-bit) by ABI contract, NOT PetscInt (64-bit in Int64 PETSc builds).
typedef int (*PhTrPluginFn)(
	size_t n,
	double *x,  double *y,  double *z,
	double *T,  double *p,
	double  time,
	double *sxx, double *syy, double *szz,
	double *sxy, double *sxz, double *syz,
	double *j2_stress_cell, double *j2_strainrate_cell,
	double *eta_cell,       double *aps_cell,
	int    *phase_in, int *phase_out);

// function pointer signature of jl_parse_opts(int *argc, char ***argv)
typedef void (*JlParseOptsFn)(int *argc, char ***argv);

// function pointer signature of jl_init_with_image_handle(void *handle)
typedef void (*JlInitWithImageHandleFn)(void *handle);

// function pointer signature of jl_atexit_hook(int status)
typedef void (*JlAtexitHookFn)(int status);

//---------------------------------------------------------------------------
// module-local state (Phase 1: single global plugin instance, one per
// process. See the header comment on why this is process-wide, not
// per-LaMEMLibSolve-call, state.)
namespace
{
	PetscBool     initTried  = PETSC_FALSE; // PhTrPluginLoad already attempted in this process?
	PetscBool     active     = PETSC_FALSE; // plugin successfully loaded?
	PetscDLHandle handle     = NULL;         // handle of the plugin library
	PhTrPluginFn  fn         = NULL;         // lamem_phase_transition
	JlAtexitHookFn atexitFn  = NULL;         // jl_atexit_hook (optional)

	// SoA scratch buffers, reused & grown across time steps
	PetscInt      bufcap = 0;
	PetscScalar  *bx = NULL, *by = NULL, *bz = NULL, *bT = NULL, *bp = NULL;
	PetscScalar  *bsxx = NULL, *bsyy = NULL, *bszz = NULL;
	PetscScalar  *bsxy = NULL, *bsxz = NULL, *bsyz = NULL;
	PetscScalar  *bj2s = NULL, *bj2e = NULL, *beta = NULL, *baps = NULL;

	// The plugin ABI uses Cint (int32_t), whereas PetscInt is 64-bit in this
	// build (double_real_Int64). Keep dedicated int32 buffers and convert
	// explicitly to/from PetscInt when reading/writing markers; phase IDs are
	// always small, so no overflow check is needed.
	int32_t      *bphase_in = NULL, *bphase_out = NULL;

	// per-time-step cell-centred J2 invariant buffers (indexed like jr->svCell[])
	PetscInt      cellcap = 0;
	PetscScalar  *cellJ2Stress = NULL, *cellJ2StrainRate = NULL;
}
//---------------------------------------------------------------------------
PetscBool PhTrPluginIsActive(void)
{
	return active;
}
//---------------------------------------------------------------------------
// PetscRegisterFinalize callback: releases plugin resources exactly once per
// process, at PetscFinalize() time (i.e. after the last LaMEMLibSolve() call,
// however many times it ran). Julia cannot be re-initialised once torn down,
// and dlclose()-ing a library holding an initialised Julia runtime is not
// supported, so this must never run more than once and must never run from
// LaMEMLibSolve() itself.
static PetscErrorCode PhTrPluginFinalize(void)
{
	PetscFunctionBeginUser;

	if(active && atexitFn)
	{
		// best-effort: some embedded Julia runtimes hang or crash on
		// finalize. See the Phase 1 report for the outcome observed here.
		atexitFn(0);
	}

	if(handle)
	{
		PetscCall(PetscDLClose(&handle));
	}

	handle   = NULL;
	fn       = NULL;
	atexitFn = NULL;
	active   = PETSC_FALSE;

	PetscCall(PetscFree(bx));         PetscCall(PetscFree(by));   PetscCall(PetscFree(bz));
	PetscCall(PetscFree(bT));         PetscCall(PetscFree(bp));
	PetscCall(PetscFree(bsxx));       PetscCall(PetscFree(bsyy)); PetscCall(PetscFree(bszz));
	PetscCall(PetscFree(bsxy));       PetscCall(PetscFree(bsxz)); PetscCall(PetscFree(bsyz));
	PetscCall(PetscFree(bj2s));       PetscCall(PetscFree(bj2e));
	PetscCall(PetscFree(beta));       PetscCall(PetscFree(baps));
	PetscCall(PetscFree(bphase_in));  PetscCall(PetscFree(bphase_out));
	PetscCall(PetscFree(cellJ2Stress));
	PetscCall(PetscFree(cellJ2StrainRate));

	bufcap = 0;
	cellcap = 0;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PhTrPluginLoad(AdvCtx *actx)
{
	char      lib[_str_len_];
	PetscBool found;
	void     *sym;

	PetscFunctionBeginUser;

	(void)actx;

	// Julia cannot be initialised twice in one process (adjoint/inversion
	// drivers call LaMEMLibSolve() repeatedly - see src/adjoint.cpp). Only
	// attempt loading once per process; later calls reuse whatever the first
	// call established (loaded plugin, or none).
	if(initTried) PetscFunctionReturn(0);

	initTried = PETSC_TRUE;

	PetscCall(PetscOptionsGetString(NULL, NULL, "-phase_transition_lib", lib, _str_len_, &found));

	if(!found) PetscFunctionReturn(0); // nothing requested -> built-in transitions only

	// Register the (single, process-wide) teardown callback before doing
	// anything else, so a failure partway through loading still gets cleaned
	// up at PetscFinalize().
	PetscCall(PetscRegisterFinalize(PhTrPluginFinalize));

	// Load the plugin. Use PETSC_DL_NOW (resolve all symbols immediately, so
	// that the Julia runtime bundled with the plugin, and dlopen'd as one of
	// its dependencies, is fully linked in before we look up any symbol).
	// PetscDLOpen(..., PETSC_DL_NOW, ...) maps to dlopen(path, RTLD_NOW |
	// RTLD_GLOBAL) (RTLD_GLOBAL is PETSc's default dlflags2 unless
	// PETSC_DL_LOCAL is requested - see PETSc's src/sys/dll/dlimpl.c).
	PetscCall(PetscDLOpen(lib, PETSC_DL_NOW, &handle));

	// Locate jl_parse_opts. This, and jl_init_with_image_handle below, live
	// in libjulia, which is a *dependency* of the plugin library, not of
	// LaMEM: LaMEM never links against libjulia. Because the juliac-built
	// plugin records libjulia (and its own dependents) as its own linked
	// dependencies, dlopen(RTLD_NOW|RTLD_GLOBAL) on the plugin pulls them
	// into the process, and PetscDLSym (= dlsym) on the plugin's handle can
	// resolve symbols pulled in transitively through it. Verified with a
	// standalone harness that performs exactly this dlopen/dlsym sequence
	// against the actual juliac-built bundle (no julia.h, no -ljulia); see
	// doc/phase_transition_plugin_PHASE1_REPORT.md for exactly what was run.
	PetscCall(PetscDLSym(handle, "jl_parse_opts", &sym));

	if(!sym)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB,
			"phase_transition_lib: could not find jl_parse_opts via the "
			"plugin handle (%s). The plugin must be a juliac-built bundle "
			"with libjulia as a resolvable dependency.", lib);
	}

	// Tell Julia's runtime, BEFORE jl_init_with_image_handle, not to install
	// its own signal handlers for SIGSEGV/SIGBUS/SIGFPE etc. This is the
	// ABI-stable, documented way to do this (jl_parse_opts is an exported
	// Julia C entry point that predates and postdates the specific layout of
	// the jl_options struct, so this does not depend on that struct's
	// layout). A prior version of this file instead called
	// PetscPushSignalHandler(PetscSignalHandlerDefault, NULL) AFTER Julia
	// init, believing that would reclaim PETSc's handlers; this was verified
	// to be a no-op (PETSc's signal.c only calls signal() when its internal
	// SignalSet flag is false, and it has been true since PetscInitialize;
	// the push just duplicates a stack entry without changing what handler
	// is actually installed) and, on macOS, ineffective for a second reason:
	// Julia additionally installs Mach exception ports for SIGSEGV/SIGBUS
	// that a POSIX signal()/sigaction() call cannot displace. Passing
	// --handle-signals=no via jl_parse_opts prevents Julia from installing
	// either mechanism in the first place, so PETSc's handlers (installed
	// earlier, at PetscInitialize) remain the only ones in effect. See
	// doc/phase_transition_plugin_PHASE1_REPORT.md for the verification.
	{
		static char argv0[] = "lamem";
		static char argv1[] = "--handle-signals=no";
		static char *jlargv[2] = { argv0, argv1 };
		char **jlargvp = jlargv;
		int    jlargc  = 2;

		((JlParseOptsFn)sym)(&jlargc, &jlargvp);
	}

	// now locate and call jl_init_with_image_handle itself
	PetscCall(PetscDLSym(handle, "jl_init_with_image_handle", &sym));

	if(!sym)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB,
			"phase_transition_lib: could not find jl_init_with_image_handle "
			"via the plugin handle (%s).", lib);
	}

	((JlInitWithImageHandleFn)sym)((void*)handle);

	// Julia's runtime (specifically its bundled LinearAlgebra/OpenBLAS_jll,
	// which is part of the Julia base image and gets initialised regardless
	// of whether the plugin itself uses linear algebra) calls libblastrampoline's
	// lbt_forward() during jl_init_with_image_handle to register ITS OWN
	// bundled OpenBLAS as the default BLAS/LAPACK target. Because LaMEM and
	// the plugin share the SAME process-wide libblastrampoline instance
	// (confirmed: only one libblastrampoline.5.dylib is ever loaded, by
	// design - it has one install name and dyld/macOS's two-level namespace
	// resolves the plugin's dependency to the same already-loaded image),
	// this OVERWRITES the forwarding table entries that PETSc's own startup
	// (via the LBT_DEFAULT_LIBS environment variable, consulted once by
	// libblastrampoline at its own load time) had set up - even though
	// LBT_DEFAULT_LIBS is still set in the environment, Julia's init does
	// not re-read it and unconditionally registers its own default. Confirmed
	// by observing "Error: no BLAS/LAPACK library loaded for idamax_()" etc.
	// on the very next PETSc linear solve after loading the plugin.
	//
	// Mitigation: call lbt_forward() ourselves, immediately after Julia
	// init, to re-register PETSc's own BLAS libraries (clear=0, i.e.
	// additive: this does not remove Julia's own registrations, it just
	// re-points the specific symbols PETSc needs back to PETSc's chosen
	// library). lbt_forward is resolved via PetscDLSym(NULL, ...), which
	// PETSc's dlimpl.c implements as a dlsym(RTLD_DEFAULT, ...)-equivalent
	// process-wide lookup (a NULL handle is NOT the plugin's own handle:
	// looking it up through the plugin handle specifically did not resolve
	// it in testing, only the process-wide lookup did - both were tried and
	// the outcome recorded in the Phase 1 report).
	{
		typedef int32_t (*LbtForwardFn)(const char*, int32_t, int32_t, const char*);
		void *lbtSym = NULL;

		PetscCall(PetscDLSym(NULL, "lbt_forward", &lbtSym));

		if(lbtSym)
		{
			char      ilp64lib[_str_len_], lp64lib[_str_len_];
			PetscBool foundIlp64, foundLp64;

			PetscCall(PetscOptionsGetString(NULL, NULL, "-phase_transition_lbt_ilp64", ilp64lib, _str_len_, &foundIlp64));
			PetscCall(PetscOptionsGetString(NULL, NULL, "-phase_transition_lbt_lp64",  lp64lib,  _str_len_, &foundLp64));

			if(foundIlp64) ((LbtForwardFn)lbtSym)(ilp64lib, 0, 0, NULL);
			if(foundLp64)  ((LbtForwardFn)lbtSym)(lp64lib,  0, 0, NULL);

			if(foundIlp64 || foundLp64)
			{
				PetscPrintf(PETSC_COMM_WORLD, "Phase transition plugin  : re-forwarded PETSc's BLAS/LAPACK "
					"via libblastrampoline after Julia init\n");
			}
			else
			{
				PetscPrintf(PETSC_COMM_WORLD, "Phase transition plugin  : WARNING: Julia's runtime may have "
					"overridden PETSc's BLAS/LAPACK forwarding (libblastrampoline). Pass "
					"-phase_transition_lbt_ilp64 <path> -phase_transition_lbt_lp64 <path> "
					"(the same two libraries as LBT_DEFAULT_LIBS) to restore it.\n");
			}
		}
	}

	// look up the actual phase-transition entry point
	PetscCall(PetscDLSym(handle, "lamem_phase_transition", &sym));

	if(!sym)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB,
			"phase_transition_lib: symbol lamem_phase_transition not found in %s", lib);
	}

	fn = (PhTrPluginFn)sym;

	// jl_atexit_hook is optional; only used (best-effort) in PhTrPluginFinalize
	PetscCall(PetscDLSym(handle, "jl_atexit_hook", &sym));
	atexitFn = (JlAtexitHookFn)sym; // NULL if not found - fine

	active = PETSC_TRUE;

	PetscPrintf(PETSC_COMM_WORLD, "Phase transition plugin  : %s\n", lib);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// (re)allocate the per-marker SoA scratch buffers if the marker count grows
static PetscErrorCode PhTrPluginEnsureMarkerCapacity(PetscInt n)
{
	PetscFunctionBeginUser;

	if(n <= bufcap) PetscFunctionReturn(0);

	PetscCall(PetscFree(bx));         PetscCall(PetscFree(by));   PetscCall(PetscFree(bz));
	PetscCall(PetscFree(bT));         PetscCall(PetscFree(bp));
	PetscCall(PetscFree(bsxx));       PetscCall(PetscFree(bsyy)); PetscCall(PetscFree(bszz));
	PetscCall(PetscFree(bsxy));       PetscCall(PetscFree(bsxz)); PetscCall(PetscFree(bsyz));
	PetscCall(PetscFree(bj2s));       PetscCall(PetscFree(bj2e));
	PetscCall(PetscFree(beta));       PetscCall(PetscFree(baps));
	PetscCall(PetscFree(bphase_in));  PetscCall(PetscFree(bphase_out));

	PetscCall(PetscMalloc((size_t)n*sizeof(PetscScalar), &bx));
	PetscCall(PetscMalloc((size_t)n*sizeof(PetscScalar), &by));
	PetscCall(PetscMalloc((size_t)n*sizeof(PetscScalar), &bz));
	PetscCall(PetscMalloc((size_t)n*sizeof(PetscScalar), &bT));
	PetscCall(PetscMalloc((size_t)n*sizeof(PetscScalar), &bp));
	PetscCall(PetscMalloc((size_t)n*sizeof(PetscScalar), &bsxx));
	PetscCall(PetscMalloc((size_t)n*sizeof(PetscScalar), &bsyy));
	PetscCall(PetscMalloc((size_t)n*sizeof(PetscScalar), &bszz));
	PetscCall(PetscMalloc((size_t)n*sizeof(PetscScalar), &bsxy));
	PetscCall(PetscMalloc((size_t)n*sizeof(PetscScalar), &bsxz));
	PetscCall(PetscMalloc((size_t)n*sizeof(PetscScalar), &bsyz));
	PetscCall(PetscMalloc((size_t)n*sizeof(PetscScalar), &bj2s));
	PetscCall(PetscMalloc((size_t)n*sizeof(PetscScalar), &bj2e));
	PetscCall(PetscMalloc((size_t)n*sizeof(PetscScalar), &beta));
	PetscCall(PetscMalloc((size_t)n*sizeof(PetscScalar), &baps));
	PetscCall(PetscMalloc((size_t)n*sizeof(int32_t),     &bphase_in));
	PetscCall(PetscMalloc((size_t)n*sizeof(int32_t),     &bphase_out));

	bufcap = n;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// (re)allocate the per-cell J2 buffers if the local cell count grows
// (it is constant for a fixed grid/partition, but Phase 1 keeps this
// defensive in case that assumption ever changes)
static PetscErrorCode PhTrPluginEnsureCellCapacity(PetscInt ncells)
{
	PetscFunctionBeginUser;

	if(ncells <= cellcap) PetscFunctionReturn(0);

	PetscCall(PetscFree(cellJ2Stress));
	PetscCall(PetscFree(cellJ2StrainRate));

	PetscCall(PetscMalloc((size_t)ncells*sizeof(PetscScalar), &cellJ2Stress));
	PetscCall(PetscMalloc((size_t)ncells*sizeof(PetscScalar), &cellJ2StrainRate));

	cellcap = ncells;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// Compute a genuinely CELL-CENTRED second invariant (J2) of the deviatoric
// stress and of the deviatoric strain rate, for every local cell, once per
// time step. This is analogous to LaMEM's own ParaView j2_dev_stress /
// j2_strain_rate fields (src/outFunct.cpp: PVOutWriteJ2DevStress /
// PVOutWriteJ2StrainRate) in that it uses the same underlying quantities
// (svCell diagonal components + the surrounding svXYEdge/svXZEdge/svYZEdge
// off-diagonal components), but it is NOT identical: ParaView's routines
// interpolate the (already squared) per-edge-node J2 contributions onto grid
// CORNERS (via InterpXYEdgeCorner etc.), i.e. they produce a node-centred
// field. This routine instead produces a genuinely CELL-centred field, by
// averaging, for each cell, the 4 surrounding edge values of each
// off-diagonal component (the same pattern src/JacResAux.cpp's
// JacResGetSHmax/JacResGetEHmax use to build a cell-centred sxy from the
// XY-edge grid) before squaring, then combining with the cell's own diagonal
// components -- i.e. it follows LaMEM's own PVOutWriteDevStress/
// PVOutWriteStrainRate cell-vs-edge geometry (see the comment above
// GET_SXY/GET_SXY etc. there), not the corner-interpolation of
// PVOutWriteJ2DevStress/PVOutWriteJ2StrainRate.
//
// Geometry (confirmed against src/JacRes.cpp's JacResGetEffStrainRate, which
// fills svXYEdge/svXZEdge/svYZEdge from velocity gradients, and
// src/fdstag.cpp's FDSTAGCreateDMDA, which sizes DA_XY/DA_XZ/DA_YZ):
//   - DA_XY has the full (Nx,Ny) node count and Nz-1 (cell count) in Z: an
//     XY-edge varies over (i,j) at fixed k, so cell (i,j,k)'s 4 surrounding
//     XY-edges are XY(i,j,k), XY(i+1,j,k), XY(i,j+1,k), XY(i+1,j+1,k).
//   - DA_XZ has the full (Nx,Nz) node count and Ny-1 in Y: cell (i,j,k)'s 4
//     surrounding XZ-edges are XZ(i,j,k), XZ(i+1,j,k), XZ(i,j,k+1), XZ(i+1,j,k+1).
//   - DA_YZ has the full (Ny,Nz) node count and Nx-1 in X: cell (i,j,k)'s 4
//     surrounding YZ-edges are YZ(i,j,k), YZ(i,j+1,k), YZ(i,j,k+1), YZ(i,j+1,k+1).
static PetscErrorCode PhTrPluginComputeJ2(AdvCtx *actx)
{
	FDSTAG      *fs;
	JacRes      *jr;
	Scaling     *scal;
	SolVarCell  *svCell;
	Vec          lxy_s, lxz_s, lyz_s;  // local edge vectors: (effective) stress
	Vec          lxy_d, lxz_d, lyz_d;  // local edge vectors: strain rate
	PetscScalar ***axy_s, ***axz_s, ***ayz_s;
	PetscScalar ***axy_d, ***axz_d, ***ayz_d;
	PetscInt     i, j, k, sx, sy, sz, nx, ny, nz, iter, ncells;
	PetscScalar  pf;

	PetscFunctionBeginUser;

	jr   = actx->jr;
	fs   = actx->fs;
	scal = jr->scal;

	ncells = fs->nCells;
	PetscCall(PhTrPluginEnsureCellCapacity(ncells));

	// stabilized-stress prefactor, exactly as PVOutWriteDevStress/
	// PVOutWriteJ2DevStress use it
	if(jr->ctrl.initGuess) pf = 0.0;
	else                   pf = 2.0;

	// --- fill local edge vectors from svXYEdge/svXZEdge/svYZEdge, and ghost-exchange them ---
	PetscCall(FDSTAGGetLocalVectorEdge(fs, &lxy_s, &lxz_s, &lyz_s));
	PetscCall(FDSTAGGetLocalVectorEdge(fs, &lxy_d, &lxz_d, &lyz_d));

	PetscCall(DMDAVecGetArray(fs->DA_XY, lxy_s, &axy_s));
	PetscCall(DMDAVecGetArray(fs->DA_XY, lxy_d, &axy_d));
	PetscCall(DMDAGetCorners(fs->DA_XY, &sx, &sy, &sz, &nx, &ny, &nz));
	iter = 0;
	START_STD_LOOP
	{
		SolVarEdge *svEdge = &jr->svXYEdge[iter++];
		axy_s[k][j][i] = svEdge->s + pf*svEdge->svDev.eta_st*svEdge->d;
		axy_d[k][j][i] = svEdge->d;
	}
	END_STD_LOOP
	PetscCall(DMDAVecRestoreArray(fs->DA_XY, lxy_s, &axy_s));
	PetscCall(DMDAVecRestoreArray(fs->DA_XY, lxy_d, &axy_d));
	LOCAL_TO_LOCAL(fs->DA_XY, lxy_s);
	LOCAL_TO_LOCAL(fs->DA_XY, lxy_d);

	PetscCall(DMDAVecGetArray(fs->DA_XZ, lxz_s, &axz_s));
	PetscCall(DMDAVecGetArray(fs->DA_XZ, lxz_d, &axz_d));
	PetscCall(DMDAGetCorners(fs->DA_XZ, &sx, &sy, &sz, &nx, &ny, &nz));
	iter = 0;
	START_STD_LOOP
	{
		SolVarEdge *svEdge = &jr->svXZEdge[iter++];
		axz_s[k][j][i] = svEdge->s + pf*svEdge->svDev.eta_st*svEdge->d;
		axz_d[k][j][i] = svEdge->d;
	}
	END_STD_LOOP
	PetscCall(DMDAVecRestoreArray(fs->DA_XZ, lxz_s, &axz_s));
	PetscCall(DMDAVecRestoreArray(fs->DA_XZ, lxz_d, &axz_d));
	LOCAL_TO_LOCAL(fs->DA_XZ, lxz_s);
	LOCAL_TO_LOCAL(fs->DA_XZ, lxz_d);

	PetscCall(DMDAVecGetArray(fs->DA_YZ, lyz_s, &ayz_s));
	PetscCall(DMDAVecGetArray(fs->DA_YZ, lyz_d, &ayz_d));
	PetscCall(DMDAGetCorners(fs->DA_YZ, &sx, &sy, &sz, &nx, &ny, &nz));
	iter = 0;
	START_STD_LOOP
	{
		SolVarEdge *svEdge = &jr->svYZEdge[iter++];
		ayz_s[k][j][i] = svEdge->s + pf*svEdge->svDev.eta_st*svEdge->d;
		ayz_d[k][j][i] = svEdge->d;
	}
	END_STD_LOOP
	PetscCall(DMDAVecRestoreArray(fs->DA_YZ, lyz_s, &ayz_s));
	PetscCall(DMDAVecRestoreArray(fs->DA_YZ, lyz_d, &ayz_d));
	LOCAL_TO_LOCAL(fs->DA_YZ, lyz_s);
	LOCAL_TO_LOCAL(fs->DA_YZ, lyz_d);

	// --- re-read (ghosted) edge arrays and average onto cells ---
	PetscCall(DMDAVecGetArray(fs->DA_XY, lxy_s, &axy_s));
	PetscCall(DMDAVecGetArray(fs->DA_XY, lxy_d, &axy_d));
	PetscCall(DMDAVecGetArray(fs->DA_XZ, lxz_s, &axz_s));
	PetscCall(DMDAVecGetArray(fs->DA_XZ, lxz_d, &axz_d));
	PetscCall(DMDAVecGetArray(fs->DA_YZ, lyz_s, &ayz_s));
	PetscCall(DMDAVecGetArray(fs->DA_YZ, lyz_d, &ayz_d));

	PetscCall(DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz));
	iter = 0;
	START_STD_LOOP
	{
		svCell = &jr->svCell[iter];

		// average the 4 surrounding edge values per off-diagonal direction
		// onto this cell (same pattern as JacResGetSHmax/JacResGetEHmax)
		PetscScalar sxy = (axy_s[k][j][i] + axy_s[k][j][i+1] + axy_s[k][j+1][i] + axy_s[k][j+1][i+1])/4.0;
		PetscScalar sxz = (axz_s[k][j][i] + axz_s[k][j][i+1] + axz_s[k+1][j][i] + axz_s[k+1][j][i+1])/4.0;
		PetscScalar syz = (ayz_s[k][j][i] + ayz_s[k][j+1][i] + ayz_s[k+1][j][i] + ayz_s[k+1][j+1][i])/4.0;

		PetscScalar dxy = (axy_d[k][j][i] + axy_d[k][j][i+1] + axy_d[k][j+1][i] + axy_d[k][j+1][i+1])/4.0;
		PetscScalar dxz = (axz_d[k][j][i] + axz_d[k][j][i+1] + axz_d[k+1][j][i] + axz_d[k+1][j][i+1])/4.0;
		PetscScalar dyz = (ayz_d[k][j][i] + ayz_d[k][j+1][i] + ayz_d[k+1][j][i] + ayz_d[k+1][j+1][i])/4.0;

		PetscScalar sxx = svCell->sxx + pf*svCell->svDev.eta_st*svCell->dxx;
		PetscScalar syy = svCell->syy + pf*svCell->svDev.eta_st*svCell->dyy;
		PetscScalar szz = svCell->szz + pf*svCell->svDev.eta_st*svCell->dzz;

		PetscScalar J2s = 0.5*(sxx*sxx + syy*syy + szz*szz) + sxy*sxy + sxz*sxz + syz*syz;
		PetscScalar J2e = 0.5*(svCell->dxx*svCell->dxx + svCell->dyy*svCell->dyy + svCell->dzz*svCell->dzz)
		                   + dxy*dxy + dxz*dxz + dyz*dyz;

		cellJ2Stress[iter]     = sqrt(J2s)*scal->stress;
		cellJ2StrainRate[iter] = sqrt(J2e)*scal->strain_rate;

		iter++;
	}
	END_STD_LOOP

	PetscCall(DMDAVecRestoreArray(fs->DA_XY, lxy_s, &axy_s));
	PetscCall(DMDAVecRestoreArray(fs->DA_XY, lxy_d, &axy_d));
	PetscCall(DMDAVecRestoreArray(fs->DA_XZ, lxz_s, &axz_s));
	PetscCall(DMDAVecRestoreArray(fs->DA_XZ, lxz_d, &axz_d));
	PetscCall(DMDAVecRestoreArray(fs->DA_YZ, lyz_s, &ayz_s));
	PetscCall(DMDAVecRestoreArray(fs->DA_YZ, lyz_d, &ayz_d));

	PetscCall(FDSTAGRestoreLocalVectorEdge(fs, &lxy_s, &lxz_s, &lyz_s));
	PetscCall(FDSTAGRestoreLocalVectorEdge(fs, &lxy_d, &lxz_d, &lyz_d));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PhTrPluginApply(AdvCtx *actx)
{
	JacRes      *jr;
	Scaling     *scal;
	Marker      *P;
	PetscInt     i, ID, n, changed_loc, changed_glob;
	PetscInt     numPhases;
	PetscScalar  time_dim, pShift;
	int          rc;

	PetscFunctionBeginUser;

	if(!active) PetscFunctionReturn(0); // no plugin loaded -> nothing to do

	jr        = actx->jr;
	scal      = jr->scal;
	n         = actx->nummark;
	numPhases = actx->dbm->numPhases;

	PetscCall(PhTrPluginEnsureMarkerCapacity(n));

	// once-per-step, genuinely cell-centred J2 invariants (stress & strain rate)
	PetscCall(PhTrPluginComputeJ2(actx));

	// dimensional simulation time
	time_dim = jr->bc->ts->time*scal->time;

	// pressure shift used by the rheology/plasticity and by the built-in
	// "Pressure" Constant transition (see Check_Constant_Phase_Transition)
	pShift = (jr->ctrl.pShift != 0.0) ? jr->ctrl.pShift : 0.0;

	// build SoA input arrays (dimensional) from the local markers
	for(i = 0; i < n; i++)
	{
		P  = &actx->markers[i];
		ID = actx->cellnum[i];

		bx[i] = P->X[0]*scal->length;
		by[i] = P->X[1]*scal->length;
		bz[i] = P->X[2]*scal->length;
		bT[i] = P->T*scal->temperature - scal->Tshift;
		bp[i] = (P->p + pShift)*scal->stress;

		bsxx[i] = P->S.xx*scal->stress;
		bsyy[i] = P->S.yy*scal->stress;
		bszz[i] = P->S.zz*scal->stress;
		bsxy[i] = P->S.xy*scal->stress;
		bsxz[i] = P->S.xz*scal->stress;
		bsyz[i] = P->S.yz*scal->stress;

		bj2s[i] = cellJ2Stress[ID];
		bj2e[i] = cellJ2StrainRate[ID];

		beta[i] = jr->svCell[ID].svDev.eta*scal->viscosity;
		baps[i] = jr->svCell[ID].svDev.APS; // dimensionless

		bphase_in[i]  = (int32_t)P->phase;
		bphase_out[i] = (int32_t)P->phase;
	}

	// call the plugin once for all local markers on this rank. A negative
	// return value is a plugin-signalled failure (e.g. an exception caught
	// on the Julia side) - see phase_transition_plugin.h.
	rc = fn((size_t)n,
		bx, by, bz, bT, bp, (double)time_dim,
		bsxx, bsyy, bszz, bsxy, bsxz, bsyz,
		bj2s, bj2e, beta, baps,
		bphase_in, bphase_out);

	if(rc < 0)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB,
			"phase_transition_lib: plugin reported failure (lamem_phase_transition returned %d)", rc);
	}

	// write back only markers whose phase actually changed, validating each
	// new phase against the material database bounds (an out-of-range phase
	// from the plugin would otherwise corrupt svCell->phRat[] in
	// ADVInterpMarkToCell via an out-of-bounds array write)
	changed_loc = 0;

	for(i = 0; i < n; i++)
	{
		if(bphase_out[i] != bphase_in[i])
		{
			if(bphase_out[i] < 0 || bphase_out[i] >= numPhases)
			{
				SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER,
					"phase_transition_lib: plugin returned out-of-range phase %d for local marker %" PetscInt_FMT
					" (valid range: 0..%" PetscInt_FMT ")", bphase_out[i], i, numPhases-1);
			}

			P = &actx->markers[i];
			P->phase = (PetscInt)bphase_out[i];
			changed_loc++;
		}
	}

	PetscCallMPI(MPI_Allreduce(&changed_loc, &changed_glob, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD));

	if(changed_glob)
	{
		// re-validate all local marker phases (cheap, collective-safe;
		// matches the pattern ADVRemap uses before ADVInterpMarkToCell) and
		// update cell phase ratios to reflect the plugin's changes
		PetscCall(ADVCheckMarkPhases(actx));
		PetscCall(ADVInterpMarkToCell(actx));
	}

	PetscPrintf(PETSC_COMM_WORLD, "Phase transition plugin  : %" PetscInt_FMT " marker(s) changed phase\n", changed_glob);

	PetscFunctionReturn(0);
}

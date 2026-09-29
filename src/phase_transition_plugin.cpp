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
#include "advect.h"
#include "scaling.h"
#include "JacRes.h"
#include "bc.h"
#include "tssolve.h"
#include "phase_transition_plugin.h"
#include <cstddef>
//---------------------------------------------------------------------------
// C ABI of the plugin function (see phase_transition_plugin.h for the
// full, documented signature)
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

// function pointer signature of jl_init_with_image_handle(void *handle)
typedef void (*JlInitWithImageHandleFn)(void *handle);

// function pointer signature of jl_atexit_hook(int status)
typedef void (*JlAtexitHookFn)(int status);

//---------------------------------------------------------------------------
// module-local state (Phase 1: single global plugin instance, one per rank)
namespace
{
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
	PetscInt     *bphase_in = NULL, *bphase_out = NULL;
}
//---------------------------------------------------------------------------
PetscBool PhTrPluginIsActive(void)
{
	return active;
}
//---------------------------------------------------------------------------
PetscErrorCode PhTrPluginLoad(AdvCtx *actx)
{
	char      lib[_str_len_];
	PetscBool found;
	void     *sym;

	PetscFunctionBeginUser;

	active = PETSC_FALSE;
	handle = NULL;
	fn     = NULL;

	PetscCall(PetscOptionsGetString(NULL, NULL, "-phase_transition_lib", lib, _str_len_, &found));

	if(!found) PetscFunctionReturn(0); // nothing requested -> built-in transitions only

	(void)actx;

	// Load the plugin. Use PETSC_DL_NOW (resolve all symbols immediately, so
	// that the Julia runtime bundled with the plugin, and dlopen'd as one of
	// its dependencies, is fully linked in before we look up any symbol).
	PetscCall(PetscDLOpen(lib, PETSC_DL_NOW, &handle));

	// Locate jl_init_with_image_handle. This symbol lives in libjulia, which
	// is a *dependency* of the plugin library, not of LaMEM: LaMEM never
	// links against libjulia. On macOS/Linux, PetscDLOpen(..., PETSC_DL_NOW)
	// resolves the plugin's dependencies (libjulia, its stdlib images, etc.)
	// into the process at load time, and PetscDLSym on the plugin's own
	// handle can see symbols pulled in transitively through it, because the
	// underlying dlopen(path, RTLD_NOW) call resolves symbols against the
	// full set of libraries the plugin depends on. If this lookup fails on a
	// platform where that is not the case, the fallback is to PetscDLOpen the
	// libjulia.dylib/.so found next to the plugin explicitly and look the
	// symbol up there instead.
	PetscCall(PetscDLSym(handle, "jl_init_with_image_handle", &sym));

	if(!sym)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB,
			"phase_transition_lib: could not find jl_init_with_image_handle "
			"via the plugin handle (%s). The plugin must be a juliac-built "
			"bundle with libjulia as a resolvable dependency.", lib);
	}

	((JlInitWithImageHandleFn)sym)((void*)handle);

	// Julia's runtime installs its own signal handlers during jl_init_*.
	// LaMEM (via PETSc/MPI) needs to keep catching SIGSEGV/SIGFPE/SIGBUS
	// etc. through PETSc's own handler for its error reporting to work.
	// Re-install PETSc's default signal handler now that Julia is up.
	PetscCall(PetscPushSignalHandler(PetscSignalHandlerDefault, NULL));

	// look up the actual phase-transition entry point
	PetscCall(PetscDLSym(handle, "lamem_phase_transition", &sym));

	if(!sym)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB,
			"phase_transition_lib: symbol lamem_phase_transition not found in %s", lib);
	}

	fn = (PhTrPluginFn)sym;

	// jl_atexit_hook is optional; only used (best-effort) in PhTrPluginDestroy
	PetscCall(PetscDLSym(handle, "jl_atexit_hook", &sym));
	atexitFn = (JlAtexitHookFn)sym; // NULL if not found - fine

	active = PETSC_TRUE;

	PetscPrintf(PETSC_COMM_WORLD, "Phase transition plugin  : %s\n", lib);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// (re)allocate the SoA scratch buffers if the marker count grows
static PetscErrorCode PhTrPluginEnsureCapacity(PetscInt n)
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
	PetscCall(PetscMalloc((size_t)n*sizeof(PetscInt),    &bphase_in));
	PetscCall(PetscMalloc((size_t)n*sizeof(PetscInt),    &bphase_out));

	bufcap = n;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PhTrPluginApply(AdvCtx *actx)
{
	JacRes      *jr;
	Scaling     *scal;
	Marker      *P;
	SolVarCell  *svCell;
	PetscInt     i, ID, n, changed_loc, changed_glob;
	PetscScalar  time_dim;

	PetscFunctionBeginUser;

	if(!active) PetscFunctionReturn(0); // no plugin loaded -> nothing to do

	jr   = actx->jr;
	scal = jr->scal;
	n    = actx->nummark;

	PetscCall(PhTrPluginEnsureCapacity(n));

	// dimensional simulation time
	time_dim = jr->bc->ts->time*scal->time;

	// build SoA input arrays (dimensional) from the local markers
	// NOTE (Phase 1 scope): the cell-centred J2 invariants below use only the
	// diagonal (cell-centered) deviatoric stress/strain-rate components
	// (svCell->sxx,syy,szz / dxx,dyy,dzz). The off-diagonal components live on
	// separate edge-centred DMDA grids (svXYEdge/svXZEdge/svYZEdge) and
	// combining them into a true cell-centred J2 the way ParaView output does
	// (src/outFunct.cpp: PVOutWriteJ2DevStress/PVOutWriteJ2StrainRate) needs a
	// corner-interpolation pass across ghost points. That is deliberately left
	// out of Phase 1 to keep this hook cheap and allocation-free per step; see
	// the Phase 1 report for the exact simplification and its effect on the
	// reported J2 values (diagonal-only J2 is a lower bound on the true J2).
	for(i = 0; i < n; i++)
	{
		P  = &actx->markers[i];
		ID = actx->cellnum[i];
		svCell = &jr->svCell[ID];

		bx[i] = P->X[0]*scal->length;
		by[i] = P->X[1]*scal->length;
		bz[i] = P->X[2]*scal->length;
		bT[i] = P->T*scal->temperature - scal->Tshift;
		bp[i] = P->p*scal->stress;

		bsxx[i] = P->S.xx*scal->stress;
		bsyy[i] = P->S.yy*scal->stress;
		bszz[i] = P->S.zz*scal->stress;
		bsxy[i] = P->S.xy*scal->stress;
		bsxz[i] = P->S.xz*scal->stress;
		bsyz[i] = P->S.yz*scal->stress;

		// cell-centred (diagonal-only) J2 invariants of stress and strain rate
		{
			PetscScalar sxx = svCell->sxx, syy = svCell->syy, szz = svCell->szz;
			PetscScalar dxx = svCell->dxx, dyy = svCell->dyy, dzz = svCell->dzz;

			bj2s[i] = sqrt(0.5*(sxx*sxx + syy*syy + szz*szz))*scal->stress;
			bj2e[i] = sqrt(0.5*(dxx*dxx + dyy*dyy + dzz*dzz))*scal->strain_rate;
		}

		beta[i] = svCell->svDev.eta*scal->viscosity;
		baps[i] = svCell->svDev.APS; // dimensionless

		bphase_in[i]  = (PetscInt)P->phase;
		bphase_out[i] = (PetscInt)P->phase;
	}

	// call the plugin once for all local markers on this rank
	changed_loc = fn((size_t)n,
		bx, by, bz, bT, bp, (double)time_dim,
		bsxx, bsyy, bszz, bsxy, bsxz, bsyz,
		bj2s, bj2e, beta, baps,
		bphase_in, bphase_out);

	(void)changed_loc; // recompute below from the actual diffs (authoritative)

	// write back only markers whose phase actually changed
	changed_loc = 0;

	for(i = 0; i < n; i++)
	{
		if(bphase_out[i] != bphase_in[i])
		{
			P = &actx->markers[i];
			P->phase = bphase_out[i];
			changed_loc++;
		}
	}

	PetscCallMPI(MPI_Allreduce(&changed_loc, &changed_glob, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD));

	if(changed_glob)
	{
		PetscCall(ADVInterpMarkToCell(actx));
	}

	PetscPrintf(PETSC_COMM_WORLD, "Phase transition plugin  : %" PetscInt_FMT " marker(s) changed phase\n", changed_glob);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PhTrPluginDestroy(void)
{
	PetscFunctionBeginUser;

	if(active && atexitFn)
	{
		// best-effort: some Julia runtimes hang or crash on finalize when
		// embedded. See the Phase 1 report for the outcome observed here.
		atexitFn(0);
	}

	if(handle)
	{
		PetscCall(PetscDLClose(&handle));
	}

	handle = NULL;
	fn     = NULL;
	atexitFn = NULL;
	active = PETSC_FALSE;

	PetscCall(PetscFree(bx));         PetscCall(PetscFree(by));   PetscCall(PetscFree(bz));
	PetscCall(PetscFree(bT));         PetscCall(PetscFree(bp));
	PetscCall(PetscFree(bsxx));       PetscCall(PetscFree(bsyy)); PetscCall(PetscFree(bszz));
	PetscCall(PetscFree(bsxy));       PetscCall(PetscFree(bsxz)); PetscCall(PetscFree(bsyz));
	PetscCall(PetscFree(bj2s));       PetscCall(PetscFree(bj2e));
	PetscCall(PetscFree(beta));       PetscCall(PetscFree(baps));
	PetscCall(PetscFree(bphase_in));  PetscCall(PetscFree(bphase_out));

	bufcap = 0;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------

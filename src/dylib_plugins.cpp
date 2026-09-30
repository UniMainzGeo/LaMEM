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
//..............   USER-DEFINED DYLIB PLUGINS (Julia)   .....................
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
#include "parsing.h"
#include "dylib_plugins.h"
#include <cstddef>
#include <cstdint>
#include <cmath>
#if defined(PETSC_HAVE_DLADDR)
#include <dlfcn.h>
#endif
//---------------------------------------------------------------------------
typedef int      (*DylibPluginAbiVersionFn)(void);
typedef void      (*JlParseOptsFn)(int *argc, char ***argv);
typedef void      (*JlInitWithImageHandleFn)(void *handle);
typedef void      (*JlAtexitHookFn)(int status);

// libblastrampoline ABI (not linked; resolved via PetscDLSym), layout from
// libblastrampoline.h, version 5.15.0 (the deployed version)
struct LbtLibraryInfo { char *libname; void *dlhandle; const char *suffix; uint8_t *active_forwards; int32_t interface, complex_retstyle, f2c, cblas; };
struct LbtConfig      { LbtLibraryInfo **loaded_libs; uint32_t build_flags; const char **exported_symbols; uint32_t num_exported_symbols; };
typedef const LbtConfig* (*LbtGetConfigFn)(void);
typedef int32_t          (*LbtForwardFn)(const char*, int32_t, int32_t, const char*);
typedef int32_t          (*LbtGetNumThreadsFn)(void);
typedef void              (*LbtSetNumThreadsFn)(int32_t);

struct LbtSnapshot { char *libname, *suffix; }; // PetscStrallocpy'd, freed after restore
//---------------------------------------------------------------------------
namespace
{
PetscBool     initTried = PETSC_FALSE, active = PETSC_FALSE;
PetscDLHandle handle    = NULL;
DylibPluginFn  fn        = NULL;
JlAtexitHookFn atexitFn = NULL;
char          loadedPath[_str_len_] = "";

PetscInt      bufcap = 0;
PetscScalar  *bx = NULL, *by = NULL, *bz = NULL, *bT = NULL, *bp = NULL, *bT_out = NULL;
PetscScalar  *bsxx = NULL, *bsyy = NULL, *bszz = NULL, *bsxy = NULL, *bsxz = NULL, *bsyz = NULL;
PetscScalar  *bj2s = NULL, *bj2e = NULL, *beta = NULL, *baps = NULL;
int32_t      *bphase_in = NULL, *bphase_out = NULL;

PetscScalar **markerBufs[] = { &bx, &by, &bz, &bT, &bp, &bT_out,
                               &bsxx, &bsyy, &bszz, &bsxy, &bsxz, &bsyz, &bj2s, &bj2e, &beta, &baps
                             };
const int nMarkerBufs = sizeof(markerBufs)/sizeof(markerBufs[0]);

PetscInt      cellcap = 0;
PetscScalar  *cellJ2Stress = NULL, *cellJ2StrainRate = NULL;

enum { LBT_MAX_SNAPSHOT = 16 };
LbtSnapshot lbtSnap[LBT_MAX_SNAPSHOT];
int         lbtSnapCount = 0;
}
//---------------------------------------------------------------------------
PetscBool DylibPluginIsActive(void) { return active; }
//---------------------------------------------------------------------------
// phase_transitions = dylib requires a loaded plugin (lamem_phase_transition
// is mandatory in any library DylibPluginLoad accepts) - checked once at
// startup so a missing/failed dylib_plugin fails immediately, not silently
// at the first step.
PetscErrorCode DylibPluginCheckPhaseTr(DBMat *dbm)
{
	PetscFunctionBeginUser;

	if(dbm->dylibPhaseTr && !active)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER,
		        "phase_transitions = dylib requires a loaded plugin (dylib_plugin = <path> in the .dat file, or -dylib_plugin on the command line)");
	}
	if(active && !dbm->dylibPhaseTr) PetscPrintf(PETSC_COMM_WORLD, "Dylib plugin  : loaded, but phase_transitions = builtin; lamem_phase_transition will not be called\n");

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
static PetscErrorCode DylibPluginFreeBuffers(void)
{
	int i;

	PetscFunctionBeginUser;

	for(i = 0; i < nMarkerBufs; i++) PetscCall(PetscFree(*markerBufs[i]));
	PetscCall(PetscFree(bphase_in)); PetscCall(PetscFree(bphase_out));
	PetscCall(PetscFree(cellJ2Stress)); PetscCall(PetscFree(cellJ2StrainRate));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// dlclose() of a library holding an initialised Julia runtime is unsupported;
// the handle is leaked for the process lifetime, only jl_atexit_hook runs
static PetscErrorCode DylibPluginFinalize(void)
{
	PetscFunctionBeginUser;

	if(active && atexitFn) atexitFn(0);

	fn = NULL; atexitFn = NULL; active = PETSC_FALSE;

	PetscCall(DylibPluginFreeBuffers());

	bufcap = 0; cellcap = 0;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// Snapshot libblastrampoline's forwarding table before jl_init_with_image_handle
// runs: Julia's LinearAlgebra.__init__ calls lbt_forward(clear=1), which both
// overwrites the table and frees the strings behind any earlier snapshot, so
// they must be copied out first. lbt_get_config's struct layout is not
// versioned; trust it only if the resolved symbol lives in a
// "libblastrampoline.5" (macOS) or "libblastrampoline.so.5" (Linux) image
// (dladdr), else skip silently.
static PetscErrorCode DylibPluginSnapshotLbt(void)
{
	void *sym = NULL;

	PetscFunctionBeginUser;

	lbtSnapCount = 0;

	PetscCall(PetscDLSym(NULL, "lbt_get_config", &sym));
	if(!sym) PetscFunctionReturn(0);

#if defined(PETSC_HAVE_DLADDR)
	Dl_info info;
	if(!dladdr(sym, &info) || !info.dli_fname ||
	   (!strstr(info.dli_fname, "libblastrampoline.5") && !strstr(info.dli_fname, "libblastrampoline.so.5")))
	{
		PetscFunctionReturn(0);
	}
#else
	PetscFunctionReturn(0);
#endif

	const LbtConfig *cfg = ((LbtGetConfigFn)sym)();
	int              i;

	if(!cfg || !cfg->loaded_libs) PetscFunctionReturn(0);

	for(i = 0; i < LBT_MAX_SNAPSHOT && cfg->loaded_libs[i] != NULL; i++)
	{
		PetscCall(PetscStrallocpy(cfg->loaded_libs[i]->libname, &lbtSnap[i].libname));
		PetscCall(PetscStrallocpy(cfg->loaded_libs[i]->suffix,  &lbtSnap[i].suffix));
	}

	lbtSnapCount = i;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// re-forward the snapshotted BLAS libs (clear=0, additive) and restore the
// thread count, both clobbered by Julia's LinearAlgebra.__init__
static PetscErrorCode DylibPluginRestoreLbt(int32_t nthreadsBefore)
{
	void *fwdSym = NULL, *getSym = NULL, *setSym = NULL;
	int   i;

	PetscFunctionBeginUser;

	PetscCall(PetscDLSym(NULL, "lbt_forward", &fwdSym));

	if(fwdSym && lbtSnapCount > 0)
	{
		for(i = 0; i < lbtSnapCount; i++)
		{
			int32_t rc = ((LbtForwardFn)fwdSym)(lbtSnap[i].libname, 0, 0, lbtSnap[i].suffix);

			if(rc <= 0)
			{
				SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB,
				        "dylib_plugin: failed to re-forward BLAS/LAPACK library '%s' after Julia init",
				        lbtSnap[i].libname);
			}

			PetscCall(PetscFree(lbtSnap[i].libname));
			PetscCall(PetscFree(lbtSnap[i].suffix));
		}

		lbtSnapCount = 0;
	}
	else if(!fwdSym)
	{
		PetscPrintf(PETSC_COMM_WORLD, "Dylib plugin  : no libblastrampoline found; PETSc's BLAS/LAPACK forwarding was not touched\n");
	}

	PetscCall(PetscDLSym(NULL, "lbt_get_num_threads", &getSym));
	PetscCall(PetscDLSym(NULL, "lbt_set_num_threads", &setSym));
	if(getSym && setSym && nthreadsBefore > 0) ((LbtSetNumThreadsFn)setSym)(nthreadsBefore);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode DylibPluginLoad(AdvCtx *actx, FB *fb)
{
	char      lib[_str_len_];
	void     *sym;
	int32_t   nthreadsBefore = -1;

	PetscFunctionBeginUser;

	(void)actx;

	// dylib_plugin = <path> in the .dat file; -dylib_plugin <path> on the
	// command line takes precedence (same override rule getStringParam
	// uses for every other top-level parameter, e.g. msetup in ADVCreate)
	PetscCall(getStringParam(fb, _OPTIONAL_, "dylib_plugin", lib, NULL));
	if(!strlen(lib)) PetscFunctionReturn(0); // later calls with the option set are still honoured

	if(initTried)
	{
		if(active && strcmp(lib, loadedPath) != 0)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP,
			        "dylib_plugin: a different plugin ('%s') was already loaded in this process ('%s')",
			        lib, loadedPath);
		}
		PetscFunctionReturn(0);
	}

	initTried = PETSC_TRUE;
	PetscCall(PetscStrncpy(loadedPath, lib, _str_len_));
	PetscCall(PetscRegisterFinalize(DylibPluginFinalize));

	// PETSC_DL_NOW -> dlopen(RTLD_NOW|RTLD_GLOBAL), pulling in libjulia as
	// the plugin's own dependency; LaMEM never links libjulia itself
	PetscCall(PetscDLOpen(lib, PETSC_DL_NOW, &handle));

	PetscCall(PetscDLSym(handle, "jl_parse_opts", &sym));
	if(!sym) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB, "dylib_plugin: jl_parse_opts not found in %s", lib);

	// --handle-signals=no keeps PETSc's SIGSEGV/SIGBUS handler in charge;
	// --threads=1 --gcthreads=1 keeps the embedded runtime single-threaded
	{
		static char a0[] = "lamem", a1[] = "--handle-signals=no", a2[] = "--threads=1", a3[] = "--gcthreads=1";
		static char *jlargv[4] = { a0, a1, a2, a3 };
		char **jlargvp = jlargv;
		int    jlargc  = 4;

		((JlParseOptsFn)sym)(&jlargc, &jlargvp);
	}

	PetscCall(DylibPluginSnapshotLbt());
	{
		void *getSym = NULL;
		PetscCall(PetscDLSym(NULL, "lbt_get_num_threads", &getSym));
		if(getSym) nthreadsBefore = ((LbtGetNumThreadsFn)getSym)();
	}

	PetscCall(PetscDLSym(handle, "jl_init_with_image_handle", &sym));
	if(!sym) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB, "dylib_plugin: jl_init_with_image_handle not found in %s", lib);
	((JlInitWithImageHandleFn)sym)((void*)handle);

	PetscCall(DylibPluginRestoreLbt(nthreadsBefore));

	PetscCall(PetscDLSym(handle, "lamem_plugin_abi_version", &sym));
	if(sym)
	{
		int v = ((DylibPluginAbiVersionFn)sym)();
		if(v != DYLIB_PLUGIN_ABI_VERSION)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB,
			        "dylib_plugin: %s reports ABI v%d, expected v%d", lib, v, DYLIB_PLUGIN_ABI_VERSION);
		}
	}

	PetscCall(PetscDLSym(handle, "lamem_phase_transition", &sym));
	if(!sym) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB, "dylib_plugin: lamem_phase_transition not found in %s", lib);
	fn = (DylibPluginFn)sym;

	PetscCall(PetscDLSym(handle, "jl_atexit_hook", &sym));
	atexitFn = (JlAtexitHookFn)sym; // optional

	active = PETSC_TRUE;

	PetscPrintf(PETSC_COMM_WORLD, "Dylib plugin  : %s\n", lib);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
static PetscErrorCode DylibPluginEnsureMarkerCapacity(PetscInt n)
{
	int i;

	PetscFunctionBeginUser;

	if(n <= bufcap) PetscFunctionReturn(0);

	for(i = 0; i < nMarkerBufs; i++) PetscCall(PetscFree(*markerBufs[i]));
	PetscCall(PetscFree(bphase_in)); PetscCall(PetscFree(bphase_out));

	for(i = 0; i < nMarkerBufs; i++) PetscCall(PetscMalloc((size_t)n*sizeof(PetscScalar), markerBufs[i]));
	PetscCall(PetscMalloc((size_t)n*sizeof(int32_t), &bphase_in));
	PetscCall(PetscMalloc((size_t)n*sizeof(int32_t), &bphase_out));

	bufcap = n;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
static PetscErrorCode DylibPluginEnsureCellCapacity(PetscInt ncells)
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
// Fill+ghost-exchange one edge DA's stress/strain-rate fields from svEdge[]
static PetscErrorCode DylibPluginFillEdge(DM da, SolVarEdge *svEdge, PetscScalar pf, Vec ls, Vec ld)
{
	PetscScalar ***as, ***ad;
	PetscInt     i, j, k, sx, sy, sz, nx, ny, nz, iter = 0;

	PetscFunctionBeginUser;

	PetscCall(DMDAVecGetArray(da, ls, &as));
	PetscCall(DMDAVecGetArray(da, ld, &ad));
	PetscCall(DMDAGetCorners(da, &sx, &sy, &sz, &nx, &ny, &nz));

	START_STD_LOOP
	{
		SolVarEdge *e = &svEdge[iter++];
		as[k][j][i] = e->s + pf*e->svDev.eta_st*e->d;
		ad[k][j][i] = e->d;
	}
	END_STD_LOOP

	PetscCall(DMDAVecRestoreArray(da, ls, &as));
	PetscCall(DMDAVecRestoreArray(da, ld, &ad));
	LOCAL_TO_LOCAL(da, ls);
	LOCAL_TO_LOCAL(da, ld);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// Cell-centred (not ParaView's corner-centred) J2 of deviatoric stress and
// strain rate: cell diagonal components plus its 4 surrounding edge values
// per off-diagonal direction, averaged (same geometry as JacResGetSHmax).
static PetscErrorCode DylibPluginComputeJ2(AdvCtx *actx)
{
	FDSTAG      *fs = actx->fs;
	JacRes      *jr = actx->jr;
	Vec          lxy_s, lxz_s, lyz_s, lxy_d, lxz_d, lyz_d;
	PetscScalar ***axy_s, ***axz_s, ***ayz_s, ***axy_d, ***axz_d, ***ayz_d;
	PetscInt     i, j, k, sx, sy, sz, nx, ny, nz, iter = 0;
	PetscScalar  pf = jr->ctrl.initGuess ? 0.0 : 2.0;

	PetscFunctionBeginUser;

	PetscCall(DylibPluginEnsureCellCapacity(fs->nCells));

	PetscCall(FDSTAGGetLocalVectorEdge(fs, &lxy_s, &lxz_s, &lyz_s));
	PetscCall(FDSTAGGetLocalVectorEdge(fs, &lxy_d, &lxz_d, &lyz_d));

	PetscCall(DylibPluginFillEdge(fs->DA_XY, jr->svXYEdge, pf, lxy_s, lxy_d));
	PetscCall(DylibPluginFillEdge(fs->DA_XZ, jr->svXZEdge, pf, lxz_s, lxz_d));
	PetscCall(DylibPluginFillEdge(fs->DA_YZ, jr->svYZEdge, pf, lyz_s, lyz_d));

	PetscCall(DMDAVecGetArray(fs->DA_XY, lxy_s, &axy_s)); PetscCall(DMDAVecGetArray(fs->DA_XY, lxy_d, &axy_d));
	PetscCall(DMDAVecGetArray(fs->DA_XZ, lxz_s, &axz_s)); PetscCall(DMDAVecGetArray(fs->DA_XZ, lxz_d, &axz_d));
	PetscCall(DMDAVecGetArray(fs->DA_YZ, lyz_s, &ayz_s)); PetscCall(DMDAVecGetArray(fs->DA_YZ, lyz_d, &ayz_d));

	PetscCall(DMDAGetCorners(fs->DA_CEN, &sx, &sy, &sz, &nx, &ny, &nz));
	START_STD_LOOP
	{
		SolVarCell *c = &jr->svCell[iter];

		PetscScalar sxy = (axy_s[k][j][i] + axy_s[k][j][i+1] + axy_s[k][j+1][i] + axy_s[k][j+1][i+1])/4.0;
		PetscScalar sxz = (axz_s[k][j][i] + axz_s[k][j][i+1] + axz_s[k+1][j][i] + axz_s[k+1][j][i+1])/4.0;
		PetscScalar syz = (ayz_s[k][j][i] + ayz_s[k][j+1][i] + ayz_s[k+1][j][i] + ayz_s[k+1][j+1][i])/4.0;
		PetscScalar dxy = (axy_d[k][j][i] + axy_d[k][j][i+1] + axy_d[k][j+1][i] + axy_d[k][j+1][i+1])/4.0;
		PetscScalar dxz = (axz_d[k][j][i] + axz_d[k][j][i+1] + axz_d[k+1][j][i] + axz_d[k+1][j][i+1])/4.0;
		PetscScalar dyz = (ayz_d[k][j][i] + ayz_d[k][j+1][i] + ayz_d[k+1][j][i] + ayz_d[k+1][j+1][i])/4.0;

		PetscScalar sxx = c->sxx + pf*c->svDev.eta_st*c->dxx;
		PetscScalar syy = c->syy + pf*c->svDev.eta_st*c->dyy;
		PetscScalar szz = c->szz + pf*c->svDev.eta_st*c->dzz;

		PetscScalar J2s = 0.5*(sxx*sxx + syy*syy + szz*szz) + sxy*sxy + sxz*sxz + syz*syz;
		PetscScalar J2e = 0.5*(c->dxx*c->dxx + c->dyy*c->dyy + c->dzz*c->dzz) + dxy*dxy + dxz*dxz + dyz*dyz;

		cellJ2Stress[iter]     = sqrt(J2s);
		cellJ2StrainRate[iter] = sqrt(J2e);

		iter++;
	}
	END_STD_LOOP

	PetscCall(DMDAVecRestoreArray(fs->DA_XY, lxy_s, &axy_s)); PetscCall(DMDAVecRestoreArray(fs->DA_XY, lxy_d, &axy_d));
	PetscCall(DMDAVecRestoreArray(fs->DA_XZ, lxz_s, &axz_s)); PetscCall(DMDAVecRestoreArray(fs->DA_XZ, lxz_d, &axz_d));
	PetscCall(DMDAVecRestoreArray(fs->DA_YZ, lyz_s, &ayz_s)); PetscCall(DMDAVecRestoreArray(fs->DA_YZ, lyz_d, &ayz_d));
	PetscCall(FDSTAGRestoreLocalVectorEdge(fs, &lxy_s, &lxz_s, &lyz_s));
	PetscCall(FDSTAGRestoreLocalVectorEdge(fs, &lxy_d, &lxz_d, &lyz_d));

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode DylibPluginPhaseTransition(AdvCtx *actx)
{
	JacRes      *jr = actx->jr;
	Scaling     *scal = jr->scal;
	Marker      *P;
	PetscInt     i, ID, n = actx->nummark, numPhases = actx->dbm->numPhases;
	int          rc;
	PetscInt     errFlagLoc = 0, err2[2], glob2[2];
	PetscInt     badIdx = -1;
	int32_t      badPhase = 0;
	PetscBool    badT = PETSC_FALSE;
	LaMEMPluginScaling scaling;

	PetscFunctionBeginUser;

	if(!active || !actx->dbm->dylibPhaseTr) PetscFunctionReturn(0);

	PetscCall(DylibPluginEnsureMarkerCapacity(n));
	PetscCall(DylibPluginComputeJ2(actx));

	scaling.abi_version = DYLIB_PLUGIN_ABI_VERSION;
	scaling.utype       = (int32_t)scal->utype;
	scaling.length      = (double)scal->length;
	scaling.time        = (double)scal->time;
	scaling.stress      = (double)scal->stress;
	scaling.temperature = (double)scal->temperature;
	scaling.viscosity   = (double)scal->viscosity;
	scaling.strain_rate = (double)scal->strain_rate;
	scaling.velocity    = (double)scal->velocity;
	scaling.density     = (double)scal->density;
	scaling.Tshift      = (double)scal->Tshift;
	scaling.pShift      = (double)jr->ctrl.pShift;
	scaling.dt          = (double)jr->bc->ts->dt;
	scaling.step        = (int64_t)jr->bc->ts->istep;

	// arrays are LaMEM's internal (non-dimensional) units, unconverted;
	// the plugin dimensionalises using `scaling` (see dylib_plugins.h)
	for(i = 0; i < n; i++)
	{
		P  = &actx->markers[i];
		ID = actx->cellnum[i];

		bx[i] = P->X[0]; by[i] = P->X[1]; bz[i] = P->X[2];
		bT[i] = P->T;    bp[i] = P->p; // raw, without pShift

		bsxx[i] = P->S.xx; bsyy[i] = P->S.yy; bszz[i] = P->S.zz;
		bsxy[i] = P->S.xy; bsxz[i] = P->S.xz; bsyz[i] = P->S.yz;

		bj2s[i] = cellJ2Stress[ID];
		bj2e[i] = cellJ2StrainRate[ID];
		beta[i] = jr->svCell[ID].svDev.eta;
		baps[i] = jr->svCell[ID].svDev.APS;

		bphase_in[i]  = (int32_t)P->phase;
		bphase_out[i] = (int32_t)P->phase;
		bT_out[i]     = bT[i];
	}

	// every rank calls fn() and joins every collective below, even if n==0
	rc = fn((size_t)n,
	        bx, by, bz, bT, bp, (double)jr->bc->ts->time,
	        bsxx, bsyy, bszz, bsxy, bsxz, bsyz,
	        bj2s, bj2e, beta, baps,
	        bphase_in, bphase_out, bT_out,
	        &scaling);

	// validate before any write-back, so a per-rank failure cannot desync
	// the collectives below (all ranks must reach the same MPI_Allreduce)
	if(rc < 0)
	{
		errFlagLoc = 1;
	}
	else for(i = 0; i < n; i++)
		{
			if(bphase_out[i] != bphase_in[i] && (bphase_out[i] < 0 || bphase_out[i] >= numPhases))
			{
				errFlagLoc = 1; badIdx = i; badPhase = bphase_out[i]; break;
			}
			if(!std::isfinite((double)bT_out[i]))
			{
				errFlagLoc = 1; badIdx = i; badT = PETSC_TRUE; break;
			}
		}

	PetscCallMPI(MPI_Allreduce(&errFlagLoc, &glob2[0], 1, MPIU_INT, MPI_MAX, PETSC_COMM_WORLD));

	if(glob2[0])
	{
		if(rc < 0)
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB, "dylib_plugin: plugin failed on at least one rank (rc=%d)", rc);
		else if(badIdx >= 0 && badT)
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "dylib_plugin: non-finite T_out at local marker %" PetscInt_FMT, badIdx);
		else if(badIdx >= 0)
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "dylib_plugin: out-of-range phase %d at local marker %" PetscInt_FMT, badPhase, badIdx);
		else
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB, "dylib_plugin: another MPI rank reported a plugin failure");
	}

	err2[0] = 0; err2[1] = 0; // reused as local changed-phase/changed-T counters

	for(i = 0; i < n; i++)
	{
		P = &actx->markers[i];

		if(bphase_out[i] != bphase_in[i]) { P->phase = (PetscInt)bphase_out[i]; err2[0]++; }
		if(bT_out[i] != bT[i])            { P->T     = bT_out[i];               err2[1]++; }
	}

	PetscCallMPI(MPI_Allreduce(err2, glob2, 2, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD));

	// ADVInterpMarkToCell also seeds svBulk.Tn from marker T, which the next
	// JacResInitTemp call reads - must run on a T-only change too, not just phase
	if(glob2[0] || glob2[1])
	{
		PetscCall(ADVCheckMarkPhases(actx));
		PetscCall(ADVInterpMarkToCell(actx));
	}

	PetscPrintf(PETSC_COMM_WORLD, "Dylib plugin  : %" PetscInt_FMT " marker(s) changed phase, "
	            "%" PetscInt_FMT " marker(s) changed temperature\n", glob2[0], glob2[1]);

	PetscFunctionReturn(0);
}

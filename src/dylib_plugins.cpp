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
#if defined(PETSC_HAVE_DLADDR) && !defined(_WIN32)
#include <dlfcn.h>
#endif
#if defined(_WIN32)
// WIN32_LEAN_AND_MEAN skips windows.h's COM/OLE headers (objidl.h,
// oaidl.h, ...), which otherwise declare a global `byte` typedef that
// collides with C++17's std::byte under this file's `using namespace std`
// (from LaMEM.h) - and those headers are not needed for GetModuleHandleA.
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h> // GetModuleHandleA, see DylibPluginOpenWinLbt
#endif
//---------------------------------------------------------------------------
// the ABI structs have no implicit padding: on 64-bit platforms their sizes
// are fixed (the plugin reports its own via lamem_plugin_struct_sizes)
static_assert(sizeof(void*) != 8 || sizeof(LaMEMPluginMarkers) == 256, "LaMEMPluginMarkers layout changed");
static_assert(sizeof(void*) != 8 || sizeof(LaMEMPluginCells)   == 344, "LaMEMPluginCells layout changed");
static_assert(sizeof(LaMEMPluginStep)    == 24,  "LaMEMPluginStep layout changed");
static_assert(sizeof(LaMEMPluginScaling) == 112, "LaMEMPluginScaling layout changed");

typedef int32_t  (*DylibPluginAbiVersionFn)(void);
typedef int32_t  (*DylibPluginStructSizesFn)(int64_t *sizes, int32_t n);
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
#if defined(_WIN32)
PetscDLHandle libjuliaHandle = NULL, lbtHandle = NULL; // sibling DLLs, see DylibPluginOpenWinJulia/DylibPluginOpenWinLbt
#define LbtSymHandle lbtHandle // GetProcAddress needs the specific DLL's own handle
#else
#define LbtSymHandle NULL // dlsym(NULL, ...) searches the whole process (RTLD_GLOBAL)
#endif

// per-marker buffers (grown, never shrunk); the writable fields have an
// input copy (*_in) and an output array (*_out), in the same order
enum
{
	MB_X, MB_Y, MB_Z, MB_P,                         // read-only
	MB_IN,                                          // first writable input
	MB_OUT = MB_IN + 12,                            // first writable output
	MB_NUM = MB_OUT + 12
};
const int   nWritable = MB_OUT - MB_IN;
const char *writableName[] = { "T", "aps", "ats", "sxx", "syy", "szz", "sxy", "sxz", "syz", "ux", "uy", "uz" };

PetscInt  bufcap = 0;
double   *mblock = NULL; // backs mb[] and the three int32 arrays, see DylibPluginAllocBlock
double   *mb[MB_NUM];
int32_t  *mcell = NULL, *mphase_in = NULL, *mphase_out = NULL;

// per-cell buffers, copied once per step
enum
{
	CB_ETA, CB_ETA_ST, CB_I2GDT, CB_HR, CB_APS, CB_PSR,
	CB_THETA, CB_RHO, CB_IKDT, CB_ALPHA, CB_TN, CB_PN, CB_RHO_PF, CB_MF, CB_PHI, CB_HA, CB_COND,
	CB_SXX, CB_SYY, CB_SZZ, CB_HXX, CB_HYY, CB_HZZ, CB_DXX, CB_DYY, CB_DZZ,
	CB_UX, CB_UY, CB_UZ, CB_ATS, CB_ETA_CR,
	CB_DIIDIF, CB_DIIDIS, CB_DIIPRL, CB_DIIFK, CB_DIIPL, CB_YIELD,
	CB_J2S, CB_J2E,
	CB_NUM
};

PetscInt  cellcap = 0, phRatcap = 0;
double   *cblock = NULL; // backs cb[] and cfreesurf
double   *cb[CB_NUM];
int32_t  *cfreesurf = NULL;
double   *cphRat = NULL;

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
	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
static PetscErrorCode DylibPluginFreeBuffers(void)
{
	PetscFunctionBeginUser;

	PetscCall(PetscFree(mblock));
	PetscCall(PetscFree(cblock));
	PetscCall(PetscFree(cphRat));

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

	bufcap = 0; cellcap = 0; phRatcap = 0;

	PetscFunctionReturn(0);
}
#if defined(_WIN32)
//---------------------------------------------------------------------------
// On POSIX, PetscDLOpen(PETSC_DL_NOW) uses RTLD_GLOBAL, so the plugin's own
// dependencies (libjulia, libblastrampoline) are reachable via PetscDLSym. Windows
// has no equivalent: GetProcAddress only searches one specific DLL's exports, so
// those DLLs need handles of their own. JuliaC places them flat next to the plugin.
//
// libjulia.dll must be opened from the plugin's directory BEFORE the plugin itself.
// Windows binds imports by base name only: if the plugin were loaded first, its
// libjulia/libjulia-internal imports could bind to an unrelated Julia found on
// PATH, and jl_init_with_image_handle would fail its image consistency check
// against that mismatched runtime.
static PetscErrorCode DylibPluginOpenWinJulia(const char *pluginPath)
{
	char  dir[_str_len_], path[_str_len_];
	char *slash;

	PetscFunctionBeginUser;

	PetscCall(PetscStrncpy(dir, pluginPath, _str_len_));
	slash = strrchr(dir, '\\');
	if(!slash) slash = strrchr(dir, '/');
	if(slash) *slash = '\0'; else dir[0] = '\0'; // no directory given: current dir

	PetscCall(PetscSNPrintf(path, _str_len_, "%s%slibjulia.dll", dir, dir[0] ? "\\" : ""));
	PetscCall(PetscDLOpen(path, PETSC_DL_NOW, &libjuliaHandle));
	if(!libjuliaHandle) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB, "dylib_plugin: %s not found next to %s", path, pluginPath);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// libblastrampoline-5.dll must NOT be opened a second time: PETSc's own BLAS calls
// are already bound to the copy loaded with LaMEM.exe, and opening the plugin's
// copy by a different path would create a distinct module - Julia would then set up
// forwarding on that one while PETSc's calls stay on the unconfigured first ("no
// BLAS/LAPACK library loaded"). GetModuleHandleA returns the already-loaded
// module's handle, which is the same void* PetscDLOpen would have produced.
static PetscErrorCode DylibPluginOpenWinLbt(void)
{
	PetscFunctionBeginUser;

	lbtHandle = (PetscDLHandle)GetModuleHandleA("libblastrampoline-5.dll"); // NULL if not loaded, not an error

	PetscFunctionReturn(0);
}
#endif
//---------------------------------------------------------------------------
// Snapshot libblastrampoline's forwarding table before jl_init_with_image_handle
// runs: Julia's LinearAlgebra.__init__ calls lbt_forward(clear=1), which both
// overwrites the table and frees the strings behind any earlier snapshot, so
// they must be copied out first. lbt_get_config's struct layout is not
// versioned; trust it only if the resolved symbol lives in a
// "libblastrampoline.5" (macOS) or "libblastrampoline.so.5" (Linux) image
// (dladdr), else skip silently. On Windows lbtHandle already IS the specific
// "libblastrampoline-5.dll" module DylibPluginOpenWinLbt looked up by that
// exact name (or NULL if none is loaded), so there is no address to verify
// - the explicit name lookup is the check.
static PetscErrorCode DylibPluginSnapshotLbt(void)
{
	void *sym = NULL;

	PetscFunctionBeginUser;

	lbtSnapCount = 0;

	PetscCall(PetscDLSym(LbtSymHandle, "lbt_get_config", &sym));
	if(!sym) PetscFunctionReturn(0);

#if !defined(_WIN32)
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

	PetscCall(PetscDLSym(LbtSymHandle, "lbt_forward", &fwdSym));

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

	PetscCall(PetscDLSym(LbtSymHandle, "lbt_get_num_threads", &getSym));
	PetscCall(PetscDLSym(LbtSymHandle, "lbt_set_num_threads", &setSym));
	if(getSym && setSym && nthreadsBefore > 0) ((LbtSetNumThreadsFn)setSym)(nthreadsBefore);

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode DylibPluginLoad(AdvCtx *actx, FB *fb)
{
	char      lib[_str_len_];
	void     *sym;
	int32_t   nthreadsBefore = -1;
	int       abiVersion     = -1;

	PetscFunctionBeginUser;

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

#if defined(_WIN32)
	PetscCall(DylibPluginOpenWinJulia(lib)); // must precede opening the plugin, see there
#endif

	// PETSC_DL_NOW -> dlopen(RTLD_NOW|RTLD_GLOBAL), pulling in libjulia as
	// the plugin's own dependency; LaMEM never links libjulia itself
	PetscCall(PetscDLOpen(lib, PETSC_DL_NOW, &handle));
	if(!handle) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN, "dylib_plugin: could not open %s", lib);

#if defined(_WIN32)
	// RTLD_GLOBAL has no Windows equivalent: libjulia's own exports (jl_parse_opts
	// among them) are not visible via the plugin DLL's handle. See
	// DylibPluginOpenWinJulia above.
	PetscCall(PetscDLSym(libjuliaHandle, "jl_parse_opts", &sym));
#else
	PetscCall(PetscDLSym(handle, "jl_parse_opts", &sym));
#endif
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

#if defined(_WIN32)
	PetscCall(DylibPluginOpenWinLbt());
#endif

	PetscCall(DylibPluginSnapshotLbt());
	{
		void *getSym = NULL;
		PetscCall(PetscDLSym(LbtSymHandle, "lbt_get_num_threads", &getSym));
		if(getSym) nthreadsBefore = ((LbtGetNumThreadsFn)getSym)();
	}

#if defined(_WIN32)
	PetscCall(PetscDLSym(libjuliaHandle, "jl_init_with_image_handle", &sym));
#else
	PetscCall(PetscDLSym(handle, "jl_init_with_image_handle", &sym));
#endif
	if(!sym) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB, "dylib_plugin: jl_init_with_image_handle not found in %s", lib);
	((JlInitWithImageHandleFn)sym)((void*)handle);

	PetscCall(DylibPluginRestoreLbt(nthreadsBefore));

	// a plugin built for another ABI would read/write the structs with a
	// different layout: both checks are mandatory, before the first call
	PetscCall(PetscDLSym(handle, "lamem_plugin_abi_version", &sym));
	if(!sym)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB,
		        "dylib_plugin: %s does not export lamem_plugin_abi_version (built for ABI v1?), expected v%d - rebuild it with this LaMEM's LaMEMPlugin.jl",
		        lib, DYLIB_PLUGIN_ABI_VERSION);
	}
	abiVersion = (int)((DylibPluginAbiVersionFn)sym)();
	if(abiVersion != DYLIB_PLUGIN_ABI_VERSION)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB,
		        "dylib_plugin: %s reports ABI v%d, expected v%d", lib, abiVersion, DYLIB_PLUGIN_ABI_VERSION);
	}

	PetscCall(PetscDLSym(handle, "lamem_plugin_struct_sizes", &sym));
	if(!sym) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB, "dylib_plugin: lamem_plugin_struct_sizes not found in %s", lib);
	{
		const char *names[DYLIB_PLUGIN_NUM_STRUCTS] = { "LaMEMPluginMarkers", "LaMEMPluginCells", "LaMEMPluginStep", "LaMEMPluginScaling" };
		int64_t     ours [DYLIB_PLUGIN_NUM_STRUCTS] = { (int64_t)sizeof(LaMEMPluginMarkers), (int64_t)sizeof(LaMEMPluginCells),
		                                                (int64_t)sizeof(LaMEMPluginStep),    (int64_t)sizeof(LaMEMPluginScaling)
		                                              };
		int64_t     theirs[DYLIB_PLUGIN_NUM_STRUCTS] = { -1, -1, -1, -1 };
		int32_t     nret;
		int         i;

		nret = ((DylibPluginStructSizesFn)sym)(theirs, DYLIB_PLUGIN_NUM_STRUCTS);
		if(nret != DYLIB_PLUGIN_NUM_STRUCTS)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB,
			        "dylib_plugin: %s describes %d ABI structs, expected %d", lib, (int)nret, DYLIB_PLUGIN_NUM_STRUCTS);
		}
		for(i = 0; i < DYLIB_PLUGIN_NUM_STRUCTS; i++)
		{
			if(theirs[i] != ours[i])
			{
				SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB,
				        "dylib_plugin: struct layout mismatch in %s: sizeof(%s) is %lld bytes in the plugin, %lld in LaMEM - rebuild it with this LaMEM's LaMEMPlugin.jl",
				        lib, names[i], (long long)theirs[i], (long long)ours[i]);
			}
		}
	}

	PetscCall(PetscDLSym(handle, "lamem_phase_transition", &sym));
	if(!sym) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB, "dylib_plugin: lamem_phase_transition not found in %s", lib);
	fn = (DylibPluginFn)sym;

	PetscCall(PetscDLSym(handle, "jl_atexit_hook", &sym));
	atexitFn = (JlAtexitHookFn)sym; // optional

	active = PETSC_TRUE;

	PetscPrintf(PETSC_COMM_WORLD, "Dylib plugin parameters:\n");
	PetscPrintf(PETSC_COMM_WORLD, "   Library                                 : %s\n", lib);
	PetscPrintf(PETSC_COMM_WORLD, "   Plugin ABI version                      : %d\n", abiVersion);
	if(actx->dbm->dylibPhaseTr) PetscPrintf(PETSC_COMM_WORLD, "   Phase transitions                       : dylib (lamem_phase_transition is called every step, after the built-in transitions)\n");
	else                        PetscPrintf(PETSC_COMM_WORLD, "   Phase transitions                       : builtin\n");
	PetscPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------\n");

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// One block for k arrays of cap 8-byte entries each, staggered by 128 bytes
// modulo 16 KiB. Separate large mallocs are page aligned, so the 30-40
// arrays the gather/scatter loops stream through in parallel would all map
// to the same cache sets; that made those loops several times slower.
static PetscErrorCode DylibPluginAllocBlock(PetscInt cap, int k, double **block, double **arr)
{
	size_t stride = ((size_t)cap + 2047)/2048*2048 + 16; // in doubles
	int    i;

	PetscFunctionBeginUser;

	PetscCall(PetscFree(*block));
	PetscCall(PetscMalloc((size_t)k*stride*sizeof(double), block));

	for(i = 0; i < k; i++) arr[i] = *block + (size_t)i*stride;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
static PetscErrorCode DylibPluginEnsureMarkerCapacity(PetscInt n)
{
	double *arr[MB_NUM + 3];
	int     i;

	PetscFunctionBeginUser;

	if(n <= bufcap) PetscFunctionReturn(0);

	PetscCall(DylibPluginAllocBlock(n, MB_NUM + 3, &mblock, arr));

	for(i = 0; i < MB_NUM; i++) mb[i] = arr[i];
	mcell      = (int32_t*)arr[MB_NUM];
	mphase_in  = (int32_t*)arr[MB_NUM + 1];
	mphase_out = (int32_t*)arr[MB_NUM + 2];

	bufcap = n;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
static PetscErrorCode DylibPluginEnsureCellCapacity(PetscInt ncells, PetscInt numPhases)
{
	double *arr[CB_NUM + 1];
	int     i;

	PetscFunctionBeginUser;

	if(ncells > cellcap)
	{
		PetscCall(DylibPluginAllocBlock(ncells, CB_NUM + 1, &cblock, arr));

		for(i = 0; i < CB_NUM; i++) cb[i] = arr[i];
		cfreesurf = (int32_t*)arr[CB_NUM];

		cellcap = ncells;
	}

	if(ncells*numPhases > phRatcap)
	{
		PetscCall(PetscFree(cphRat));
		PetscCall(PetscMalloc((size_t)(ncells*numPhases)*sizeof(double), &cphRat));

		phRatcap = ncells*numPhases;
	}

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
// Writes cb[CB_J2S]/cb[CB_J2E]; the caller ensures their capacity.
static PetscErrorCode DylibPluginComputeJ2(AdvCtx *actx)
{
	FDSTAG      *fs = actx->fs;
	JacRes      *jr = actx->jr;
	Vec          lxy_s, lxz_s, lyz_s, lxy_d, lxz_d, lyz_d;
	PetscScalar ***axy_s, ***axz_s, ***ayz_s, ***axy_d, ***axz_d, ***ayz_d;
	PetscInt     i, j, k, sx, sy, sz, nx, ny, nz, iter = 0;
	PetscScalar  pf = jr->ctrl.initGuess ? 0.0 : 2.0;

	PetscFunctionBeginUser;

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

		cb[CB_J2S][iter] = sqrt(J2s);
		cb[CB_J2E][iter] = sqrt(J2e);

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
// Read-only structure-of-arrays copy of all local cells, once per step
static PetscErrorCode DylibPluginGatherCells(AdvCtx *actx, LaMEMPluginCells *cells)
{
	JacRes   *jr        = actx->jr;
	PetscInt  nCells    = actx->fs->nCells;
	PetscInt  numPhases = actx->dbm->numPhases;
	PetscInt  c, ph;

	PetscFunctionBeginUser;

	PetscCall(DylibPluginEnsureCellCapacity(nCells, numPhases));
	PetscCall(DylibPluginComputeJ2(actx)); // fills cb[CB_J2S], cb[CB_J2E]

	for(c = 0; c < nCells; c++)
	{
		const SolVarCell *sv = &jr->svCell[c];

		cb[CB_ETA]   [c] = sv->svDev.eta;
		cb[CB_ETA_ST][c] = sv->svDev.eta_st;
		cb[CB_I2GDT] [c] = sv->svDev.I2Gdt;
		cb[CB_HR]    [c] = sv->svDev.Hr;
		cb[CB_APS]   [c] = sv->svDev.APS;
		cb[CB_PSR]   [c] = sv->svDev.PSR;

		cb[CB_THETA] [c] = sv->svBulk.theta;
		cb[CB_RHO]   [c] = sv->svBulk.rho;
		cb[CB_IKDT]  [c] = sv->svBulk.IKdt;
		cb[CB_ALPHA] [c] = sv->svBulk.alpha;
		cb[CB_TN]    [c] = sv->svBulk.Tn;
		cb[CB_PN]    [c] = sv->svBulk.pn;
		cb[CB_RHO_PF][c] = sv->svBulk.rho_pf;
		cb[CB_MF]    [c] = sv->svBulk.mf;
		cb[CB_PHI]   [c] = sv->svBulk.phi;
		cb[CB_HA]    [c] = sv->svBulk.Ha;
		cb[CB_COND]  [c] = sv->svBulk.cond;

		cb[CB_SXX][c] = sv->sxx; cb[CB_SYY][c] = sv->syy; cb[CB_SZZ][c] = sv->szz;
		cb[CB_HXX][c] = sv->hxx; cb[CB_HYY][c] = sv->hyy; cb[CB_HZZ][c] = sv->hzz;
		cb[CB_DXX][c] = sv->dxx; cb[CB_DYY][c] = sv->dyy; cb[CB_DZZ][c] = sv->dzz;
		cb[CB_UX] [c] = sv->U[0]; cb[CB_UY][c] = sv->U[1]; cb[CB_UZ][c] = sv->U[2];

		cb[CB_ATS]   [c] = sv->ATS;
		cb[CB_ETA_CR][c] = sv->eta_cr;
		cb[CB_DIIDIF][c] = sv->DIIdif;
		cb[CB_DIIDIS][c] = sv->DIIdis;
		cb[CB_DIIPRL][c] = sv->DIIprl;
		cb[CB_DIIFK] [c] = sv->DIIfk;
		cb[CB_DIIPL] [c] = sv->DIIpl;
		cb[CB_YIELD] [c] = sv->yield;

		cfreesurf[c] = (int32_t)sv->FreeSurf;

		for(ph = 0; ph < numPhases; ph++) cphRat[c*numPhases + ph] = sv->phRat[ph];
	}

	cells->ncells    = (size_t)nCells;
	cells->numPhases = (int32_t)numPhases;
	cells->reserved  = 0;

	cells->eta    = cb[CB_ETA];    cells->eta_st = cb[CB_ETA_ST]; cells->I2Gdt = cb[CB_I2GDT];
	cells->Hr     = cb[CB_HR];     cells->aps    = cb[CB_APS];    cells->psr   = cb[CB_PSR];

	cells->theta  = cb[CB_THETA];  cells->rho    = cb[CB_RHO];    cells->IKdt  = cb[CB_IKDT];
	cells->alpha  = cb[CB_ALPHA];  cells->Tn     = cb[CB_TN];     cells->pn    = cb[CB_PN];
	cells->rho_pf = cb[CB_RHO_PF]; cells->mf     = cb[CB_MF];     cells->phi   = cb[CB_PHI];
	cells->Ha     = cb[CB_HA];     cells->cond   = cb[CB_COND];

	cells->sxx = cb[CB_SXX]; cells->syy = cb[CB_SYY]; cells->szz = cb[CB_SZZ];
	cells->hxx = cb[CB_HXX]; cells->hyy = cb[CB_HYY]; cells->hzz = cb[CB_HZZ];
	cells->dxx = cb[CB_DXX]; cells->dyy = cb[CB_DYY]; cells->dzz = cb[CB_DZZ];

	cells->free_surf = cfreesurf;
	cells->ux        = cb[CB_UX]; cells->uy = cb[CB_UY]; cells->uz = cb[CB_UZ];
	cells->ats       = cb[CB_ATS];
	cells->eta_cr    = cb[CB_ETA_CR];
	cells->DIIdif    = cb[CB_DIIDIF]; cells->DIIdis = cb[CB_DIIDIS]; cells->DIIprl = cb[CB_DIIPRL];
	cells->DIIfk     = cb[CB_DIIFK];  cells->DIIpl  = cb[CB_DIIPL];
	cells->yield     = cb[CB_YIELD];

	cells->j2_stress     = cb[CB_J2S];
	cells->j2_strainrate = cb[CB_J2E];

	cells->phRat = cphRat;

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode DylibPluginPhaseTransition(AdvCtx *actx)
{
	JacRes             *jr   = actx->jr;
	Scaling            *scal = jr->scal;
	Marker             *P;
	PetscInt            i, n = actx->nummark, numPhases = actx->dbm->numPhases;
	int                 f;
	int32_t             rc;
	PetscInt            errFlagLoc = 0, errFlagGlob = 0, cnt[3], glob[3];
	PetscInt            badIdx = -1;
	int                 badField = -1; // index into writableName, -1: phase
	int32_t             badPhase = 0;
	LaMEMPluginMarkers  markers;
	LaMEMPluginCells    cells;
	LaMEMPluginStep     step;
	LaMEMPluginScaling  scaling;

	PetscFunctionBeginUser;

	if(!active || !actx->dbm->dylibPhaseTr) PetscFunctionReturn(0);

	PetscCall(DylibPluginEnsureMarkerCapacity(n));
	PetscCall(DylibPluginGatherCells(actx, &cells));

	scaling.abi_version      = DYLIB_PLUGIN_ABI_VERSION;
	scaling.utype            = (int32_t)scal->utype;
	scaling.length           = (double)scal->length;
	scaling.time             = (double)scal->time;
	scaling.stress           = (double)scal->stress;
	scaling.temperature      = (double)scal->temperature;
	scaling.viscosity        = (double)scal->viscosity;
	scaling.strain_rate      = (double)scal->strain_rate;
	scaling.velocity         = (double)scal->velocity;
	scaling.density          = (double)scal->density;
	scaling.conductivity     = (double)scal->conductivity;
	scaling.expansivity      = (double)scal->expansivity;
	scaling.dissipation_rate = (double)scal->dissipation_rate;
	scaling.Tshift           = (double)scal->Tshift;
	scaling.pShift           = (double)jr->ctrl.pShift;

	step.time = (double)jr->bc->ts->time;
	step.dt   = (double)jr->bc->ts->dt;
	step.step = (int64_t)jr->bc->ts->istep;

	// one pass over the markers: copy every field (internal units, unconverted)
	// and pre-fill the outputs with the inputs
	{
		double         *x = mb[MB_X], *y = mb[MB_Y], *z = mb[MB_Z], *p = mb[MB_P];
		double         *in[MB_OUT - MB_IN], *out[MB_OUT - MB_IN];
		const PetscInt *cellnum = actx->cellnum;

		for(f = 0; f < nWritable; f++) { in[f] = mb[MB_IN + f]; out[f] = mb[MB_OUT + f]; }

		for(i = 0; i < n; i++)
		{
			P = &actx->markers[i];

			x[i] = P->X[0]; y[i] = P->X[1]; z[i] = P->X[2];
			p[i] = P->p; // raw, without pShift

			mcell[i]     = (int32_t)cellnum[i];
			mphase_in[i] = mphase_out[i] = (int32_t)P->phase;

			in[0] [i] = out[0] [i] = P->T;
			in[1] [i] = out[1] [i] = P->APS;
			in[2] [i] = out[2] [i] = P->ATS;
			in[3] [i] = out[3] [i] = P->S.xx;
			in[4] [i] = out[4] [i] = P->S.yy;
			in[5] [i] = out[5] [i] = P->S.zz;
			in[6] [i] = out[6] [i] = P->S.xy;
			in[7] [i] = out[7] [i] = P->S.xz;
			in[8] [i] = out[8] [i] = P->S.yz;
			in[9] [i] = out[9] [i] = P->U[0];
			in[10][i] = out[10][i] = P->U[1];
			in[11][i] = out[11][i] = P->U[2];
		}

		markers.n          = (size_t)n;
		markers.cell_index = mcell;
		markers.x          = x; markers.y = y; markers.z = z; markers.p = p;
		markers.phase_in   = mphase_in;
		markers.phase_out  = mphase_out;

		markers.T_in   = in[0];  markers.aps_in = in[1];  markers.ats_in = in[2];
		markers.sxx_in = in[3];  markers.syy_in = in[4];  markers.szz_in = in[5];
		markers.sxy_in = in[6];  markers.sxz_in = in[7];  markers.syz_in = in[8];
		markers.ux_in  = in[9];  markers.uy_in  = in[10]; markers.uz_in  = in[11];

		markers.T_out   = out[0]; markers.aps_out = out[1];  markers.ats_out = out[2];
		markers.sxx_out = out[3]; markers.syy_out = out[4];  markers.szz_out = out[5];
		markers.sxy_out = out[6]; markers.sxz_out = out[7];  markers.syz_out = out[8];
		markers.ux_out  = out[9]; markers.uy_out  = out[10]; markers.uz_out  = out[11];
	}

	// every rank calls fn() and joins every collective below, even if n==0
	rc = fn(&markers, &cells, &step, &scaling);

	// validate before any write-back, so a per-rank failure cannot desync
	// the collectives below (all ranks must reach the same MPI_Allreduce)
	if(rc < 0)
	{
		errFlagLoc = 1;
	}
	else
	{
		for(i = 0; i < n; i++)
		{
			if(mphase_out[i] != mphase_in[i] && (mphase_out[i] < 0 || mphase_out[i] >= numPhases))
			{
				errFlagLoc = 1; badIdx = i; badPhase = mphase_out[i]; break;
			}
		}
		for(f = 0; f < nWritable && !errFlagLoc; f++)
		{
			const double *out = mb[MB_OUT + f];

			for(i = 0; i < n; i++)
			{
				if(!std::isfinite(out[i])) { errFlagLoc = 1; badIdx = i; badField = f; break; }
			}
		}
	}

	PetscCallMPI(MPI_Allreduce(&errFlagLoc, &errFlagGlob, 1, MPIU_INT, MPI_MAX, PETSC_COMM_WORLD));

	if(errFlagGlob)
	{
		if(rc < 0)
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB, "dylib_plugin: plugin failed on at least one rank (rc=%d)", (int)rc);
		else if(badIdx >= 0 && badField >= 0)
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "dylib_plugin: non-finite %s_out at local marker %" PetscInt_FMT, writableName[badField], badIdx);
		else if(badIdx >= 0)
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "dylib_plugin: out-of-range phase %d at local marker %" PetscInt_FMT, badPhase, badIdx);
		else
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB, "dylib_plugin: another MPI rank reported a plugin failure");
	}

	// write back only what changed; count markers that changed phase,
	// temperature, and any other field
	cnt[0] = cnt[1] = cnt[2] = 0;
	{
		const double *in[MB_OUT - MB_IN], *out[MB_OUT - MB_IN];

		for(f = 0; f < nWritable; f++) { in[f] = mb[MB_IN + f]; out[f] = mb[MB_OUT + f]; }

		for(i = 0; i < n; i++)
		{
			PetscInt other = 0;

			P = &actx->markers[i];

			if(mphase_out[i] != mphase_in[i]) { P->phase = (PetscInt)mphase_out[i]; cnt[0]++; }
			if(out[0][i]     != in[0][i])     { P->T     = out[0][i];               cnt[1]++; }

			if(out[1] [i] != in[1] [i]) { P->APS  = out[1] [i]; other = 1; }
			if(out[2] [i] != in[2] [i]) { P->ATS  = out[2] [i]; other = 1; }
			if(out[3] [i] != in[3] [i]) { P->S.xx = out[3] [i]; other = 1; }
			if(out[4] [i] != in[4] [i]) { P->S.yy = out[4] [i]; other = 1; }
			if(out[5] [i] != in[5] [i]) { P->S.zz = out[5] [i]; other = 1; }
			if(out[6] [i] != in[6] [i]) { P->S.xy = out[6] [i]; other = 1; }
			if(out[7] [i] != in[7] [i]) { P->S.xz = out[7] [i]; other = 1; }
			if(out[8] [i] != in[8] [i]) { P->S.yz = out[8] [i]; other = 1; }
			if(out[9] [i] != in[9] [i]) { P->U[0] = out[9] [i]; other = 1; }
			if(out[10][i] != in[10][i]) { P->U[1] = out[10][i]; other = 1; }
			if(out[11][i] != in[11][i]) { P->U[2] = out[11][i]; other = 1; }

			cnt[2] += other;
		}
	}

	PetscCallMPI(MPI_Allreduce(cnt, glob, 3, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD));

	// ADVInterpMarkToCell maps every changed marker field back to the cells
	// (phRat, svBulk.Tn - read by the next JacResInitTemp -, APS, ATS, the
	// history stress and the displacement), so it must also run when only T
	// or another non-phase field changed
	if(glob[0] || glob[1] || glob[2])
	{
		PetscCall(ADVCheckMarkPhases(actx));
		PetscCall(ADVInterpMarkToCell(actx));
	}

	PetscPrintf(PETSC_COMM_WORLD, "Dylib plugin  : %" PetscInt_FMT " marker(s) changed phase, "
	            "%" PetscInt_FMT " marker(s) changed temperature, "
	            "%" PetscInt_FMT " marker(s) changed other fields\n", glob[0], glob[1], glob[2]);

	PetscFunctionReturn(0);
}

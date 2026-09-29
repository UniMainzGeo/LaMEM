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
#if defined(PETSC_HAVE_DLADDR)
#include <dlfcn.h>
#endif
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
// libblastrampoline's public C ABI (struct layout copied verbatim from its
// public header, e.g. <libblastrampoline_jll>/include/libblastrampoline.h or
// the libblastrampoline artifact's include/libblastrampoline.h - NOT
// guessed: field order/types were read directly from that header before
// writing this). LaMEM does not link libblastrampoline; these are only used
// to interpret the pointer lbt_get_config() returns, itself resolved at
// runtime via PetscDLSym.
struct LbtLibraryInfo
{
	char       *libname;
	void       *dlhandle;
	const char *suffix;
	uint8_t    *active_forwards;
	int32_t     interface;
	int32_t     complex_retstyle;
	int32_t     f2c;
	int32_t     cblas;
};

struct LbtConfig
{
	LbtLibraryInfo **loaded_libs;
	uint32_t          build_flags;
	const char      **exported_symbols;
	uint32_t          num_exported_symbols;
};

typedef const LbtConfig* (*LbtGetConfigFn)(void);
typedef int32_t          (*LbtForwardFn)(const char*, int32_t, int32_t, const char*);
typedef int32_t          (*LbtGetNumThreadsFn)(void);
typedef void              (*LbtSetNumThreadsFn)(int32_t);

// A snapshot of one libblastrampoline-forwarded library, taken BEFORE Julia
// init (Julia's own LinearAlgebra.__init__ calls lbt_forward(..., clear=1),
// which both overwrites the forwarding table AND frees the strings behind
// the previous lbt_get_config() snapshot - so the libname/suffix must be
// copied out, not just pointer-saved, before jl_init_with_image_handle runs).
struct LbtSnapshot
{
	char *libname; // PetscStrallocpy'd copy
	char *suffix;  // PetscStrallocpy'd copy, or NULL
};

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
	char          loadedPath[_str_len_] = ""; // path passed to the first successful PhTrPluginLoad

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

	// snapshot of PETSc's libblastrampoline forwarding table, taken just
	// before jl_init_with_image_handle, used to restore it afterwards (see
	// PhTrPluginLoad); freed once restored
	enum { LBT_MAX_SNAPSHOT = 16 };
	LbtSnapshot lbtSnap[LBT_MAX_SNAPSHOT];
	int         lbtSnapCount = 0;
}
//---------------------------------------------------------------------------
PetscBool PhTrPluginIsActive(void)
{
	return active;
}
//---------------------------------------------------------------------------
// PetscRegisterFinalize callback: releases plugin resources exactly once per
// process, at PetscFinalize() time (i.e. after the last LaMEMLibSolve() call,
// however many times it ran). Julia cannot be re-initialised once torn down.
//
// IMPORTANT: this deliberately does NOT call PetscDLClose(&handle). Once
// jl_init_with_image_handle has run, dlclose()-ing the library that holds
// the initialised Julia runtime is not supported by Julia (it does not
// expect its own image to be unloaded independently of process exit), so
// the handle is intentionally leaked for the life of the process; only
// jl_atexit_hook (best-effort) is called, once.
static PetscErrorCode PhTrPluginFinalize(void)
{
	PetscFunctionBeginUser;

	if(active && atexitFn)
	{
		// best-effort: some embedded Julia runtimes hang or crash on
		// finalize. See the Phase 1 report for the outcome observed here.
		atexitFn(0);
	}

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
// Snapshot PETSc's current libblastrampoline forwarding table (the
// libraries it points at, per-library suffix), BEFORE Julia touches
// anything. Must be called before jl_init_with_image_handle: Julia's
// LinearAlgebra.__init__ calls lbt_forward(..., clear=1), which both
// overwrites the table AND frees the very strings lbt_get_config() would
// otherwise still be pointing at, so the strings are copied out here.
//
// VERSION SAFETY: lbt_config_t/lbt_library_info_t (defined above, LbtConfig/
// LbtLibraryInfo) have no version field, and libblastrampoline exports no
// version-query symbol either, so there is no direct way to ask "is this
// struct layout the one I coded against?". The struct layout used here was
// copied field-for-field from libblastrampoline's own public header at
// version 5.15.0 (the version actually deployed here: the loaded image is
// /workspace/destdir/lib/libblastrampoline.5.dylib, a copy of the header
// exists in several Yggdrasil build-artifact trees on this machine, e.g.
// .../aarch64-apple-darwin20-libgfortran5-cxx11-mpi+mpitrampoline/destdir/
// include/libblastrampoline.h, matching the upstream
// github.com/JuliaLinearAlgebra/libblastrampoline "include/libblastrampoline.h"
// for that release). As the closest available guard, resolve lbt_get_config's
// OWN address with dladdr() and only trust the struct layout if the
// containing image's path names "libblastrampoline.5" (i.e. is recognisably
// an LBT 5.x build) -- this cannot detect a struct layout change WITHIN the
// 5.x series, only guards against a hypothetical future LBT 6+ that
// re-orders/extends the struct under the same-looking "libblastrampoline.5"
// name never being mistaken for 5.15.0's layout by a differently-tagged
// image name. If dladdr is unavailable (PETSC_HAVE_DLADDR undefined) or the
// image name doesn't match, this snapshot is skipped and
// PhTrPluginRestoreLbt falls back to -phase_transition_lbt_ilp64/-lp64 or
// LBT_DEFAULT_LIBS. A future libblastrampoline 6 (or any release that
// changes this struct) would need this code updated together with it.
static PetscErrorCode PhTrPluginSnapshotLbt(void)
{
	void *sym = NULL;

	PetscFunctionBeginUser;

	lbtSnapCount = 0;

	// lbt_get_config is exported by libblastrampoline (confirmed with `nm`
	// on the actual libblastrampoline.5.dylib in this deployment); resolved
	// process-wide via PetscDLSym(NULL, ...) rather than through the plugin
	// handle, for the same reason lbt_forward is (see PhTrPluginRestoreLbt).
	PetscCall(PetscDLSym(NULL, "lbt_get_config", &sym));

	if(!sym) PetscFunctionReturn(0); // no libblastrampoline visible - nothing to snapshot

#if defined(PETSC_HAVE_DLADDR)
	{
		Dl_info info;

		if(!dladdr(sym, &info) || !info.dli_fname || !strstr(info.dli_fname, "libblastrampoline.5"))
		{
			PetscPrintf(PETSC_COMM_WORLD, "Phase transition plugin  : WARNING: lbt_get_config resolved to an "
				"unrecognised libblastrampoline image (%s); not trusting the struct layout this code was "
				"written against (libblastrampoline 5.15.0). Falling back to "
				"-phase_transition_lbt_ilp64/-phase_transition_lbt_lp64 or LBT_DEFAULT_LIBS.\n",
				(info.dli_fname ? info.dli_fname : "(unknown)"));
			PetscFunctionReturn(0);
		}
	}
#else
	// No dladdr available on this platform: cannot even perform the
	// name-based sanity check above, so do not risk misinterpreting an
	// unknown libblastrampoline layout - skip the snapshot and rely on the
	// options/LBT_DEFAULT_LIBS fallback in PhTrPluginRestoreLbt.
	PetscPrintf(PETSC_COMM_WORLD, "Phase transition plugin  : WARNING: dladdr unavailable on this platform; "
		"cannot verify the libblastrampoline version before trusting its struct layout. Falling back to "
		"-phase_transition_lbt_ilp64/-phase_transition_lbt_lp64 or LBT_DEFAULT_LIBS.\n");
	PetscFunctionReturn(0);
#endif

	{
		const LbtConfig *cfg = ((LbtGetConfigFn)sym)();
		int              i;

		if(!cfg || !cfg->loaded_libs) PetscFunctionReturn(0);

		for(i = 0; i < LBT_MAX_SNAPSHOT && cfg->loaded_libs[i] != NULL; i++)
		{
			PetscCall(PetscStrallocpy(cfg->loaded_libs[i]->libname, &lbtSnap[i].libname));
			PetscCall(PetscStrallocpy(cfg->loaded_libs[i]->suffix,  &lbtSnap[i].suffix));
		}

		lbtSnapCount = i;
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
// Restore PETSc's BLAS/LAPACK forwarding (and thread count) after Julia's
// own LinearAlgebra.__init__ has clobbered it during jl_init_with_image_handle.
//
// Root cause (confirmed): Julia's base runtime includes LinearAlgebra and
// OpenBLAS_jll regardless of whether the plugin code itself uses linear
// algebra. Its __init__ calls libblastrampoline's lbt_forward(libopenblas,
// clear=1, ...), which (a) clears every existing forward and installs its
// own, and (b) resets the BLAS thread count via its own default (typically
// CPU threads / 2, unless OPENBLAS_NUM_THREADS/OMP_NUM_THREADS pin it).
// Because LaMEM and the plugin's bundled libblastrampoline end up being the
// SAME process-wide instance -- confirmed with DYLD_PRINT_LIBRARIES: only
// one libblastrampoline.5.dylib is ever loaded -- this is dyld's ordinary
// install-name de-duplication (both the LaMEM process and the plugin
// bundle reference the library by the same two-level-namespace install
// name, "@rpath/libblastrampoline.5.dylib", so dyld resolves the plugin's
// copy to the identical already-loaded image; this is NOT some special
// property of two-level namespace symbol *resolution*, just ordinary image
// de-duplication by install name), the clobbering is real and was observed
// to break the very next PETSc linear solve ("no BLAS/LAPACK library
// loaded for idamax_()" / "dgemm_()"). It was also observed to silently
// change PETSc's BLAS thread count (LinearAlgebra.__init__ sets it to
	// something based on CPU_THREADS unless OPENBLAS_NUM_THREADS/OMP_NUM_THREADS
	// is set in the environment) -- since the plugin's OpenBLAS IS PETSc's
	// OpenBLAS (same de-duplicated image), this would silently multithread
	// PETSc's own BLAS calls under MPI, competing with MPI ranks for cores.
	//
	// Fix: restore both, from a snapshot taken BEFORE jl_init_with_image_handle
	// (see PhTrPluginSnapshotLbt): re-forward each previously-registered
	// library with clear=0 (additive: does not remove Julia's own
	// registrations, only re-points the specific symbols PETSc needs back
	// to PETSc's chosen library) and check the return value (>0 symbols
	// forwarded); restore the thread count captured before init.
	//
	// -phase_transition_lbt_ilp64/-phase_transition_lbt_lp64 remain as an
	// explicit override for the rare case lbt_get_config is not available
	// (falls back to parsing LBT_DEFAULT_LIBS on ';', matching what PETSc's
	// own libblastrampoline consulted at its own load time) or the
	// snapshot found nothing.
static PetscErrorCode PhTrPluginRestoreLbt(int32_t nthreadsBefore)
{
	void *fwdSym = NULL;

	PetscFunctionBeginUser;

	PetscCall(PetscDLSym(NULL, "lbt_forward", &fwdSym));

	if(!fwdSym)
	{
		PetscPrintf(PETSC_COMM_WORLD, "Phase transition plugin  : WARNING: lbt_forward not found; "
			"cannot verify/restore PETSc's BLAS/LAPACK forwarding after Julia init.\n");
		PetscFunctionReturn(0);
	}

	if(lbtSnapCount > 0)
	{
		int i;

		for(i = 0; i < lbtSnapCount; i++)
		{
			int32_t rc = ((LbtForwardFn)fwdSym)(lbtSnap[i].libname, 0, 0, lbtSnap[i].suffix);

			if(rc <= 0)
			{
				SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB,
					"phase_transition_lib: failed to re-forward BLAS/LAPACK library '%s' "
					"(suffix '%s') via libblastrampoline after Julia init (lbt_forward "
					"returned %d symbols forwarded); PETSc's linear solves would silently "
					"break from here on.", lbtSnap[i].libname, lbtSnap[i].suffix ? lbtSnap[i].suffix : "(none)", (int)rc);
			}

			PetscCall(PetscFree(lbtSnap[i].libname));
			PetscCall(PetscFree(lbtSnap[i].suffix));
		}

		lbtSnapCount = 0;

		PetscPrintf(PETSC_COMM_WORLD, "Phase transition plugin  : re-forwarded %d BLAS/LAPACK "
			"librar%s via libblastrampoline after Julia init\n", (int)i, i == 1 ? "y" : "ies");
	}
	else
	{
		// fallback: no snapshot (lbt_get_config unavailable, or nothing was
		// registered yet when we looked) - try the explicit override
		// options, then LBT_DEFAULT_LIBS itself, parsed on ';'
		char      ilp64lib[_str_len_], lp64lib[_str_len_];
		PetscBool foundIlp64, foundLp64;
		int       nForwarded = 0;

		PetscCall(PetscOptionsGetString(NULL, NULL, "-phase_transition_lbt_ilp64", ilp64lib, _str_len_, &foundIlp64));
		PetscCall(PetscOptionsGetString(NULL, NULL, "-phase_transition_lbt_lp64",  lp64lib,  _str_len_, &foundLp64));

		if(foundIlp64) { if(((LbtForwardFn)fwdSym)(ilp64lib, 0, 0, NULL) > 0) nForwarded++; }
		if(foundLp64)  { if(((LbtForwardFn)fwdSym)(lp64lib,  0, 0, NULL) > 0) nForwarded++; }

		if(!foundIlp64 && !foundLp64)
		{
			const char *envLibs = getenv("LBT_DEFAULT_LIBS");

			if(envLibs)
			{
				char buf[2*_str_len_], *tok, *saveptr;

				PetscCall(PetscStrncpy(buf, envLibs, sizeof(buf)));
				for(tok = strtok_r(buf, ";", &saveptr); tok; tok = strtok_r(NULL, ";", &saveptr))
				{
					if(((LbtForwardFn)fwdSym)(tok, 0, 0, NULL) > 0) nForwarded++;
				}
			}
		}

		if(nForwarded == 0)
		{
			PetscPrintf(PETSC_COMM_WORLD, "Phase transition plugin  : WARNING: could not snapshot PETSc's "
				"BLAS/LAPACK libraries before Julia init, and no fallback (-phase_transition_lbt_ilp64/"
				"-phase_transition_lbt_lp64/LBT_DEFAULT_LIBS) restored any. PETSc's linear solves may "
				"break from here on.\n");
		}
		else
		{
			PetscPrintf(PETSC_COMM_WORLD, "Phase transition plugin  : re-forwarded %d BLAS/LAPACK "
				"librar%s via the fallback path (snapshot was unavailable)\n", nForwarded, nForwarded == 1 ? "y" : "ies");
		}
	}

	// restore the BLAS thread count LinearAlgebra.__init__ silently changed
	{
		void *getSym = NULL, *setSym = NULL;

		PetscCall(PetscDLSym(NULL, "lbt_get_num_threads", &getSym));
		PetscCall(PetscDLSym(NULL, "lbt_set_num_threads", &setSym));

		if(getSym && setSym && nthreadsBefore > 0)
		{
			((LbtSetNumThreadsFn)setSym)(nthreadsBefore);
		}
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
PetscErrorCode PhTrPluginLoad(AdvCtx *actx)
{
	char      lib[_str_len_];
	PetscBool found;
	void     *sym;
	int32_t   nthreadsBefore = -1;

	PetscFunctionBeginUser;

	(void)actx;

	PetscCall(PetscOptionsGetString(NULL, NULL, "-phase_transition_lib", lib, _str_len_, &found));

	if(!found)
	{
		// nothing requested on this call -> built-in transitions only.
		// NOTE: this intentionally does NOT set initTried, so that a LATER
		// call in the same process (e.g. a second LaMEMLibSolve() under an
		// inversion driver) that DOES pass -phase_transition_lib is still
		// honoured, consistent with "only once per process, and only once
		// a plugin is actually requested" rather than "only the first call
		// counts no matter what".
		PetscFunctionReturn(0);
	}

	if(initTried)
	{
		// A plugin was already loaded (or a load was already attempted) in
		// this process. Julia cannot be re-initialised, so a second,
		// different plugin cannot be swapped in - fail loudly rather than
		// silently keep using whatever was loaded first.
		if(active && strcmp(lib, loadedPath) != 0)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP,
				"phase_transition_lib: a different plugin ('%s') was already loaded earlier in this "
				"process ('%s'); Julia cannot be re-initialised with a different plugin. This can "
				"happen under adjoint/inversion drivers that call LaMEMLibSolve() repeatedly.",
				lib, loadedPath);
		}

		PetscFunctionReturn(0); // same path (or a prior attempt failed to load anything) -> no-op
	}

	initTried = PETSC_TRUE;

	PetscCall(PetscStrncpy(loadedPath, lib, _str_len_));

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

	// Tell Julia's runtime, BEFORE jl_init_with_image_handle:
	//   --handle-signals=no : do not install Julia's own signal handlers /
	//     Mach exception ports for SIGSEGV/SIGBUS/SIGFPE etc. This is the
	//     ABI-stable, documented way to do this (jl_parse_opts is an
	//     exported Julia C entry point that predates and postdates the
	//     specific layout of the jl_options struct, so this does not
	//     depend on that struct's layout). A prior version of this file
	//     instead called PetscPushSignalHandler(PetscSignalHandlerDefault,
	//     NULL) AFTER Julia init, believing that would reclaim PETSc's
	//     handlers; this was verified to be a no-op (PETSc's signal.c only
	//     calls signal() when its internal SignalSet flag is false, and it
	//     has been true since PetscInitialize; the push just duplicates a
	//     stack entry without changing what handler is actually installed)
	//     and, on macOS, ineffective for a second reason: Julia
	//     additionally installs Mach exception ports for SIGSEGV/SIGBUS
	//     that a POSIX signal()/sigaction() call cannot displace. See
	//     doc/phase_transition_plugin_PHASE1_REPORT.md for the live
	//     verification (a real SIGSEGV correctly reaches the HOST's
	//     handler, not Julia's, with this fix in place).
	//   --threads=1 --gcthreads=1 : pin the embedded runtime to a single
	//     thread. JULIA_NUM_THREADS/JULIA_NUM_GC_THREADS in the calling
	//     user's environment would otherwise start additional Julia
	//     threads whose GC safepoint mechanism relies on the very signal
	//     handling that --handle-signals=no just disabled - see the ABI
	//     header for the consequences (deep recursion becomes a hard
	//     SEGV rather than a catchable StackOverflowError; no Julia SIGINT
	//     handling; jl_parse_opts itself calls exit() on an unparseable
	//     option, so a typo here would kill the whole LaMEM process).
	{
		static char argv0[] = "lamem";
		static char argv1[] = "--handle-signals=no";
		static char argv2[] = "--threads=1";
		static char argv3[] = "--gcthreads=1";
		static char *jlargv[4] = { argv0, argv1, argv2, argv3 };
		char **jlargvp = jlargv;
		int    jlargc  = 4;

		((JlParseOptsFn)sym)(&jlargc, &jlargvp);
	}

	// Snapshot PETSc's libblastrampoline forwarding table BEFORE Julia's
	// own init has a chance to clobber it (see PhTrPluginSnapshotLbt).
	PetscCall(PhTrPluginSnapshotLbt());
	{
		void *getSym = NULL;
		PetscCall(PetscDLSym(NULL, "lbt_get_num_threads", &getSym));
		if(getSym) nthreadsBefore = ((LbtGetNumThreadsFn)getSym)();
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

	// Restore PETSc's BLAS/LAPACK forwarding and thread count, both
	// silently clobbered by Julia's own LinearAlgebra.__init__ during the
	// call above (see PhTrPluginRestoreLbt for the full explanation).
	PetscCall(PhTrPluginRestoreLbt(nthreadsBefore));

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
	PetscInt     errFlagLoc, errFlagGlob;
	PetscInt     badIdx = -1;
	int32_t      badPhase = 0;

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

	// build SoA input arrays (dimensional) from the local markers. n==0 is
	// handled the same way as n>0 below: every rank still calls fn() and
	// takes part in every collective (MPI_Allreduce for the error flag and
	// for the changed-marker count), it just does so with n=0 - a rank
	// with zero local markers must never skip these collectives, or ranks
	// that DO have markers would deadlock waiting for it.
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

	// call the plugin once for all local markers on this rank (even if
	// n==0). A negative return value is a plugin-signalled failure (e.g.
	// an exception caught on the Julia side) - see phase_transition_plugin.h.
	rc = fn((size_t)n,
		bx, by, bz, bT, bp, (double)time_dim,
		bsxx, bsyy, bszz, bsxy, bsxz, bsyz,
		bj2s, bj2e, beta, baps,
		bphase_in, bphase_out);

	// --- Pass 1: VALIDATE ALL returned phases before writing anything back ---
	// A plugin can fail (rc<0) or return an out-of-range phase on just SOME
	// ranks. Calling SETERRQ directly from inside a per-rank check would
	// make only the failing rank(s) abort while the others carry on into
	// the collectives below (MPI_Allreduce, PetscPrintf) - a classic
	// collective-mismatch hang. Instead: compute a local error flag first,
	// MPI_Allreduce it (MAX) across all ranks, and only THEN have every
	// rank call SETERRQ together if any rank detected a problem - keeping
	// the error path collective, exactly like the success path.
	errFlagLoc = 0;

	if(rc < 0)
	{
		errFlagLoc = 1;
	}
	else
	{
		for(i = 0; i < n; i++)
		{
			if(bphase_out[i] != bphase_in[i] && (bphase_out[i] < 0 || bphase_out[i] >= numPhases))
			{
				errFlagLoc = 1;
				badIdx     = i;
				badPhase   = bphase_out[i];
				break;
			}
		}
	}

	PetscCallMPI(MPI_Allreduce(&errFlagLoc, &errFlagGlob, 1, MPIU_INT, MPI_MAX, PETSC_COMM_WORLD));

	if(errFlagGlob)
	{
		if(rc < 0)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB,
				"phase_transition_lib: plugin reported failure on at least one rank "
				"(this rank's lamem_phase_transition returned %d)", rc);
		}
		else if(badIdx >= 0)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER,
				"phase_transition_lib: plugin returned out-of-range phase on at least one rank "
				"(this rank: phase %d for local marker %" PetscInt_FMT ", valid range: 0..%" PetscInt_FMT ")",
				badPhase, badIdx, numPhases-1);
		}
		else
		{
			// this rank saw no problem itself, but another rank did -
			// still abort collectively rather than silently continuing
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB,
				"phase_transition_lib: another MPI rank reported a plugin failure "
				"(bad return value or out-of-range phase); aborting collectively.");
		}
	}

	// --- Pass 2: only now write back, once every rank is known-good ---
	changed_loc = 0;

	for(i = 0; i < n; i++)
	{
		if(bphase_out[i] != bphase_in[i])
		{
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

// Standalone harness used by the t40_PhaseTransitionPlugin testset to
// verify, against a built plugin (normally ptlib_constant):
//   - lamem_plugin_abi_version() and lamem_plugin_struct_sizes() agree with
//     src/dylib_plugins.h, which this file includes directly (as plain C,
//     LAMEM_PLUGIN_ABI_ONLY), so the header's ABI part is also checked to
//     compile without C++ or PETSc
//   - the "loud failure" guards in LaMEMPlugin.jl's lamem_pt_wrapper: a
//     LaMEMPluginScaling whose length/time/stress are not strictly positive
//     returns -2, one with a different abi_version returns -3, instead of
//     silently reading garbage
// With a second argument "demo" (plugin = ptlib_demo_fields) it instead
// checks the write path of lamem_pt_wrapper on hand-made cell data: the
// rule's APS and pressure writes arrive in internal units (the pressure via
// the exact inverse of (p + pShift)*stress), untouched entries stay
// bit-for-bit as pre-filled.
//
// This mirrors, in a separate process, what the real LaMEM binary does when
// loading a plugin: dlopen (RTLD_NOW|RTLD_GLOBAL), jl_parse_opts with
// --handle-signals=no, jl_init_with_image_handle, then call
// lamem_phase_transition. It must run as its own process rather than being
// exercised via ccall from the Julia test harness itself, because
// initialising a second Julia runtime inside a process that is already one
// (the Julia process running `runtests.jl`) is unsupported.
#include <dlfcn.h>
#include <stdio.h>
#include <string.h>

#define LAMEM_PLUGIN_ABI_ONLY
#include "dylib_plugins.h"

typedef int32_t (*abi_fn)(void);
typedef int32_t (*sizes_fn)(int64_t*, int32_t);

static int same(double a, double b) { return !memcmp(&a, &b, sizeof(double)); }

// ptlib_demo_fields: markers of phase 4 in a cell with Tn >= 1350 C and
// >= 50 % phase 4 get aps = 1 and p *= 0.99 (dimensional), unless aps >= 1
static int run_demo(DylibPluginFn f)
{
	// marker 0 gets seeded, marker 1 already is (left untouched)
	double  x[2] = {0.0, 0.0}, T[2] = {1.2, 1.2}, pin[2] = {0.37, 0.37}, aps[2] = {0.0, 1.0}, o[2] = {0.0, 0.0};
	double  Tout[2], pout[2], apsout[2], w[10][2];
	double  Tn[1] = {1.7}, c0[1] = {0.0}, phr[6] = {0.0, 0.0, 0.0, 0.0, 1.0, 0.0};
	int32_t cell[2] = {0, 0}, phin[2] = {4, 4}, phout[2] = {4, 4}, fs[1] = {0};
	int     i, j;

	for (i = 0; i < 2; i++) { Tout[i] = T[i]; pout[i] = pin[i]; apsout[i] = aps[i]; }
	for (j = 0; j < 10; j++) for (i = 0; i < 2; i++) w[j][i] = 0.0;

	LaMEMPluginMarkers m = {
		2, cell, x, x, x,
		phin, T, pin, aps, o, o, o, o, o, o, o, o, o, o,
		phout, Tout, pout, apsout, w[0], w[1], w[2], w[3], w[4], w[5], w[6], w[7], w[8], w[9]
	};
	LaMEMPluginCells c = {
		1, 6, 0,
		c0, c0, c0, c0, c0, c0,
		c0, c0, c0, c0, Tn, c0, c0, c0, c0, c0, c0,
		c0, c0, c0, c0, c0, c0, c0, c0, c0, fs, c0, c0, c0, c0, c0,
		c0, c0, c0, c0, c0, c0,
		c0, c0, phr
	};
	LaMEMPluginStep s = { 0.0, 0.01, 5 };
	LaMEMPluginScaling sc = {
		DYLIB_PLUGIN_ABI_VERSION, 2,
		1.0, 1.0, 1000.0, 1000.0, 1.0, 1.0, 1.0, 1.0,
		1.0, 1.0, 1.0,
		273.15, 0.05
	};

	int32_t rc = f(&m, &c, &s, &sc);

	// the same floating-point operations lamem_pt_wrapper performs
	double pdim = (pin[0] + sc.pShift)*sc.stress;
	double pnew = (1.0 - 0.01)*pdim;
	double pexp = pnew/sc.stress - sc.pShift;

	int ok_aps   = apsout[0] == 1.0 && same(apsout[1], aps[1]);
	int ok_p     = same(pout[0], pexp) && same(pout[1], pin[1]);
	int ok_other = same(Tout[0], T[0]) && same(Tout[1], T[1]) && phout[0] == 4 && phout[1] == 4;

	for (j = 0; j < 10; j++) for (i = 0; i < 2; i++) ok_other = ok_other && same(w[j][i], 0.0);

	printf("demo: rc=%d, aps %s, p %s (p_out %.17g, expected %.17g), untouched %s\n", (int)rc,
	       ok_aps ? "ok" : "WRONG", ok_p ? "ok" : "WRONG", pout[0], pexp, ok_other ? "ok" : "WRONG");

	return 0;
}

int main(int argc, char **argv)
{
	if (argc < 2) { fprintf(stderr, "usage: %s <plugin.dylib/.so>\n", argv[0]); return 1; }

	void *h = dlopen(argv[1], RTLD_NOW | RTLD_GLOBAL);
	if (!h) { fprintf(stderr, "dlopen: %s\n", dlerror()); return 1; }

	void (*parse)(int*, char***) = (void(*)(int*,char***)) dlsym(h, "jl_parse_opts");
	char a0[] = "lamem", a1[] = "--handle-signals=no";
	char *jlargv[2] = { a0, a1 };
	char **jlargvp = jlargv;
	int   jlargc   = 2;
	if (parse) parse(&jlargc, &jlargvp);

	void (*init)(void*) = (void(*)(void*)) dlsym(h, "jl_init_with_image_handle");
	if (!init) { fprintf(stderr, "no jl_init_with_image_handle\n"); return 2; }
	init(h);

	DylibPluginFn f = (DylibPluginFn) dlsym(h, "lamem_phase_transition");
	if (!f) { fprintf(stderr, "no lamem_phase_transition\n"); return 3; }

	if (argc > 2 && !strcmp(argv[2], "demo")) return run_demo(f);

	abi_fn abiver = (abi_fn) dlsym(h, "lamem_plugin_abi_version");
	printf("abi_version=%d (LaMEM: %d)\n", abiver ? (int)abiver() : -999, DYLIB_PLUGIN_ABI_VERSION);

	sizes_fn sizes = (sizes_fn) dlsym(h, "lamem_plugin_struct_sizes");
	if (sizes)
	{
		int64_t ours[DYLIB_PLUGIN_NUM_STRUCTS] = { (int64_t)sizeof(LaMEMPluginMarkers), (int64_t)sizeof(LaMEMPluginCells),
		                                           (int64_t)sizeof(LaMEMPluginStep),    (int64_t)sizeof(LaMEMPluginScaling) };
		int64_t theirs[DYLIB_PLUGIN_NUM_STRUCTS] = { -1, -1, -1, -1 };
		int32_t nret = sizes(theirs, DYLIB_PLUGIN_NUM_STRUCTS);
		printf("struct sizes: plugin %lld %lld %lld %lld, header %lld %lld %lld %lld -> %s\n",
		       (long long)theirs[0], (long long)theirs[1], (long long)theirs[2], (long long)theirs[3],
		       (long long)ours[0], (long long)ours[1], (long long)ours[2], (long long)ours[3],
		       (nret == DYLIB_PLUGIN_NUM_STRUCTS && !memcmp(ours, theirs, sizeof(ours))) ? "ok" : "MISMATCH");
	}
	else printf("struct sizes: lamem_plugin_struct_sizes missing\n");

	// one marker in cell 0: T_internal=1.5, temperature=1000, Tshift=273.15
	// -> T_dim=1226.85C >= 1200 -> the Constant-transition plugin flips
	// phase 2->3, rc=1
	double  x[1] = {0.0}, T[1] = {1.5}, o[1] = {0.0}, Tout[1] = {1.5}, w[13][1];
	double  phr[6] = {0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
	int32_t cell[1] = {0}, pin[1] = {2}, pout[1] = {2}, fs[1] = {0};
	int     i;

	for (i = 0; i < 13; i++) w[i][0] = 0.0;

	LaMEMPluginMarkers m = {
		1, cell, x, x, x,
		pin, T, o, o, o, o, o, o, o, o, o, o, o, o,
		pout, Tout, w[0], w[1], w[2], w[3], w[4], w[5], w[6], w[7], w[8], w[9], w[10], w[11]
	};
	LaMEMPluginCells c = {
		1, 6, 0,
		o, o, o, o, o, o,
		o, o, o, o, o, o, o, o, o, o, o,
		o, o, o, o, o, o, o, o, o, fs, o, o, o, o, o,
		o, o, o, o, o, o,
		o, o, phr
	};
	LaMEMPluginStep s = { 0.0, 0.01, 5 };
	LaMEMPluginScaling good = {
		DYLIB_PLUGIN_ABI_VERSION, 2,
		1.0, 1.0, 1.0, 1000.0, 1.0, 1.0, 1.0, 1.0,
		1.0, 1.0, 1.0,
		273.15, 0.0
	};
	int32_t c1 = f(&m, &c, &s, &good);
	printf("good struct: rc=%d, phase %d -> %d\n", (int)c1, (int)pin[0], (int)pout[0]);

	// SAME struct but length corrupted to 0 - must trigger the guard (-2)
	LaMEMPluginScaling bad = good;
	bad.length = 0.0;
	printf("bad struct (length=0): rc=%d\n", (int)f(&m, &c, &s, &bad));

	// wrong ABI version in the scaling struct - must trigger -3
	LaMEMPluginScaling old = good;
	old.abi_version = 2;
	printf("bad struct (abi_version=2): rc=%d\n", (int)f(&m, &c, &s, &old));

	return 0;
}

// Standalone harness used by the t40_PhaseTransitionPlugin testset to
// verify the ABI v2 scaling-struct "loud failure" guard implemented in
// LaMEMPlugin.jl's lamem_pt_wrapper (see src/phase_transition_plugin.h for
// the authoritative struct layout, and LaMEMPlugin.jl for the guard
// itself: it rejects a LaMEMPluginScaling whose length/time/stress are not
// strictly positive, returning -2, instead of silently reading garbage).
//
// This mirrors, in a separate process, exactly what the real LaMEM binary
// does when loading a plugin: dlopen (RTLD_NOW|RTLD_GLOBAL), jl_parse_opts
// with --handle-signals=no, jl_init_with_image_handle, then call
// lamem_phase_transition. It must run as its own process rather than being
// exercised via ccall from the Julia test harness itself, because
// initialising a second Julia runtime inside a process that is already one
// (the Julia process running `runtests.jl`) is unsupported - see
// src/phase_transition_plugin.h's "this design cannot be used in-process
// from a Julia host" note.
#include <dlfcn.h>
#include <stdio.h>
#include <stdint.h>

// Mirrors src/phase_transition_plugin.h's `struct LaMEMPluginScaling`
// field-for-field (NOT guessed - copied from that header).
struct LaMEMPluginScaling
{
	int32_t     abi_version;
	int32_t     utype;
	double      length;
	double      time;
	double      stress;
	double      temperature;
	double      viscosity;
	double      strain_rate;
	double      velocity;
	double      density;
	double      Tshift;
	double      pShift;
	double      dt;
	int64_t     step;
	const char *lbl_length;
	const char *lbl_time;
	const char *lbl_stress;
	const char *lbl_temperature;
	const char *lbl_viscosity;
	const char *lbl_strain_rate;
	const char *lbl_velocity;
	const char *lbl_density;
};

typedef int (*pt_fn)(size_t,
	double*, double*, double*, double*, double*, double,
	double*, double*, double*, double*, double*, double*,
	double*, double*, double*, double*,
	int*, int*, double*, const struct LaMEMPluginScaling*);

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

	pt_fn f = (pt_fn) dlsym(h, "lamem_phase_transition");
	if (!f) { fprintf(stderr, "no lamem_phase_transition\n"); return 3; }

	double x[1] = {0.0}, T[1] = {1.5}, o[1] = {0.0};
	int pin[1] = {2}, pout[1];
	double Tout[1];

	// Test 1: a plausible, valid scaling struct.
	struct LaMEMPluginScaling good = {
		2, 2,
		1.0, 1.0, 1.0, 1000.0, 1.0, 1.0, 1.0, 1.0,
		273.15, 0.0, 0.01, 5,
		0, 0, 0, 0, 0, 0, 0, 0
	};
	int c1 = f(1, x, x, x, T, o, 0.0, o, o, o, o, o, o, o, o, o, o, pin, pout, Tout, &good);
	printf("good struct: rc=%d\n", c1);

	// Test 2: the SAME struct but with length corrupted to 0 - must trigger
	// the guard (rc == -2), not silently misbehave.
	struct LaMEMPluginScaling bad = good;
	bad.length = 0.0;
	int c2 = f(1, x, x, x, T, o, 0.0, o, o, o, o, o, o, o, o, o, o, pin, pout, Tout, &bad);
	printf("bad struct (length=0): rc=%d\n", c2);

	return 0;
}

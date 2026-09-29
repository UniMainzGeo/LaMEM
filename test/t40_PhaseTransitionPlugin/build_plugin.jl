# Build the t40 phase-transition plugin bundle (ptlib_constant.jl) with
# juliac, using the JuliaC.jl package. Run this manually (or via CI, once
# juliac is available there) from this directory:
#
#   julia +1.13.0 --startup-file=no --project=<a project with JuliaC installed> build_plugin.jl
#
# It produces build_constant/lib/libptlib_constant.{dylib,so} (plus the
# bundled Julia runtime alongside it), which test/runtests.jl's
# "t40_PhaseTransitionPlugin" testset looks for and skips (with a clear
# message) if not present -- this lets the testset run wherever a
# juliac-capable Julia is available while not requiring one in ordinary CI
# runs that only build LaMEM's C/C++ code.
#
# Built and verified with Julia 1.13.0 (juliaup channel "+1.13.0") against
# JuliaC.jl as installed in a scratch project; any reasonably recent Julia
# with a working `juliac`/JuliaC.jl setup should work, but 1.13.0 is what
# this file's own bundle was actually produced and tested with (see
# doc/phase_transition_plugin_PHASE1_REPORT.md).
using JuliaC

JuliaC.main([
    "--output-lib", "libptlib_constant.dylib", # ".so" is substituted automatically on Linux
    "--bundle", "build_constant",
    "--trim=safe",
    "--compile-ccallable",
    "--experimental",
    "ptlib_constant.jl",
])

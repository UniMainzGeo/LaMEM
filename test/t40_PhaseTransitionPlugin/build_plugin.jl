# Build a t40 phase-transition plugin bundle with juliac (JuliaC.jl). Run:
#   julia +1.13.0 --startup-file=no --project=<env with JuliaC> build_plugin.jl [src.jl] [name]
# src.jl defaults to ptlib_constant.jl, name to its stem; produces
# build_<name>/lib/libptlib_<name>.{dylib,so}. --output-lib takes NO
# extension (JuliaC appends the platform's own dlext; a wrong explicit one
# is a hard error in JuliaC's link_products, not silently substituted).
using JuliaC

src  = length(ARGS) >= 1 ? ARGS[1] : "ptlib_constant.jl"
name = length(ARGS) >= 2 ? ARGS[2] : replace(splitext(src)[1], "ptlib_" => "")

JuliaC.main([
    "--output-lib", "libptlib_$name",
    "--bundle", "build_$name",
    "--project", dirname(abspath(src)),
    "--trim=safe",
    "--compile-ccallable",
    "--experimental",
    src,
])

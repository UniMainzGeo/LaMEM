# Build a dylib_plugin bundle with juliac (JuliaC.jl). Run from anywhere:
#   julia --project=@juliac scripts/dylib_plugins/build_plugin.jl <rule.jl> [output_dir]
# output_dir defaults to next to rule.jl; produces
# <output_dir>/build_<name>/lib/libptlib_<name>.{dylib,so}, name = rule's
# stem with a leading "ptlib_" stripped. --output-lib takes NO extension
# (JuliaC appends the platform's own dlext; a wrong explicit one is a hard
# error in JuliaC's link_products, not silently substituted). No --project
# is passed to JuliaC.main: rule files have no package deps of their own, so
# JuliaC falls back to the caller's own active project (@juliac) instead of
# instantiating a separate one - this is also why no Project.toml/Manifest.toml
# ever appears next to a rule file or in output_dir.
using JuliaC

length(ARGS) >= 1 || error("usage: build_plugin.jl <rule.jl> [output_dir]")

src        = abspath(ARGS[1])
name       = replace(splitext(basename(src))[1], "ptlib_" => "")
output_dir = length(ARGS) >= 2 ? abspath(ARGS[2]) : dirname(src)

JuliaC.main([
    "--output-lib", "libptlib_$name",
    "--bundle", joinpath(output_dir, "build_$name"),
    "--trim=safe",
    "--compile-ccallable",
    "--experimental",
    src,
])

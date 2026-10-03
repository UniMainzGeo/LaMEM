# this compiles LaMEM using the PETSc_jll libraries
using PETSc_jll
using MPICH_jll

cd("../src")

# read command-line arguments
is64bit  = any(contains.(ARGS, "int64"))
do_check = any(contains.(ARGS, "check"))

# Where setup_packages.jl staged PETSc. The Makefile reads PETSC_OPT/PETSC_DEB/
# FASTSCAPE_LIB as POSIX paths. On Windows the real destdir is a native path, but
# the CI job mounts it at /workspace/destdir inside MSYS2, so that is what make sees.
destdir       = get(ENV, "LAMEM_CI_DESTDIR", "/workspace/destdir")
petsc_destdir = Sys.iswindows() ? "/workspace/destdir" : destdir
pjoin(parts...) = join(parts, "/") # not joinpath: must stay POSIX on Windows too

# the Int32 build has no separate debug variant
opt, deb = is64bit ? ("double_real_Int64", "double_real_Int64_deb") :
                     ("double_real_Int32", "double_real_Int32")
println("Using PETSc that has $(is64bit ? 64 : 32)bit integers")

# Take the environment (dynamic libraries etc.) from the PETSc
cmd = addenv(PETSc_jll.ex42(),
             "PETSC_OPT"     => pjoin(petsc_destdir, "lib", "petsc", opt),
             "PETSC_DEB"     => pjoin(petsc_destdir, "lib", "petsc", deb),
             "FASTSCAPE_LIB" => pjoin(petsc_destdir, "lib", "fastscape"))

# PETSc_jll.ex42() carries its own baked-in env and does not inherit ours, so
# anything make's child processes need from the calling shell must be forwarded.
if Sys.iswindows()
    # make itself is found via the real shell PATH, but its children inherit the
    # Cmd's env, which lacks MSYS2's /mingw64/bin where c++ lives
    i = findfirst(startswith("PATH="), cmd.env)
    jll_path = isnothing(i) ? "" : cmd.env[i][6:end]
    cmd = addenv(cmd, "PATH" => ENV["PATH"] * ";" * jll_path)
end
if haskey(ENV, "LIBRARY_PATH") # e.g. macOS CI pointing the linker at libemutls_w
    cmd = addenv(cmd, "LIBRARY_PATH" => ENV["LIBRARY_PATH"])
end

@show pkgversion(PETSc_jll)
#@show pkgversion(MPICH_jll)

make(args) = run(Cmd(`make $args`, env = cmd.env))

if do_check
    # only check source formatting, don't compile
    println("---- Checking LaMEM source formatting ----")
    make(`mode=opt check`)
else
    println("Compiling LaMEM")

    println("---- Compiling LaMEM opt version ----")
    make(`mode=opt surf=scape all`)

    println("---- Compiling LaMEM deb version ----")
    make(`mode=deb surf=scape all`)
end

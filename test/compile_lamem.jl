# this compiles LaMEM using the PETSc_jll libraries
using PETSc_jll
using MPICH_jll

cd("../src")

# read command-line arguments
is64bit  = any(contains.(ARGS, "int64"))
do_check = any(contains.(ARGS, "check"))

# see setup_packages.jl: on Linux CI this is the fixed /workspace/destdir; on other
# platforms (no such writable root-level directory) it is a user-writable one, given
# as a plain native path (e.g. "C:\lamem_destdir" on Windows - the same form Julia's
# own filesystem calls use in setup_packages.jl).
destdir = get(ENV, "LAMEM_CI_DESTDIR", "/workspace/destdir")

# src/Makefile reads PETSC_OPT/PETSC_DEB/FASTSCAPE_LIB itself (`include
# ${PETSC_DIR}/lib/petsc/conf/variables` and friends), which on Windows runs under
# MSYS2's `make` and needs a POSIX-style path - but PETSc_jll's own
# lib/petsc/.../conf/variables has a further `include
# /workspace/destdir/.../petscvariables` line baked in absolutely at Yggdrasil
# build time (BinaryBuilder's own sandbox path, not relocatable). The Windows CI
# job's "Mount /workspace/destdir..." step gives MSYS2 an /etc/fstab entry so that
# POSIX path resolves to the real destdir - so PETSC_OPT/PETSC_DEB/FASTSCAPE_LIB
# must use literal /workspace/destdir on Windows too (matching Linux, where
# destdir already *is* /workspace/destdir), not a drive-letter-derived path.
function msys2_path(p::AbstractString)
    Sys.iswindows() || return p
    p == "/workspace/destdir" && return p
    p = replace(p, "\\" => "/")
    m = match(r"^([A-Za-z]):(.*)$", p)
    isnothing(m) && return p
    return "/" * lowercase(m.captures[1]) * m.captures[2]
end
psep(parts...) = msys2_path(join(parts, "/"))
petsc_destdir = Sys.iswindows() ? "/workspace/destdir" : destdir

# Take the environment (dynamic libraries etc.) from the PETSc
if is64bit
    println("Using PETSc that has 64bit integers")
    cmd = addenv(PETSc_jll.ex42(),
                    "PETSC_OPT"=>psep(petsc_destdir, "lib", "petsc", "double_real_Int64"),
                    "PETSC_DEB"=>psep(petsc_destdir, "lib", "petsc", "double_real_Int64_deb"),
                )

else
    println("Using PETSc that has 32bit integers")
    cmd = addenv(PETSc_jll.ex42(),
                    "PETSC_OPT"=>psep(petsc_destdir, "lib", "petsc", "double_real_Int32"),
                    "PETSC_DEB"=>psep(petsc_destdir, "lib", "petsc", "double_real_Int32"),
                )
end

cmd = addenv(cmd, "FASTSCAPE_LIB"=>psep(petsc_destdir, "lib", "fastscape"))

@show pkgversion(PETSc_jll)
#@show pkgversion(MPICH_jll)

# src/Makefile is a POSIX Makefile (uname, readlink, $(shell ...), etc.). On Linux it
# runs directly; on Windows this script itself must be invoked from inside an MSYS2
# shell (see the CI workflow's Windows job), which provides a real `make` and the
# rest of the POSIX toolchain the Makefile needs - Julia's own `run` still just
# execs the `make` found on PATH, whichever one that is.
if do_check
    # only check source formatting, don't compile
    println("---- Checking LaMEM source formatting ----")
    check_format = Cmd(`make mode=opt check`, env = cmd.env)
    run(check_format)
else
    println("Compiling LaMEM")

    println("---- Compiling LaMEM opt version ----")
    compile_lamem = Cmd(`make mode=opt surf=scape all`, env = cmd.env)
    run(compile_lamem)

    println("---- Compiling LaMEM deb version ----")
    compile_lamem = Cmd(`make mode=deb surf=scape all`, env = cmd.env)
    run(compile_lamem)


end

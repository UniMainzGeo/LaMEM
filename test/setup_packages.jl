# this downloads the required packages

# Add PETSc with required version
using Pkg
Pkg.add(name="PETSc_jll", version="3.25.4")
Pkg.add(name="MPICH_jll", version="5.0.1")
Pkg.add(name="Fastscapelib_jll")

# Copy the relevant directories over
using PETSc_jll, MPICH_jll, Fastscapelib_jll

# Destination prefix that compile_lamem.jl's PETSC_OPT/PETSC_DEB/FASTSCAPE_LIB point
# at. Linux CI relies on the fixed /workspace/destdir path (root-owned, needs sudo,
# and is what the real Yggdrasil recipe's cross-compilation sandbox also uses); on
# other platforms there is no such fixed, writable root-level directory, so an
# ordinary user-writable one is used instead and sudo is skipped, since it does not
# exist on Windows runners and is not needed for a directory the user already owns.
destdir = get(ENV, "LAMEM_CI_DESTDIR", "/workspace/destdir")
use_sudo = !Sys.iswindows() && destdir == "/workspace/destdir"

maybe_sudo(cmd::Cmd) = use_sudo ? `sudo -E $cmd` : cmd

mkpath(destdir)

# copy the contents of all directories in a single one
for path in PETSc_jll.PATH_list
    cur_dir = path[1:end-3]

    # copy mpi directories - we somehow have to do that one by one
    dirs = ["bin", "lib", "include", "share"]
    for d in dirs
        if isdir(joinpath(cur_dir, d))
            run(maybe_sudo(`cp -rf $(joinpath(cur_dir, d)) $destdir/`))
        end
    end

    # On Windows, PETSc_jll nests each precision/index variant's own tree
    # (lib/petsc/conf/variables and friends, read directly by LaMEM's
    # Makefile) under bin/petsc/<variant>/ instead of lib/petsc/<variant>/
    # like on Unix - so the generic copy above lands it at
    # destdir/bin/petsc/... instead of the destdir/lib/petsc/... path
    # compile_lamem.jl's PETSC_OPT/PETSC_DEB point at. Mirror it there too.
    if Sys.iswindows() && isdir(joinpath(cur_dir, "bin", "petsc"))
        mkpath(joinpath(destdir, "lib"))
        run(maybe_sudo(`cp -rf $(joinpath(cur_dir, "bin", "petsc")) $(joinpath(destdir, "lib"))/`))
    end
end

"""
    copy all files
"""
function cp_files(srcdir, destdir; force=true)
    for f in readdir(srcdir)
        if isfile(joinpath(srcdir, f))
            src = joinpath(srcdir, f)
            dst = joinpath(destdir, f)
            #cp(src, dst, force=force)
            run(maybe_sudo(`cp -rf $src $dst`))

        end
    end
    return nothing
end

# And all required dynamic libraries (except petsc)
for srcdir in PETSc_jll.LIBPATH_list
    if !contains(srcdir, "petsc")
        dest_dir = joinpath(destdir, "lib")
        cp_files(srcdir, dest_dir)
    end
end

# copy PETSc directories
#run(`sudo -E cp -rf $petsc_dir/lib /workspace/destdir`)

# print
run(`ls $(joinpath(destdir, "lib"))`)

# Stage the FastScape library where the LaMEM Makefile expects it (FASTSCAPE_LIB)
fastscape_dir = joinpath(destdir, "lib", "fastscape")
run(maybe_sudo(`mkdir -p $fastscape_dir`))
for srcdir in Fastscapelib_jll.LIBPATH_list
    cp_files(srcdir, fastscape_dir)
end
run(`ls $fastscape_dir`)

# PETSc_jll's own lib/petsc/<variant>/lib/petsc/conf/variables has a further
# `include /workspace/destdir/.../petscvariables` line baked in absolutely at
# Yggdrasil build time (BinaryBuilder's own cross-compilation sandbox path,
# not relocatable). On Linux this is harmless, since destdir there already
# IS /workspace/destdir. On the GH-hosted macOS runner, no such fixed,
# writable root-level directory exists at all - even creating one via sudo
# mkdir/ln -s fails there with "Read-only file system" - so that absolute
# include can never resolve on its own; rewrite it, in our own copy of the
# file only, to point at wherever destdir actually is. Windows is handled
# differently (an MSYS2 /etc/fstab mount makes the literal POSIX path
# /workspace/destdir itself resolve to the real destdir - see the CI
# workflow's Windows job), so it must keep the original /workspace/destdir
# string in this file, not get it rewritten to a native Windows path here.
if !Sys.iswindows() && destdir != "/workspace/destdir"
    petsc_lib_dir = joinpath(destdir, "lib", "petsc")
    if isdir(petsc_lib_dir)
        for variant in readdir(petsc_lib_dir)
            variables_file = joinpath(petsc_lib_dir, variant, "lib", "petsc", "conf", "variables")
            if isfile(variables_file)
                contents = read(variables_file, String)
                fixed = replace(contents, "/workspace/destdir" => destdir)
                if fixed != contents
                    # cp -rf preserves the source Julia artifact's read-only
                    # permissions (artifacts are stored read-only on disk),
                    # so our own copy under destdir is read-only too - make
                    # it writable before overwriting its contents in place.
                    chmod(variables_file, 0o644)
                    write(variables_file, fixed)
                end
            end
        end
    end
end

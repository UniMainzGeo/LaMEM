# this downloads the required packages

# Add PETSc
using Pkg
Pkg.add(name="PETSc_jll", version="3.18.6")

# Copy the relevant directories over
using PETSc_jll

mpi_dir   = PETSc_jll.PATH_list[1][1:end-3]
petsc_dir = PETSc_jll.PATH_list[2][1:end-3]
@show mpi_dir
@show petsc_dir

# copy mpi directories - we somehow have to do that one by one
dirs = ["bin","lib","include","share"]
for d in dirs
    run(`sudo -E cp -r $mpi_dir/$d /workspace/destdir/`)
end

# Same with petsc
dirs = ["bin","lib","share"]
for d in dirs
    run(`sudo -E cp -r $mpi_dir/$d /workspace/destdir/`)
end


"""
    copy all files 
"""
function cp_files(srcdir, destdir; force=true)
    for f in readdir(srcdir)
        if isfile(joinpath(srcdir,f))
            src = joinpath(srcdir,f)
            dst = joinpath(destdir,f)
            #cp(src, dst, force=force)
            run(`sudo -E cp -r $src $dst`)

        end
    end
    return nothing
end

# And all required dynamic libraries (except petsc)
for srcdir in PETSc_jll.LIBPATH_list
    if !contains(srcdir, "/julia") && !contains(srcdir,"petsc")
        #run(`sudo -E cp -r $d/"*".so /workspace/destdir/lib/`)
        #run(`sudo -E cp -r $d/ /workspace/destdir/lib/`)
        
        #dest_dir = "/Users/kausb/Downloads/workspace/destdir/lib"
        dest_dir = "/workspace/destdir/lib/"

        cp_files(srcdir, dest_dir)
    end
end

# copy PETSc directories
run(`sudo -E cp -rf $petsc_dir/lib /workspace/destdir`)

# print
run(`ls /workspace/destdir/lib`);

# Compile LaMEM
println("Compiling LaMEM")
cd("../src")
@show pwd()

# Take the environment (dynamic libraries etc.) from the PETSc
cmd = addenv(PETSc_jll.ex42(), "PETSC_OPT"=>"/workspace/destdir/lib/petsc/double_real_Int64")
compile_lamem = Cmd(`make mode=opt all`, env = cmd.env)
run(compile_lamem)

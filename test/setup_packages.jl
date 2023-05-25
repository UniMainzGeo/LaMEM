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

# And all required dynamic libraries (except petsc)
for d in PETSc_jll.LIBPATH_list
    if !contains(d, "/julia") && !contains(d,"petsc")
        run(`sudo -E cp -r $d /workspace/destdir/lib`)
    end
end

# print
run(`ls /workspace/destdir/lib`);

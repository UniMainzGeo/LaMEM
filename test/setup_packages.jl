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

# print content of directories
# copy mpi directories - we somehow have to do that one by one
dirs = ["bin","lib","include","share"]
for d in dirs
    run(`sudo -E cp -r $mpi_dir/$d /workspace/destdir/`)
end

run(`ls /workspace/destdir/`);

# print
run(`ls /workspace/destdir/lib`);

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
run(`ls $(mpi_dir)`);

# copy mpi directories
#run(`cp -r $mpi_dir /Users/kausb/Downloads/workspace/destdir/`)

mpi_dir = mpi_dir*"*"
run(`sudo -E cp -r $mpi_dir/bin /workspace/destdir/`)

run(`ls /workspace/destdir/`);

# print
run(`ls /workspace/destdir/lib`);

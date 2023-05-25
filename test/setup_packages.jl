# this downloads the required packages

# Add PETSc
using Pkg
Pkg.add("PETSc_jll", version="3.18.6")

# Copy the relevant directories over
using PETSc_jll

mpi_dir   = PETSc_jll.PATH_list[1][1:end-3]
petsc_dir = PETSc_jll.PATH_list[2][1:end-3]

# copy mpi directories
run(`sudo -E cp -r $mpi_dir /workspace/destdir/`)

# print
run(`ls /workspace/destdir/`);

# print
run(`ls /workspace/destdir/lib`);

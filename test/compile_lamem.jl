# this compiles LaMEM using the PETSc_jll libraries
using PETSc_jll
using MPICH_jll

# Compile LaMEM
println("Compiling LaMEM")
cd("../src")

# read command-line argument
if any(contains.(ARGS,"int64"))
    is64bit = true
else
    is64bit = false
end

# Take the environment (dynamic libraries etc.) from the PETSc
if is64bit
    println("Compiling LaMEM with PETSc that has 64bit integers")
    cmd = addenv(PETSc_jll.ex42(), 
                    "PETSC_OPT"=>"/workspace/destdir/lib/petsc/double_real_Int64",
                    "PETSC_DEB"=>"/workspace/destdir/lib/petsc/double_real_Int64_deb",
                )
    
else
    println("Compiling LaMEM with PETSc that has 32bit integers")
    cmd = addenv(PETSc_jll.ex42(), 
                    "PETSC_OPT"=>"/workspace/destdir/lib/petsc/double_real_Int32",
                    "PETSC_DEB"=>"/workspace/destdir/lib/petsc/double_real_Int32",
                )
end

@show pkgversion(PETSc_jll)
#@show pkgversion(MPICH_jll)

println("---- Compiling LaMEM opt version ----")
compile_lamem = Cmd(`make mode=opt all`, env = cmd.env)
run(compile_lamem)

println("---- Compiling LaMEM deb version ----")
compile_lamem = Cmd(`make mode=deb all`, env = cmd.env)
run(compile_lamem)




# this compiles LaMEM using the PETSc_jll libraries
using PETSc_jll
using MPICH_jll

cd("../src")

# read command-line arguments
is64bit  = any(contains.(ARGS, "int64"))
do_check = any(contains.(ARGS, "check"))

# Take the environment (dynamic libraries etc.) from the PETSc
if is64bit
    println("Using PETSc that has 64bit integers")
    cmd = addenv(PETSc_jll.ex42(), 
                    "PETSC_OPT"=>"/workspace/destdir/lib/petsc/double_real_Int64",
                    "PETSC_DEB"=>"/workspace/destdir/lib/petsc/double_real_Int64_deb",
                )
    
else
    println("Using PETSc that has 32bit integers")
    cmd = addenv(PETSc_jll.ex42(), 
                    "PETSC_OPT"=>"/workspace/destdir/lib/petsc/double_real_Int32",
                    "PETSC_DEB"=>"/workspace/destdir/lib/petsc/double_real_Int32",
                )
end

@show pkgversion(PETSc_jll)
#@show pkgversion(MPICH_jll)

if do_check
    # only check source formatting, don't compile
    println("---- Checking LaMEM source formatting ----")
    check_format = Cmd(`make mode=opt check`, env = cmd.env)
    run(check_format)
else
    println("Compiling LaMEM")

    println("---- Compiling LaMEM opt version ----")
    compile_lamem = Cmd(`make mode=opt all`, env = cmd.env)
    run(compile_lamem)

    println("---- Compiling LaMEM deb version ----")
    compile_lamem = Cmd(`make mode=deb all`, env = cmd.env)
    run(compile_lamem)
end

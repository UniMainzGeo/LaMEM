# this starts the testing suite from the CI
# the difference is that we use LaMEM compiled vs dynamic libraries so we need to add the environment
using Pkg,LaMEM_C, Test


if "is64bit" in ARGS
    args_local = ["use_dynamic_lib","is64bit","test_superlu"]
else
    # 32bit PETSc installation
    args_local = ["use_dynamic_lib"]
end

# run test suite
Pkg.test("LaMEM_C", test_args=args_local)


exit()

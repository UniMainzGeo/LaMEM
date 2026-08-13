# this starts the testing suite from the CI
# the difference is that we use LaMEM compiled vs dynamic libraries so we need to add the environment
using Pkg,LaMEM_C, Test

if "is64bit" in ARGS
    args_local = ["use_dynamic_lib","is64bit"]
else
    # 32bit PETSc installation
    args_local = ["use_dynamic_lib"]
end

# pass through test selectors (bare numbers or ranges, e.g. "37" or "03-07")
# so that the CI can run only a subset of the testsuite
for a in ARGS
    if occursin(r"^\d+(-\d+)?$", a)
        args_local = push!(args_local, a)
    end
end

# run test suite
Pkg.test("LaMEM_C", test_args=args_local)

exit()

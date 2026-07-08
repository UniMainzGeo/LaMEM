using Pkg

args_local = String[]
if "is64bit" in ARGS
    push!(args_local, "is64bit")
end
if "valgrind" in ARGS
    push!(args_local, "valgrind")
end
if isempty(args_local)
    # 32bit PETSc installation, no valgrind
    args_local = [""]
end

# compile LaMEM if required
cur_dir = pwd()
cd("../src")
run(`make mode=opt all`)
run(`make mode=deb all`)
cd(cur_dir)

# run test suite
Pkg.test("LaMEM_C", test_args=args_local)

exit()

using Pkg

if "is64bit" in ARGS
    args_local = ["is64bit",]
else
    # 32bit PETSc installation
    args_local = ["",]
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

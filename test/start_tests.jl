using Pkg

# compile LaMEM if required
cur_dir = pwd()
cd("../src")
run(`make mode=opt all`)
run(`make mode=deb all`)
cd(cur_dir)

# run test suite
Pkg.test("LaMEM_C")

exit()
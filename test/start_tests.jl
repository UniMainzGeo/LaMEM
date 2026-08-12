using Pkg

args_local = String[]
if "is64bit" in ARGS
    push!(args_local, "is64bit")
end
if "valgrind" in ARGS
    push!(args_local, "valgrind")
end
if "check" in ARGS
    push!(args_local, "check")
end

# Any remaining arguments are treated as test selectors (e.g. "01", "05",
# "03-07") and forwarded so a subset of tests can be run, e.g.:
#   make test 01 05 32
#   make test 03-07 11 12-17

known_flags = ("is64bit", "valgrind", "check", "create_plots")

test_selectors = [a for a in ARGS if !(a in known_flags)]

append!(args_local, test_selectors)

# compile LaMEM if required
cur_dir = pwd()
cd("../src")
run(`make mode=opt all`)
run(`make mode=deb all`)
cd(cur_dir)

# run test suite
Pkg.test("LaMEM_C", test_args=args_local)

exit()

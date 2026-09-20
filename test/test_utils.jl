# These are tools that help perform the LaMEM tests, which run LaMEM locally
using LinearAlgebra, Glob
#import LaMEM.Run: deactivate_multithreading
#include("julia/IO_functions.jl")

if use_dynamic_lib
    #using LaMEM.LaMEM_jll.PETSc_jll
    using PETSc_jll
    using OpenBLAS_jll   # Julia's bundled ILP64 OpenBLAS, forwarded to PETSc via LBT (see add_dylibs)
end

export run_lamem_local_test, perform_lamem_test, clean_test_directory, run_lamem_save_grid_local, mpiexec
export CreatePartitioningFile_local, LaMEM_has_fastscape


if use_dynamic_lib
    mpiexec = if PETSc_jll.MPICH_jll.is_available()
        PETSc_jll.MPICH_jll.mpiexec()
    elseif PETSc_jll.MicrosoftMPI_jll.is_available()
        PETSc_jll.MicrosoftMPI_jll.mpiexec()
    elseif PETSc_jll.OpenMPI_jll.is_available()
        PETSc_jll.OpenMPI_jll.mpiexec()
    else
        warning("")
        nothing
    end

else
    mpiexec = "mpiexec"
end



"""
    run_lamem_local_test(ParamFile::String, cores::Int64=1, args::String=""; 
                        outfile="test.out", bin_dir="../../bin", deb=false,
                        mpiexec="mpiexec", valgrind=false)

This runs a LaMEM simulation with given `ParamFile` on 1 or more cores, while writing the output to a local log file.

If `valgrind=true`, the run is instead performed under Valgrind (memcheck), invoked directly - for
cores==1 as `valgrind ... exec`, for cores>1 as `mpiexec -n cores valgrind ... exec`, using the same
`mpiexec` as a normal parallel run - and always uses the debug (`deb`) build regardless of the
`deb` kwarg passed in, since Valgrind needs debug info to produce useful reports. Valgrind's XML
report(s) are written next to the log file as `outfile_<pid>.xml`, and its combined stdout/stderr as
`outfile.out`.

"""
function run_lamem_local_test(ParamFile::String, cores::Int64=1, args::String=""; 
                outfile="test.out", bin_dir="../../bin", deb=false,
                mpiexec="mpiexec", valgrind::Bool=false)
    
    cur_dir = pwd()
    if valgrind
        # Valgrind needs debug info (and is far more useful against an unoptimized build) - always use
        # the debug binary here, regardless of which one the test itself asked for.
        exec = joinpath(cur_dir,bin_dir,"deb","LaMEM")
    elseif deb
        exec=joinpath(cur_dir,bin_dir,"deb","LaMEM")
    else
        exec=joinpath(cur_dir,bin_dir,"opt","LaMEM")
    end

    success = true
    dylibs, mpipath = get_dylibs()
    args = split(args)

    if valgrind
        # Run LaMEM directly under Valgrind, using the same mpiexec as a normal parallel run (for
        # cores==1 we skip mpiexec entirely, just like the non-valgrind path below does).
        # Valgrind's XML report is written per-rank as outbase_<pid>.xml; combined stdout/stderr as
        # outbase.out.
        outbase = isempty(outfile) ? "valgrind_out" : first(splitext(outfile))

        valgrind_cmd = `valgrind -v --leak-check=full --track-origins=yes --show-reachable=yes --xml=yes --xml-file=$(outbase)_%p.xml --child-silent-after-fork=yes -q`

        if cores == 1
            perform_run = Cmd(`$(valgrind_cmd) $(exec) -ParamFile $(ParamFile) $args`)
        else
            perform_run = Cmd(`$(mpiexec) -n $(cores) $(valgrind_cmd) $(exec) -ParamFile $(ParamFile) $args`)
        end

        perform_run = addenv(add_dylibs(perform_run, dylibs), "MPIWRAP_DEBUG"=>"quiet")

        try
            # Open once and share the same IOStream for stdout & stderr, so they interleave into a
            # single file correctly (matching bash's `> file 2>&1`) instead of two independent writers.
            open("$(outbase).out", "w") do io
                run(pipeline(perform_run, stdout=io, stderr=io))
            end
        catch
            println("An error occured running Valgrind in directory: $(cur_dir) ")
            println("while running the command:")
            println(perform_run)
            success = false
        end

        return success
    end

    try
        if cores==1
            perform_run = Cmd(`$(exec) -ParamFile $(ParamFile) $args`);
            
            # add dynamic libraries to the path (if specified)
            perform_run = add_dylibs(perform_run, dylibs)

           ## perform_run = deactivate_multithreading(perform_run)

            #perform_run = addenv(perform_run,"PATH"=>mpipath)

            # Run LaMEM on a single core, which does not require a working MPI
            try 
                if !isempty(outfile)
                    run(pipeline(perform_run, stdout=outfile));
                else
                    run(perform_run);
                end
            catch
                println("An error occured in directory: $(cur_dir) ")
                println("while running the script:")
                println(perform_run)
                success = false
            end
        else

            perform_run = Cmd(`$(mpiexec) -n $(cores) $(exec) -ParamFile $(ParamFile) $args`);

            # add dynamic libraries to the path (if specified)
            perform_run = add_dylibs(perform_run, dylibs)

       ##     perform_run = deactivate_multithreading(perform_run)

           # perform_run = addenv(perform_run,"PATH"=>mpipath)
  
            # set correct environment
            #mpirun = setenv(mpiexec, LaMEM_jll.JLLWrappers.JLLWrappers.LIBPATH_env=>LaMEM_jll.LIBPATH[]);
            # Run LaMEM in parallel
            try 
                if !isempty(outfile)
                    run(pipeline(perform_run, stdout=outfile));
                else
                    run(perform_run);
                end
            catch
                println(perform_run)
                success = false
            end
        end
    catch
        success = false
    end
  
    return success
end

"""
    LaMEM_has_fastscape(; bin_dir="../bin", deb=false)

Checks whether the LaMEM binary currently installed in `bin_dir` was compiled
with FastScape support (`make surf=scape`), by running it with the
`-fastscape_info` flag and inspecting its output. Used to skip FastScape-only
tests when running against a binary that wasn't built with FastScape enabled.
"""
function LaMEM_has_fastscape(; bin_dir="../bin", deb=false)

    cur_dir = pwd()
    if deb
        exec = joinpath(cur_dir, bin_dir, "deb", "LaMEM")
    else
        exec = joinpath(cur_dir, bin_dir, "opt", "LaMEM")
    end

    if !isfile(exec)
        return false
    end

    dylibs, _ = get_dylibs()
    has_fastscape = false
    try
        perform_run = add_dylibs(Cmd(`$(exec) -fastscape_info`), dylibs)
        out = read(perform_run, String)
        has_fastscape = contains(out, "FASTSCAPE_ENABLED")
    catch
        has_fastscape = false
    end

    return has_fastscape
end

function get_line_containing(stringarray::Vector{SubString{String}}, lookfor::String)

	for line in stringarray
		   if contains(line, lookfor)
		   foundline=line
		   return foundline
		   end
	end
end

# Matches a single "Object Type   Creations   Destructions [  Memory  Descendants' Mem.]"
# row from PETSc's `-log_view` memory-usage table. The trailing Memory/Descendants' Mem.
# columns are only present for some PETSc builds (e.g. debug builds track allocation sizes;
# optimized builds often only track Creations/Destructions counts), so they are optional:
#   "              Vector   258            258     18944816     0."   (deb-style, with sizes)
#   "              Vector   140            140"                       (opt-style, counts only)
# The object-type name is matched non-greedily so it can contain spaces.
const MEMORY_ROW_RE = r"^\s*([A-Za-z][A-Za-z0-9 /'\-]*?)\s+(\d+)\s+(\d+)(?:\s+[\d.]+\s+[\d.]+\.?)?\s*$"

"""
    counts = parse_memory_usage(file::String)

Parses a `-log_view` log `file` and returns a `Dict{String,Tuple{Int,Int}}` mapping each
PETSc object type (Vector, Matrix, Index Set, ...) to its summed `(creations,
destructions)` across the whole file (all event stages combined). Returns an empty `Dict`
if the file contains no `-log_view` object-tracking table.

The table's header line ("Object Type   Creations   Destructions ...") is used as the
trigger to start parsing, rather than a fixed marker string, since its exact wording
(and whether Memory/Descendants' Mem. columns are present) varies between PETSc builds.
"""
function parse_memory_usage(file::String)
    counts = Dict{String,Tuple{Int,Int}}()
    in_section = false

    open(file) do io
        for line in eachline(io)
            if !in_section
                if occursin("Creations", line) && occursin("Destructions", line)
                    in_section = true
                end
                continue
            end

            m = match(MEMORY_ROW_RE, line)
            m === nothing && continue

            name = strip(m.captures[1])
            creations    = parse(Int, m.captures[2])
            destructions = parse(Int, m.captures[3])

            prev = get(counts, name, (0, 0))
            counts[name] = (prev[1] + creations, prev[2] + destructions)
        end
    end

    return counts
end

"""
    success = check_memory_usage(file::String)

Checks a `-log_view` log `file` (see `parse_memory_usage`) and reports whether every
PETSc object type was destroyed as many times as it was created. On a mismatch (a likely
missing `*Destroy()` call) this prints the offending object types and returns `false`.
If the file contains no `-log_view` table at all, a warning is printed and `false`
is returned. Otherwise returns `true`.
"""
function check_memory_usage(file::String)
    counts = parse_memory_usage(file)

    if isempty(counts)
        println("WARNING: no -log_view memory usage table found in $file; nothing to check")
        return false
    end

    mismatches = [(name, c, d) for (name, (c, d)) in counts if c != d]
    if isempty(mismatches)
        return true
    end

    println("LEAK SUSPECTED in $file (Creations != Destructions):")
    println("  $(rpad("Object Type", 24)) | $(rpad("Creations", 10)) | Destructions")
    for (name, c, d) in sort(mismatches, by = x -> x[1])
        printstyled("  $(rpad(name, 24)) | $(rpad(c, 10)) | $d\n", color = :red)
    end

    return false
end


"""
    Procpartname = CreatePartitioningFile_local(ParamFile::String, cores::Int64=1, args::String=""; bin_dir="../../bin", deb=false,mpiexec="mpiexec", dylibs="")

Create a processor partitioning file with a locally build version of LaMEM (potentially compiled vs. dynamic libraries)
"""
function CreatePartitioningFile_local(ParamFile::String, cores::Int64=1, args::String=""; 
                LaMEM_dir="../../bin", deb=false,
                mpiexec="mpiexec", verbose=false)
    

	if cores==1	& verbose==true
		return print("No partitioning file required for 1 core model setup \n")	
	end
    
    cur_dir = pwd()
	ParamFile    = abspath(ParamFile)
	args         = args*"-mode save_grid"

    # run local lamem & save output to file. This takes care of locally build LaMEM vs 
    run_lamem_local_test(ParamFile, cores, args; 
            outfile="savegrid.log", bin_dir=LaMEM_dir, deb=deb,
            mpiexec=mpiexec)
            
    logoutput = String(read("savegrid.log"))

    arr          = split(logoutput,"\n")
	foundline    = get_line_containing(arr,"Processor grid  [nx, ny, nz]         : ")
	foundline    = join(map(x -> isspace(foundline[x]) ? "" : foundline[x], 1:length(foundline)))
	
	sprtlftbrkt  = split(foundline,"[")
	sprtrghtbrkt = split(sprtlftbrkt[3],"]")
	separatecoma = split(sprtrghtbrkt[1],",")
	procnumbers  = parse.(Int, separatecoma)
	Procpartname = "ProcessorPartitioning_$(cores)cpu_$(procnumbers[1]).$(procnumbers[2]).$(procnumbers[3]).bin" 
	if !isfile(joinpath((splitdir(ParamFile)[1]),Procpartname))
		Procpartname = nothing
	end
	cd(cur_dir)
    return Procpartname

    #=
	arr          = split(logoutput,"\n")
	foundline    = get_line_containing(arr,"Processor grid  [nx, ny, nz]         : ")
	foundline    = join(map(x -> isspace(foundline[x]) ? "" : foundline[x], 1:length(foundline)))
	sprtlftbrkt  = split(foundline,"[")
	sprtrghtbrkt = split(sprtlftbrkt[3],"]")
	separatecoma = split(sprtrghtbrkt[1],",")
	procnumbers  = parse.(Int, separatecoma)
	Procpartname = "ProcessorPartitioning_$(cores)cpu_$(procnumbers[1]).$(procnumbers[2]).$(procnumbers[3]).bin" 
	if isfile(joinpath((splitdir(ParamFile)[1]),Procpartname))
		return Procpartname
	else
	return Nothing
	end
    =#
end


"""
    out = extract_info_logfiles(file::String, keywords::NTuple{N,String}=("|Div|_inf","|Div|_2","|mRes|_2"), split_sign="=", remove_substrings="")

Extracts values from the logfile `file` specified after `keywords` (optionally defining a `split_sign`).
This will generally return a NamedTuple with Vectors 

"""
function extract_info_logfiles(file::String, keywords::NTuple{N,String}=("|Div|_inf","|Div|_2","|mRes|_2"), split_signs="=") where N

    split_sign = split_signs;
    out=()
    for i=1:N        
        d =  []
        if isa(split_signs,Tuple)
            split_sign = split_signs[i]
        end
        
        open(file) do f
            while ! eof(f) 
                # read a new / next line for every iteration          
                line = readline(f)
                if contains(line, keywords[i])
                    num = extract_value_from_string(line, keywords[i], split_sign);
                    push!(d,num)    # add value to vector
                end
            end
        end
        
        out =  (out..., Float64.(d))    # add vector to tuple
    end

    out_NT = NamedTuple{Symbol.(keywords)}(out)

    return out_NT
end


"""
    extract_value_from_string(line_origin::String, keyword::String, split_sign="=", type=Float64)

This extracts a numerical value from `line_origin` after the `keyword`. Optionally, there can be a `split_sign` in between and we can transfer it to another `type` than `Float64`.`
"""
function extract_value_from_string(line_origin::String, keyword::String, split_sign="=", type=Float64)
    
    if length(split_sign)==0
        num = extract_lastvalue_from_string(line_origin, keyword, type)
    else
        # find keyword
        ind  = findlast(keyword,line_origin)
        line = line_origin[ind[end]+1:end]

        # remove split signs
        line = replace(line, split_sign=>"") 
        line = replace(line, ","=>"")           # comma

        # Extract value
        num=NaN
        if !isempty(ind)
            line_vec = split(line)
            try
                num = parse(type,line_vec[1])
            catch
                println("Problem parsing line: $line_origin")
            end
        end
    end

    return num
end

"""
    extract_lastvalue_from_string(line_origin::String, keyword::String, type=Float64)

This extracts the last numerical value from `line_origin` after the `keyword`.
"""
function extract_lastvalue_from_string(line_origin::String, keyword::String, type=Float64)
    
    # find keyword
    ind  = findlast(keyword,line_origin)
    line = line_origin[ind[end]+1:end]

    # Extract value
    num=[];
    line_vec = split(line)
    for i=1:lastindex(line_vec)
        val=nothing
        try
            val = parse(type,line_vec[i])
        catch
        end
        if !isnothing(val)
            push!(num,val)
        end
    end
    
    return num[1]
end


"""
    success = compare_logfiles(new::String, expected::String, 
                        keywords::NTuple{N,String}=("|Div|_inf","|Div|_2","|mRes|_2"), 
                        accuracy=((rtol=1e-6,), (rtol=1e-6,), (rtol=1e-6,)),
                        split_sign="=",
                        remove_substring="")

This compares two logfiles (different parameters which can be indicated). If the length of the vectors is not the same, or the accuracy criteria are not met, `success=false` and info is displayed, to help track down the matter.
We scan the file for lines with the given `keywords`, and extract numerical values from it.

Arguments:
    `split_sign`        : Can be a `Tuple`, containing the sign after which we split the string
    
"""
function compare_logfiles(new::String, expected::String, 
                        keywords::NTuple{N,String}=("|Div|_inf","|Div|_2","|mRes|_2"), 
                        accuracy=((rtol=1e-6,), (rtol=1e-6,), (rtol=1e-6,));
                        split_sign="=") where N

    new_out = extract_info_logfiles(new, keywords, split_sign)
    exp_out = extract_info_logfiles(expected, keywords, split_sign)
    
    test_status = true
    for i=1:N 
        rtol, atol = 0,0
        if  haskey(accuracy[i], :rtol)
            rtol = accuracy[i].rtol;
        end
        if  haskey(accuracy[i], :atol)
            atol = accuracy[i].atol;
        end
        if length(new_out[i])==length(exp_out[i])
            te =  isapprox(new_out[i], exp_out[i], rtol=rtol, atol=atol)
            if te==false
                println("Problem with comparison of $(keywords[i]):")
                print_differences( new_out[i], exp_out[i], accuracy[i])
                test_status = false
            end
        else
            println("Problem with comparison of $(keywords[i]):")
            println("length of vectors not the same (new: $(length(new_out[i])), expected: $(length(exp_out[i]))")
            test_status = false
        end
       
    end
  

    return test_status
end

# Pretty formatting of errors
function print_differences(new, expected, accuracy)
    n = 24;
    atol_expected = 0.0;
    if haskey(accuracy, :atol)
        atol_expected = accuracy.atol;
    end
    rtol_expected = 0.0;
    if haskey(accuracy, :rtol)
        rtol_expected = accuracy.rtol;
    end

    println("      $(rpad("New",n)) | $(rpad("Expected",n)) | $(rpad("rtol (<$(rtol_expected))",n)) | $(rpad("atol  (<$(atol_expected))",n))")

    for i in eachindex(new)
        atol = norm(new[i] - expected[i])
        rtol = atol/max(norm(new[i]), norm(expected[i]))
        col = :normal
        if atol >  max(atol_expected, rtol_expected*max(norm(new[i]), norm(expected[i])))
            col = :red
        end
        printstyled("$(rpad(i,4))  $(rpad(new[i],n)) | $(rpad(expected[i],n)) | $(rpad(rtol,n)) | $(rpad(atol,n)) \n", color=col)
    end

    return nothing
end



"""
    clean_test_directory(dir)

This cleans a certain (test) directory of all that LaMEM and the tests may have generated.
"""
function clean_test_directory(dir)
    
    cur_dir = pwd();

    clean_directory(dir)
        
    cd(dir)
    for f in glob("*.out")
        rm(f)
    end
    for f in glob("*.log")
        rm(f)
    end
    for f in glob("*.bin")
        rm(f)
    end
    for f in glob("markers*")
        rm(f, force=true, recursive=true)
    end
    for f in glob("restart")
        rm(f, force=true, recursive=true)
    end
    cd(cur_dir)  # return to directory       

end


"""
    get_dylibs()
This retrieves dynamic libraries, required to run LaMEM. It assumes that the global variable `use_dynamic_lib` is present
"""
function get_dylibs()
    if use_dynamic_lib
        # NOTE: LIBPATH is a Ref{String}; it must be dereferenced. Passing the Ref itself
        # stringifies to `Base.RefValue{String}("...")`, which no loader can parse.
        dylibs = PETSc_jll.LIBPATH[];

        mpi_path = if PETSc_jll.MPICH_jll.is_available()
            PETSc_jll.MPICH_jll.PATH_list[1]
        elseif PETSc_jll.MicrosoftMPI_jll.is_available()
            PETSc_jll.MicrosoftMPI_jll.PATH_list[1]
        elseif PETSc_jll.OpenMPI_jll.is_available()
            PETSc_jll.OpenMPI_jll.PATH_list[1]
        else
            nothing
        end

    else
        dylibs = ""
        mpi_path = ""
    end
    
    return dylibs, mpi_path
end


"""
    add_dylibs(cmd, dylibs)
Adds the runtime environment a LaMEM binary built against the PETSc_jll libraries needs:

- The dynamic library search path, via the loader variable that is correct for the current
  platform (`LD_LIBRARY_PATH` on Linux, `DYLD_FALLBACK_LIBRARY_PATH` on macOS, `PATH` on
  Windows); `JLLWrappers.LIBPATH_env` names that variable for us.

- `LBT_DEFAULT_LIBS`. PETSc_jll >= 3.25 no longer links OpenBLAS directly but libblastrampoline
  (LBT). LBT ships no BLAS of its own: it forwards to whatever library it is told to load. Inside
  Julia, LinearAlgebra does that forwarding; a standalone executable gets an empty trampoline and
  either segfaults on its first BLAS call (VecNorm -> BLASdot) or floods stderr with
  "no BLAS/LAPACK library loaded for dgemm_()". `LBT_DEFAULT_LIBS` tells LBT what to forward to
  at load time, as a `;`-separated list.

  BOTH interfaces have to be loaded. The stack mixes them: libpetsc, UMFPACK and CHOLMOD call the
  ILP64, `_64_`-suffixed symbols (`dgemm_64_`), while SuperLU_DIST and MUMPS call the unsuffixed
  LP64 ones (`dgemm_`). Loading only the ILP64 half leaves the direct solvers with an empty
  trampoline - a MUMPS run emits ~78k "no BLAS/LAPACK library loaded" errors and aborts. So
  forward Julia's stdlib OpenBLAS (ILP64 `libopenblas64_`) *and* OpenBLAS32_jll (LP64), which
  PETSc_jll depends on for exactly these libraries.

- Single-threaded BLAS and OpenMP. SuperLU_DIST is built with OpenMP and calls the LP64 BLAS
  from inside its parallel regions. With threading left at its default, those calls
  intermittently reach OpenBLAS with a corrupted argument and the run prints

      ** On entry to DGEMM  parameter number  8 had an illegal value

  dozens of times, after which the solve produces garbage. It is nondeterministic: the same
  binary on the same input fails most runs, which is why a test can pass in one CI run and
  fail in the next on an identical tree.

  This is not the whole story for t23_Permeable, which still diverges on Linux/Int64 with
  superlu_dist even with the threads pinned (see PR #84); that test stays on MUMPS.

  `OMP_NUM_THREADS` is the variable that matters - pinning `OPENBLAS_NUM_THREADS` alone does not
  help, since the threads in question are SuperLU_DIST's, not OpenBLAS's. `OPENBLAS_NUM_THREADS`
  and `VECLIB_MAXIMUM_THREADS` are set alongside it to keep the BLAS itself single-threaded, the
  same set LaMEM.jl pins in `deactivate_multithreading`. The test suite compares residuals
  against stored reference values, so it wants reproducible arithmetic rather than threaded
  speed anyway.
"""
function add_dylibs(cmd::Cmd, dylibs)
    # Thread pinning applies to every run, including a locally built PETSc (no dylibs to add):
    # the OpenMP race below is a property of SuperLU_DIST, not of the PETSc_jll packaging.
    cmd = addenv(cmd, "OMP_NUM_THREADS"        => "1",
                      "OPENBLAS_NUM_THREADS"   => "1",
                      "VECLIB_MAXIMUM_THREADS" => "1")

    if isempty(dylibs)
        return cmd
    end

    # ILP64 (libpetsc, UMFPACK, CHOLMOD) and LP64 (SuperLU_DIST, MUMPS); LBT needs both
    blas_libs = join([OpenBLAS_jll.libopenblas_path,
                      PETSc_jll.OpenBLAS32_jll.libopenblas_path], ";")

    return addenv(cmd, PETSc_jll.JLLWrappers.LIBPATH_env => dylibs,
                       "LBT_DEFAULT_LIBS"               => blas_libs)
end


"""
    perform_lamem_test( dir::String, 
                        ParamFile::String, 
                        expectedFile::String; 
                        keywords=("|Div|_inf","|Div|_2","|mRes|_2"), 
                        accuracy=((rtol=1e-6,atol=1e-10), (rtol=1e-6,), (rtol=1e-6,)), 
                        cores::Int64=1, 
                        args::String="",
                        bin_dir="../bin",  
                        deb=false, 
                        mpiexec="mpiexec",
                        split_sign="=", 
                        debug::Bool=false, 
                        create_expected_file::Bool=false, 
                        clean_dir::Bool=true,
                        valgrind::Bool = (@isdefined(use_valgrind) ? use_valgrind : false),
                        memcheck::Bool = (@isdefined(use_memcheck) ? use_memcheck : false))

This performs a LaMEM simulation and compares certain keywords of the logfile with results of a previous simulation        

Parameters:
- `dir`: directory in which the LaMEM `*.dat` ParamFile is located
- `ParamFile`: name of the LaMEM input file
- `expectedFile`: name of the file with earlier results, versus which we compare (WARNING! WITHOUT EXTENSION)
- `keywords`: Tuple with keywords, which contain numerical values
- `accuracy`: Tuple that contains relative (`rtol`), and (optionally) absolute `atol` tolerances
- `cores`: Number of cores on which to perform the test
- `args`: Optional LaMEM command line Arguments
- `bin_dir`: directory where the LaMEM binaries are, relative to the current one
- `deb`: run with debug version of LaMEM?
- `mpiexec`: mpi executable
- `split_sign`: split sign (or Tuple of it)
- `debug`: set to true if you simply want to see the output of the simulation (no test done)
- `create_expected_file`: create an expected file
- `clean_dir`: delete all timestep & pvd files at the end?
- `valgrind`: run this test through Valgrind instead of running LaMEM directly.
   If not specified explicitly, this defaults to the global `use_valgrind` flag (set by `runtests.jl`
   when `julia start_tests.jl ... valgrind` is used), so existing test definitions don't need to change.
- `memcheck`: append `-log_view` to the run and, instead of the usual keyword-based numeric comparison,
   check that every PETSc object type reported by `-log_view` was destroyed as many times as it was
   created. The test fails if there is a mismatch (a likely missing `*Destroy()` call).
   If not specified explicitly, this defaults to the global `use_memcheck`.

"""
function perform_lamem_test(dir::String, ParamFile::String, expectedFile::String; 
                keywords=("|Div|_inf","|Div|_2","|mRes|_2"), accuracy=((rtol=1e-6,), (rtol=1e-6,), (rtol=1e-6,)), 
                cores::Int64=1, args::String="",
                bin_dir="../bin", deb=false, mpiexec="mpiexec",
                split_sign="=", 
                debug::Bool=false, create_expected_file::Bool=false, clean_dir::Bool=true,
                valgrind::Bool = (@isdefined(use_valgrind) ? use_valgrind : false),
                memcheck::Bool = (@isdefined(use_memcheck) ? use_memcheck : false))
   
    if valgrind == true
        valgrind_flag = "[valgrind]"
    elseif memcheck == true
        valgrind_flag = "[memcheck]"
    else
        valgrind_flag = ""
    end
    
    if deb
        opt_mod="deb"
    else
        opt_mod="opt"
    end 

    # print info about running tests                
    @info "Performing test $ParamFile in directory $dir on $cores cores in $opt_mod mode $valgrind_flag"
    
    cur_dir = pwd();
    cd(dir)

    bin_dir = joinpath(cur_dir,bin_dir);
    if debug==true
        outfile = "";
    else
        outfile = "$expectedFile.out";
    end
    if create_expected_file==true
        outfile = "$expectedFile.expected";
        debug = true;
    end
	
	expectedFile = "$expectedFile.expected";

    if memcheck
        # -log_view makes PETSc print, per object type, how many objects were created vs.
        # destroyed. It is appended to whatever command-line args the test already uses,
        # and is written to the same stdout log used below.
        args = isempty(args) ? "-log_view" : args * " -log_view"
    end

    # perform simulation 
    success = run_lamem_local_test(ParamFile, cores, args, outfile=outfile, bin_dir=bin_dir, deb=deb, mpiexec=mpiexec, valgrind=valgrind);

    # Under Valgrind we always run the debug binary (see run_lamem_local_test),
    # Since the point of a valgrind run is leak detection, not numeric validation, skip that comparison here.
    if success==true && debug==false && memcheck
        # Memory-leak check instead of the usual numeric comparison - see check_memory_usage.
        success = check_memory_usage(outfile)
    elseif success==true && debug==false && !valgrind
        # compare logfiles 
        success = compare_logfiles(outfile, expectedFile, keywords, accuracy, split_sign=split_sign)
    end
    if create_expected_file==true
        println("Created expected file: $expectedFile in directory $dir")
    end
    if !success
        # something went wrong with executing the file (likely @ PETSc error)
        # Display some useful info here that helps debugging
        println("Problem detected with test; see this on commandline with: ")
        println("  dir=$(joinpath(cur_dir,dir)) ")
        println("  ParamFile=$(ParamFile) ")
        println("  cores=$(cores) ")
        println("  args=$(args) ")
        println("  outfile=$(outfile) ")
        println("  bindir=$(bin_dir) ")
        println("  deb=$(deb) ")
        println("  mpiexec=$(mpiexec) ")
        println("  success = run_lamem_local_test(ParamFile, cores, args, outfile=nothing, bin_dir=bin_dir, deb=deb, mpiexec=mpiexec);")
    end

    cd(cur_dir)  # return to directory       

    if clean_dir
       clean_test_directory(dir)
    end 
    
    return success
end 




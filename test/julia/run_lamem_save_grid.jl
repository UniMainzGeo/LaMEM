# This contains routines to saving processor partitioning file to run LaMEM from julia.
# Returns the name of the processor partitioning file
#
# Note: This downloads the BinaryBuilder version of LaMEM, which is not necessarily the latest version of LaMEM 
#       (or the same as the current repository), since we have to manually update the builds.
using Base.Sys   


function run_lamem_with_log(ParamFile::String, cores::Int64=1, args::String=""; wait=true, deactivate_multithreads=true)
    if iswindows() && cores>1
        cores=1;
        println("LaMEM_jll does not support parallel runs on windows; using 1 core instead")
    end
	out = Pipe()
    if cores==1
        # Run LaMEM on a single core, which does not require a working MPI
        cmd = `$(LaMEM_jll.LaMEM()) -ParamFile $(ParamFile) $args`
        if deactivate_multithreads
            cmd = deactivate_multithreading(cmd)
        end
		cmd1 = pipeline(ignorestatus(cmd),stdout=out)
		run(cmd1, wait=wait);
		close(out.in)
	else
	
        # set correct environment
       # mpirun = setenv(mpiexec, LaMEM_jll.JLLWrappers.JLLWrappers.LIBPATH_env=>LaMEM_jll.LIBPATH[]);
		mpirun = mpiexec

        # create command-line object
		cmd = `$(mpirun) -n $cores --map-by -ParamFile $(ParamFile) $args `
        if deactivate_multithreads
            cmd = deactivate_multithreading(cmd)
        end

        # Run LaMEM in parallel
		cmd1 = pipeline(ignorestatus(cmd),stdout=out)
        run(cmd1, wait=wait);
		close(out.in)
    end
	stdout = String(read(out))

    return stdout
end

"""
    deactivate_multithreading(cmd)

This deactivates multithreading
"""
function deactivate_multithreading(cmd::Cmd)
    # multithreading of the BLAS libraries that is installed by default with the julia BLAS
    # does not work well. Switch that off:
    cmd = addenv(cmd,"OMP_NUM_THREADS"=>1)
    cmd = addenv(cmd,"VECLIB_MAXIMUM_THREADS"=>1)
    return cmd
end


function JuliaStringToArray(input)


    arr = split(input,"\n")
	return arr
end

#=
function get_line_containing(stringarray::Vector{SubString{String}}, lookfor::String)


	for line in stringarray
		   if contains(line, lookfor)
		   foundline=line
		   return foundline
		   end
	end
end
=#

""" 
	ProcessorPartFile = run_lamem_save_grid(ParamFile::String, cores::Int64=1; verbose=true, directory=pwd())
This calls LaMEM simulation, for using the parameter file `ParamFile` 
and creates processor partitioning file `"ProcessorPartitioning_Xcpu_X.Y.Z.bin"` for `{X}` number of cores. 

# Example:
```julia
julia> using LaMEM
julia> ParamFile="../../input_models/BuildInSetups/FallingBlock_Multigrid.dat";
julia> ProcessorPartFile = run_lamem_save_grid(ParamFile, 2)
```
"""
function run_lamem_save_grid(ParamFile::String, cores::Int64=1; verbose=true, directory=pwd())
	if cores==1	& verbose==true
		return print("No partitioning file required for 1 core model setup \n")	
	end
	if iswindows() && cores>1
        cores=1;
        println("LaMEM_jll does not support parallel runs on windows; using 1 core instead")
    end
	cur_dir = pwd();
	cd(directory)

	ParamFile    = abspath(ParamFile)
	logoutput    = run_lamem_with_log(ParamFile, cores,"-mode save_grid" )
	
	arr          = JuliaStringToArray(logoutput)
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
end

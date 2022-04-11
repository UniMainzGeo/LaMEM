# This contains routines to run LaMEM from julia.
#
# Note: This downloads the BinaryBuilder version of LaMEM, which is not necessarily the latest version of LaMEM 
#       (or the same as the current repository), since we have to manually update the builds.




using LaMEM_jll

# load the correct mpi
const mpiexec = if isdefined(LaMEM_jll,:MPICH_jll)
    LaMEM_jll.MPICH_jll.mpiexec()
elseif isdefined(LaMEM_jll,:MicrosoftMPI_jll) 
    LaMEM_jll.MicrosoftMPI_jll.mpiexec()
else
    nothing
end

mpirun = addenv(mpiexec, LaMEM_jll.JLLWrappers.JLLWrappers.LIBPATH_env=>LaMEM_jll.LIBPATH[]);
    
""" 
    run_lamem(ParamFile::String, cores::Int64=1, args:String="")

This starts a LaMEM simulation, for using the parameter file `ParamFile` on `cores` number of cores. 
Optional additional command-line parameters can be specified with `args`.

# Example:
The first step is to ensure that `LaMEM_jll` is installed on your system. You only need to do this once, or once LaMEM_jll is updated. 
```julia
julia> import Pkg
julia> Pkg.add("LaMEM_jll")
```

Next you can call LaMEM with:
```julia
julia> ParamFile="../../input_models/BuildInSetups/FallingBlock_Multigrid.dat";
julia> run_lamem(ParamFile)
```

Do the same on 2 cores with
```julia
julia> ParamFile="../../input_models/BuildInSetups/FallingBlock_Multigrid.dat";
julia> run_lamem(ParamFile, 2)
```

"""
function run_lamem(ParamFile::String, cores::Int64=1, args::String="")

    if cores==1
        # Run LaMEM on a single core, which does not require a working MPI
        run(`$(LaMEM_jll.LaMEM()) -ParamFile $(ParamFile) $(args)`);
    else
        # Run LaMEM in parallel
        run(`$(mpirun) -n $cores $(LaMEM_jll.LaMEM_path) -ParamFile $(ParamFile) $(args)`);
    end

    return nothing
end
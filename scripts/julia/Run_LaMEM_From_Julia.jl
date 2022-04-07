# This contains routines to run LaMEM from julia.
# Note that this downloads the BinaryBuilder version of LaMEM, which is not necessarylu 

using LaMEM_jll

# load the correct mpi
const mpiexec = if LaMEM_jll.MPICH_jll.is_available()
    LaMEM_jll.MPICH_jll.mpiexec()
elseif MAGEMin_jll.MicrosoftMPI_jll.is_available()
    LaMEM_jll.MicrosoftMPI_jll.mpiexec()
else
    nothing
end

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

    
    mpirun = addenv(mpiexec, LaMEM_jll.JLLWrappers.JLLWrappers.LIBPATH_env=>LaMEM_jll.LIBPATH[]);
    
    # Run LaMEM in parallel
    run(`$(mpirun) -n $cores $(LaMEM_jll.LaMEM_path) -ParamFile $(ParamFile) $(args)`);

    return nothing
end
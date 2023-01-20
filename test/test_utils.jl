# These are tools that help perform the LaMEM tests run LaMEM locally

"""
    run_lamem_local_test(ParamFile::String, cores::Int64=1, args::String=""; 
                        outfile="test.out", bin_dir="../../bin", opt=true, deb=false,
                        mpiexec="mpiexec")

This runs a LaMEM simulation with given `ParamFile` on 1 or more cores, while writing the output to a local log file.

"""
function run_lamem_local_test(ParamFile::String, cores::Int64=1, args::String=""; 
                outfile="test.out", bin_dir="../../bin", opt=true, deb=false,
                mpiexec="mpiexec")
    
    cur_dir = pwd()
    if opt
        exec=joinpath(cur_dir,bin_dir,"opt","LaMEM")
    elseif deb
        exec=joinpath(cur_dir,bin_dir,"deb","LaMEM")
    end

   
    if cores==1
        perform_run = `$(exec) -ParamFile $(ParamFile) $args`;
        
        # Run LaMEM on a single core, which does not require a working MPI
        if !isempty(outfile)
            run(pipeline(perform_run, stdout=outfile));
        else
            run(perform_run);
        end
    else
        perform_run = `$(mpiexec) -n $(cores) $(exec) -ParamFile $(ParamFile) $args`;
        # set correct environment
        #mpirun = setenv(mpiexec, LaMEM_jll.JLLWrappers.JLLWrappers.LIBPATH_env=>LaMEM_jll.LIBPATH[]);

        # Run LaMEM in parallel
        if !isempty(outfile)
            run(pipeline(perform_run, stdout=outfile));
        else
            run(perform_run);
        end

    end

    return nothing
end


"""
    out = extract_info_logfiles(file::String, keyword::NTuple{N,String}=("|Div|_inf","|Div|_2","|mRes|_2"), split_sign="=")

Extracts values from the logfile `file` specified after `keywords` (optionally defining a `split_sign`).
This will generally return a NamedTuple with Vectors 
"""
function extract_info_logfiles(file::String, keyword::NTuple{N,String}=("|Div|_inf","|Div|_2","|mRes|_2"), split_sign="=") where N

    out=()
    for i=1:N        
        d =  []
        open(file) do f
            while ! eof(f) 
                # read a new / next line for every iteration          
                line = readline(f)
                    if contains(line, keyword[i])
                        num = parse(Float64,split(line,split_sign)[end])
                        push!(d,num)    # add value to vector
                    end
            end
        end
        
        out =  (out..., Float64.(d))    # add vector to tuple
    end

    out_NT = NamedTuple{Symbol.(keyword)}(out)

    return out_NT
end
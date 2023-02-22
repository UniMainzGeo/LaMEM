# These are tools that help perform the LaMEM tests, which run LaMEM locally
using LinearAlgebra
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

    success = true
    if cores==1
        perform_run = `$(exec) -ParamFile $(ParamFile) $args`;
        
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
        perform_run = `$(mpiexec) -n $(cores) $(exec) -ParamFile $(ParamFile) $args`;
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
  

    return success
end


"""
    out = extract_info_logfiles(file::String, keywords::NTuple{N,String}=("|Div|_inf","|Div|_2","|mRes|_2"), split_sign="=", remove_substrings="")

Extracts values from the logfile `file` specified after `keywords` (optionally defining a `split_sign`).
This will generally return a NamedTuple with Vectors 

"""
function extract_info_logfiles(file::String, keywords::NTuple{N,String}=("|Div|_inf","|Div|_2","|mRes|_2"), split_signs="=", remove_substrings="") where N

    split_sign = split_signs;
    remove_substring = remove_substrings;
    out=()
    for i=1:N        
        d =  []
        if isa(split_signs,Tuple)
            split_sign = split_signs[i]
        end
        if isa(remove_substrings,Tuple)
            remove_substring = remove_substrings[i]
        end
        
        open(file) do f
            while ! eof(f) 
                # read a new / next line for every iteration          
                line = readline(f)
                if contains(line, keywords[i])
                    if contains(line, "[")
                        # remove everything in between [  ] (units etc.)
                        ind = (findfirst("[", line)[1],findlast("]", line)[1])
                        line = line[1:ind[1]-1]*line[ind[end]+1:end]
                    end 
                    if !isempty(remove_substring)
                        if contains(line,remove_substring)
                            line = replace(line, remove_substring=>"")
                        end
                    end
                    num=NaN
                    try
                        num = parse(Float64,split(line,split_sign)[end])
                    catch
                        error("Problem parsing line: $line")
                    end
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
    success = compare_logfiles(new::String, expected::String, 
                        keywords::NTuple{N,String}=("|Div|_inf","|Div|_2","|mRes|_2"), 
                        accuracy=((rtol=1e-6,), (rtol=1e-6,), (rtol=1e-6,)),
                        split_sign="=",
                        remove_substring="")

This compares two logfiles (different parameters which can be indicated). If the length of the vectors is not the same, or the accuracy criteria are not met, `success=false` and info is displayed, to help track down the matter.
We scan the file for lines with the given `keywords`, and extract numerical values from it.

Arguments:
    `split_sign`        : Can be a `Tuple`, containing the sign after which we split the string
    `remove_substring`  : Optional `Tuple` with substrings to be stripped from the line, before numerical value is extracted
    

"""
function compare_logfiles(new::String, expected::String, 
                        keywords::NTuple{N,String}=("|Div|_inf","|Div|_2","|mRes|_2"), 
                        accuracy=((rtol=1e-6,), (rtol=1e-6,), (rtol=1e-6,));
                        split_sign="=", remove_substring="") where N

    new_out = extract_info_logfiles(new, keywords, split_sign, remove_substring)
    exp_out = extract_info_logfiles(expected, keywords, split_sign, remove_substring)

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
    println("      $(rpad("New",n)) | $(rpad("Expected",n)) | $(rpad("rtol (<$(accuracy.rtol))",n)) | $(rpad("atol  (<$(accuracy.atol))",n))")

    for i=1:length(new)
        atol = norm(new[i] - expected[i])
        rtol = atol/max(norm(new[i]), norm(expected[i]))
        col = :normal
        if atol>  max(accuracy.atol, accuracy.rtol*max(norm(new[i]), norm(expected[i])))
            col = :red
        end
        printstyled("$(rpad(i,4))  $(rpad(new[i],n)) | $(rpad(expected[i],n)) | $(rpad(rtol,n)) | $(rpad(atol,n)) \n", color=col)
    end

    return nothing
end


function perform_lamem_test(dir::String, ParamFile::String, expectedFile::String; 
                keywords=("|Div|_inf","|Div|_2","|mRes|_2"), accuracy=((rtol=1e-6,), (rtol=1e-6,), (rtol=1e-6,)), 
                cores::Int64=1, args::String="",
                bin_dir="../bin",  opt=true, deb=false, mpiexec="mpiexec",
                debug=false, split_sign="=", remove_substring="")

    cur_dir = pwd();
    cd(dir)

    bin_dir = joinpath(cur_dir,bin_dir);
    if debug==true
        outfile = "";
    else
        outfile = "test_$(cores).out";
    end

    # perform simulation 
    success = run_lamem_local_test(ParamFile, cores, args, outfile=outfile, bin_dir=bin_dir, opt=opt, deb=deb, mpiexec=mpiexec);

    if success && debug==false
        # compare logfiles 
        success = compare_logfiles(outfile, expectedFile, keywords, accuracy, split_sign=split_sign, remove_substring=remove_substring)
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
        println("  opt=$(opt) ")
        println("  deb=$(deb) ")
        println("  mpiexec=$(mpiexec) ")
        println("  success = run_lamem_local_test(ParamFile, cores, args, outfile=nothing, bin_dir=bin_dir, opt=opt, deb=deb, mpiexec=mpiexec);")
    end

    cd(cur_dir)  # return to directory       
    
    return success
end 


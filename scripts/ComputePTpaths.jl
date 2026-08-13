# This is a small wrapper to process all passive tracer files in the current directory 
#  and generate a matlab file from it
#
# Make sure that the path below is correct and run it as:
#   julia ComputePTpaths.jl

push!(LOAD_PATH, "../../scripts/julia/");   # path
using PassiveTracers                        # load julia module
ExtractPTpaths()                            # transfer


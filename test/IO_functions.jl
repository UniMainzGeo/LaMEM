module IO_functions
# this contains I/O routines of LaMEM, which don't require LaMEM_jll

include("read_timestep.jl")
export read_LaMEM_PVTR_file, read_LaMEM_PVTS_file, field_names, readPVD, read_LaMEM_PVTU_file


#include("utils_IO.jl")
#export IO_functions.clean_directory, changefolder, project_onto_crosssection


end
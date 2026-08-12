# This is a small wrapper to process all passive tracer files in the current directory 
#  and generate a matlab file from it
#
# Make sure that the path below is correct and run it as:
#   julia TransferVTK_To_Matlab.jl   <FileName>
#
# Hete <FileName> is the filename of the simulation (as specified in the *.pvd file)


push!(LOAD_PATH, "../../scripts/julia/");   # path
using ReadLaMEM_Timestep                    # load julia module

if length(ARGS)>0
    FileName    =   ARGS[1];                # read filename from command (or change that if you wish)
else
    FileName    =   "PlumeLithosphereInteraction";  
end

TransferVTK2MAT(FileName)                   # transfer to matlab


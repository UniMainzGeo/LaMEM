# This processes all Timestep directories in the current LaMEM 
# directory and extracts the minimum location of phase 1

# Read LaMEM routines
include("../../scripts/julia/ReadLaMEM_Timestep.jl")



" Reads all Timesteps in the current directory and extracts the minimum coordinate of phase 0

"
function ReadMinZ(FileName)

    # define a macro that finds all directories with a certain pattern
    searchdir(path,key) = filter(x->occursin(key,x), readdir(path))

    if length(ARGS)>0
        FileName    =   ARGS[1];            # command line argument
        FileName    = "$FileName.pvtr"
    else
        FileName    =   "RT_FSSA.pvtr";
    end

    print("Processing file name $FileName \n")

    # Create a new directory that will contain the 2D plots
    if ~ispath("Visualization_2D")
        mkdir("Visualization_2D")
    end

    # Go over all Timestep subdirectories in the current directory
    Directories         = searchdir(pwd(),"Timestep_")

    Time_vec        = zeros(length(Directories))
    Zmin_vec        = zeros(length(Directories))
    for iStep=1:length(Directories)

        DirName     =   Directories[iStep];
        print("Processing directory $DirName \n")

        # Extract the timestep info from the directory directory name 
        id      =   findlast("_",DirName);
        Time    =   parse(Float64,DirName[id[1]+1:length(DirName)]);

        # ------------------------------------------------------------------------------------- 
        # Read data from timestep
        data        =   Read_VTR_File(DirName, FileName);

        # Extract 2D fields from file
        phase,x,z   =   ReadField_2D_pVTR(data, "phase [ ]"     ,       "Scalar");     # phase
            
        # ------------------------------------------------------------------------------------- 

        # extract all points that have phase 1
        ind       =   findall( ((x->x<0.5)), phase)
        
        z_min = 1e6;
        for id in ind
        # global z_min;
        z_val = z[id[1]]
        x_val = x[id[2]]
        if z_val< z_min
            z_min = z_val;
        end 

        end
    #    print("Time=$Time,  z_min=$z_min\n")
        

        Time_vec[iStep] = Time;
        Zmin_vec[iStep] = z_min;
    
    end

    return Time_vec, Zmin_vec
end
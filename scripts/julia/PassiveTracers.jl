module PassiveTracers
# This module contains various routines to deal with processing data on LaMEM passive tracers 
# Usage:
#   1) Make sure you are within the LaMEM directory that contains the passive tracer files
#   2) Add include path with Julia scripts:
#        push!(LOAD_PATH, "../../scripts/julia/")
#   3) Run script, e.g.:
#          
#       
# 

# for easier usage outside the module 
export ExtractPTpaths, PostProcess_PassiveTracers

# Read required modules
using ReadLaMEM_Timestep        # Read LaMEM routines
using Printf        
using MAT                       # write matlab files



#-------------------------------------------------------------------------------------------
"
    Example that reads LaMEM/VTK passive tracers files into Julia, combined P/T evolution 
    into a single file and save that as a matlab file.
    
"
function ExtractPTpaths()
    # This processes all Timestep directories in the current LaMEM 
    # directory, extracts Passive Tracers and creates a matlab file
    # of it.
    #
    # Copy this file to your directory & change the initial path here

    # define a macro that finds all directories with a certain pattern
    searchdir(path,key) = filter(x->occursin(key,x), readdir(path))
    
    if length(ARGS)>0
        FileName    =   ARGS[1];            # command line argument
        FileName    = "$FileName.pvtu"
    else
        print("No passive tracer filename given on command-line; trying to determine one... \n")
        PT_filenames = searchdir(pwd(),"_passive_tracers.pvd");
        if length(PT_filenames)==0
            error("Cannot determime a passive tracer file. Please add it on the command-line")
        end

        PT_filenames = PT_filenames[1]
        PT_filenames = PT_filenames[1:end-4];
        FileName     = "$PT_filenames.pvtu"
    end

    print("Postprocessing passive tracer files $FileName \n")

    # Go over all Timestep subdirectories in the current directory
    Directories         =   searchdir(pwd(),"Timestep_")

    nMax                =   length(Directories)
    #nMax                =   20;
    minZ                =   [];
    for iStep=1:nMax
        local DirName, data, points_coord, ID, Active, P,T,zCoord, xCoord, New, id, Time
        local OutFileName;
        global minZ, OutNames, Time_vec, P_mat, T_mat, x_mat, z_mat;
        global maxP, maxT;

        DirName     =   Directories[iStep];
        print("Processing directory $DirName \n")

        # Extract the timestep info from the directory directory name 
        id      =   findlast("_",DirName);
        Time    =   parse(Float64,DirName[id[1]+1:length(DirName)]);

        # ------------------------------------------------------------------------------------- 
        # Read data from timestep
        data        =   Read_VTU_File(DirName, FileName);
        
        points_coord =  ReadField_VTU(data,"coords");                   # coordinate array
        ID          =   ReadField_VTU(data,"ID");                       # ID of all points
        Active      =   ReadField_VTU(data,"Active");                   # Active or not?
        P           =   ReadField_VTU(data,"Pressure [MPa]");           # Pressure
        T           =   ReadField_VTU(data,"Temperature [C]");          # Temperature
        MeltFrac    =   ReadField_VTU(data,"Mf_Grid [ ]");              # Melt fraction
        zCoord      =   points_coord[:,3];   
        xCoord      =   points_coord[:,1];   

        numP        =   length(P);
        if iStep==1
            Time_vec        =   Array{Float64}(undef,   nMax);
            P_mat           =   Array{Float32}(undef,   numP, nMax);
            T_mat           =   Array{Float32}(undef,   numP, nMax);
            x_mat           =   Array{Float32}(undef,   numP, nMax);
            z_mat           =   Array{Float32}(undef,   numP, nMax);
        end 

        # Max P & T
        id          =   findall(x->x==1, Active);
        if sizeof(id)>0
            maxP[id]   =   max.(maxP[id],P[id]);                  
            maxT[id]   =   max.(maxT[id],T[id]);                  
        end

        # Save Data in matrix
        numP                =   length(P);
        Time_vec[iStep]     =   Time;
        P_mat[:,iStep]      =   P[:];
        T_mat[:,iStep]      =   T[:];
        x_mat[:,iStep]      =   xCoord[:]; 
        z_mat[:,iStep]      =   zCoord[:]; 
    end

    # Generate MATLAB file that saves all the info on the particles
    file = matopen("PTpaths_PassiveTracers.mat", "w")
    write(file, "Time_vec",     Time_vec);
    write(file, "Pressure",     P_mat);
    write(file, "Temperature",  T_mat);
    write(file, "xCoord",       x_mat);
    write(file, "zCoord",       z_mat);
    close(file)
    print("Saved passive tracer data to matlab file: PTpaths_PassiveTracers.mat \n") 

end
#-------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------
"
    Example that shows how we can read in LaMEM/VTK passive tracers files into Julia, 
    do computations based on the fields that are on the passive tracers and write new 
    VTK passive tracers files to disk.

    You can modify this according to your needs

"
function PostProcess_PassiveTracers()

    # define a macro that finds all directories with a certain pattern
    searchdir(path,key) = filter(x->occursin(key,x), readdir(path))

    if length(ARGS)>0
        FileName    =   ARGS[1];            # command line argument
        FileName    = "$FileName.pvtr"
    else
        FileName    =   "PlumeLithosphereInteraction_passive_tracers.pvtu";
    end

    print("Postprocessing passive tracer files $FileName \n")

    # Go over all Timestep subdirectories in the current directory
    Directories         =   searchdir(pwd(),"Timestep_")

    nMax                =   length(Directories)
    #nMax                =   20;
    minZ                =   [];
    for iStep=1:nMax
        local DirName, data, points_coord, ID, Active, P,T,zCoord, New, id, Time
        local OutFileName;
        global minZ, OutNames, Time_vec;

        DirName     =   Directories[iStep];
        print("Processing directory $DirName \n")

        # Extract the timestep info from the directory directory name 
        id      =   findlast("_",DirName);
        Time    =   parse(Float64,DirName[id[1]+1:length(DirName)]);

        # ------------------------------------------------------------------------------------- 
        # Read data from timestep
        data        =   Read_VTU_File(DirName, FileName);
        
        points_coord =  ReadField_VTU(data,"coords");                   # coordinate array
        ID          =   ReadField_VTU(data,"ID");                       # ID of all points
        Active      =   ReadField_VTU(data,"Active");                   # Active or not?
        P           =   ReadField_VTU(data,"Pressure [MPa]");           # Pressure
        T           =   ReadField_VTU(data,"Temperature [C]");          # Temperature
        zCoord      =   points_coord[:,3];   

        if iStep==1
            minZ        =   zCoord;
            OutNames    =   Array{String}(undef,    nMax);
            Time_vec    =   Array{Float64}(undef,   nMax);
        end 

        # filter out data with -Inf (that was kind of a LaMEM hack to mark of )
        id          =   findall(x->x<-1e4, zCoord);
        if sizeof(id)>0
            Active[id]  .=   0;
        end

        # Determine minimum coordinate of active particles
        id          =   findall(x->x==1, Active);
        if sizeof(id)>0
            minZ[id]   =   min.(minZ[id],zCoord[id]);                   # obtain minimum of evert point in array
        end

        # Save data to new files --------------------------------------------------------------
    
        # Add data to object
        New = [];
        New = AddNewField_VTU(New, data, zCoord, "zCoord");        
        New = AddNewField_VTU(New, data, minZ,   "minZ");
        New = AddNewField_VTU(New, data, Active, "Active");             # will overwrite the previous "Active" field in file
    
        # write file back to disk
        OutFileName     =   "PassiveTracers_modified.vtu"
        WriteNewPoints_VTU(DirName, OutFileName, data, New);
        
        # Store filename and timestep info (for PVD files) ------------------------------------
        Time_vec[iStep]     =   Time;
        OutNames[iStep]     =   "$DirName/$OutFileName";

    end

    # Generate PVD file from the directories & filenames ----------------------------------
    WritePVD_File("PassiveTracers_modified.pvd", OutNames, Time_vec);

end # end of postprocess passive tracers
#-------------------------------------------------------------------------------------------


end # end of module PassiveTracers

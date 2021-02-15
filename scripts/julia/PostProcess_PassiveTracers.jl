# This processes all Timestep directories in the current LaMEM 
# directory and creates visualizations of them all
using Printf

# Read LaMEM routines
include("../../../scripts/julia/ReadLaMEM_Timestep.jl")

# define a macro that finds all directories with a certain pattern
searchdir(path,key) = filter(x->occursin(key,x), readdir(path))

if length(ARGS)>0
    FileName    =   ARGS[1];            # command line argument
    FileName    = "$FileName.pvtu"
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
    global minZ, OutNames, Time_vec, ActivationAge, MeltingEvents, prevMeltFrac, MeltFrac;
    global maxP, maxT, ZirconAge;

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

    if iStep==1
        minZ            =   zCoord;
        maxP            =   P;
        maxT            =   T;
        ActivationAge   =   zCoord.*0;                              # will store the time that a particle is activated
        ZirconAge       =   zCoord.*0;
        MeltingEvents   =   zCoord.*0;                              # how many melting events did this particle see?  
        prevMeltFrac    =   zCoord.*0;                              # stores the melt fraction of the previous timestep
        OutNames        =   Array{String}(undef,    nMax);
        Time_vec        =   Array{Float64}(undef,   nMax);
    end 

    # Determine minimum coordinate of active particles
    id          =   findall(x->x==1, Active);
    if sizeof(id)>0
        minZ[id]   =   min.(minZ[id],zCoord[id]);                   # obtain minimum of evert point in array
    end

    # Max P & T
    id          =   findall(x->x==1, Active);
    if sizeof(id)>0
        maxP[id]   =   max.(maxP[id],P[id]);                  
        maxT[id]   =   max.(maxT[id],T[id]);                  
    end

    # Determine the time that a particle is activated
    id          =   findall(x->x==0, Active);
    if sizeof(id)>0
       ActivationAge[id]   .=   Time;                              # this will be updated until it becomes active
    end
    ZirconAge .= Time .- ActivationAge;
 
    for i=1:length(MeltFrac)
        if ( (prevMeltFrac[i]==0.0) && (MeltFrac[i]>0.02) && ( (Time - ActivationAge[i]) > 10 ))
            MeltingEvents[i] = MeltingEvents[i] + 1.0;  
        end
    end

    #minZ_coord = minimum(minZ);
    #print("minimum Z coordinate = $minZ_coord \n")

    # Save data to new files --------------------------------------------------------------
   
    # Add data to object
    New = [];
    New = AddNewField_VTU(New, data, zCoord,        "zCoord");        
    New = AddNewField_VTU(New, data, minZ,          "minZ");
    New = AddNewField_VTU(New, data, maxP,          "maxP");
    New = AddNewField_VTU(New, data, maxT,          "maxT");
    New = AddNewField_VTU(New, data, Active,        "Active");              # will overwrite the previous "Active" field in file
    New = AddNewField_VTU(New, data, ActivationAge, "ActivationAge");       # Shows when a particle became active, so we see how long they have been around
    New = AddNewField_VTU(New, data, ZirconAge,     "ZirconAge");           # How long ago did the first zircon form?
    New = AddNewField_VTU(New, data, MeltingEvents, "MeltingEvents");       # Shows the # of melting episodes a particle had
    New = AddNewField_VTU(New, data, MeltFrac,      "MeltFrac");            # 
    New = AddNewField_VTU(New, data, prevMeltFrac,  "prevMeltFrac");        # 
    


    # write file back to disk
    OutFileName     =   "PassiveTracers_modified.vtu"
    WriteNewPoints_VTU(DirName, OutFileName, data, New);
    
    # Store filename and timestep info (for PVD files) ------------------------------------
    Time_vec[iStep]     =   Time;
    OutNames[iStep]     =   "$DirName/$OutFileName";


    prevMeltFrac        =   MeltFrac;                                            # store melt fraction of previous timestep
end

# Generate PVD file from the directories & filenames ----------------------------------
WritePVD_File("PassiveTracers_modified.pvd", OutNames, Time_vec);



